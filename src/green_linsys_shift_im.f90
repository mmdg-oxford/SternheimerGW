  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
SUBROUTINE green_linsys_shift_im (green, xk1, iw0, mu, iq, nwgreen)
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, iunigk
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc
  USE constants,            ONLY : degspin, pi, tpi, RYTOEV, eps8
  USE cell_base,            ONLY : tpiba2
  USE ener,                 ONLY : ef
  USE klist,                ONLY : xk, wk, nkstot
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, g2kin, et, ecutwfc
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE control_gw,           ONLY : rec_code, niter_gw, nmix_gw, tr2_gw, &
                                   alpha_pv, lgamma, lgamma_gamma, convt, &
                                   nbnd_occ, alpha_mix, ldisp, rec_code_read, &
                                   where_rec, current_iq, ext_recover, &
                                   eta, tr2_green, maxter_green, prec_shift,&
                                   multishift
  USE nlcc_gw,              ONLY : nlcc_any
  USE units_gw,             ONLY : iuwfc, lrwfc, iuwfcna, iungreen, lrgrn
  USE eqv,                  ONLY : evq, eprectot
  USE qpoint,               ONLY : xq, npwq, igkq, nksq, ikks, ikqs
  USE disp,                 ONLY : nqs, x_q
  USE freq_gw,              ONLY : fpol, fiu, nfs, nfsmax, wgreen, deltaw, w0pmw
  USE gwsigma,              ONLY : sigma_c_st, ecutsco, ecutprec
  USE gvect,                ONLY : g, ngm
  USE mp,                   ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_images,            ONLY : nimage, my_image_id, intra_image_comm,   &
                                   me_image, nproc_image, inter_image_comm
  USE mp_global,            ONLY : nproc_pool_file, &
                                   nproc_bgrp_file, nproc_image_file
  USE mp_bands,             ONLY : nproc_bgrp, ntask_groups
  USE mp_world,             ONLY : nproc, mpime
  USE mp_pools,             ONLY : inter_pool_comm

  IMPLICIT NONE 

 !should be freq blocks...
 !COMPLEX(DP) :: gr_A_shift(npwx, nwgreen)
  COMPLEX(DP), ALLOCATABLE :: gr_A_shift(:,:)

  COMPLEX(DP) :: gr_A(npwx, 1), rhs(npwx , 1)
  COMPLEX(DP) :: gr(npwx, 1), ci, cw 
  COMPLEX(DP) :: green(sigma_c_st%ngmt, sigma_c_st%ngmt, nwgreen)
  COMPLEX(DP), ALLOCATABLE :: etc(:,:)

  REAL(DP) :: dirac, x, delta, support
  REAL(DP) :: k0mq(3) 
  REAL(DP) :: w_ryd(nwgreen)
  REAL(DP) , allocatable :: h_diag (:,:)
  REAL(DP)               :: eprec_gamma
  REAL(DP) :: thresh, anorm, averlt, dr2, sqrtpi
  REAL(DP) :: tr_cgsolve = 1.0d-4
  REAL(DP) :: ehomo, elumo, mu

  INTEGER :: nwgreen
  INTEGER :: iw, igp, iw0
  INTEGER :: iq, ik0
  INTEGER :: rec0, n1, gveccount
  INTEGER, ALLOCATABLE      :: niters(:)
  INTEGER :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             ltaver,     & ! average counter
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ig,         & ! counter on G vectors
             ndim,       & ! integer actual row dimension of dpsi
             is,         & ! counter on spin polarizations
             nt,         & ! counter on types
             na,         & ! counter on atoms
             nrec, nrec1,& ! the record number for dvpsi and dpsi
             ios,        & ! integer variable for I/O control
             mode          ! mode index
    INTEGER  :: igkq_ig(npwx) 
    INTEGER  :: igkq_tmp(npwx) 
    INTEGER  :: counter
    INTEGER  :: igstart, igstop, ngpool, ngr, igs, ngvecs

  REAL(DP) :: gam(3)
  REAL(DP) :: xk1(3)

  LOGICAL :: conv_root
  EXTERNAL cg_psi, ch_psi_all_green

    ALLOCATE  (h_diag (npwx, 1))
    ALLOCATE  (etc(nbnd_occ(1), nkstot))
    if(multishift) ALLOCATE(gr_A_shift(npwx, nwgreen))

    ci = (0.0d0, 1.0d0)
!Convert freq array generated in freqbins into rydbergs.
    do iw =1, nwgreen
       w_ryd(iw) = w0pmw(1,iw)/RYTOEV
    enddo
    call start_clock('greenlinsys')
    where_rec='no_recover'
!This should ensure the Green's fxn has the correct -\delta for \omega <!\epsilon_{F}:
!This smooths out variations and I think makes sense
   ikq = iq
   call gk_sort(x_q(1,ikq), ngm, g, ( ecutwfc / tpiba2 ),&
                 npw, igk, g2kin )
   npwq = npw
!Need a loop to find all plane waves below ecutsco when igkq takes us outside of this sphere.
!igkq_tmp is gamma centered index up to ngmsco,
!igkq_ig  is the linear index for looping up to npwq.
    counter = 0
    igkq_tmp(:) = 0
    igkq_ig(:)  = 0 
    do ig = 1, npwx
       if((igkq(ig).le.sigma_c_st%ngmt).and.((igkq(ig)).gt.0)) then
           counter = counter + 1
           igkq_tmp (counter) = igkq(ig)
           igkq_ig  (counter) = ig
       endif
    enddo
!Now the G-vecs up to the correlation cutoff have been divided between pools.
!Calculates beta functions (Kleinman-Bylander projectors), with
!structure factor, for all atoms, in reciprocal space
    call init_us_2 (npwq, igkq, x_q (1, ikq), vkb)
    do ig = 1, npwq
       g2kin (ig) = ((x_q (1,ikq) + g (1, igkq(ig) ) ) **2 + &
                     (x_q (2,ikq) + g (2, igkq(ig) ) ) **2 + &
                     (x_q (3,ikq) + g (3, igkq(ig) ) ) **2 ) * tpiba2
    enddo
    green  = (0.0d0, 0.0d0)
    h_diag = 0.d0
    if(multishift) then
      do ig = 1, npwq
         !if(g2kin(ig).le.ecutprec) then
           h_diag(ig,1) =  1.0d0
         !else
         !if(prec_shift) then
         !h_diag(ig,1)= 1.d0/max(1.0d0, g2kin(ig)/(eprectot(nbnd_occ(1),ikq)))
         !else
         !h_diag(ig,1) =  1.0d0
         !endif
         !endif
      enddo
    else
      do ig = 1, npwq
         h_diag(ig,1)= 1.d0/max(1.0d0, g2kin(ig)/(eprectot(nbnd_occ(1),ikq)))
      enddo
    endif
!On first frequency block we do the seed system with BiCG:
    gveccount = 1
    if(multishift) gr_A_shift = (0.0d0, 0.d0)
    call para_img(counter, igstart, igstop)
!allocate list to keep track of the number of residuals for each G-vector:
    ngvecs = igstop-igstart + 1
    if(.not.allocated(niters)) ALLOCATE(niters(ngvecs))
    niters(:) = 0
    do ig = igstart, igstop
!Doing Linear System with Wavefunction cutoff (full density) for each perturbation. 
          if(multishift) then
             rhs(:,:)  = (0.0d0, 0.0d0)
             rhs(igkq_ig(ig), 1) = -(1.0d0, 0.0d0)
             gr_A(:,:) = (0.0d0, 0.0d0)
             lter = 0
             etc(:, :) = CMPLX( 0.0d0, 0.0d0, kind=DP)
             cw = CMPLX( 0.0d0, 0.0d0, kind=DP) 
             conv_root = .true.
             anorm = 0.0d0
             call cbcg_solve_green(ch_psi_all_green, cg_psi, etc(1,ikq), rhs, gr_A, h_diag,   &
                                   npwx, npwq, tr2_green, ikq, lter, conv_root, anorm, 1, npol, &
                                   cw , niters(gveccount), .true.)
             call green_multishift_im(npwx, npwq, nwgreen, niters(gveccount), 1, w_ryd(1), mu, gr_A_shift)
             if (niters(gveccount).ge.maxter_green) then
                   WRITE(1000+mpime, '(5x,"Gvec: ", i4)') ig
                   gr_A_shift(:,:) = dcmplx(0.0d0,0.0d0)
             endif
             do iw = 1, nwgreen
               do  igp = 1, counter
                   green (igkq_tmp(ig), igkq_tmp(igp),iw) = green (igkq_tmp(ig), igkq_tmp(igp),iw) + &
                                                            gr_A_shift(igkq_ig(igp),iw)
               enddo
             enddo
             gveccount = gveccount + 1
          else if(.not.multishift) then
             do iw = 1, nwgreen
                rhs(:,:)  = (0.0d0, 0.0d0)
                rhs(igkq_ig(ig), 1) = -(1.0d0, 0.0d0)
                gr_A(:,:) = (0.0d0, 0.0d0)
                lter = 0
                etc(:, :) = CMPLX( 0.0d0, 0.0d0, kind=DP)
                cw = CMPLX( mu, w_ryd(iw), kind=DP) 
                conv_root = .true.
                anorm = 0.0d0
                call cbcg_solve(ch_psi_all_green, cg_psi, etc(1,1), rhs, gr_A, h_diag,   &
                                npwx, npwq, tr2_green, ikq, lter, conv_root, anorm, 1, npol, &
                                cw, maxter_green, .true.)
                if (lter.ge.maxter_green) then
                   WRITE(1000+mpime, '(5x,"Gvec: ", i4)') ig
                   gr_A = dcmplx(0.0d0,0.0d0)
                endif
                do igp = 1, counter
                   green (igkq_tmp(ig), igkq_tmp(igp),iw) = green (igkq_tmp(ig), igkq_tmp(igp),iw) + &
                                                            gr_A(igkq_ig(igp),1)
                enddo
             enddo
             gveccount = gveccount + 1
          endif
    enddo !igstart
#ifdef __PARA
    call mp_barrier(inter_image_comm)
    call mp_sum(green, inter_image_comm)
#endif __PARA

if(allocated(niters))     deallocate(niters)
if(allocated(h_diag))     deallocate(h_diag)
if(allocated(etc))        deallocate(etc)
if(allocated(gr_A_shift)) deallocate(gr_A_shift)

call stop_clock('greenlinsys')
return
END SUBROUTINE green_linsys_shift_im
