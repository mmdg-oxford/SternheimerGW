!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
!> G TIMES W PRODUCT
!! On the real frequency axis the Green's function
!! Requires access to the occupied staets \psi_{n\k}(\r)
!! Therefore we use the same tricks as for sigma_exch.f90
SUBROUTINE sigma_c_re(ik0) 

  USE cell_base,         ONLY : tpiba2, tpiba, omega, alat, at,bg
  USE constants,         ONLY : e2, fpi, RYTOEV, tpi, eps8, pi
  USE control_flags,     ONLY : noinv
  USE control_gw,        ONLY : lgamma, eta, godbyneeds, padecont, cohsex, modielec, &
                                trunc_2d, tmp_dir_coul, output
  USE disp,              ONLY : nqs, nq1, nq2, nq3, wq, x_q, xk_kpoints
  USE eqv,               ONLY : evq
  USE expand_igk_module, ONLY : expand_igk
  USE freq_gw,           ONLY : fpol, fiu, nfs, nfsmax, &
                                nwcoul, nwgreen, nwalloc, nwsigma, wtmp, wcoul, &
                                wgreen, wsigma, wsigmamin, wsigmamax, &
                                deltaw, wcoulmax, ind_w0mw, ind_w0pw, &
                                w0pmw, wgtcoul
  USE gvect,             ONLY : g, ngm, nl
  USE gwsigma,           ONLY : sigma_c_st
  USE io_files,          ONLY : prefix, tmp_dir
  USE io_global,         ONLY : stdout, ionode_id, ionode, meta_ionode
  USE kinds,             ONLY : DP
  USE kinds_gw,          ONLY : i8b
  USE klist,             ONLY : wk, xk, nkstot, nks
  USE lr_symm_base,      ONLY : nsymq, invsymq, gi, gimq, irgq, irotmq, minus_q
  USE lsda_mod,          ONLY : nspin
  USE mp,                ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_global,         ONLY : mp_global_end
  USE mp_images,         ONLY : nimage, my_image_id, intra_image_comm,   &
                                me_image, nproc_image, inter_image_comm
  USE mp_pools,          ONLY : inter_pool_comm
  USE mp_world,          ONLY : nproc, mpime
  USE output_mod,        ONLY : filcoul
  USE qpoint,            ONLY : xq, npwq, nksq, ikks, ikqs
  USE sigma_io_module,   ONLY : sigma_io_write_c
  USE symm_base,         ONLY : nsym, s, time_reversal, t_rev, ftau, invs, nrot
  USE timing_module,     ONLY : time_sigma_c, time_GW_product, time_sigma_setup, &
                                time_sigma_comm, time_sigma_io
  USE units_gw,          ONLY : iuncoul, iungreen, iunsigma, lrsigma, lrcoul, lrgrn, iuwfc, lrwfc
  USE wvfct,             ONLY : nbnd, npw, npwx, g2kin

  IMPLICIT NONE

  COMPLEX(DP)         :: ci, czero
  COMPLEX(DP)         :: phase
  COMPLEX(DP)         :: aux (sigma_c_st%dfftt%nnr)
!Sigma arrays
  COMPLEX (DP), ALLOCATABLE :: sigma(:,:,:)
  COMPLEX (DP), ALLOCATABLE :: sigma_g(:,:,:)
!Pade arrays
  COMPLEX(DP), ALLOCATABLE :: z(:), u(:), a(:)
!W arrays 
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_g_R (:,:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul_pade_g (:,:)
  COMPLEX(DP), ALLOCATABLE :: scrcoul(:,:)
!G arrays:
  COMPLEX(DP), ALLOCATABLE :: greenf_g(:,:,:), greenfr(:,:)
  COMPLEX(DP) :: cprefac
  COMPLEX(DP), ALLOCATABLE  ::  eigv(:,:)
!v array
!COMPLEX(DP), ALLOCATABLE ::  barcoul(:,:), barcoulr(:,:), barcoul_R(:,:)
  REAL(DP) :: qg2, qg, qxy, qz
  REAL(DP) :: w_ryd(nwcoul), w_rydsig(nwsigma)
  REAL(DP) :: xq_ibk(3), xq_ibz(3)
!q-vector of coulomb potential:
  REAL(DP) :: xq_coul(3)
  REAL(DP) :: rcut, spal
!CHECK FOR NAN's
  REAL(DP)     :: ar, ai
!For dirac delta fxn.
  REAL(DP)     :: dirac, x, support, zcut
  REAL(DP) :: ehomo, elumo, mu
!FREQUENCY GRIDS/COUNTERS
  INTEGER  :: iwim, iw, ikq
  INTEGER  :: iw0, iw0mw, iw0pw
!COUNTERS
  INTEGER :: ig, igp, irr, icounter, ir, irp
  INTEGER :: iqs, nkr
  INTEGER :: iq, ipol, iqrec
  INTEGER :: ikmq, ik0, ik, nkpool
  INTEGER :: rec0, ios
  INTEGER :: counter, ierr
  INTEGER :: inversym, screening
!SYMMETRY
  INTEGER               :: isym, jsym, isymop, nig0
  INTEGER, ALLOCATABLE  :: gmapsym(:,:)
!For G^NA
  INTEGER     :: igkq_ig(npwx) 
  INTEGER     :: igkq_tmp(npwx) 
  INTEGER     :: ss(3,3)
  INTEGER     :: ibnd
  INTEGER     :: iw0start, iw0stop
!For running PWSCF need some variables 
  LOGICAL             :: pade_catch
  LOGICAL             :: found_q
  LOGICAL             :: limq, inv_q, found
!File related:
  character(len=256) :: tempfile, filename
!Complete file name
  integer(i8b) :: unf_recl
  REAL (DP)   :: xk1(3), aq(3)
  REAL(DP)    :: sxq(3,48), xqs(3,48)
  INTEGER     :: imq, isq(48), nqstar, nkpts
  INTEGER     :: i
  REAL(DP)    :: wgt(nsym), xk_un(3,nsym)
  INTEGER     :: ixk1, iqrec_old
  INTEGER     :: isym_k(nsym), nig0_k(nsym), iqrec_k(nsym)
  LOGICAL, EXTERNAL :: eqvect
  LOGICAL     :: invq_k(nsym), found_qc
  LOGICAL     :: k_cycle(nsym)
  INTEGER     :: iqcoul
  real(DP), parameter :: eps=1.e-5_dp
  logical :: found_k
  INTEGER :: ikstar, ik1

#define DIRECT_IO_FACTOR 8 

   CALL start_clock(time_sigma_c)
   CALL start_clock(time_sigma_setup)

   CALL expand_igk

! iG(W-v)
   ALLOCATE ( scrcoul_g       (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)    )
   ALLOCATE ( scrcoul_g_R     (sigma_c_st%ngmt, sigma_c_st%ngmt, nfs)    )
   ALLOCATE ( scrcoul_pade_g  (sigma_c_st%ngmt, sigma_c_st%ngmt)         )
   ALLOCATE ( greenf_g        (sigma_c_st%ngmt, sigma_c_st%ngmt, nwgreen) )
!These go on the big grid...
   ALLOCATE ( scrcoul        (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr))
   ALLOCATE ( greenfr        (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr))
!Technically only need gmapsym up to sigma_c_st%ngmt or ngmgrn...
   ALLOCATE ( gmapsym  (ngm, nrot)   )
   ALLOCATE ( eigv     (ngm, nrot)   )
!This is a memory hog...
   ALLOCATE (sigma  (sigma_c_st%dfftt%nnr, sigma_c_st%dfftt%nnr, nwsigma))
   ALLOCATE  (z(nfs), a(nfs), u(nfs))
   w_ryd(:) = wcoul(:nwcoul)/RYTOEV
   w_rydsig(:) = wsigma(:nwsigma)/RYTOEV
   WRITE(6,"( )")
   WRITE(6,'(4x,"Direct product GW for k0(",i3," ) = (",3f12.7," )")') ik0, (xk_kpoints(ipol, ik0), ipol=1,3)
   WRITE(6,"( )")
   WRITE(6,'(4x, "ngmsco, ", i4, " nwsigma, ", i4)') sigma_c_st%ngmt, nwsigma
   WRITE(6,'(4x, "nrsco, ", i4, " nfs, ", i4)') sigma_c_st%dfftt%nnr, nfs
   zcut = 0.50d0*sqrt(at(1,3)**2 + at(2,3)**2 + at(3,3)**2)*alat
   ci = (0.0d0, 1.d0)
   czero = (0.0d0, 0.0d0)
   sigma(:,:,:) = (0.0d0, 0.0d0)
   CALL gmap_sym(nrot, s, ftau, gmapsym, eigv, invs)
   IF(allocated(sigma)) THEN
     WRITE(6,'(4x,"Sigma allocated")')
   ELSE
     WRITE(6,'(4x,"Sigma too large!")')
     CALL mp_global_end()
     STOP
   ENDIF
   WRITE(6,'(4x, "nsym, nsymq, nsymbrav ", 3i4)') nsym, nsymq, nrot 
!Set appropriate weights for points in the brillouin zone.
!Weights of all the k-points are in odd positions in list.
!nksq is number of k points not including k+q.
!Every processor needs access to the files: _gw0si.coul1 and _gw0si.green1
   call mp_barrier(inter_image_comm)
#ifdef __PARA
!OPEN coulomb file (only written to by head node).
   filename = trim(prefix)//"."//trim(filcoul)//"1"
   tempfile = trim(tmp_dir_coul) // trim(filename)
   unf_recl = DIRECT_IO_FACTOR * int(lrcoul, kind=kind(unf_recl))
   open(iuncoul, file = trim(adjustl(tempfile)), iostat = ios, &
   form = 'unformatted', status = 'OLD', access = 'direct', recl = unf_recl)
#endif
  CALL para_img(nwsigma, iw0start, iw0stop)
  WRITE(6, '(5x, "nwsigma ",i4, " iw0start ", i4, " iw0stop ", i4)') nwsigma, iw0start, iw0stop
  WRITE(1000+mpime, '(5x, "nwsigma ",i4, " iw0start ", i4, " iw0stop ", i4)') nwsigma, iw0start, iw0stop
!ONLY PROCESSORS WITH K points to process: 
  WRITE(6,'(4x,"Starting Frequency Integration")')
!kpoints split between pools
  CALL get_homo_lumo (ehomo, elumo)
  mu = ehomo + 0.5d0*(elumo-ehomo)
  call mp_barrier(inter_pool_comm)
  call mp_bcast(mu, ionode_id ,inter_pool_comm)
  call mp_barrier(inter_pool_comm)
  WRITE(6,'("mu", f12.7)') mu*RYTOEV
  WRITE(1000+mpime,'("mu", f12.7)') mu*RYTOEV

  CALL stop_clock(time_sigma_setup)

DO iq = 1, nqs
   iqcoul = 1
   found_qc = .false.
   scrcoul_g(:,:,:)   = dcmplx(0.0d0, 0.0d0)

   cprefac = (deltaw/RYTOEV)*wq(iq)*(0.0d0, 1.0d0)/tpi

   if(.not.modielec) CALL davcio(scrcoul_g, lrcoul, iuncoul, iq, -1)
   CALL coulpade(scrcoul_g(1,1,1), xq(1))
   iqrec_k = 0
   isym_k  = 0
   nig0_k  = 0
   invq_k  = .false.
   DO isymop = 1, nsym
      CALL rotate(xq, aq, s, nsym, invs(isymop))
      xk1 = xk_kpoints(:,ik0) - aq(:)
      nig0 = 1
      inv_q=.false.
      call find_qG_ibz(xk1, s, iqrec, isym, nig0, found_q, inv_q)
      found_k = .false.
      do ikstar = 1, nks 
         found_k  = (abs(xk(1,ikstar) - x_q(1,iqrec)).le.eps).and. &
                    (abs(xk(2,ikstar) - x_q(2,iqrec)).le.eps).and. & 
                    (abs(xk(3,ikstar) - x_q(3,iqrec)).le.eps) 
         if (found_k) then
            ik1 = ikstar
            exit
         endif
      enddo
      k_cycle(isymop) = found_k
      iqrec_k(isymop) = ik1
      isym_k(isymop)  = isym
      nig0_k(isymop)  = nig0
      invq_k(isymop)  = inv_q
   ENDDO
   iqrec_old = 0
   DO isymop = 1, nsym
  !Same approach as in exchange this means we can
  !read wave.
     if(.not.k_cycle(isymop)) CYCLE 
     iqrec = iqrec_k (isymop) 

     ! Green's function
     CALL stop_clock(time_sigma_c)
     IF (iqrec_old /= iqrec) WRITE(1000+mpime,'("RUNNING THROUGH GREENS FUNCTION")')
     IF (iqrec_old /= iqrec) CALL green_linsys_shift_re(greenf_g(1,1,1), mu, iqrec)
     CALL stop_clock(time_sigma_c)

     iqrec_old = iqrec
     isym   = isym_k(isymop)
     nig0   = nig0_k(isymop)
     inv_q  = invq_k(isymop)
     if(inv_q) write(1000+mpime, '("Need to use time reversal")')
     write(1000+mpime, '("xk point, isym, iqrec, nig0")')
     write(1000+mpime, '(3f11.7, 3i4)') x_q(:, iqrec), isym, iqrec, nig0
!Start integration over iw +/- wcoul.
!Rotate W and initialize necessary 
!quantities for pade_continuation or godby needs.
  IF(iw0stop-iw0start+1.gt.0) THEN
     DO iw0 = iw0start, iw0stop
        DO iw = 1, nwcoul
           IF ( iw/2*2.eq.iw ) THEN
                cprefac = cprefac * 4.d0/3.d0
           ELSE
                cprefac = cprefac * 2.d0/3.d0
           ENDIF
           IF (iw.eq.1) THEN
               cprefac = (1.0d0/3.0d0)*(deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi
           ELSE IF (iw.eq.nwcoul) THEN
               cprefac = (1.0d0/3.0d0)*(deltaw/RYTOEV) * wq(iq) * (0.0d0, 1.0d0)/ tpi
           ENDIF
           CALL construct_w(scrcoul_g(1,1,1), scrcoul_pade_g(1,1), (w_ryd(iw)-w_rydsig(iw0)))

           CALL start_clock(time_GW_product)
           scrcoul = czero
           CALL fft6_c(scrcoul_pade_g(1,1), scrcoul(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isymop, +1)
           greenfr(:,:) = czero
           if(.not.inv_q) then
             CALL fft6_g(greenf_g(1,1,iw), greenfr(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isym, nig0, +1)
             sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*cprefac*greenfr(:,:)*scrcoul(:,:)
           else
             CALL fft6_g(greenf_g(1,1,iw+nwcoul), greenfr(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isym, nig0, +1)
             sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*cprefac*conjg(greenfr(:,:))*scrcoul(:,:)
           endif
           CALL stop_clock(time_GW_product)

           CALL construct_w(scrcoul_g(1,1,1), scrcoul_pade_g(1,1), (-w_rydsig(iw0)-w_ryd(iw)))

           CALL start_clock(time_GW_product)
           scrcoul = czero
           CALL fft6_c(scrcoul_pade_g(1,1), scrcoul(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isymop, +1)
           greenfr(:,:) = czero
           if(.not.inv_q) then
             CALL fft6_g(greenf_g(1,1,iw+nwcoul), greenfr(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isym, nig0, +1)
             sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*cprefac*greenfr(:,:)*scrcoul(:,:)
           else
             CALL fft6_g(greenf_g(1,1,iw), greenfr(1,1), sigma_c_st, gmapsym(1,1), eigv(1,1), isym, nig0, +1)
             sigma (:,:,iw0) = sigma (:,:,iw0) + (1.0d0/dble(nsym))*cprefac*conjg(greenfr(:,:))*scrcoul(:,:)
           endif
           CALL stop_clock(time_GW_product)

        ENDDO !on frequency convolution over w'
     ENDDO !on iw0  
  ENDIF
 ENDDO!ISYMOP
ENDDO!iq
DEALLOCATE ( gmapsym          )
DEALLOCATE ( greenfr          )
DEALLOCATE ( greenf_g         )
DEALLOCATE ( scrcoul          )
DEALLOCATE ( scrcoul_pade_g   )
DEALLOCATE ( scrcoul_g, scrcoul_g_R )
DEALLOCATE ( z,a,u )

  CALL start_clock(time_sigma_comm)

#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
  CALL mp_sum(sigma, inter_pool_comm)
  CALL mp_barrier(inter_image_comm)
  CALL mp_sum(sigma, inter_image_comm)
#endif __PARA

  CALL stop_clock(time_sigma_comm)
  CALL start_clock(time_sigma_io)

  IF (meta_ionode) THEN
    ALLOCATE ( sigma_g (sigma_c_st%ngmt, sigma_c_st%ngmt, nwsigma))
    IF(allocated(sigma_g)) THEN
       WRITE(6,'(4x,"Sigma_g allocated")')
    ELSE
       WRITE(6,'(4x,"Sigma_g too large!")')
       CALL mp_global_end()
       STOP
    ENDIF
    WRITE(6,'(4x,"Sigma in G-Space")')
    sigma_g = (0.0d0,0.0d0)
    DO iw = 1, nwsigma
       CALL fft6(sigma_g(1,1,iw), sigma(1,1,iw), sigma_c_st, -1)
    ENDDO
!Now write Sigma in G space to file. 
    CALL davcio (sigma_g, lrsigma, iunsigma, ik0, 1)
    CALL sigma_io_write_c(output%unit_sigma, ik0, sigma_g)
    WRITE(6,'(4x,"Sigma Written to File")')
    DEALLOCATE ( sigma_g  )
  ENDIF !ionode
  CALL mp_barrier(inter_image_comm)
  DEALLOCATE ( sigma  )

  CALL stop_clock(time_sigma_io)
  CALL stop_clock(time_sigma_c)

END SUBROUTINE sigma_c_re
