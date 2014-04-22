!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE gwq_init()
!----------------------------------------------------------------------------
  !
  !   THIS ROUTINE INITIALIZES K.E. = (k+G)**2 and VKB (kleinman-bylander projectors).
  !     "This subroutine computes the quantities necessary to describe the
  !     local and nonlocal pseudopotential in the GW program."
  !
  !     In detail it computes:
  !     a) The local potential at G-G'. Needed for the part of the dynamic
  !        matrix independent of deltapsi.

  !     b) The local potential at q+G-G'. Needed for the second
  !        second part of the dynamical matrix.

  !     c) The D coefficients for the US pseudopotential or the E_l parame
  !        of the KB pseudo. In the US case it prepares also the integrals
  !        qrad and qradq which are needed for computing Q_nm(G) and
  !        Q_nm(q+G)

  !     d) The functions vkb(k+G) needed for the part of the dynamical matrix
  !        independent of deltapsi.

  !     e) The becp functions for the k points

  !     e') The derivative of the becp term with respect to a displacement

  !     HL the only part of this routine that is still relevant for GW calculations
  !     is the following:

  !     0) initializes the structure factors

  !     a0) compute rhocore for each atomic-type if needed for nlcc
   
  !     In particular we need to make sure that this vkb (k+q+G) is appropriate
  !     f) The functions vkb(k+q+G), needed for the linear system and the
  !        second part of the dynamical matrix.
  !
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : bg, tpiba, tpiba2, omega
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE becmod,               ONLY : calbec
  USE constants,            ONLY : eps8, tpi
  USE gvect,                ONLY : g, ngm, ecutwfc
  USE klist,                ONLY : xk
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunigk
  USE atom,                 ONLY : msh, rgrid
  USE vlocal,               ONLY : strf
  USE spin_orb,             ONLY : lspinorb
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : noncolin, npol
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf
  USE eqv,                  ONLY : vlocq, evq, eprec
  USE gwus,                 ONLY : becp1, alphap, dpqq, dpqq_so
  USE nlcc_gw,              ONLY : nlcc_any
  USE control_gw,           ONLY : trans, zue, epsil, lgamma, all_done, nbnd_occ
  USE units_gw,             ONLY : lrwfc, iuwfc
  USE qpoint,               ONLY : xq, igkq, npwq, nksq, eigqts, ikks, ikqs
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum

  !
  IMPLICIT NONE
  !
  ! ... local variables

  INTEGER :: nt, ik, ikq, ipol, ibnd, ikk, na, ig, irr, imode0
    ! counter on atom types
    ! counter on k points
    ! counter on k+q points
    ! counter on polarizations
    ! counter on bands
    ! index for wavefunctions at k
    ! counter on atoms
    ! counter on G vectors
  REAL(DP) :: arg
    ! the argument of the phase
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)
    ! used to compute alphap
  COMPLEX(DP), EXTERNAL :: zdotc
  !
  !
  IF (all_done) RETURN
  !
  CALL start_clock( 'gwq_init' )
  !
  ALLOCATE( aux1( npwx*npol, nbnd ) )    

!HL structure factors

  DO na = 1, nat
     arg = ( xq(1) * tau(1,na) + &
             xq(2) * tau(2,na) + &
             xq(3) * tau(3,na) ) * tpi
     eigqts(na) = CMPLX( COS( arg ), - SIN( arg ) ,kind=DP)
  END DO

  !                 
  !
  ! ... a0) compute rhocore for each atomic-type if needed for nlcc
  !

  IF ( nlcc_any ) WRITE(6, '(" NLCC ")')
  IF ( nlcc_any ) CALL set_drhoc( xq )
  !
  ! ... b) the fourier components of the local potential at q+G
  !
  vlocq(:,:) = 0.D0
  !
  DO nt = 1, ntyp
     !
     IF (upf(nt)%tcoulombp) then
        CALL setlocq_coul ( xq, upf(nt)%zp, tpiba2, ngm, g, omega, vlocq(1,nt) )
     ELSE
        CALL setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
                   upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, &
                   vlocq(1,nt) )
     END IF
     !
  END DO
  !


  IF ( nksq > 1 ) REWIND( iunigk )
  !
  DO ik = 1, nksq
     !
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     !
     IF ( lsda ) current_spin = isk( ikk )
     !
     ! ... g2kin is used here as work space
     ! ... sorts k+g in order of increasing magnitude, up to ecut
     ! ... NB: this version should yield the same ordering for different ecut
     ! ...     and the same ordering in all machines
     !HL

     ! write(6,*) nksq, ikk, ikq 
     ! write(6,*) ngm, npw

     CALL gk_sort( xk(1,ikk), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )

     !
     ! ... if there is only one k-point evc, evq, npw, igk stay in memory
     !

     IF ( nksq > 1 ) WRITE( iunigk ) npw, igk

     !

     IF ( lgamma ) THEN

     !
        npwq = npw
     !
     ELSE   
     !
        CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ), &
                      npwq, igkq, g2kin )
     !
        IF ( nksq > 1 ) WRITE( iunigk ) npwq, igkq
     !
        IF ( ABS( xq(1) - ( xk(1,ikq) - xk(1,ikk) ) ) > eps8 .OR. &
             ABS( xq(2) - ( xk(2,ikq) - xk(2,ikk) ) ) > eps8 .OR. &
             ABS( xq(3) - ( xk(3,ikq) - xk(3,ikk) ) ) > eps8 ) THEN
           WRITE( stdout,'(/,5x,"k points #",i6," and ", &
                  & i6,5x," total number ",i6)') ikk, ikq, nksq
           WRITE( stdout, '(  5x,"Expected q ",3f10.7)')(xq(ipol), ipol=1,3)
           WRITE( stdout, '(  5x,"Found      ",3f10.7)')((xk(ipol,ikq) &
                                                -xk(ipol,ikk)), ipol = 1, 3)
           CALL errore( 'gwq_init', 'wrong order of k points', 1 )
        END IF
        !
     END IF

     ! ... d) The functions vkb(k+G)
     !
     CALL init_us_2( npw, igk, xk(1,ikk), vkb )
     !
     ! ... read the wavefunctions at k
     !
     CALL davcio( evc, lrwfc, iuwfc, ikk, -1 )

     !
     ! ... e) we compute the becp terms which are used in the rest of
     ! ...    the code
     !

     CALL calbec (npw, vkb, evc, becp1(ik) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!ALL of this between !!!-!!! can be deleted, has to do with phonon perturbations.!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     ! ... e') we compute the derivative of the becp term with respect to an
     !         atomic displacement i.e. equation B9
     !
     ! DO ipol = 1, 3
     !    aux1=(0.d0,0.d0)
     !    DO ibnd = 1, nbnd
     !       DO ig = 1, npw
     !          aux1(ig,ibnd) = evc(ig,ibnd) * tpiba * ( 0.D0, 1.D0 ) * & 
     !                          ( xk(ipol,ikk) + g(ipol,igk(ig)) )
     !       END DO
     !       IF (noncolin) THEN
     !          DO ig = 1, npw
     !             aux1(ig+npwx,ibnd)=evc(ig+npwx,ibnd)*tpiba*(0.D0,1.D0)*& 
     !                       ( xk(ipol,ikk) + g(ipol,igk(ig)) )
     !          END DO
     !       END IF
     !    END DO
     !    CALL calbec (npw, vkb, aux1, alphap(ipol,ik) )
     ! END DO

     !

     IF ( .NOT. lgamma ) &
        CALL davcio( evq, lrwfc, iuwfc, ikq, -1 )

     !
     ! diagonal elements of the unperturbed Hamiltonian, 
     ! needed for preconditioning
     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     do ig = 1, npwq
        g2kin (ig) = ( (xk (1,ikq) + g (1, igkq(ig)) ) **2 + &
                       (xk (2,ikq) + g (2, igkq(ig)) ) **2 + &
                       (xk (3,ikq) + g (3, igkq(ig)) ) **2 ) * tpiba2
     enddo

     aux1 = (0.d0,0.d0)

     DO ig = 1, npwq
        aux1 (ig,1:nbnd_occ(ikk)) = g2kin (ig) * evq (ig, 1:nbnd_occ(ikk))
     END DO


     IF (noncolin) THEN
        DO ig = 1, npwq
           aux1 (ig+npwx,1:nbnd_occ(ikk)) = g2kin (ig)* &
                                  evq (ig+npwx, 1:nbnd_occ(ikk))
        END DO
     END IF

     DO ibnd=1,nbnd_occ(ikk)
        eprec (ibnd,ik) = 1.35d0 * zdotc(npwx*npol,evq(1,ibnd),1,aux1(1,ibnd),1)
     END DO
     !
  END DO

!#ifdef __PARA
!     CALL mp_sum ( eprec, intra_pool_comm )
!#endif

  !
  DEALLOCATE( aux1 )

  CALL stop_clock( 'gwq_init' )
  !
  RETURN
  !
END SUBROUTINE gwq_init
