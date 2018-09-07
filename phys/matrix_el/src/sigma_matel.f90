!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! 
! Copyright (C) 2010 - 2018
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! SternheimerGW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! SternheimerGW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with SternheimerGW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
MODULE sigma_matel

  IMPLICIT NONE

CONTAINS

SUBROUTINE matrix_element(ikpt, calc)

  USE buffers,              ONLY : get_buffer, close_buffer
  USE buiol,                ONLY : buiol_check_unit
  USE cell_base,            ONLY : tpiba2
  USE constants,            ONLY : RYTOEV, eps14
  USE container_interface,  ONLY : element_type, no_error
  USE control_gw,           ONLY : do_imag
  USE control_lr,           ONLY : lgamma
  USE disp,                 ONLY : xk_kpoints, num_k_pts
  USE driver,               ONLY : calculation
  USE ener,                 ONLY : ef
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE freq_gw,              ONLY : nwsigma, nwsigwin
  USE freqbins_module,      ONLY : freqbins_type
  USE gvect,                ONLY : ngm, g, gl, igtongl
  USE gvecw,                ONLY : ecutwfc
  USE gw_data,              ONLY : var_exch, var_corr
  USE gwsigma,              ONLY : nbnd_sig, corr_conv, exch_conv, ecutsco, ecutsex
  USE io_files,             ONLY : diropn 
  USE io_global,            ONLY : stdout, meta_ionode
  USE kinds_gw,             ONLY : i8b
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk, lgauss
  USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
  USE qpoint,               ONLY : npwq
  USE scf,                  ONLY : rho, rho_core, rhog_core, scf_type, v
  USE sigma_expect_mod,     ONLY : sigma_expect_after_wavef_ordering
  USE timing_module,        ONLY : time_matel
  USE units_gw,             ONLY : iuwfc, lrwfc
  USE wavefunctions,        ONLY : evc
  USE wvfct,                ONLY : nbnd, npw, npwx, g2kin

IMPLICIT NONE

  !> the data of the GW calculation
  INTEGER, INTENT(IN) :: ikpt
  TYPE(calculation), INTENT(INOUT) :: calc
  TYPE(element_type) element

  COMPLEX(DP), ALLOCATABLE  :: sigma_band_con(:,:,:), sigma_band_x(:,:,:), sigma_band_c(:,:,:)
  COMPLEX(DP)               :: psic(dffts%nnr), vpsi(ngm)
  COMPLEX(DP)               :: ZdoTC, vxc(nbnd_sig,nbnd_sig)
  REAL(DP)                  :: vtxc, etxc
  INTEGER                   :: igk(npwx), ikq
  INTEGER                   :: ig, ibnd, jbnd, ir, ierr
  INTEGER                   :: ng
  INTEGER                   :: sigma_c_ngm, sigma_x_ngm

  !> energy of the highest occupied state
  REAL(dp) ehomo

  !> energy of the lowest unoccupied state
  REAL(dp) elumo

  !> the chemical potential
  REAL(dp) mu

  !> real constant of 0.5
  REAL(dp),    PARAMETER :: half = 0.5_dp

  !> complex constant of 0
  COMPLEX(dp), PARAMETER :: zero = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

  CALL start_clock(time_matel)

  IF ( .NOT. meta_ionode ) RETURN

  nbnd   = nbnd_sig 
  lgamma = .true.

  !
  ! define the chemical potential
  !
  IF (.NOT.lgauss) THEN
    !
    ! for semiconductors choose the middle of the gap
    CALL get_homo_lumo(ehomo, elumo)
    mu = half * (ehomo + elumo)
    !
  ELSE
    !
    ! for metals set it to the Fermi energy
    mu = ef
    !
  END IF

  ! check for gamma point
  IF(ALL(ABS(xk_kpoints(:,ikpt)) <= eps14)) THEN
     ikq = 1
  ELSE
     ikq = 2
  END IF

  WRITE(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikpt, xk_kpoints(:,ikpt)
  WRITE(stdout,'(/4x,"k0(",i3," ) = (", 3f7.3, " )")') ikq, xk(:,ikq)
  ! set matrix elements to 0
  vxc          = 0
  ! create map to G ordering at current k-point
  CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ),&
                npw, igk, g2kin )
  npwq = npw
  ! read wave functions at current k-point
  CALL get_buffer (evc, lrwfc, iuwfc, ikq)
  evc(npw + 1:, :) = zero

  !
  ! expectation value of V_xc
  !

  ! generate v_xc(r) in real space:
  v%of_r = 0
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  ! loop over all bands
  DO jbnd = 1, nbnd_sig

    ! extract right wave function
    psic = 0
    DO ig = 1, npwq
      ! map G-vectors according to smooth igk
      psic(dffts%nl(igk(ig))) = evc(ig, jbnd)
    END DO

    ! FFT wave function to real space
    CALL invfft ('Wave', psic(:), dffts)

    ! product of wave function and potential
    DO ir = 1, dfftp%nnr
      psic (ir) = psic(ir) * v%of_r (ir,1)
    END DO

    ! FFT back to reciprocal space
    CALL fwfft ('Wave', psic(:), dffts)

    ! reorder back to original order
    DO ig = 1, npwq
      vpsi(ig) = psic(dffts%nl(igk(ig)))
    END DO
 
    ! multiply with left wave function 
    DO ibnd = 1, nbnd_sig
      vxc(ibnd,jbnd) = ZdoTC (npwq, evc (1, ibnd), 1, vpsi, 1)
    END DO

  END DO ! bands

  !
  ! expectation value of Sigma_x
  !
  ! read exchange self energy
  element%variable = var_exch
  element%access_index = ikpt
  CALL calc%data%read_element(element, ierr)

  IF (ierr == no_error) THEN

    ! sanity check
    IF (ABS(exch_conv - ecutsex) < eps14 .OR. &
        ABS(exch_conv) < eps14) THEN
      sigma_x_ngm = calc%grid%exch_fft%ngm
    ELSE IF((exch_conv < ecutsex) .AND. (exch_conv > 0.0)) THEN
      DO ng = 1, ngm
         IF ( gl( igtongl (ng) ) <= (exch_conv/tpiba2)) sigma_x_ngm = ng
      END DO
    ELSE
      CALL errore("sigma_matel", "Exch Conv must be greater than zero and less than ecutsex", 1)
    END IF

    ! evaluate matrix elements for exchange
    CALL sigma_expect_after_wavef_ordering(calc%data%exch,evc,igk,sigma_band_x)

  ELSE

    ALLOCATE(sigma_band_x(nbnd_sig, nbnd_sig, 1))
    sigma_band_x = zero

  END IF

  !
  ! expectation value of Sigma_c:
  !

  ! open file containing correlation part of sigma
  element%variable = var_corr
  element%access_index = ikpt
  CALL calc%data%read_element(element, ierr)

  IF (ierr == no_error) THEN

    ! For convergence tests corr_conv can be set at input lower than ecutsco.
    ! This allows you to calculate the correlation energy at lower energy cutoffs
    IF (ABS(corr_conv - ecutsco) < eps14) THEN
      sigma_c_ngm = calc%grid%corr_fft%ngm
    ELSE IF(corr_conv < ecutsco .AND. corr_conv > 0.0) THEN
      DO ng = 1, ngm
        IF (gl( igtongl (ng) ) <= (corr_conv/tpiba2)) sigma_c_ngm = ng
      END DO
    ELSE
      CALL errore("sigma_matel", "Corr Conv must be greater than zero and less than ecutsco", 1)
    END IF

    ! evaluate expectation value of wave function
    CALL sigma_expect_after_wavef_ordering(calc%data%corr(:,:,:,1),evc,igk,sigma_band_c)

  ELSE

    ALLOCATE(sigma_band_c(nbnd_sig, nbnd_sig, nwsigma))
    sigma_band_c = zero

  END IF

  !
  ! analytic continuation from imaginary frequencies to real ones
  ! 
  IF (do_imag) THEN
    ! We can set arbitrary \Sigma(\omega) energy windows with analytic continuation:
    ALLOCATE (sigma_band_con(nbnd_sig, nbnd_sig, nwsigwin))
    ! print selfenergy on the imaginary axis.
    CALL print_matel_im(ikq, vxc, sigma_band_x, sigma_band_c, AIMAG(calc%freq%sigma), nwsigma)
    ! do analytic continuation and print selfenergy on the real axis.
    sigma_band_con(:,:,:) = zero
    CALL sigma_pade(sigma_band_c, sigma_band_con, mu, AIMAG(calc%freq%sigma), calc%freq%window, nwsigwin)
    CALL print_matel(ikq, vxc, sigma_band_x, sigma_band_con, calc%freq%window, nwsigwin)
   deallocate(sigma_band_con)
  ELSE
    ! print sigma on real axis
    CALL print_matel(ikq, vxc, sigma_band_x, sigma_band_c, REAL(calc%freq%sigma) + mu, nwsigma)
  END IF

  IF (ALLOCATED(calc%data%exch)) DEALLOCATE(calc%data%exch)
  IF (ALLOCATED(calc%data%corr)) DEALLOCATE(calc%data%corr)
  CALL stop_clock(time_matel)

END SUBROUTINE matrix_element

END MODULE sigma_matel
