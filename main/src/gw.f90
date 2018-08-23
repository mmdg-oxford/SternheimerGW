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
PROGRAM gw
!-----------------------------------------------------------------------
!... This is the main driver of the SternheimerGW code.
!-----------------------------------------------------------------------
  USE close_gwq_mod,        ONLY : close_gwq
  USE control_gw,           ONLY : do_sigma_exx, do_sigma_matel, do_coulomb, &
                                   do_sigma_c, do_q0_only, output
  USE disp,                 ONLY : num_k_pts, w_of_k_start, w_of_k_stop
  USE driver,               ONLY : calculation
  USE exchange_module,      ONLY : exchange_wrapper
  USE freq_gw,              ONLY : nwsigma, nwsigwin
  USE gwsigma,              ONLY : nbnd_sig
  USE io_global,            ONLY : meta_ionode
  USE pp_output_mod,        ONLY : pp_output_open_all
  USE run_nscf_module,      ONLY : run_nscf
  USE setup,                ONLY : setup_calculation
  USE sigma_matel,          ONLY : matrix_element
  USE sigma_module,         ONLY : sigma_wrapper
  USE timing_module,        ONLY : time_setup

  IMPLICIT NONE

  INTEGER             :: ik
  LOGICAL             :: do_band, do_matel

  TYPE(calculation) calc

  CALL setup_calculation(calc)

! Calculation W
  IF(do_coulomb) CALL do_stern(calc%config_coul, calc%grid, calc%freq)
  ik = 1
  do_band  = .TRUE.
  do_matel = .TRUE.
! Calculation of CORRELATION energy \Sigma^{c}_{k}=\sum_{q}G_{k-q}{W_{q}-v_{q}}:
  IF (.NOT.do_q0_only) THEN
      DO ik = w_of_k_start, w_of_k_stop
         CALL start_clock(time_setup)
         CALL run_nscf(do_band, do_matel, ik, calc%config)
         CALL initialize_gw(.FALSE.)
         CALL stop_clock(time_setup)
         IF (do_sigma_c) CALL sigma_wrapper(ik, calc, calc%grid, calc%config_green, &
           calc%freq, calc%vcut, calc%config, calc%debug)
! Calculation of EXCHANGE energy \Sigma^{x}_{k}= \sum_{q}G_{k}{v_{k-S^{-1}q}}:
         IF (do_sigma_exx) call exchange_wrapper(ik, calc, calc%grid, calc%vcut)
! Calculation of Matrix Elements <n\k| V^{xc}, \Sigma^{x}, \Sigma^{c}(iw) |n\k>:
         IF (do_sigma_matel) then
           IF (meta_ionode .AND. ik == w_of_k_start) then         
             call pp_output_open_all(num_k_pts, nbnd_sig, nwsigwin, nwsigma, output)
           END IF
           CALL matrix_element(ik, calc)
         END IF
         CALL clean_pw_gw(.TRUE.)
      END DO
  END IF
  CALL close_gwq(.TRUE., calc%data)
  CALL stop_gw(.TRUE.)

END PROGRAM gw
