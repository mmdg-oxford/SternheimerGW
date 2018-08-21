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
MODULE setup

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_calculation(calc)
    !
    USE check_stop,        ONLY: check_stop_init
    USE control_gw,        ONLY: do_imag
    USE environment,       ONLY: environment_start
    USE freq_gw,           ONLY: nwsigma, nwsigwin, wsigmamin, wsigmamax, wcoulmax, nwcoul, &
                                 wsig_wind_min, wsig_wind_max, nwsigwin
    USE freqbins_module,   ONLY: freqbins
    USE gw_driver,         ONLY: calculation
    USE gw_opening,        ONLY: gw_opening_logo, gw_opening_message
    USE gwsigma,           ONLY: ecutsco, ecutsex
    USE mp_global,         ONLY: mp_startup
    USE sigma_grid_module, ONLY: sigma_grid
    USE timing_module,     ONLY: time_setup
    !
    TYPE(calculation), INTENT(OUT) :: calc
    CHARACTER(*), PARAMETER :: code = 'SternheimerGW'
    !
    ! Initialize MPI, clocks, print initial messages
    CALL mp_startup(start_images = .TRUE.)
    CALL gw_opening_logo()
    CALL environment_start(code)
    CALL gw_opening_message()
    ! Initialize GW calculation, read ground state information.
    ! Initialize frequency grids, FFT grids for correlation
    ! and exchange operators, open relevant GW-files.
    CALL start_clock(time_setup)
    CALL gwq_readin(calc%config_coul, calc%config_green, calc%freq, calc%vcut, calc%debug)
    CALL check_stop_init()
    CALL check_initial_status()
    CALL freqbins(do_imag, wsigmamin, wsigmamax, nwsigma, wcoulmax, nwcoul, &
                  wsig_wind_min, wsig_wind_max, nwsigwin, calc%freq)
    CALL sigma_grid(calc%freq, ecutsex, ecutsco, calc%grid)
    CALL opengwfil(calc%grid)
    CALL stop_clock(time_setup)
    !
  END SUBROUTINE setup_calculation

END MODULE setup
