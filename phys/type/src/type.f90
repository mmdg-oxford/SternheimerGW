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
MODULE gw_type

  USE debug_module,         ONLY: debug_type
  USE freqbins_module,      ONLY: freqbins_type
  USE gw_data,              ONLY: gw_data_container
  USE select_solver_module, ONLY: select_solver_type
  USE setup_nscf_module,    ONLY: sigma_config_type
  USE sigma_grid_module,    ONLY: sigma_grid_type
  USE truncation_module,    ONLY: vcut_type

  PRIVATE
  PUBLIC calculation

  TYPE calculation
    !> the debug configuration of the calculation
    TYPE(debug_type) debug
    !> stores the frequencies uses for the calculation
    TYPE(freqbins_type) freq
    !> interface the larger arrays with a parallel file IO
    TYPE(gw_data_container) data
    !> stores the configuration of the linear solver for the screened Coulomb interaction
    TYPE(select_solver_type) config_coul
    !> stores the configuration of the linear solver for the Green's function
    TYPE(select_solver_type) config_green
    !> stores the configuration of the self-energy calculation
    TYPE(sigma_config_type), ALLOCATABLE :: config(:)
    !> stores the FFT grids used in the calculation
    TYPE(sigma_grid_type) grid
    !> stores the truncated Coulomb potential
    TYPE(vcut_type) vcut
  END TYPE calculation

END MODULE gw_type
