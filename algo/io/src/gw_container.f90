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
MODULE gw_container

  IMPLICIT NONE

CONTAINS

  SUBROUTINE open_container(data)
    !
    USE container_interface, ONLY: configuration
    USE control_gw,          ONLY: output
    USE disp,                ONLY: xk_kpoints, num_k_pts
    USE gw_data,             ONLY: gw_data_container, var_k_point
    !
    TYPE(gw_data_container), INTENT(OUT) :: data
    TYPE(configuration) config
    INTEGER ierr
    ! 
    config%filename = output%file_data
    CALL data%open(config, ierr)
    CALL errore(__FILE__, "Error opening container for GW data", ierr)
    ALLOCATE(data%k_point(3, num_k_pts))
    data%k_point = xk_kpoints(:,:num_k_pts)
    CALL data%write_variable(var_k_point, ierr)
    CALL errore(__FILE__, "Error writing k points", ierr)
    ! 
  END SUBROUTINE open_container

  SUBROUTINE close_container(data)
    !
    USE gw_data, ONLY: gw_data_container
    TYPE(gw_data_container), INTENT(INOUT) :: data
    INTEGER ierr
    !
    CALL data%close(ierr)
    CALL errore(__FILE__, "Error closing container of GW data", ierr)
    !
  END SUBROUTINE close_container

  SUBROUTINE write_exch_dim(data, exch_dim)
    !
    USE gw_data, ONLY: gw_data_container, var_exch
    TYPE(gw_data_container), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: exch_dim(:)
    INTEGER ierr
    !
    CALL data%write_dimension(var_exch, exch_dim, ierr)
    CALL errore(__FILE__, "Error writing dimension of exchange", ierr)
    !
  END SUBROUTINE write_exch_dim

  SUBROUTINE write_corr_dim(data, corr_dim)
    !
    USE gw_data, ONLY: gw_data_container, var_corr
    TYPE(gw_data_container), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: corr_dim(:)
    INTEGER ierr
    !
    CALL data%write_dimension(var_corr, corr_dim, ierr)
    CALL errore(__FILE__, "Error writing dimension of correlation", ierr)
    !
  END SUBROUTINE write_corr_dim

END MODULE gw_container
