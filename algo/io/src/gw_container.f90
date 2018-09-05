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

  TYPE gw_dimension
    INTEGER, ALLOCATABLE :: exch(:), corr(:)
  END TYPE gw_dimension

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
    config%filename = TRIM(output%file_data)
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

  LOGICAL FUNCTION consistent_dimension(data, dims)
    !
    USE gw_data, ONLY: gw_data_container, var_exch, var_corr
    TYPE(gw_data_container), INTENT(INOUT) :: data
    TYPE(gw_dimension), INTENT(IN) :: dims
    !
    consistent_dimension = .FALSE.
    IF (ALLOCATED(dims%exch)) THEN
      IF (.NOT.dimension_correct(data, var_exch, dims%exch)) RETURN
    END IF
    IF (ALLOCATED(dims%corr)) THEN
      IF (.NOT.dimension_correct(data, var_corr, dims%corr)) RETURN
    END IF
    consistent_dimension = .TRUE.
    !
  END FUNCTION consistent_dimension

  LOGICAL FUNCTION dimension_correct(data, var, ref_dim)
    !
    USE container_interface, ONLY: no_error
    USE gw_data,             ONLY: gw_data_container
    TYPE(gw_data_container), INTENT(INOUT) :: data
    INTEGER, INTENT(IN) :: var, ref_dim(:)
    INTEGER ierr
    INTEGER, ALLOCATABLE :: test_dim(:)
    !
    dimension_correct = .TRUE.
    CALL data%read_dimension(var, test_dim, ierr)
    CALL errore(__FILE__, "Error reading dimension", ierr)
    IF (.NOT.ALLOCATED(test_dim)) RETURN
    CALL data%check_dimension(ref_dim, test_dim, ierr)
    dimension_correct = (ierr == no_error)
    !
  END FUNCTION dimension_correct 

  SUBROUTINE write_dimension(data, dims)
    !
    USE gw_data, ONLY: gw_data_container, var_exch, var_corr
    TYPE(gw_data_container), INTENT(INOUT) :: data
    TYPE(gw_dimension), INTENT(IN) :: dims
    INTEGER ierr
    !
    IF (ALLOCATED(dims%exch)) THEN
      CALL data%write_dimension(var_exch, dims%exch, ierr)
      CALL errore(__FILE__, "Error writing dimension of exchange", ierr)
    END IF
    IF (ALLOCATED(dims%corr)) THEN
      CALL data%write_dimension(var_corr, dims%corr, ierr)
      CALL errore(__FILE__, "Error writing dimension of correlation", ierr)
    END IF
    !
  END SUBROUTINE write_dimension

END MODULE gw_container
