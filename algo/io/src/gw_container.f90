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

  PUBLIC
  PRIVATE copy_coul, copy_corr, copy_exch

CONTAINS

  SUBROUTINE open_container(filename, data)
    !
    USE container_interface, ONLY: configuration
    USE gw_data,             ONLY: gw_data_container
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(gw_data_container), INTENT(OUT) :: data
    TYPE(configuration) config
    INTEGER ierr
    ! 
    config%filename = TRIM(filename)
    CALL data%open(config, ierr)
    CALL errore(__FILE__, "Error opening container for GW data", ierr)
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

  SUBROUTINE write_k_point(data)
    !
    USE gw_data, ONLY: gw_data_container, var_k_point
    USE disp,    ONLY: xk_kpoints, num_k_pts
    TYPE(gw_data_container), INTENT(INOUT) :: data
    INTEGER ierr
    !
    data%k_point = xk_kpoints(:,:num_k_pts)
    CALL data%write_variable(var_k_point, ierr)
    CALL errore(__FILE__, "Error writing k points", ierr)
    !
  END SUBROUTINE write_k_point

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

  SUBROUTINE backup(data)
    !
    USE container_interface, ONLY: configuration
    USE control_gw,          ONLY: output, do_coulomb, do_sigma_exx, do_sigma_c
    USE gw_data,             ONLY: gw_data_container
    USE io_global,           ONLY: meta_ionode, stdout
    USE mp,                  ONLY: mp_barrier
    USE mp_world,            ONLY: world_comm
    TYPE(gw_data_container), INTENT(INOUT) :: data
    TYPE(gw_data_container) backup_data
    CHARACTER(LEN=:), ALLOCATABLE :: filename, backup_file
    !
    WRITE(stdout, '(5x,a)') 'Inconsistent old datafile found, creating backup'
    filename = TRIM(output%file_data)
    backup_file = filename // '.bak'
    CALL close_container(data)
    CALL mp_barrier(world_comm)
    IF (meta_ionode) CALL RENAME(filename, backup_file)
    CALL mp_barrier(world_comm)
    !
    CALL open_container(filename, data)
    CALL open_container(backup_file, backup_data)
    IF (.NOT.do_coulomb) CALL copy_coul(data, backup_data)
    IF (.NOT.do_sigma_exx) CALL copy_exch(data, backup_data)
    IF (.NOT.do_sigma_c) CALL copy_corr(data, backup_data)
    !
  END SUBROUTINE backup

  SUBROUTINE copy_coul(data, backup_data)
    !
    USE container_interface, ONLY: element_type
    USE gw_data,             ONLY: gw_data_container, var_coul
    USE mp_world,            ONLY: mpime, nproc
    TYPE(gw_data_container), INTENT(INOUT) :: data, backup_data
    TYPE(element_type) element
    INTEGER, ALLOCATABLE :: dims(:)
    INTEGER ierr, indx
    !
    CALL backup_data%read_dimension(var_coul, dims, ierr)
    CALL errore(__FILE__, "Error reading Coulomb dimension", ierr)
    IF (.NOT.ALLOCATED(dims)) CALL errore(__FILE__, "Coulomb matrix not found", 1) 
    CALL data%write_dimension(var_coul, dims, ierr)
    DO indx = mpime + 1, dims(SIZE(dims)), nproc
      element%variable = var_coul
      element%access_index = indx
      CALL backup_data%read_element(element, ierr)
      CALL errore(__FILE__, "Error reading Coulomb data from backup file", ierr)
      data%coul = backup_data%coul
      CALL data%write_element(element, ierr)
    END DO
    !
  END SUBROUTINE copy_coul

  SUBROUTINE copy_corr(data, backup_data)
    !
    USE container_interface, ONLY: element_type
    USE gw_data,             ONLY: gw_data_container, var_corr
    USE mp_world,            ONLY: mpime, nproc
    TYPE(gw_data_container), INTENT(INOUT) :: data, backup_data
    TYPE(element_type) element
    INTEGER, ALLOCATABLE :: dims(:)
    INTEGER ierr, indx
    !
    CALL backup_data%read_dimension(var_corr, dims, ierr)
    CALL errore(__FILE__, "Error reading Correlation dimension", ierr)
    IF (.NOT.ALLOCATED(dims)) CALL errore(__FILE__, "Correlation matrix not found", 1)
    CALL data%write_dimension(var_corr, dims, ierr)
    DO indx = mpime + 1, dims(SIZE(dims)), nproc
      element%variable = var_corr
      element%access_index = indx
      CALL backup_data%read_element(element, ierr)
      CALL errore(__FILE__, "Error reading Correlation data from backup file", ierr)
      data%corr = backup_data%corr
      CALL data%write_element(element, ierr)
    END DO
    !
  END SUBROUTINE copy_corr

  SUBROUTINE copy_exch(data, backup_data)
    !
    USE container_interface, ONLY: element_type
    USE gw_data,             ONLY: gw_data_container, var_exch
    USE mp_world,            ONLY: mpime, nproc
    TYPE(gw_data_container), INTENT(INOUT) :: data, backup_data
    TYPE(element_type) element
    INTEGER, ALLOCATABLE :: dims(:)
    INTEGER ierr, indx
    !
    CALL backup_data%read_dimension(var_exch, dims, ierr)
    CALL errore(__FILE__, "Error reading Correlation dimension", ierr)
    IF (.NOT.ALLOCATED(dims)) CALL errore(__FILE__, "Correlation matrix not found", 1)
    CALL data%write_dimension(var_exch, dims, ierr)
    DO indx = mpime + 1, dims(SIZE(dims)), nproc
      element%variable = var_exch
      element%access_index = indx
      CALL backup_data%read_element(element, ierr)
      CALL errore(__FILE__, "Error reading Correlation data from backup file", ierr)
      data%exch = backup_data%exch
      CALL data%write_element(element, ierr)
    END DO
    !
  END SUBROUTINE copy_exch

END MODULE gw_container
