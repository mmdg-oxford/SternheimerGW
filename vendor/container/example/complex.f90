! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
MODULE complex_array

  USE container_interface

  IMPLICIT NONE

  PRIVATE
  PUBLIC var_double, var_single, complex_array_container

  ENUM, BIND(C)
    ENUMERATOR :: var_double = 1
    ENUMERATOR var_single
  END ENUM

  TYPE, EXTENDS(container_type) :: complex_array_container
    COMPLEX(DP), ALLOCATABLE :: double(:,:,:)
    COMPLEX, ALLOCATABLE :: single(:,:)
  CONTAINS
    PROCEDURE :: internal_init => init
    PROCEDURE :: internal_read_variable => read_variable
    PROCEDURE :: internal_read_element => read_element
    PROCEDURE :: internal_write_variable => write_variable
    PROCEDURE :: internal_write_element => write_element
    PROCEDURE :: internal_update_offset => update_offset
  END TYPE complex_array_container

CONTAINS

  SUBROUTINE init(container)
    !
    CLASS(complex_array_container), INTENT(OUT) :: container
    container%num_dim = [3, 2]
    !
  END SUBROUTINE init

  SUBROUTINE write_variable(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(complex_array_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    !
    SELECT CASE (config%variable)
    CASE (var_double)
      IF (ALLOCATED(config%dimension)) THEN
        CALL container%check_dimension(config%dimension, SHAPE(container%double), ierr)
      ELSE
        CALL MPI_FILE_WRITE(container%filehandle, SHAPE(container%double), 3, &
          MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      END IF
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%double, SIZE(container%double), &
        MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
    CASE (var_single)
      IF (ALLOCATED(config%dimension)) THEN
        CALL container%check_dimension(config%dimension, SHAPE(container%single), ierr)
      ELSE
        CALL MPI_FILE_WRITE(container%filehandle, SHAPE(container%single), 2, &
          MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      END IF
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%single, SIZE(container%single), &
        MPI_COMPLEX, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE write_variable

  SUBROUTINE write_element(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(complex_array_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable)
    CASE (var_double)
      CALL container%increase_offset(dims, MPI_DOUBLE_COMPLEX, ierr)
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%double, SIZE(container%double), &
        MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
    CASE (var_single)
      CALL container%increase_offset(dims, MPI_COMPLEX, ierr)
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%single, SIZE(container%single), &
        MPI_COMPLEX, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE write_element

  SUBROUTINE read_variable(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(complex_array_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    SELECT CASE (config%variable)
    CASE (var_double)
      IF (ALLOCATED(container%double)) DEALLOCATE(container%double)
      ALLOCATE(container%double(dims(1), dims(2), dims(3)))
      CALL mpi_func(container%filehandle, container%double, SIZE(container%double), &
        MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
    CASE (var_single)
      IF (ALLOCATED(container%single)) DEALLOCATE(container%single)
      ALLOCATE(container%single(dims(1), dims(2)))
      CALL mpi_func(container%filehandle, container%single, SIZE(container%single), &
        MPI_COMPLEX, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE read_variable

  SUBROUTINE read_element(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(complex_array_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable)
    CASE (var_double)
      CALL container%increase_offset(dims, MPI_DOUBLE_COMPLEX, ierr)
      IF (ierr /= no_error) RETURN
      IF (ALLOCATED(container%double)) DEALLOCATE(container%double)
      ALLOCATE(container%double(dims(1), dims(2), 1))
      CALL mpi_func(container%filehandle, container%double, SIZE(container%double), &
        MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
    CASE (var_single)
      CALL container%increase_offset(dims, MPI_COMPLEX, ierr)
      IF (ierr /= no_error) RETURN
      IF (ALLOCATED(container%single)) DEALLOCATE(container%single)
      ALLOCATE(container%single(dims(1), 1))
      CALL mpi_func(container%filehandle, container%single, SIZE(container%single), &
        MPI_COMPLEX, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE read_element

  SUBROUTINE update_offset(container, variable, dimension, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(complex_array_container), INTENT(INOUT) :: container
    INTEGER, INTENT(IN) :: variable, dimension(:)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER(KIND=MPI_OFFSET_KIND) offset
    !
    SELECT CASE (variable)
    CASE (var_double)
      CALL container%increase_offset(dimension, MPI_DOUBLE_COMPLEX, ierr)
    CASE (var_single)
      CALL container%increase_offset(dimension, MPI_COMPLEX, ierr)
    END SELECT
    IF (ierr /= no_error) RETURN
    CALL container%update_offset(ierr)
    !
  END SUBROUTINE update_offset

END MODULE complex_array
