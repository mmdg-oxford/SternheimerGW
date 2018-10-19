! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
MODULE two_array

  USE container_interface

  IMPLICIT NONE

  PRIVATE
  PUBLIC var_second, var_first, two_array_container

  ENUM, BIND(C)
    ENUMERATOR :: var_second = 1
    ENUMERATOR var_first
  END ENUM

  TYPE, EXTENDS(container_type) :: two_array_container
    REAL, ALLOCATABLE :: second(:)
    REAL(DP), ALLOCATABLE :: first(:,:)
  CONTAINS
    PROCEDURE :: internal_init => init
    PROCEDURE :: internal_read_variable => read_variable
    PROCEDURE :: internal_read_element => read_element
    PROCEDURE :: internal_write_variable => write_variable
    PROCEDURE :: internal_write_element => write_element
    PROCEDURE :: internal_update_offset => update_offset
  END TYPE two_array_container

CONTAINS

  SUBROUTINE init(container)
    !
    CLASS(two_array_container), INTENT(OUT) :: container
    CALL allocate_copy_from_to([1, 2], container%num_dim)
    !
  END SUBROUTINE init

  SUBROUTINE write_variable(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(two_array_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    !
    SELECT CASE (config%variable)
    CASE (var_second)
      IF (ALLOCATED(config%dimension)) THEN
        CALL container%check_dimension(config%dimension, SHAPE(container%second), ierr)
      ELSE
        CALL MPI_FILE_WRITE(container%filehandle, SHAPE(container%second), 1, &
          MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      END IF
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%second, SIZE(container%second), &
        MPI_REAL, MPI_STATUS_IGNORE, ierr)
    CASE (var_first)
      IF (ALLOCATED(config%dimension)) THEN
        CALL container%check_dimension(config%dimension, SHAPE(container%first), ierr)
      ELSE
        CALL MPI_FILE_WRITE(container%filehandle, SHAPE(container%first), 2, &
          MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      END IF
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%first, SIZE(container%first), &
        MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE write_variable

  SUBROUTINE write_element(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(two_array_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    CALL allocate_copy_from_to(config%dimension, dims)
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable)
    CASE (var_second)
      CALL container%increase_offset(dims, MPI_REAL, ierr)
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%second, SIZE(container%second), &
        MPI_REAL, MPI_STATUS_IGNORE, ierr)
    CASE (var_first)
      CALL container%increase_offset(dims, MPI_DOUBLE_PRECISION, ierr)
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%first, SIZE(container%first), &
        MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE write_element

  SUBROUTINE read_variable(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(two_array_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    CALL allocate_copy_from_to(config%dimension, dims)
    SELECT CASE (config%variable)
    CASE (var_second)
      IF (ALLOCATED(container%second)) DEALLOCATE(container%second)
      ALLOCATE(container%second(dims(1)))
      CALL mpi_func(container%filehandle, container%second, SIZE(container%second), &
        MPI_REAL, MPI_STATUS_IGNORE, ierr)
    CASE (var_first)
      IF (ALLOCATED(container%first)) DEALLOCATE(container%first)
      ALLOCATE(container%first(dims(1), dims(2)))
      CALL mpi_func(container%filehandle, container%first, SIZE(container%first), &
        MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE read_variable

  SUBROUTINE read_element(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(two_array_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    CALL allocate_copy_from_to(config%dimension, dims)
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable)
    CASE (var_second)
      CALL container%increase_offset(dims, MPI_REAL, ierr)
      IF (ierr /= no_error) RETURN
      IF (ALLOCATED(container%second)) DEALLOCATE(container%second)
      ALLOCATE(container%second(1))
      CALL mpi_func(container%filehandle, container%second, SIZE(container%second), &
        MPI_REAL, MPI_STATUS_IGNORE, ierr)
    CASE (var_first)
      CALL container%increase_offset(dims, MPI_DOUBLE_PRECISION, ierr)
      IF (ierr /= no_error) RETURN
      IF (ALLOCATED(container%first)) DEALLOCATE(container%first)
      ALLOCATE(container%first(dims(1), 1))
      CALL mpi_func(container%filehandle, container%first, SIZE(container%first), &
        MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE read_element

  SUBROUTINE update_offset(container, variable, dimension, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(two_array_container), INTENT(INOUT) :: container
    INTEGER, INTENT(IN) :: variable, dimension(:)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER(KIND=MPI_OFFSET_KIND) offset
    !
    SELECT CASE (variable)
    CASE (var_second)
      CALL container%increase_offset(dimension, MPI_REAL, ierr)
    CASE (var_first)
      CALL container%increase_offset(dimension, MPI_DOUBLE_PRECISION, ierr)
    END SELECT
    IF (ierr /= no_error) RETURN
    CALL container%update_offset(ierr)
    !
  END SUBROUTINE update_offset

END MODULE two_array
