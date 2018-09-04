! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
MODULE simple

  USE container_interface

  IMPLICIT NONE

  PRIVATE
  PUBLIC var_array, simple_container

  ENUM, BIND(C)
    ENUMERATOR :: var_array = 2
  END ENUM

  TYPE, EXTENDS(container_type) :: simple_container
    INTEGER, ALLOCATABLE :: array(:,:)
  CONTAINS
    PROCEDURE :: internal_init => init
    PROCEDURE :: internal_read_variable => read_variable
    PROCEDURE :: internal_read_element => read_element
    PROCEDURE :: internal_write_variable => write_variable
    PROCEDURE :: internal_write_element => write_element
    PROCEDURE :: internal_update_offset => update_offset
  END TYPE simple_container

CONTAINS

  SUBROUTINE init(container)
    !
    CLASS(simple_container), INTENT(OUT) :: container
    container%num_dim = [1, 2]
    !
  END SUBROUTINE init

  SUBROUTINE write_variable(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(simple_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    !
    SELECT CASE (config%variable)
    CASE (var_array)
      IF (ALLOCATED(config%dimension)) THEN
        CALL container%check_dimension(config%dimension, SHAPE(container%array), ierr)
      ELSE
        CALL MPI_FILE_WRITE(container%filehandle, SHAPE(container%array), 2, &
          MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      END IF
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%array, SIZE(container%array), &
        MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE write_variable

  SUBROUTINE write_element(container, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(simple_container), INTENT(INOUT) :: container
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable)
    CASE (var_array)
      CALL container%increase_offset(dims, MPI_INTEGER, ierr)
      IF (ierr /= no_error) RETURN
      CALL MPI_FILE_WRITE(container%filehandle, container%array, SIZE(container%array), &
        MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE write_element

  SUBROUTINE read_variable(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(simple_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    SELECT CASE (config%variable)
    CASE (var_array)
      IF (ALLOCATED(container%array)) DEALLOCATE(container%array)
      ALLOCATE(container%array(dims(1), dims(2)))
      CALL mpi_func(container%filehandle, container%array, SIZE(container%array), &
        MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE read_variable

  SUBROUTINE read_element(container, mpi_func, config, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(simple_container), INTENT(INOUT) :: container
    EXTERNAL mpi_func
    TYPE(internal_config), INTENT(IN) :: config
    INTEGER, INTENT(OUT) :: ierr
    INTEGER, ALLOCATABLE :: dims(:)
    !
    dims = config%dimension
    dims(SIZE(dims)) = config%access_index - 1
    SELECT CASE (config%variable)
    CASE (var_array)
      CALL container%increase_offset(dims, MPI_INTEGER, ierr)
      IF (ierr /= no_error) RETURN
      IF (ALLOCATED(container%array)) DEALLOCATE(container%array)
      ALLOCATE(container%array(dims(1), 1))
      CALL mpi_func(container%filehandle, container%array, SIZE(container%array), &
        MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
    END SELECT
    !
  END SUBROUTINE read_element

  SUBROUTINE update_offset(container, variable, dimension, ierr)
    !
    INCLUDE 'mpif.h'
    CLASS(simple_container), INTENT(INOUT) :: container
    INTEGER, INTENT(IN) :: variable, dimension(:)
    INTEGER, INTENT(OUT) :: ierr
    INTEGER(KIND=MPI_OFFSET_KIND) offset
    !
    SELECT CASE (variable)
    CASE (var_array)
      CALL container%increase_offset(dimension, MPI_INTEGER, ierr)
    END SELECT
    IF (ierr /= no_error) RETURN
    CALL container%update_offset(ierr)
    !
  END SUBROUTINE update_offset

END MODULE simple
