!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! 
! Copyright (C) 2010 - 2016 
! Henry Lambert, Martin Schlipf, and Feliciano Giustino
!
! Sternheimer-GW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Sternheimer-GW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Sternheimer-GW. If not, see
! http://www.gnu.org/licenses/gpl.html .
!
!------------------------------------------------------------------------------ 
!> This module provides the subroutines to parallelize a task over a communicator.
MODULE parallel_module

  IMPLICIT NONE

  !> wrapper for the MPI_ALLGATHERV routine
  !!
  !! The subroutines are called with an array that is distributed according
  !! across the communicator. We must know which part of the array is done by
  !! which process.
  !!
  !! It is recommended that you use the parallel_task subroutine to generate
  !! the distribution. Then you can just pass the num_task array generated by
  !! that subroutine.
  INTERFACE mp_allgatherv
    MODULE PROCEDURE mp_allgatherv_rv, mp_allgatherv_rm, &
                     mp_allgatherv_cv, mp_allgatherv_cm
  END INTERFACE mp_allgatherv

  PRIVATE

  PUBLIC parallel_task, mp_allgatherv

CONTAINS

  !> Parallelize a given number of tasks over given communicator.
  !!
  !! Distribution of tasks over processes in communicator according to the
  !! following strategy
  !!
  SUBROUTINE parallel_task(comm, num_task_total, first_task, last_task, num_task)

    USE mp_global, ONLY: mp_rank, mp_size

    !> the communicator over which the tasks are distributed
    INTEGER, INTENT(IN)  :: comm

    !> total number of tasks to assign
    INTEGER, INTENT(IN)  :: num_task_total

    !> first task assigned to this process
    INTEGER, INTENT(OUT) :: first_task

    !> last task assigned to this process
    INTEGER, INTENT(OUT) :: last_task

    !> the number of tasks assigned to each process
    !! @note you can use this for the mp_gather wrapper in this module
    INTEGER, INTENT(OUT), ALLOCATABLE :: num_task(:)

    !> number of processes in communicator
    INTEGER num_proc

    !> rank of this process in communicator
    INTEGER my_rank

    !> minimal number of task per process
    INTEGER num_task_min

    !> remainder after assigning an equal amount of processes
    INTEGER num_remain

    !> at this point we add one additional task
    INTEGER last_proc

    ! determine rank and size
    ! note: add 1 to rank, because Fortran indices start counting at 1
    my_rank = mp_rank(comm) + 1
    num_proc = mp_size(comm)

    ! allocate array for assigned tasks
    ALLOCATE(num_task(num_proc))

    !! 1. Distribute num_task / size(comm) on every process.
    num_task_min = num_task_total / num_proc

    !! 2. Determine the remainder of num_remain of unassigned tasks.
    num_remain = MOD(num_task_total, num_proc)

    !! 3. Assign the last processes in the list an extra task.
    last_proc = num_proc - num_remain
    num_task(:last_proc) = num_task_min
    num_task(last_proc+1:) = num_task_min + 1

    !! 4. Determine first and last task for current process.
    last_task = SUM(num_task(:my_rank))
    first_task = last_task - num_task(my_rank) + 1

  END SUBROUTINE parallel_task

  !> Gather a vector of reals.
  SUBROUTINE mp_allgatherv_rv(comm, num_task, array)

    USE kinds, ONLY: dp
    USE mp,    ONLY: mp_size
    USE parallel_include

    !> The communicator across which the tasks are distributed.
    INTEGER,  INTENT(IN)    :: comm

    !> The distribution of the tasks across the processes.
    INTEGER,  INTENT(IN)    :: num_task(:)

    !> The vector that is distributed across the processes.
    REAL(dp), INTENT(INOUT) :: array(:)

    !> displacement array - sum the number of tasks x dimension of the array
    INTEGER,  ALLOCATABLE   :: displ(:)

    ! number of processes
    INTEGER num_proc

    ! loop variable for processes
    INTEGER iproc
 
    ! MPI error code
    INTEGER ierr

    ! we only need to gather the vector if it is distributed
#ifdef __MPI
  
    ! determine number of processes
    num_proc = SIZE(num_task)

    ! calculate displacement
    ALLOCATE(displ(num_proc))

    ! note: C array storage - starting at 0
    displ(1) = 0
    DO iproc = 1, num_proc - 1
      displ(iproc + 1) = displ(iproc) + num_task(iproc)
    END DO ! iproc

    ! sanity test of the input
    IF (mp_size(comm) /= num_proc) &
      CALL errore(__FILE__, "error in allgatherv - communicator size and distribution don't match", 1)
    IF (displ(num_proc) + num_task(num_proc) /= SIZE(array)) &
      CALL errore(__FILE__, "array size does not match with assigned tasks", 1)

    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, array, num_task, displ, MPI_DOUBLE_PRECISION, comm, ierr)
    CALL errore(__FILE__, "error in mpi_allgatherv call", ierr)

    DEALLOCATE(displ)

#endif

  END SUBROUTINE mp_allgatherv_rv

  !> Gather a vector of complex.
  SUBROUTINE mp_allgatherv_cv(comm, num_task, array)

    USE kinds, ONLY: dp
    USE mp,    ONLY: mp_size
    USE parallel_include

    !> The communicator across which the tasks are distributed.
    INTEGER,     INTENT(IN)    :: comm

    !> The distribution of the tasks across the processes.
    INTEGER,     INTENT(IN)    :: num_task(:)

    !> The vector that is distributed across the processes.
    COMPLEX(dp), INTENT(INOUT) :: array(:)

    !> displacement array - sum the number of tasks x dimension of the array
    INTEGER,     ALLOCATABLE   :: displ(:)

    ! number of processes
    INTEGER num_proc

    ! loop variable for processes
    INTEGER iproc
 
    ! MPI error code
    INTEGER ierr

    ! we only need to gather the vector if it is distributed
#ifdef __MPI
  
    ! determine number of processes
    num_proc = SIZE(num_task)

    ! calculate displacement
    ALLOCATE(displ(num_proc))

    ! note: C array storage - starting at 0
    displ(1) = 0
    DO iproc = 1, num_proc - 1
      displ(iproc + 1) = displ(iproc) + num_task(iproc)
    END DO ! iproc

    ! sanity test of the input
    IF (mp_size(comm) /= num_proc) &
      CALL errore(__FILE__, "error in allgatherv - communicator size and distribution don't match", 1)
    IF (displ(num_proc) + num_task(num_proc) /= SIZE(array)) &
      CALL errore(__FILE__, "array size does not match with assigned tasks", 1)

    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, array, num_task, displ, MPI_DOUBLE_COMPLEX, comm, ierr)
    CALL errore(__FILE__, "error in mpi_allgatherv call", ierr)

    DEALLOCATE(displ)

#endif

  END SUBROUTINE mp_allgatherv_cv

  !> Gather a matrix of reals.
  SUBROUTINE mp_allgatherv_rm(comm, num_task, array)

    USE kinds, ONLY: dp
    USE mp,    ONLY: mp_size
    USE parallel_include

    !> The communicator across which the tasks are distributed.
    INTEGER,  INTENT(IN)    :: comm

    !> The distribution of the tasks across the processes.
    INTEGER,  INTENT(IN)    :: num_task(:)

    !> The vector that is distributed across the processes.
    REAL(dp), INTENT(INOUT) :: array(:,:)

    !> actual size of sent data - number of task x dimension of the array
    INTEGER,  ALLOCATABLE   :: receive(:)

    !> displacement array - sum the number of tasks x dimension of the array
    INTEGER,  ALLOCATABLE   :: displ(:)

    ! number of processes
    INTEGER num_proc

    ! loop variable for processes
    INTEGER iproc
 
    ! MPI error code
    INTEGER ierr

    ! we only need to gather the vector if it is distributed
#ifdef __MPI
  
    ! determine number of processes
    num_proc = SIZE(num_task)

    ! calculate number of sent information
    ALLOCATE(receive(num_proc))
    receive = num_task * SIZE(array, 1)

    ! calculate displacement
    ALLOCATE(displ(num_proc))

    ! note: C array storage - starting at 0
    displ(1) = 0
    DO iproc = 1, num_proc - 1
      displ(iproc + 1) = displ(iproc) + receive(iproc)
    END DO ! iproc

    ! sanity test of the input
    IF (mp_size(comm) /= num_proc) &
      CALL errore(__FILE__, "error in allgatherv - communicator size and distribution don't match", 1)
    IF (displ(num_proc) + receive(num_proc) /= SIZE(array)) &
      CALL errore(__FILE__, "array size does not match with assigned tasks", 1)

    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, array, receive, displ, MPI_DOUBLE_PRECISION, comm, ierr)
    CALL errore(__FILE__, "error in mpi_allgatherv call", ierr)

    DEALLOCATE(displ)
    DEALLOCATE(receive)

#endif

  END SUBROUTINE mp_allgatherv_rm

  !> Gather a matrix of complex.
  SUBROUTINE mp_allgatherv_cm(comm, num_task, array)

    USE kinds, ONLY: dp
    USE mp,    ONLY: mp_size
    USE parallel_include

    !> The communicator across which the tasks are distributed.
    INTEGER,     INTENT(IN)    :: comm

    !> The distribution of the tasks across the processes.
    INTEGER,     INTENT(IN)    :: num_task(:)

    !> The vector that is distributed across the processes.
    COMPLEX(dp), INTENT(INOUT) :: array(:,:)

    !> actual size of sent data - number of task x dimension of the array
    INTEGER,     ALLOCATABLE   :: receive(:)

    !> displacement array - sum the number of tasks x dimension of the array
    INTEGER,     ALLOCATABLE   :: displ(:)

    ! number of processes
    INTEGER num_proc

    ! loop variable for processes
    INTEGER iproc
 
    ! MPI error code
    INTEGER ierr

    ! we only need to gather the vector if it is distributed
#ifdef __MPI
  
    ! determine number of processes
    num_proc = SIZE(num_task)

    ! calculate number of sent information
    ALLOCATE(receive(num_proc))
    receive = num_task * SIZE(array, 1)

    ! calculate displacement
    ALLOCATE(displ(num_proc))

    ! note: C array storage - starting at 0
    displ(1) = 0
    DO iproc = 1, num_proc - 1
      displ(iproc + 1) = displ(iproc) + receive(iproc)
    END DO ! iproc

    ! sanity test of the input
    IF (mp_size(comm) /= num_proc) &
      CALL errore(__FILE__, "error in allgatherv - communicator size and distribution don't match", 1)
    IF (displ(num_proc) + receive(num_proc) /= SIZE(array)) &
      CALL errore(__FILE__, "array size does not match with assigned tasks", 1)

    CALL MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, array, receive, displ, MPI_DOUBLE_COMPLEX, comm, ierr)
    CALL errore(__FILE__, "error in mpi_allgatherv call", ierr)

    DEALLOCATE(displ)
    DEALLOCATE(receive)

#endif

  END SUBROUTINE mp_allgatherv_cm

END MODULE parallel_module
