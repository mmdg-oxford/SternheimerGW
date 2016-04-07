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
!> this module provides functionality to reorder an array
!!
!! given an array and a map, the array is modified in-place so that the
!! new order is compatible with the map
MODULE reorder_mod

  USE kinds, ONLY : dp

  IMPLICIT NONE

  INTERFACE reorder
    MODULE PROCEDURE reorder_r1d, reorder_r2d, &
                     reorder_c1d, reorder_c2d
  END INTERFACE reorder

CONTAINS

  !> reorder a one dimensional real array
  !!
  !! \param array on entry: array with old order; on exit: array with order
  !!              according to map
  !! \param map new order for array
  !! \param max_indx skip elements larger that max_indx (optional default: no limit)
  SUBROUTINE reorder_r1d(array, map, max_indx)

    REAL(dp), INTENT(INOUT) :: array(:)
    INTEGER,  INTENT(IN)    :: map(:)
    INTEGER,  INTENT(IN), OPTIONAL :: max_indx

    INTEGER upper

    ! array and map must have same size
    CALL errore("reorder", "array shape mismatch", size(array) /= size(map))

    ! determine the upper limit
    IF (present(max_indx)) THEN
      upper = max_indx
    ELSE
      upper = size(array)
    END IF

    ! if the map points somewhere copy it
    WHERE ( map /= 0 .AND. map <= upper )
      array = array(map)
    ELSEWHERE
      array = 0
    END WHERE

  END SUBROUTINE reorder_r1d

  !> reorder a two dimensional real array (first index)
  !!
  !! \param array on entry: array with old order; on exit: array with order
  !!              according to map
  !! \param map new order for array
  !! \param max_indx skip elements larger that max_indx (optional default: no limit)
  SUBROUTINE reorder_r2d(array, map, max_indx)

    REAL(dp), INTENT(INOUT) :: array(:,:)
    INTEGER,  INTENT(IN)    :: map(:)
    INTEGER,  INTENT(IN), OPTIONAL :: max_indx

    INTEGER ii

    DO ii = 1, size(array,2)
      CALL reorder_r1d(array(:,ii), map, max_indx)
    END DO

  END SUBROUTINE reorder_r2d

  !> reorder a one dimensional complex array
  !!
  !! \param array on entry: array with old order; on exit: array with order
  !!              according to map
  !! \param map new order for array
  !! \param max_indx skip elements larger that max_indx (optional default: no limit)
  SUBROUTINE reorder_c1d(array, map, max_indx)

    COMPLEX(dp), INTENT(INOUT) :: array(:)
    INTEGER,     INTENT(IN)    :: map(:)
    INTEGER,  INTENT(IN), OPTIONAL :: max_indx

    INTEGER upper

    ! array and map must have same size
    CALL errore("reorder", "array shape mismatch", size(array) /= size(map))

    ! determine the upper limit
    IF (present(max_indx)) THEN
      upper = max_indx
    ELSE
      upper = size(array)
    END IF

    ! if the map points somewhere copy it
    WHERE ( map /= 0 .AND. map <= upper )
      array = array(map)
    ELSEWHERE
      array = 0
    END WHERE

  END SUBROUTINE reorder_c1d

  !> reorder a two dimensional complex array (first index)
  !!
  !! \param array on entry: array with old order; on exit: array with order
  !!              according to map
  !! \param map new order for array
  !! \param max_indx skip elements larger that max_indx (optional default: no limit)
  SUBROUTINE reorder_c2d(array, map, max_indx)

    COMPLEX(dp), INTENT(INOUT) :: array(:,:)
    INTEGER,     INTENT(IN)    :: map(:)
    INTEGER,  INTENT(IN), OPTIONAL :: max_indx

    INTEGER ii

    DO ii = 1, size(array,2)
      CALL reorder_c1d(array(:,ii), map, max_indx)
    END DO

  END SUBROUTINE reorder_c2d

END MODULE reorder_mod