!------------------------------------------------------------------------------
!
! This file is part of the SternheimerGW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2018 Quantum ESPRESSO group,
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
SUBROUTINE clean_pw_gw()
  !-----------------------------------------------------------------------
  !
  ! This routine deallocates all the variables of pwscf and of the
  ! GW code and resets the same variables as after reading input in
  ! gwq_readin, so that it is possible to start a calculation at
  ! a new q.
  !
  USE buffers,      ONLY: close_buffer
  USE lr_symm_base, ONLY: nsymq
  USE units_gw,     ONLY: iuwfc, iubar
  !
  IMPLICIT NONE
  !
  CALL clean_pw(.FALSE.)
  CALL deallocate_gwq()
  nsymq = 0
  CALL close_buffer(iuwfc, 'delete')
  CALL close_buffer(iubar, 'delete')
  !
END SUBROUTINE clean_pw_gw
