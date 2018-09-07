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
!> Provides routines that generate output for QE's plotband program.
MODULE pp_output_mod

  USE kinds,       ONLY : dp
  USE gw_type_mod, ONLY : pp_output_type

  IMPLICIT NONE

CONTAINS

  !> Open all files that are use for PP printing
  !!
  !! \param nks number of k-points the data will contain
  !! \param nbnd number of bands for GW quasiparticle energies
  !! \param nw_re number of frequency points on real axis
  !! \param nw_im number of frequency point on imaginary axis
  !! \param output type that contains all files for PP
  !!
  SUBROUTINE pp_output_open_all(nks, nbnd, nw_re, nw_im, output)

    USE gw_type_mod, ONLY : name_length, output_type

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(IN) :: nw_re
    INTEGER, INTENT(IN) :: nw_im
    TYPE(output_type), INTENT(INOUT) :: output

    CHARACTER(LEN=name_length) prefix

    ! bands for band structures
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_dft)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_gw)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_vxc)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_exchange)
    CALL pp_output_open(nks, nbnd, output%directory, output%pp_renorm)

    ! prefix contains directory and file prefix
    prefix = TRIM(output%directory)//"/"//TRIM(output%prefix)//"."

  END SUBROUTINE pp_output_open_all

  !> Open a file to print the data for QE's plotband program
  !!
  !! If the filename is not set, this routine will just clear the to_file
  !! flag, so that later parts of the code can test whether data is meant
  !! to be written. If the filename is present, the file is opened and the
  !! unit is stored in the type. Then the header for the data is written
  !! into the file.
  !!
  !! \param nks number of k-points the data will contain
  !! \param nbnd number of bands the data will have
  !! \param directory directory in which files are created
  !! \param output type that contains the filename on input and the unit
  !! and some metadata after the return of the function
  !!
  SUBROUTINE pp_output_open(nks, nbnd, directory, output)

    USE io_files, ONLY : seqopn

    INTEGER, INTENT(IN) :: nks
    INTEGER, INTENT(IN) :: nbnd
    CHARACTER(*), INTENT(IN) :: directory
    TYPE(pp_output_type), INTENT(INOUT) :: output

    INTEGER, EXTERNAL :: find_free_unit
    LOGICAL exst

    NAMELIST /plot/ nks, nbnd

    ! if no filename is present, clear to_file flag and exit
    output%to_file = (output%filename /= '')
    IF (.NOT.output%to_file) RETURN

    ! set metadata
    output%num_band = nbnd
    output%num_freq = 1
    output%num_kpoint = nks

    ! open the file
    output%iunit = find_free_unit()
    CALL seqopn(output%iunit, output%filename, "FORMATTED", exst, TRIM(directory)//"/")

    ! write namelist to file
    WRITE(output%iunit, NML=plot)

  END SUBROUTINE pp_output_open

  !> Close a file associated to a particular output type
  !!
  !! Depending on the format of the file either the xml or the
  !! regular file closing routine are used.
  !!
  !! \param output file (as pp_output_type) which is closed
  !!
  SUBROUTINE pp_output_close(output)

    USE iotk_module, ONLY: iotk_write_end, iotk_close_write

    TYPE(pp_output_type), INTENT(IN) :: output

    LOGICAL opnd

    ! don't do anything if file is not open
    IF (.NOT.output%to_file) RETURN

    INQUIRE(UNIT = output%iunit, OPENED = opnd)
    IF (.NOT.opnd) CALL errore(__FILE__, output%filename//' not opened', 1)
    CLOSE(output%iunit)

  END SUBROUTINE pp_output_close

  !> Print an 1D array for QE's plotband program
  !!
  !! The array has the dimension of the number of bands. It is assumed that
  !! this routine is called nks times, where nks is the number of k-points
  !! specified when opening the file.
  !!
  !! \param output file (as pp_output_type) to which the data is printed
  !! \param kpt vector containing the current k-point
  !! \param data array containing the data
  SUBROUTINE pp_output(output, kpt, data)

    TYPE(pp_output_type), INTENT(IN) :: output
    REAL(dp), INTENT(IN) :: kpt(3)
    REAL(dp), INTENT(IN) :: data(:)

    LOGICAL opnd

    ! don't do anything if file is not required
    IF (.NOT.output%to_file) RETURN

    !
    ! sanity test of the input
    !
    CALL errore(__FILE__, 'data array size inconsistent (band)', output%num_band - SIZE(data))
    CALL errore(__FILE__, 'data array size inconsistent (freq)', output%num_freq - 1)
    INQUIRE(UNIT = output%iunit, OPENED = opnd)
    IF (.NOT.opnd) CALL errore(__FILE__, output%filename//' not opened', 1)

    !
    ! write the data to the file
    !
    WRITE(output%iunit, '(5x,3f10.6)') kpt
    WRITE(output%iunit, '(10f15.8)') data
    ! add an empty line at the end of one data set
    WRITE(output%iunit,*)

  END SUBROUTINE pp_output

END MODULE pp_output_mod
