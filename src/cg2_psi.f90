  !-----------------------------------------------------------------------
  ! Copyright (C) 2010-2015 Henry Lambert, Feliciano Giustino
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !-----------------------------------------------------------------------
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------
subroutine cg2_psi (lda, n, m, psi, h_diag)
  !-----------------------------------------------------------------
  !
  !    This routine gives a preconditioning to the linear system solver.
  !    The preconditioning is diagonal in reciprocal space
  !
  !
  USE kinds, only : DP
  USE noncollin_module, only : noncolin, npol
  implicit none

  integer :: lda, n, m
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors

  complex(DP) :: psi (lda*npol, m)
  ! inp/out: the vector to be preconditioned

  real(DP) :: h_diag (lda*npol, m)
  !COMPLEX(DP) :: h_diag (lda*npol, m)
  ! input: the preconditioning vector

  integer :: k, i
  ! counter on bands
  ! counter on the elements of the vector
  !
  !transform solution vector x = E^{-T}x':
  do k = 1, m
     do i = 1, n
        psi (i, k) = psi (i, k) * sqrt(h_diag (i, k))
     enddo
  enddo

  IF (noncolin) THEN
     do k = 1, m
        do i = 1, n
           psi (i+lda, k) = psi (i+lda, k) * h_diag (i+lda, k)
        enddo
     enddo
  END IF
  return
end subroutine cg2_psi
