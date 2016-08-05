!------------------------------------------------------------------------------
!
! This file is part of the Sternheimer-GW code.
! Parts of this file are taken from the Quantum ESPRESSO software
! P. Giannozzi, et al, J. Phys.: Condens. Matter, 21, 395502 (2009)
!
! Copyright (C) 2010 - 2016 Quantum ESPRESSO group,
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
subroutine allocate_gwq
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
  USE kinds,    ONLY : DP
  USE klist,    ONLY : nks, nkstot
  USE wvfct,    ONLY : nbnd, npwx
  USE gvect,    ONLY : ngm
  USE lsda_mod, ONLY : nspin
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag
  USE wavefunctions_module,  ONLY : evc
  USE spin_orb,       ONLY : lspinorb
  USE becmod,         ONLY : bec_type, becp, allocate_bec_type
  USE uspp,           ONLY : okvan, nkb
  USE paw_variables,  ONLY : okpaw
  USE uspp_param,     ONLY : nhm
  USE freq_gw,        ONLY : fiu, nfs, nfsmax
  USE ions_base,      ONLY : nat, ntyp => nsp
  USE lrus,          ONLY : becp1
  USE qpoint,        ONLY : nksq, eigqts, igkq
  USE eqv_gw,        ONLY : dpsi, evq, vlocq, dmuxc, dvpsi, &
                            dvbare, dpsim, dpsip, eprectot
  USE units_gw,      ONLY : this_pcxpsi_is_on_file, this_dvkb3_is_on_file
  USE control_gw,    ONLY : lgamma  
  USE fft_base,      ONLY : dfftp, dffts
  USE disp,          ONLY : gmap, eval_occ

  IMPLICIT NONE

  INTEGER :: ik, ipol
  !
  !   FOR LGAMMA
  if (lgamma) then
    !
    ! q=0  : evq is a pointer to evc
    !
    evq  => evc
  else
    !
    ! q!=0 : evq is and calculated at point k+q
    !
    allocate (evq ( npwx*npol , nbnd))    
    allocate (igkq ( npwx))    
  endif

  allocate (dvpsi ( npwx*npol , nbnd))    
  allocate (vlocq ( ngm , ntyp))    
  allocate (eprectot ( nbnd, nkstot) )
  allocate (eigqts ( nat))
  allocate (dmuxc (dfftp%nnr , nspin_mag , nspin_mag))    
  allocate (dvbare(dffts%nnr))    

  ALLOCATE (becp1(nksq))
  DO ik=1,nksq
     call allocate_bec_type ( nkb, nbnd, becp1(ik) )
  END DO
  CALL allocate_bec_type ( nkb, nbnd, becp )

  return
end subroutine allocate_gwq
