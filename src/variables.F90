! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY May 2018
!     Module to define all of the variables.
! --------------------------------------------------------------------------------------------------

module variables
  use double
  use domain, only : ntl, nx, ny, olx, oly

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 11/05/18 - Begun line 1 re-write of reduced gravity model
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  ! Prognostic variables.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)   :: h1, h2
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: u1, u2
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: v1, v2

  ! Diagnostic variables.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: div, ke
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: psi, xi
 
  ! Tendency terms.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly,1:ntl)   :: Gh
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly,1:ntl) :: Gu
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly,1:ntl) :: Gv

  ! Masks for boundary conditions.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: tmask, tland
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: umask, uland
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: vmask, vland
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: maskavNS, maskBC3, maskBC4, maskBC6
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: maskavEW, maskBC1, maskBC2, maskBC5
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: zmask, zland

  ! Model parameters & forcing.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: f
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: Kgm, Kv
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: taux
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: tauy
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: AhD
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: AhZ

! --------------------------------------------------------------------------------------------------

end module variables

! --------------------------------------------------------------------------------------------------
