! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module grid
  use double
  use domain, only : nx, ny, olx, oly

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 11/05/18 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none

! --

  ! Grid point locations.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: xc, yc     ! h points/cell centres.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: xf, yf     ! u points.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: xg, yg     ! v points.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: xz, yz     ! xi points/cell corners.

  ! Grid spacings for different types of cell.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: dxf, dyf   ! distance across h cell.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: dxc, dyg   ! distance across u cell.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: dxg, dyc   ! distance across v cell.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: dxv, dyu   ! distance across xi cell.

  ! Grid cell surface areas.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: Ac         ! dxg*dyg, h cells.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: Aw         ! dxv*dyf, u cells.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: As         ! dxf*dyu, v cells.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: Az         ! dxc*dyc, xi cells.

  ! Reciprocals of grid spacings.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: rdxf, rdyf ! reciprocal of dxf, dyf, h cell.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: rdxc, rdyg ! reciprocal of dxc, dyg, u cell.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: rdxg, rdyc ! reciprocal of dxg, dyc, v cell.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: rdxv, rdyu ! reciprocal of dxv, dyu, xi cell.

  ! Reciprocals of surface areas.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: rAc        ! reciprocal of dxf*dyf, h cells.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)   :: rAw        ! reciprocal of dxg*dyg, u cells.
  real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)   :: rAs        ! reciprocal of dxg*dyc, v cells.
  real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: rAz        ! reciprocal of dxv*dyu, xi cells.

! --------------------------------------------------------------------------------------------------

end module grid

! --------------------------------------------------------------------------------------------------
