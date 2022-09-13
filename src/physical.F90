! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
!     Specification of physical parameters.
! --------------------------------------------------------------------------------------------------

module physical
  use double

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------
!  Variables :     (In SI units unless otherwise stated)
!                   radius          - radius of planet.
!                   g               - gravitational acceleration.
!                   pi              - 3.1415... to machine precision.
!                   omega           - rotation rate of planet.
!                   deg2rad/rad2deg - conversion between degrees and radians.
! --------------------------------------------------------------------------------------------------

  implicit none

! --

  ! Physical parameters.
  real(kind=dp), parameter :: radius=6.373d6, g=0.02d0, pi=4.0d0*ATAN(1.0d0)
  real(kind=dp), parameter :: omega=2.0d0*pi/86400.0d0
  real(kind=dp), parameter :: rho0=1000.0d0

  ! Degree/radian conversion parameters.
  real(kind=dp), parameter :: deg2rad=pi/180.0d0, rad2deg=180.0d0/pi

! --------------------------------------------------------------------------------------------------

end module physical

! --------------------------------------------------------------------------------------------------
