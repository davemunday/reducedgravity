! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
!     Specification of computational model domain, i.e. grid size, spacing, timestep, etc.
! --------------------------------------------------------------------------------------------------

module domain
  use double

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  ! Number of grid boxes in computational domain.
  integer(kind=dp), parameter :: nx=212_dp, ny=146_dp ! 1o
!  integer(kind=dp), parameter :: nx=424_dp, ny=288_dp ! 1/2o

  ! Overlap extent.
  integer(kind=dp), parameter :: olx=2_dp, oly=2_dp

  ! Number of time levels to use for tendency terms.
  integer(kind=dp), parameter :: ntl=3_dp

  ! Grid spacing.
  real(kind=dp), dimension(1-olx:nx+olx), parameter :: dx=1.
  real(kind=dp), dimension(1-oly:ny+oly), parameter :: dy=1.
  real(kind=dp), parameter :: xzorigin=0.0, yzorigin=-68.0d0-dy(1)

  ! Initial timestep + integration length in no. of timesteps.
  integer(kind=dp), parameter :: dtpd=144_dp, nyr=250_dp
  integer(kind=dp), parameter :: nt0=0_dp*nyr*12_dp*30_dp*dtpd, nt=nyr*12_dp*30_dp*dtpd
  integer(kind=dp), parameter :: ntd=30_dp*dtpd, ntp=nyr*12_dp*30_dp*dtpd, nts=5_dp*12_dp*30_dp*dtpd

  ! Timestep length.
  real(kind=dp), parameter :: dt=600.0d0 ! dt per day = 86400/600 = 144

! --------------------------------------------------------------------------------------------------

end module domain

! --------------------------------------------------------------------------------------------------
