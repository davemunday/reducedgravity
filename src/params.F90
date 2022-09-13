! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module params
  use double

#include "CPP_OPTIONS.h"

! Specification of model parameters.

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of SWEBE to enable specification of numerical domain in terms of
!            fluid region plus overlaps.
! --------------------------------------------------------------------------------------------------
!  Variables :
!              Axy0           - amplitude of viscosity.
!              beta           - meridional gradient in Coriolis parameter.
!              cp             - specific heat capacity (at constant pressure).
!              delta_s/t      - amplitude of cross-channel salinity/potential temperature difference.
!              f0             - reference Coriolis parameter (@ x=0, y=0).
!              FWF0           - amplitude of freshwater flux forcing.
!              h0             - initial height of the free surface.
!              Kxy0           - amplitude of tracer diffusivity.
!              lambda_s/t     - reciprocal of restoring timescale, 0=no restoring.
!              Q0             - amplitude of heat flux forcing.
!              r0             - amplitude of bottom friction.
!              rho0           - Boussinesq (& EOS) reference density.
!              salty0/theta0  - reference salinity/potential temperature for EOS.
!              sbeta/talpha   - haline contraction/thermal expansion coefficient for EOS.
!              slip           - lateral momentum boundary conditions, 0=freeslip, 1=noslip.
!              taux0/tauy0    - amplitude of zonal/meridional wind stress.
!              z0             - decay scale for (exponential) structure function.
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  ! Friction parameters + wind forcing.
  real(kind=dp), parameter :: AhD0=0.025d0, A4D0=0.0d0
  real(kind=dp), parameter :: AhZ0=0.025d0, A4Z0=0.0d0
  real(kind=dp), parameter :: Kgm0=1000.0d0
  real(kind=dp), parameter :: Kv0=1.d-5, hstar=10., lambda=1.E-6, yrelax=-67.5
  real(kind=dp), parameter :: slip=1.0d0
  real(kind=dp), parameter :: taux0=0.15d0, tauy0=0.0d0

! --

  ! Input directory for pickup runs.
  character(len=5), parameter :: input_dir='data/'

  ! Output directory for data.
  character(len=5), parameter :: output_dir='data/'

! --------------------------------------------------------------------------------------------------

end module params

! --------------------------------------------------------------------------------------------------
