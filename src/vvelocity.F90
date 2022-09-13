! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module vvelocity
  use double
  use domain
  use grid
  use overlaps
  use diagnostics

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  ! Interface blocks to overload common operations.

! --------------------------------------------------------------------------------------------------

  contains

! --------------------------------------------------------------------------------------------------
! Boundary conditions for no normal and/or no tangential flow on U velocity.

  subroutine apply_vvelocity_boundary_conditions( slip, newvvelocity, maskavEW, maskBC1, maskBC2, maskBC5 )

    real(kind=dp), intent(in) :: slip
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: maskavEW, maskBC1, maskBC2, maskBC5

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(inout) :: newvvelocity

    ! No normal flow BCs.
    newvvelocity = maskavEW * newvvelocity

    ! Tangential flow BCs - walls to the west.
    newvvelocity(1:nx-1,1:ny+1) = newvvelocity(1:nx-1,1:ny+1) +                                    &
         (1.0d0-2.0d0*slip)*newvvelocity(2:nx,1:ny+1)*maskBC1(2:nx,1:ny+1)

    ! Walls to the east.
    newvvelocity(2:nx,1:ny+1) = newvvelocity(2:nx,1:ny+1) +                                        &
         (1.0d0-2.0d0*slip)*newvvelocity(1:nx-1,1:ny+1)*maskBC2(1:nx-1,1:ny+1)

    ! Thin walls - overwrites the ghostpoints!
    !newvvelocity%field(2:nx-1,1:ny+1) = newvvelocity%field(2:nx-1,1:ny+1)*maskavEW%field(2:nx-1,1:ny+1) +&
    !                0.5d0*( newvvelocity%field(2:nx-1,1:ny+1)+newvvelocity%field(3:nx,1:ny+1) ) *&
    !                maskBC5%field(2:nx-1,1:ny+1)

    ! Update the overlap regions.
    call update_overlaps_v( newvvelocity )

  end subroutine apply_vvelocity_boundary_conditions

! --------------------------------------------------------------------------------------------------
! Tendency due to viscous dissipation of meridional momentum.

  function calculate_gvah( divergence_viscosity, vorticity_viscosity, recip_area, dx, dy, maskavEW,&
                           divergence, vorticity, thickness ) &
           result( tendency )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: divergence_viscosity, dx
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: divergence, thickness
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly),   intent(in) :: recip_area
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: vorticity_viscosity, dy, vorticity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly),   intent(in) :: maskavEW

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: Gvau
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: Gvav

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: tendency

    ! Initialise the tendency to prevent seg faults.
    tendency = 0.0d0
    Gvau = 0.0d0
    Gvav = 0.0d0

    ! Calculate the viscous flux of zonal momentum due to horizontal divergence.
    Gvau = divergence_viscosity * thickness * divergence
    call update_overlaps_t( Gvau )

    ! Calculate the viscous flux of zonal momentum due to horizontal vorticity.
    Gvav(1:nx+1,1:ny+1) = vorticity_viscosity(1:nx+1,1:ny+1)                                       &
                        * 0.25d0*( thickness(0:nx,0:ny) + thickness(0:nx,1:ny+1)                   &
                                 + thickness(1:nx+1,0:ny) + thickness(1:nx+1,1:ny+1) )             &
                        * vorticity(1:nx+1,1:ny+1)
    call update_overlaps_z( Gvav )
    
    ! Calculate the tendency as the divergence of the viscous fluxes.
    tendency = calculate_divergence_on_vvelocity( recip_area, dx, dy, Gvav, Gvau )

    ! Divide by thickness.
    !tendency(1:nx,1:ny+1) = 2.0d0 * maskavEW(1:nx,1:ny+1) * tendency(1:nx,1:ny+1)                  &
    !                      / ( thickness(1:nx,1:ny+1) + thickness(1:nx,0:ny) )

    ! Update the overlaps.
    call update_overlaps_v( tendency )

  end function calculate_gvah

! --------------------------------------------------------------------------------------------------
! Tendency due to Coriolis acceleration and/or nonlinear vorticity flux.

  function calculate_gvf( coriolis, thickness, uvelocity, maskavNS ) result( tendency )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: coriolis
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly),   intent(in) :: uvelocity, maskavNS

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    tendency = 0.0d0

    ! Calculate the tendency due to Coriolis acceleration - KE conserving scheme.
    !tendency(1:nx,1:ny+1) = -0.5d0 * rdyc(1:nx,1:ny+1) * (                                         &
    !                                 coriolis(1:nx,1:ny+1) *                                       &
    !                                 0.5d0*( dyg(1:nx,0:ny)*uvelocity(1:nx,0:ny)                   &
    !                                       + dyg(1:nx,1:ny+1)*uvelocity(1:nx,1:ny+1) ) *           &
    !                                 0.25d0*( thickness(0:nx-1,0:ny) + thickness(1:nx,0:ny)        &
    !                                        + thickness(0:nx-1,1:ny+1) + thickness(1:nx,1:ny+1) )  &
    !                               + coriolis(2:nx+1,1:ny+1) *                                     &
    !                                 0.5d0*( dyg(2:nx+1,0:ny)*uvelocity(2:nx+1,0:ny)               &
    !                                       + dyg(2:nx+1,1:ny+1)*uvelocity(2:nx+1,1:ny+1) ) *       &
    !                                 0.25d0 * ( thickness(1:nx,0:ny) + thickness(2:nx+1,0:ny)      &
    !                                          + thickness(1:nx,1:ny+1) + thickness(2:nx+1,1:ny+1) ) )

    tendency(1:nx,1:ny+1) = - 0.50d0 * rdyc(1:nx,1:ny+1) *                                         &
                                       ( coriolis(1:nx,1:ny+1)+coriolis(2:nx+1,1:ny+1) ) *         &
                                       0.25d0*( dyg(1:nx,1:ny+1)*uvelocity(1:nx,1:ny+1)*maskavNS(1:nx,1:ny+1)            &
                                              + dyg(2:nx+1,1:ny+1)*uvelocity(2:nx+1,1:ny+1)*maskavNS(2:nx+1,1:ny+1)        &
                                              + dyg(2:nx+1,0:ny)*uvelocity(2:nx+1,0:ny)*maskavNS(2:nx+1,0:ny)            &
                                              + dyg(1:nx,0:ny)*uvelocity(1:nx,0:ny)*maskavNS(1:nx,0:ny) )
    
    ! Divide by thickness.
    !tendency(1:nx,1:ny+1) = 2.0d0 * tendency(1:nx,1:ny+1) / ( thickness(1:nx,0:ny) + thickness(1:nx,1:ny+1) )

    ! Update the overlaps.
    call update_overlaps_v( tendency )

  end function calculate_gvf

! --------------------------------------------------------------------------------------------------
! Tendency due to gradient in geopotential height.

  function calculate_gvm( gravity, thickness ) result( tendency )

    real(kind=dp) :: gravity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: thickness

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    tendency = 0.0d0

    tendency(1:nx,1:ny+1) = -gravity*rdyc(1:nx,1:ny+1)*( thickness(1:nx,1:ny+1) - thickness(1:nx,0:ny) )

    ! Update the overlaps.
    call update_overlaps_v( tendency )

  end function calculate_gvm

! --------------------------------------------------------------------------------------------------
! Tendency due to surface wind stress.

  function calculate_gvt( wind_stress, reference_density, thickness ) result( tendency)

    real(kind=dp), intent(in) :: reference_density
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: thickness
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: wind_stress

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    tendency = 0.0d0

    ! Calculate the tendency due to wind stress.
    tendency(1:nx,1:ny+1) = 2.0d0 * wind_stress(1:nx,1:ny+1) / reference_density                           &
                          / ( thickness(1:nx,0:ny) + thickness(1:nx,1:ny+1) )

    ! Update the overlaps.
    call update_overlaps_v( tendency )

  end function calculate_gvt

! --------------------------------------------------------------------------------------------------
! Calculate the V velocity tendencies.

  function calculate_vvelocity_tendency( divergence_viscosity, vorticity_viscosity,                &
                                         coriolis, gravity, wind_stress, reference_density,        &
                                         thickness, uvelocity, vvelocity,                          &
                                         vorticity, kinetic_energy, divergence,                    &
                                         ocean, maskavNS, maskavEW )                               &
           result( vvelocity_tendency )
 
    ! Variables for viscous dissipation.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: divergence_viscosity, divergence
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: vorticity_viscosity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: Gvah

    ! Variables for Coriolis acceleration.
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: coriolis
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: Gvf

    ! Variables for gradient in kinetic energy,
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: Gvke

    ! Variables for gradient in geopotential.
    real(kind=dp), intent(in) :: gravity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: Gvm

    ! Variables for wind_stress forcing.
    real(kind=dp), intent(in) :: reference_density
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: wind_stress
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: Gvt

    ! Variables for relative vorticity flux (nonlinear advection).
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: Gvxi

    ! Model fields.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: thickness, kinetic_energy, ocean
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly),   intent(in) :: uvelocity, maskavNS
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly),   intent(in) :: vvelocity, maskavEW
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: vorticity

    ! Output variable.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: vvelocity_tendency

    ! Initialise tendencies to prevent seg faults.
    Gvah = 0.0d0
    Gvf = 0.0d0
    Gvke = 0.0d0
    Gvm = 0.0d0
    Gvt = 0.0d0
    Gvxi = 0.0d0

    ! Calculate the tendency due to viscous dissipation.
    Gvah = calculate_gvah( divergence_viscosity, vorticity_viscosity, rAs, dxf, dyu, maskavEW,     &
                           divergence, vorticity, ocean*thickness )
    Gvah(1:nx,1:ny+1) = 2.0d0 * Gvah(1:nx,1:ny+1)*maskavEW(1:nx,1:ny+1)                            &
                      / ( thickness(1:nx,1:ny+1) + thickness(1:nx,0:ny) )
    call update_overlaps_v( Gvah )

    ! Calculate the tendency due to the Coriolis acceleration
    Gvf = calculate_gvf( coriolis, thickness, uvelocity, maskavNS )
    Gvf = Gvf*maskavEW

    ! Calculate the tendency due to gradients in the kinetic energy,
    Gvke = calculate_gvm( 1.0d0, ocean*kinetic_energy )
    Gvke = Gvke*maskavEW

    ! Calculate the tendency due to gradients in the geopotential height.
    Gvm = calculate_gvm( gravity, ocean*thickness )
    Gvm = Gvm*maskavEW

    ! Calculate the tendency due to wind stress.
    Gvt = calculate_gvt( wind_stress, reference_density, thickness )
    Gvt = Gvt*maskavEW

    ! Calculate the tendency due to nonlinear vorticity flux.
    Gvxi = calculate_gvf( vorticity, thickness, uvelocity, maskavNS )
    Gvxi = Gvxi*maskavEW

    ! Sum all of the tendencies together.
    vvelocity_tendency = 0.0d0
    vvelocity_tendency = ( Gvah + Gvf + Gvke + Gvm + Gvt + Gvxi )

! --

  end function calculate_vvelocity_tendency

! --------------------------------------------------------------------------------------------------
! Timestep V velocity forwards in time.

  function step_vvelocity_forward( ab_coeffs, timestep, oldvvelocity, vvelocity_tendency )         &
           result( newvvelocity )

    real(kind=dp), dimension(1:ntl), intent(in) :: ab_coeffs
    real(kind=dp), intent(in) :: timestep
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: oldvvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly,1:ntl), intent(in) :: vvelocity_tendency

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: total_tendency
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: newvvelocity

    ! Initialise newuvelocity to prevent seg faults.
    newvvelocity = oldvvelocity

    ! Sum the different time level together weighted by ab_coeffs.
    total_tendency = ab_coeffs(1) * vvelocity_tendency(:,:,1)                                      &
         + ab_coeffs(2) * vvelocity_tendency(:,:,2)                                                &
         + ab_coeffs(3) * vvelocity_tendency(:,:,3)

    ! Extrapolate forwards in time.
    newvvelocity(1:nx,1:ny+1) = oldvvelocity(1:nx,1:ny+1) + timestep*total_tendency(1:nx,1:ny+1)

    ! Update the overlap regions.
    call update_overlaps_v( newvvelocity )

  end function step_vvelocity_forward

! --------------------------------------------------------------------------------------------------

end module vvelocity

! --------------------------------------------------------------------------------------------------
