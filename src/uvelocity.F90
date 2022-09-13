! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module uvelocity
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

  contains

! --------------------------------------------------------------------------------------------------
! Boundary conditions for no normal and/or no tangential flow on U velocity.

  subroutine apply_uvelocity_boundary_conditions( slip, newuvelocity, maskavNS, maskBC3, maskBC4, maskBC6 )

    real(kind=dp), intent(in) :: slip
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: maskavNS, maskBC3, maskBC4, maskBC6

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(inout) :: newuvelocity

    ! No normal flow BCs.
    newuvelocity = maskavNS * newuvelocity

    ! Tangential flow BCs - walls to the south.
    newuvelocity(1:nx+1,1:ny-1) = newuvelocity(1:nx+1,1:ny-1) +                                    &
         (1.0d0-2.0d0*slip)*newuvelocity(1:nx+1,2:ny)*maskBC3(1:nx+1,2:ny)

    ! Walls to the north.
    newuvelocity(1:nx+1,2:ny) = newuvelocity(1:nx+1,2:ny) +                                        &
         (1.0d0-2.0d0*slip)*newuvelocity(1:nx+1,1:ny-1)*maskBC4(1:nx+1,1:ny-1)

    ! Thin walls overwrites ghostpoints!
    !newuvelocity%field(1:nx+1,2:ny-1) = newuvelocity%field(1:nx+1,2:ny-1)*maskavNS%field(1:nx+1,2:ny-1) +&
    !                0.5d0*( newuvelocity%field(1:nx+1,2:ny-1)+newuvelocity%field(1:nx+1,3:ny) ) *&
    !                maskBC6%field(1:nx+1,2:ny-1)

    ! Update the overlap regions.
    call update_overlaps_u( newuvelocity )

  end subroutine apply_uvelocity_boundary_conditions

! --------------------------------------------------------------------------------------------------
! Tendency due to viscous dissipation of zonal momentum.

  function calculate_guah( divergence_viscosity, vorticity_viscosity, recip_area, dx, dy, maskavNS,&
                           divergence, vorticity, thickness ) &
           result( tendency )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: divergence_viscosity, dy
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: divergence, thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly),   intent(in) :: recip_area
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: vorticity_viscosity, dx, vorticity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly),   intent(in) :: maskavNS

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: Guau
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: Guav

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: tendency

    ! Initialise the tendency to prevent seg faults.
    tendency = 0.0d0
    Guau = 0.0d0
    Guav = 0.0d0

    ! Calculate the viscous flux of zonal momentum due to horizontal divergence.
    Guau = divergence_viscosity * thickness * divergence
    call update_overlaps_t( Guau )

    ! Calculate the viscous flux of zonal momentum due to horizontal vorticity.
    Guav(1:nx+1,1:ny+1) = -vorticity_viscosity(1:nx+1,1:ny+1)                                      &
                        * 0.25d0*( thickness(0:nx,0:ny) + thickness(0:nx,1:ny+1)                   &
                                 + thickness(1:nx+1,0:ny) + thickness(1:nx+1,1:ny+1) )             &
                        * vorticity(1:nx+1,1:ny+1)
    call update_overlaps_z( Guav )
    
    ! Calculate the tendency as the divergence of the viscous fluxes.
    tendency = calculate_divergence_on_uvelocity( recip_area, dx, dy, Guau, Guav )

    ! Divide by thickness.
    !tendency(1:nx+1,1:ny) = 2.0d0 * maskavNS(1:nx+1,1:ny) * tendency(1:nx+1,1:ny)                  &
    !                      / ( thickness(0:nx,1:ny) + thickness(1:nx+1,1:ny) )

    ! Update the overlaps.
    call update_overlaps_u( tendency )

  end function calculate_guah

! --------------------------------------------------------------------------------------------------
! Tendency due to Coriolis acceleration and/or nonlinear vorticity flux.

  function calculate_guf( coriolis, thickness, vvelocity, maskavEW ) result( tendency )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: coriolis
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: thickness
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly),   intent(in) :: vvelocity, maskavEW

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    tendency = 0.0d0

    ! Calculate the tendency due to Coriolis acceleration - KE conserving scheme.
    !tendency(1:nx+1,1:ny) = 0.5d0 * rdxc(1:nx+1,1:ny) * (                                          &
    !                                coriolis(1:nx+1,1:ny) *                                        &
    !                                0.5d0*( dxg(0:nx,1:ny)*vvelocity(0:nx,1:ny)                    &
    !                                      + dxg(1:nx+1,1:ny)*vvelocity(1:nx+1,1:ny) ) *            &
    !                                0.25d0*( thickness(0:nx,0:ny-1) + thickness(1:nx+1,0:ny-1)     &
    !                                       + thickness(0:nx,1:ny) + thickness(1:nx+1,1:ny) )       &
    !                              + coriolis(1:nx+1,2:ny+1) *                                      &
    !                                0.5d0*( dxg(0:nx,2:ny+1) * vvelocity(0:nx,2:ny+1)              &
    !                                      + dxg(1:nx+1,2:ny+1) * vvelocity(1:nx+1,2:ny+1) ) *      &
    !                                0.25d0 * ( thickness(0:nx,1:ny) + thickness(1:nx+1,1:ny)       &
    !                                         + thickness(0:nx,2:ny+1) + thickness(1:nx+1,2:ny+1) ) )

    tendency(1:nx+1,1:ny) =   0.50d0 * rdxc(1:nx+1,1:ny) *                                         &
                                       ( coriolis(1:nx+1,1:ny)+coriolis(1:nx+1,2:ny+1) ) *         &
                                       0.25d0*( dxg(1:nx+1,1:ny)*vvelocity(1:nx+1,1:ny)*maskavEW(1:nx+1,1:ny)            &
                                              + dxg(0:nx,1:ny)*vvelocity(0:nx,1:ny)*maskavEW(0:nx,1:ny)                &
                                              + dxg(1:nx+1,2:ny+1)*vvelocity(1:nx+1,2:ny+1)*maskavEW(1:nx+1,2:ny+1)        &
                                              + dxg(0:nx,2:ny+1)*vvelocity(0:nx,2:ny+1)*maskavEW(0:nx,2:ny+1) )
 
    ! Divide by thickness.
    !tendency(1:nx+1,1:ny) = 2.0d0 * tendency(1:nx+1,1:ny) / ( thickness(0:nx,1:ny) + thickness(1:nx+1,1:ny) )

    ! Update the overlaps.
    call update_overlaps_u( tendency )

  end function calculate_guf

! --------------------------------------------------------------------------------------------------
! Tendency due to gradient in geopotential height.

  function calculate_gum( gravity, thickness ) result( tendency )

    real(kind=dp) :: gravity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: thickness

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    tendency = 0.0d0

    tendency(1:nx+1,1:ny) = -gravity*rdxc(1:nx+1,1:ny)*( thickness(1:nx+1,1:ny) - thickness(0:nx,1:ny) )

    ! Update the overlaps.
    call update_overlaps_u( tendency )

  end function calculate_gum

! --------------------------------------------------------------------------------------------------
! Tendency due to surface wind stress.

  function calculate_gut( wind_stress, reference_density, thickness ) result( tendency)

    real(kind=dp), intent(in) :: reference_density
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: wind_stress

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    tendency = 0.0d0

    ! Calculate the tendency due to wind stress.
    tendency(1:nx+1,1:ny) = 2.0d0 * wind_stress(1:nx+1,1:ny) / reference_density                           &
                          / ( thickness(0:nx,1:ny) + thickness(1:nx+1,1:ny) )

    ! Update the overlaps.
    call update_overlaps_u( tendency )

  end function calculate_gut

! --------------------------------------------------------------------------------------------------
! Calculate the U velocity tendencies.

  function calculate_uvelocity_tendency( divergence_viscosity, vorticity_viscosity,                &
                                         coriolis, gravity, wind_stress, reference_density,        &
                                         thickness, uvelocity, vvelocity,                          &
                                         vorticity, kinetic_energy, divergence,                    &
                                         ocean, maskavNS, maskavEW )                               &
           result( uvelocity_tendency )

    ! Variables for viscous dissipation.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: divergence_viscosity, divergence
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: vorticity_viscosity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: Guah

    ! Variables for Coriolis acceleration.
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: coriolis
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: Guf

    ! Variables for gradient in kinetic energy.
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: Guke

    ! Variables for gradient in geopotential.
    real(kind=dp), intent(in) :: gravity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: Gum

    ! Variables for wind_stress forcing.
    real(kind=dp), intent(in) :: reference_density
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: wind_stress
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: Gut

    ! Variables for relative vorticity flux (nonlinear advection).
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: Guxi

    ! Model fields.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: thickness, kinetic_energy, ocean
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly),   intent(in) :: uvelocity, maskavNS
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly),   intent(in) :: vvelocity, maskavEW
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: vorticity

    ! Output variable.
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: uvelocity_tendency

    ! Initialise tendencies to prevent seg faults.
    Guah = 0.0d0
    Guf = 0.0d0
    Guke = 0.0d0
    Gum = 0.0d0
    Gut = 0.0d0
    Guxi = 0.0d0

    ! Calculate the tendency due to viscous dissipation.
    Guah = calculate_guah( divergence_viscosity, vorticity_viscosity, rAw, dxv, dyf, maskavNS,     &
                           divergence, vorticity, ocean*thickness )
    Guah(1:nx+1,1:ny) = 2.0d0 * Guah(1:nx+1,1:ny)*maskavNS(1:nx+1,1:ny)                            &
                      / ( thickness(0:nx,1:ny) + thickness(1:nx+1,1:ny) )
    call update_overlaps_u( Guah )

    ! Calculate the tendency due to the Coriols acceleration.
    Guf = calculate_guf( coriolis, thickness, vvelocity, maskavEW )
    Guf = Guf*maskavNS

    ! Calculate the tendency due to gradients in kinetic energy.
    Guke = calculate_gum( 1.0d0, ocean*kinetic_energy )
    Guke = Guke*maskavNS

    ! Calculate the tendency due to gradients in the geopotential height.
    Gum = calculate_gum( gravity, ocean*thickness )
    Gum = Gum*maskavNS

    ! Calculate the tendency due to wind stress.
    Gut = calculate_gut( wind_stress, reference_density, thickness )
    Gut = Gut*maskavNS

    ! Calculate the tendency due to nonlinear vorticity flux.
    Guxi = calculate_guf( vorticity, thickness, vvelocity, maskavEW )
    Guxi = Guxi*maskavNS

! --

    ! Sum all of the tendencies together.
    uvelocity_tendency = 0.0d0
    uvelocity_tendency = ( Guah + Guf + Guke + Gum + Gut + Guxi )

! --

  end function calculate_uvelocity_tendency

! --------------------------------------------------------------------------------------------------
! Timestep U velocity forwards in time.

  function step_uvelocity_forward( ab_coeffs, timestep, olduvelocity, uvelocity_tendency )         &
           result( newuvelocity )

    real(kind=dp), dimension(1:ntl), intent(in) :: ab_coeffs
    real(kind=dp), intent(in) :: timestep
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: olduvelocity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly,1:ntl), intent(in) :: uvelocity_tendency

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: total_tendency
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: newuvelocity

    ! Initialise newuvelocity to prevent seg faults.
    newuvelocity = olduvelocity

    ! Sum the different time level together weighted by ab_coeffs.
    total_tendency = ab_coeffs(1) * uvelocity_tendency(:,:,1)                                      &
         + ab_coeffs(2) * uvelocity_tendency(:,:,2)                                                &
         + ab_coeffs(3) * uvelocity_tendency(:,:,3)

    ! Extrapolate forwards in time.
    newuvelocity(1:nx+1,1:ny) = olduvelocity(1:nx+1,1:ny) + timestep * total_tendency(1:nx+1,1:ny)

    ! Update the overlap regions.
    call update_overlaps_u( newuvelocity )

  end function step_uvelocity_forward

! --------------------------------------------------------------------------------------------------

end module uvelocity

! --------------------------------------------------------------------------------------------------
