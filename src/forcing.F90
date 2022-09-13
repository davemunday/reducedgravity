! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! This module contains subroutines that will usually be problem specific, and thus user-defined. In 
! other words; initialisation of friction and wind forcing. The coefficients are all defined on the 
! thickness points. In the timestepping routines, care is taken to ensure that they are correctly
! averaged to the velocity points, when needed. This allows for simple implementation of spatially 
! varying friction, etc. In order to make user modification of this routines to create 
! spatially-varying fields, this module already uses grid, although it is not strictly necessary for
! the default of constant fields. Any analytical formulas used to calculate variable coefficients 
! should use xc & yc (the x-/y- coordinates of cell centres/thickness points).
! --------------------------------------------------------------------------------------------------

module forcing
  use double
  use domain
  use grid
  use overlaps
  use physical

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 11/05/18 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  contains

! --------------------------------------------------------------------------------------------------

  function set_bottom_friction( bottom_friction_magnitude ) result( coefficient_of_friction )

    real(kind=dp), intent(in) :: bottom_friction_magnitude

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: coefficient_of_friction

    integer(kind=dp) :: k

    ! Calculate the bottom/interfacial friction.
    coefficient_of_friction = bottom_friction_magnitude

  end function set_bottom_friction

! --------------------------------------------------------------------------------------------------

  function set_lateral_diffusion( diffusivity_magnitude ) result( diffusivity )

    real(kind=dp), intent(in) :: diffusivity_magnitude

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: diffusivity

    ! Calculate the diffusivity.
    diffusivity = diffusivity_magnitude
    
  end function set_lateral_diffusion

! --------------------------------------------------------------------------------------------------

  function set_biharmonic_diffusion( diffusivity_magnitude ) result( diffusivity )

    real(kind=dp), intent(in) :: diffusivity_magnitude

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: diffusivity

    ! Calculate the diffusivity.
    diffusivity = diffusivity_magnitude
    
  end function set_biharmonic_diffusion

! --------------------------------------------------------------------------------------------------

  subroutine set_lateral_friction( divergence_viscosity_magnitude, vorticity_viscosity_magnitude,  &
                                   divergence_viscosity, vorticity_viscosity )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)     :: LD
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: LZ

    real(kind=dp), intent(in) :: divergence_viscosity_magnitude
    real(kind=dp), intent(in) :: vorticity_viscosity_magnitude

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out)     :: divergence_viscosity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(out) :: vorticity_viscosity

    ! Initialise the coefficient_of_friction to prevent seg faults.
    divergence_viscosity = divergence_viscosity_magnitude
    vorticity_viscosity = vorticity_viscosity_magnitude
   
    ! Calculate the length scales.
    LD = dxf*dyf/( dxf+dyf )
    LZ = dxv*dyu/( dxv+dyu )

    ! Calculate the viscosity.
    divergence_viscosity(1:nx,1:ny) = 0.25*divergence_viscosity_magnitude*Ac(1:nx,1:ny)/dt
    vorticity_viscosity(1:nx+1,1:ny+1) = 0.25*vorticity_viscosity_magnitude*Az(1:nx+1,1:ny+1)/dt
    
    ! Update the overlaps.
    call update_overlaps_t( divergence_viscosity )
    call update_overlaps_z( vorticity_viscosity )
    
  end subroutine set_lateral_friction

! --------------------------------------------------------------------------------------------------

  subroutine set_biharmonic_friction( divergence_viscosity_magnitude, vorticity_viscosity_magnitude,&
                                      divergence_viscosity, vorticity_viscosity )

    real(kind=dp), intent(in) :: divergence_viscosity_magnitude
    real(kind=dp), intent(in) :: vorticity_viscosity_magnitude

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out)     :: divergence_viscosity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(out) :: vorticity_viscosity

    ! Initialise the coefficient_of_friction to prevent seg faults.
    divergence_viscosity = divergence_viscosity_magnitude
    vorticity_viscosity = vorticity_viscosity_magnitude

  end subroutine set_biharmonic_friction

! --------------------------------------------------------------------------------------------------

  subroutine set_wind( taux_magnitude, tauy_magnitude, taux, tauy )

    real(kind=dp), intent(in) :: taux_magnitude, tauy_magnitude

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(out) :: taux
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(out) :: tauy

    integer(kind=dp) :: j

    ! Initialise the wind stresses to prevent seg faults.
    taux = 0.0d0
    tauy = 0.0d0

    DO j=1,ny
      IF( yf(0,j).LE.-33. )THEN
        taux(:,j) = taux_magnitude*cos(pi*((yf(:,j)+50.5)/35.))*cos(pi*((yf(:,j)+50.5)/35.))
      ELSE
        taux(:,j) = 0.
      ENDIF
    END DO
!    PRINT*,taux(1,:)

  end subroutine set_wind

! --------------------------------------------------------------------------------------------------

end module forcing

! --------------------------------------------------------------------------------------------------
