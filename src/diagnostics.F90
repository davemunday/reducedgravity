! --------------------------------------------------------------------------------------------------
!   D.R. MUNDAY November 2014
!   This contains the functions required to calculate any diagnostic fields from the prognostic
!   fields. Some of these fields are used to timestep the model.
! --------------------------------------------------------------------------------------------------

module diagnostics
  use double
  use overlaps

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  contains

! --------------------------------------------------------------------------------------------------
! Calculates the divergence, as per MITgcm manual. This gives viscous dissipation that is more 
! consistent with the vector invariant form of the equations of motion.

  function calculate_divergence_on_thickness( recip_area, dx, dy, uvelocity, vvelocity ) result( divergence )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)  :: divergence
    
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in)   :: recip_area
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: dy, uvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: dx, vvelocity

    ! Initialise divergence properly to prevent seg faults.
    divergence = 0.0d0

    ! Calculate divergence of the velocity as the difference between fluxes.
    divergence(1:nx,1:ny) = recip_area(1:nx,1:ny) * ( dy(2:nx+1,1:ny)*uvelocity(2:nx+1,1:ny)       &
                                                    - dy(1:nx,1:ny)*uvelocity(1:nx,1:ny)           &
                                                    + dx(1:nx,2:ny+1)*vvelocity(1:nx,2:ny+1)       &
                                                    - dx(1:nx,1:ny)*vvelocity(1:nx,1:ny) )

    ! Update the overlap regions.
    call update_overlaps_t( divergence )

  end function calculate_divergence_on_thickness

! --------------------------------------------------------------------------------------------------
! Calculates divergence of a vorticity and thickness field onto the uvelocity points.

  function calculate_divergence_on_uvelocity( recip_area, dx, dy, thickness, vorticity ) result( divergence )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly)  :: divergence

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in)   :: recip_area
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in)     :: dy, thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: dx, vorticity

    ! Initialise divergence properly to prevent seg faults.
    divergence = 0.0d0

    ! Calculate divergence as the difference between fluxes.
    divergence(1:nx+1,1:ny) = recip_area(1:nx+1,1:ny) * ( dy(1:nx+1,1:ny)*thickness(1:nx+1,1:ny)     &
                                                        - dy(0:nx,1:ny)*thickness(0:nx,1:ny)         &
                                                        + dx(1:nx+1,2:ny+1)*vorticity(1:nx+1,2:ny+1) &
                                                        - dx(1:nx+1,1:ny)*vorticity(1:nx+1,1:ny) )

    ! Update the overlap regions.
    call update_overlaps_u( divergence )

  end function calculate_divergence_on_uvelocity

! --------------------------------------------------------------------------------------------------
! Calculate divergence of a voritcity and thickness field onto the vvelocity points.

  function calculate_divergence_on_vvelocity( recip_area, dx, dy, vorticity, thickness ) result( divergence )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly)  :: divergence

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in)   :: recip_area
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in)     :: dx, thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: dy, vorticity

    ! Initialise divergence properly to prevent seg faults.
    divergence = 0.0d0

    ! Calculate divergence as the difference between fluxes.
    divergence(1:nx,1:ny+1) = recip_area(1:nx,1:ny+1) * ( dy(2:nx+1,1:ny+1)*vorticity(2:nx+1,1:ny+1) &
                                                        - dy(1:nx,1:ny+1) * vorticity(1:nx,1:ny+1)   &
                                                        + dx(1:nx,1:ny+1) * thickness(1:nx,1:ny+1)   &
                                                        - dx(1:nx,0:ny) * thickness(1:nx,0:ny) )

    ! Update the overlap regions.
    call update_overlaps_v( divergence )

  end function calculate_divergence_on_vvelocity

! --------------------------------------------------------------------------------------------------
! Calculate the kinetic energy from the provided velocity fields.
  
  function calculate_kinetic_energy( uvelocity, vvelocity ) result( kinetic_energy )
    
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly)  :: kinetic_energy

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: uvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: vvelocity

    ! Initialise speed_squared properly to prevent seg faults.
    kinetic_energy = 0.0d0

    ! Calculate speed_squared as average of squares to conserve KE.
    kinetic_energy(1:nx,1:ny) = 0.25d0*( uvelocity(1:nx,1:ny)**2 + uvelocity(2:nx+1,1:ny)**2       &
                                       + vvelocity(1:nx,1:ny)**2 + vvelocity(1:nx,2:ny+1)**2 )

    ! Update the overlap regions.
    call update_overlaps_t( kinetic_energy )

  end function calculate_kinetic_energy

! --------------------------------------------------------------------------------------------------
! Calculates transport streamfunction by integrating the zonal velocity field.

  function calculate_streamfunction( dy, thickness, uvelocity, maskavNS ) result( streamfunction )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly)   :: streamfunction

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in)   :: thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: dy, uvelocity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: maskavNS

    ! Internal arguments.
    integer(kind=dp) :: i, j

    ! Initialise streamfunction properly to prevent seg faults.
    streamfunction = 0.0d0
    
    ! Transport streamfunction on vorticity points from the u velocity.

    ! Set the boundary condition on the northern boundary.
    do i = 1,nx+1
      streamfunction(i,ny+1) = 0.0d0
    end do

    ! Integrate southwards.
    do j = ny,1,-1
      do i = 1,nx+1
        streamfunction(i,j) = streamfunction(i,j+1) + 0.5d0*( thickness(i,j) + thickness(i-1,j) )  &
                     * dy(i,j)*maskavNS(i,j)*uvelocity(i,j)
      end do
    end do

    ! Set the overlaps.
    call update_overlaps_z( streamfunction )

  end function calculate_streamfunction

! --------------------------------------------------------------------------------------------------
! Calculates the vorticity using the circulation around the grid point.

  function calculate_vorticity( recip_area, dx, dy, uvelocity, vvelocity ) result( vorticity )
    
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: vorticity

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: recip_area
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in)   :: dx, uvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in)   :: dy, vvelocity

    ! Initialise the vorticity field properly so everything is allocated properly.
    vorticity = 0.0d0

    ! Calculate the vorticity as circulation per unit area.
    vorticity(1:nx+1,1:ny+1) = recip_area(1:nx+1,1:ny+1) * (                                       &
                               uvelocity(1:nx+1,0:ny)*dx(1:nx+1,0:ny)                              &
                             + vvelocity(1:nx+1,1:ny+1)*dy(1:nx+1,1:ny+1)                          &
                             - uvelocity(1:nx+1,1:ny+1)*dx(1:nx+1,1:ny+1)                          &
                             - vvelocity(0:nx,1:ny+1)*dy(0:nx,1:ny+1) )

    ! Set the overlaps
    call update_overlaps_z( vorticity )

  end function calculate_vorticity

! --------------------------------------------------------------------------------------------------

end module diagnostics

! --------------------------------------------------------------------------------------------------
