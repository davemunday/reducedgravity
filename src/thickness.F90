! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module thickness
  use double
  use domain
  use grid
  use physical, only : g, pi
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
! Apply boundary conditions to the thickness field by updating the overlaps and fixing any land
! points to their initial values.

  subroutine apply_thickness_boundary_conditions( oldthickness, newthickness, ocean_mask, land_mask )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: oldthickness, ocean_mask, land_mask

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(inout) :: newthickness

    ! Thickness outside domain fixed at initial value, to prevent blow up when ghost point velocities 
    ! (from no-slip conditions) lead to non-zero advection terms!
    newthickness(1:nx,1:ny) = newthickness(1:nx,1:ny)*ocean_mask(1:nx,1:ny)                        &
         + oldthickness(1:nx,1:ny)*land_mask(1:nx,1:ny)

    ! Update the overlap region.
    call update_overlaps_t( newthickness )

  end subroutine apply_thickness_boundary_conditions

! --------------------------------------------------------------------------------------------------
! Calcuate the zonal thickness gradient.

  function calculate_dhdx( thickness ) result( dhdx )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: thickness

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: dhdx

    ! Initialise the gradient.
    dhdx = 0.0d0

    ! Take the zonal gradient.
    dhdx(1:nx+1,1:ny) = ( thickness(1:nx+1,1:ny) - thickness(0:nx,1:ny) ) / dxc(1:nx+1,1:ny)

    ! Update the overlaps.
    call update_overlaps_u( dhdx )
    
  end function calculate_dhdx
  
! --------------------------------------------------------------------------------------------------
! Calcuate the meridional thickness gradient.

  function calculate_dhdy( thickness ) result( dhdy )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: thickness

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: dhdy

    ! Initialise the gradient.
    dhdy = 0.0d0

    ! Take the zonal gradient.
    dhdy(1:nx,1:ny+1) = ( thickness(1:nx,1:ny+1) - thickness(1:nx,1:ny+1) ) / dyc(1:nx,1:ny+1)

    ! Update the overlaps.
    call update_overlaps_v( dhdy )
    
  end function calculate_dhdy
  
! --------------------------------------------------------------------------------------------------
!  Calculate the thickness tendency due to thickness diffusion.

  function calculate_Ghkgm( kappa_gm, thickness, maskavNS, maskavEW ) result( tendency )
  
!  --

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: kappa_gm, thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: maskavNS
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: maskavEW

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: ghu
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: ghv

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: tendency

!  --

    ! Initialise tendency term.
    tendency = 0.0

!  --

    ! Calculate the zonal thickness flux due to GM.
    ghu = calculate_dhdx( thickness )
    ghu(1:nx+1,1:ny) = 0.5d0 * ( kappa_gm(1:nx+1,1:ny)+kappa_gm(0:nx,1:ny) ) *                     &
                       maskavNS(1:nx+1,1:ny) * ghu(1:nx+1,1:ny)
    call update_overlaps_u( ghu )

    ! Calculate the meridional thickness flux due to GM.
    ghv = calculate_dhdy( thickness )
    ghv(1:nx,1:ny+1) = 0.5d0 * ( kappa_gm(1:nx,1:ny+1)+kappa_gm(1:nx,0:ny) ) *                     &
                       maskavEW(1:nx,1:ny+1) * ghv(1:nx,1:ny+1)
    call update_overlaps_v( ghv )

    ! Take the divergence to form the tendency.
    tendency = calculate_divergence_on_thickness( rAc, dxg, dyg, ghu, ghv )
    
!  --

  end function calculate_Ghkgm

! --------------------------------------------------------------------------------------------------
!  Calculate the thickness tendency due to diapycnal diffusion.

  function calculate_Ghkv( kappa_v, thickness ) result( tendency )
  
!  --

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: kappa_v, thickness

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: tendency

!  --

    ! Initialise tendency term.
    tendency = 0.0

!  --

    tendency(1:nx,1:ny) = kappa_v(1:nx,1:ny) / thickness(1:nx,1:ny)

    ! Update the overlap regions.
    call update_overlaps_t( tendency )

!  --

  end function calculate_Ghkv

! --------------------------------------------------------------------------------------------------
!  Calculate the thickness tendency due to diabatic forcing.

  function calculate_Ghrelax( relax_boundary, relax_thickness, relax_timescale,                    &
                              thickness )                                                          &
                              result( tendency )

!  --

    real(kind=dp), intent(in) :: relax_boundary, relax_thickness, relax_timescale
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: thickness

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: tendency
 
!  --

    ! Initialise tendency term.
    tendency = 0.0
 
!  --

    tendency(1:nx,1:ny) = 0.5d0 * relax_timescale *                                                &
                          ( 1.0d0-tanh( 0.75d0*pi*(yc(1:nx,1:ny)-relax_boundary) ) ) *             &
                          ( relax_thickness-thickness(1:nx,1:ny) )

    ! Update overlaps.
    call update_overlaps_t( tendency )
    !print*, 'tanh =', 0.5d0*( 1.0d0 - tanh( 0.75d0*pi*(yc(1,1:ny)-relax_boundary ) ) )
    
!  --

  end function calculate_Ghrelax

! --------------------------------------------------------------------------------------------------
! Calculate the thickness flux at U velocity points.

  function calculate_ghu( oldthickness, olduvelocity, maskavNS ) result( thickness_flux )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: oldthickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: olduvelocity, maskavNS

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: thickness_flux

#ifdef FLUX_LIMIT_ADVECTION
    integer(kind=dp) :: i, j
    real(kind=dp) :: Cr, Rjp, Rj, Rjm, upstream_bias
#endif

    ! Initialise the thickness_flux to prevent seg faults.
    thickness_flux = 0.0d0

    ! Calculate the thickness_flux using 2nd order centred differences.
    thickness_flux(1:nx+1,1:ny) = maskavNS(1:nx+1,1:ny) *                                          &
                        0.5d0*( oldthickness(0:nx,1:ny) + oldthickness(1:nx+1,1:ny) ) *            &
                        olduvelocity(1:nx+1,1:ny)

#ifdef FLUX_LIMIT_ADVECTION
    ! Try to flux limit the above expression - see MITgcm/pkg/generic_advdiff/gad_fluxlimit_adv_x.F
    do j = 1,ny
      do i = 1,nx+1

        ! Calculate the slope in different directions.
        Rjp = maskavNS(i+1,j) * ( oldthickness(i+1,j) - oldthickness( i ,j) ) ! to the east
        Rj  = maskavNS( i ,j) * ( oldthickness( i ,j) - oldthickness(i-1,j) ) ! at the point
        Rjm = maskavNS(i-1,j) * ( oldthickness(i-1,j) - oldthickness(i-2,j) ) ! to the west

        ! Calculate the slope ratio, taking into account the sign of the advection.
        if( Rj/=0.0d0 )then
          if( olduvelocity(i,j) > 0. )then
            Cr = Rjm/Rj
          else
            Cr = Rjp/Rj
          end if
        else
          if( olduvelocity(i,j) > 0. )then
            Cr = Rjm*1.0d20
          else
            Cr = Rjp*1.0d20
          end if
        end if

        ! Apply the "Superbee" limiter of Roe (1985).
        Cr = max( 0.0d0, max( min(1.0d0,2.0d0*Cr), min(2.0d0,Cr) ) )

        ! Calculate the upstream-biased flux.
        upstream_bias = thickness_flux(i,j) - 0.5d0*abs(olduvelocity(i,j))*Rj

        ! Flux limit the initial centred-in-space flux estimate.
        thickness_flux(i,j) = upstream_bias + Cr*( thickness_flux(i,j) - upstream_bias )

      end do
    end do
#endif

    ! Update the overlap regions.
    call update_overlaps_u( thickness_flux )

  end function calculate_ghu

! --------------------------------------------------------------------------------------------------
! Calculate the thickness flux at V velocity points.

  function calculate_ghv( oldthickness, oldvvelocity, maskavEW ) result( thickness_flux )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: oldthickness
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: oldvvelocity, maskavEW

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: thickness_flux

#ifdef FLUX_LIMIT_ADVECTION
    integer(kind=dp) :: i, j
    real(kind=dp) :: Cr, Rjp, Rj, Rjm, upstream_bias
#endif

    ! Initialise the thickness_flux to prevent seg faults.
    thickness_flux = 0.0d0

    ! Calculate the thickness_flux using 2nd order centred differences.
    thickness_flux(1:nx,1:ny+1) = maskavEW(1:nx,1:ny+1) *                                          &
                        0.5d0*( oldthickness(1:nx,0:ny) + oldthickness(1:nx,1:ny+1) ) *            &
                        oldvvelocity(1:nx,1:ny+1)

#ifdef FLUX_LIMIT_ADVECTION
    ! Try to flux limit the above expression - see MITgcm/pkg/generic_advdiff/gad_fluxlimit_adv_y.F
    do j = 1,ny+1
      do i = 1,nx

        ! Calculate the slope in different direction.
        Rjp = maskavEW(i,j+1) * ( oldthickness(i,j+1) - oldthickness(i, j ) ) ! to the north
        Rj  = maskavEW(i, j ) * ( oldthickness(i, j ) - oldthickness(i,j-1) ) ! at the point
        Rjm = maskavEW(i,j-1) * ( oldthickness(i,j-1) - oldthickness(i,j-2) ) ! to the south

        ! Calculate the slope ratio, taking into account the sign of the advection.
        if( Rj/=0.0d0 )then
          if( oldvvelocity(i,j) > 0. )then
            Cr = Rjm/Rj
          else
            Cr = Rjp/Rj
          end if
        else
          if( oldvvelocity(i,j) > 0. )then
            Cr = Rjm*1.0d20
          else
            Cr = Rjp*1.0d20
          end if
        end if

        ! Apply the "Superbee" limiter of Roe (1985).
        Cr = max( 0.0d0, max( min(1.0d0,2.0d0*Cr), min(2.0d0,Cr) ) )

        ! Calculate the upstream-biased flux.
        upstream_bias = thickness_flux(i,j) - 0.5d0*abs(oldvvelocity(i,j))*Rj

        ! Flux limit the initial centred-in-space flux estimate.
        thickness_flux(i,j) = upstream_bias + Cr*( thickness_flux(i,j) - upstream_bias )

      end do
    end do
#endif

    ! Update the overlap regions.
    call update_overlaps_v( thickness_flux )

  end function calculate_ghv

! --------------------------------------------------------------------------------------------------

  function calculate_thickness_advection( oldthickness, olduvelocity, oldvvelocity, maskavNS, maskavEW ) &
           result( tendency )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: oldthickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: olduvelocity, maskavNS
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: oldvvelocity, maskavEW

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly) :: Ghu
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly) :: Ghv

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    tendency = 0.0d0

    ! Calculate the thickness fluxes at the U/V velocity points.
    Ghu = calculate_ghu( oldthickness, olduvelocity, maskavNS )
    Ghv = calculate_ghv( oldthickness, oldvvelocity, maskavEW )

    ! Take the divergence of the thickness fluxes.
    tendency = calculate_divergence_on_thickness( rAc, dxg, dyg, Ghu, Ghv )
    tendency = -tendency ! flip sign for tendency.

  end function calculate_thickness_advection

! --------------------------------------------------------------------------------------------------

  function calculate_thickness_tendency( diapycnal_diffusion, kappa_gm,                            &
                                         relax_boundary, relax_thickness, relax_timescale, ocean,  &
                                         thickness, uvelocity, vvelocity, maskavNS, maskavEW )     &
           result( tendency )

    ! Variables for thickness advection.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: thickness, ocean
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: uvelocity, maskavNS
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: vvelocity, maskavEW
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: Ghuv

    ! Variables for diapycnal diffusion.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: diapycnal_diffusion
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: Ghkv

    ! Variables for thickness relaxtion.
    real(kind=dp), intent(in) :: relax_boundary, relax_thickness, relax_timescale
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: Ghrelax
   
    ! Variables for GM thickness advection.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: kappa_gm
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: Ghgm
   
    ! Output variable.
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: tendency

    ! Initialise tendency to prevent seg faults.
    Ghuv = 0.0d0
    Ghkv = 0.0d0
    Ghgm = 0.0d0

    ! Calculate the tendency due to the advection of thickness.
    Ghuv = calculate_thickness_advection( thickness, uvelocity, vvelocity, maskavNS, maskavEW )
    Ghuv = ocean*Ghuv

    ! Calculate the tendency due to diapycnal diffusion.
    Ghkv = calculate_Ghkv( diapycnal_diffusion, thickness )

    ! Calculate the tendency due to thickness relaxtion.
    Ghrelax = calculate_Ghrelax( relax_boundary, relax_thickness, relax_timescale, thickness )

    ! Calculate the tendency due to thickness diffusion ala GM.
    Ghgm = calculate_Ghkgm( kappa_gm, thickness, maskavNS, maskavEW )

    tendency = ( Ghuv + Ghkv + Ghrelax + Ghgm )

    ! Update the overlap region.
    call update_overlaps_t( tendency )

  end function calculate_thickness_tendency

! --------------------------------------------------------------------------------------------------
! Timestep layer thickness forwards in time.

  function step_thickness_forward( ab_coeffs, timestep, oldthickness, thickness_tendency,          &
                                   relax_thickness ) result( newthickness )

    real(kind=dp), dimension(1:ntl), intent(in) :: ab_coeffs
    real(kind=dp), intent(in) :: timestep, relax_thickness
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: oldthickness
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly,1:ntl) :: thickness_tendency

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: total_tendency
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly) :: newthickness

    ! Initialise newthickness to prevent seg faults.
    newthickness = oldthickness
    total_tendency = 0.0d0

    ! Sum the different time level together weighted by ab_coeffs.
    total_tendency = ab_coeffs(1) * thickness_tendency(:,:,1)                                      &
         + ab_coeffs(2) * thickness_tendency(:,:,2)                                                &
         + ab_coeffs(3) * thickness_tendency(:,:,3)

    ! Extrapolate forwards in time.
    newthickness(1:nx,1:ny) = oldthickness(1:nx,1:ny) + timestep * total_tendency(1:nx,1:ny)

    ! Enforce a minimum thickness to prevent NaN's.
    newthickness = max( newthickness, relax_thickness )

    ! Update the overlap region.
    call update_overlaps_t( newthickness )

  end function step_thickness_forward

! --------------------------------------------------------------------------------------------------

end module thickness

! --------------------------------------------------------------------------------------------------
