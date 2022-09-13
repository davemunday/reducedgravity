! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
!     This module contains the subroutines required to initialise the model.
! --------------------------------------------------------------------------------------------------

module initial
  use double
  use physical
  use domain
  use grid
  use overlaps

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none
  
! --------------------------------------------------------------------------------------------------

  contains

! --------------------------------------------------------------------------------------------------
! Initialises the shape of the ocean basins according to some idealised layout.

    subroutine init_basin( ocean_mask, land_mask )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: ocean_mask, land_mask

    ! Use idealised bathymetry.

    ! Make the whole thing ocean.
    ocean_mask = 1.0d0

    ! Put land on all four boundaries
    ocean_mask(1:nx,1) = 0.0d0
    ocean_mask(1:nx,ny) = 0.0d0
    !ocean_mask(1,1:ny) = 0.0d0
    !ocean_mask(nx,1:ny) = 0.0d0

    ! Use an if loop, which tests the number of grid boxes, to place the partial walls.
    if( ny == 146_dp )then
       print*, 'ny =', ny
       ! Add partial walls to the east and west - Cape Horn at 32oS.
       ocean_mask(1:8,1+37:) = 0
       ocean_mask(205:nx,1+37:) = 0
    
       ! Divide into Pacific and Atlantic basins - narrow Drake Passage.
       ocean_mask(129:144,1+11:) = 0

       ! Add the Antipodean landmass to the Indo-Pacific basin.
       !ocean_mask(59:78,1+1:1+24) = 0
    elseif( ny == 288_dp )then
       print*, 'initial.F90 : ny =', ny
       stop
    else
    endif
    
    ! Set the overlaps.
    call update_overlaps_t( ocean_mask )

    ! Specify the land_mask.
    land_mask = 1.0 - ocean_mask

  end subroutine init_basin

! --------------------------------------------------------------------------------------------------
!  Initialise no normal flow boundary conditions on u and v

  subroutine init_bc( olduvelocity, newuvelocity, oldvvelocity, newvvelocity, maskavNS, maskavEW )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: maskavNS
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: maskavEW

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(inout) :: olduvelocity, newuvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(inout) :: oldvvelocity, newvvelocity

    ! Apply no normal flow boundary condition to the old and new U velocity fields.
    olduvelocity = maskavNS * olduvelocity
    newuvelocity = maskavNS * newuvelocity

    ! Apply no normal flow boundary condition to the old and new V velocity fields.
    oldvvelocity = maskavEW * oldvvelocity
    newvvelocity = maskavEW * newvvelocity

  end subroutine init_bc

! --------------------------------------------------------------------------------------------------
! Calculates the Coriolis parameter at velocity points.

  function init_coriolis() result( coriolis )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly) :: coriolis

    ! Initialise the coriolis output.
    coriolis = 0.0d0

    ! Set the Coriolis parameter for a spherical polar grid.
    coriolis = 2.0d0*omega*SIN( deg2rad*yz )

  end function init_coriolis
  
! --------------------------------------------------------------------------------------------------
! Initialises diagnostic fields to sensible things.

  subroutine init_diagnostic_fields( divergence, kinetic, streamfunction, vorticity )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: divergence, kinetic
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(out) :: streamfunction, vorticity

    divergence = 0.0d0
    kinetic = 0.0d0
    streamfunction = 0.0d0
    vorticity = 0.0d0
    
  end subroutine init_diagnostic_fields
    
! --------------------------------------------------------------------------------------------------
! Sets up grid spacing and defines x and y locations at thickness and velocity points. Defines 
! commonly used arrays (reciprocals of average box sizes, etc).

  subroutine init_grid

    integer(kind=dp) :: i, j
    real(kind=dp), dimension(1:nx+1,1:ny+1) :: dlat, dlon

! --

    ! Coordinates of xi points/cell corners.

    ! Initialise to zero.
    xz = 0.0d0
    yz = 0.0d0

    ! x coordinates of cell corners.
    xz(1-olx,1-oly:ny+1+oly) = xzorigin - sum( dx(1-olx:0) )
    do i = 1-olx+1,nx+1+olx
      xz(i,1-oly:ny+1+oly) = xz(i-1,1-oly:ny+1+oly) + dx(i-1)
    end do

    ! y coordinates of cell corners.
    yz(1-olx:nx+1+olx,1-oly) = yzorigin - sum( dy(1-olx:0) )
    do j = 1-oly+1,ny+1+oly
       yz(1-olx:nx+1+olx,j) = yz(1-olx:nx+1+olx,j-1) + dy(j-1)
    end do

! --

    ! Coordinates of th points/cell centres.

    ! Initialise to zero.
    xc = 0.0d0
    yc = 0.0d0

    ! x coordinates of cell centres.
    xc(1-olx,1-oly:ny+oly) = xzorigin - sum( dx(1-olx:0) ) + 0.5d0*dx(1-olx)
    do i = 1-olx+1,nx+olx
      xc(i,1-oly:ny+oly) = xc(i-1,1-oly:ny+oly) + 0.5d0*( dx(i-1)+dx(i) )
    end do

    ! y coordinates of cell centres.
    yc(1-olx:nx+olx,1-oly) = yzorigin - sum( dy(1-oly:0) ) + 0.5d0*dy(1-oly)
    do j = 1-oly+1,ny+oly
      yc(1-olx:nx+olx,j) = yc(1-olx:nx+olx,j-1) + 0.5d0*( dy(j-1)+dy(j) )
    end do

! --

    ! Coordinates of u velocity points.

    ! Initialise to zero.
    xf = 0.0d0
    yf = 0.0d0

    ! x coordinates of u velocity points.
    xf(1-olx,1-oly:ny+oly) = xzorigin - sum( dx(1-olx:0) )
    do i = 1-olx+1,nx+1+olx
      xf(i,1-oly:ny+oly) = xf(i-1,1-oly:ny+oly) + dx(i-1)
    end do

    ! y coordinates of u velocity points.
    yf(1-olx:nx+olx,1-oly) = yzorigin - sum( dy(1-oly:0) ) + 0.5d0*dy(1-oly)
    do j = 1-oly+1,ny+oly
      yf(1-olx:nx+1+olx,j) = yf(1-olx:nx+1+olx,j-1) + 0.5d0*( dy(j-1)+dy(j) )
    end do

! --

    ! Coordinates of v velocity points.

    ! Initialise to zero.
    xg = 0.0d0
    yg = 0.0d0

    ! x coordinates of v velocity points.
    xg(1-olx,1-oly:ny+1+oly) = xzorigin - sum( dx(1-olx:0) ) + 0.5d0*dx(1-olx)
    do i = 1-olx+1,nx+olx
      xg(i,1-oly:ny+1+oly) = xg(i-1,1-oly:ny+1+oly) + 0.5d0*( dx(i-1)+dx(i) )
    end do

    ! y coordinates of v velocity points.
    yg(1-olx:nx+olx,1-oly) = yzorigin - sum( dy(1-oly:0) )
    do j = 1-oly+1,ny+1+oly
      yg(1-olx:nx+olx,j) = yg(1-olx:nx+olx,j-1) + dy(j-1)
    end do

! --

    ! Initialise to zero.
    dxf = 0.0d0
    dxg = 0.0d0
    dyf = 0.0d0
    dyg = 0.0d0
    
    ! Grid spacings in x.
    do j = 1,ny,1
      dxf(1:nx,j) = radius*deg2rad*dx(1:nx)*cos( deg2rad*yc(1:nx,j) )
    end do
    do j = 1,ny+1,1
      dxg(1:nx,j) = radius*deg2rad*dx(1:nx)*cos( deg2rad*yg(1:nx,j) )
    end do

    ! Grid spacings in y.
    do i = 1,nx,1
      dyf(i,1:ny) = radius*deg2rad*dy(1:ny)
    end do
    do i = 1,nx+1,1
      dyg(i,1:ny) = radius*deg2rad*dy(1:ny)
    end do

! --

    ! Set the overlaps.
    call update_overlaps_t( dxf )
    call update_overlaps_v( dxg )
    call update_overlaps_t( dyf )
    call update_overlaps_u( dyg )

! --

    ! Initialise to zero.
    dxc = 0.0d0
    dxv = 0.0d0
    dyc = 0.0d0
    dyu = 0.0d0
    
    ! Average the known grid spacings onto the unknown grid spacings.
    dxc(1:nx+1,1:ny) = 0.5d0*( dxf(0:nx,1:ny)+dxf(1:nx+1,1:ny) )
    dxv(1:nx+1,1:ny+1) = 0.5d0*( dxg(0:nx,1:ny+1)+dxg(1:nx+1,1:ny+1) )

    dyc(1:nx,1:ny+1) = 0.5d0*( dyf(1:nx,0:ny)+dyf(1:nx,1:ny+1) )
    dyu(1:nx+1,1:ny+1) = 0.5d0*( dyg(1:nx+1,0:ny)+dyg(1:nx+1,1:ny+1) )

! --

    ! Set the overlaps.
    call update_overlaps_u( dxc )
    call update_overlaps_z( dxv )
    call update_overlaps_v( dyc )
    call update_overlaps_z( dyu )

! --

    ! Surface area of tracer cells.
    dlon(1:nx,1:ny) = rad2deg*dxg(1:nx,1:ny)/radius
    dlat(1:nx,1:ny) = rad2deg*dyg(1:nx,1:ny)/radius
    Ac(1:nx,1:ny) = deg2rad*dlon(1:nx,1:ny)*radius*radius*                                          &
      abs( sin( deg2rad*(yz(1:nx,1:ny)+dlat(1:nx,1:ny)) )-sin( deg2rad*yz(1:nx,1:ny) ) )

    ! Surface area of vorticity cells.
    dlon(1:nx+1,1:ny+1) = rad2deg*dxc(0:nx,0:ny)/radius
    dlat(1:nx+1,1:ny+1) = rad2deg*dyc(0:nx,0:ny)/radius
    Az(1:nx+1,1:ny+1) = deg2rad*dlon(1:nx+1,1:ny+1)*radius*radius*                                 &
      abs( sin( deg2rad*(yc(0:nx,0:ny)+dlat) )-sin( deg2rad*yc(0:nx,0:ny) ) )

    ! Surface area of u velocity cells.
    dlon(1:nx+1,1:ny) = rad2deg*dxv(0:nx,1:ny)/radius
    dlat(1:nx+1,1:ny) = rad2deg*dyf(0:nx,1:ny)/radius
    Aw(1:nx+1,1:ny) = deg2rad*dlon(1:nx+1,1:ny)*radius*radius*                                     &
      abs( sin( deg2rad*(yg(0:nx,1:ny)+dlat(1:nx+1,1:ny)) )-sin( deg2rad*yg(0:nx,1:ny) ) )

    ! Surface area of v velocity cells.
    dlon(1:nx,1:ny+1) = rad2deg*dxf(1:nx,0:ny)/radius
    dlat(1:nx,1:ny+1) = rad2deg*dyu(1:nx,0:ny)/radius
    As(1:nx,1:ny+1) = deg2rad*dlon(1:nx,1:ny+1)*radius*radius*                                     &
      abs( sin( deg2rad*(yf(1:nx,0:ny)+dlat(1:nx,1:ny+1)) )-sin( deg2rad*yf(1:nx,0:ny)) )

! --

    ! Set the overlaps.
    call update_overlaps_t( Ac )
    call update_overlaps_u( Aw )
    call update_overlaps_v( As )
    call update_overlaps_z( Az )

! --

    ! The way the reciprocals of the grid geometry parameters is calculated is
    ! independent of the type of grid (Cartesian vs. spherical polar, etc).

    ! Reciprocal of grid spacings in x.
    rdxg = 1.0d0 / dxg
    rdxc = 1.0d0 / dxc
    rdxv = 1.0d0 / dxv
    rdxf = 1.0d0 / dxf

    ! Reciprocal of grid spacings in y.
    rdyg = 1.0d0 / dyg
    rdyc = 1.0d0 / dyc
    rdyf = 1.0d0 / dyf
    rdyu = 1.0d0 / dyu

    ! Reciprocal of Cell surface areas.
    rAc = 1.0d0 / Ac
    rAz = 1.0d0 / Az
    rAw = 1.0d0 / Aw
    rAs = 1.0d0 / As

! --

  end subroutine init_grid

! --------------------------------------------------------------------------------------------------
! Uses mask, which is loaded in init_bathy, to construct the arrays required by the ghost-point 
! treatment of boundary conditions.

  subroutine init_masks( ocean_mask, land_mask, maskavNS, maskavEW,                                &
                         umask, uland, vmask, vland, zmask, zland,                                 &
                         maskBC1, maskBC2, maskBC3, maskBC4, maskBC5, maskBC6 )

! --

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in)      :: ocean_mask

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out)     :: land_mask
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(out)   :: umask, uland
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(out)   :: maskavNS, maskBC3, maskBC4, maskBC6
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(out)   :: vmask, vland
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(out)   :: maskavEW, maskBC1, maskBC2, maskBC5
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(out) :: zmask, zland

! --

    ! Mask for no normal flow on U velocity.
    maskavNS = 0.0d0
    maskavNS(1:nx+1,1:ny) = ocean_mask(1:nx+1,1:ny)*ocean_mask(0:nx,1:ny)
    call update_overlaps_u( maskavNS )

    ! Mask for no normal flow on V velocity.
    maskavEW = 0.0d0
    maskavEW(1:nx,1:ny+1) = ocean_mask(1:nx,1:ny+1)*ocean_mask(1:nx,0:ny)
    call update_overlaps_v( maskavEW )

! --

    ! Mask for tangential flow boundary conditions on U velocity.
    maskBC3 = 0.0d0
    maskBC3(1:nx+1,1:ny) = ocean_mask(1:nx+1,1:ny)*ocean_mask(0:nx,1:ny)*        &
                                 (1.0d0-ocean_mask(1:nx+1,0:ny-1))*(1.0d0-ocean_mask(0:nx,0:ny-1))
    call update_overlaps_u( maskBC3 )

    maskBC4 = 0.0d0
    maskBC4(1:nx+1,1:ny) = ocean_mask(1:nx+1,1:ny)*ocean_mask(0:nx,1:ny)*        &
                                 (1.0d0-ocean_mask(1:nx+1,2:ny+1))*(1.0d0-ocean_mask(0:nx,2:ny+1))
    call update_overlaps_u( maskBC4 )

! --

    ! Mask for tangential flow boundary conditions on U velocity at thin walls.
    maskBC6 = 0.0d0
    maskBC6(1:nx+1,1:ny) = (1.0d0-ocean_mask(1:nx+1,1:ny))* &
                                 (1.0d0-ocean_mask(0:nx,0:ny-1))*&
                                 ocean_mask(1:nx+1,0:ny-1)*ocean_mask(1:nx+1,2:ny+1)*  &
                                 ocean_mask(0:nx,0:ny-1)*ocean_mask(0:nx,2:ny+1)
    call update_overlaps_u( maskBC6 )

! --

    ! Masks for tangential boundary conditions on V velocity.
    maskBC1 = 0.0d0
    maskBC1(1:nx,1:ny+1) = ocean_mask(1:nx,1:ny+1) *ocean_mask(1:nx,0:ny)*  &
                              (1.d0-ocean_mask(0:nx-1,1:ny+1))*(1.0d0-ocean_mask(0:nx-1,0:ny))
    call update_overlaps_v( maskBC1 )

    maskBC2 = 0.0d0
    maskBC2(1:nx,1:ny+1) = ocean_mask(1:nx,1:ny+1)*ocean_mask(1:nx,0:ny)*        &
                                (1.0d0-ocean_mask(2:nx+1,1:ny+1))*(1.0d0-ocean_mask(2:nx+1,0:ny))
    call update_overlaps_v( maskBC2 )

! --

    ! Mask for tangential flow boundary conditions on V velocity at thin walls.
    maskBC5 = 0.0d0
    maskBC5(1:nx,1:ny+1) = (1.0d0-ocean_mask(1:nx,1:ny+1))*(1.0d0-ocean_mask(1:nx,0:ny))*&
                                 ocean_mask(0:nx-1,1:ny+1)*ocean_mask(2:nx+1,1:ny+1)*  &
                                 ocean_mask(0:nx-1,0:ny)*ocean_mask(2:nx+1,0:ny)
    call update_overlaps_v( maskBC5 )

! --

    ! Mask for land at thickness points.
    land_mask = 1.0d0 - ocean_mask

! --

    ! Mask for U points at land.
    umask = 0.0d0
    umask(1:nx+1,1:ny) = ocean_mask(1:nx+1,1:ny)                                       &
                             + (1.0d0-ocean_mask(1:nx+1,1:ny))*ocean_mask(0:nx,1:ny)
    call update_overlaps_u( umask )

    ! Mask for V points at land.
    vmask = 0.0d0
    vmask(1:nx,1:ny+1) = ocean_mask(1:nx,1:ny+1)                                       &
                             + (1.0d0-ocean_mask(1:nx,1:ny+1))*ocean_mask(1:nx,0:ny)
    call update_overlaps_v( vmask )

    ! Mask for Z points at land.
    zmask = 0.0d0
    zmask(1:nx+1,1:ny+1) = 1.0d0 -                                                               &
                             (1.0d0-ocean_mask(1:nx+1,1:ny+1))*(1.0d0-ocean_mask(0:nx,1:ny+1))*&
                             (1.0d0-ocean_mask(1:nx+1,0:ny))*(1.0d0-ocean_mask(0:nx,0:ny))
    call update_overlaps_z( zmask )

! --

    ! Mask for land at U points.
    uland = 1.0d0 - umask

    ! Mask for ladn at V points.
    vland = 1.0d0 - vmask

    ! Mask for land at Z points.
    zland = 1.0d0 - zmask

! --

  end subroutine init_masks

! --------------------------------------------------------------------------------------------------
! Reads in pickup data from a series of binary files and puts it in the right place.

  subroutine init_pickup( input_dir, thickness, uvelocity, vvelocity,                              &
                          thickness_tendency, uvelocity_tendency, vvelocity_tendency )


    character(len=10) :: timestamp

    character(len=*) :: input_dir

    integer(kind=dp) :: l

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(out) :: thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(out) :: uvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(out) :: vvelocity

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly,1:ntl),   intent(out) :: thickness_tendency
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly,1:ntl), intent(out) :: uvelocity_tendency
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly,1:ntl), intent(out) :: vvelocity_tendency

    ! Write the total number of timesteps to the timestamp for the filename.
    write( timestamp, '(i0.10)' ) nt0
    write(6,*) timestamp
    write(6,*) input_dir//'Gh.'//timestamp//'.pick'

    ! Read the thickness field from binary file.
    open(unit=19,file=input_dir//'h1.'//timestamp//'.pick',status='unknown')
    read(19,*) thickness
    close(10)

    ! Read the uvelocity field from binary file.
    open(unit=20,file=input_dir//'u1.'//timestamp//'.pick',status='unknown')
    read(20,*) uvelocity
    close(20)

    ! Read the vvelocity field from binary file.
    open(unit=21,file=input_dir//'v1.'//timestamp//'.pick',status='unknown')
    read(21,*) vvelocity
    close(21)

    ! Read the thickness tendency field from binary file.
    open(unit=22,file=input_dir//'Gh.'//timestamp//'.pick',status='unknown')
    read (22,*) thickness_tendency
    close(22)

    ! Read the uvelocity tendency field from binary file.
    open(unit=23,file=input_dir//'Gu.'//timestamp//'.pick',status='unknown')
    read(23,*) uvelocity_tendency
    close(23)

    ! Read the vvelocity tendency field from binary file.
    open(unit=24,file=input_dir//'Gv.'//timestamp//'.pick',status='unknown')
    read(24,*) vvelocity_tendency
    close(24)

  end subroutine init_pickup

! --------------------------------------------------------------------------------------------------
! Initialises prognostic fields to sensible things.

  subroutine init_prognostic_fields( oldthickness, newthickness, olduvelocity, newuvelocity,       &
                                     oldvvelocity, newvvelocity,                                   &
                                     thickness_tendency, uvelocity_tendency, vvelocity_tendency )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: oldthickness, newthickness
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: olduvelocity, newuvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: oldvvelocity, newvvelocity

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: thickness_tendency
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: uvelocity_tendency
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(out) :: vvelocity_tendency

    integer(kind=dp) :: l

    ! Initialise u velocity to zero.
    olduvelocity = 0.0d0
    newuvelocity = olduvelocity

    ! Initialise v velocity to zero.
    oldvvelocity = 0.0d0
    newvvelocity = oldvvelocity

    ! Initialise with flat free-surface.
    oldthickness = 500.0d0
    newthickness = oldthickness

    ! Initialise tendency terms.
    thickness_tendency = 0.0d0
    uvelocity_tendency = 0.0d0
    vvelocity_tendency = 0.0d0
    
  end subroutine init_prognostic_fields

! --------------------------------------------------------------------------------------------------

  subroutine swap_time_levels( newthickness, newuvelocity, newvvelocity,                           &
                               oldthickness, olduvelocity, oldvvelocity,                           &
                               thickness_tendency, uvelocity_tendency, vvelocity_tendency )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(inout)   :: newthickness, oldthickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(inout) :: newuvelocity, olduvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(inout) :: newvvelocity, oldvvelocity

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly,1:ntl), intent(inout)   :: thickness_tendency
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly,1:ntl), intent(inout) :: uvelocity_tendency
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly,1:ntl), intent(inout) :: vvelocity_tendency

    integer(kind=dp) :: l

    oldthickness = newthickness
    olduvelocity = newuvelocity
    oldvvelocity = newvvelocity

    do l = ntl,2,-1
      thickness_tendency(:,:,l) = thickness_tendency(:,:,l-1)
      uvelocity_tendency(:,:,l) = uvelocity_tendency(:,:,l-1)
      vvelocity_tendency(:,:,l) = vvelocity_tendency(:,:,l-1)
    end do

  end subroutine swap_time_levels

! --------------------------------------------------------------------------------------------------

end module initial

! --------------------------------------------------------------------------------------------------
