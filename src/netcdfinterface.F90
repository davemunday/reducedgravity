! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module netcdfinterface
  use double
  use domain
  use grid
  use netcdf
  use numerical

#include "CPP_OPTIONS.h"

  implicit none

! --------------------------------------------------------------------------------------------------
! This is a netcdf interface for the writing of dump files. Each dump file contains the same fields:
!    1) time;
!    2) thickness;
!    3) transport streamfunction;
!    4) u velocity;
!    5) v velocity;
!    6) vorticity;
!    7) speed squared;
!    8) divergence of velocity;
! Variables containing the x/y coordinates of all these variables are also included. These routines 
! can be used to write multiple netcdf files (as long as each file has a unique ncid).
! --------------------------------------------------------------------------------------------------

  contains

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of SWEBE to enable specification of numerical domain in terms of
!            fluid region plus overlaps.
! 24/11/14 - Introduced ioflag, which is set by macro in CPP_OPTIONS.h. This macro decides whether
!            or not to include overlaps in the output.
! --------------------------------------------------------------------------------------------------

  subroutine create_dump_file( filename, ncid, varids, ocean_mask, land_mask )

!  --

    integer(kind=4) :: status
    integer(kind=4) :: timid, xcid, ycid, xfid, yfid, xgid, ygid, xzid, yzid
    integer(kind=4) :: xc_id, yc_id, xf_id, yf_id, xg_id, yg_id, xz_id, yz_id
    integer(kind=4) :: tid, thid, uid, vid, xiid, psiid
    integer(kind=4) :: mid, lid, usqid
    integer(kind=4) :: divid
    integer(kind=dp) :: k
    real(kind=dp), dimension(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag) :: thoutput

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: ocean_mask, land_mask

    integer(kind=4), intent(out) :: ncid

    integer(kind=4), dimension(8), intent(out) :: varids

    character(len=*), intent(in) :: filename

!  --

    ! Try and open a netcdf file without overwriting it if it already exists.
    status = nf90_create( filename, nf90_share, ncid )
    if( status /= nf90_noerr ) call handle_err(status)

!  --

    ! Create the dimensions.
    status = nf90_def_dim( ncid, 'TIME', nf90_unlimited, timid )
    status = nf90_def_dim( ncid, 'XC', int(nx+2*olx*ioflag), xcid )
    status = nf90_def_dim( ncid, 'YC', int(ny+2*oly*ioflag), ycid )
    status = nf90_def_dim( ncid, 'XF', int(nx+2*olx*ioflag+1), xfid )
    status = nf90_def_dim( ncid, 'YF', int(ny+2*oly*ioflag), yfid )
    status = nf90_def_dim( ncid, 'XG', int(nx+2*olx*ioflag), xgid )
    status = nf90_def_dim( ncid, 'YG', int(ny+2*oly*ioflag+1), ygid )
    status = nf90_def_dim( ncid, 'XZ', int(nx+2*olx*ioflag+1), xzid )
    status = nf90_def_dim( ncid, 'YZ', int(ny+2*oly*ioflag+1), yzid )

!  --

    ! Create the variables.
    status = nf90_def_var( ncid, 'TIME', nf90_double, (/ timid /), tid )
    status = nf90_put_att( ncid, tid, 'long_name', 'model time' )
    status = nf90_put_att( ncid, tid, 'units', 's' )

!  --

    status = nf90_def_var( ncid, 'XC', nf90_double, (/ xcid /), xc_id )
    status = nf90_put_att( ncid, xc_id, 'long_name', 'x coordinate of thickness points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xc_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xc_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xc_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YC', nf90_double, (/ ycid /), yc_id )
    status = nf90_put_att( ncid, yc_id, 'long_name', 'y coordinate of thickness points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yc_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yc_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yc_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'H', nf90_double, (/ xcid, ycid, timid /), thid )
    status = nf90_put_att( ncid, thid, 'long_name', 'fluid layer thickness' )
    status = nf90_put_att( ncid, thid, 'units', 'm' )

    status = nf90_def_var( ncid, 'MASK', NF90_int, (/ xcid, ycid /), mid )
    status = nf90_put_att( ncid, mid, 'long_name', 'ocean mask' )

    status = nf90_def_var( ncid, 'LAND', NF90_int, (/ xcid, ycid /), lid )
    status = nf90_put_att( ncid, lid, 'long_name', 'land mask' )

    status = nf90_def_var( ncid, 'USQ', nf90_double, (/ xcid, ycid, timid /), usqid )
    status = nf90_put_att( ncid, usqid, 'long_name', 'speed squared' )
    status = nf90_put_att( ncid, usqid, 'units', 'm2/s2' )

    status = nf90_def_var( ncid, 'DIV', nf90_double, (/ xcid, ycid, timid /), divid )
    status = nf90_put_att( ncid, divid, 'long_name', 'divergence of velocity' )
    status = nf90_put_att( ncid, divid, 'units', '/s' )

!  --

    status = nf90_def_var( ncid, 'XF', nf90_double, (/ xfid /), xf_id )
    status = nf90_put_att( ncid, xf_id, 'long_name', 'x coordinate of u velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xf_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xf_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xf_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YF', nf90_double, (/ yfid /), yf_id )
    status = nf90_put_att( ncid, yf_id, 'long_name', 'y coordinate of u velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yf_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yf_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yf_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'U', nf90_double, (/ xfid, ycid, timid /), uid )
    status = nf90_put_att( ncid, uid, 'long_name', 'u velocity' )
    status = nf90_put_att( ncid, uid, 'units', 'm/s' )

!  --

    status = nf90_def_var( ncid, 'XG', nf90_double, (/ xgid /), xg_id )
    status = nf90_put_att( ncid, xg_id, 'long_name', 'y coordinate of v velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xg_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xg_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xg_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YG', nf90_double, (/ ygid /), yg_id )
    status = nf90_put_att( ncid, yg_id, 'long_name', 'y coordinate of v velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yg_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yg_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yg_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'V', nf90_double, (/ xcid, ygid, timid /), vid )
    status = nf90_put_att( ncid, vid, 'long_name', 'v velocity' )
    status = nf90_put_att( ncid, vid, 'units', 'm/s' )

!  --

    status = nf90_def_var( ncid, 'XZ', nf90_double, (/ xzid /), xz_id )
    status = nf90_put_att( ncid, xz_id, 'long_name', 'x coordinate of vorticity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xz_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xz_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xz_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YZ', nf90_double, (/ yzid /), yz_id )
    status = nf90_put_att( ncid, yz_id, 'long_name', 'y coordinate of vorticity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yz_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yz_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yz_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'XI', nf90_double, (/ xzid, yzid, timid /), xiid )
    status = nf90_put_att( ncid, xiid, 'long_name', 'vorticity' )
    status = nf90_put_att( ncid, xiid, 'units', '/s' )

    status = nf90_def_var( ncid, 'PSI', nf90_double, (/ xzid, yzid, timid /), psiid )
    status = nf90_put_att( ncid, psiid, 'long_name', 'transport streamfunction' )
    status = nf90_put_att( ncid, psiid, 'units', 'm3/s' )

!  --

    ! Create attributes with the number of grid boxes and overlap size.
    status = nf90_put_att( ncid, nf90_global, 'nx', nx )
    status = nf90_put_att( ncid, nf90_global, 'ny', ny )
    status = nf90_put_att( ncid, nf90_global, 'olx', olx )
    status = nf90_put_att( ncid, nf90_global, 'oly', oly )

!  --

    ! End define mode.
    status = NF90_ENDDEF( ncid )
    if( status.NE.nf90_noerr ) call handle_err( status )

!  --

    ! Write the coordinate variables to file.
    status = nf90_put_var( ncid, xc_id, xc(1-olx*ioflag:nx+olx*ioflag,1), (/ 1 /) )
    status = nf90_put_var( ncid, yc_id, yc(1,1-oly*ioflag:ny+oly*ioflag), (/ 1 /) )

    status = nf90_put_var( ncid, xf_id, xf(1-olx*ioflag:nx+olx*ioflag+1,1), (/ 1 /) )
    status = nf90_put_var( ncid, yf_id, yf(1,1-oly*ioflag:ny+oly*ioflag), (/ 1 /) )

    status = nf90_put_var( ncid, xg_id, xg(1-olx*ioflag:nx+olx*ioflag,1), (/ 1 /) )
    status = nf90_put_var( ncid, yg_id, yg(1,1-oly*ioflag:ny+oly*ioflag+1), (/ 1 /) )

    status = nf90_put_var( ncid, xz_id, xz(1-olx*ioflag:nx+olx*ioflag+1,1), (/ 1 /) )
    status = nf90_put_var( ncid, yz_id, yz(1,1-oly*ioflag:ny+oly*ioflag+1), (/ 1 /) )


!  --

    ! Write the land/ocean masks to file.
    status = nf90_put_var( ncid, mid, ocean_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )
    status = nf90_put_var( ncid, lid, land_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )

!  --

    ! Create array of variable ids.
    varids(1) = tid
    varids(2) = thid
    varids(3) = uid
    varids(4) = vid
    varids(5) = xiid
    varids(6) = psiid
    varids(7) = usqid
    varids(8) = divid

!  --

  end subroutine create_dump_file

! --------------------------------------------------------------------------------------------------

  subroutine create_grid_file( filename )

!  --

    integer(kind=4) :: status
    integer(kind=4) :: xcid, ycid, xfid, yfid, xgid,ygid, xzid, yzid
    integer(kind=4) :: xc_id, yc_id, xf_id, yf_id, xg_id, yg_id, xz_id, yz_id
    integer(kind=4) :: dxfid, dyfid, dxgid, dygid, dxcid, dycid, dxvid, dyuid
    integer(kind=4) :: acid, asid, awid, azid

    integer(kind=4) :: ncid

    character(len=* ), intent(in) :: filename

!  --

    ! Try and open a netcdf file without overwriting it if it already exists.
    status = nf90_create( filename, nf90_share, ncid )
    if( status /= nf90_noerr ) call handle_err(status)

!  --

    ! Create the dimensions.
    status = nf90_def_dim( ncid, 'XC', int(nx+2*olx*ioflag), xcid )
    status = nf90_def_dim( ncid, 'YC', int(ny+2*oly*ioflag), ycid )
    status = nf90_def_dim( ncid, 'XF', int(nx+2*olx*ioflag+1), xfid )
    status = nf90_def_dim( ncid, 'YF', int(ny+2*oly*ioflag), yfid )
    status = nf90_def_dim( ncid, 'XG', int(nx+2*olx*ioflag), xgid )
    status = nf90_def_dim( ncid, 'YG', int(ny+2*oly*ioflag+1), ygid )
    status = nf90_def_dim( ncid, 'XZ', int(nx+2*olx*ioflag+1), xzid )
    status = nf90_def_dim( ncid, 'YZ', int(ny+2*oly*ioflag+1), yzid )

!  --

    ! Create the variables.
    status = nf90_def_var( ncid, 'XC', nf90_double, (/ xcid /), xc_id )
    status = nf90_put_att( ncid, xc_id, 'long_name', 'x coordinate of thickness points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xc_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xc_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xc_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YC', nf90_double, (/ ycid /), yc_id )
    status = nf90_put_att( ncid, yc_id, 'long_name', 'y coordinate of thickness points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yc_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yc_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yc_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'DXF', nf90_double, (/ xcid, ycid /), dxfid )
    status = nf90_put_att( ncid, dxfid, 'long_name', 'width of thickness cell' )
    status = nf90_put_att( ncid, dxfid, 'units', 'm' )

    status = nf90_def_var( ncid, 'DYF', nf90_double, (/ xcid, ycid /), dyfid )
    status = nf90_put_att( ncid, dyfid, 'long_name', 'height of thickness cell' )
    status = nf90_put_att( ncid, dyfid, 'units', 'm' )

    status = nf90_def_var( ncid, 'AC', nf90_double, (/ xcid, ycid /), acid )
    status = nf90_put_att( ncid, acid, 'long_name', 'area of thickness cell' )
    status = nf90_put_att( ncid, acid, 'units', 'm2' )

!  --

    status = nf90_def_var( ncid, 'XF', nf90_double, (/ xfid /), xf_id )
    status = nf90_put_att( ncid, xf_id, 'long_name', 'x coordinate of u velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xf_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xf_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xf_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YF', nf90_double, (/ yfid /), yf_id )
    status = nf90_put_att( ncid, yf_id, 'long_name', 'y coordinate of u velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yf_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yf_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yf_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'DXC', nf90_double, (/ xfid, yfid /), dxcid )
    status = nf90_put_att( ncid, dxcid, 'long_name', 'width of u velocity cell' )
    status = nf90_put_att( ncid, dxcid, 'units', 'm' )

    status = nf90_def_var( ncid, 'DYG', nf90_double, (/ xfid, yfid /), dygid )
    status = nf90_put_att( ncid, dygid, 'long_name', 'height of u velocity cell' )
    status = nf90_put_att( ncid, dygid, 'units', 'm' )

    status = nf90_def_var( ncid, 'AW', nf90_double, (/ xfid, yfid /), awid )
    status = nf90_put_att( ncid, awid, 'long_name', 'area of u velocity cell' )
    status = nf90_put_att( ncid, awid, 'units', 'm2' )

!  --

    status = nf90_def_var( ncid, 'XG', nf90_double, (/ xgid /), xg_id )
    status = nf90_put_att( ncid, xg_id, 'long_name', 'y coordinate of v velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xg_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xg_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xg_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YG', nf90_double, (/ ygid /), yg_id )
    status = nf90_put_att( ncid, yg_id, 'long_name', 'y coordinate of v velocity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yg_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yg_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yg_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'DXG', nf90_double, (/ xgid, ygid /), dxgid )
    status = nf90_put_att( ncid, dxgid, 'long_name', 'width of v velocity cell' )
    status = nf90_put_att( ncid, dxgid, 'units', 'm' )

    status = nf90_def_var( ncid, 'DYC', nf90_double, (/ xgid, ygid /), dycid )
    status = nf90_put_att( ncid, dycid, 'long_name', 'height of v velocity cell' )
    status = nf90_put_att( ncid, dycid, 'units', 'm' )

    status = nf90_def_var( ncid, 'AS', nf90_double, (/ xgid, ygid /), asid )
    status = nf90_put_att( ncid, asid, 'long_name', 'area of v velocity cell' )
    status = nf90_put_att( ncid, asid, 'units', 'm2' )

!  --

    status = nf90_def_var( ncid, 'XZ', nf90_double, (/ xzid /), xz_id )
    status = nf90_put_att( ncid, xz_id, 'long_name', 'x coordinate of vorticity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, xz_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, xz_id, 'units', 'degrees longitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, xz_id, 'units', 'degrees' )
#endif

    status = nf90_def_var( ncid, 'YZ', nf90_double, (/ yzid /), yz_id )
    status = nf90_put_att( ncid, yz_id, 'long_name', 'y coordinate of vorticity points' )
#ifdef USE_CARTESIAN_GRID
    status = nf90_put_att( ncid, yz_id, 'units', 'm' )
#elif defined USE_SPHERICAL_POLAR_GRID
    status = nf90_put_att( ncid, yz_id, 'units', 'degrees latitude' )
#elif defined USE_CYLINDRICAL_GRID
    status = nf90_put_att( ncid, yz_id, 'units', 'm' )
#endif

    status = nf90_def_var( ncid, 'DXV', nf90_double, (/ xzid, yzid /), dxvid )
    status = nf90_put_att( ncid, dxvid, 'long_name', 'width of vorticity cell' )
    status = nf90_put_att( ncid, dxvid, 'units', 'm' )

    status = nf90_def_var( ncid, 'DYU', nf90_double, (/ xzid, yzid /), dyuid )
    status = nf90_put_att( ncid, dyuid, 'long_name', 'height of vorticity cell' )
    status = nf90_put_att( ncid, dyuid, 'units', 'm' )

    status = nf90_def_var( ncid, 'AZ', nf90_double, (/ xzid, yzid /), azid )
    status = nf90_put_att( ncid, azid, 'long_name', 'area of vorticity cell' )
    status = nf90_put_att( ncid, azid, 'units', 'm2' )

!  --

    ! Create attributes with the number of grid boxes and overlap size.
    status = nf90_put_att( ncid, nf90_global, 'nx', nx )
    status = nf90_put_att( ncid, nf90_global, 'ny', ny )
    status = nf90_put_att( ncid, nf90_global, 'olx', olx )
    status = nf90_put_att( ncid, nf90_global, 'oly', oly )

!  --

    ! End define mode.
    status = NF90_ENDDEF( ncid )
    if( status.NE.nf90_noerr ) call handle_err( status )

!  --

    ! Write the coordinate variables to file.
    status = nf90_put_var( ncid, xc_id, xc(1-olx*ioflag:nx+olx*ioflag,1), (/ 1 /) )
    status = nf90_put_var( ncid, yc_id, yc(1,1-oly*ioflag:ny+oly*ioflag), (/ 1 /) )

    status = nf90_put_var( ncid, xf_id, xf(1-olx*ioflag:nx+olx*ioflag+1,1), (/ 1 /) )
    status = nf90_put_var( ncid, yf_id, yf(1,1-oly*ioflag:ny+oly*ioflag), (/ 1 /) )

    status = nf90_put_var( ncid, xg_id, xg(1-olx*ioflag:nx+olx*ioflag,1), (/ 1 /) )
    status = nf90_put_var( ncid, yg_id, yg(1,1-oly*ioflag:ny+oly*ioflag+1), (/ 1 /) )

    status = nf90_put_var( ncid, xz_id, xz(1-olx*ioflag:nx+olx*ioflag+1,1), (/ 1 /) )
    status = nf90_put_var( ncid, yz_id, yz(1,1-oly*ioflag:ny+oly*ioflag+1), (/ 1 /) )

!  --

    ! Write the grid spacings to file
    status = nf90_put_var( ncid, dxfid, dxf(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )
    status = nf90_put_var( ncid, dyfid, dyf(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )
    status = nf90_put_var( ncid, dxcid, dxc(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )
    status = nf90_put_var( ncid, dygid, dyg(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )
    status = nf90_put_var( ncid, dxgid, dxg(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1), (/ 1, 1 /) )
    status = nf90_put_var( ncid, dycid, dyc(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1), (/ 1, 1 /) )
    status = nf90_put_var( ncid, dxvid, dxv(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1), (/ 1, 1 /) )
    status = nf90_put_var( ncid, dyuid, dyu(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1), (/ 1, 1 /) )
    status = nf90_put_var( ncid, acid, Ac(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )
    status = nf90_put_var( ncid, awid, Aw(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag), (/ 1, 1 /) )
    status = nf90_put_var( ncid, asid, As(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1), (/ 1, 1 /) )
    status = nf90_put_var( ncid, azid, Az(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1), (/ 1, 1 /) )

!  --

    ! Close the grid file.
    call close_dump_file( ncid )

!  --

  end subroutine create_grid_file

! --------------------------------------------------------------------------------------------------

  subroutine write_dump_file( ncid, varids, t, ndump, thickness, streamfunction, uvelocity,        &
                              vvelocity, vorticity, speed_squared, divergence,                     &
                              ocean_mask, land_mask, umask, uland, vmask, vland, zmask, zland )

!  --

    integer(kind=4) :: status
    integer(kind=dp) :: k

    real(kind=dp), dimension(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)     :: thoutput
    real(kind=dp), dimension(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag)   :: uoutput
    real(kind=dp), dimension(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1)   :: voutput
    real(kind=dp), dimension(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1) :: zoutput

    integer(kind=4), intent(in) :: ncid, t, ndump
    integer(kind=4), dimension(8), intent(in) :: varids
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in)     :: ocean_mask, land_mask
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in)   :: umask, uland
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in)   :: vmask, vland
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: zmask, zland

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),     intent(in) :: thickness, divergence, speed_squared
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly),   intent(in) :: uvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly),   intent(in) :: vvelocity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(in) :: vorticity, streamfunction

! --

    ! Increment the time axis.
    status = nf90_put_var( ncid, varids(1), t*dt, start=(/ ndump /) )

! --

    ! Mask thickness.
    thoutput(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag) =                            &
         thickness(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)!                          &
!         *ocean_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)                        &
!         - 000000.0d0*dble(1-ioflag)*land_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)

    ! Dump thickness to file.
    status = nf90_put_var( ncid, varids(2), thoutput, (/ 1, 1, ndump /) )

! --

    ! Mask the u velocity.
    uoutput(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag) =                             &
         uvelocity(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag)!                        &
!         *umask(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag)                           &
!         - 000000.0d0*dble(1-ioflag)*uland(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag)

    ! Dump the u velocity to file.
    status = nf90_put_var( ncid, varids(3), uoutput, (/ 1, 1, ndump /) )

! --

    ! Mask the v velocity.
    voutput(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1) =                             &
         vvelocity(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1)!                        &
!         *vmask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1)                           &
!         - 000000.0d0*dble(1-ioflag)*vland(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag+1)

    ! Dump the v velocity to file.
    status = nf90_put_var( ncid, varids(4), voutput, (/ 1, 1, ndump /) )

! --

    ! Mask the vorticity.
    zoutput(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1) =                           &
         vorticity(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1)!                      &
!         *zmask(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1)                         &
!         - 000000.0d0*dble(1-ioflag)*zland(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1)

    ! Dump vorticity to filed
    status = nf90_put_var( ncid, varids(5), zoutput, (/ 1, 1, ndump /) )

! --

    ! Mask streamfunction.
    zoutput(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1) =                           &
         streamfunction(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1)!                 &
!         *zmask(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1)                         &
!         - 000000.0d0*dble(1-ioflag)*zland(1-olx*ioflag:nx+olx*ioflag+1,1-oly*ioflag:ny+oly*ioflag+1)

    ! Dump streamfunction to file.
    status = nf90_put_var( ncid, varids(6), zoutput, (/ 1, 1, ndump /) )

! --

    ! Mask speed squared.
    thoutput(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag) =                              &
         speed_squared(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)!                      &
!         *ocean_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)                        &
!         - 000000.0d0*dble(1-ioflag)*land_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)

    ! Dump speed squared to file.
    status = nf90_put_var( ncid, varids(7), thoutput, (/ 1, 1, ndump /) )

! --

    ! Mask divergence of velocity.
    thoutput(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag) =                              &
         divergence(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)!                         &
!         *ocean_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)                        &
!         - 000000.0d0*dble(1-ioflag)*land_mask(1-olx*ioflag:nx+olx*ioflag,1-oly*ioflag:ny+oly*ioflag)

    ! Dump divergence to file.
    status = nf90_put_var( ncid, varids(8), thoutput, (/ 1, 1, ndump /) )

! --

  end subroutine write_dump_file

! --------------------------------------------------------------------------------------------------

  subroutine close_dump_file( ncid )

    integer(kind=4) :: status

    integer(kind=4), intent(in) :: ncid

    ! Close the netcdf file.
    status = NF90_close( ncid )

    ! Stop if there is an error.
    if( status /= nf90_noerr ) call handle_err(status)

  end subroutine close_dump_file

! --------------------------------------------------------------------------------------------------

  subroutine handle_err( status )

    integer, intent(in) :: status

    if( status.NE.nf90_noerr )then
      PRint*, trim( nf90_strerror( status ) )
      stop "Stopped"
    ENDif

  end subroutine handle_err

! --------------------------------------------------------------------------------------------------

end module netcdfinterface

! --------------------------------------------------------------------------------------------------
