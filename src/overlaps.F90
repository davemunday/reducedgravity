! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
!     Subroutines that deal with the overlap regions, largely to update them.
! --------------------------------------------------------------------------------------------------

module overlaps
  use double
  use domain, only : nx, ny, olx, oly

#include "CPP_OPTIONS.h"

  implicit none


! --------------------------------------------------------------------------------------------------
! 11/05/18 - Begun line 1 re-write of reduced gravity model.
! --------------------------------------------------------------------------------------------------
! Interface to overload the subroutines for use on integers and real numbers.

  interface update_overlaps_t
    module procedure update_overlaps_tfield, update_overlaps_tmask
  end interface update_overlaps_t

  interface update_overlaps_u
    module procedure update_overlaps_ufield, update_overlaps_umask
  end interface update_overlaps_u

  interface update_overlaps_v
    module procedure update_overlaps_vfield, update_overlaps_vmask
  end interface update_overlaps_v

  interface update_overlaps_z
    module procedure update_overlaps_zfield, update_overlaps_zmask
  end interface update_overlaps_z

contains

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for a real field on T points.

  subroutine update_overlaps_tfield( tfield )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(inout) :: tfield
    
    ! Reentrance in X.
    tfield(1-olx:0,:) = tfield(nx-olx+1:nx,:)
    tfield(nx+1:nx+olx,:) = tfield(1:olx,:)

    ! Reentrance in Y.
    tfield(:,1-oly:0) = tfield(:,ny-oly+1:ny)
    tfield(:,ny+1:ny+oly) = tfield(:,1:oly)

  end subroutine update_overlaps_tfield

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for an integer field on T points.

  subroutine update_overlaps_tmask( tfield )

    integer(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(inout) :: tfield
    
    ! Reentrance in X.
    tfield(1-olx:0,:) = tfield(nx-olx+1:nx,:)
    tfield(nx+1:nx+olx,:) = tfield(1:olx,:)

    ! Reentrance in Y.
    tfield(:,1-oly:0) = tfield(:,ny-oly+1:ny)
    tfield(:,ny+1:ny+oly) = tfield(:,1:oly)

  end subroutine update_overlaps_tmask

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for a real field on U points.

  subroutine update_overlaps_ufield( ufield )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(inout) :: ufield

    ! Reentrance in X.
    ufield(nx+1:nx+olx+1,:) = ufield(1:olx+1,:)
    ufield(1-olx:1,:) = ufield(nx-olx+1:nx+1,:)

    ! Reentrance in Y.
    ufield(:,ny+1:ny+oly) = ufield(:,1:oly)
    ufield(:,1-oly:0) = ufield(:,ny-oly+1:ny)

  end subroutine update_overlaps_ufield

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for an integer field on U points.

  subroutine update_overlaps_umask( ufield )

    integer(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(inout) :: ufield

    ! Reentrance in X.
    ufield(nx+1:nx+olx+1,:) = ufield(1:olx+1,:)
    ufield(1-olx:1,:) = ufield(nx-olx+1:nx+1,:)

    ! Reentrance in Y.
    ufield(:,ny+1:ny+oly) = ufield(:,1:oly)
    ufield(:,1-oly:0) = ufield(:,ny-oly+1:ny)

  end subroutine update_overlaps_umask

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for a real field on V points.

  subroutine update_overlaps_vfield( vfield )

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(inout) :: vfield

    ! Reentrance in X.
    vfield(1-olx:0,:) = vfield(nx-olx+1:nx,:)
    vfield(nx+1:nx+olx,:) = vfield(1:olx,:)

    ! Reentrance in Y.
    vfield(:,ny+1:ny+oly+1) = vfield(:,1:oly+1)
    vfield(:,1-oly:1) = vfield(:,ny-oly+1:ny+1)

  end subroutine update_overlaps_vfield

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for an integer field on V points.

  subroutine update_overlaps_vmask( vfield )

    integer(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(inout) :: vfield

    ! Reentrance in X.
    vfield(1-olx:0,:) = vfield(nx-olx+1:nx,:)
    vfield(nx+1:nx+olx,:) = vfield(1:olx,:)

    ! Reentrance in Y.
    vfield(:,ny+1:ny+oly+1) = vfield(:,1:oly+1)
    vfield(:,1-oly:1) = vfield(:,ny-oly+1:ny+1)

  end subroutine update_overlaps_vmask

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for a real field on Z points.

  subroutine update_overlaps_zfield( zfield )

    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(inout) :: zfield

    ! Reentrance in X.
    zfield(1-olx:1,:) = zfield(nx-olx+1:nx+1,:)
    zfield(nx+1:nx+olx+1,:) = zfield(1:olx+1,:)

    ! Reentrance in Y.
    zfield(:,1-oly:1) = zfield(:,ny-oly+1:ny+1)
    zfield(:,ny+1:ny+oly+1) = zfield(:,1:oly+1)

  end subroutine update_overlaps_zfield

! --------------------------------------------------------------------------------------------------
! This subroutine updates the overlap region for an integer field on Z points.

  subroutine update_overlaps_zmask( zfield )

    integer(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+1+oly), intent(inout) :: zfield

    ! Reentrance in X.
    zfield(1-olx:1,:) = zfield(nx-olx+1:nx+1,:)
    zfield(nx+1:nx+olx+1,:) = zfield(1:olx+1,:)

    ! Reentrance in Y.
    zfield(:,1-oly:1) = zfield(:,ny-oly+1:ny+1)
    zfield(:,ny+1:ny+oly+1) = zfield(:,1:oly+1)

  end subroutine update_overlaps_zmask

! --------------------------------------------------------------------------------------------------

end module overlaps

! --------------------------------------------------------------------------------------------------
