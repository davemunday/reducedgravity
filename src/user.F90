!  ----------------------------------------------------------------------------
!     D.R. MUNDAY January 2014
!  ----------------------------------------------------------------------------

module user
  use double
  use params
  use grid

#include "CPP_OPTIONS.h"

  implicit none

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

!  ----------------------------------------------------------------------------

  contains

!  ----------------------------------------------------------------------------

  subroutine calculate_circumpolar_transport( output_directory, no_iterations, timestep,           &
                                              uvelocity, thickness, maskavNS )

    integer(kind=dp) :: j, k

    character(len=*), intent(in) :: output_directory

    integer(kind=dp), intent(in) :: no_iterations
    real(kind=dp), intent(in) :: timestep

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly), intent(in) :: thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: uvelocity
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: maskavNS

    real(kind=dp) :: circumpolar_transport

    ! Initialise the circumpolar_transport.
    circumpolar_transport = 0.0d0

    ! Calculate the transport in the middle of the domain.
    do j = 1,ny
       circumpolar_transport = circumpolar_transport + 0.5d0*( thickness(138,j)+thickness(137,j) ) &
                               * dyg(137,j)*maskavNS(137,j)*uvelocity(137,j)
    end do

    open(29, file=output_directory//'transport.dat',access='append', status='unknown' )
    write(29,*) timestep*no_iterations, circumpolar_transport/1.0d6
    close(29)

  end subroutine calculate_circumpolar_transport

!  ----------------------------------------------------------------------------

end module user

! -----------------------------------------------------------------------------
