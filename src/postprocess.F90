! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module postprocess
  use double
  use grid
  use domain

#include "CPP_OPTIONS.h"

! This module contains post-processing routines for the when the model has finished running.

! --------------------------------------------------------------------------------------------------
! 13/11/14 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  contains

! --------------------------------------------------------------------------------------------------

  subroutine dump_pickup( iteration_no, output_dir, thickness, uvelocity, vvelocity,               &
                          thickness_tendency, uvelocity_tendency, vvelocity_tendency )

    character(len=10) :: timestamp

    character(len=*), intent(in) :: output_dir

    integer(kind=dp) :: l
    integer(kind=dp), intent(in) :: iteration_no

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly),   intent(in) :: thickness
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly), intent(in) :: uvelocity
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly), intent(in) :: vvelocity

    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+oly,1:ntl),   intent(in) :: thickness_tendency
    real(kind=dp), dimension(1-olx:nx+1+olx,1-oly:ny+oly,1:ntl), intent(in) :: uvelocity_tendency
    real(kind=dp), dimension(1-olx:nx+olx,1-oly:ny+1+oly,1:ntl), intent(in) :: vvelocity_tendency

    ! Write the total number of timesteps to the timestamp for the filename.
    write( timestamp, '(i0.10)' ) nt0 + iteration_no

    ! Write the thickness field to binary file.
    open(unit=19,file=output_dir//'h1.'//timestamp//'.pick',status='unknown')
    write(19,*) thickness
    close(19)

    ! Write the uvelocity field to binary file.
    open(unit=20,file=output_dir//'u1.'//timestamp//'.pick',status='unknown')
    write(20,*) uvelocity
    close(20)

    ! Write the vvelocity field to binary file.
    open(unit=21,file=output_dir//'v1.'//timestamp//'.pick',status='unknown')
    write(21,*) vvelocity
    close(21)

    ! Write the thickness tendency field to binary file.
    open(unit=22,file=output_dir//'Gh.'//timestamp//'.pick',status='unknown')
    write(22,*) thickness_tendency
    close(22)

    ! Write the uvelocity tendency field to binary file.
    open(unit=23,file=output_dir//'Gu.'//timestamp//'.pick',status='unknown')
    write(23,*) uvelocity_tendency
    close(23)

    ! Write the vvelocity tendency field to binary file.
    open(unit=24,file=output_dir//'Gv.'//timestamp//'.pick',status='unknown')
    write(24,*) vvelocity_tendency
    close(24)

  end subroutine dump_pickup

! --------------------------------------------------------------------------------------------------

end module postprocess

! --------------------------------------------------------------------------------------------------
