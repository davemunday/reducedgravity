! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
! --------------------------------------------------------------------------------------------------

module numerical
  use double
  use domain, only : ntl

#include "CPP_OPTIONS.h"

! Specification of numerical parameters.

! --------------------------------------------------------------------------------------------------
! 11/05/18 - Begun line 1 re-write of reduced gravity model.
! --------------------------------------------------------------------------------------------------

#include "CPP_OPTIONS.h"

  implicit none

! --------------------------------------------------------------------------------------------------

  ! Coefficients for Adams-Bashforth 2nd/3rd order timestepping.
  real(kind=dp), dimension(1:ntl) :: a

! --

#if defined INCLUDE_OVERLAPS
  ! Flag to output overlaps and unmask ghostpoints.
  integer(kind=dp), parameter :: ioflag=1_dp
#else
  integer(kind=dp), parameter :: ioflag=0_dp
#endif

! --------------------------------------------------------------------------------------------------

end module numerical

! --------------------------------------------------------------------------------------------------
