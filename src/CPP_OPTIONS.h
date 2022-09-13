!  ----------------------------------------------------------------------------
!     D.R. MUNDAY May 2018
!  ----------------------------------------------------------------------------

#ifndef CPP_OPTIONS_H
#define CPP_OPTIONS_H

! --------------------------------------------------------------------------------------------------
! Numerical solution options.

! Select a timestepping routine - AB3 (Adams-Bashforth 3rd order) is HIGHLY RECOMMENDED.
#undef USE_AB2_DT
#define USE_AB3_DT

!  ----------------------------------------------------------------------------
!  Equation format options.

! Activate flux limitation of thickness advection
#undef FLUX_LIMIT_ADVECTION

! Activate time averaging code.
#undef INCLUDE_TAVE

! Choose whether to include overlaps and unmasked ghostpoints in output.
#undef INCLUDE_OVERLAPS

!  ----------------------------------------------------------------------------
!  Numerical method choies.
!  ----------------------------------------------------------------------------

#endif /* CPP_OPTIONS_H */

!  ----------------------------------------------------------------------------
