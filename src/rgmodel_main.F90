! --------------------------------------------------------------------------------------------------
!     D.R. MUNDAY November 2014
!     Main model routine.
! --------------------------------------------------------------------------------------------------

program rgmodel_main
  use double
  use variables
  use initial
  use diagnostics
  use forcing
  use params
  use netcdfinterface
  use thickness
  use uvelocity
  use vvelocity
  use user
  use postprocess

#include "CPP_OPTIONS.h"

! --------------------------------------------------------------------------------------------------
! 11/05/18 - Begun line 1 re-write of rg model.
! --------------------------------------------------------------------------------------------------

  implicit none

! --------------------------------------------------------------------------------------------------

  ! Timestepping variables.
  integer(kind=dp) :: k, t

  ! String to contain the number of iterations/timesteps.
  character(len=10) :: timestamp

  ! Variables for instantaneous state and time average netCDF output.
  integer(kind=4) :: state_ncid
  integer(kind=4), dimension(8) :: state_varids

! --------------------------------------------------------------------------------------------------
! Starting RG Model.

  write(6,*) ''
  write(6,*) '-- INITIALISING RGMODEL v0 --'
  
! --------------------------------------------------------------------------------------------------
! Set some non-global variables.

  ! Set timestep counter to zero.
  t = 0_dp

  ! Write the initial iteration number to the timestamp
  write(timestamp,'(I0.10)') nt0

! --------------------------------------------------------------------------------------------------

  ! Initialise the shape of the basin.
  call init_basin( tmask, tland )

  ! Initialise the grid.
  call init_grid
  
  ! Initialise the ocean/land mask and boundary condition masks.
  call init_masks( tmask, tland, maskavNS, maskavEW, umask, uland, vmask, vland, zmask, zland,     &
                   maskBC1, maskBC2, maskBC3, maskBC4, maskBC5, maskBC6 )
  
  ! Calculate the Coriolis parameter at vorticity points.
  f = init_coriolis( )
  
  ! Initialise the friction and forcing fields.
  Kv = set_lateral_diffusion( Kv0 )
  Kgm = set_lateral_diffusion( Kgm0 )
  call set_lateral_friction( AhD0, AhZ0, AhD, AhZ )
  call set_wind( taux0, tauy0, taux, tauy )

  ! Initialise prognostic fields to sensible things.
  call init_prognostic_fields( h1, h2, u1, u2, v1, v2, Gh, Gu, Gv )

  ! If starting later than the beginning of time, read in data from pickup files.
  if( nt0 > 0 )then
     call init_pickup( input_dir, h1, u1, v1, Gh, Gu, Gv )
     ! Initialise the update thickness to be the same as the current level to
     ! ensure smooth pickup.
     h2 = h1
     u2 = u1
     v2 = v1
  end if

  ! Initialise diagnostic fields.
  call init_diagnostic_fields( div, ke, psi, xi )

  ! Calculate diagnostic variables using the initial thickness/velocity fields.
  ke = tmask*calculate_kinetic_energy( u1, v1 )
  psi = calculate_streamfunction( dyg, tmask*h1, u1, maskavNS )
  xi = zmask*calculate_vorticity( rAz, dxc, dyc, u1, v1 )
  div = calculate_divergence_on_thickness( rAc, dxg, dyg, u1, v1 )

! --------------------------------------------------------------------------------------------------

  ! Initialise boundary conditions if starting from the beginning of time.
  if( nt0.eq.0 ) then
    call init_bc( u1, u2, v1, v2, maskavNS, maskavEW )
  ENDif

! --------------------------------------------------------------------------------------------------
! Initialise netcdf output.

  ! Initialise the netcdf output file for the instantaneous state.
  call create_dump_file( output_dir//'state.'//timestamp//'.nc', state_ncid, state_varids, tmask, tland )

! --------------------------------------------------------------------------------------------------
! Dump initial state.

  ! Write the initial state to netcdf file.
  write(6,*) 'Writing initial state to nc file'
  call write_dump_file( state_ncid, state_varids, int(t), int(t/nts+1), h1, psi, u1, v1, xi,       &
       ke, div, tmask, tland, umask, uland, vmask, vland, zmask, zland )
  
! --------------------------------------------------------------------------------------------------
! Write the grid information to file.

  call create_grid_file( output_dir//'grid.'//timestamp//'.nc' )

! --------------------------------------------------------------------------------------------------

  ! Choose whether to dump the initial circumpolar transport to file, depends on whether picing up 
  ! and/or writing to a new directory.
  if( (nt0==0).or.(input_dir/=output_dir) )then
    call calculate_circumpolar_transport( output_dir, t+nt0, dt, u1, h1, maskavNS )
  endif

! --------------------------------------------------------------------------------------------------

  ! Warn that we are about to start.
  write(6,*) ' '
  write(6,*) 'Entering timestepping loop'

!!                                     MAIN TIMESTEPPING LOOP                                     !!

  ! Begin loop over time.
  do t = 1,nt

! --

#if defined USE_AB2_DT
    ! Use a select case to correctly set timestepping coefficients.
    select case (nt0+t)
      case(1)
        a(1) = 1.0d0
        a(2) = 0.0d0
      case(2:)
        a(1) = 1.5d0
        a(2) = -0.5d0
      case default
        write(6,*) 'ERROR: TIMESTEPPING COEFFICIENTS INCORRECTLY SPECIFIED.'
        stop
    end select
#elif defined USE_AB3_DT
    ! Use a select case to correctly set timestepping coefficients.
    select case (nt0+t)
      case(1)
        a(1) = 1.0d0
        a(2) = 0.0d0
        a(3) = 0.0d0
      case(2)
        a(1) = 1.5d0
        a(2) = -0.5d0
        a(3) = 0.0d0
      case(3:)
        a(1) = 23.0d0/12.0d0
        a(2) = -16.0d0/12.0d0
        a(3) = 5.0d0/12.0d0
      case default
        write(6,*) 'ERROR: TIMESTEPPING COEFFICIENTS INCORRECTLY SPECIFIED.'
        stop
    end select
#endif

! --

    ! Prognostication of thickness.

    ! Evaluate the tendency due to (the divergence of) thickness advection.
    Gh(:,:,1) = calculate_thickness_tendency( Kv, Kgm, yrelax, hstar, lambda, tmask,               &
                                              h1, u1, v1, maskavNS, maskavEW )

    ! Step fluid layer thickness forwards in time.
    h2 = step_thickness_forward( a, dt, h1, Gh, hstar )

    ! Apply boundary conditions to the new thickness field.
    call apply_thickness_boundary_conditions( h1, h2, tmask, tland )
    !print*, 'h2(1,1:3) =', h2(1,1:3)

! --
    
    ! Prognostication of u velocity.

    ! Evalute the uvelocity tendency due to forcing and dissipation.
    Gu(:,:,1) = calculate_uvelocity_tendency( AhD, AhZ, f, g, taux, rho0,                          &
                                              h1, u1, v1, xi, ke, div,                             &
                                              tmask, maskavNS, maskavEW )

    ! Step fluid layer uvelocity forwards in time.
    u2 = step_uvelocity_forward( a, dt, u1, Gu )

    ! Apply boundary conditions to the new uvelocity field.
    call apply_uvelocity_boundary_conditions( slip, u2, maskavNS, maskBC3, maskBC4, maskBC6 )
    !print*, 'u2(1,1:3) =', u2(1,1:3)
    
! --

    ! Prognostication of v velocity.

    ! Evalute the uvelocity tendency due to forcing and dissipation.
    Gv(:,:,1) = calculate_vvelocity_tendency( AhD, AhZ, f, g, tauy, rho0,                          &
                                              h1, u1, v1, xi, ke, div,                             &
                                              tmask, maskavNS, maskavEW )

    ! Step fluid layer vvelocity forwards in time.
    v2 = step_vvelocity_forward( a, dt, v1, Gv )

    ! Apply boundary conditions to the new uvelocity field.
    call apply_vvelocity_boundary_conditions( slip, v2, maskavEW, maskBC1, maskBC2, maskBC5 )
    !print*, 'v2(1,1:3) =', v2(1,1:3)

! --

    ! Update diagnostic fields to the latest time level.
    ke = tmask*calculate_kinetic_energy( u2, v2 )
    psi = calculate_streamfunction( dyg, tmask*h2, u2, maskavNS )
    xi = zmask*calculate_vorticity( rAz, dxc, dyc, u2, v2 )
    div = calculate_divergence_on_thickness( rAc, dxg, dyg, u2, v2 )

! --

    ! Shuffle time levels around.
    call swap_time_levels( h2, u2, v2, h1, u1, v1, Gh, Gu, Gv )

! --

    ! If enough timesteps have elapsed, dump the state variables to netcdf file.
    if( mod(t,nts) == 0 )then
      write(6,*) ' '
      write(6,*) 'Writing state at t =', int(t), ', t/nts =', int(t/nts)
      call write_dump_file( state_ncid, state_varids, int(t), int(t/nts+1), h1, psi, u1, v1, xi,   &
                            ke, div, tmask, tland, umask, uland, vmask, vland, zmask, zland )
      endif

! --

    ! If enough timesteps have elapsed, dump pickup files for sufficient variables to allow for a smooth pickup.
    if( mod(t,ntp) == 0 )then
      call dump_pickup( t, output_dir, h1, u1, v1, Gh, Gu, Gv )
    end if

! --

    ! If enough timesteps have elapsed, do some diagnostics.
    if( mod(t,ntd) == 0 )then
      call calculate_circumpolar_transport( output_dir, t+nt0, dt, u1, h1, maskavNS )
    endif

! --

  end do
  ! End loop over time.

!!                                   END MAIN TIMESTEPPING LOOP                                   !!

  write(6,*) ''
  write(6,*) 'Exiting timestepping loop'


! --------------------------------------------------------------------------------------------------
!  Post-processing.

!  --

  ! Close the instantaneous state netcdf file.
  call close_dump_file( state_ncid )

!  --

  write(6,*) ''
  write(6,*) 'Writing final state to pickup'

  ! Dump pickup files for sufficient variables to allow for a smooth pickup.
  call dump_pickup( t-1, output_dir, h1, u1, v1, Gh, Gu, Gv )
  
! --------------------------------------------------------------------------------------------------
! Ending RG Model.

  write(6,*) ''
  write(6,*) '-- ENDING RGMODEL v0 --'

! --------------------------------------------------------------------------------------------------

end program rgmodel_main

! --------------------------------------------------------------------------------------------------
