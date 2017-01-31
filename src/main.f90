! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Driver for 1D Kinematic Driver (KiD) model 
!
! Author: Ben Shipway
!
! For version details see header_data.f90
!

Module main
 
  Use typeKind
  Use parameters, only : dt, dg_dt, nx
  Use namelists, only : read_namelist
  Use runtime, only : time, time_step, n_times
  Use switches
  Use set_profiles,  only : read_profiles
  Use interpolation, only : interpolate_input, interpolate_forcing
  Use diagnostics,   only : save_diagnostics_1d, save_diagnostics_2d, &
       write_diagnostics, query_dgstep
  Use derived_fields, only : calc_derived_fields
  Use advection_interface, only : advect_column
  Use mphys_interface, only : mphys_column
  Use stepfields, only : step_column
  Use divergence, only : diverge_column

  Implicit none

contains

  subroutine main_loop
    
    integer :: itime      ! loop counter for time
    !
    ! Start by reading in namelists
    !

    if (l_namelists) call read_namelist

    ! Set up the initial fields and forcing
    if (l_input_file)then
       call read_profiles(input_file)
    else
       call read_profiles(icase)
    end if

    call interpolate_input(ifiletype)

    call interpolate_forcing

    call calc_derived_fields

    ! Do we want to do diagnostics on this timestep?
    call query_dgstep

    if ( nx == 1 ) then 
       call save_diagnostics_1d
    else 
       call save_diagnostics_2d
    endif

    do itime=1,n_times

       time=time+dt
       time_step=time_step+1
  
       ! Do we want to do diagnostics on this timestep?
       call query_dgstep

       call interpolate_forcing

       call calc_derived_fields

       if (l_advect)then
          call advect_column(scheme_id=0)
       end if

       if (l_diverge)then
          call diverge_column
       end if

       if (l_mphys)then
          call mphys_column(scheme_id=imphys)
       end if

       call step_column
       
       if ( nx == 1 ) then
          call save_diagnostics_1d
       else
          call save_diagnostics_2d
       endif

    end do

    if (l_write_dgs) call write_diagnostics

  end subroutine main_loop

End Module main
