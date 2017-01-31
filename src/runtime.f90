! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Runtime variables
!
!
Module runtime

  Use typeKind
  Implicit None

  ! For time control
  real(wp) :: time = 0.
  integer  :: time_step = 0
  integer  :: n_times   = 0  
  integer  :: n_force_times = 0 ! number of times for forcing
                                ! can differ from n_times (see test cases)
  
  ! For diagostics...
  integer  :: i_dgtime  = 0   ! diagnostic time index
  integer  :: n_dgtimes = 0   ! final number of diagnostic times

  logical  :: l_dgstep=.true. ! True if diagnostics calculated
                              ! on this timestep

end Module runtime
