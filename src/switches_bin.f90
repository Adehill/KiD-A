! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Switches specific to bin microphysics
!
Module switches_bin

  Implicit None

  logical ::                      &
       l_act=.false.       ! Allow activation
  
  logical ::                      &
       l_cond_evap=.true. ! Allow condensation and evaporatio (bin model)

  logical ::                      &
       l_coll_coal=.false. ! Allow collision-coalescence (bin model)

  logical ::                      &
       l_break=.False.     ! Allow collisional breakup (bin model)

  logical ::                      &
       l_sed_ult=.False.   ! switch for ultimate sedimentation

  logical ::                      &
       l_two_mom=.True.    ! switch to permit single moment bin micro
                           ! default is true for standard two moment
                           ! bin micro

end Module switches_bin
