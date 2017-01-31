! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Some parameters and arrays needed by different schemes
! Note that some of these are no longer parameters, but can be 
! passed in through a namelist (see namelists.f90)
!
!

#ifndef DEF_NX
#define DEF_NX 1
#endif

#ifndef DEF_NZ
#define DEF_NZ 120
#endif

module parameters
  

  Use typeKind
  Implicit None

  ! control parameters
  real               :: dt=1.        ! Time step 
  integer, parameter :: nz=DEF_NZ    ! Number of height levels
  integer, parameter :: nx=DEF_NX    ! Number of horizontal gridpoints

  ! microphysics specific parameters
  integer, parameter :: max_nmoments=3
  integer, parameter :: max_nbins=34
  integer, parameter :: nspecies=5
#if SHIPWAY_MICRO == 1
  integer, parameter :: naerosol=9 ! number of aerosol species
#else
  integer, parameter :: naerosol=3 ! number of aerosol species
#endif

  ! diagnostic parameters
  integer, parameter ::     &
       max_dgs       = 300  & ! Maximum number of dg variables in
                              ! a diagnostic array
      ,max_char_len  = 200    ! Maximum variable name length

  real(wp)           :: dg_dt  = 60.
                        ! time interval between diagnostic calcs
                        ! should be an integer multiple of dt

  integer            :: diaglevel = 5 ! switch for level of diagnostic output 
  
  real(wp)           :: dgstart = 0.0 ! start time for diagnostics 

 ! other parameters
  real(wp),parameter :: unset_real=-999
  integer,parameter  :: unset_integer=-999

  ! Microphysics parameters
  integer :: num_aero_moments(naerosol) = 1
  integer :: num_aero_bins(naerosol)= 1
  real(wp) :: aero_mom_init(max_nmoments)=(/  &   ! Obsolete
         100. &  ! Initialise number              ! Obsolete
      ,  0. &  ! Initialise mass                  ! Obsolete
      ,  0. &  ! Initialise 3rd moment?           ! Obsolete
        /)

  real(wp) :: aero_N_init(naerosol)=0.0
  real(wp) :: aero_rd_init(naerosol)=0.0
  real(wp) :: aero_sig_init(naerosol)=0.0

  integer :: num_h_moments(nspecies)= (/  &
        1 & ! cloud
       ,1 & ! rain
       ,2 & ! ice
       ,1 & ! snow
       ,1 & ! graupel
       /)

  ! If using bulk microphysics, simply set these to 1.
  integer :: num_h_bins(nspecies)= 1

  real(wp) :: mom_init(max_nmoments)= (/  &
        0. &  ! Initialise mass
       ,1. &  ! Initialise number
       ,0. &  ! Initialise volume
       /)

  character(10) :: h_names(nspecies)= &
       (/ 'cloud     '   &
       ,  'rain      '   &
       ,  'ice       '   &
       ,  'snow      '   &
       ,  'graupel   '  /)

  character(10) :: mom_names(max_nmoments)= &
       (/ 'mass      '   &
       ,  'number    '   &
       ,  'volume    '  /)

  character(10) :: mom_units(max_nmoments)= &
       (/ 'kg/kg     '   &
       ,  '/kg       '   &
       ,  'm3/kg     '  /)

  character(10) :: aero_names(naerosol)= & 
     (/ 'aitken    '   & 
     ,  'accum     '   & 
     ,  'coarse    '   & 
#if SHIPWAY_MICRO == 1
     ,  'active    '   & 
     ,  'active_r  '   & 
     ,  'dust      '   & 
     ,  'dust_ice  '   & 
     ,  'dust_cloud'   & 
     ,  'sol_ice   '   & 
#endif
     /) 
  
  character(10) :: aero_mom_names(max_nmoments)= & 
     (/ 'number    '   & 
     ,  'mass      '   & 
     ,  'unknown   '  /) 
  
  character(10) :: aero_mom_units(max_nmoments)= & 
     (/ 'kg/kg     '   & 
     ,  '/kg       '   & 
     ,  'unknown   '  /) 

   ! if using bin microphysics set variable below to 
   ! determine the last bin of cloud
   ! if not using bin, this value is ignored
   !integer :: split_bins = 15

end module parameters
