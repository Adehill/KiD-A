! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Input variables
!
!
module input_variables

  Use typeKind
  Implicit none

  ! For ecmwf input file
  real(wp), allocatable ::  &
       time_in(:)       & ! input times
      ,nlev_in(:)       &
      ,nlevp1_in(:)     &
      ,coor_par_a(:)    &
      ,coor_par_b(:)    &
      ,p_surf_in(:)     & ! surface pressure(Pa)
      ,T_in(:,:)        & ! temperature (K)
      ,q_in(:,:)        & ! vapour (kg/kg)
      ,Tforce_in(:,:,:)     & ! Advective tendency: Temperature in K/s
      ,qforce_in(:,:,:)     & ! Advective tendency: vapour in kg/kg/s
      ,p_in(:,:)        & ! pressure in Pa on full (physics) levels
      ,p_half_in(:,:)   & ! pressure in Pa on half (velocity) levels
      ,omega_in(:,:)    & ! vertical pressure velocity Pa/s
      ,rho_in(:,:)      & ! density kg/m3
      ,z_in(:,:)        & ! height in m on full (physics) levels
      ,z_half_in(:,:)   & ! height in m on half (velocity) levels
      ,w_t(:,:,:)         & ! time varying subsidence velocity (m/s)
      ,v_t(:,:,:)         & ! time varying horizontal velocity (m/s) 	  
      ,w_t_half(:,:,:)         & 
      ,v_t_half(:,:,:)         & 
      ,wth_surf_t(:,:)    & ! time varying surface theta flux (Km/s)
      ,wqv_surf_t(:,:)      ! time varying subsidence velocity (kg/kg m/s)

  integer ::       &
       nt_in            &  !Number of input times
      ,nk_in            &  !Number of input levels
      ,nk_half_in          !Number of input half levels (for staggered
  ! grids)


  ! Forcing data timestep (for standard cases)
  real(wp) :: dt_for = 0


end module input_variables
