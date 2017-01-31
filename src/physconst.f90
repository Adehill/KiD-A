! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Contains values of all physical constants as parameters                          
!      CP, R, RATIO_MOL_WTS, G, RLVAP, RLSUB,                                      
! and functions of the above                                                       
!      RLVAP_ON_CP, RLSUB_ON_CP, R_ON_CP, G_ON_2                                   
Module physconst

  Use typeKind

  Implicit None


  REAL(WP), PARAMETER :: CP = 1005.                                                
  ! Specific heat of gas at constant pressure                     
  REAL(WP), PARAMETER :: R = 287.05                                                
  ! Universal gas constant                                        
  REAL(WP), PARAMETER :: RATIO_MOL_WTS = 1.608                                     
  ! Ratio of molecular weights of air and water vapour            
  REAL(WP), PARAMETER :: G = 9.81                                                  
  ! Acceleration due to gravity                                   
  REAL(WP), PARAMETER :: RLVAP = 2.501e6                                           
  ! Latent heat of vapourisation                                  
  REAL(WP), PARAMETER :: RLSUB = 2.834e6                                           
  ! Latent heat of sublimation                                    
  REAL(WP), PARAMETER :: RLFUS = RLSUB - RLVAP                                     
  ! Latent heat of fusion                                         
  REAL(WP), PARAMETER :: RLVAP_ON_CP = RLVAP/CP                                    
  REAL(WP), PARAMETER :: RLSUB_ON_CP = RLSUB/CP                                    
  REAL(WP), PARAMETER :: RLFUS_ON_CP = RLSUB_ON_CP-RLVAP_ON_CP                     
  REAL(WP), PARAMETER :: R_ON_CP = R/CP                                            
  REAL(WP), PARAMETER :: G_ON_2 = 0.5*G                                            
  ! Time data which is useful to keep constant
  REAL(WP), PARAMETER ::      &                                                     
  &   secinhr   = 3600.   &   ! second in hour                                   
  & , hrinday   = 24.     &   ! hours in day                                     
  & , secsinday = 86400.  &   ! seconds in day                                   
  & , rsecsinday = 1./secsinday      

  real(wp), parameter :: p0=100000. ! reference surface pressure
  real(wp), parameter :: pi=3.141592653589793239_wp ! pi

  real(wp), parameter :: tk0c=273.15 ! temperature in Kelvin at
                                     ! 0 degrees C  


  real(wp), parameter :: rhoW   = 1000.   ! density of water
  real(wp), parameter :: Rw     = 461.5   ! gas conbstant for water vapour
  real(wp), parameter :: THcond = 0.0243  ! thermal conductivity of air 
  real(wp), parameter :: diffWV = 2.26E-5 ! diffusivity of water vapour in air 
  real(wp), parameter :: visair = 1.44E-5 ! kinematic viscosity of air      
  real(wp), parameter :: Cwater = 4187.   ! specific heat capacity of liquid water 
  real(wp), parameter :: Cice   = 2093.   ! specific heat capacity of ice

end Module physconst
