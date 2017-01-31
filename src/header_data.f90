! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module header_data
  !
  ! Module to place header data to go into diagnostic output
  !

  Use parameters, only : max_char_len
  Implicit None

  character(100) ::    &
        modelName=     & ! driver model name - Do not change this!!
        'KiD'          &  
       ,version=       & ! driver model version number - Do not change this!!
       '!Version!'     &  
       ,revision=      & ! driver model version number - Do not change this!!
       '!Revision!'    &    
       ,title=         & ! title - Do not change this!!
       '1D Kinematic microphysics tests' & 
       ,references=    & ! references - Do not change this!!
       ''              &
       ,username=      & ! User name
       'Put Your Name Here'   &
       ,institution=   & ! Institution
       ''          &
       ,mphys_id=      & ! default microphysics scheme identifier  
       'LEM2.4'        &
       ,advection_id=  & ! default advection scheme identifier
       'LEM2.4'        &
       ,comments=      & ! user comments
       ''              

  ! arrays to hold current date and time
  integer       :: dt_values(8)
  character(5)  :: dt_zone
  character(39) :: dateString

end module header_data
