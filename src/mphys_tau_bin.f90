! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!
module mphys_tau_bin

contains

  subroutine mphys_tau_bin_interface

    write(*,'(8(/, A))') & 
         '===========================================================' &
         ,'This is a dummy stub for the Tel-Aviv University bin ' &
         , 'microphysics code. This code is not part of the current'  &
         ,'distribution.  To obtain the code, either email '& 
         ,'ben.ship'//&
         'way@metoffice.gov.uk or visit the downloads' &
         ,'page at http://www.convection.info/microphysics' &
         ,'The program will now exit...' &
         ,'==========================================================' 
    stop

  end subroutine mphys_tau_bin_interface

end module
