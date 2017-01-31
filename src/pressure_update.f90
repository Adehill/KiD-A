module pressure_update

  ! The initial exner term is calculated in subroutine z2exner (in test_cases.f90).
  ! This is assumed to constant in the horizontal, which is OK for fixed_theta, 
  ! but not good when theta responds to advection and microphysics. This routine
  ! updates the pressure term depending on T so that density can vary depending on 
  ! theta
  Use column_variables
  Use physconst, only : Ru=>R, r_on_cp, g, p0, pi, cp
  Use parameters, only : nz, nx
  Use interpolation, only : interpolate, smooth1D

  Implicit none

contains

  subroutine exner_timestep_upd(theta_local)
    
    ! Update the exner term depending on theta
    
    real(wp), intent(in) :: theta_local(:,0:) 
    ! theta_local is theta when called from stepfields.f90 and theta + 
    ! (dtheta_adv(k,i)+dtheta_div(k,i))*dt when called before mphys
    !
    integer :: k, j

    do j=1,nx
       ! call smooth1D to be consistent with the initial z2exner
       call smooth1D(theta_local(:,j), theta_ref(:))
       exner(:,j)=(rho(:)*Ru*theta_ref(:)/p0)**(r_on_cp/(1.-r_on_cp))
    end do
    
  end subroutine exner_timestep_upd
  
end module pressure_update
