module aerosols
  ! This module replaces the more comprehensive version from the Shipway 
  ! microphysics scheme - to be used in setting up aerosol distributions

  ! For setting aerosol
  Use class_species, only: species_allocate
  Use column_variables, only: aerosol_init
  Use parameters, only: num_aero_moments, num_aero_bins, nz, naerosol
  Use switches, only: l_ainit
  Use physconst, only : pi

  Implicit None

  integer, parameter ::                                     &
       FourByteReal = selected_real_kind(P =  6, R =  37)   &
       ,EightByteReal = selected_real_kind(P = 13, R =  307)

  integer, parameter ::                 &
     !       wp=FourByteReal                  & ! real working precision 
     wp=EightByteReal                    ! real working precision 


contains

  function moment_logn(N, rm, sigma, p)
    !
    ! Calculate the moments of a lognormal distribution
    !
    ! int_0^infinty r^p*n(r)dr
    !
    ! where 
    !
    !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
    !            *exp(-.5*(log(r/rm)/log(sigma))^2)
    !
    ! rm = median radius
    ! log(rm) = mean of log(n(r)/N)
    ! log(sigma) = s.d. of log(n(r)/N)
    !
    
    real(wp) :: N, rm, sigma
    real(wp) :: p ! calculate pth moment
    real(wp) :: moment_logn

    moment_logn=N*rm**p*exp(.5*p*p*log(sigma)**2)

  end function moment_logn


  subroutine set_aerosol(Ninit, indices, lainits_in, Nds_in, sigmas_in, rds_in, densitys_in, fscale)

    integer, intent(in) :: Ninit
    integer, intent(in) :: indices(Ninit)
    logical, intent(in) :: lainits_in(Ninit)
    real(wp), intent(in) :: Nds_in(Ninit)
    real(wp), intent(in) :: sigmas_in(Ninit)
    real(wp), intent(in) :: rds_in(Ninit)
    real(wp), intent(in) :: densitys_in(Ninit)
    real(wp), intent(in) :: fscale(nz)

    real(wp) :: nd, sigma, rd, density
    integer :: in, ih, k

    do in=1, Ninit
      ih=indices(in)
      if (lainits_in(ih).and.ih<=naerosol)then
        l_ainit(ih)=.True.
        Nd=Nds_in(ih)
        sigma=sigmas_in(ih)
        rd=rds_in(ih)
        density=densitys_in(ih)
        do k=1,nz
          call species_allocate(aerosol_init(k,ih) &
             , num_aero_moments(ih), num_aero_bins(ih), ih, 0.)
          if (num_aero_moments(ih) > 0) &
             aerosol_init(k,ih)%moments(1,1)=Nd
          if (num_aero_moments(ih) > 1) &
             aerosol_init(k,ih)%moments(1,2)=(4.*pi*density/3.)*moment_logn(Nd, rd, sigma, 3._wp)*fscale(k)
        end do
      end if
    end do

  end subroutine set_aerosol


end module aerosols
