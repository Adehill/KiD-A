! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Interpolate routines
!
!
module interpolation

  Use typeKind

  Implicit none

contains

  subroutine interpolate_input(input_type)

    Use parameters, only : dt, nspecies, num_h_moments, &
         num_h_bins, mom_init, num_aero_moments, &
         num_aero_bins,aero_mom_init, nz, nx
    Use runtime, only : n_times, n_force_times
    Use column_variables
    Use physconst, only : Ru=>R, g, r_on_cp, p0
    Use switches
    Use input_variables
    Use class_species, only: species_allocate

    integer, optional :: input_type   ! specifies format
    ! for input fields

    !local variables
    integer :: itype, k, j, ih,imom
    
    real :: tend_init=0.  ! initial value for tendency terms

! Initialize aerosols and their tendencies
    do ih=1,naerosol
       do j=0,nx+1
          do k=1,nz
             call species_allocate(aerosol(k,j,ih) &
                  , num_aero_moments(ih), num_aero_bins(ih), ih)
             call species_allocate(daerosol_mphys(k,j,ih) &
                  , num_aero_moments(ih), num_aero_bins(ih), ih &
                  , tend_init)
             call species_allocate(daerosol_adv(k,j,ih) &
                  , num_aero_moments(ih), num_aero_bins(ih), ih &
               , tend_init)
             call species_allocate(daerosol_div(k,j,ih) &
                  , num_aero_moments(ih), num_aero_bins(ih), ih &
                  , tend_init)
             call species_allocate(daerosol_force(k,j,ih) &
                  , num_aero_moments(ih), num_aero_bins(ih), ih &
                  , tend_init)
             do imom=1,num_aero_moments(ih)
                aerosol(k,j,ih)%moments(:,imom)=aero_mom_init(imom) 
             end do
          end do
       end do
    end do

    ! Initialize hydrometeors and their tendencies
    do ih=1,nspecies
       do j=0,nx+1
          do k=1,nz
             call species_allocate(hydrometeors(k,j,ih) &
                  , num_h_moments(ih), num_h_bins(ih), ih)
             call species_allocate(dhydrometeors_mphys(k,j,ih) &
                  , num_h_moments(ih), num_h_bins(ih), ih &
                  , tend_init)
             call species_allocate(dhydrometeors_adv(k,j,ih) &
                  , num_h_moments(ih), num_h_bins(ih), ih &
                  , tend_init)
             call species_allocate(dhydrometeors_div(k,j,ih) &
                  , num_h_moments(ih), num_h_bins(ih), ih &
                  , tend_init)
             call species_allocate(dhydrometeors_force(k,j,ih) &
                  , num_h_moments(ih), num_h_bins(ih), ih &
                  , tend_init)
             do imom=1,num_h_moments(ih) 
                hydrometeors(k,j,ih)%moments(:,imom)=mom_init(imom)
             end do
          end do
       end do
    end do


    if (present(input_type)) then
       itype=input_type
    else
       itype=iukmo_lem
    end if


    select case (itype)
    case(itest_case)
       ! Initialize some species
       do ih=1,nspecies
          if (l_hinit(ih))then
             do j=0,nx+1 
                do k=1,nz
                   do imom=1,num_h_moments(ih)
                      hydrometeors(k,j,ih)%moments(:,imom) = &
                           hydrometeors_init(k,ih)%moments(:,imom)
                   end do
                end do
             end do
          end if
        end do
        ! Initialize some aerosol species 
        do ih=1,naerosol 
          if (l_ainit(ih))then 
            do j=0,nx+1 
              do k=1,nz 
                do imom=1,num_aero_moments(ih) 
                  aerosol(k,j,ih)%moments(:,imom) = & 
                     aerosol_init(k,ih)%moments(:,imom) 
                end do
              end do
            end do
          end if
        end do
      end select


  end subroutine interpolate_input

  
  subroutine interpolate_forcing

    Use parameters, only : nz, nx, dt
    Use column_variables, only : w, v, w_half, v_half,  z, z_half, x, x_half, & 
       Tforce, qforce, wth_surf, wqv_surf &
       ,Trelax, qv, theta, exner, thinit, qvinit
    Use input_variables, only : w_t, time_in, Tforce_in, qforce_in &
         , wth_surf_t,  wqv_surf_t, v_t, w_t_half, v_t_half
    Use runtime, only : time, n_times, n_force_times
    Use switches, only: iforce_method

    integer :: itime, loc_ntimes, k

    itime=interval(time_in, time)

    if ( n_force_times > 0 .and. n_force_times < n_times ) then 
       loc_ntimes = n_force_times
    else
       loc_ntimes = n_times
    endif

    if (itime>=1.and. itime < loc_ntimes ) then
       w(1:nz,0:nx+1)= w_t(1:nz,0:nx+1,itime) +               &
            ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  &
            *(w_t(1:nz,0:nx+1,itime+1)-w_t(1:nz,0:nx+1,itime))
       w_half(1:nz,0:nx+1)= w_t_half(1:nz,0:nx+1,itime) +               & 
          ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  & 
          *(w_t_half(1:nz,0:nx+1,itime+1)-w_t_half(1:nz,0:nx+1,itime))  

       v(1:nz,0:nx+1)= v_t(1:nz,0:nx+1,itime) +               &
            ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  &
            *(v_t(1:nz,0:nx+1,itime+1)-v_t(1:nz,0:nx+1,itime))
       if (nx > 1) then 
          v_half(1:nz,0:nx+1)= v_t_half(1:nz,0:nx+1,itime) +               & 
             ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  & 
             *(v_t_half(1:nz,0:nx+1,itime+1)-v_t_half(1:nz,0:nx+1,itime))          
        else 
          v_half = v 
        endif

        select case (iforce_method)
        case(0)
          Tforce(1:nz,0:nx+1)=Tforce_in(1:nz,0:nx+1,itime) +               &
            ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  &
            *(Tforce_in(1:nz,0:nx+1,itime+1)-Tforce_in(1:nz,0:nx+1,itime))
         
         qforce(1:nz,0:nx+1)=qforce_in(1:nz,0:nx+1,itime) +               &
            ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  &
            *(qforce_in(1:nz,0:nx+1,itime+1)-qforce_in(1:nz,0:nx+1,itime))
         
         wth_surf(0:nx+1)=wth_surf_t(0:nx+1,itime) +               &
            ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  &
            *(wth_surf_t(0:nx+1,itime+1)-wth_surf_t(0:nx+1,itime))
         
         wqv_surf(0:nx+1)=wqv_surf_t(0:nx+1,itime) +               &
            ((time-time_in(itime))/(time_in(itime+1)-time_in(itime)))  &
            *(wqv_surf_t(0:nx+1,itime+1)-wqv_surf_t(0:nx+1,itime))
       case(1)
         do k=1,nz
           Tforce(k,0:nx+1)= (thinit(k,0:nx+1) - sum(theta(k,0:nx+1))/(nx+2)) &
              /exner(k,0:nx+1)/Trelax(k)
           qforce(k,0:nx+1)= (qvinit(k,0:nx+1) - sum(qv(k,0:nx+1))/(nx+2))/Trelax(k)
         end do
       case(2)
         do k=1,nz
           Tforce(k,0:nx+1)= (thinit(k,0:nx+1) - theta(k,0:nx+1)) &
              /exner(k,0:nx+1)/Trelax(k)
           qforce(k,0:nx+1)= (qvinit(k,0:nx+1) - qv(k,0:nx+1))/Trelax(k)
         end do
       end select

    end if



  end subroutine interpolate_forcing

  integer function interval(x,x0)
    ! Find interval in 1D monotonic array x in which value x0 lies
    ! return the lower index of this interval

    real(wp) :: x(:)
    real(wp) :: x0
    integer :: i,nx
    
    nx=size(x) 
    if ((x0-x(1))*(x0-x(nx))>0.)then
       if (x(1)<x(nx))then
          if (x0<x(1))i=0
          if (x0>x(nx))i=nx
       else
          if (x0>x(1))i=0
          if (x0<x(nx))i=nx
       end if
    else
       do i=1,nx-1
          if (x(i) <= x0 .and. x(i+1) > x0) exit
       end do
    end if

    interval=i

  end function interval

  integer function nearest_loc(x,x0)
    !Find nearest value to x0 in a 1D monotonic array x

    real(wp) :: x(:)
    real(wp) :: x0
    integer :: i

    i=interval(x,x0)

    if (abs(x(i+1)-x0) <= abs(x(i)-x0))i=i+1

    nearest_loc=i

  end function nearest_loc

  subroutine make_wgrid(z,znew)
    !
    ! Subroutine to make a staggered w grid
    !
    real(wp), intent(in)  :: z(:)
    real(wp), intent(out) :: znew(:)
    
    !local variables
    integer :: k, nz

    integer :: kpm=1 !or 0
    
    nz=size(z(:))
    
    do k=2-kpm,nz-kpm
       znew(k)=.5*(z(k+kpm)+z(k+kpm-1))
    end do
    k=1+(nz-1)*kpm    
    znew(k)=2.*z(k)-znew(k-1)
    
  end subroutine make_wgrid

  subroutine make_vgrid(x,xnew)
    !
    ! Subroutine to make a staggered v grid
    !
    real(wp), intent(in)  :: x(0:)
    real(wp), intent(out) :: xnew(0:)
    
    !local variables
    integer :: j, nx

    integer :: jpm=1 !or 0
    
    nx=size(x(:))-1
    
    ! this assumes the 0 and nx+1 array bounds
    do j=1-jpm,nx-jpm
       xnew(j)=.5*(x(j)+x(j+jpm))
    end do

    j = nx
    xnew(j)=2.*x(j)-xnew(j-1)

  end subroutine make_vgrid
 

  subroutine interpolate(z,f,znew,fnew,scheme_id)
    !
    ! Interpolate f from z points to znew points
    !
    ! Currently only linear interpolation coded up
    ! as option 1.  More to go in later....
    !
    real(wp), intent(in)  :: z(:),f(:),znew(:)
    real(wp), intent(out) :: fnew(:)
    integer, intent(in), optional :: scheme_id

    !local variables
    integer :: cscheme_id, k,l,nz,nznew

    nz=size(z(:))
    nznew=size(znew(:))

    if (present(scheme_id))then
       cscheme_id=scheme_id
    else
       cscheme_id=1
    end if

    select case(cscheme_id)
    case(1)
       !linear interpolation
       do k=1,nznew
          l=interval(z,znew(k))
          if (l==0) l=l+1    !extrapolate
          if (l==nz) l=l-1 !extrapolate
          fnew(k)=f(l)+(f(l+1)-f(l))*(znew(k)-z(l))/(z(l+1)-z(l)) 
       end do
    end select

  end subroutine interpolate

  subroutine interpolate_x(x,f,xnew,fnew,scheme_id)
    !
    ! Interpolate f from z points to znew points
    !
    ! Currently only linear interpolation coded up
    ! as option 1.  More to go in later....
    !
    real(wp), intent(in)  :: x(:),f(:),xnew(:)
    real(wp), intent(out) :: fnew(:)
    integer, intent(in), optional :: scheme_id

    !local variables
    integer :: cscheme_id, k,l,nx,nxnew

    nx=size(x(:))
    nxnew=size(xnew(:))

    if (present(scheme_id))then
       cscheme_id=scheme_id
    else
       cscheme_id=1
    end if

    select case(cscheme_id)
    case(1)
       !linear interpolation
       do k=1,nxnew
          l=interval(x(1:nx-1),xnew(k))
          if (l==0) l=l+1    !extrapolate
          if (l>=(nxnew)) l=l-1 !extrapolate
          fnew(k)=f(l)+(f(l+1)-f(l))*(xnew(k)-x(l))/(x(l+1)-x(l))
       end do
       
    end select

  end subroutine interpolate_x

  subroutine smooth1D(field, field_out, n)
    ! Smooth field
    real(wp), intent(in) :: field(:)
    real(wp), intent(out) :: field_out(:)    
    integer, intent(in), optional :: n

    !local variables
    integer :: cn, i, k, nf

    if (present(n))then
       cn=n
    else
       cn=10
    end if

    nf=size(field)
    field_out(:)=field(:)
    do i=1,cn
       do k=2,nf
          field_out(k)=.5*(field_out(k)+field_out(k-1))
       end do
    end do

  end subroutine smooth1D

end module interpolation

