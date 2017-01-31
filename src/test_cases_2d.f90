! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Some 2D test cases
!
! This module follows the same format as test_cases
! but the cases below are 2-D cases. 

module test_cases_2d

  Use typeKind
  Use parameters, only : dt, num_h_moments, num_h_bins &
     , aero_N_init, aero_rd_init, aero_sig_init
  Use runtime, only : n_times, n_force_times
  Use common_physics
  Use column_variables
  Use input_variables
  Use derived_fields, only : calc_derived_fields
  Use diagnostics, only : save_dg, i_dgtime
  Use physconst, only : Ru=>R, r_on_cp, g, p0, pi, cp
  Use interpolation, only : interpolate, smooth1D, make_wgrid, &
     make_vgrid, interpolate_x 
  Use switches
  Use class_species, only: species_allocate
  Use test_cases, only: allocate_forcing, z2exner, set_standard_profile 
  Use namelists, only: amp_fact, set_Nc

  Use aerosols, only: moment_logn, awp=>wp, set_aerosol

  Implicit none

  !local variables
  integer, private :: k, itime, ih, j, i
  
  ! Levels for interpolation
  real(wp), allocatable :: &
        pHeight(:)         & ! height
       ,pTheta(:)          & ! theta
       ,pqv(:)             & ! qv
       ,pRH(:)               ! RH

  ! 1-D column for interpolation. This is assigned to 
  ! 2-D arrays 
  real(wp), allocatable :: &
        Theta_1d(:)          & ! theta
       ,qv_1d(:)               ! qv  

  real(wp), allocatable :: RH(:) ! RH
  real(wp) :: &
        maxT  & ! max time
       ,maxZ  & ! max height
       ,maxW  & ! max updraught velocity
       ,maxX    ! max horizontal extent
  real(wp) :: force_dt

  
contains

  subroutine standard_2d_cases(icase)

    integer, intent(in) :: icase

    !local variables
    real(wp) :: t     ! temporary local time variable

    ! for setting aerosol
    integer, parameter :: Ninit=naerosol ! Number of aerosol we can intialize
    integer :: indices(Ninit)
    logical :: lainits(Ninit)
    real(wp) :: Nds(Ninit)
    real(wp) :: sigmas(Ninit)
    real(wp) :: rds(Ninit)
    real(wp) :: densitys(Ninit)
    real(wp) :: fscale(nz)

    print *, '2d cases called', icase

    select case (icase)

    case(itest_2d)
      
      print*, 'This test case has now been removed - Sorry.'
      stop

    case(igcss_2d_Cu)
       !=====================================================
       ! 2d cumulus case based on morrison & grabowski (2006)
       !=====================================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=3000.
       if (all(tctrl==0.))tctrl(1)=3600.
       if (all(xctrl==0.))xctrl(1)=9000.
       if (all(wctrl==0.))wctrl(1)=1.0
       if (ipctrl==0)ipctrl=1

       maxZ=zctrl(1)
       maxX=xctrl(1)
       maxT=tctrl(1)
       maxW=wctrl(1)
       n_times=int(maxT/dt)

       n_force_times=int(maxT/100.0)

       call set_2D_Cu_thermo_field(maxZ, maxX)

       call allocate_forcing(nz,nx,n_force_times)

       call set_2D_Cu_wind_field(maxW, maxZ, maxX, n_force_times)

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 

    case(igcss_2d_ISDAC)
       !================================
       ! 2d Stratocumulus case based 
       ! on morrison & grabowski (2007)
       ! Adapted for mixed-phase 
       !================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=825. - 300
       if (all(tctrl==0.))tctrl(1)=7200.
       if (all(xctrl==0.))xctrl(1)=800.
       if (all(wctrl==0.))wctrl(1)=1.0
       if (all(lhf_ctrl==0.))lhf_ctrl(1)=0.
       if (all(lhf_ctrl==0.))shf_ctrl(1)=0.

       if (ipctrl==0)ipctrl=4

       maxZ=zctrl(1)
       maxX=xctrl(1)
       maxT=tctrl(1)
       maxW=wctrl(1)
       n_times= int(maxT/dt)
       
       n_force_times = 2

       call set_2D_ISDAC_thermo_field(maxZ, maxX)

       call allocate_forcing(nz,nx,n_force_times)

       call set_2D_Sc_wind_field(maxW, maxZ, maxX, n_force_times)
       
       ! set the thermodynamic forcing
       do k=1,nz
          qforce_in(k,:,:)=0.0
          Tforce_in(k,:,:)=0.0
       enddo

       field_mask(nz,:)=0.0

    case(igcss_2d_Sc)
       !================================
       ! 2d Stratocumulus case based 
       ! on morrison & grabowski (2007)
       !================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=1040.
       if (all(tctrl==0.))tctrl(1:2)=(/12600.,1800./)
       ! if (all(tctrl==0.))tctrl(1)=10800.
       if (all(xctrl==0.))xctrl(1)=2000.
       if (all(wctrl==0.))wctrl(1)=1.0
       if (all(lhf_ctrl==0.))lhf_ctrl(1:2)=(/0., 3./)
       if (all(lhf_ctrl==0.))shf_ctrl(1:2)=(/0., -3./)
       !if (all(lhf_ctrl==0.))lhf_ctrl(1)=3.
       !if (all(lhf_ctrl==0.))shf_ctrl(1)=-3.
       
       if (ipctrl==0)ipctrl=2

       maxZ=zctrl(1)
       maxX=xctrl(1)
       maxT=tctrl(1)
       maxW=wctrl(1)
       n_times= int(maxT/dt)

       n_force_times = 3
       force_dt = 7200.0
       !n_force_times = n_times

       call set_2D_Sc_thermo_field(maxZ, maxX)

       call allocate_forcing(nz,nx,n_force_times)

       call set_2D_Sc_wind_field(maxW, maxZ, maxX, n_force_times)
       
       ! set the thermodynamic forcing
       do itime=1,n_force_times
          t=(itime-1)*force_dt
          if (t<tctrl(2))then
             do k=1,nz
                qforce_in(k,:,itime)=lhf_ctrl(1)/(RLVAP*rho(k)*maxZ)
                Tforce_in(k,:,itime)=shf_ctrl(1)/(CP*rho(k)*maxZ)
             enddo
          else
            do k=1,nz
                qforce_in(k,:,itime)=lhf_ctrl(2)/(RLVAP*rho(k)*maxZ)
                Tforce_in(k,:,itime)=shf_ctrl(2)/(CP*rho(k)*maxZ)
             enddo
          endif
       enddo

        !do k=1,nz
        !   qforce_in(k,:,:)=lhf_ctrl(1)/(RLVAP*rho(k)*maxZ)
        !   Tforce_in(k,:,:)=shf_ctrl(1)/(CP*rho(k)*maxZ)
        !enddo

       field_mask(nz,:)=0.0

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 

    case(iwmo_case1)
       !================================
       ! 2d Stratocumulus case based 
       ! for wmo workshop 2012
       !================================
       ! Set default control values
       if (all(zctrl==0.))zctrl(1)=1525.
       if (all(tctrl==0.))tctrl(1)=1800.
       if (all(xctrl==0.))xctrl(1)=1500.
       if (all(wctrl==0.))wctrl(1)=1.0
       if (all(lhf_ctrl==0.))lhf_ctrl(1)=0.
       if (all(lhf_ctrl==0.))shf_ctrl(1)=0.

       if (ipctrl==0)ipctrl=2

       maxZ=zctrl(1)
       maxX=xctrl(1)
       maxT=tctrl(1)
       maxW=wctrl(1)
       n_times= int(maxT/dt)

       force_dt = 600.0
      
       n_force_times = int(maxT/force_dt)

       call set_WMO_Case1_thermo_field(maxZ, maxX)

       call allocate_forcing(nz,nx,n_force_times)

       call set_WMO_Case1_wind_field(maxW, maxZ, maxX, n_force_times, force_dt)
       
       iforce_method=2 ! relax back to local initial profile
       qforce_in(:,:,:)=0.0
       Tforce_in(:,:,:)=0.0
       ! set the relaxation timescale
       do k=1,nz
         Trelax(k)=300.*exp(z(k)/200.)
       enddo

       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 

       !field_mask(nz,:)=0.0

    case (igcss_2d_squall)
       !================================
       ! 2d Squall case
       !================================
       ! Set default control values
       if (all(zctrl==0.))then
         zctrl(1)=11000.
         zctrl(2)=4000.
         zctrl(3)=500. 
         zctrl(4)=12000. 
       end if
       if (all(tctrl==0.))tctrl(1)=10.
       if (all(xctrl==0.))then
         xctrl(1)=10000.
         xctrl(2)=40000.
         xctrl(3)=40000.
         xctrl(4)=120000.
         xctrl(5)=240000.
       end if
       if (all(wctrl==0.))then
         wctrl(1)=10.0
         wctrl(2)=1.3
         wctrl(3)=1.5
         wctrl(4)=4.0
         wctrl(5)=.1e-3
       end if
       if (ipctrl==0)ipctrl=2

       maxZ=zctrl(4)
       maxX=xctrl(5)
       maxT=tctrl(1)
       n_times= int(maxT/dt)
       
       n_force_times = 2

       call set_gate_thermo_profile(maxZ, maxX)

       call allocate_forcing(nz,nx,n_force_times)

       call set_squall(wctrl(1), wctrl(2), wctrl(3), wctrl(4), wctrl(5), & 
          zctrl(1), zctrl(2), zctrl(3),                                  &
          xctrl(1), xctrl(2), xctrl(3), xctrl(4)                         & 
          ,n_force_times)
       
       do ih=1,naerosol
         indices(ih)=ih
         if (aero_N_init(ih) > 0)lainits(ih) = .true.
         Nds(ih)    = aero_N_init(ih)
         sigmas(ih) = aero_sig_init(ih)
         rds(ih)    = aero_rd_init(ih)
       end do
       densitys(:) = 1777.
       fscale(:) = 1.

       call set_aerosol(Ninit, indices, lainits, Nds, sigmas, rds, densitys, fscale) 

     case default
       print*, ''
       print*, '================================='
       print*, 'Warning: 2D case chosen'
       print*, 'Test case value not recognized - ',icase
       print*, 'Was a 1D case intended?'
       print*, 'Exitting...'
       print*, '================================='
       print*, ''
       stop

    end select

    ! Save initial profiles
    thinit(:,:)=theta(:,:)
    qvinit(:,:)=qv(:,:)

  end subroutine standard_2d_cases

 subroutine set_2D_Cu_thermo_field(maxZ, maxX)
    !
    ! Set up the 2D field for cumulus based on 
    ! Morrison and Grabowski (2007)  
    !

    real(wp), intent(in) :: maxZ, maxX

    ! local allocatable arrays for temperature and presssure
    real(wp), allocatable :: &
         press_cu(:)  & ! pressure for Cu case
         ,temp_cu(:)  & ! temperature for Cu case
         ,rh_cu(:)     ! RH for Cu case

    real(wp) :: tempk, tempkm, delz, delt, tavi

    integer :: nlevs, km1

    nlevs = 26

    allocate(pHeight(nlevs))
    allocate(pTheta(nlevs))
    allocate(pqv(nlevs))
    allocate(theta_1d(nz))
    allocate(qv_1d(nz))

    allocate(press_cu(nlevs))
    allocate(temp_cu(nlevs))
    allocate(rh_cu(nlevs))

    ! pqv in g/kg
    pqv=(/14.5,  14.5,  14.5,  14.0,  13.7,  13.9,  13.9,   &
         10.3,  10.3,  10.0,   9.9,   8.9,   7.9,   4.0,   2.3, &
         1.2,   1.2,   0.9,   0.6,   2.0,   1.6,   0.4,   1.5, &
         0.9,   0.5,   0.4/)

    press_cu=(/1014., 1010., 1000.,  990.,  975.,  960.,  950., &
          925.,  900.,  875.,  850.,  825.,  800.,  790.,  775., & 
          765.,  755.,  745.,  730.,  715.,  700.,  650.,  600., &
          550.,  500.,  450./)

    temp_cu=(/298.4, 298.0, 296.8, 295.7, 295.0, 293.7, 293.1, &
         291.4, 290.0, 288.0, 286.5, 285.1, 284.2, 284.5, 284.1, &
         284.4, 283.4, 284.2, 284.4, 283.2, 282.0, 280.0, 275.7, &
         272.5, 269.5, 266./)

    if (l_cu_cold)then
      do k=1,nlevs
        rh_cu(k) =  pqv(k)/qsaturation(temp_cu(k),press_cu(k))
      end do
      
      temp_cu=temp_cu-20 ! reduce temperature for testing ice
      
      do k=1,nlevs
        pqv(k) =  rh_cu(k)*qsaturation(temp_cu(k),press_cu(k))
      end do
    end if

    ! set qv kg/kg
    pqv(:) = pqv(:)*1.e-3
    
    ! calculate theta from temp_cu
    do k = 1, nlevs
       ptheta(k) = temp_cu(k)*(1.e3/press_cu(k))**r_on_cp
    enddo
    
    ! calculate approximate height from pressure 
    pheight(1) = 0.0
    do k = 2, nlevs
       km1 = k-1
       tempk = ptheta(k) * (1.e3/press_cu(k))**(-r_on_cp) &
            * (1. + .6*pqv(k))
       tempkm = ptheta(km1) * (1.e3/press_cu(km1))**(-r_on_cp) &
            * (1. + .6*pqv(km1))
       
       delt=tempk-tempkm
       if(delt.gt.1.e-4) then
          tavi=log(tempk/tempkm)/delt
       else
          tavi=1./tempk
       endif
       
       delz=-ru/(tavi*g)*log(press_cu(k)/press_cu(km1))
       pheight(k) = pheight(km1) + delz
    enddo
 
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    do j=0, nx+1
       x(j)=maxX*j/float(nx)
    enddo
    
    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)   
    
    do j = 1,nx
       theta(:,j) = theta_1d(:)
       qv(:,j) = qv_1d(:)
    enddo

    p_surf=press_cu(1)*100.
    call z2exner
    
    call calc_derived_fields
    
    deallocate(qv_1d)
    deallocate(theta_1d)
    deallocate(pqv)
    deallocate(rh_cu)
    deallocate(pTheta)
    deallocate(pHeight)

    deallocate(press_cu)
    deallocate(temp_cu)

  end subroutine set_2D_Cu_thermo_field

 subroutine set_2D_Sc_thermo_field(maxZ, maxX)
    !
    ! Set up the 2D field for Stratocumulus based on 
    ! Morrison and Grabowski (2007)   
    !

    real(wp), intent(in) :: maxZ, maxX

    ! local allocatable arrays for temperature and presssure
    real(wp), allocatable :: &
         press_Sc(:)  & ! pressure for Cu case
         ,temp_Sc(:)     ! temperature for Cu case

    real(wp) :: tempk, tempkm, delz, delt, tavi, scal_height

    integer :: nlevs, km1

    nlevs = 2

    allocate(pHeight(nlevs))
    allocate(pTheta(nlevs))
    allocate(pqv(nlevs))
    allocate(theta_1d(nz))
    allocate(qv_1d(nz))

    allocate(press_Sc(nlevs))
    allocate(temp_Sc(nlevs))

    ! pqv in g/kg
    pqv=(/8.5, 8.5/)

    press_Sc=(/1000., 850./)
    !srf_pressure = 1015.
    !press_Sc(1) = srf_pressure

    temp_Sc=(/288., 288./)

    ! set qv kg/kg
    pqv(:) = pqv(:)*1.e-3
    
    ! theta and temp_Sc are the same
    ptheta(:) = temp_Sc(:)
    
    ! calculate approximate height from pressure 
    pheight(1) = 0.0
    do k = 2, nlevs
       km1 = k-1
       tempk = ptheta(k) * (1.e3/press_Sc(k))**(-r_on_cp) &
            * (1. + .6*pqv(k))
       tempkm = ptheta(km1) * (1.e3/press_Sc(km1))**(-r_on_cp) &
            * (1. + .6*pqv(km1))

       delt=tempk-tempkm
       if(delt.gt.1.e-4) then
          tavi=log(tempk/tempkm)/delt
       else
          tavi=1./tempk
       endif
       
       delz=-Ru/(tavi*g)*log(press_Sc(k)/press_Sc(km1))
       pheight(k) = pheight(km1) + delz
    enddo     
    
    do k=1,nz
       z(k)=(maxZ*(k-1)/float(nz))-(10.)
    end do
  
    do j=0, nx+1
       x(j)=maxX*j/float(nx)
    enddo
    

    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)   
    
    do j = 1,nx
       theta(:,j) = theta_1d(:)
       qv(:,j) = qv_1d(:)
    enddo

    theta(:,0) = theta(:,1)
    theta(:,nx+1) = theta(:,nx)

    qv(:,0) = qv(:,1)
    qv(:,nx+1) = qv(:,nx)   

    p_surf=p0
    call z2exner

    call calc_derived_fields

    deallocate(qv_1d)
    deallocate(theta_1d)
    deallocate(pqv)
    deallocate(pTheta)
    deallocate(pHeight)

    deallocate(press_Sc)
    deallocate(temp_Sc)

  end subroutine set_2D_Sc_thermo_field
  
 subroutine set_WMO_Case1_thermo_field(maxZ, maxX)
    !
    ! Set up the 2D field for Stratocumulus based on 
    ! WMO 2012 modelling Case1
    !

    real(wp), intent(in) :: maxZ, maxX
    real(wp) :: z0=1520. 
    
    ! z and z_half declared in column variables as 1:nz. Redefine here 
    ! so that that are declared as 0:nz for this case
    !real(wp) :: z(0:nz), z_half(0:nz)
    real(wp) :: full_zlevel_offset, full_xlevel_offset

    full_zlevel_offset = ((maxZ/float(nz)))/2.
    full_xlevel_offset = ((maxX/float(nx)))/2.

    do k=1,nz
       z(k)=(maxZ*(k-1)/float(nz))-full_zlevel_offset
    end do
    call make_wgrid(z, z_half)
    
    print *, z
    print *, z_half

    do j=0, nx+1
       x(j)=maxX*(j-1)/float(nx)-full_xlevel_offset
    end do    
    call make_vgrid(x,x_half)
 
    do k=1,nz
      if (z(k) <= z0)then
        theta(k,:) = 289.
        qv(k,:) = 7.5e-3
      else
        theta(k,:) = 303. + (z(k)-z0)**(1./3.)
        qv(k,:) = 0.5e-3
      end if
    end do

    theta(:,0) = theta(:,1)
    theta(:,nx+1) = theta(:,nx)

    qv(:,0) = qv(:,1)
    qv(:,nx+1) = qv(:,nx)   

    p_surf=p0
    call z2exner

    call calc_derived_fields

  end subroutine set_WMO_Case1_thermo_field
  

 subroutine set_2D_ISDAC_thermo_field(maxZ, maxX)
    !
    ! Set up the 2D field for Stratocumulus based on 
    ! Morrison and Grabowski (2007)   
    !

    real(wp), intent(in) :: maxZ, maxX

    ! local allocatable arrays for temperature and presssure
    real(wp), allocatable :: &
         press_Sc(:)  & ! pressure for Cu case
         ,temp_Sc(:)     ! temperature for Cu case

    real(wp) :: tempk, tempkm, delz, delt, tavi, scal_height

    integer :: nlevs, km1

    nlevs = 2

    allocate(pHeight(nlevs))
    allocate(pTheta(nlevs))
    allocate(pqv(nlevs))
    allocate(theta_1d(nz))
    allocate(qv_1d(nz))

    allocate(press_Sc(nlevs))
    allocate(temp_Sc(nlevs))

    ! pqv in g/kg
    pqv=(/1.5, 1.5/)

    press_Sc=(/980., 916./)

    temp_Sc=(/264.5, 264.5/)

    ! set qv kg/kg
    pqv(:) = pqv(:)*1.e-3
    
    ! theta and temp_Sc are the same
    ptheta(:) = temp_Sc(:)
    
    ! calculate approximate height from pressure 
    pheight(1) = 0.0
    do k = 2, nlevs
       km1 = k-1
       tempk = ptheta(k) * (1.e3/press_Sc(k))**(-r_on_cp) &
            * (1. + .6*pqv(k))
       tempkm = ptheta(km1) * (1.e3/press_Sc(km1))**(-r_on_cp) &
            * (1. + .6*pqv(km1))
       
       delt=tempk-tempkm
       if(delt.gt.1.e-4) then
          tavi=log(tempk/tempkm)/delt
       else
          tavi=1./tempk
       endif
       
       delz=-Ru/(tavi*g)*log(press_Sc(k)/press_Sc(km1))
       pheight(k) = pheight(km1) + delz
    enddo
 
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do
    
    do j=0, nx+1
       x(j)=maxX*j/float(nx)
    enddo
    
    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)   
    
    do j = 1,nx
       theta(:,j) = theta_1d(:)
       qv(:,j) = qv_1d(:)
    enddo

    theta(:,0) = theta(:,1)
    theta(:,nx+1) = theta(:,nx)

    qv(:,0) = qv(:,1)
    qv(:,nx+1) = qv(:,nx)   

    p_surf=press_Sc(1)*100.
    call z2exner

    call calc_derived_fields
               
    deallocate(qv_1d)
    deallocate(theta_1d)
    deallocate(pqv)
    deallocate(pTheta)
    deallocate(pHeight)

    deallocate(press_Sc)
    deallocate(temp_Sc)

  end subroutine set_2D_ISDAC_thermo_field

  subroutine set_2D_wind_field(ip, maxZ, maxX, n_times)

    integer, intent(in) :: ip, n_times ! case switch
    real(wp), intent(in) :: maxZ, maxX  

    select case(ip)
    case (1) 
      call set_2D_Cu_wind_field(1._wp,maxZ, maxX, n_times)
    case (2)
      call set_2D_Sc_wind_field(1._wp, maxZ, maxX, n_times)
    case (3)
      call orographic1(maxZ, maxX, n_times)
    case (4)
      call set_2D_Sc_wind_field(1._wp, maxZ, maxX, n_times)
    case default
       call set_2D_Cu_wind_field(1._wp, maxZ, maxX, n_times)
    end select

  end subroutine set_2D_wind_field

  subroutine set_2D_Cu_wind_field(maxW, maxZ, maxX, n_times)
    !
    ! Set up the 2D wind field for cumulus based on 
    ! Morrison and Grabowski (2007) 
    ! 
    ! Using formulation described in Appendix of
    ! MG07

    real(wp), intent(in) :: maxZ, maxX, maxW
    integer, intent(in) :: n_times

    !local variables
    real(wp) :: t     ! temporary local time variable
    real(wp) :: Lz
    real(wp) :: zp(nz+1), dzp(nz+1)
    real(wp) :: xp(nx+1), dxp(nx+1)
    real(wp) :: phi(0:nx+1,nz+1) 
                      ! streamfunction for cumulus 
                      ! convection
    real(wp) :: ux(nx+1,nz+1), uz(nx+1,nz+1)
    real(wp) :: zscale1, zscale2
                      ! depth of the inflow and 
                      ! outflow (respectively)
    real(wp) :: ampl0, ampl20, xscale0, ampa, ampb &
         , amp2a, amp2b, tscale1, tscale2, ampl, ampl2 &
         , xscale, t1, ztop, x0, zz1, zz2 &
         , xl, xx, zl, xcen, zsh
    real(wp) :: scal_fac & ! factor to scale w and v to 0 above 
                           ! zscale2
         , dw0 & ! change in w from zscale2 to 0
         , dv0   ! change in v from zscale2 to 0
    
    integer :: nxmid

    ! DETERMINE THE DEPTH OF THE INFLOW (ZSCALE1) 
    ! AND OUTFLOW (ZSCALE2)(in METERS)
    zscale1=1.7*1.e3
    zscale2=2.7*1.e3
 
    ! CALCULATE X AND Z DISTANCES (IN METERS)    
    dzp(1:nz) = dz(1:nz)
    dzp(nz+1) = dz(nz)
    
    dxp(1:nx) = dx(1:nx)
    dxp(nx+1) = dx(nx)

    do k=1,nz+1
       zp(k)=(k-1)*dzp(k)
    enddo
    do i=1,nx+1
       xp(i)=(i-1)*dxp(i)
    enddo
    
    !
    ! INITIAL DATA FOR THE STREAMFUNCTION. AMPL IS WMAX IN M/S,
    ! AMPL2 IS THE VALUE OF LINEAR SHEAR OVER THE ZSCALE2 DEPTH (M/S)
    ! XSCALE IS WIDTH OF THE UPDRAFT IN M
    ! AMPA AND AMPB CONTROL THE TEMPORAL FLUCTUATIONS OF WMAX 
    ! AMP2A AND AMP2B CONTROL THE TEMPORAL FLUCTUATIONS OF SHEAR (TILT)
    ! TSCALE1 AND TSCALE2 ARE PERIODS OF COSINE FLUCTUATIONS OF THE ABOVE

    ! AMPL0=1.0 original value
    AMPL0=maxW
    AMPL20=0.
    XSCALE0=1.8*1.e3
    AMPA=3.5
    AMPB=3.0
    AMP2A=0.6*1.e3
    AMP2B=0.5*1.e3
    tscale1=600.
    tscale2=900.

    do itime=1,n_times
       t=(itime-1)*100.
       time_in(itime)=t
       ! set parameters for each time period
       !
       !  AMPL AND AMPL2 VARYING IN TIME
       !
       if (t.lt.300.) then
          AMPL=AMPL0
          AMPL2=AMPL20
          XSCALE=XSCALE0
       elseif (t.ge.300..and.t.lt.900.) then
          AMPL=AMPL0
          AMPL2=AMP2A*(cos(pi*((t-300.)/tscale1 - 1.)) + 1.)
          XSCALE=XSCALE0
       elseif (t.ge.900..and.t.le.1500.) then
          AMPL=AMPA*(cos(pi*((t-900.)/tscale1 + 1.)) +1.) + AMPL0
          AMPL2=AMP2A*(cos(pi*((t-300.)/tscale1 - 1.)) + 1.)
          XSCALE=XSCALE0
       elseif (t.ge.1500..and.t.lt.2100.) then
          AMPL=AMPB*(cos(pi*(t-1500.)/tscale2) +1.) + 2.*AMPL0
          AMPL2=AMP2B*(cos(pi*((t-1500.)/tscale2 - 1.)) + 1.)
          XSCALE=XSCALE0
       elseif (t.ge.2100..and.t.lt.2400.) then
          AMPL=AMPB*(cos(pi*(t-1500.)/tscale2) +1.) + 2.*AMPL0
          AMPL2=AMP2B*(cos(pi*((t-1500.)/tscale2 - 1.)) + 1.)
          XSCALE=XSCALE0
       else
          t1=2400.
          AMPL=AMPB*(cos(pi*(t1-1500.)/tscale2) +1.) + 2.*AMPL0
          AMPL2=AMP2B*(cos(pi*((t1-1500.)/tscale2 - 1.)) + 1.)
          XSCALE=XSCALE0
       endif
       !
       ! ADOPT THE AMPL FOR STREAMFUNCTION CALCULATION TO BE W_MAX/K_x
       AMPL=AMPL/pi*XSCALE
       !
       ! DEFINE STREAMFUNCTION AS A FUNCTION OF HEIGHT
       ! ALSO, CENTRALIZE THE UPDRAFT TO OCCUPY ONLY THE INNER XSCALE
       ! OF THE DOMAIN
       !
       ZTOP=zp(nz+1)/ZSCALE1
       XCEN=.5*xp(nx+1)
       NXMID=NX/2+1
       X0=(xp(nx+1)-XSCALE)/2.
       DO I=1,NXMID
          DO K=1,NZ+1
             ZZ1=zp(k)/ZSCALE1
             ZZ2=zp(k)/ZSCALE2
             XL=2.*SQRT((xp(i)-XCEN)**2)
             XX=XL/XSCALE
             ZL=2.*ZSCALE1
             PHI(i,k)=0.
             IF (XX.LE.1.) THEN
                IF(ZZ1.LT.1.) THEN
                   PHI(i,k)=-cos(pi*(xp(i)-X0)/XSCALE)*sin(pi*zp(k)/ZL)
                   PHI(i,k)=PHI(i,k)*AMPL
                ELSEIF (ZZ1.GE.1..AND.ZZ2.LE.1.) THEN
                   ZL=2.*(ZSCALE2-ZSCALE1)
                   PHI(i,k)=-cos(pi*(xp(i)-X0)/XSCALE) &
                        *sin(pi*(.5+(zp(k)-ZSCALE1)/ZL))
                   PHI(i,k)=PHI(i,k)*AMPL
                ENDIF
             ELSE
                IF(ZZ1.LT.1.) THEN
                   PHI(i,k)=-sin(pi*zp(k)/ZL)
                   PHI(i,k)=PHI(i,k)*AMPL
                ELSEIF (ZZ1.GE.1..AND.ZZ2.LE.1.) THEN
                   ZL=2.*(ZSCALE2-ZSCALE1)
                   PHI(i,k)=-sin(pi*(.5+(zp(k)-ZSCALE1)/ZL))
                   PHI(i,k)=PHI(i,k)*AMPL 
                ENDIF
             ENDIF
          end do
       end do
       !
       ! APPLY THE SYMMETRY CONDITION
       DO I=NXMID,NX+1
          DO K=1,NZ+1
             PHI(i,k)=-PHI(NX+1-i+1,k)
          ENDDO
       ENDDO
       !
       !     ADD LINEAR SHEAR TO PRODUCE A WEAK TILT OF THE UPDRAFT
       DO K=1,NZ+1
          ZSH=zp(k)/ZSCALE2
          DO I=1,NX+1
             PHI(i,k)=PHI(i,k) - AMPL2*.5*ZSH**2.
          Enddo
       Enddo
       !  calculate rho*vel by derivation of streamfunction and normalize
       !  rho*ux velocity:
       do i=1,nx+1
          do k=1,nz
             ux(i,k)=-(phi(i,k+1)-phi(i,k))/dzp(k)  *dt/dxp(i)
          enddo
       enddo
       !  rho*uz velocity
       do k=1,nz+1
          do i=1,nx
             uz(i,k)=(phi(i+1,k)-phi(i,k))/dxp(i)  *dt/dzp(k)
          enddo
       enddo
       
       ! velocity fields
       do i=1,nx
          do k=1,nz
            v_t(k,i,itime)=0.5*(ux(i,k)+ux(i+1,k))/dt*dxp(i) /rho(k)
            w_t(k,i,itime)=0.5*(uz(i,k)+uz(i,k+1))/dt*dzp(k) /rho(k)
          end do
       end do

       v_t(:,0,itime) =  v_t(:,1,itime) 
       v_t(:,nx+1,itime) =  v_t(:,nx,itime) 
       w_t(:,0,itime) =  w_t(:,1,itime) 
       w_t(:,nx+1,itime) =  w_t(:,nx,itime)           
       do j = 0, nx+1 
         call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
       enddo
       do k = 1, nz 
         call interpolate_x(x,v_t(k,:,itime),x_half,v_t_half(k,:,itime),scheme_id=1) 
       enddo
       
       call wind_chk(w_t_half(:,:,itime), v_t_half(:,:,itime)) 
     end do
       
  end subroutine set_2D_Cu_wind_field  

  subroutine set_2D_Sc_wind_field(w_peak, maxZ, maxX, n_times)
    !
    ! Set up the 2D wind field for Stratocumulus based on 
    ! Morrison and Grabowski (2007) 
    ! 
    real(wp), intent(in) :: w_peak  ! peak updraught strenght
    real(wp), intent(in) :: maxZ, maxX
    integer, intent(in) :: n_times

    !local variables
    real(wp) :: t     ! temporary local time variable
    !real(wp) :: zp(nz+1), dzp(nz+1)
    real(wp) :: zp(nz), dzp(nz)
    real(wp) :: xp(nx+1), dxp(nx+1)
    real(wp) :: phi(1:nx+1,nz+1)
                      ! streamfunction for cumulus 
                      ! convection
    !real(wp) :: ux(nx+1,nz+1), uz(nx+1,nz+1)
    real(wp) :: ux(nx+1,nz), uz(nx,nz+1)

    real(wp) :: zscale, xscale
                      ! depth of the inflow and 
                      ! outflow (respectively)
    real(wp) :: ampl, ztop, x0, xcen
    real(wp) :: field_2d(nz,1:nx)
    
    dzp(1:nz) = dz(1:nz)
    !dzp(nz+1) = dz(nz)
    
    dxp(1:nx) = dx(1:nx)
    dxp(nx+1) = dx(nx)

    !do k=1,nz+1
    do k=1,nz
       zp(k)=((k-1)*dzp(k))-(dz_half(1))
    enddo
    do i=1,nx+1
       xp(i)=(i-1)*dxp(i)
    enddo

    XSCALE=xp(nx+1)
    ZSCALE=1000.0 ! set ZSCALE = 1020 m, while Z(NZ) = 1020 m
   ! print *, zscale
   ! print *, z_half

!
!! INITIAL DATA FOR THE STREAMFUNCTION. AMPL IS WMAX IN M/S,
!
    AMPL=w_peak
!
! ADOPT THE AMPL FOR STREAMFUNCTION CALCULATION TO BE W_MAX/K_x
    AMPL=0.5*AMPL/pi*XSCALE
!
! DEFINE STREAMFUNCTION AS A FUNCTION OF HEIGHT
! ALSO, CENTRALIZE THE UPDRAFT TO OCCUPY ONLY THE INNER XSCALE
! OF THE DOMAIN
!
    ZTOP=zp(nz)/ZSCALE
    XCEN=.5*xp(nx)
    X0=(xp(nx+1)-XSCALE)/2.
    phi = 0.0
    DO I=1,NX+1
       DO K=1,NZ+1
          PHI(i,k)=-cos(2.*pi*(xp(i)-X0)/XSCALE)* &
               sin(pi*zp(k)/ZSCALE)
          PHI(i,k)=PHI(i,k)*AMPL
       ENDDO
    ENDDO

    do itime=1,n_times
       t=(itime-1)*maxT
       time_in(itime)=t

       !  calculate rho*vel by derivation of streamfunction and normalize
       !  rho*ux velocity:
       do i=1,nx+1
          do k=1,nz
             ux(i,k)=-(phi(i,k+1)-phi(i,k))/dzp(k)*dt/dxp(i)
          enddo
       enddo
       !  rho*uz velocity
       do k=1,nz+1
          do i=1,nx
             uz(i,k)=(phi(i+1,k)-phi(i,k))/dxp(i)*dt/dzp(k)
          enddo
       enddo
       ! ****CHECK RHO CALC****
       ! velocity fields
       
       do i=1,nx
          !w_t(1,i,itime) = 0.0 ! ensure 0.0 velocity at the surface
          do k=1,nz
             v_t(k,i,itime)=0.5*(ux(i,k)+ux(i+1,k))/dt*dxp(i) /rho(k)
             w_t(k,i,itime)=0.5*(uz(i,k)+uz(i,k+1))/dt*dzp(k) /rho(k)
          end do
       end do

       v_t(:,0,itime) = v_t(:,nx,itime)
       v_t(:,nx+1,itime) = v_t(:,1,itime)
       
       w_t_half(:,0,itime) = w_t_half(:,nx,itime)
       w_t_half(:,nx+1,itime) = w_t_half(:,1,itime)
        

        do j = 0, nx+1 
          call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
       enddo

        do k = 1, nz 
          call interpolate_x(x,v_t(k,:,itime),x_half,v_t_half(k,:,itime),scheme_id=1) 
        enddo
        
        v_t_half(:,0,itime) = v_t_half(:,1,itime)
        v_t_half(:,nx+1,itime) = v_t_half(:,nx,itime)

        w_t_half(:,0,itime) = w_t_half(:,1,itime)
        w_t_half(:,nx+1,itime) = w_t_half(:,nx,itime) 
        

!!$        if (itime < 2) then 
!!$           field_2d(:,1:nx)= w_t(:,1:nx,itime)
!!$           call save_dg( field_2d, 'w before wind_chk', itime, units='m/s', dim='z,x')
!!$           field_2d(:,1:nx)= w_t_half(:,1:nx,itime)
!!$           call save_dg( field_2d, 'w_half before wind_chk', itime, units='m/s', dim='z,x')
!!$           field_2d(:,1:nx)= v_t(:,1:nx,itime)
!!$           call save_dg( field_2d, 'v before wind_chk', itime, units='m/s', dim='z,x')
!!$           field_2d(:,1:nx)= v_t_half(:,1:nx,itime)
!!$           call save_dg( field_2d, 'v_half before wind_chk', itime, units='m/s', dim='z,x')
!!$         
!!$           do k = 1,nz
!!$              do i = 1,nx
!!$                 field_2d(k,i)= uz(i,k)
!!$              enddo
!!$           enddo
!!$           call save_dg( field_2d, 'uz', itime, units='m/s', dim='z,x')
!!$           do k = 1,nz
!!$              do i = 1,nx
!!$                 field_2d(k,i)= phi(i,k)
!!$              enddo
!!$           enddo
!!$           call save_dg( field_2d, 'phi', itime, units='m/s', dim='z,x')
!!$           
!!$        endif

        call wind_chk(w_t_half(:,:,itime), v_t_half(:,:,itime))
        
!!$        if (itime < 2) then 
!!$           field_2d(:,1:nx)= w_t(:,1:nx,itime)
!!$           call save_dg( field_2d, 'w after wind_chk', itime, units='m/s', dim='z,x')
!!$           field_2d(:,1:nx)= w_t_half(:,1:nx,itime)
!!$           call save_dg( field_2d, 'w_half after wind_chk', itime, units='m/s', dim='z_half,x')
!!$           field_2d(:,1:nx)= v_t(:,1:nx,itime)
!!$           call save_dg( field_2d, 'v after wind_chk', itime, units='m/s', dim='z,x')
!!$           field_2d(:,1:nx)= v_t_half(:,1:nx,itime)
!!$           call save_dg( field_2d, 'v_half after wind_chk', itime, units='m/s', dim='z,x')
!!$        endif
     end do
!!$     print *, w_t_half(nz, :, 1), z_half(nz)

  end subroutine set_2D_Sc_wind_field

  subroutine set_WMO_Case1_wind_field(w_peak, maxZ, maxX, n_times, force_dt)
    !
    ! Set up the 2D wind field for Stratocumulus based on 
    ! WMO 2012 Case 1 
    ! 
    real(wp), intent(in) :: w_peak  ! peak updraught strenght
    real(wp), intent(in) :: maxZ, maxX, force_dt
    integer, intent(in) :: n_times 

    !local variables
    real(wp) :: t     ! temporary local time variable
    real(wp) :: zp(0:nz+1), dzp(0:nz+1)
    real(wp) :: xp(0:nx+1), dxp(0:nx+1)
    real(wp) :: phi(0:nx+1,0:nz+1)
                      ! streamfunction for cumulus 
                      ! convection
    real(wp) :: ux(0:nx+1,1:nz+1), uz(0:nx+1,1:nz+1)
    real(wp) :: zscale, xscale
                      ! depth of the inflow and 
                      ! outflow (respectively)
    real(wp) :: ampl(n_times)
    real(wp) :: ztop, x0, xcen
    
    ! CALCULATE X AND Z DISTANCES (IN METERS)
    dzp(1:nz) = dz(1:nz)
    dzp(nz+1) = dz(nz)
    ! Set up the z grid for defining phi so that phi is 
    ! defined on equivalent of KiD full levels. This 
    ! grid adds a 0 level and an extra level at the top
    zp(0) = z_half(1)-dzp(1)
    zp(1:nz) = z_half(1:nz)
    zp(nz+1) = z_half(nz)+dzp(nz)
    
    dxp(0:nx) = dx(0:nx)
    dxp(nx+1) = dx(nx)
    
    ! KiD x-grid is set-up correctly for phi so just 
    ! equate xp and x
    xp(0:nx+1) = x_half(0:nx+1)

    XSCALE=xp(nx+1)
    ZSCALE=zp(nz)
!
!! INITIAL DATA FOR THE STREAMFUNCTION. AMPL IS WMAX IN M/S,
!
    do itime=1,n_times
       t=(itime-1)*force_dt
       time_in(itime)=t
       print *, t
       if (t <= 600.0) then
          AMPL(itime)=0.0
       else
          AMPL(itime)=w_peak
       endif
       ! ADOPT THE AMPL FOR STREAMFUNCTION CALCULATION TO BE W_MAX/K_x
       AMPL(itime)=0.5*AMPL(itime)/pi*XSCALE
    enddo
          
!
! ADOPT THE AMPL FOR STREAMFUNCTION CALCULATION TO BE W_MAX/K_x
!    AMPL(itime)=0.5*AMPL/pi*XSCALE
!    AMPL=AMPL/pi*XSCALE
!
! DEFINE STREAMFUNCTION AS A FUNCTION OF HEIGHT
! ALSO, CENTRALIZE THE UPDRAFT TO OCCUPY ONLY THE INNER XSCALE
! OF THE DOMAIN
!
    do itime=1,n_times
       !t=(itime-1)*maxT
       !time_in(itime)=t
       !print *, time_in
       XCEN=.5*xp(nx+1)
       X0=0.0!(xp(nx+1)-XSCALE)/2.
       phi = 0.0
       DO I=0,NX+1
          DO K=0,NZ+1
             PHI(i,k)=-cos(2.*pi*(xp(i)-X0)/XSCALE)* &
                  sin(pi*zp(k)/ZSCALE)
             PHI(i,k)=PHI(i,k)*AMPL(itime)
          ENDDO
       ENDDO


       !  calculate rho*vel by derivation of streamfunction and normalize
       !  rho*ux velocity:
       do i=0,nx+1
          do k=1,nz+1
             ux(i,k)=-(phi(i,k)-phi(i,k-1))/dzp(k)
          enddo
       enddo
       !  rho*uz velocity
       do k=1,nz
          do i=1,nx+1
             uz(i,k)=(phi(i,k)-phi(i-1,k))/dxp(i)
          enddo
       enddo
      
       do i=1,nx+1
          do k=1,nz
             v_t_half(k,i,itime)=(ux(i,k))/rho(k)
             w_t_half(k,i,itime)=(uz(i,k))/rho_half(k)
          end do
       end do  
       
       !!AH - setting the boundaries like below causes a divergence 
       !!     in the first column
       !set the extra 0 and nx+1 columns
        v_t_half(:,0,itime) = v_t_half(:,nx,itime)
        v_t_half(:,nx+1,itime) = v_t_half(:,1,itime)

        w_t_half(:,0,itime) = w_t_half(:,nx,itime)
        w_t_half(:,nx+1,itime) = w_t_half(:,1,itime)

        call wind_chk(w_t_half(:,:,itime), v_t_half(:,:,itime))

        !!AH - interpolate w_half onto w so w exists on full level
        !!     This may be needed by some schemes
        do j = 0, nx+1 
          call interpolate(z_half,w_t_half(:,j,itime),z,w_t(:,j,itime),scheme_id=1) 
       enddo

       do k = 1, nz 
          call interpolate_x(x_half,v_t_half(k,:,itime),x,v_t(k,:,itime),scheme_id=1) 
       enddo

    end do    
    
    !k = 10 
    !j = 1
    !print *, v_t_half(k,j-1,1),v_t_half(k,j,1),phi(j-1,k),phi(j,k) 
    !j = 2
    !print *, v_t_half(k,j-1,1),v_t_half(k,j,1),phi(j,k)
    !print *, 'X-grid'
    !do j = 0, nx+1
    !   print *, x_half(j), xp(j),  v_t_half(1,j, 1) 
    !enddo
    !print *, x0, xcen
    
  end subroutine set_WMO_Case1_wind_field

  subroutine orographic1(maxZ, maxX, n_times)
    !
    ! Set up the 2D wind field for Stratocumulus based on 
    ! Morrison and Grabowski (2007) 
    ! 
    real(wp), intent(in) :: maxZ, maxX
    integer, intent(in) :: n_times

    !local variables
    real(wp) :: t     ! temporary local time variable
    real(wp) :: zp(nz+1), dzp(nz+1)
    real(wp) :: xp(nx+1), dxp(nx+1)
    real(wp) :: phi(0:nx+1,nz+1) 
                      ! streamfunction for flow
                      ! around two dipoles (taken
                      ! as proxy for 2 simple hills)
    real(wp) :: ux(nx+1,nz+1), uz(nx+1,nz+1)
    real(wp) :: wmax, xcentre
    real(wp) :: xa, xb, pratio, ha, hb, ca, cb, wa, wb
    real(wp) :: Ra, Rb, thresh, zoff
    real(wp) :: znd, xnd

    
    !=========================================
    ! Examine grid
    !=========================================
    dzp(1:nz) = dz(:)
    dzp(nz+1) = dz(nz)
    
    dxp(1:nx) = dx(1:nx)
    dxp(nx+1) = dx(nx)

    do k=1,nz+1
       zp(k)=(k-1)*dzp(k)
    enddo
    do i=1,nx+1
       xp(i)=(i-1)*dxp(i)
    enddo

    xcentre=.5*xp(nx+1)

    !=========================================
    ! Set parameters defining the hills
    !=========================================
    pratio = 1.05  ! ratio of circleA radius to hillA height

    xa = -2000.  ! location of hill a 
    xb = 2000. ! location of hill b 
    xa = xa + xcentre
    xb = xb + xcentre

    ha = 1800.  ! height of hill a
    hb = 1300.  ! height of hill b

    ca = ha/pratio  ! height of circle a
    cb = hb/pratio  ! height of circle b

    wa = 2000.  ! width of hill a
    wb = 2000.  ! width of hill a

    !=========================================
    ! Determine streamline to use as hills
    !=========================================

    Ra=pratio
    Rb=sqrt((xa-xb)*(xa-xb)/(wb*wb)+ha*ha/(cb*cb))
    thresh = ha*(1.-1./(Ra*Ra)) +  ha*(1.-1./(Rb*Rb))

    zoff=thresh/2.

    !=========================================
    ! Set amplitude  
    !=========================================
    wmax = 4.
    
    !=========================================
    ! Set up potential
    !=========================================
    DO I=1,nx+1
      xnd=xp(i)
       DO K=1,nz+1
         znd=(zp(k)+zoff)
         Ra=sqrt((xnd-xa)*(xnd-xa)/(wa*wa)+znd*znd/(ca*ca))
         Rb=sqrt((xnd-xb)*(xnd-xb)/(wb*wb)+znd*znd/(cb*cb))
         phi(i,k) = znd*(1.-1./(Ra*Ra)) +  znd*(1.-1./(Rb*Rb))

         if (phi(i,k) < .5*thresh)then
           field_mask(k,i)=0.0
         end if
         phi(i,k)=-phi(i,k)*wmax
       ENDDO
    ENDDO

    call save_dg(field_mask(1:nz,1:nx)*transpose(phi(1:nx, 1:nz)), 'phi', 1)
    
    !=========================================
    ! Calculate velocities
    !=========================================

    do itime=1,n_times
       t=itime*dt
       time_in(itime)=t

       !  calculate rho*vel by derivation of streamfunction and normalize
       !  rho*ux velocity:
       do i=1,nx+1
          do k=1,nz
             ux(i,k)=-(phi(i,k+1)-phi(i,k))/dzp(k)*dt/dxp(i)
          enddo
       enddo
       !  rho*uz velocity
       do k=1,nz+1
          do i=1,nx
             uz(i,k)=(phi(i+1,k)-phi(i,k))/dxp(i)*dt/dzp(k)
          enddo
       enddo
       ! ****CHECK RHO CALC****
       ! velocity fields
       do i=1,nx
          do k=1,nz
             v_t(k,i,itime)=0.5*(ux(i,k)+ux(i+1,k))/dt*dxp(i) /rho(k)
             w_t(k,i,itime)=0.5*(uz(i,k)+uz(i,k+1))/dt*dzp(k) /rho(k)
          end do
       end do

       do j = 0, nx+1 
         call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
       enddo
       do k = 1, nz 
         call interpolate_x(x,v_t(k,:,itime),x_half,v_t_half(k,:,itime),scheme_id=1) 
       enddo
       
       call wind_chk(w_t_half(:,:,itime), v_t_half(:,:,itime)) 
     end do
       
  end subroutine orographic1


  subroutine set_squall(A1, A2, A3, Q, S, & 
     Z1, Z2, Z3, X1, X2, xc1, xc2, n_times)
    !
    ! Set up the 2D squall line case
    ! c.f. Slawinska et al 2009, QJRMS
    ! 
    real(wp), intent(inout) :: A1, A2, A3, S, Q, & 
     Z1, Z2, Z3, X1, X2, xc1, xc2
    integer, intent(in) :: n_times

    !local variables
    real(wp) :: t     ! temporary local time variable
    real(wp) :: zp(nz), dzp(nz)
    real(wp) :: xp(nx+1), dxp(nx+1)
    real(wp) :: zm1(nz+1),zm2(nz+1),zm3(nz+1)
    real(wp) :: xm1(nx+1),xm2(nx+1)
    real(wp) :: phi(0:nx+1,nz+1)  ! streamfunction 
    real(wp) :: ux(nx+1,nz+1), uz(nx+1,nz+1)
    real(wp) :: zscale, xscale
                      ! depth of the inflow and 
                      ! outflow (respectively)
    real(wp) :: AX, phi_tmp
    
    ! CALCULATE X AND Z DISTANCES (IN METERS)

    dzp(1:nz) = dz(1:nz)
    
    dxp(1:nx) = dx(1:nx)
    dxp(nx+1) = dx(nx)

    do k=1,nz
       zp(k)=(k-1)*dzp(k)
       zm1(k) = min(zp(k), Z1)
       zm2(k) = min(Z2, max(0., zp(k)-Z3))
       zm3(k) = min(Z2, max(0., zp(k)-Z3-Z2))
    enddo
    do i=1,nx+1
       xp(i)=(i-1)*dxp(i)
       xm1(i) = max(-X1, min(X1, xp(i)-xc1))
       xm2(i) = max(-X2, min(X2, xp(i)-xc2))
    enddo

!
    phi = 0.0

    DO I=1,NX+1
      DO K=1,NZ
        ! Convective part
        PHI(i,k)=(A1*2*X1/pi)*sin(0.5*pi*xm1(i)/X1)*sin(pi*zm1(k)/Z1)
!         ! Convective part is additionally skewed by rho
        PHI(i,k) = PHI(i,k)*rho(k)
      ENDDO
    END DO
    ! renormalize convective part given extra rho factor
    ! so that A1 still give peak amplitude
    AX = (A1*2*X1/pi)/maxval(phi(:,:))
    phi(:,:) = phi(:,:)*AX
    DO I=1,NX+1
       DO K=1,NZ
        ! Large scale constant + shear
         phi_tmp = (-0.5*S*zp(k)*zp(k) - Q*zp(k))/rho(k)
         PHI(i,k)=Phi(i,k) + phi_tmp
        ! Lower strat downdraft
         phi_tmp= -(A2*2.0*X2/pi)*sin(0.5*pi*xm2(i)/X2)*sin(pi*zm2(k)/Z2)**2 
         PHI(i,k)=Phi(i,k) + phi_tmp
        ! Upper strat updraft
         phi_tmp= (A3*2.0*X2/pi)*sin(0.5*pi*xm2(i)/X2)*sin(pi*zm3(k)/Z2)**2 
         PHI(i,k)=Phi(i,k) + phi_tmp
         !multiply by rho so that amplitudes reflect velocity not mass flux
         PHI(i,k) = rho(k)*PHI(i,k)
       ENDDO
    ENDDO

    !  calculate rho*vel by derivation of streamfunction and normalize
    !  rho*ux velocity:
    do k=1,nz-1
      do i=1,nx+1
        ux(i,k)=-(phi(i,k+1)-phi(i,k))/dzp(k)*dt/dxp(i)
      enddo
    enddo
    !  rho*uz velocity
    do k=1,nz-1
      do i=1,nx
        uz(i,k)=(phi(i+1,k)-phi(i,k))/dxp(i)*dt/dzp(k)
      enddo
    enddo

       ! velocity fields
    do itime=1,n_times
       t=(itime-1)*maxT
       time_in(itime)=t
       do i=1,nx
          do k=1,nz-1
             v_t(k,i,itime)=0.5*(ux(i,k)+ux(i+1,k))/dt*dxp(i) /rho(k)
             w_t(k,i,itime)=0.5*(uz(i,k)+uz(i,k+1))/dt*dzp(k) /rho(k)
          end do
       end do

        v_t(:,0,itime) = v_t(:,nx,itime)
        v_t(:,nx+1,itime) = v_t(:,1,itime)

        w_t(:,0,itime) = w_t(:,nx,itime)
        w_t(:,nx+1,itime) = w_t(:,1,itime) 
        	 
        do j = 0, nx+1 
          call interpolate(z,w_t(:,j,itime),z_half,w_t_half(:,j,itime),scheme_id=1) 
        enddo
        do k = 1, nz 
          call interpolate_x(x,v_t(k,:,itime),x_half,v_t_half(k,:,itime),scheme_id=1) 
        enddo
        
        call wind_chk(w_t_half(:,:,itime), v_t_half(:,:,itime)) 
      end do
       
  end subroutine set_squall

  subroutine set_gate_thermo_profile(maxZ, maxX)
    !
    ! Set up the 2D field for GATE 
    !

    real(wp), intent(in) :: maxZ, maxX

    ! local allocatable arrays for temperature and presssure
    real(wp), allocatable :: &
         press_cu(:)  & ! pressure (mb)
         ,temp_cu(:)    ! temperature  (C)

    real(wp) :: tempk, tempkm, delz, delt, tavi

    integer :: nlevs, km1

    nlevs = 23

    allocate(pHeight(nlevs))
    allocate(pTheta(nlevs))
    allocate(pqv(nlevs))
    allocate(theta_1d(nz))
    allocate(qv_1d(nz))

    allocate(press_cu(nlevs))
    allocate(temp_cu(nlevs))

    ! pqv in g/kg
    pqv=(/.178E+02, 0.172E+02, 0.156E+02, 0.134E+02, 0.111E+02, &
       .888E+01, 0.631E+01, 0.487E+01, 0.396E+01, 0.200E+01,    &
       .984E+00, 0.806E+00, 0.370E+00, 0.135E+00, 0.599E-01,    &
       .258E-01, 0.123E-01, 0.582E-02, 0.367E-02, 0.589E-02,    &
       .104E-02, 0.247E-02, 0.585E-02/)

    press_cu=(/1008.00, 991.25, 945.50, 893.79, 836.06, 772.82, 705.22, &
       635.05, 564.48, 495.73, 430.71, 370.78, 316.72, 268.82,          &
       226.98, 190.82, 159.87, 133.55, 111.29,  92.56,  52.31,          &
       22.08,   9.32/)

    temp_cu=(/25.26,  24.13,  21.04,  18.66,  16.50,  13.41, 9.06, &
       3.73,  -1.51,  -6.97, -14.09, -22.44, -30.57, -39.60,       &
       -48.69, -57.40, -65.21, -72.58, -76.71, -74.98, -74.98,     &
       -74.98, -74.98/)

    ! set temp K
    temp_cu=temp_cu + 273.15 
    
    ! set qv kg/kg
    pqv(:) = pqv(:)*1.e-3
    
    ! calculate theta from temp_cu
    do k = 1, nlevs
       ptheta(k) = temp_cu(k)*(1.e3/press_cu(k))**r_on_cp
    enddo
    
    ! calculate approximate height from pressure 
    pheight(1) = 0.0
    do k = 2, nlevs
       km1 = k-1
       tempk = ptheta(k) * (1.e3/press_cu(k))**(-r_on_cp) &
            * (1. + .6*pqv(k))
       tempkm = ptheta(km1) * (1.e3/press_cu(km1))**(-r_on_cp) &
            * (1. + .6*pqv(km1))
       
       delt=tempk-tempkm
       if(delt.gt.1.e-4) then
          tavi=log(tempk/tempkm)/delt
       else
          tavi=1./tempk
       endif
       
       delz=-Ru/(tavi*g)*log(press_cu(k)/press_cu(km1))
       pheight(k) = pheight(km1) + delz
    enddo
 
    do k=1,nz
       z(k)=maxZ*k/float(nz)
    end do

    do j=0, nx+1
       x(j)=maxX*j/float(nx)
    enddo
    
    call interpolate(pHeight, ptheta, z, theta_1d)
    call interpolate(pHeight, pqv, z, qv_1d)   
    
    do j = 1,nx
       theta(:,j) = theta_1d(:)
       qv(:,j) = qv_1d(:)
    enddo

    p_surf=press_cu(1)*100.
    call z2exner
    
    call calc_derived_fields
    
    deallocate(qv_1d)
    deallocate(theta_1d)
    deallocate(pqv)
    deallocate(pTheta)
    deallocate(pHeight)

    deallocate(press_cu)
    deallocate(temp_cu)

  end subroutine set_gate_thermo_profile
          
 subroutine wind_chk(w_tmp, v_tmp)
    
   ! This subrourine takes in w and v fields and modifys them 
   ! to ensure non-divergence.
   ! First, the v winds are corrected while the w winds are kept 
   ! constant. Starting at j=0 (v kept constant here) then 
   ! working across the grid.
   ! Finally if periodic bounds are used and there is a global
   ! divergence, then the w field at j=1 is modified to
   ! compensate.  
   !
   ! Note that w and v are assumed to be on half levels
   !
   !
   ! Care should be taken when using this routine, that the flow
   ! field is not grossly modified.

    real(wp), intent(inout) :: w_tmp(nz,0:nx+1)
    real(wp), intent(inout) :: v_tmp(nz,0:nx+1)
    integer :: j, k

    ! local variables
    real(wp) :: delta_w(nz,0:nx+1), delta_v(nz,0:nx+1)
    real(wp) :: div(nz,0:nx+1), dv, dw

    ! Diagnose initial divergence
    do j = 1, nx+1
      do k = 2,nz
        
        delta_w(k,j) = w_tmp(k,j)*rho_half(k) &
           - w_tmp(k-1,j)*rho_half(k-1)
        delta_v(k,j) = (v_tmp(k,j) - v_tmp(k,j-1))*rho(k)

        div(k,j) = delta_w(k,j)/dz(k) + delta_v(k,j)/dx(j)
        if (j<nx+1)call save_dg(k,j,div(k,j), 'div1', 1)
        
      end do
    end do

    do j = 2, nx+1
      do k = 2,nz
        
          dv = - delta_w(k,j)*dx_half(j)/dz_half(k) 
          !v_tmp(k,j) = v_tmp(k,j-1) + dv/rho(k)
      enddo

     ! v_tmp(1,j) = v_tmp(2,j) 
    enddo

    ! periodic condition...
    !v_tmp(:,0) = v_tmp(:,nx)

    if (l_periodic_bound)then ! must ensure global non-divergence so alter w in column 1
      j=1
      do k=2,nz
        dw = - (v_tmp(k,j) - v_tmp(k,j-1))*rho(k)*dz_half(k)/dx_half(j)
        !w_tmp(k,j) = (w_tmp(k-1,j)*rho_half(k-1) + dw)/rho_half(k)
      end do
    end if

    do j = 1, nx+1
       do k = 2,nz

        delta_w(k,j) = w_tmp(k,j)*rho_half(k) &
           - w_tmp(k-1,j)*rho_half(k-1)
        delta_v(k,j) = (v_tmp(k,j) - v_tmp(k,j-1))*rho(k)
        
        div(k,j) = delta_w(k,j)/dz(k) + delta_v(k,j)/dx(j)
        if (j<nx+1)call save_dg(k,j,div(k,j), 'div2', 1)
         
      end do
    end do
    
  end subroutine wind_chk


end module test_cases_2d
