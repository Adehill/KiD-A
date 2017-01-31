! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Physics which is used by the driver, e.g. a standard 
! saturation calculation
!
!
Module common_physics

  Use typeKind
  Use physconst

  implicit none

  Interface qsaturation
     module procedure qsaturation_sp, qsaturation_dp
  end Interface

  Interface qisaturation
     module procedure qisaturation_sp, qisaturation_dp
  end Interface
  
contains

  REAL(SP) FUNCTION Qsaturation_sp (T, p)  
    IMPLICIT NONE  
    !
    ! Function to return the saturation mr over water
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
    !
    real(sp), INTENT (IN) :: T, p  
    ! Temperature in Kelvin
    ! Pressure in mb
    real(sp), PARAMETER::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882, &
         qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
    if (T > qsa3 .and. p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            > qsa4)then
       Qsaturation_sp = qsa1 / (p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            - qsa4)
    else
       qsaturation_sp=0.
    end if

  END FUNCTION Qsaturation_sp

  REAL(DP) FUNCTION Qsaturation_dp (T, p)  
    IMPLICIT NONE  
    !
    ! Function to return the saturation mr over water
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
    !
    real(dp), INTENT (IN) :: T, p  
    ! Temperature in Kelvin
    ! Pressure in mb
    real(dp), PARAMETER::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882, &
         qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
    if (T > qsa3 .and. p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            > qsa4)then
       Qsaturation_dp = qsa1 / (p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            - qsa4)
    else
       qsaturation_dp=0.
    end if

  END FUNCTION Qsaturation_dp

  REAL(SP) FUNCTION QIsaturation_sp (T, p)  
    IMPLICIT NONE  
    !
    ! Function to return the saturation mr over ice
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-21.8745584*(T-273.15)/(T-7.66))-6.109)
    !
    real(sp), INTENT (IN) :: T, p  
    ! Temperature in Kelvin
    ! Pressure in mb
    real(sp), PARAMETER::tk0c = 273.15, qsa1 = 3.8, qsa2 = -21.8745584, &
         qsa3 = 7.66, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
    if (T > qsa3 .and. p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            > qsa4)then
       Qisaturation_sp = qsa1 / (p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            - qsa4)
    else
       Qisaturation_sp=0.
    end if

  End FUNCTION QIsaturation_sp


  REAL(DP) FUNCTION QIsaturation_dp (T, p)  
    IMPLICIT NONE  
    !
    ! Function to return the saturation mr over ice
    ! Based on tetans formular
    ! QS=3.8/(P*EXP(-21.8745584*(T-273.15)/(T-7.66))-6.109)
    !
    real(dp), INTENT (IN) :: T, p  
    ! Temperature in Kelvin
    ! Pressure in mb
    real(dp), PARAMETER::tk0c = 273.15, qsa1 = 3.8, qsa2 = -21.8745584, &
         qsa3 = 7.66, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
    if (T > qsa3 .and. p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            > qsa4)then
       Qisaturation_dp = qsa1 / (p * EXP (qsa2 * (t - tk0c) / (T - qsa3) ) &
            - qsa4)
    else
       Qisaturation_dp=0.
    end if

  End FUNCTION QIsaturation_dp

  Subroutine salr(z,tsurf,psurf,qsurf,theta_salr)

    Implicit None

    real(wp), intent(in)  :: z(:)          ! height grid (m)
    real(wp), intent(in)  :: tsurf         ! surface temperature (K)
    real(wp), intent(in)  :: psurf         ! surface pressure (Pa)
    real(wp), intent(in)  :: qsurf         ! surface vapour mixing
                                           ! ratio (kg/kg)
    real(wp), intent(out) :: theta_salr(:) ! theta on saturated adiabat

    ! local variables
    integer, parameter :: nfinelevs=1000
    integer :: k, nz, k2
    real(wp) :: dz, qs, tbar, p_mb, l, xi, exner
    real(wp), dimension(0:nfinelevs) :: t, p, q, z_fine

    ! parameters for stratosphere
    real(wp), parameter :: z_pause=16.0e3 ! height tropopause
    real(wp), parameter :: n=0.02         ! bv freq
    real(wp) :: strat_const, n2g

    real(wp), allocatable :: tsalr(:), psalr(:)

    nz=size(z)
    allocate(tsalr(nz))
    allocate(psalr(nz))

    ! calculate refernce profile on fine grid
    dz=(z(nz)+1.)/float(nfinelevs)
    do k=0,nfinelevs
       z_fine(k)=k*dz ! could add zsurf
    end do

    strat_const=g*g/(cp*n*n)
    n2g=n*n/g

    t(0)=tsurf
    p(0)=psurf

    do k=0,nfinelevs-1

       p_mb=0.01*p(k)
       if (t(k) .ge. tk0c) then
          qs=qsaturation(t(k),p_mb)
          l=rlvap
       else
          qs=qsaturation(t(k),p_mb) ! qisaturation(t(k),p_mb)
          l=rlvap                   ! rlsub
       endif

       if (qsurf .le. qs) then
          q(k)=qsurf
          t(k+1)=t(k)-g*dz/cp
       else
          ! constant bv=0.02 in stratosphere
          if (z_fine(k+1) .ge. z_pause) then
             q(k)=0.0
             t(k+1)=exp(n2g*(z_fine(k+1)-z_fine(k)))*(t(k)-strat_const)+strat_const
          else
             q(k)=qs
             t(k+1)=t(k)-g*dz*(1.0+l*qs/(r*t(k))) /       &
                  (cp+l**2.0*qs/(Rw*t(k)**2.0))
          endif
       endif

       tbar=2.0*t(k+1)*t(k)/(t(k)+t(k+1))
       p(k+1)=p(k)*exp(-g*dz/(r*tbar))

    end do
    q(nfinelevs)=0.0

    ! interpolate profile to model grid
    do k=1,nz
       do k2=0,nfinelevs-1
          if (z_fine(k2+1) .gt. z(k)) exit
       end do
       xi=(z(k)-z_fine(k2))/(z_fine(k2+1)-z_fine(k2))
       tsalr(k)=(1.0-xi)*t(k2)+xi*t(k2+1)
       psalr(k)=(1.0-xi)*p(k2)+xi*p(k2+1)
    end do

    ! obtain potential temperature
    do k=1,nz
       exner=(psalr(k)/p0)**r_on_cp
       theta_salr(k)=tsalr(k)/exner
    end do
    theta_salr(1)=theta_salr(2)

    deallocate(psalr)
    deallocate(tsalr)

  end subroutine salr


end Module common_physics
