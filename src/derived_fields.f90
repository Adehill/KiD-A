! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set and calculate derived fields here
!
!

Module derived_fields
  
  Use column_variables
  Use interpolation, only: interpolate, make_wgrid, make_vgrid &
       ,interpolate_x
  Use physconst, only : r_on_cp, g, p0, pi
  Use switches, only : icase

  Implicit none

contains

  subroutine calc_derived_fields
    
    ! Fields derived from model values, but kept 
    ! in column_variables for use in some 
    ! microphysics schemes
    integer :: k,j

    TdegK(:,:)=exner(:,:)*theta(:,:)
    pmb(:,:)=.01*p0*exner(:,:)**(1./r_on_cp)

    ! If we need a staggered grid in vertical...

    if (icase /= 650) then 
       call make_wgrid(z,z_half)
       if ( nx > 1 ) then 
       ! If we need a staggered grid in horizontal
          call make_vgrid(x,x_half)
       else 
          x_half = x
       endif
    else
       !print *, 'make_wgrid called in test_case_2d'
    endif
     
    ! density on half levels
    do k=1,nz-1
       rho_half(k)=sqrt(rho(k)*rho(k+1))
    end do
    rho_half(nz)=rho(nz)*rho(nz)/rho(nz-1)
    
    !dz
    do k=2, nz
       dz(k)=z(k)-z(k-1)
       dz_half(k)=z_half(k)-z_half(k-1)
    end do
    dz(1)=dz(2)
    dz_half(1)=dz_half(2)

    !dx 
    if (nx > 1) then 
       do j=1,nx
          dx(j)=x(j)-x(j-1)
          dx_half(j)=x_half(j)-x_half(j-1)
       end do
       dx(0)=dx(1)
       dx(nx+1)=dx(nx)
       dx_half(0)=dx_half(1)
       dx_half(nx+1)=dx_half(nx)
    else
       dx(:) = 1.
       dx_half(:) = 1.
    endif
       
    
  end subroutine calc_derived_fields


End Module derived_fields
