! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to calculate divergence/convergence terms
!

module divergence
  
  Use parameters, only : num_h_moments, num_h_bins, nz, nspecies &
     ,num_aero_moments,num_aero_bins, aero_mom_init
  Use column_variables
  Use diagnostics, only : dz_half
  Use switches, only : isurface, isurface_fixed
  

  Implicit None

contains

  subroutine diverge_column

    real(wp) :: dudx(nz), dwrhodz(nz)
    integer :: ih, imom, ibin, k, j
    
    do k=1,nz-1
       dwrhodz(k)=(w_half(k+1,1)*rho_half(k+1)-w_half(k,1)*rho_half(k))/dz_half(k)
    end do
    dwrhodz(nz)=0.
    dudx(:)=-dwrhodz(:)/rho(:)

    dudx(1)=0 ! This should be taken care of in the surface options

    do j = 1, nx
       !theta
       dtheta_div(:,j)=-theta(:,j)*dudx(:)
       !qv
       dqv_div(:,j)=-qv(:,j)*dudx(:)
       !ss
       dss_div(:,j)=-ss(:,j)*dudx(:)
    end do
       
    !aerosol
    do ih=1,naerosol
       do ibin=1,num_aero_bins(ih)
          do imom=1,num_aero_moments(ih)
             do j = 1,nx
                do k=1,nz
                   daerosol_div(k,j,ih)%moments(ibin,imom)= &
                        -aerosol(k,j,ih)%moments(ibin,imom)*dudx(k)
                end do
             end do
          end do
       end do
    end do
    !hydrometeors
    do ih=1,nspecies
       do ibin=1,num_h_bins(ih)
          do imom=1,num_h_moments(ih)
             do j = 1,nx
                do k=1,nz
                   dhydrometeors_div(k,j,ih)%moments(ibin,imom)= &
                        -hydrometeors(k,j,ih)%moments(ibin,imom)*dudx(k)
                end do
             end do
          end do
       end do
    end do
    

  end subroutine diverge_column

end module divergence
