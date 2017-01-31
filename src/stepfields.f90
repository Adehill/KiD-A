! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Increment the fields
!
!
module stepfields

  Use runtime
  Use typeKind

  Implicit None

  integer:: k, j

 contains 
   
   subroutine step_column
     Use parameters, only : dt, nspecies, &
          num_h_moments, num_h_bins, num_aero_moments, num_aero_bins, nz, &
          nx
     Use physconst, only : Ru => R, p0, r_on_cp
     Use switches

     Use column_variables
     Use pressure_update, only : exner_timestep_upd

     !local variables
     integer :: ih, imom, ibin
     real(wp) :: pos_tmp, srf_tmp(0:nx+1) ! temporary storage
     logical :: sp_test ! test for rounding error
     character(100) :: fmt ! format string
     real :: field(nz,0:nx+1)

     !------------------------------
     ! Manipulate mask if necessary
     !------------------------------
     if (.not. l_periodic_bound .and. nx>1)then
       field_mask(:,0)=0.
       field_mask(:,1)=0.
       field_mask(:,nx)=0.
       field_mask(:,nx+1)=0.
     end if
       
     !-----
     !theta
     !-----
     if (.not. l_fix_theta)then
        srf_tmp(:)=theta(1,:)
        ! Advective part
        if (.not. l_noadv_theta) then
           theta(:,:)=theta(:,:)+dtheta_adv(:,:)*field_mask(:,:)*dt
        endif
        ! mphys part
        if (.not. l_nomphys_theta) then
           theta(:,:)=theta(:,:)+dtheta_mphys(:,:)*field_mask(:,:)*dt
        endif
        ! divergence part
        if (.not. l_noadv_theta) then
           theta(:,:)=theta(:,:)+dtheta_div(:,:)*field_mask(:,:)*dt
        endif
        ! forced part
        theta(:,:)=theta(:,:)+Tforce(:,:)*exner(:,:)*field_mask(:,:)*dt
        ! surface value
        if (isurface_fixed==isurface_fixed)theta(1,:)=srf_tmp(:)
     end if
     !-----
     !qv
     !-----
     if (.not. l_fix_qv)then
        srf_tmp(:)=qv(1,:)
        ! Advective part
        if (.not. l_noadv_qv)then
           if (l_posadv_qv)then
              qv(:,:)=qv(:,:)+max(0.,dqv_adv(:,:))*field_mask(:,:)*dt
           else
              qv(:,:)=qv(:,:)+dqv_adv(:,:)*field_mask(:,:)*dt
           end if
        end if
        ! mphys part
        if (.not. l_nomphys_qv) then
           qv(:,:)=qv(:,:)+dqv_mphys(:,:)*field_mask(:,:)*dt
        endif
        ! divergence part
        qv(:,:)=qv(:,:)+dqv_div(:,:)*field_mask(:,:)*dt
        ! forced part
        qv(:,:)=qv(:,:)+qforce(:,:)*field_mask(:,:)*dt
        ! surface value
        if (isurface_fixed==isurface_fixed)qv(1,:)=srf_tmp(:)

     end if
     !------------
     !aerosol
     !------------
     ! Advective part
     if (.not. l_fix_aerosols)then
        do ih=1,naerosol
           do imom=1,num_aero_moments(ih)
              do ibin=1,num_aero_bins(ih)
                 do j=1,nx
                 do k=1,nz
                    aerosol(k,j,ih)%moments(ibin,imom)=           &
                         aerosol(k,j,ih)%moments(ibin,imom) + (   &
                         daerosol_adv(k,j,ih)%moments(ibin,imom)  &
                         )*field_mask(k,j)*dt
                 end do
                 end do
              end do
           end do
        end do
     end if
     ! mphys part
     do ih=1,naerosol
        do imom=1,num_aero_moments(ih)
           do ibin=1,num_aero_bins(ih)
              do j=1,nx
              do k=1,nz
                 aerosol(k,j,ih)%moments(ibin,imom)=              &
                      aerosol(k,j,ih)%moments(ibin,imom) + (      &
                      + daerosol_mphys(k,j,ih)%moments(ibin,imom) &
                      )*field_mask(k,j)*dt
              end do
              end do
           end do
        end do
     end do
     ! divergence part
     if (.not. l_fix_aerosols)then
        do ih=1,naerosol
           do imom=1,num_aero_moments(ih)
              do ibin=1,num_aero_bins(ih)
                 do j=1,nx
                 do k=1,nz
                    aerosol(k,j,ih)%moments(ibin,imom)=            &
                         aerosol(k,j,ih)%moments(ibin,imom) + (    &
                         + daerosol_div(k,j,ih)%moments(ibin,imom) &
                         )*field_mask(k,j)*dt
                 end do
                 end do
              end do
           end do
        end do
     end if
     !------------
     !hydrometeors
     !------------
     ! Advective part
     if (.not. l_noadv_hydrometeors)then
        do ih=1,nspecies
           do imom=1,num_h_moments(ih)
              do ibin=1,num_h_bins(ih)
                 do j=1,nx
                 do k=1,nz
                    hydrometeors(k,j,ih)%moments(ibin,imom)=           &
                         hydrometeors(k,j,ih)%moments(ibin,imom) + (   &
                         dhydrometeors_adv(k,j,ih)%moments(ibin,imom)  &
                         )*field_mask(k,j)*dt
                 end do
                 end do
              end do
           end do
        end do
     end if
     ! mphys part
     do ih=1,nspecies
        do imom=1,num_h_moments(ih)
           do ibin=1,num_h_bins(ih)
              do j=1,nx
              do k=1,nz
                 hydrometeors(k,j,ih)%moments(ibin,imom)=              &
                      hydrometeors(k,j,ih)%moments(ibin,imom) + (      &
                      + dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                      )*field_mask(k,j)*dt
              end do
              end do
           end do
        end do
     end do
     ! divergence part
     if (.not. l_nodiv_hydrometeors)then
        do ih=1,nspecies
           do imom=1,num_h_moments(ih)
              do ibin=1,num_h_bins(ih)
                 do j=1,nx
                 do k=1,nz
                    hydrometeors(k,j,ih)%moments(ibin,imom)=            &
                         hydrometeors(k,j,ih)%moments(ibin,imom) + (    &
                         + dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                         )*field_mask(k,j)*dt 
                 end do
                 end do
              end do
           end do
        end do
     end if
     !      ! forced part - not yet implemented
     !         do ih=1,nspecies
     !            do imom=1,num_h_moments(ih)
     !               do ibin=1,num_h_bins(ih)
     !                  hydrometeors(k,ih)%moments(ibin,imom)=              &
     !                       hydrometeors(k,ih)%moments(ibin,imom) + (      &
     !                       + dhydrometeors_force(k,ih)%moments(ibin,imom) &
     !                       )*field_mask(k,j)*dt
     !               end do
     !            end do
     !         end do

     ! Positivity check - to prevent rounding errors giving negative
     ! values - !!! NB this violates conservation !!!
     if (l_force_positive)then
     do ih=1,nspecies
        do imom=1,num_h_moments(ih)
           do ibin=1,num_h_bins(ih)
              do j=1,nx
              do k=1,nz
                 pos_tmp=max(0., hydrometeors(k,j,ih)%moments(ibin&
                      &,imom))
                 if (hydrometeors(k,j,ih)%moments(ibin,imom)<0.) then
                    sp_test=1.e3*SPACING(MAX( &
                         dt*dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                         ,dt*dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                         ,dt*dhydrometeors_adv(k,j,ih)%moments(ibin,imom) &
                         ,hydrometeors(k,j,ih)%moments(ibin,imom) &
                         - dt*(                                   &
                         dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                         +dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                         +dhydrometeors_adv(k,j,ih)%moments(ibin,imom)))) &
                         < ABS(hydrometeors(k,j,ih)%moments(ibin,imom))
                    if (sp_test)then
                       fmt='( A, /, A, T33, E11.3, /, &
                            &  A, T33, E11.3, /, &
                            &  A, T18, I2,  /, &
                            &  A, T18, I2,  /, &
                            &  A, T18, I2 )'
                       write(*,fmt)  'Warning: some &
                            &negative numbers have been generated which do&
                            & not look like rounding error.', 'Estimat&
                            &ed rounding error: ', 10.*SPACING(MAX( &
                            dt*dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                            ,dt*dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                            ,dt*dhydrometeors_adv(k,j,ih)%moments(ibin,imom) &
                            ,hydrometeors(k,j,ih)%moments(ibin,imom) &
                            - dt*(                                   &
                            dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                            +dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                            +dhydrometeors_adv(k,j,ih)%moments(ibin&
                            &,imom)))),   &
                            'Negative value of hydrometeor: ', &
                            hydrometeors(k,j,ih)%moments(ibin,imom), &
                            'species: ',ih, &
                            'moment: ',imom, &
                            'bin: ',ibin
                       write(*, '(''Timestep: '', I4, '' z level: '', &
                            &I4)') time_step, k
                       write(*,*) ''
                       write(*,*)   &
                             dt*dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                             ,dt*dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                             ,dt*dhydrometeors_adv(k,j,ih)%moments(ibin,imom) &
                             ,hydrometeors(k,j,ih)%moments(ibin,imom) &
                             - dt*(                                   &
                             dhydrometeors_mphys(k,j,ih)%moments(ibin,imom) &
                             +dhydrometeors_div(k,j,ih)%moments(ibin,imom) &
                             +dhydrometeors_adv(k,j,ih)%moments(ibin&
                             &,imom))
                    end if
                 end if
                 hydrometeors(k,j,ih)%moments(ibin,imom)=pos_tmp
              end do
              end do
           end do
        end do
     end do
     end if
     !----------------
     ! update pressure
     !----------------
     if (l_pupdate) then
        call exner_timestep_upd(theta)
     endif

   end subroutine step_column

 end module stepfields

