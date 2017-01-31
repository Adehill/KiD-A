! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing interface to hugh morrison's microphysics code 
!
module mphys_morr_two_moment

  Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, mom_names, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use module_mp_morr_two_moment
  Use diagnostics, only: save_dg, i_dgtime

  Implicit None

  !Logical switches 
  logical :: micro_unset=.True.
  integer:: ih, imom
  character(max_char_len) :: name, units

contains
  
  Subroutine mphys_morrison_interface
    
    real :: t1d(nz), p1d(nz), dz1d(nz),qv1d(nz),qc1d(nz) &
         ,qr1d(nz), qi1d(nz), ni1d(nz), qs1d(nz)         &
         , qg1d(nz), ns1d(nz), nr1d(nz), ng1d(nz)        &
         , w1d(nz)

    real :: wvar1d(nz)

    real :: qc_tend1d(nz), qi_tend1d(nz), qni_tend1d(nz), & 
         qr_tend1d(nz), ni_tend1d(nz), ns_tend1d(nz), & 
         nr_tend1d(nz), t_tend1d(nz), qv_tend1d(nz), &
         qg_tend1d(nz), ng_tend1d(nz)

    real :: precprt1d, snowrt1d

    real :: effc1d(nz), effi1d(nz), effs1d(nz),        &
         effr1d(nz), effg1d(nz)

    real :: qrcu1d(nz), qscu1d(nz), qicu1d(nz)

    real :: qgsten(nz), qrsten(nz), qisten(nz),  &
         qnisten(nz), qcsten(nz)

    ! KiD_2D diag arrays
    real :: qrsten_2d(nz,nx),precprt2d(nx), snowrt2d(nx) 

    integer :: kts, kte, i, j, k
    
    kts=1
    kte=nz
    j=1

    precprt1d = 0.
    snowrt1d  = 0.

    do i=1,nx
       do k=1,nz
       ! zero some of these for safety
          qc_tend1d(k)  = 0.
          qi_tend1d(k)  = 0.
          qni_tend1d(k) = 0.
          qr_tend1d(k)  = 0.
          ni_tend1d(k)  = 0.
          ns_tend1d(k)  = 0.
          nr_tend1d(k)  = 0.
          t_tend1d(k)   = 0.
          qv_tend1d(k)  = 0.
          qg_tend1d(k)  = 0.
          ng_tend1d(k)  = 0.
          effc1d(k)     = 0.
          effi1d(k)     = 0.
          effs1d(k)     = 0.
          effr1d(k)     = 0.
          effg1d(k)     = 0.
          qrcu1d(k)     = 0.
          qscu1d(k)     = 0.
          qicu1d(k)     = 0.
          qgsten(k)     = 0.
          qrsten(k)     = 0.
          qisten(k)     = 0.
          qnisten(k)    = 0.
          qcsten(k)     = 0.

          t1d(k) = (theta(k,i) + (dtheta_adv(k,i)+dtheta_div(k,i))*dt )*exner(k,i)
          p1d(k) = p0*exner(k,i)**(1./r_on_cp)
          dz1d(k) = dz(k)
          qv1d(k) = qv(k,i)+ (dqv_adv(k,i)+dqv_div(k,i))*dt
          qc1d(k) = hydrometeors(k,i,1)%moments(1,1) & 
               + (dhydrometeors_adv(k,i,1)%moments(1,1) &
               + dhydrometeors_div(k,i,1)%moments(1,1))*dt


          if (num_h_moments(2) >= 1) &
               qr1d(k) = hydrometeors(k,i,2)%moments(1,1)& 
               + (dhydrometeors_adv(k,i,2)%moments(1,1) &
               + dhydrometeors_div(k,i,2)%moments(1,1))*dt
          if (num_h_moments(2) >= 2) &
               nr1d(k) = hydrometeors(k,i,2)%moments(1,2)& 
               + (dhydrometeors_adv(k,i,2)%moments(1,2) &
               + dhydrometeors_div(k,i,2)%moments(1,2))*dt
          if (num_h_moments(3) >= 1) &
               qi1d(k) = hydrometeors(k,i,3)%moments(1,1)& 
               + (dhydrometeors_adv(k,i,3)%moments(1,1) &
               + dhydrometeors_div(k,i,3)%moments(1,1))*dt
          if (num_h_moments(3) >= 2) &
               ni1d(k) = hydrometeors(k,i,3)%moments(1,2)& 
               + (dhydrometeors_adv(k,i,3)%moments(1,2) &
               + dhydrometeors_div(k,i,3)%moments(1,2))*dt
          if (num_h_moments(4) >= 1) &
               qs1d(k) = hydrometeors(k,i,4)%moments(1,1)& 
               + (dhydrometeors_adv(k,i,4)%moments(1,1) &
               + dhydrometeors_div(k,i,4)%moments(1,1))*dt
          if (num_h_moments(4) >= 2) &
               ns1d(k) = hydrometeors(k,i,4)%moments(1,2)& 
               + (dhydrometeors_adv(k,i,4)%moments(1,2) &
               + dhydrometeors_div(k,i,4)%moments(1,2))*dt
          if (num_h_moments(5) >= 1) &
               qg1d(k) = hydrometeors(k,i,5)%moments(1,1)& 
               + (dhydrometeors_adv(k,i,5)%moments(1,1) &
               + dhydrometeors_div(k,i,5)%moments(1,1))*dt
          if (num_h_moments(5) >= 2) &
               ng1d(k) = hydrometeors(k,i,5)%moments(1,2)& 
               + (dhydrometeors_adv(k,i,5)%moments(1,2) &
               + dhydrometeors_div(k,i,5)%moments(1,2))*dt
          
          
          wvar1d(k) = 0.5 ! hard-wired not coupled to forcing!
          w1d(k) = w_half(k,i)
       end do

       ! Initialise microphysics 
       if (micro_unset)then
          call morr_two_moment_init
          micro_unset=.False.
       end if


       call morr_two_moment_micro(qc_tend1d, qi_tend1d, qni_tend1d, &
            qr_tend1d, ni_tend1d, ns_tend1d, nr_tend1d,             &
            qc1d, qi1d, qs1d, qr1d, ni1d, ns1d, nr1d,               &
            t_tend1d, qv_tend1d, t1d, qv1d, p1d, dz1d, w1d, wvar1d, &
            precprt1d, snowrt1d,                                    &
            effc1d, effi1d, effs1d, effr1d, dt,                     &
            i,i,j,j,kts,kte,                                        &
            i,i,j,j,kts,kte,                                        &
            qg_tend1d, ng_tend1d, qg1d, ng1d, effg1d,               &
            qrcu1d, qscu1d, qicu1d,                                 &
            qgsten, qrsten, qisten, qnisten, qcsten)
       
       qrsten_2d(:,i)=qrsten(:)
       precprt2d(i) = precprt1d
       snowrt2d(i) = snowrt1d

       
       ! save tendencies
       do k=1,nz
          dtheta_mphys(k,i)=t_tend1d(k)/exner(k,i)
          dqv_mphys(k,i)=qv_tend1d(k)
          dhydrometeors_mphys(k,i,1)%moments(1,1)= qc_tend1d(k)
          dhydrometeors_mphys(k,i,2)%moments(1,1)= qr_tend1d(k)
          dhydrometeors_mphys(k,i,2)%moments(1,2)= nr_tend1d(k)
          dhydrometeors_mphys(k,i,3)%moments(1,1)= qi_tend1d(k)
          dhydrometeors_mphys(k,i,3)%moments(1,2)= ni_tend1d(k)
          dhydrometeors_mphys(k,i,4)%moments(1,1)= qni_tend1d(k)
          dhydrometeors_mphys(k,i,4)%moments(1,2)= ns_tend1d(k)
          dhydrometeors_mphys(k,i,5)%moments(1,1)= qg_tend1d(k)
          dhydrometeors_mphys(k,i,5)%moments(1,2)= ng_tend1d(k)
       end do
       
    end do

    ! Save some diagnostics
    imom=1
    !rain sed
    ih=2
    imom=1
    name=trim(h_names(ih))//'_'//trim(mom_names(imom))//'_sed'
    units=trim(mom_units(imom))//'/s'
    if ( nx == 1 ) then 
       call save_dg(qrsten(:), name, i_dgtime, units, dim='z')
    else
       call save_dg(qrsten_2d(1:nz, 1:nx), name, i_dgtime, units, dim='z')
    endif

    !rain ppt
    ih=2
    name='total_surface_ppt'
    units='kg/kg m'
    if (nx == 1) then
       call save_dg(precprt1d, name, i_dgtime,  units, dim='time')
    else
       call save_dg(precprt2d(1:nx), name, i_dgtime,  units, dim='time')
    endif
    !ice ppt
    ih=3
    name='surface_ppt_for_snow'
    units='kg/kg m'
    if (nx == 1) then
       call save_dg(snowrt1d, name, i_dgtime,  units, dim='time')
    else
       call save_dg(snowrt2d(1:nx), name, i_dgtime,  units, dim='time')
    endif

  end Subroutine mphys_morrison_interface

end module mphys_morr_two_moment
