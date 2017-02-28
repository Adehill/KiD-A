!
! This Module contains the code required to call 
! Tel-aviv university (TAU) bin microphysics subroutines
! within the 1-D framework
!
!
module mphys_tau_bin

  Use parameters, only : num_h_moments, num_h_bins, &
       nspecies, mom_names, h_names, mom_units, max_char_len, &
       num_aero_moments,num_aero_bins, aero_mom_init
  Use column_variables
  Use physconst, only : p0, this_r_on_cp=>r_on_cp, pi
  Use mphys_tau_bin_declare
  
  Use switches_bin
  Use module_mp_tau_bin 
  Use module_bin_init

 implicit none

  real :: q_lem(JMINP:JMAXP, KKP, NQP)
  real :: th_lem(JMINP:JMAXP, KKP)
  real :: sq_lem(JMINP:JMAXP, KKP, NQP)
  real :: sth_lem(JMINP:JMAXP, KKP)
  real :: w_lem(JMINP:JMAXP, KKP)

  real :: tbase(JMINP:JMAXP, KKP)

  integer i

  !Logical switches 
  logical :: l_ice=.False.
  logical :: micro_unset=.True.

 type qindex
     integer :: ispecies ! Species index
     integer :: imoment  ! moment index
  end type qindex
     
  type(qindex), allocatable :: qindices(:)
  integer :: nqs ! number of qindices (possibly different from nqp)


contains

  subroutine set_micro

    rdt=1./dt

    ! vapour
    
    iq=1
    iqv=iq

    if (num_h_bins(1) >= 11) then
      IMICROBIN=1
      IRAINBIN=1
       
      call qcount(iqss, iq)   ! advected supersat.

!      call qcount(iql, iq)   ! total water for diag
      
      if (num_aero_moments(1) >= 1) then
         do i = 1,ln2
            call qcount(IAERO_BIN(i),iq) ! aerosol bins
         enddo
      end if
      do i = 1,lk
        call qcount(ICDKG_BIN(i),iq) ! cloud mass bins
      enddo
      do i = 1,lk
        call qcount(ICDNC_BIN(i),iq) ! cloud number bins
      enddo
      
      nqs=aero_bin+ln2+lk+lk+1

      allocate(qindices(nqs))      

      ! Set qindices for later use
      if (num_aero_moments(1) >= 1) then
         do i = 1,ln2
            qindices(IAERO_BIN(i))%ispecies=1
            qindices(IAERO_BIN(i))%imoment=1
         enddo
      end if
      if (num_h_bins(1) >= 4)then ! cloud mass
         do i = 1,lk
            qindices(ICDKG_BIN(i))%ispecies=1
            qindices(ICDKG_BIN(i))%imoment=1
            qindices(ICDNC_BIN(i))%ispecies=1
            qindices(ICDNC_BIN(i))%imoment=2
         enddo
      end if

    end if ! bin model selected
  
  end subroutine set_micro


 
  subroutine qcount(var, count)
    integer, intent(out) :: var
    integer, intent(inout) ::  count

    count=count+1
    var=count
  end subroutine qcount

  subroutine mphys_tau_bin_interface

    ! Warm bin scheme
    integer ::j, k, iq

    ! base and top of the z grid
    integer :: kts, kte

    kts = 1
    kte = kkp

    ! Set up input arrays...
    dzn(:)=dz_half(:)

    rhon(:)=rho(:)
    rdz_on_rhon(:)=1./(dz(:)*rhon(:))

    ! Set up microphysics species
    if (micro_unset)then ! See later call to microset
       call set_micro
    end if


    do j=jminp,jmaxp
       q_lem (j, :, iqv) = qv(:, j)
       q_lem (j, :, iqss) = ss(:, j)
       rprefrcp(j,:)=exner(:,j)
       prefrcp(j,:)=1./rprefrcp(j,:)
       prefn(j,:)=p0*exner(:,j)**(1./this_r_on_cp)
       ! Reference temperature (this is fixed in lem, but 
       ! shouldn't make a difference for microphysics if we 
       ! just set it to be the current profile (i.e. th'=0)
       tbase(j,:)=theta(:,j)*exner(:,j)

       do iq=1,ln2
         ih=qindices(IAERO_BIN(iq))%ispecies
         imom=qindices(IAERO_BIN(iq))%imoment
         do k = kts,kte
            q_lem (j, k, IAERO_BIN(iq)) = aerosol(k,j,ih)%moments(iq,imom)
         enddo
      enddo
       do iq=1,lk
         ! mass bins
         ih=qindices(ICDKG_BIN(iq))%ispecies
         imom=qindices(ICDKG_BIN(iq))%imoment
         do k = kts,kte
            q_lem (j, k, ICDKG_BIN(iq)) = hydrometeors(k,j,ih)%moments(iq,imom)
         enddo
         ! number bins
         ih=qindices(ICDNC_BIN(iq))%ispecies
         imom=qindices(ICDNC_BIN(iq))%imoment
         do k = kts,kte
            q_lem (j, k, ICDNC_BIN(iq)) = hydrometeors(k,j,ih)%moments(iq,imom)
         enddo
      enddo
       th_lem (j, :) = 0.0
       w_lem(j,:)=w_half(:,j)
    end do

    if (micro_unset)then
       tref(:) =  tbase(1,:)
       call bin_init !initialises the cloud bin categories
       call data     !reads in and sets the coll-coal kernal

       DO IQ = 1,LN2
         ih=qindices(IAERO_BIN(iq))%ispecies
         imom=qindices(IAERO_BIN(iq))%imoment
         DO J = JMINP , JMAXP
            do k = kts,kte
               CCNORIG(J,k,IQ) = aerosol(k,j,ih)%moments(iq,imom)
            enddo
          ENDDO
       ENDDO

       DO K = 1, KKP
         DO J = JMINP, JMAXP
           TOTCCNORIG(J,K) = 0.0
           DO IQ = 1, LN2
             TOTCCNORIG(J,K) = TOTCCNORIG(J,K) + CCNORIG(J,K,IQ)
           ENDDO
         ENDDO
       ENDDO
       DO IQ = 1, Ln2
          CCNORIGTOT(IQ) = 0.0
          CCNORIGAVG(IQ) = 0.0
           DO K = 1, KKP
             DO J = JMINP, JMAXP
               CCNORIGTOT(IQ) = CCNORIGTOT(IQ) + CCNORIG(J,K,IQ)
             ENDDO
           ENDDO
           CCNORIGAVG(IQ) = CCNORIGTOT(IQ)/(JJP*KKP)
       ENDDO    
        micro_unset=.False.
    end if

    ! This bit doesn't yet use the switches on advection...
    do j=jminp,jmaxp
       sth_lem(j,:)=dtheta_adv(:,j)+dtheta_div(:,j)
       sq_lem(j,:,iqv)=dqv_adv(:,j)+dqv_div(:,j)

       sq_lem(j,:,iqss)=dss_adv(:,j)+dss_div(:,j)
       
       do iq=1,ln2
          ih=qindices(iaero_bin(iq))%ispecies
          imom=qindices(iaero_bin(iq))%imoment
          do k=1,nz
             sq_lem(j,k,iaero_bin(iq))=(daerosol_adv(k,j,ih)%moments(iq,imom) &
                  + daerosol_div(k,j,ih)%moments(iq,imom))
          end do
       enddo
       do iq=1,lk
          ih=qindices(icdkg_bin(iq))%ispecies
          imom=qindices(icdkg_bin(iq))%imoment
          do k=1,nz
             sq_lem(j,k,icdkg_bin(iq))=dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                  + dhydrometeors_div(k,j,ih)%moments(iq,imom)
          end do
          ih=qindices(icdnc_bin(iq))%ispecies
          imom=qindices(icdnc_bin(iq))%imoment
          do k=1,nz
             sq_lem(j,k,icdnc_bin(iq))=dhydrometeors_adv(k,j,ih)%moments(iq,imom) &
                  + dhydrometeors_div(k,j,ih)%moments(iq,imom)
          end do
       end do
     end do
! test if the transport has moved mass and number around
     DO K = 1, nz
        DO J = jminp,jmaxp
           DO IQ = 1, LK
              CALL ADVECTcheck(J,K,iq,DT,Q_lem(J,K,ICDKG_BIN(iq)),              &
    &                 Q_lem(J,K,ICDNC_BIN(iq)),SQ_lem(J,K,ICDKG_BIN(iq)),  &
    &                 SQ_lem(J,K,ICDNC_BIN(iq)))
           ENDDO
        ENDDO
     ENDDO

     call tau_bin(1,kts, kte, tbase, q_lem, sth_lem, sq_lem, dt, rdt )


     do j=jminp,jmaxp
        sth_lem(j,:)=sth_lem(j,:)-(dtheta_adv(:,j)+dtheta_div(:,j))
        sq_lem(j,:,iqv)=sq_lem(j,:,iqv)-(dqv_adv(:,j)+dqv_div(:,j))

        do iq=1,ln2
           ih=qindices(iaero_bin(iq))%ispecies
           imom=qindices(iaero_bin(iq))%imoment
           do k=1,nz
              sq_lem(j,k,iaero_bin(iq))=sq_lem(j,k,iaero_bin(iq))       &
                   - (daerosol_adv(k,j,ih)%moments(iq,imom)                &
                   + daerosol_div(k,j,ih)%moments(iq,imom))
           end do
        enddo
        
        do iq=1,lk
           ih=qindices(icdkg_bin(iq))%ispecies
           imom=qindices(icdkg_bin(iq))%imoment
           do k=1,nz
              sq_lem(j,k,icdkg_bin(iq))= sq_lem(j,k,icdkg_bin(iq))      &
                   - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                   + dhydrometeors_div(k,j,ih)%moments(iq,imom))
           end do
           ih=qindices(icdnc_bin(iq))%ispecies
           imom=qindices(icdnc_bin(iq))%imoment
           do k=1,nz
              sq_lem(j,k,icdnc_bin(iq))= sq_lem(j,k,icdnc_bin(iq))      &
                   - (dhydrometeors_adv(k,j,ih)%moments(iq,imom)           &
                   + dhydrometeors_div(k,j,ih)%moments(iq,imom))
           end do
        end do
        
        ! For now set no microphysics on the bottom level - this would be 
       ! better done by having a subterranian level 0 in column variables
        sq_lem(j,1,1)=0
        
        dtheta_mphys(:,j)=sth_lem(j,:)

        dqv_mphys(:,j)=sq_lem(j,:,iqv)
        
       !
       ! update supersaturation field here (not in step fields)
       !
        ss(:,j) = q_lem(j,:,iqss)
    
        do iq=1,ln2
           ih=qindices(iaero_bin(iq))%ispecies
           imom=qindices(iaero_bin(iq))%imoment
           do k=1,nz
              daerosol_mphys(k,j,ih)%moments(iq,imom) =              &
                   sq_lem(j,k,iaero_bin(iq))
           end do
        enddo

        do iq=1,lk
           ih=qindices(icdkg_bin(iq))%ispecies
           imom=qindices(icdkg_bin(iq))%imoment
           do k=1,nz
              dhydrometeors_mphys(k,j,ih)%moments(iq,imom) =         &
                   sq_lem(j,k,icdkg_bin(iq))
           end do
           ih=qindices(icdnc_bin(iq))%ispecies
           imom=qindices(icdnc_bin(iq))%imoment
           do k=1,nz
              dhydrometeors_mphys(k,j,ih)%moments(iq,imom) =        &
                   sq_lem(j,k,icdnc_bin(iq))
           end do
        end do
     end do
     
   end subroutine mphys_tau_bin_interface

!DECK ADVcheck
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
      SUBROUTINE  ADVECTcheck(J,K,IQ,DT,ZQmass,ZQnum,Sourcemass,&
           Sourcenum)
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS      
      IMPLICIT NONE
      
!CALL PRAMETR
!CALL RI
!CALL XD      
!CALL GRID1     
!local variable
      REAL :: QmassFLD, ZQmass, Qmass
      REAL :: Qnumfield, ZQnum, Qnum
      REAL :: Sourcemass
      REAL :: Sourcenum
      REAL :: AVG, AVGinit
      REAL :: DT,RDT
!loop counters
      INTEGER J, K, IQ      
          
      RDT = 1.0/DT

!First calculate the new field      
      QmassFLD=(ZQmass+(DT*Sourcemass))
      Qnumfield = (ZQnum+(DT*Sourcenum))
!Change units to microphys units (just for consistency
 
      QmassFLD = (QmassFLD*rhon(k))/1.e3
      Qnumfield = (Qnumfield*rhon(k))/1.e6
      Sourcemass = (Sourcemass*rhon(k))/1.e3
      Sourcenum = (Sourcenum*rhon(k))/1.e6

!calculate the average particle size for the bin

      IF(Qnumfield  >  0.0)THEN
        AVG =QmassFLD/Qnumfield
        IF(AVG >  (2.*X_BIN(IQ))) THEN
          sourcenum=((QmassFLD/(2.*X_BIN(IQ)-1.e-20))-                        &
     &                             ((ZQnum*rhon(k))/1.e6))*RDT
        ENDIF
        IF(AVG <  X_BIN(IQ).AND.AVG >  0.0) THEN
          sourcenum=((QmassFLD/(X_BIN(IQ)+1.e-20))-                           &
     &                             ((ZQnum*rhon(k))/1.e6))*RDT
        ENDIF

      ENDIF
!
      Sourcemass = (Sourcemass*1.e3)/rhon(k)
      Sourcenum  = (Sourcenum*1.e6)/rhon(k)
      
!do positivity check after normalisation as changing SQ by normalisation 
!can lead to negatives      
      IF((ZQmass+(DT*Sourcemass) <  0.0).or.                            &
     &                    (ZQnum+(DT*Sourcenum) <  0.0)) THEN
        Sourcemass = -0.99999 * ZQmass/DT
        Sourcenum =  -0.99999 * ZQnum/DT
      ENDIF

      IF (ABS(Sourcemass) <  1.e-20.OR. ABS(Sourcenum) <  1.e-20) THEN
          Sourcemass = 0.0
          Sourcenum = 0.0
      ENDIF

      END subroutine ADVECTcheck



end module
