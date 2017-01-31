! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing (slightly modified) ultf2d subroutine.
!
!
module ultf2d_mod

  Use typeKind
  Use parameters
  Use runtime
  Use switches
  Use column_variables, only: dx,field_mask, x_half, x

  ! Parameters needed by ultf2d
  integer, parameter, private :: &
       &  KKP=nz  &
       & , NPROC=1 &
       & , NPES=1
  ! JJP, JMINP and JMAXP set dynamically depending on 
  ! 1D or 2D run
  integer, private :: JJP, JMINP, JMAXP, JCOL
           
  REAL, private :: CX=100,CY=100
  INTEGER, private :: INFO, IRIGHT, ILEFT

  integer, private :: kmin=1, kmax=nz, kdzof=0  
  real,private :: adz(KKP), ardz(KKP)
  real,private :: adzn(KKP), ardzn(KKP)
  real,private :: dt_ultf2d

  real, allocatable :: &
        fflxl(:,:), fflxb(:,:), w_ultf2d(:,:),v_ultf2d(:,:) &
        , field_ultf2d(:,:)

  character(max_char_len) :: name, units

 contains

  subroutine ultf2d_interface(field, v, w, z, zw, rho, rhow, &
       field_adv, wf_surf)

    real(wp), intent(in) :: field(:,0:)
    real(wp), intent(out):: field_adv(:,0:)
    real(wp), intent(in) :: w(:,0:), v(:,0:), z(:), zw(:), rho(:), rhow(:)    

    real(wp), intent(in), optional :: wf_surf(:)
    ! local variables
    integer :: k, j
    real(wp) :: dx_arb=1., dwrhodz(KKP), dudx(KKP)

    do k=2,nz
       adz(k)=zw(k)-zw(k-1)
       adzn(k)=z(k)-z(k-1)
       ardz(k)=1./adz(k)
       ardzn(k)=1./adzn(k)
    end do

    do k=2,KKP-1
       dwrhodz(k)=(w(k+1,1)*rhow(k+1)-w(k,1)*rhow(k))/adz(k)
    end do
    dwrhodz(KKP)=0.

    dudx=0. ! If l_diverge, then use the separate 
            ! divergence calculation, else if l_diverge_advection
            ! then assume a divergent field in the advection
            ! calculation.
    if (.not. l_diverge .and. l_diverge_advection)  dudx(:)=-dwrhodz(:)/rho(:)

    ! set up the jjp (y-dirn) grid depending on 1D or 2D
    if (nx == 1) then 
       JJP = 3
       JMINP = JJP-1
       JMAXP = JJP+1
       JCOL = JJP-2
    else
       JJP = nx
       JMINP = 0
       JMAXP = JJP+1
       ! in 2D JCOL not required
       JCOL = JJP-2
    endif
    ! allocate the fluxes, velocities and prognostic field used in ULTIMATE
    allocate(fflxl(0:JJP+1,KKP))
    allocate(fflxb(0:JJP+1,KKP))
    allocate(w_ultf2d(0:JJP+1,KKP))
    allocate(v_ultf2d(0:JJP+1,KKP))
    allocate(field_ultf2d(0:JJP+1,KKP))

    ! initialise variables to 0
    fflxl=0.0
    fflxb=0.0
    w_ultf2d=0.0
    v_ultf2d=0.0
    field_ultf2d=0.0

    ! need to flip arrays for advection
    if (nx == 1 ) then
       do j = 0,jjp+1
          w_ultf2d(j,:) = w(:,1)
          v_ultf2d(j,:) = v(:,1)
          field_ultf2d(j,:) = field(:,1)
       end do
       v_ultf2d(1,:)=dx_arb*dudx(:)
       v_ultf2d(3,:)=-dx_arb*dudx(:)
    else
       do k = 1,kkp 
          do j = 1,jjp
             w_ultf2d(j,k) = w(k,j)
             v_ultf2d(j,k) = v(k,j)
             !if (j == 1) then
             !   print *, v_ultf2d(j,k), j,k
             !   print *, v_ultf2d(0,k), '0',k
             !endif
             field_ultf2d(j,k) = field(k,j)
          end do
       end do
       if (l_periodic_bound) then
            ! wrap winds for periodic boundaries
          ! for some reason field and wind in 0
          ! column do not pass properly 
          w_ultf2d(0,:) = w(:,nx)
          v_ultf2d(0,:) = v(:,nx)
          field_ultf2d(0,:) = field(:,nx)

          w_ultf2d(jjp+1,:) = w(:,1)
          v_ultf2d(jjp+1,:) = v(:,1)
          field_ultf2d(jjp+1,:) = field(:,1)
       endif
    endif
    
    dt_ultf2d=dt

    call ultf2d(dt, field_ultf2d, v_ultf2d, w_ultf2d, fflxl, fflxb &
         , ardz, ardzn, adzn ,kmin,kmax,kdzof                      & 
         , v_ultf2d(nx,:), field_ultf2d(nx,:), JJP)


    field_adv=0.

    ! field_adv calc is same as LEM but w source is first term
    !

    if (nx == 1) then 
       do k=2,size(field(:,1))-1
          field_adv(k,jcol)=(ardz(k)/rho(k))*(w(k-1,jcol)*FFLXB(jcol,K)*rhow(k-1) &
               - W(K,jcol)*FFLXB(jcol,K+1)*rhow(k)) &
               + (( v_ultf2d(jcol-1,k)*fflxl(jcol,k) - &
               v_ultf2d(jcol,k)*fflxl(jcol+1,k)))/dx(jcol)
       end do
    else
       do j=1,jjp
          !print *, 'flux', x_half(j), x(j) , fflxl(j,18), w(k-1, j), fflxb(j,18)
          do k=2,size(field(:,1))-1
             field_adv(k,j)=(ardz(k)/rho(k))*(w(k-1,j)*FFLXB(j,K)*rhow(k-1) &
                  - W(K,j)*FFLXB(j,K+1)*rhow(k)) &
                  + (( v_ultf2d(j-1,k)*fflxl(j,k) - &
                  v_ultf2d(j,k)*fflxl(j+1,k)))/dx(j)
          enddo       
          ! deal with top level k = kkp when nx > 1
          ! if 1-D no fluxes into the top level
          k = kkp
          field_adv(k,j) = (ardz(k)/rho(k))*(w(k-1,j)*FFLXB(j,K)*rhow(k-1)) &
               + ( v_ultf2d(j-1,k)*fflxl(j,k) - &
               v_ultf2d(j,k)*fflxl(j+1,k))/dx(j)
             
          field_adv(1,j)=-(ardz(1)/rho(k))* &
               (W(1,j)*FFLXB(j,2)*rhow(k))  &
               + ( v_ultf2d(j-1,1)*fflxl(j,1) - &
               v_ultf2d(j,1)*fflxl(j+1,1))/dx(j)  
       end do
    endif

    if (l_periodic_bound) then 
       field_adv(:,0) = field_adv(:,nx)
       field_adv(:,nx+1) = field_adv(:,1)
    endif

    ! Mask those points we want to keep fixed 
    ! (e.g. subterrainean in orographic cases)
    field_adv=field_adv*field_mask

    k=1
    select case(isurface)
    case(isurface_fixed)
       field_adv(k,:)=0
    case(isurface_flux)
       do j = 1,jjp
          if (present(wf_surf))then
             field_adv(k,j)=(ardz(k)/rho(k))*(wf_surf(j)*rhow(1) &
                  &       - W(K,j)*FFLXB(j,K+1)*rhow(k))
          else
             field_adv(k,j)=-(ardz(k)/rho(k))* &
                  &        W(K,j)*FFLXB(j,K+1)*rhow(k)
          endif
       enddo
    end select

    deallocate(fflxl)
    deallocate(fflxb)
    deallocate(w_ultf2d)
    deallocate(v_ultf2d)
    deallocate(field_ultf2d)

  end subroutine ultf2d_interface

      SUBROUTINE ULTF2D(DT,ZF,V,W,FFLXL,FFLXB,ARDZ,ARDZN,ADZN, &
     & KMIN,KMAX,KDZOF,VHALO,ZFHALO, JJP)

!..CALCULATES THE EFFECTIVE FACE VALUES FOR ADVECTION USING                        
!..LEONARD'S ULTIMATE QUICKEST SCHEME WITH FIRST MULTI-DIM LIMITER                 
!..AND GRADT TERMS ADDED TO QUICKEST SCHEME.                                       
!..DEALS WITH X-INDEPENDENT 2D CASE                                                
!..NOTE: ALL FACE COURANT NUMBERS (EG. VFACE) HAVE FACTOR OF 0.5                   
!..INCLUDED. REMEMBER THIS ROUTINE CALCULATES THE RIGHT/TOP FACE FLUX              
!..ALTHOUGH IT SAVES IT IN FFLXL(J+1), ETC.                                        
!------------------------------------------------------------------                
       IMPLICIT NONE                                                               
                                                                                   
!*CALL PRAMETR                                                                      
!*CALL GRID1                                                                        
!*CALL ME                                                                           
                                                                                   
!                                                                                  
! Subroutine arguments                                                             
!                                                                                  
! Intent (IN)                                                                      
!
      INTEGER, INTENT(IN) ::  &
     &    JJP                    ! From KiD model depending on 1-D or 2-D  
      INTEGER, INTENT(IN) ::  &                                                     
     &    KMIN,KMAX           &   ! bottom and top grid-lev for advection           
     &   ,KDZOF                  ! =1 for advection of W, 0 otherwise,             
!                                !       to shift vertical grid index              
      REAL, INTENT(IN) ::     &                                                     
     &    DT                     ! timestep (s)                                    
      REAL, INTENT(IN), DIMENSION(0:JJP+1,KKP) ::     &                             
     &    ZF                  &   ! advected field                                  
     &   , V                  &   ! y comp.t of velocity                            
     &   , W                     ! z comp.t of velocity                            
      REAL, INTENT(IN), DIMENSION(KKP) ::             &                             
     &    ARDZ                &   ! vertical grid arrays, equal to                  
     &   ,ARDZN               &   ! standard model variables (without `A')          
     &   ,ADZN                   ! for all fields except W for which               
!                                ! they are staggered.                             
      REAL, INTENT(IN), DIMENSION(KKP) ::  &                                        
     &    ZFHALO,  & ! J-1 values for ZF for multi PE case                          
     &     VHALO     ! J-1 values for ZV for multi PE case                          
!                                                                                  
! Intent (OUT)                                                                     
!                                                                                  
      REAL,  INTENT(OUT), DIMENSION(0:JJP+1,KKP) ::     &                           
     &    FFLXB       &  ! flux through bottom cell face                            
     &   ,FFLXL         ! flux through left (y-dirn) cell face                     
!                                                                                  
! Local variables:                                                                 
!                                                                                  
      INTEGER         &                                                             
     &   K,J          &  ! loop counters                                            
     &  ,JSU          &  ! index of upwinded node in y-direction                    
     &  ,JOFSET,KOFSET & ! 0/1 off-set for J,K-index, for upwinding                 
     &  ,JPERIOD       & ! wrap-around for J-index                                  
     &  ,KP0,KP1       & ! K and K+1, but restricted to KKP (to                     
!                       ! extrapolate gradient at upper boundary)                  
     &  ,KPOSWF        & ! 0/1 off-set for K-index, for upwinding                   
     &  ,KSU,KSL        ! upper and lower vertical grid indexes for                
!                       ! GRADT term for horizonal advection to                    
!                       ! incorporate upper b.c.                                   
                                                                                   
      REAL :: &                                                                        
     &   RDZU,RDZC,RDZD &! upwinded 1/(grid size)                                   
     &  ,DZD            &! downwind grid spacing                                    
     &  ,FXU,FXC,FXD    &! upwinded nodal points in x-direction                     
     &  ,FYU,FYC,FYD    &! upwinded nodal points in y-direction                     
     &  ,FZU,FZC,FZD    &! upwinded nodal points in z-direction                     
     &  ,SUM_CFL_OUT    &! sum of absolute value of out-flowing CFL nos             
     &  ,FGT1,FGT2      &! `GRADT' terms                                            
     &  ,R6             &! 1/6                                                      
     &  ,VFACE,WFACE    &! advecting velocity on cell face                          
     &  ,ADVPOS,ADVNEG  &! sign indicators (0 or 1)                                 
     &  ,YCFLMOD,ZCFLMOD &! absolute value of face CFL number                        
     &  ,HRCFL          &! 1/SUM_CFL_OUT                                            
     &  ,HCFL           &! 0.5*CFLMOD                                               
     &  ,CFLCOEF        &! coefficient in QUICKEST calculation                      
     &  ,SDIFFF         &! 0/1 sign indicator for non-monotonicity                  
     &  ,SIGNV,SIGNW    ! 0.5*(sign of V and W)                                    
      ! AJEQ0 no longer required at V2.3                                           
!--------------------------------------------------------------                    
! scalars used in calculation of flux-limited cell face value                      
!  (same set for X and Y directions)                                               
      REAL ::           &                                                           
     &  FZDELS        & ! downwind-upwind difference of F                           
     &  ,RFZDELS      &  ! 1/FZDELS                                                 
     &  ,FZCURVS      &  ! curvature of F                                           
     &  ,FZFS         &  ! 1D QUICKEST face-value of F - (1st GRADT term)           
     &  ,FZNFS        &  ! normalised face value                                    
     &  ,FZNCS        &  ! normalised central node value                            
     &  ,FZNREFS      &  ! normalised reference value                               
     &  ,FZTEMP       &  ! limited face value                                       
     &  ,FYDELS,RFYDELS,FYCURVS,FYFS,FYNFS,FYNCS,FYNREFS,FYTEMP                    
! Arrays to pass the haloes for the wrapping of FFLXL at end of routine            
      REAL, DIMENSION(KKP) :: HALOSEND, HALORECV                                   
                                                                                   
!..VERTICAL FLUX, SET UP STENCIL                                                   
      R6=1.0/6.0                                                                   
!                                                                                  
!.....VERTICAL FLUXES                                                              
!                                                                                  
      DO 100 J=1,JJP                                                               
!---FIRST CASE: K = KMIN = 1 ---------------------------                           
! THIS SECTION ONLY FOR VERTICAL ADVECTION OF W                                    
        IF(KMIN.EQ.1)THEN                                                          
          K=1                                                                      
         IF(W(J,K).GE.0.0)THEN                                                     
           FZU = ZF(J,1)       ! K-1 MAPPED ONTO 1                                 
           FZC = ZF(J,K)                                                           
           FZD = ZF(J,K+1)                                                         
         ELSE                                                                      
           FZU = ZF(J,K+2)                                                         
           FZC = ZF(J,K+1)                                                         
           FZD = ZF(J,K)                                                           
         ENDIF                                                                     
         FZDELS = FZD - FZU                                                        
         FZCURVS = FZD - 2.0*FZC + FZU                                             
         IF( ABS(FZCURVS) .GE. ABS(FZDELS) )THEN                                   
           FFLXB(J,K+1) = FZC                                                      
         ELSE                                                                      
          IF(W(J,K).GE.0.0)THEN                                                    
           SUM_CFL_OUT=DT*( W(J,K)*ARDZ(K+KDZOF)+ &                                 
     &      CY*(MAX(0.,V(J,K))+ABS(MIN(0.,V(J-1,K)))) )                            
!..1ST GRADT TERM, IN Y DIRN                                                       
           VFACE=0.125*(V(J-1,K)+V(J,K)+V(J-1,K+1)+V(J,K+1))                       
           JSU=J-NINT(SIGN(1.0,VFACE))                                             
           FGT1=ABS(VFACE)*CY*DT*(FZC-ZF(JSU,K))                                   
!..2ND GRADT TERM, IN X DIRN, IS ZERO                                              
           RDZU=ARDZN(2)  ! K MAPPED ONTO 2                                        
           RDZC=ARDZ(K+KDZOF)                                                      
         ELSE                                                                      
           SUM_CFL_OUT=DT*(                        &                                
     &       ARDZ(K+1+KDZOF)*(MAX(0.,W(J,K+1))-W(J,K))+  &                          
     &       CY*(MAX(0.,V(J,K+1))+ABS(MIN(0.,V(J-1,K+1)))) )                       
!..1ST GRADT TERM, IN Y DIRN                                                       
           VFACE=0.125*(V(J-1,K)+V(J,K)+V(J-1,K+1)+V(J,K+1))                       
           JSU=J-NINT(SIGN(1.0,VFACE))                                             
           FGT1=ABS(VFACE)*CY*DT*(FZC-ZF(JSU,K+1))                                 
!..2ND GRADT TERM, IN X DIRN, IS ZERO                                              
           RDZU=ARDZN(K+2)                                                         
           RDZC=ARDZ(K+1+KDZOF)                                                    
         ENDIF                                                                     
!..CALCULATE FLUXES                                                                
         DZD=ADZN(K+1)                                                             
         RDZD=ARDZN(K+1)                                                           
         ZCFLMOD=ABS(W(J,K)*ARDZN(K+1)*DT)                                         
         HRCFL=1./(SUM_CFL_OUT+epsilon(1.d0))                                             
         HCFL=0.5*ZCFLMOD                                                          
         CFLCOEF=R6*(1.0-ZCFLMOD*ZCFLMOD)*(DZD*DZD*RDZC)                           
!
         if (FZDELS < 0.0) then
            RFZDELS=1.0/(FZDELS-sqrt(epsilon(1.0_wp)))
         else
            RFZDELS=1.0/(FZDELS+sqrt(epsilon(1.0_wp)))
         endif
         FZFS=0.5*(FZD+FZC)-HCFL*(FZD-FZC)-         &                               
     &        CFLCOEF*((FZD-FZC)*RDZD-(FZC-FZU)*RDZU) - FGT1                          
!                           !   GRADT TERM, FGT1, ADDED BEFORE LIMITING            
         FZNFS=(FZFS-FZU)*RFZDELS                                                  
         FZNCS=(FZC-FZU)*RFZDELS                                                   
         FZNREFS=FZNCS*HRCFL                                                       
         FZTEMP=MAX(FZNCS,MIN(FZNFS,FZNREFS,1.0))                                  
         FFLXB(J,K+1) = ( FZTEMP*FZDELS + FZU )                                    
         ENDIF    ! first-order upwind                                             
!                                                                                  
       ENDIF  ! KMIN EQ 1                                                          
!---SECOND CASE: K = KKP-1 (BOTTOM OF KKP'TH BOX)-------------------               
       K=KKP-1                                                                     
       IF(W(J,K).GE.0.0)THEN                                                       
         FZU = ZF(J,K-1)                                                           
         FZC = ZF(J,K)                                                             
         FZD = ZF(J,K+1)                                                           
       ELSE                                                                        
         FZU = 2.0*ZF(J,K+1)-ZF(J,K)                                               
         FZC = ZF(J,K+1) !top boundary condition corrected at V2.1                 
         FZD = ZF(J,K)                                                             
       ENDIF                                                                       
       FZDELS = FZD - FZU                                                          
       FZCURVS = FZD - 2.0*FZC + FZU                                               
       IF( ABS(FZCURVS) .GE. ABS(FZDELS) )THEN                                     
         FFLXB(J,K+1) = FZC                                                        
       ELSE                                                                        
       IF(W(J,K).GE.0.0)THEN                                                       
           SUM_CFL_OUT=DT*(                  &                                      
     &      ARDZ(K+KDZOF)*(W(J,K)+ABS(MIN(0.,W(J,K-1))))+   &
     &           CY*(MAX(0.,V(J,K))+ABS(MIN(0.,V(J-1,K)))) )                       
!..1ST GRADT TERM, IN Y DIRN                                                       
           VFACE=0.125*(V(J-1,K)+V(J,K)+V(J-1,K+1)+V(J,K+1))                       
           JSU=J-NINT(SIGN(1.0,VFACE))                                             
           FGT1=ABS(VFACE)*CY*DT*(FZC-ZF(JSU,K))                                   
!..2ND GRADT TERM, IN X DIRN, IS ZERO                                              
           RDZU=ARDZN(K)                                                           
           RDZC=ARDZ(K+KDZOF)                                                      
       ELSE                                                                        
           SUM_CFL_OUT=DT*(                                  &                     
     &       -W(J,K)*ARDZ(MIN(KKP,K+1+KDZOF))+               &                      
     &       CY*(MAX(0.,V(J,K+1))+ABS(MIN(0.,V(J-1,K+1)))) )                       
!_Extrapolate linearly to level equidistant above upper boundary                   
           FZD=ZF(J,K)   !                                                         
!..1ST GRADT TERM, IN Y DIRN                                                       
           VFACE=0.125*(V(J-1,K)+V(J,K)+V(J-1,K+1)+V(J,K+1))                       
           JSU=J-NINT(SIGN(1.0,VFACE))                                             
           FGT1=ABS(VFACE)*CY*DT*(FZC-ZF(JSU,K+1))                                 
!..2ND GRADT TERM, IN X DIRN, IS ZERO                                              
           RDZU=ARDZN(KKP)        !K+2 MAPPED TO KKP                               
           RDZC=ARDZ(MIN(KKP,K+1+KDZOF))                                           
       ENDIF                                                                       
!..CALCULATE FLUXES...................                                             
       DZD=ADZN(K+1)                                                             
       RDZD=ARDZN(K+1)                                                             
       ZCFLMOD=ABS(W(J,K)*ARDZN(K+1)*DT)                                           
       HRCFL=1./(SUM_CFL_OUT+epsilon(1.d0))                                               
       HCFL=0.5*ZCFLMOD                                                            
       CFLCOEF=R6*(1.0-ZCFLMOD*ZCFLMOD)*(DZD*DZD*RDZC)                             
!
       if (FZDELS < 0.0) then
          RFZDELS=1.0/(FZDELS-sqrt(epsilon(1.0_wp)))
       else
          RFZDELS=1.0/(FZDELS+sqrt(epsilon(1.0_wp)))
       endif
       FZFS=0.5*(FZD+FZC)-HCFL*(FZD-FZC)-           &
     &     CFLCOEF*((FZD-FZC)*RDZD-(FZC-FZU)*RDZU) - FGT1                          
       FZNFS=(FZFS-FZU)*RFZDELS                                                    
       FZNCS=(FZC-FZU)*RFZDELS                                                     
       FZNREFS=FZNCS*HRCFL                                                         
       FZTEMP=MAX(FZNCS,MIN(FZNFS,FZNREFS,1.0))                                    
       FFLXB(J,K+1) = ( FZTEMP*FZDELS + FZU )                                      
       ENDIF  ! first-order upwind                                                 
!                                                                                  
!-----CODE FOR K=2,KKP-2 ------------------------------------                      
       DO 110 K=2,KKP-2                                                            
        SIGNW=SIGN(.5,W(J,K))                                                      
        ADVPOS=0.5+SIGNW                                                           
        ADVNEG=0.5-SIGNW                                                           
        KOFSET=NINT(ADVNEG)                                                        
        FZU = ZF(J,K-1)*ADVPOS + ZF(J,K+2)*ADVNEG                                  
        FZC = ZF(J,K+KOFSET)                                                       
        FZD = ZF(J,K+1)*ADVPOS + ZF(J,K)  *ADVNEG                                  
        FZDELS = FZD - FZU                                                         
        FZCURVS = FZD - 2.0*FZC + FZU                                              
        IF( ABS(FZCURVS) .GE. ABS(FZDELS) )THEN                                    
          FFLXB(J,K+1) = FZC
        ELSE                                                                       
        SUM_CFL_OUT=DT*( ARDZ(K+KOFSET)*                           &
     &   (MAX(0.,W(J,K+KOFSET))+ABS(MIN(0.,W(J,K-1+KOFSET))))+     &                
     &    CY*(MAX(0.,V(J,K+KOFSET))+ABS(MIN(0.,V(J-1,K+KOFSET)))) )                
!..1ST GRADT TERM, IN Y DIRN                                                       
           VFACE=0.125*(V(J-1,K)+V(J,K)+V(J-1,K+1)+V(J,K+1))                       
           JSU=J-NINT(SIGN(1.0,VFACE))                                             
           FGT1=ABS(VFACE)*CY*DT*(FZC-ZF(JSU,K+KOFSET))                            
!..2ND GRADT TERM, IN X DIRN, IS ZERO                                              
        RDZU=ARDZN(K) *ADVPOS + ARDZN(K+2)*ADVNEG                                  
        RDZC=ARDZ(K+KDZOF)*ADVPOS + ARDZ(K+1+KDZOF)*ADVNEG                         
!        !corrected at V1.5 (previously wrong for non-uniform DZ)                  
!                                                                                  
!..CALCULATE FLUXES                                                                
        DZD=ADZN(K+1)                                                              
        RDZD=ARDZN(K+1)                                                            
        ZCFLMOD=ABS(W(J,K)*ARDZN(K+1)*DT)                                          
        HRCFL=1./(SUM_CFL_OUT+epsilon(1.d0))                                              
        HCFL=0.5*ZCFLMOD                                                           
        CFLCOEF=R6*(1.0-ZCFLMOD*ZCFLMOD)*(DZD*DZD*RDZC)                            
!
        if (FZDELS < 0.0) then
           RFZDELS=1.0/(FZDELS-sqrt(epsilon(1.0_wp)))
        else
           RFZDELS=1.0/(FZDELS+sqrt(epsilon(1.0_wp)))
        endif
        
        FZFS=0.5*(FZD+FZC)-HCFL*(FZD-FZC)-                   &                      
     &   CFLCOEF*((FZD-FZC)*RDZD-(FZC-FZU)*RDZU) - FGT1  
        FZNFS=(FZFS-FZU)*RFZDELS                                                   
        FZNCS=(FZC-FZU)*RFZDELS                                                    
        FZNREFS=FZNCS*HRCFL                                                        
        FZTEMP=MAX(FZNCS,MIN(FZNFS,FZNREFS,1.0))                                   
        FFLXB(J,K+1) = ( FZTEMP*FZDELS + FZU ) 
        !print *, 'flux fzu', fflxb(j, K+1), j, k
        ENDIF ! first-order upwind                                                 
!                                                                                  
 110   CONTINUE                                                                    
 100  CONTINUE                                                                     
!                                                                                  
!.....Y-DIRECTION FLUXES  (ALL IN ONE GO)                                          
!                                                                                  
      DO 111 K=2,KMAX                                                              
        KP1=MIN(K+1,KKP)                                                           
        KP0=KP1-1                                                                  
        ! Occasionally in the J=0 loop, J-1 values are required                    
        ! which are not held in the ZV, ZF arrays for the multi                    
        ! PE case. This requires the use of additional halo                        
        ! arrays and the reinstatement of an IF test for J=0                       
        ! (removed at V1.4). To prevent this test slowing the code                 
        ! when not required the J=0 loop has been split from the                   
        ! main J loop for the FFLXL fluxes at V2.3                                 
        ! Thus, JPERIOD and AJEQ0 can also be removed from the                     
        ! main J loop.                                                             
       J=0                                                                         
!..SET UP STENCIL                                                                  
        SIGNV=SIGN(.5,V(J,K))                                                      
        ADVPOS=0.5+SIGNV                                                           
        ADVNEG=0.5-SIGNV                                                           
        JOFSET=NINT(ADVNEG)                                                        
        JPERIOD=JJP                                                                
        IF(NPES.GT.1)THEN                                                          
          SUM_CFL_OUT=DT*( ARDZ(K)*                                &                
     &       (MAX(0.,W(J+JOFSET,K))+ABS(MIN(0.,W(J+JOFSET,K-1))))  &                
     &      +CY*( MAX(0.,V(J+JOFSET,K))                            &                
     &         +JOFSET*ABS(MIN(0.,V(J,K)))                         &                
     &         +(1-JOFSET)*ABS(MIN(0., VHALO(K)))                  &                
     &         ) )                                                                 
          FYU=ZFHALO(K)*ADVPOS + ZF(J+2,K)*ADVNEG                                  
        ELSE                                                                       
          SUM_CFL_OUT=DT*( ARDZ(K)*                                &                
     &      (MAX(0.,W(J+JOFSET,K))+ABS(MIN(0.,W(J+JOFSET,K-1))))   &                
     &      +CY*(MAX(0.,V(J+JOFSET,K))                             &                
     &          +ABS(MIN(0.,V(J-1+JOFSET+JPERIOD,K)))) )                           
          FYU=ZF(J-1+JPERIOD,K)*ADVPOS + ZF(J+2,K)*ADVNEG                          
        ENDIF                                                                      
        FYC=ZF(J,K)  *ADVPOS + ZF(J+1,K)*ADVNEG                                    
        FYD=ZF(J+1,K)*ADVPOS + ZF(J,K)  *ADVNEG                                    
!..1ST GRADT TERM, IN Z DIRN                                                       
        WFACE=0.125*(W(J,K)+W(J,K-1)+W(J+1,K)+W(J+1,K-1))                          
        KPOSWF=NINT(0.5+SIGN(0.5,WFACE))                                           
        KSU=KP1*(1-KPOSWF)+(K-1)*KPOSWF                                            
        KSL=KP0*(1-KPOSWF)+K*KPOSWF                                                
        FGT1=ABS(WFACE)*ARDZ(K+KDZOF)*DT*                           &               
     &    (ZF(J+JOFSET,KSL)-ZF(J+JOFSET,KSU))                                      
!..2ND GRADT TERM, IN X DIRN, IS ZERO                                              
!..CALCULATE FLUX                                                                  
        YCFLMOD=ABS(V(J,K)*CY*DT  )                                                
        HRCFL=1./(SUM_CFL_OUT+epsilon(1.d0))                                              
        HCFL=0.5*YCFLMOD                                                           
        CFLCOEF=R6*(1.0-YCFLMOD*YCFLMOD)                                           
        FYDELS=FYD-FYU
        if (fydels > 0.0) then
           RFYDELS=1.0/(FYDELS+sqrt(epsilon(1.0_wp)))
        else
           RFYDELS=1.0/(FYDELS-sqrt(epsilon(1.0_wp)))
        endif
        FYCURVS=FYD-2.0*FYC+FYU
        FYFS=0.5*(FYD+FYC)-HCFL*(FYD-FYC)-CFLCOEF*FYCURVS            &              
          - FGT1             !   GRADT TERM (ADDED BEFORE LIMITING)
        FYNFS=(FYFS-FYU)*RFYDELS                                                   
        FYNCS=(FYC-FYU)*RFYDELS
        FYNREFS=FYNCS*HRCFL                                                        
        SDIFFF=SIGN(0.5,ABS(FYCURVS)-ABS(FYDELS))                                  
        FYTEMP=MAX(FYNCS,MIN(FYNFS,FYNREFS,1.0))                                   
        FFLXL(J+1,K)=FYC*(0.5+SDIFFF)                                &              
     &     +(FYTEMP*FYDELS+FYU)*(0.5-SDIFFF) 
! END OF J=0                                                                       
       DO J=1,JJP-1                                                                
!..SET UP STENCIL                                                                  
        SIGNV=SIGN(.5,V(J,K))                                                      
        ADVPOS=0.5+SIGNV                                                           
        ADVNEG=0.5-SIGNV                                                           
        JOFSET=NINT(ADVNEG)                                                        
        FYU = ZF(J-1,K)*ADVPOS + ZF(J+2,K)*ADVNEG                                  
        FYC = ZF(J,K)  *ADVPOS + ZF(J+1,K)*ADVNEG                                  
        FYD = ZF(J+1,K)*ADVPOS + ZF(J,K)  *ADVNEG                                  
        FYDELS = FYD - FYU                                                         
        FYCURVS = FYD - 2.0*FYC + FYU                                              
        IF( ABS(FYCURVS) .GE. ABS(FYDELS) )THEN                                    
          FFLXL(J+1,K) = FYC 
        ELSE                                                                       
        SUM_CFL_OUT=DT*( ARDZ(K)*                                     &             
     &   (MAX(0.,W(J+JOFSET,K))+ABS(MIN(0.,W(J+JOFSET,K-1))))         &             
     &  +CY*(MAX(0.,V(J+JOFSET,K))                                    &             
     &          +ABS(MIN(0.,V(J-1+JOFSET,K)))) )                                   
!..1ST GRADT TERM, IN Z DIRN                                                       
        WFACE=0.125*(W(J,K)+W(J,K-1)+W(J+1,K)+W(J+1,K-1))                          
        KPOSWF=NINT(0.5+SIGN(0.5,WFACE))                                           
        KSU=KP1*(1-KPOSWF)+(K-1)*KPOSWF                                            
        KSL=KP0*(1-KPOSWF)+K*KPOSWF                                                
        FGT1=ABS(WFACE)*ARDZ(K+KDZOF)*DT*                             &             
     &    (ZF(J+JOFSET,KSL)-ZF(J+JOFSET,KSU))                                      
!..2ND GRADT TERM, IN X DIRN, IS ZERO                                              
!..CALCULATE FLUX                                                                  
        YCFLMOD=ABS(V(J,K)*CY*DT  )                                                
        HRCFL=1./(SUM_CFL_OUT+epsilon(1.d0))                                              
        HCFL=0.5*YCFLMOD                                                           
        CFLCOEF=R6*(1.0-YCFLMOD*YCFLMOD)
        if (fydels > 0.0) then
           RFYDELS=1.0/(FYDELS+sqrt(epsilon(1.0_wp)))
        else
           RFYDELS=1.0/(FYDELS-sqrt(epsilon(1.0_wp)))
        endif

        FYFS=0.5*(FYD+FYC)-HCFL*(FYD-FYC)-CFLCOEF*FYCURVS             &             
     &     - FGT1             !   GRADT TERM (ADDED BEFORE LIMITING)               
        FYNFS=(FYFS-FYU)*RFYDELS                                                   
        FYNCS=(FYC-FYU)*RFYDELS                                                    
        FYNREFS=FYNCS*HRCFL                                                        
        FYTEMP=MAX(FYNCS,MIN(FYNFS,FYNREFS,1.0))                                   
        FFLXL(J+1,K) = ( FYTEMP*FYDELS + FYU ) 
        ENDIF  ! first-order upwind                                                
       ENDDO                                                                       
 111  CONTINUE                                                                     
                                                                                   
!..DEAL WITH VARIOUS BOUNDARY VALUES                                               
                                                                                   
      IF(KMIN.EQ.2)THEN     ! ALL FIELDS EXCEPT W                                  
!                                                                                  
!   SET LOWER BOUNDARY TO ARBITRARY VALUE SINCE IT WILL BE MULTIPLIED              
!   BY W(K=1)=0 WHEN USED.                                                         
!                                                                                  
       DO 402 J=JMINP,JMAXP                                                        
        FFLXB(J,2)=0.                                                              
 402   CONTINUE                                                                    
      ENDIF                                                                        
      ! Wrap FFLXL for calculation in calling routine                              
!       IF(NPES.GT.1)THEN                                                            
!         HALOSEND(1)=0.0                                                            
!         DO K=2,KKP                                                                 
!           HALOSEND(K)=FFLXL(1,K)                                                   
!         ENDDO                                                                      
!         CALL GC_GSYNC(NPROC,INFO)                                                  
!         CALL GC_RSEND(4113,KKP,ILEFT,INFO,HALORECV,HALOSEND)                       
!         CALL GC_RRECV(4113,KKP,IRIGHT,INFO,HALORECV,HALOSEND)                      
!         CALL GC_GSYNC(NPROC,INFO)                                                  
!         DO K=1,KKP                                                                 
!            FFLXL(JJP+1,K)=HALORECV(K)                                              
!         ENDDO                                                                      
!       ELSE  ! single PE case                                                       
      DO K=2,KKP                                                                 
         FFLXL(JJP+1,K)=FFLXL(1,K)
      ENDDO

!      ENDIF                                                                        
!                                                                                  
      IF(JMINP.EQ.0)THEN                                                           
       DO 400 K=2,KKP                                                              
        FFLXB(0,K)=0.                                                              
        FFLXB(JJP+1,K)=0.                                                          
        FFLXL(0,K)=0.                                                              
 400   CONTINUE                                                                    
      ENDIF                                                                        
      RETURN                                                                       
      END  SUBROUTINE
                                                                                   
!SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS           

end module ultf2d_mod
