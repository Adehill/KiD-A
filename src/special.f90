Module special

  Use typeKind
  
  Implicit None

contains

  !==================================!
  ! Gamma function                   !
  !==================================!

  real function gamma(x, y)  
    implicit none    
    
    real(wp), intent(in) :: x
    real(wp), intent(in), optional :: y
    real(wp) :: f,g,z
    REAL(WP), PARAMETER :: pi=3.14159
    
    if (.not. present(y))then
    
      f=huge(x)
      g=1
      z=x
      if ((z+int(abs(z))).ne.0)then
        do while(z<3)
          ! Lets use a recursion relation for Gamma functions
          ! to get a large argument and use Stirlings formula
          g=g*z
          z=z+1
        end do
        
        ! This is just stirlings formula...
        f=(1.-2.*(1-2./(3.*z*z))/(7.*z*z))/(30.*z*z)
        f=(1.-f)/(12.*z)+z*(log(z)-1)
        f=(exp(f)/g)*sqrt(2.*pi/z)
        
      endif
      gamma=f
      
    else

!!$       gamma=max(0., gammaup(x,y))

    end if
    
  end function gamma


  !==================================!
  ! Upper incomplete gamma function  !
  !==================================!
!!$
!!$      real(wp) function gammaup(a,x)
!!$      
!!$      implicit none
!!$
!!$      real(wp) :: a,x,gammafuncnew,eigamma
!!$   
!!$! integral(x to infty) t^a-1 exp(-t) dt
!!$   
!!$      gammaup=(gamma(a))*(1-x**a*eigamma(a,x))
!!$   
!!$      end function gammaup

  !==================================!
  ! Entire Incomplete gamma function !
  !==================================!

!!$      real(wp) function eigamma(a,x)
!!$      
!!$      implicit none
!!$      
!!$      real(wp), intent(in) :: a,x
!!$      real(wp) :: p,g,b,z,f
!!$      integer :: j
!!$      real(wp), parameter :: pi=3.14159
!!$
!!$      if (x < 100.0) then     
!!$ 
!!$      g=1
!!$      p=1
!!$      b=a
!!$      do while(b<=1.)
!!$         b=b+1
!!$         g=g*x
!!$         p=p*b+g
!!$      end do
!!$      b=b+1
!!$
!!$      j=int(5.*(3.+abs(x))/2.)
!!$      f=1./(j+b-x)
!!$
!!$      do while(j/=0)
!!$         j=j-1
!!$         f=(f*x+1.)/(j+b)
!!$      end do
!!$
!!$      p=p+f*g*x
!!$      g=(1.-(2./(7.*b*b))*(1.-2./(3.*b*b)))/(30.*b*b)
!!$      g=(g-1.)/(12.*b)-b*(log(b)-1.)
!!$
!!$      f=p*exp(g-x)*sqrt(b/(2.*pi))
!!$
!!$      eigamma=f
!!$   
!!$      else
!!$     
!!$      eigamma=0.0
!!$
!!$      endif
!!$   
!!$      end function eigamma
   

end Module special
