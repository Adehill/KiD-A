! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Class for the microphysics mixing ratio/ number concentration, etc. 
! for each species 
!
! 
!
!

module class_species

  Use typeKind
  Use parameters, only :max_nmoments, max_nbins, unset_real

  implicit none

  type, public :: species
     integer :: h_id
     integer :: nmoments
     integer :: nbins
     real(wp) , pointer :: moments(:,:)
  end type species

  type, public :: species_hw
     ! As species, but with hard-wired sizes for moments
     integer :: h_id
     integer :: nmoments
     integer :: nbins
     real(wp) :: moments(max_nbins,max_nmoments)
  end type species_hw

  interface species_allocate
     module procedure species_allocate_pnt, species_allocate_hw
  end interface

contains

  subroutine species_allocate_pnt(var, nmoments, nbins, id &
       , initval)
    !
    ! initialize and allocate space for species
    !
    
    type(species), intent(inout) :: var
    integer, intent(in) :: nmoments, nbins, id
    real, intent(in), optional :: initval
    real :: cinitval

    cinitval=unset_real
    if (present(initval))cinitval=initval

    allocate(var%moments(max(1,nbins), max(1,nmoments)))
    var%moments=cinitval
    var%nmoments=nmoments
    var%nbins=nbins
    var%h_id=id
    
  end subroutine species_allocate_pnt


  subroutine species_allocate_hw(var, nmoments, nbins, id &
       , initval)
    !
    ! initialize and allocate space for species
    ! This is for hardwired species, so only don't
    ! do any allocation
    !

    type(species_hw), intent(inout) :: var
    integer, intent(in) :: nmoments, nbins, id
    real, intent(in), optional :: initval
    real :: cinitval

    cinitval=unset_real
    if (present(initval))cinitval=initval

    ! don't do any allocation
    var%moments=cinitval
    var%nmoments=nmoments
    var%nbins=nbins
    var%h_id=id
    
  end subroutine species_allocate_hw


end module class_species
