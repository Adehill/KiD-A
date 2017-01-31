! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Sets the initial profiles
!
!

module set_profiles

  Use typeKind
  Use netcdf
  Use input_variables
  Use switches
  Use parameters, only : nz, nx
  Implicit none

  interface read_profiles
     module procedure read_profiles_file, read_profiles_standard
  end interface

contains

  subroutine read_profiles_file(file, input_type)

    character(len=*), optional :: file
    integer, optional :: input_type

    ! Local variables    
    integer :: cinput_type
    integer(incdfp) :: ncid, status            &
         , timeid, nlevid, nlevp1id,p_surfid   &
         , timedimid, nlevdimid, nlevp1dimid   &
         , coor_lena, coor_lenb, xtype, attnum &
         , Tid, qid, Tforceid, qforceid, omegaid   &
         , len

    if (present(input_type))then
       cinput_type=input_type
    else
       cinput_type=iecmwf
    end if

    select case(cinput_type)
!!$    case(iecmwf)
!!$       status=NF90_OPEN(trim(file), nf90_noWrite, ncid)
!!$       status=NF90_INQ_DIMID(ncid,'time',timedimid)
!!$       status=NF90_INQ_DIMID(ncid,'nlev',nlevdimid)
!!$       status=NF90_INQ_DIMID(ncid,'nlevp1',nlevp1dimid)
!!$       status=NF90_INQ_VARID(ncid,'time',timeid)
!!$       status=NF90_INQ_VARID(ncid,'nlev',nlevid)
!!$       status=NF90_INQ_VARID(ncid,'nlevp1',nlevp1id)
!!$       status=NF90_INQ_VARID(ncid,'ps',p_surfid)
!!$       status=NF90_INQ_VARID(ncid,'t',Tid)
!!$       status=NF90_INQ_VARID(ncid,'q',qid)
!!$       status=NF90_INQ_VARID(ncid,'tadv',Tforceid)
!!$       status=NF90_INQ_VARID(ncid,'qadv',qforceid)
!!$       status=NF90_INQ_VARID(ncid,'omega',omegaid)
!!$       status=NF90_Inquire_Dimension(ncid, timedimid, len=len)
!!$       nt_in=len
!!$       status=NF90_Inquire_Dimension(ncid, nlevdimid, len=len)
!!$       nk_in=len
!!$       status=NF90_Inquire_Dimension(ncid, nlevp1dimid, len=len)
!!$       nk_half_in=len
!!$       status=NF90_Inquire_Attribute(ncid, NF90_Global, 'coor_par_a',&
!!$            & xtype, &
!!$            & coor_lena, attnum)
!!$       status=NF90_Inquire_Attribute(ncid, NF90_Global, 'coor_par_a',&
!!$            & xtype, &
!!$            & coor_lenb, attnum)
!!$       if (nz /= nk_half_in)then
!!$          print*, 'Incorrect number of levels in file'
!!$          print*, 'nk_half_in=',nk_half_in
!!$          print*, 'nk_in=',nk_in
!!$          print*, 'coor_lena=',coor_lena
!!$          print*, 'coor_lenb=',coor_lenb
!!$          print*, 'nz=', nz
!!$       else
!!$          allocate(coor_par_a(coor_lena))
!!$          allocate(coor_par_b(coor_lenb))
!!$          status=NF90_Get_Att(ncid, NF90_Global, 'coor_par_a',&
!!$               & coor_par_a)
!!$          status=NF90_Get_Att(ncid, NF90_Global, 'coor_par_b',&
!!$               & coor_par_b)
!!$       end if
!!$       allocate(time_in(nt_in))
!!$       allocate(nlev_in(nk_in))
!!$       allocate(nlevp1_in(nk_half_in))
!!$       allocate(p_surf_in(nt_in))
!!$       allocate(T_in(nk_in,nt_in))
!!$       allocate(q_in(nk_in,nt_in))
!!$       allocate(Tforce_in(nk_in,nt_in))
!!$       allocate(qforce_in(nk_in,nt_in))
!!$       allocate(omega_in(nk_in,nt_in))
!!$       status=NF90_GET_VAR(ncid, timeid, time_in)
!!$       status=NF90_GET_VAR(ncid, nlevid, nlev_in)
!!$       status=NF90_GET_VAR(ncid, nlevp1id, nlevp1_in)
!!$       status=NF90_GET_VAR(ncid, p_surfid, p_surf_in)
!!$       status=NF90_GET_VAR(ncid, Tid, T_in)
!!$       status=NF90_GET_VAR(ncid, qid, q_in)
!!$       status=NF90_GET_VAR(ncid, Tforceid, Tforce_in)
!!$       status=NF90_GET_VAR(ncid, qforceid, qforce_in)
!!$       status=NF90_GET_VAR(ncid, omegaid, omega_in)
!!$       status=NF90_CLOSE(ncid)
    end select

  end subroutine read_profiles_file

  subroutine read_profiles_standard(icase)

    Use test_cases, only : standard_cases
    Use test_cases_2d, only : standard_2d_cases

    integer, intent(in) :: icase

    if (nx == 1) then 

       call standard_cases(icase)

    else
       
       call standard_2d_cases(icase)
       
    endif


  end subroutine read_profiles_standard



end module set_profiles
