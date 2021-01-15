!=======================================================================
! A New Algorithm for Calculation of Consine Solar Zenith Angle
! Author: Linjiong Zhou
! Verified: Minghua Zhang
! E-mail: linjiongzhou@hotmail.com
! Date  : 2015.02.22
! Ref.  : Zhou et al., GRL, 2015
!
! Usage:
! put this module file into cesm?_?_?/models/csm_share/shr
! replace "shr_orb_cosz = ..." in shr_orb_mod.F90 by:
! "call new_cosz(lat, lon, declin, jday, shr_orb_cosz)"
!=======================================================================

module new_cosz_mod

    use shr_kind_mod
    use shr_file_mod
    use shr_sys_mod
    use shr_log_mod

    implicit none
    private
    save

    real(shr_kind_r8) :: rdt   ! radiation time step (second)

    logical :: initialized = .false.

    public :: new_cosz

contains
!=======================================================================
subroutine new_cosz(lat, lon, dec, jday, cosz)

    implicit none

!-----------------------------------------------------------------------
! In/Out Arguements

    real(shr_kind_r8), intent(in) :: lat   ! latitude (radian)
    real(shr_kind_r8), intent(in) :: lon   ! longitude (radian)
    real(shr_kind_r8), intent(in) :: dec   ! solar declination (radian)
    real(shr_kind_r8), intent(in) :: jday  ! Julian calendar day (1.xx to 365.xx)

    real(shr_kind_r8), intent(out) :: cosz ! cosine solar zenith angle (1)

!-----------------------------------------------------------------------
! Local Arguments

    real(shr_kind_r8) :: aa, bb
    real(shr_kind_r8) :: pi
    real(shr_kind_r8) :: del, phi
    real(shr_kind_r8) :: cos_h, h
    real(shr_kind_r8) :: t1, t2, dt
    real(shr_kind_r8) :: tt1, tt2, tt3, tt4

    real(shr_kind_r8) :: dtime ! dtime = nnnn, Model time step in seconds. Default is dycore dependent.
    integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                               ! or hours (negative).  Default: -1

    integer :: ierr            ! error code
    integer :: unitn           ! namelist unit number

!-----------------------------------------------------------------------
! All variable in the namelist group: cam_inparm, useless variables

    integer            :: phys_alltoall, phys_chnk_per_thd, phys_loadbalance, phys_twin_algorithm, nlvdry, iradae,&
                          iradlw, irad_always, rayk0
    integer            :: mfilt(6), lcltod_start(6), lcltod_stop(6), ndens(6), nhtfrq(6)
    real(shr_kind_r8)  :: pertlim, raykrange, raytau0
    character(len=1)   :: avgflag_pertape(6)
    character(len=8)   :: inithist
    character(len=24)  :: fexcl1(750), fexcl2(750), fexcl3(750), fexcl4(750), fexcl5(750), fexcl6(750)
    character(len=26)  :: fincl1(750), fincl2(750), fincl3(750), fincl4(750), fincl5(750), fincl6(750),           &
                          fwrtpr1(750), fwrtpr2(750), fwrtpr3(750), fwrtpr4(750), fwrtpr5(750), fwrtpr6(750)
    character(len=128) :: iopfile
    character(len=128) :: fincl1lonlat(750), fincl2lonlat(750), fincl3lonlat(750), fincl4lonlat(750),             &
                          fincl5lonlat(750), fincl6lonlat(750)
    character(len=200) :: scm_clubb_iop_name
    character(len=256) :: ncdata, absems_data, cam_branch_file, bnd_topo, bnd_topo2, efield_hflux_file, efield_lflux_file,   &
                          efield_wei96_file, qbo_forcing_file
    character(len=256) :: hfilename_spec(6)
    logical            :: print_energy_errors, empty_htapes, inithist_all, readtrace, use_64bit_nc,               &
                          print_step_cost, pbuf_global_allocate, spectralflux, scm_crm_mode, scm_diurnal_avg,     &
                          scm_iop_srf_prop, scm_relaxation, tracers_flag, nlte_use_mo, qbo_cyclic, qbo_use_forcing
    logical            :: collect_column_output(6)
    namelist /cam_inparm/ phys_alltoall, phys_chnk_per_thd, phys_loadbalance, &
                          phys_twin_algorithm, nlvdry, iradae, &
                          iradlw, irad_always, rayk0
    namelist /cam_inparm/ mfilt, lcltod_start, lcltod_stop, ndens, nhtfrq
    namelist /cam_inparm/ pertlim, raykrange, raytau0
    namelist /cam_inparm/ avgflag_pertape
    namelist /cam_inparm/ inithist
    namelist /cam_inparm/ fexcl1, fexcl2, fexcl3, fexcl4, fexcl5, fexcl6
    namelist /cam_inparm/ fincl1, fincl2, fincl3, fincl4, fincl5, fincl6, &
                          fwrtpr1, fwrtpr2, fwrtpr3, fwrtpr4, &
                          fwrtpr5, fwrtpr6
    namelist /cam_inparm/ iopfile
    namelist /cam_inparm/ fincl1lonlat, fincl2lonlat, fincl3lonlat, fincl4lonlat, fincl5lonlat, fincl6lonlat
    namelist /cam_inparm/ scm_clubb_iop_name
    namelist /cam_inparm/ ncdata, absems_data, cam_branch_file, bnd_topo, bnd_topo2, &
                          efield_hflux_file, efield_lflux_file, &
                          efield_wei96_file, qbo_forcing_file
    namelist /cam_inparm/ hfilename_spec
    namelist /cam_inparm/ print_energy_errors, empty_htapes, inithist_all, &
                          readtrace, use_64bit_nc, &
                          print_step_cost, pbuf_global_allocate, spectralflux, &
                          scm_crm_mode, scm_diurnal_avg, &
                          scm_iop_srf_prop, scm_relaxation, tracers_flag, &
                          nlte_use_mo, qbo_cyclic, qbo_use_forcing
    namelist /cam_inparm/ collect_column_output

!-----------------------------------------------------------------------
! Get Radiation Step from Namelist
  
    namelist /cam_inparm/ dtime, iradsw

    if (.not. initialized) then

        unitn = shr_file_getUnit() ! get an unused unit number
        open(unitn, file = 'atm_in', status = 'old', iostat = ierr)
        ! added by fkc, find cam_inparm position
        call find_group_name(unitn, 'cam_inparm', status=ierr)
        if (ierr .ne. 0) then
            write(shr_log_unit, *) 'File atm_in not found, setting default values.'
        else
            ierr = 1
            do while (ierr .ne. 0)
                read(unitn, cam_inparm, iostat = ierr)
                if (ierr .ne. 0) then
                    call shr_sys_abort('new_cosz: namelist read returns an error condition for cam_inparm')
                end if
            end do
        end if
        close(unitn)
        call shr_file_freeUnit(unitn)

        if (iradsw .lt. 0) then
            rdt = - iradsw * 3600.0_shr_kind_r8
        else
            rdt = iradsw * dtime
        end if

        initialized = .true.

    end if

!-----------------------------------------------------------------------
! Compute Half-day Length

    pi = acos(- 1.0_shr_kind_r8)

    ! adjust latitude so that its tangent will be defined
    del = lat
    if (lat .eq.   pi / 2.0_shr_kind_r8) then
        del = lat - 1.0e-05_shr_kind_r8
    end if
    if (lat .eq. - pi / 2.0_shr_kind_r8) then
        del = lat + 1.0e-05_shr_kind_r8
    end if

    ! adjust declination so that its tangent will be defined
    phi = dec
    if (dec .eq.   pi / 2.0_shr_kind_r8) then
        phi = dec - 1.0e-05_shr_kind_r8
    end if
    if (dec .eq. - pi / 2.0_shr_kind_r8) then
        phi = dec + 1.0e-05_shr_kind_r8
    end if

    ! define the cosine of the half-day length
    ! adjust for cases of all daylight or all night
    cos_h = - tan(del) * tan(phi)
    if (cos_h .le. - 1.0_shr_kind_r8) then
        h = pi
    end if
    if (cos_h .ge.   1.0_shr_kind_r8) then
        h = 0.0_shr_kind_r8
    end if
    if (cos_h .gt. - 1.0_shr_kind_r8 .and. cos_h .lt. 1.0_shr_kind_r8) then
        h = acos(cos_h)
    end if

!-----------------------------------------------------------------------
! Define Local Time t and t + dt

    ! adjust t to be between -pi and pi
    t1 = (jday - int(jday)) * 2.0_shr_kind_r8 * pi + lon - pi
    if (t1 .ge.   pi) t1 = t1 - 2.0_shr_kind_r8 * pi
    if (t1 .lt. - pi) t1 = t1 + 2.0_shr_kind_r8 * pi

    dt = rdt / 86400.0_shr_kind_r8 * 2.0_shr_kind_r8 * pi
    t2 = t1 + dt

!-----------------------------------------------------------------------
! Comput Cosine Solar Zenith angle

    ! define terms needed in the cosine zenith angle equation
    aa = sin(lat) * sin(dec)
    bb = cos(lat) * cos(dec)

    ! define the hour angle
    ! force it to be between -h and h
    ! consider the situation when the night period is too short
    if (t2 .ge. pi .and. t1 .le. pi .and. pi - h .le. dt) then
        tt2 = h
        tt1 = min(max(t1,                      - h),                        h)
        tt4 = min(max(t2, 2.0_shr_kind_r8 * pi - h), 2.0_shr_kind_r8 * pi + h)
        tt3 = 2.0_shr_kind_r8 * pi - h
    else if (t2 .ge. - pi .and. t1 .le. - pi .and. pi - h .le. dt) then
        tt2 = - 2.0_shr_kind_r8 * pi + h
        tt1 = min(max(t1, - 2.0_shr_kind_r8 * pi - h), - 2.0_shr_kind_r8 * pi + h)
        tt4 = min(max(t2,                        - h),                          h)
        tt3 = - h
    else
        if (t2 .gt. pi) then
            tt2 = min(max(t2 - 2.0_shr_kind_r8 * pi, - h), h)
        else if (t2 .lt. - pi) then
            tt2 = min(max(t2 + 2.0_shr_kind_r8 * pi, - h), h)
        else
            tt2 = min(max(t2                       , - h), h)
        end if
        if (t1 .gt. pi) then
            tt1 = min(max(t1 - 2.0_shr_kind_r8 * pi, - h), h)
        else if (t1 .lt. - pi) then
            tt1 = min(max(t1 + 2.0_shr_kind_r8 * pi, - h), h)
        else
            tt1 = min(max(t1                       , - h), h)
        end if
        tt4 = 0.0_shr_kind_r8
        tt3 = 0.0_shr_kind_r8
    end if

    ! perform a time integration to obtain cosz if desired
    ! output is valid over the period from t to t + dt
    if (tt2 .gt. tt1 .or. tt4 .gt. tt3) then
        cosz = (aa * (tt2 - tt1) + bb * (sin(tt2) - sin(tt1))) / dt + &
               (aa * (tt4 - tt3) + bb * (sin(tt4) - sin(tt3))) / dt
    else
        cosz = 0.0_shr_kind_r8
    end if

end subroutine new_cosz
!! added by fkc, this find_group_name and to_lower function is copied from
!! models/atm/src/utils/namelist_utils.F90 
subroutine find_group_name(unit, group, status)

!---------------------------------------------------------------------------------------
! Purpose: 
! Search a file that contains namelist input for the specified namelist group
! name.
! Leave the file positioned so that the current record is the first record of
! the
! input for the specified group.
! 
! Method: 
! Read the file line by line.  Each line is searched for an '&' which may only
! be preceded by blanks, immediately followed by the group name which is case
! insensitive.  If found then backspace the file so the current record is the
! one containing the group name and return success.  Otherwise return -1.
!
! Author:  B. Eaton, August 2007
!---------------------------------------------------------------------------------------

   integer,          intent(in)  :: unit     ! fortran unit attached to file
   character(len=*), intent(in)  :: group    ! namelist group name
   integer,          intent(out) :: status   ! 0 for success, -1 if group name not found

   ! Local variables

   integer           :: len_grp
   integer           :: ios    ! io status
   character(len=80) :: inrec  ! first 80 characters of input record
   character(len=80) :: inrec2 ! left adjusted input record
   character(len=len(group)) :: lc_group

   !---------------------------------------------------------------------------

   len_grp = len_trim(group)
   lc_group = to_lower(group)

   ios = 0
   do while (ios <= 0)

      read(unit, '(a)', iostat=ios, end=100) inrec

      if (ios <= 0) then  ! ios < 0  indicates an end of record condition

         ! look for group name in this record

         ! remove leading blanks
         inrec2 = adjustl(inrec)

         ! check for leading '&'
         if (inrec2(1:1) == '&') then

            ! check for case insensitive group name
            if (trim(lc_group) == to_lower(inrec2(2:len_grp+1))) then

               ! found group name.  backspace to leave file position at this
               ! record
               backspace(unit)
               status = 0
               return

            end if
         end if
      end if

   end do

   100 continue  ! end of file processing
   status = -1

end subroutine find_group_name

function to_lower(str)

!----------------------------------------------------------------------- 
! Purpose: 
! Convert character string to lower case.
! 
! Method: 
! Use achar and iachar intrinsics to ensure use of ascii collating sequence.
!
! Author:  B. Eaton, July 2001
!     
! $Id$
!----------------------------------------------------------------------- 
   implicit none

   character(len=*), intent(in) :: str      ! String to convert to lower case
   character(len=len(str))      :: to_lower

! Local variables

   integer :: i                ! Index
   integer :: aseq             ! ascii collating sequence
   integer :: upper_to_lower   ! integer to convert case
   character(len=1) :: ctmp    ! Character temporary
!-----------------------------------------------------------------------
   upper_to_lower = iachar("a") - iachar("A")

   do i = 1, len(str)
      ctmp = str(i:i)
      aseq = iachar(ctmp)
      if ( aseq >= iachar("A") .and. aseq <= iachar("Z") ) &
           ctmp = achar(aseq + upper_to_lower)  
      to_lower(i:i) = ctmp
   end do

end function to_lower
                                                                                                     




end module new_cosz_mod
