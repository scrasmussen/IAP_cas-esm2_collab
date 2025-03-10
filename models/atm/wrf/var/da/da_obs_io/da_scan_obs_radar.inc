subroutine da_scan_obs_radar (iv, filename)

   !---------------------------------------------------------------------------
   ! Purpose: Scan the radar observation file
   !---------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(inout) :: iv
   character(len=*),  intent(in)    :: filename

   integer                       :: i, j, n, iost, nlevels, fm
   integer                       :: file_radar
   integer                       :: iunit

   type (radar_multi_level_type) :: platform

   character (LEN = 120)         :: char_file_radar
   character (LEN = 120)         :: char_ned

   logical                       :: outside
   integer                       :: n_dup, ndup

   real                          :: obs_time
   integer                       :: iyear, imonth, iday, ihour, imin

   if (trace_use) call da_trace_entry("da_scan_obs_radar")

   ! 1. open file
   ! ============

   call da_get_unit(iunit)
   open(unit   = iunit,     &
        FILE   = trim(filename), &
        FORM   = 'FORMATTED',  &
        ACCESS = 'SEQUENTIAL', &
        iostat =  iost,     &
        STATUS = 'OLD')

   if (iost /= 0) then
      ! Does not matter of radar file missing
      call da_warning(__FILE__,__LINE__, &
         (/"Cannot open radar file "//filename/))
      call da_free_unit(iunit) 
      if (trace_use) call da_trace_exit("da_scan_obs_radar")
      return
   end if

   ! 2. read total radar
   ! ===================

   ! 2.1 read first line
   !     ---------------

   read (unit=iunit, fmt = '(A)', iostat = iost) char_file_radar
   if (iost /= 0) then
      ! Does matter if present and unreadable
      call da_error(__FILE__,__LINE__, &
         (/"Cannot read radar file"/))
   end if

   ! 2.3 total radar number

   read (unit=char_file_radar (15:17),fmt='(I3)', iostat = iost) file_radar

   ! 2.4 skip one lines

   read (unit=iunit, fmt = '(A)', iostat = iost)

   ! 3. read radar data

   do n = 1, file_radar

      ! 3.1 skip one blank line

      read (unit=iunit, fmt = '(A)', iostat = iost)

      ! 3.2 read header

      read (unit=iunit, fmt = '(A)', iostat = iost) char_ned

      ! 3.3 read header information

      read (unit=char_ned (69:74), fmt='(I6)', iostat = iost) platform % stn % numobs

      ! 3.4 skip two lines

      read (unit=iunit, fmt = '(A)', iostat = iost)
      read (unit=iunit, fmt = '(A)', iostat = iost)

      ! 3.5 loop over records

      reports: do j = 1, platform % stn % numobs

         ! 3.5.1 read station general info

         read (unit = iunit, iostat = iost, &
                      fmt = '(A12,3X,A19,2X,2(F12.3,2X),F8.1,2X,I6)') &
                      platform % info % platform,  &
                      platform % info % date_char, &
                      platform % info % lat,       &
                      platform % info % lon,       &
                      platform % info % elv,       &
                      platform % info % levels

         read(unit=platform % info % platform (4:6), fmt='(I3)') fm

         !     3.5.2 read each level

         do i = 1, platform % info % levels
            ! height
            platform%each (i) = radar_each_level_type(missing_r, missing, -1.0,&
               field_type(missing_r, missing, missing_r, missing, missing_r), & ! rv
               field_type(missing_r, missing, missing_r, missing, missing_r))   ! rf

            read (unit = iunit, fmt = '(3X, F12.1, 2(F12.3,I4,F12.3,2X))') &
                             platform % each (i) % height,           &
                             platform % each (i) % rv % inv,         &
                             platform % each (i) % rv % qc,          &
                             platform % each (i) % rv % error,       &
                             platform % each (i) % rf % inv,         &
                             platform % each (i) % rf % qc,          &
                             platform % each (i) % rf % error
         end do

        ! Check if outside of the time range:

         read (platform%info%date_char,'(i4,4(1x,i2))') &
               iyear, imonth, iday, ihour, imin
         call da_get_julian_time (iyear,imonth,iday,ihour,imin,obs_time)
         if ( obs_time < time_slots(0) .or. &
              obs_time >= time_slots(num_fgat_time) ) then
            cycle
         endif

         call da_llxy (platform%info, platform%loc, outside)

         nlevels = platform%info%levels

         if (nlevels > max_ob_levels) then
             write(unit=message(1),fmt='(A,2I8)') &
                ' radar=> nlevels > max_ob_levels:',nlevels, max_ob_levels
             call da_warning(__FILE__,__LINE__,message(1:1))

             nlevels = max_ob_levels
             platform%info%levels = nlevels
         else if (nlevels < 1) then
            cycle reports
         end if

         iv%info(radar)%ntotal = iv%info(radar)%ntotal + 1
         if (outside) then
            cycle reports
         end if

         ! Loop over duplicating obs for global
         n_dup = 1
         if (global .and. &
            (platform%loc%i == ids .or. platform%loc%i == ide)) n_dup= 2
   
         do ndup = 1, n_dup
            select case (fm)

            case (128)
               iv%info(radar)%nlocal = iv%info(radar)%nlocal + 1

               if (iv%info(radar)%nlocal > max_radar_input) then
                  write(unit=message(1),fmt='(A,I6,A,I6)') &
                     ' radar #= ',iv%info(radar)%nlocal, ' > max_radar_input = ', max_radar_input
                  call da_error(__FILE__,__LINE__,message(1:1))
               end if

            case default;
               write(unit=stdout, fmt='(a)') 'Warning: unsaved obs found:'

               write(unit=stdout, fmt='(2a)') &
                  'platform % info % platform=', platform % info % platform

               write(unit=stdout, fmt='(a, i3)') &
                  'platform % info % levels=', platform % info % levels
            end select

         iv%info(radar)%max_lev = max(iv%info(radar)%max_lev, platform%info%levels)
         end do        !  loop over duplicate
      end do reports
   end do

   close (iunit)
   call da_free_unit(iunit)

   if (trace_use) call da_trace_exit("da_scan_obs_radar")


end subroutine da_scan_obs_radar


