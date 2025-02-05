#include <define.h>

module timemgr

   use precision
   use colm_varctl
   use abortutils, only: endrun
   implicit none

   real(r8):: dtime                           ! time step (senconds)
   integer :: mstep                           ! model step for simulation [-]
   integer :: istep                           ! current model step

   integer :: idate(3)                        ! calendar (year, julian day, seconds)
   integer :: idate_p(3)                      ! current model calendar 

contains

   subroutine ticktime 

!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
!     Notes: Greenwich time
!
!=======================================================================

      idate_p(:) = idate(:)

      idate(3) = idate(3) + nint(dtime)

      if(idate(3)>=86400)then
         idate(3) = idate(3) - 86400
         idate(2) = idate(2) + 1

         if(idate(2)>get_ndays_of_year()) then
            idate(1) = idate(1) + 1
            idate(2) = 1
         endif
      endif

   end subroutine ticktime

   subroutine set_curr_date(idate_c)

      integer, intent(in) :: idate_c(3)

      idate = idate_c

   end subroutine set_curr_date

   subroutine get_curr_date(yr, mon, day, sec)

      integer, intent(out) :: yr
      integer, intent(out) :: mon
      integer, intent(out) :: day
      integer, intent(out) :: sec

      integer k

      do k = 1, 12
         if(idate(2) .le. get_lastday_of_month(k)) then
            mon = k

            if(k.gt.1) then
               day = idate(2) - get_lastday_of_month(k-1)
            else
               day = idate(2)
            end if

            exit
         end if
      end do

      yr = idate(1)
      sec = idate(3)

   end subroutine get_curr_date

   function is_leapyear(yr)

      integer, intent(in), optional :: yr
      logical is_leapyear

      integer year

      if (present(yr)) then
         year = yr
      else
         year = idate_p(1)
      end if

      if(enable_leapyear) then
         is_leapyear = (mod(year,4)==0 .and. mod(year,100)/=0) .or. mod(year,400)==0
      else
         is_leapyear = .false.
      end if

   end function is_leapyear

   function get_lastday_of_month(mon,yr)

      integer, intent(in) :: mon
      integer, intent(in), optional :: yr
      integer get_lastday_of_month

      integer months(0:12)
      logical leapyr

      if (present(yr)) then
         leapyr = is_leapyear(yr)
      else
         leapyr = is_leapyear()
      end if

      if (leapyr) then
         months = (/0,31,60,91,121,152,182,213,244,274,305,335,366/)
      else
         months = (/0,31,59,90,120,151,181,212,243,273,304,334,365/)
      end if

      get_lastday_of_month = months(mon)

   end function get_lastday_of_month

   function get_ndays_of_year(yr)

      integer, intent(in), optional :: yr
      integer get_ndays_of_year

      logical leapyr

      if (present(yr)) then
         leapyr = is_leapyear(yr)
      else
         leapyr = is_leapyear()
      end if

      if (leapyr) then
         get_ndays_of_year = 366
      else
         get_ndays_of_year = 365
      end if

   end function get_ndays_of_year

   function get_ndays_of_month(mon)

      integer, intent(in)  :: mon
      integer get_ndays_of_month

      get_ndays_of_month = get_lastday_of_month(mon) - get_lastday_of_month(mon-1)

   end function get_ndays_of_month

   function get_curr_calday(offset)

! Return calendar day at end of current timestep with optional offset.
! Calendar day 1.0 = 0Z on Jan 1.

! Arguments
      integer, optional, intent(in) :: offset  ! Offset from current time in seconds.
                                            ! Positive for future times, negative 
                                            ! for previous times.
! Return value
      real(r8) :: get_curr_calday

      if(present(offset)) then
         get_curr_calday = idate(2) + float(idate(3))/86400 + float(offset)/86400
      else
         get_curr_calday = idate(2) + float(idate(3))/86400
      end if

      if (get_curr_calday.gt.366) then
          get_curr_calday = get_curr_calday-365
      end if

      if ( (get_curr_calday < 1.0) .or. (get_curr_calday > 366.0) )then
         call endrun( 'get_curr_calday'//': error get_curr_calday out of bounds' )
      end if
   
   end function get_curr_calday

   subroutine get_prev_date(yr, mon, day, sec)

      integer, intent(out) :: yr
      integer, intent(out) :: mon
      integer, intent(out) :: day
      integer, intent(out) :: sec

      integer k

      do k = 1, 12
         if(idate_p(2) .le. get_lastday_of_month(k)) then
            mon = k

            if(k.gt.1) then
               day = idate_p(2) - get_lastday_of_month(k-1)
            else
               day = idate_p(2)
            end if

            exit
         end if
      end do

      yr = idate_p(1)
      sec = idate_p(3)

   end subroutine get_prev_date

   function get_curr_year()

      integer get_curr_year

      get_curr_year = idate(1) 

   end function get_curr_year

   function get_step_size()

      real(r8) get_step_size

      get_step_size = dtime

   end function get_step_size

   function get_nstep()

      integer get_nstep

      get_nstep = istep

   end function get_nstep

   subroutine get_dates_range(year1,month1,day1,second1,year2,month2,day2,second2,nsec)

      integer, intent(in) :: year1
      integer, intent(in) :: month1     ! month of the year1
      integer, intent(in) :: day1       ! day of the month1
      integer, intent(in) :: second1    ! second of the day1
      integer, intent(in) :: year2
      integer, intent(in) :: month2     ! month of the year2
      integer, intent(in) :: day2       ! day of the month2
      integer, intent(in) :: second2    ! second of the day2
      integer, intent(out):: nsec       

      integer nday, yr
      integer begyr, endyr, begmn, endmn, begdy, enddy, begsec, endsec
      logical reversed

   !* -----------------------

      reversed = .false.

      begyr = year1
      endyr = year2
      begmn = month1
      endmn = month2
      begdy = day1
      enddy = day2
      begsec = second1
      endsec = second2

   !* nsec = [begyr:begmn:begdy:begsec] ~ [endyr:endmn:enddy-1:86400] + endsec

      if ((year1 >  year2) .or. & 
          (year1 == year2 .and. month1 >  month2) .or. &
          (year1 == year2 .and. month1 == month2 .and. day1 > day2) .or. &
          (year1 == year2 .and. month1 == month2 .and. day1 == day2 .and. second1 > second2)) then
         reversed = .true.

         begyr = year2
         endyr = year1
         begmn = month2
         endmn = month1
         begdy = day2
         enddy = day1
         begsec = second2
         endsec = second1
      end if

      nday = 0
      nsec = 0

      do yr = begyr, endyr
         if (begyr.eq.endyr) then
            nday = (get_lastday_of_month(endmn-1,yr)+(enddy-1)) - (get_lastday_of_month(begmn-1,yr)+begdy) + 1
         else if (yr.eq.begyr) then
            nday = nday + (get_ndays_of_year(yr)-(get_lastday_of_month(begmn-1,yr)+begdy)+1)
         else if (yr.eq.endyr) then
            nday = nday + (get_lastday_of_month(endmn-1,yr)+enddy-1)
         else
            nday = nday + get_ndays_of_year(yr)
         end if
      end do

      nsec = nday*86400 - begsec + endsec

      if (reversed) nsec = -nsec

   end subroutine get_dates_range

end module timemgr
