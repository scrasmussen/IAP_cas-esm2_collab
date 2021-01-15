#include <define.h>
 
subroutine lpwrite(idate_p,idate,fout,frestart,cdate)

  use precision
  use colm_varctl
  use nchistMod
  use timemgr, only : dtime, get_lastday_of_month

  implicit none

  integer, intent(in) :: idate_p(3)
  integer, intent(in) :: idate(3)
  character(LEN=255), intent(in)  :: fout 
  character(LEN=255), intent(inout) :: frestart
  character(LEN=255), intent(out) :: cdate

  integer year  , month  , day  , hour  , minute  , second
  integer yearp, monthp, dayp, hourp, minutep, secondp
  integer m, n, nyear

  monthp = 0
  do m = 1, 12
     if(idate_p(2).le.get_lastday_of_month(m))then
        monthp = m
        exit
     end if
  end do

  yearp   = idate_p(1)
  dayp    = idate_p(2) - get_lastday_of_month(monthp-1)
  hourp   = idate_p(3)/3600
  minutep = (idate_p(3)-3600*hourp)/60
  secondp = mod(idate_p(3),60)

  month = 0
  do m = 1, 12
     if(idate(2).le.get_lastday_of_month(m))then
        month = m
        exit
     end if
  end do

  year   = idate(1)
  day    = idate(2) - get_lastday_of_month(month-1)
  hour   = idate(3)/3600
  minute = (idate(3)-3600*hour)/60
  second = mod(idate(3),60)

  do n = 1, nMaxHist

     if(histArray(n)%is_valid) then
        histArray(n)%is_ready = .false.
        histArray(n)%is_newyear = .false.

        if(yearp.ne.year) histArray(n)%is_newyear = .true.

        if(histArray(n)%y_interval.gt.0 .and. &
           mod(yearp,histArray(n)%y_interval).eq.0 .and. yearp.ne.year) then
           write(cdate,'(i4.4)') yearp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%m_interval.gt.0 .and. &
           mod(monthp,histArray(n)%m_interval).eq.0 .and. monthp.ne.month) then
           write(cdate,'(i4.4,"-",i2.2)') yearp,monthp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%d_interval.gt.0 .and. &
           mod(dayp,histArray(n)%d_interval).eq.0 .and. dayp.ne.day) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') yearp,monthp,dayp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%h_interval.gt.0 .and. &
           mod(hourp,histArray(n)%h_interval).eq.0 .and. hourp.ne.hour) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i2.2)') yearp,monthp,dayp,hourp
           histArray(n)%fhistory = trim(fout)//'-'//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if
     end if

  end do

  write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1),idate(2),idate(3)

  if((trim(restart_freq).eq."yearly" .and. (yearp.ne.year)) .or. &
     (trim(restart_freq).eq."monthly" .and. (monthp.ne.month)) .or. &
     (trim(restart_freq).eq."daily" .and. (dayp.ne.day))) then
     frestart = trim(fout)//'-restart-'//trim(cdate)
  endif

! Used as metadata in netcdf file
  write(cdate,'(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') &
                year,   month,     day,    hour,  minute,  second

end subroutine lpwrite
! ------------------------------------------------------------------------
! EOP
