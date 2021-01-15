#include <define.h>
 
subroutine lpwrite(foutdir,frestart,cdate)

  use precision
  use colm_varctl
  use timemgr
  use nchistMod

  implicit none

  character(LEN=255), intent(in)  :: foutdir
  character(LEN=255), intent(inout) :: frestart
  character(LEN=255), intent(out) :: cdate

  character(LEN=255) fout

  integer year, month, day, hour, minute, second
  integer n

  logical newyear, newmonth, newday, newhour

  newyear  = is_newyear()
  newmonth = is_newmonth()
  newday   = is_newday()
  newhour  = is_newhour()

  fout = trim(foutdir)//'/'//trim(caseid)//'-colm-'

  call get_yrmndyhhmmss(idate_p,year,month,day,hour,minute,second)

  do n = 1, nMaxHist

     if(histArray(n)%is_valid) then
        histArray(n)%is_ready = .false.
        histArray(n)%is_newyear = newyear

        if(histArray(n)%yearly .and. newyear) then
           write(cdate,'(i4.4)') year
           histArray(n)%fhistory = trim(fout)//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%monthly .and. newmonth) then
           write(cdate,'(i4.4,"-",i2.2)') year,month
           histArray(n)%fhistory = trim(fout)//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%daily .and. newday) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') year,month,day
           histArray(n)%fhistory = trim(fout)//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%sixhourly .and. mod(hour+1,6).eq.0 .and. newhour) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i2.2)') year,month,day,hour
           histArray(n)%fhistory = trim(fout)//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%threehourly .and. mod(hour+1,3).eq.0 .and. newhour) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i2.2)') year,month,day,hour
           histArray(n)%fhistory = trim(fout)//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if

        if(histArray(n)%onehourly .and. newhour) then
           write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i2.2)') year,month,day,hour
           histArray(n)%fhistory = trim(fout)//trim(cdate)//'.nc'
           histArray(n)%is_ready = .true.
        end if
     end if

  end do

  write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1),idate(2),idate(3)

  if((trim(restart_freq).eq."yearly" .and. newyear) .or. &
     (trim(restart_freq).eq."monthly" .and. newmonth) .or. &
     (trim(restart_freq).eq."daily" .and. newday)) then
     frestart = trim(fout)//'restart-'//trim(cdate)
  endif

! Used as metadata in netcdf file
  write(cdate,'(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)') &
                year,    month,   day,     hour,    minute,  second

end subroutine lpwrite
! ------------------------------------------------------------------------
! EOP
