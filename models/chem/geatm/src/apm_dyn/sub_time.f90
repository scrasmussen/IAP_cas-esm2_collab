subroutine int22char( flag,year,month,day,hour,minute,second,cyear,&
                   &  cmonth,cday,chour,cminute,csecond,&
                   &  timestr )
integer :: year,month,day,hour,minute,second
character*2 :: cyear*4,cmonth,cday,chour,cminute,csecond
character*8 :: flag
character*19 :: timestr
if(flag.eq.'int2char') then
 write(cyear,'(i4.4)') year
 write(cmonth,'(i2.2)') month
 write(cday,'(i2.2)') day
 write(chour,'(i2.2)') hour
 write(cminute,'(i2.2)') minute
 write(csecond,'(i2.2)') second
 timestr=cyear//'-'//cmonth//'-'//cday//'_'//chour//':'//cminute//':'//csecond
elseif(flag.eq.'char2int') then
 read(cyear,*) year
 read(cmonth,*) month
 read(cday,*) day
 read(chour,*) hour
 read(cminute,*) minute
 read(csecond,*) second
 timestr=cyear//'-'//cmonth//'-'//cday//'_'//chour//':'//cminute//':'//csecond
endif
end subroutine int22char


subroutine plus8h( delt_d,delt_h,delt_mini,delt_sec,  year,month,day,hour,mini,sec )
implicit none
integer :: delt_d,delt_h,delt_mini,year,month,day,hour,mini
integer :: delt_sec,sec
!real :: delt_d,delt_h,delt_mini,year,month,day,hour,mini
!real :: delt_sec,sec
logical :: leap_year
integer :: mdays(12),mdays1(12),mdays2(12)
data mdays1/31,28,31,30,31,30,31,31,30,31,30,31/
data mdays2/31,29,31,30,31,30,31,31,30,31,30,31/
call if_leap(year,leap_year)
if(leap_year) then
 mdays=mdays2
else
 mdays=mdays1
endif
if(delt_d.gt.28) then
 print*,'delt_d is too large'
 call abort()
endif

day=day+delt_d
if(day.gt.mdays(month)) then
 day=day-mdays(month)
 month=month+1
 if(month.gt.12) then
  month=month-1
  day=day+mdays(month)-mdays(12)
  month=1
  year=year+1
 endif
endif

hour=hour+delt_h
if(hour.ge.24) then
 hour=hour-24
 day=day+1
 if(day.gt.mdays(month)) then
  day=1
  month=month+1
  if(month.gt.12) then
   month=1
   year=year+1
  endif
 endif
endif

mini=mini+delt_mini
if(mini.ge.60) then
 hour=hour+1
 mini=mini-60
 if(hour.ge.24) then
  hour=hour-24
  day=day+1
  if(day.gt.mdays(month)) then
   day=1
   month=month+1
   if(month.gt.12) then
    month=1
    year=year+1
   endif
  endif
 endif
endif

sec=sec+delt_sec
IF(sec.ge.60) THEN
 mini=mini+1
 sec=sec-60
if(mini.ge.60) then
 hour=hour+1
 mini=mini-60
 if(hour.ge.24) then
  hour=hour-24
  day=day+1
  if(day.gt.mdays(month)) then
   day=1
   month=month+1
   if(month.gt.12) then
    month=1
    year=year+1
   endif
  endif
 endif
endif
ENDIF

end

subroutine checkinput(year,month,day,mdays)
implicit none
integer :: year,month,day
integer :: mdays(12)
if(day.gt.mdays(month)) then
 print*,'input erro'
 call abort()
endif
end

subroutine if_leap(year,leap_year)
integer :: year
logical :: leap_year
if((mod(year,400).eq.0).or.(mod(year,4).eq.0.and.mod(year,100).ne.0)) then
 leap_year=.true.
else
 leap_year=.false.
endif
end


