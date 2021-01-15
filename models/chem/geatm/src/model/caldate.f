      subroutine caldate(idate)
c  
c----CAMx v5.40 111010
c  
c     CALDATE converts date from Julian (YYJJJ) format to calender
c     (YYMMDD) format
c                            
c     Copyright 1996 - 2010
c     ENVIRON International Corporation
c            
c     Modifications:  
c        none
c  
c     Input arguments:  
c        idate               julian date (YYJJJ)  
c              
c     Output arguments:  
c        idate               calender date (YYMMDD)  
c              
c     Routines Called:  
c        none
c              
c     Called by:  
c        DRYDEP
c        DRYDEPRT
c        CHRTIME
c        PIGDRIVE
c
      integer nday(12)
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----Entry point
c
c-----If it is already in calender date, return
c
      if (idate.gt.100000) goto 9999
      iyear = idate/1000
      jday = idate - iyear*1000
c
      nday(2) = 28
      if (mod(iyear,4).eq.0) nday(2) = 29
      mday = 0
      do 10 imonth = 1,12
        mday = mday + nday(imonth)
        if (mday.ge.jday) go to 20
 10   continue
 20   iday = jday - (mday - nday(imonth))
      idate = iyear*10000 + imonth*100 + iday
c
 9999 return
      end
