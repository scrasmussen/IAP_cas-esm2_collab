
      real time,dt
      integer date
     
      date = 13365
      time = 2359

      read(*,*) dt

      call uptime(time,date,dt)
     
      print*,time
      print*,date

      end





      subroutine uptime(time,date,deltat)
      implicit none
c
c----CAMx v5.40 111010
c
c     UPTIME increments simulation time/date for the given time step
c
c     Copyright 1996 - 2010
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        time            current simulation time (HHMMM)
c        date            current simulation date (YYJJJ)
c        deltat          current model coarse grid timestep (s)
c
c     Output arguments:
c        time            new simulation time (HHMMM)
c        date            new simulation date (YYJJJ)
c
c     Routines called:
c        none
c
c     Called by:
c        CAMx
c
      integer date
      real    tmin
      real    deltat
c
      real time, hour
      real fuzz
c
c-----Entry point
c
      hour = aint(time/100.)
      tmin = time - 100.*hour
      tmin = tmin + deltat/60.
      fuzz = max(0.02,0.1*(deltat/60.))
      if (tmin.ge.60.-fuzz .and. tmin.le.60.+fuzz) then
        tmin = 0.
        hour = hour + 1.
        if (hour.ge.24.) then
          hour = hour - 24.
          date = date + 1
          if( MOD(date,1000) .GT. 365 ) then
             if( MOD(INT(date/1000),4) .EQ. 0 ) then
                if( MOD(date,1000) .EQ. 367 )
     &                     date = (INT(date/1000)+1)*1000 + 1
             else
                date = (INT(date/1000)+1)*1000 + 1
             endif
          endif
        endif
      endif
c
      time = 100.*hour + tmin
c
      return
      end
