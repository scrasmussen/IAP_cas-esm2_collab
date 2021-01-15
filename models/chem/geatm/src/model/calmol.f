      subroutine calmol(temp,temp0,press,press0,deltaz,wind,landuse
     &                    ,el,i,j)
c 
c-----CAMx v4.11s 041206
c 
c     MICROMET calculates surface layer micro-meteorological flux-gradient
c     relationships and variables based on Louis (1979)
c  
c     Copyright 1996-2004
c     ENVIRON International Corporation
c  
c     Modifications:
c        8/15/02    Specified limits for denominators from 1e-10 to 1e-20
c        3/21/03    Modified Z/L calculation to use input midpoint height
c                   (instead of default 10 m) and to increase range from 
c                   (-1 to +1) to (-2.5 to +1.5) 
c             
c     Input arguments:  
c        temp                Layer 1 temperature (K) 
c        temp0               surface temperature (K)
c        press               Layer 1 pressure (mb)
c        press0              Surface layer pressure (mb)
c        deltaz              Layer 1  midpoint height (m)
c        wind                Layer 1 total wind speed (m/s)
c        z0                  Surface roughness length (m)
c             
c     Output arguments:  
c        ustar               friction velocity (m/s)
c        psih                Similarity stability correction term
c             
c     Routines called:  
c        none  
c             
c     Called by:  
c        DRYDEP 
c 
      real z0lu(24)  !lijie add from camx surface routhness(m) in landuse
      real landuse
      data vk/0.4/, g/9.8/, gamma/0.286/
      data z0lu/3., 0.25, 0.25, 0.25, 0.25, 0.1,
     &        0.05, 0.05, 0.05, 0.1, 1., 1.,
     &        1.,     1., 1., 0.0001,0.15,1.,
     &        0.02, 0.1, 1., 0.1,  0.03,0.1/
c
c-----Entry point
c-----to translate the 
      ilanduse=landuse
c
c-----Calculate potential temperature and richardson number
c
      z0    = z0lu(ilanduse) ! to get surface roughness (m) lijie
      theta = temp*(1000./press)**gamma
      theta0 = temp0*(1000./press0)**gamma
      dtheta = theta - theta0
      thetabar = (theta + theta0)/2.
      ri = (g/thetabar)*deltaz*dtheta/(wind**2 + 1.e-20)
c
c-----Determine stability functions
c
      zscale = vk/alog(deltaz/z0)
      if (ri.lt.0.) then
        cm    = 69.56*sqrt(deltaz/z0)*zscale**2
        ch    = 49.82*sqrt(deltaz/z0)*zscale**2
        fm    = 1. - 9.4*ri/(1. + cm*sqrt(abs(ri)))
        fh    = 1. - 9.4*ri/(1. + ch*sqrt(abs(ri)))
      else
        fm = 1./((1. + 4.7*ri)**2)
        fh = fm
      endif
c
c-----Calculate micromet variables
c
      ustar2 = fm*(wind*zscale)**2
      ustar2 = amax1(1.e-20,ustar2)
      ustar = sqrt(ustar2)
      thstar = 1.35*zscale**2*wind*dtheta*fh/ustar
      el = ustar2*temp/(vk*g*thstar + 1.e-20)
      elabs = abs(el) 
      if (elabs.ge.1.e4) then 
        psih = 0.0 
      elseif (el.lt.0.) then 
        zoverl = amin1(2.5,deltaz/elabs)
        zmel = alog(zoverl) 
        psih = exp(0.598 + 0.39*zmel - 0.090*zmel*zmel) 
      else    
        zoverl = amin1(1.5,deltaz/elabs)
        psih = -5.*zoverl
      endif
c
      return
      end
