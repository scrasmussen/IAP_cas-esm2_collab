      subroutine wrf_micromet(temp,temp0,press,press0,deltaz,wind,z0,pbl, &
                          ustar,eli,wstar,ri)
! 
!     MICROMET calculates surface layer micro-meteorological flux-gradient
!     relationships and variables based on Louis (1979)
!  
!     Input arguments:  
!        temp                Layer 1 temperature (K) 
!        temp0               Surface temperature (K)
!        press               Layer 1 pressure (mb)
!        press0              Surface pressure (mb)
!        deltaz              Layer 1 midpoint height (m)
!        wind                Layer 1 total wind speed (m/s)
!        z0                  Surface roughness length (m)
!        pbl                 PBL depth (m)
!             
!     Output arguments:  
!        ustar               friction velocity (m/s)
!        eli                 inverse M-O length (1/m)
!        wstar               convective velocity scale (m/s)
!        ri                  richardson number
!             
      implicit none
!
      real temp,temp0,press,press0,deltaz,wind,z0,pbl,ustar,eli,wstar
      real theta,theta0,dtheta,thetabar,ri,zscale,cm,ch,fm,fh, &
           ustar2,thstar,el,elabs,zz,xl
      real vk,g,gamma,xinf,elimin
      data vk/0.4/, g/9.8/, gamma/0.286/, xinf/100./, elimin/1.e-4/
!
!-----Entry point
!
!-----Calculate potential temperature and richardson number
!
      theta = temp*(1000./press)**gamma
      theta0 = temp0*(1000./press0)**gamma
      dtheta = theta - theta0
      thetabar = (theta + theta0)/2.
      ri = (g/thetabar)*deltaz*dtheta/(wind**2 + 1.e-20)
!
!-----Determine stability functions
!
      zscale = vk/alog(deltaz/z0)
      if (ri.le.0.) then
        cm    = 69.56*sqrt(deltaz/z0)*zscale**2
        ch    = 49.82*sqrt(deltaz/z0)*zscale**2
        fm    = 1. - 9.4*ri/(1. + cm*sqrt(abs(ri)))
        fh    = 1. - 9.4*ri/(1. + ch*sqrt(abs(ri)))
      else
        fm = 1./((1. + 4.7*ri)**2)
        fh = fm
      endif
!
!-----Calculate micromet variables
!
      ustar2 = fm*(wind*zscale)**2
      ustar2 = amax1(1.e-20,ustar2)
      ustar = sqrt(ustar2)
      thstar = 1.35*zscale**2*wind*dtheta*fh/ustar
      el = ustar2*temp/(vk*g*thstar + 1.e-20)
      eli = 1./el
      eli = sign(max(abs(eli),elimin),eli)
      el = 1./eli
      elabs = abs(el) 

      wstar = 0.
      if (el.lt.0.) wstar = (pbl*ustar**3./(vk*elabs))**(1./3.)

      return
      end
