      subroutine vdiffimp_mark(nn,dt,vdep,depth,rho,rkv,rr,nparddm,
     &                    fcup,fcdn,ldoipts,ismMax,TMPsmconv,kktop,i,j)
c
c----CAMx v4.40 061025
c
c     VDIFFIMP performs vertical diffusion of concentrations using 
c     an implicit method, where a tri-diagonal matrix is solved.
c     This version also performs vertical diffusion of sensitivities
c     if DDM is enabled.
c
c     Copyright 1996-2006
c     ENVIRON International Corporation
c          
c     Modifications:
c        4/17/00   Revised diffusion equations to weight fluxes by density
c
c     Input arguments:
c        nn                number of layers
c        dt                time step (s)
c        vdep              deposition velocity (m/s)
c        depth             layer depth (m)
c        rho               atmospheric density (mb/K)
c        rkv               vertical diffusion coefficient (mb m2/s/K)
c        rr                species concentrations (umol/m3) followed by 
c                          sensitivities (umol/m3/parameter unit).
c        nparddm           number of parameters for which sensitivities
c                          are calculated.  Should be zero if DDM is not
c                          enabled.
c        ldoipts           flag to calculate data for Process Analysis
c
c     Output arguments:
c        rr                species concentrations (umol/m3) followed by
c                          sensitivities (umol/m3/parameter unit)
c        fcup              change in layer concentration due to flux across
c                          upper interface (umol/m3) -- FOR Process Analysis
c        fcdn              change in layer concentration due to flux across
c                          lower interface (umol/m3) -- FOR Process Analysis
c
c     Routines Called:
c        TRDIAG
c
c     Called by:
c        DIFFUS
c
c
c======================== Process Analysis Begin ====================================
c
c
      real aa_old(nn),cc_old(nn)
      real fcup(nn),fcdn(nn)
      logical ldoipts
c
c========================= Process Analysis End =====================================
c
      integer :: nparddm
      dimension rr(nn+nn*nparddm),depth(nn),rkv(nn),rho(nn)
      dimension aa(nn),bb(nn),cc(nn),rrtmp(nn+nn*nparddm)
c
c==============source appointment
        integer,parameter :: ismMaxSet=100
        real smthis(ismMaxSet)
        real smother(ismMaxSet)
      integer :: ismMax,kktop
      real  ::  TMPsmconv(ismMax,nn)
c===============================      
c-----Entry point
c
c-----Lower boundary condition
c
      aa(1) = 0.
      bb(1) = 1. + dt/depth(1)*
     &             (vdep + 2.*rkv(1)/(depth(2)+depth(1))/rho(1))
c
c-----Upper boundary condition
c
      cc(nn) = 0.
      bb(nn) = 1. + dt/depth(nn)*
     &              2.*rkv(nn-1)/(depth(nn-1)+depth(nn))/rho(nn)
c
      do k = 2,nn
        aa(k) = -dt/depth(k)*2.*rkv(k-1)/(depth(k-1)+depth(k))/rho(k-1)
      enddo
      do k = 1,nn-1
        cc(k) = -dt/depth(k)*2.*rkv(k)/(depth(k+1)+depth(k))/rho(k+1)
      enddo
      do k = 2,nn-1
        bb(k) = 1.
     &        + dt/depth(k)*2.*rkv(k-1)/(depth(k-1)+depth(k))/rho(k)
     &        + dt/depth(k)*2.*rkv(k)/(depth(k+1)+depth(k))/rho(k)
      enddo
c
c======================== Process Analysis Begin ====================================
c
      if (ldoipts) then
        do k = 1,nn
          aa_old(k) = 0.
          cc_old(k) = 0.
          if (k.gt.1)  aa_old(k) = aa(k)*rho(k-1)
          if (k.lt.nn) cc_old(k) = cc(k)*rho(k+1)
        enddo
      endif
c 
c========================= Process Analysis End =====================================
c=======================source appoint=======
       do k=1,nn
         rrtmp(k)=rr(k)
       enddo
c=========================       
c   
c
c-----Solve the equations
c
      call trdiag_mark(aa,bb,cc,rr,nn,1+nparddm)
c
c
c======================== Process Analysis Begin ====================================
c
      if (ldoipts) then
        do k = 2,nn-1
           fcup(k) = (rr(k)/rho(k) - rr(k+1)/rho(k+1))*cc_old(k)
           fcdn(k) = (rr(k)/rho(k) - rr(k-1)/rho(k-1))*aa_old(k)
        enddo
c
c-----Lower boundary
c
        fcup(1) = (rr(1)/rho(1) - rr(2)/rho(2))*cc_old(1)
        fcdn(1) = -rr(1)*dt/depth(1)*vdep
c
c-----Upper boundary
c
        fcup(nn) =  0.0
        fcdn(nn) =  (rr(nn)/rho(nn) - rr(nn-1)/rho(nn-1))*aa_old(nn)
      endif
c 
c========================= Process Analysis End =====================================
c
c========================source appointment========      
      do k=1,nn
       IF(fcup(k)>=0) then
        do is=1,ismMax
          smthis(is)=TMPsmconv(is,k)
         if(k==nn) then
            smother(is)=TMPsmconv(is,k)
         else        
            smother(is)=TMPsmconv(is,k+1)
         endif   
        enddo 
        deltc1=fcup(k)
        deltc1=max(fcup(k),1.e-20)
        rrtmp(k) = max(rrtmp(k),1.e-20)       


        if(k==kktop)then
         call GetBounChange(deltc1,rrtmp(k),smthis,2,ismMax)
        else if(k>kktop) then
         do is =1 ,ismMax
         if(is==2) then
          smthis(is)=1.0
         else
          smthis(is) =0.0
         endif
         enddo
        else
         call SMmixing(rrtmp(k),smthis,deltc1,smother,ismMax)
        endif
        
        do is=1,ismMax
         TMPsmconv(is,k)=smthis(is)
        enddo

       ENDIF !fcup

      IF(fcdn(k)>=0)then
        do is=1,ismMax
         smthis(is)=TMPsmconv(is,k)
         if(k==1) then
        smother(is)=TMPsmconv(is,k)
         else
         smother(is)=TMPsmconv(is,k-1)
         endif
        enddo      
        deltc2=fcdn(k)
        deltc2=max(fcdn(k), 1.e-20)
        rrtmp(k) = max (rrtmp(k), 1.e-20)
     
        
       if(k>kktop) then
         do is = 1 , ismMax
          if(is==2) then
           smthis(2)=1.0
          else
           smthis(2) = 0.0
          endif
         enddo 
       else
         call SMmixing(rrtmp(k)+fcup(k),smthis,deltc2,smother,ismMax)
       endif       
 
       do is=1,ismMax
         TMPsmconv(is,k)=smthis(is)
       enddo
      ENDIF !fcdn
      enddo !k
      return
      end
