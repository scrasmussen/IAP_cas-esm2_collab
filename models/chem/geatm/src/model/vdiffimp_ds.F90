      subroutine vdiffimp_ds(nn,dt,vdep,depth,rho,rkv,rr,nparddm, &
                         fcup,fcdn,ldoipts)
!
!----CAMx v4.40 061025
!
!     VDIFFIMP performs vertical diffusion of concentrations using 
!     an implicit method, where a tri-diagonal matrix is solved.
!     This version also performs vertical diffusion of sensitivities
!     if DDM is enabled.
!
!     Copyright 1996-2006
!     ENVIRON International Corporation
!          
!     Modifications:
!        4/17/00   Revised diffusion equations to weight fluxes by density
!
!     Input arguments:
!        nn                number of layers
!        dt                time step (s)
!        vdep              deposition velocity (m/s)
!        depth             layer depth (m)
!        rho               atmospheric density (mb/K)
!        rkv               vertical diffusion coefficient (mb m2/s/K)
!        rr                species concentrations (umol/m3) followed by 
!                          sensitivities (umol/m3/parameter unit).
!        nparddm           number of parameters for which sensitivities
!                          are calculated.  Should be zero if DDM is not
!                          enabled.
!        ldoipts           flag to calculate data for Process Analysis
!
!     Output arguments:
!        rr                species concentrations (umol/m3) followed by
!                          sensitivities (umol/m3/parameter unit)
!        fcup              change in layer concentration due to flux across
!                          upper interface (umol/m3) -- FOR Process Analysis
!        fcdn              change in layer concentration due to flux across
!                          lower interface (umol/m3) -- FOR Process Analysis
!
!     Routines Called:
!        TRDIAG
!
!     Called by:
!        DIFFUS
!
!
!======================== Process Analysis Begin ====================================
!
!
      real aa_old(nn),cc_old(nn)
      real fcup(nn),fcdn(nn)
      logical ldoipts
!
!========================= Process Analysis End =====================================
!
      integer :: nparddm
      dimension rr(nn+nn*nparddm),depth(nn),rkv(nn),rho(nn)
      dimension aa(nn),bb(nn),cc(nn)
!
!-----Entry point
!
!-----Lower boundary condition
!
      aa(1) = 0.
      bb(1) = 1. + dt/depth(1)* &
                  (vdep + 2.*rkv(1)/(depth(2)+depth(1))/rho(1))
!
!-----Upper boundary condition
!
      cc(nn) = 0.
      bb(nn) = 1. + dt/depth(nn)*&
                   2.*rkv(nn-1)/(depth(nn-1)+depth(nn))/rho(nn)
!
      do k = 2,nn
        aa(k) = -dt/depth(k)*2.*rkv(k-1)/(depth(k-1)+depth(k))/rho(k-1)
      enddo
      do k = 1,nn-1
        cc(k) = -dt/depth(k)*2.*rkv(k)/(depth(k+1)+depth(k))/rho(k+1)
      enddo
      do k = 2,nn-1
        bb(k) = 1. &
             + dt/depth(k)*2.*rkv(k-1)/(depth(k-1)+depth(k))/rho(k) &
             + dt/depth(k)*2.*rkv(k)/(depth(k+1)+depth(k))/rho(k)
      enddo
!
!======================== Process Analysis Begin ====================================
!
      if (ldoipts) then
        do k = 1,nn
          aa_old(k) = 0.
          cc_old(k) = 0.
          if (k.gt.1)  aa_old(k) = aa(k)*rho(k-1)
          if (k.lt.nn) cc_old(k) = cc(k)*rho(k+1)
        enddo
      endif
! 
!========================= Process Analysis End =====================================
!
!
!-----Solve the equations
!
      call trdiag_ds(aa,bb,cc,rr,nn,1+nparddm)
             
!
!
!======================== Process Analysis Begin ====================================
!
      if (ldoipts) then
        do k = 2,nn-1
           fcup(k) = (rr(k)/rho(k) - rr(k+1)/rho(k+1))*cc_old(k)
           fcdn(k) = (rr(k)/rho(k) - rr(k-1)/rho(k-1))*aa_old(k)
        enddo
!
!-----Lower boundary
!
        fcup(1) = (rr(1)/rho(1) - rr(2)/rho(2))*cc_old(1)
        fcdn(1) = -rr(1)*dt/depth(1)*vdep
!
!-----Upper boundary
!
        fcup(nn) =  0.0
        fcdn(nn) =  (rr(nn)/rho(nn) - rr(nn-1)/rho(nn-1))*aa_old(nn)
      endif
! 
!========================= Process Analysis End =====================================
!
      return
      end
