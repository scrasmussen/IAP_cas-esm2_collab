      subroutine diffus_mark(myid,ig,sx,ex,sy,ey,nlay,deltat,
     &                  rkv,depth,tempk,press,
     &                  concmark,atm,ismMax,smconv,kktop)
c
c----CAMx v4.40 061025
c
c     DIFFUS drives 3-D diffusion of concentrations.  This version also
c     performs diffusion of sensitivities if DDM is enabled.
c
c     Copyright 1996-2006
c     ENVIRON International Corporation
c          
c     Modifications:
c        4/17/00   Revised diffusion equations to weight fluxes by density
c       12/07/01   added instructions for OMP
c        1/13/03   added deposited mass array
c       10/12/04   Multiple substeps as f(Kv) applied for vertical diffusion 
c
c     Input arguments:
c        igrd              grid index
c        ncol              number of columns
c        nrow              number of rows
c        nlay              number of layers
c        nspc              number of species
c        nsen              number of species times number of DDM parameters
c                          or 1, whichever is larger
c        deltat            time step (s)
c        dx                cell size in x-direction (m)
c        dy                cell size in y-direction (m)
c        idfin             map of nested grids in this grid
c        vdep              deposition velocity (m/s)
c        rkx/y             horizontal diffusion coefficient (m2/s)
c        rkv               vertical diffusion coefficient (m2/s)
c        depth             layer depth (m)
c        tempk             temperature (K)
c        press             pressure (mb)
c        mapscl            map-scale factor at cell centroid
c        conc              species concentrations (umol/m3)
c        sens              DDM sensitivities (umol/m3/parameter unit)
c        tarray2           CPU timing arguments (s)
c        strz              string for labeling the z diffusion process
c        strxy             string for labeling the x/y diffusion process
c        ipa_xy            2-D gridded array to identify if cell is
c                          in a IPRM sub-domain
c        ipa_lay           3-D gridded array to identify which IPRM sub-domain
c                          each layer is in
cc        nparddm           number of parameters for which sensitivities
c                          are calculated.  Should be zero if DDM is not
c                          enabled.
c       smconv              the fraction of every tracer
c
c     Output arguments:
c        concmark              species concentrations (umol/m3)
c        sens              DDM sensitivities (umol/m3/parameter unit)
c        fluxes            fluxes across the boundaries (umol)
c        depfld            2-D array of dry deposited mass (mol/ha, g/ha)
c
c     Routines Called:
c        VDIFFIMP
c
c     Called by:
c        EMISTRNS
c
c
c
                    
c======================== Process Analysis Begin ====================================
c
c
      real fcup(nlay),fcdn(nlay)
      logical ldoipts
      integer :: nddmsp
c
c========================= Process Analysis End =====================================
c
      integer :: sx,ex,sy,ey
      integer :: MXTRSP      
      parameter(MXTRSP=0) ! the sensitity here no use
      dimension concmark(sx-1:ex+1,sy-1:ey+1,nlay),
     &    vdep(sx-1:ex+1,sy-1:ey+1),kktop(sx-1:ex+1,sy-1:ey+1)
      real rkx(sx-1:ex+1,sy-1:ey+1,nlay),rky(sx-1:ex+1,sy-1:ey+1,nlay),
     &    rkv(sx-1:ex+1,sy-1:ey+1,nlay),depth(sx-1:ex+1,sy-1:ey+1,nlay),
     &  tempk(sx-1:ex+1,sy-1:ey+1,nlay),press(sx-1:ex+1,sy-1:ey+1,nlay),
     &     mapscl(sx-1:ex+1,sy-1:ey+1),dx(sx-1:ex+1,sy-1:ey+1,nlay),
     &     atm(sx-1:ex+1,sy-1:ey+1,nlay),dair(sx-1:ex+1,sy-1:ey+1,nlay)
      dimension c1d(nlay+nlay*MXTRSP),d1d(nlay),rk1d(nlay),
     &          ro1d(nlay)
      integer nstepv(sx:ex,sy:ey)
      real dtv,dtmp
      real convfac
c=======================source appointment
      integer :: ismMax,myid
      real  ::  smconv(ismMax,sx-1:ex+1,sy-1:ey+1,nlay)
      real  ::  TMPsmconv(ismMax,nlay)
c============================
      mapscl=1. !lijie
      vdep=0.0  ! here no dry depositions
      nddmsp =0
!-----------convert ppbv to umol/m3---
      do i=sx,ex
       do j=sy,ey
        do k=1,nlay
!          if(ig==11.and.i==57.and.j==37.and.k==8)
!     &   print*,concmark(57,37,8),smconv(2,57,37,8),'mark begin'
               
         convfac=44.9*(273./tempk(i,j,k))*(press(i,j,k)/1013.)
          concmark(i,j,k)=concmark(i,j,k)*convfac/1e+03
          enddo
       enddo
      enddo  
c-----Entry point
c
c-----Vertical diffusion
c
c
c-----Determine max vertical diffusion time step and
c     apply within a sub-step loop
c
      do 601 j = sy,ey
        i1 = sx
        i2 = ex
        do 501 i = i1,i2
          dtv = deltat
          if (deltat.lt.300.) goto 401
          do k = 1,nlay-1
            dtmp = 0.75*depth(i,j,k)*depth(i,j,k+1)/rkv(i,j,k)
            dtv = amin1(dtv,dtmp)
          enddo
          dtv = amax1(dtv,300.)
 401      nstepv(i,j) = INT( 0.999*deltat/dtv ) + 1
 501    continue
 601  continue
c
c$omp parallel default(shared)
c$omp&  private(ispc,i1,i2,i,j,k,ro1d,d1d,c1d,rk1d,lk,
c$omp&          ioff,fluxbot,isen,rokp1,ldoipts,fcup,
c$omp&          fcdn,ipa_idx,n,dtv)
c
c$omp do schedule(dynamic)
c
        do 60 j = sy,ey
          i1 = sx
          i2 = ex
          do 50 i = i1,i2
c
c-----Skip cells occupied by child grids; load 1-D arrays
c
            dtv = deltat/FLOAT( nstepv(i,j) )
            do n = 1,nstepv(i,j)
            
c=======================source appointent
            DO  is=1,ismMax
            DO  k=1,nlay
              TMPsmconv(is,k)=smconv(is,i,j,k)
            ENDDO
            ENDDO
c==========================================
            do k = 1,nlay
              ro1d(k) = press(i,j,k)/tempk(i,j,k)
              d1d(k) = depth(i,j,k)
              c1d(k) = concmark(i,j,k)
            enddo
            do k = 1,nlay-1
              rokp1 = press(i,j,k+1)/tempk(i,j,k+1)
              rk1d(k) = rkv(i,j,k)*(ro1d(k) + rokp1)/2.
            enddo
            rk1d(nlay) = 0.
c
c
c======================== Process Analysis Begin ====================================
c
            ldoipts = .TRUE.
c
c========================= Process Analysis End =====================================
c
            call vdiffimp_mark(nlay,dtv,vdep(i,j),d1d,ro1d,rk1d,c1d,
     &     nddmsp,fcup,fcdn,ldoipts,ismMax,TMPsmconv,kktop(i,j),i,j)
c
c
c
            do k = 1,nlay
              concmark(i,j,k) = c1d(k)
             enddo

            DO  is=1,ismMax
            DO  k=1,nlay
              smconv(is,i,j,k)=TMPsmconv(is,k)
            ENDDO         
            ENDDO 
            enddo ! nstepv
  50      continue
  60    continue
c
!-----------convert umol/m3 to ppb---
      do i=sx,ex
       do j=sy,ey
         do k=1,nlay
         convfac=44.9*(273./tempk(i,j,k))*(press(i,j,k)/1013.)
         concmark(i,j,k)=1e03*concmark(i,j,k)/convfac
!        if(ig==11.and.i==57.and.j==37.and.k==8)
!     & print*,concmark(57,37,8),smconv(2,57,37,8),'mark after'
                       
         enddo
        enddo
      enddo
      return
      end
