      subroutine diffus_ds(myid,ig,sx,ex,sy,ey,nlay,deltat,&
                       rkv,depth,tempk,press,&
                       conc,concom,atm,kktop,nspecom,&
                       IS,IA)
!
!----CAMx v4.40 061025
!
!     DIFFUS drives 3-D diffusion of concentrations.  This version also
!     performs diffusion of sensitivities if DDM is enabled.
!
!     Copyright 1996-2006
!     ENVIRON International Corporation
!          
!     Modifications:
!        4/17/00   Revised diffusion equations to weight fluxes by density
!       12/07/01   added instructions for OMP
!        1/13/03   added deposited mass array
!       10/12/04   Multiple substeps as f(Kv) applied for vertical diffusion 
!
!     Input arguments:
!        igrd              grid index
!        ncol              number of columns
!        nrow              number of rows
!        nlay              number of layers
!        nspc              number of species
!        nsen              number of species times number of DDM parameters
!                          or 1, whichever is larger
!        deltat            time step (s)
!        dx                cell size in x-direction (m)
!        dy                cell size in y-direction (m)
!        idfin             map of nested grids in this grid
!        vdep              deposition velocity (m/s)
!        rkx/y             horizontal diffusion coefficient (m2/s)
!        rkv               vertical diffusion coefficient (m2/s)
!        depth             layer depth (m)
!        tempk             temperature (K)
!        press             pressure (mb)
!        mapscl            map-scale factor at cell centroid
!        conc              species concentrations (umol/m3)
!        sens              DDM sensitivities (umol/m3/parameter unit)
!        tarray2           CPU timing arguments (s)
!        strz              string for labeling the z diffusion process
!        strxy             string for labeling the x/y diffusion process
!        ipa_xy            2-D gridded array to identify if cell is
!                          in a IPRM sub-domain
!        ipa_lay           3-D gridded array to identify which IPRM sub-domain
!                          each layer is in
!        nparddm           number of parameters for which sensitivities
!                          are calculated.  Should be zero if DDM is not
!                          enabled.
!
!     Output arguments:
!        conc              species concentrations (umol/m3)
!        sens              DDM sensitivities (umol/m3/parameter unit)
!        fluxes            fluxes across the boundaries (umol)
!        depfld            2-D array of dry deposited mass (mol/ha, g/ha)
!
!     Routines Called:
!        VDIFFIMP
!
!     Called by:
!        EMISTRNS
!
!
!
                    
!======================== Process Analysis Begin ====================================
!
!
      real fcup(nlay),fcdn(nlay)
      logical ldoipts
      integer :: nddmsp

!========================= Process Analysis End =====================================
!
      integer :: sx,ex,sy,ey,nspecom
      integer :: MXTRSP      
      parameter(MXTRSP=0) ! the sensitity here no use
      dimension conc(sx-1:ex+1,sy-1:ey+1,nlay),vdep(sx-1:ex+1,sy-1:ey+1) &
      ,kktop(sx-1:ex+1,sy-1:ey+1),concom(sx-1:ex+1,sy-1:ey+1,nlay,nspecom)
      real rkx(sx-1:ex+1,sy-1:ey+1,nlay),rky(sx-1:ex+1,sy-1:ey+1,nlay),&
      rkv(sx-1:ex+1,sy-1:ey+1,nlay),depth(sx-1:ex+1,sy-1:ey+1,nlay),&
       tempk(sx-1:ex+1,sy-1:ey+1,nlay),press(sx-1:ex+1,sy-1:ey+1,nlay),&
          mapscl(sx-1:ex+1,sy-1:ey+1),dx(sx-1:ex+1,sy-1:ey+1,nlay),&
          atm(sx-1:ex+1,sy-1:ey+1,nlay),dair(sx-1:ex+1,sy-1:ey+1,nlay)
     
      dimension c1d(nlay+nlay*MXTRSP),d1d(nlay),rk1d(nlay),ro1d(nlay)
      integer nstepv(sx:ex,sy:ey)
      real dtv,dtmp
      real convfac
      real ratio
      mapscl=1. !lijie
      vdep=0.0  ! here no dry depositions
      nddmsp =0
!-----Entry point
!
!-----Vertical diffusion
!
!
!-----Determine max vertical diffusion time step and
!     apply within a sub-step loop
!
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
 401      nstepv(i,j) = INT( 0.999*deltat/dtv ) + 10
 501    continue
 601  continue
!
!
        do 60 j = sy,ey
          i1 = sx
          i2 = ex
          do 50 i = i1,i2
!
!-----Skip cells occupied by child grids; load 1-D arrays
!
            dtv = deltat/FLOAT( nstepv(i,j) )
            do n = 1,nstepv(i,j)

            do k = 1,nlay
              ro1d(k) = press(i,j,k)/tempk(i,j,k)
              d1d(k) = depth(i,j,k)
              c1d(k) = conc(i,j,k)
            enddo
            do k = 1,nlay-1
              rokp1 = press(i,j,k+1)/tempk(i,j,k+1)
              rk1d(k) = rkv(i,j,k)*(ro1d(k) + rokp1)/2.
            enddo
            rk1d(nlay) = 0.

!
!======================== Process Analysis Begin ====================================
!
            ldoipts = .TRUE.
!
!========================= Process Analysis End =====================================
!
            call vdiffimp_ds(nlay,dtv,vdep(i,j),d1d,ro1d,rk1d,c1d,   &
                    nddmsp,fcup,fcdn,ldoipts)

            do k = 1,nlay-1

              if(c1d(k)-conc(i,j,k).ne.0) then
               ratio = (fcup(k)+fcdn(k))/(c1d(k)-conc(i,j,k))          
              else
               ratio = 1.0
              endif

              IF(ratio==0.and. (c1d(k)-conc(i,j,k)).ne.0.) THEN
                fcup(k) = c1d(k)-conc(i,j,k)
                fcdn(k) = 0.0
                ratio = 1.
              ENDIF                 

              fcup(k)= fcup(k) /ratio
              fcdn(k)= fcdn(k) / ratio

              conc_tmp0 = conc(i,j,k)
              conc_tmp1 = conc(i,j,k) + fcup(k)
              conc_tmp2 = conc(i,j,k) + fcup(k) + fcdn(k)

              conc(i,j,k) = amax1(conc(i,j,k),1.E-20)


              do iduc = 1,nspecom
                concom(i,j,k,iduc) = amax1(concom(i,j,k,iduc),1.E-20)


               IF(fcup(k).GT.0.) THEN
                 IF(fcup(k)>conc(i,j,k+1)) GOTO 100
                 conc(i,j,k+1) = amax1(conc(i,j,k+1),1.E-20)
                 concom(i,j,k+1,iduc)=amax1(concom(i,j,k+1,iduc),1.E-20)

                concom(i,j,k,iduc) = concom(i,j,k,iduc) +   &
                        fcup(k)*concom(i,j,k+1,iduc)/conc(i,j,k+1)

               ELSE
                  IF(ABS(fcup(k))>conc_tmp0 ) GOTO 100
                concom(i,j,k,iduc) = concom(i,j,k,iduc) +  &
                         fcup(k)*concom(i,j,k,iduc)/conc_tmp0


               ENDIF

100   CONTINUE
               IF(fcdn(k).GT.0.) THEN
                if(K==1) THEN
                  conc_1=conc(i,j,1)
                  concom_1 = concom(i,j,1,iduc)
                else
                  conc_1 = conc(i,j,k-1)
                  concom_1 = concom(i,j,k-1,iduc)
                endif
               
                IF(fcdn(k)>conc_1 ) GOTO 101
                 conc_1 = amax1(conc_1,1.E-20)           
                 concom_1 = amax1(concom_1,1.E-20)
                 concom(i,j,k,iduc) = concom(i,j,k,iduc) +&
                      fcdn(k)*concom_1/conc_1
               ELSE
                 
                 IF(ABS(fcdn(k))>conc_tmp0) GOTO 101
                 concom(i,j,k,iduc) = concom(i,j,k,iduc) +&
                              fcdn(k)*concom(i,j,k,iduc)/conc_tmp1
               ENDIF
                 concom(i,j,k,iduc) = amax1(concom(i,j,k,iduc),1.E-21)
101    CONTINUE

              enddo ! iduc

            enddo ! k                                      

              do k = 1,nlay
               conc(i,j,k) = c1d(k)
               conc(i,j,k) = amax1(conc(i,j,k),1.E-20)
              enddo ! k

            enddo ! nstepv
  50      continue
  60    continue
!
!-----------convert umol/m3 to ppb---
      do i=sx,ex
       do j=sy,ey
         do k=1,nlay
          conc(i,j,k)=amax1(conc(i,j,k),0.0)
         enddo
        enddo
      enddo
      return
      end
