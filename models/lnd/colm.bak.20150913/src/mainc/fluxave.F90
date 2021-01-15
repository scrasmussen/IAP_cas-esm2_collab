
#include <define.h>

! ----------------------------------------------------------------------
!                             FLUX TABLE                               !
! ----------------------------------------------------------------------
! perfrom grid-average from subgrid 1d vector
! subgrid to grid average mapping: average a subgrid input vector [fldv] 
! of length numpatch to a output array [fldv] of length numgrid
!
! Created by Yongjiu Dai
!--------------!-------------------------------------------------------
! pft level fluxes
!--------------!-------------------------------------------------------
! 01: taux     ! wind stress: E-W [kg/m/s2]
! 02: tauy     ! wind stress: N-S [kg/m/s2]
! 03: fsena    ! sensible heat from canopy height to atmosphere [W/m2]
! 04: lfevpa   ! latent heat flux from canopy height to atmosphere [W/m2]
! 05: fevpa    ! evapotranspiration from canopy to atmosphere [mm/s]
! 06: fsenl    ! sensible heat from leaves [W/m2]
! 07: fevpl    ! evaporation+transpiration from leaves [mm/s]
! 08: etr      ! transpiration rate [mm/s]
! 09: fseng    ! sensible heat flux from ground [W/m2]
! 10: fevpg    ! evaporation heat flux from ground [mm/s]
! 11: fgrnd    ! ground heat flux [W/m2]
! 12: sabvsun  ! solar absorbed by sunlit canopy [W/m2]
! 13: sabvsha  ! solar absorbed by shaded [W/m2]
! 14: sabg     ! solar absorbed by ground [W/m2]
! 15: olrg     ! outgoing long-wave radiation from ground+canopy [W/m2]
! 16: rnet     ! net radiation [W/m2]
! 17: zerr     ! the error of energy balance [W/m2]
! 18: assim    ! canopy assimilation rate [mol m-2 s-1]
! 19: respc    ! respiration (plant+soil) [mol m-2 s-1]
! 20: fmicr    ! microbial respiration    [mol m-2 s-1]
! 21: tlsun    ! sunlit leaf temperature [K]
! 22: tlsha    ! shaded leaf temperature [K]
! 23: ldew     ! depth of water on foliage [mm]
! 24: sigf     ! fraction of veg cover, excluding snow-covered veg [-]
! 25: green    ! leaf greenness
! 26: lai      ! leaf area index
! 27: sai      ! stem area index
! 28: alb(1,1) ! averaged albedo [visible, direct]
! 29: alb(1,2) ! averaged albedo [visible, diffuse]
! 30: alb(2,1) ! averaged albedo [near-infrared, direct]
! 31: alb(2,2) ! averaged albedo [near-infrared,diffuse]
! 32: emis     ! averaged bulk surface emissivity
! 33: z0ma     ! effective roughness [m]
!--------------!-------------------------------------------------------
! 34: trad     ! radiative temperature of surface [K]
! 35: ustar    ! u* in similarity theory [m/s]
! 36: tstar    ! t* in similarity theory [kg/kg]
! 37: qstar    ! q* in similarity theory [kg/kg]
! 38: zol      ! dimensionless height (z/L) used in Monin-Obukhov theory
! 39: rib      ! bulk Richardson number in surface layer
! 40: fm       ! integral of profile function for momentum
! 41: fh       ! integral of profile function for heat
! 42: fq       ! integral of profile function for moisture
!--------------!-------------------------------------------------------
! 43: tref     ! 2 m height air temperature [kelvin]
! 44: qref     ! 2 m height air specific humidity [kg/kg]
! 45: u10m     ! 10m u-velocity [m/s]
! 46: v10m     ! 10m v-velocity [m/s]
! 47: f10m     ! integral of profile function for momentum at 10m [-]
!--------------!-------------------------------------------------------
! column level fluxes
!--------------!-------------------------------------------------------
! 48: xerr     ! the error of water banace [mm/s]
! 49: rsur     ! surface runoff [mm/s]
! 50: rnof     ! total runoff [mm/s]
!--------------!-------------------------------------------------------
! 51:60: tss   ! soil temperature [K]
! 61:70: wliq  ! liquid water in soil layers [kg/m2]
! 71:80: wice  ! ice lens in soil layers [kg/m2]
!--------------!-------------------------------------------------------
! 81: tg       ! ground surface temperature [K]
! 82: scv      ! snow cover, water equivalent [mm]
! 83: snowdp   ! snow depth [meter]
! 84: fsno     ! fraction of snow cover on ground
!----------------------------------------------------------------------
! 85: us       ! wind in eastward direction [m/s]
! 86: vs       ! wind in northward direction [m/s]
! 87: tm       ! temperature at reference height [kelvin]
! 88: qm       ! specific humidity at reference height [kg/kg]
! 89: prc      ! convective precipitation [mm/s]
! 90: prl      ! large scale precipitation [mm/s]
! 91: pbot     ! atmospheric pressure at the surface [pa]
! 92: frl      ! atmospheric infrared (longwave) radiation [W/m2]
! 93: solar    ! downward solar radiation at surface [W/m2]
!----------------------------------------------------------------------
! pft level fluxes
!----------------------------------------------------------------------
! 94:  qsubl        ! sublimation rate from snow pack [kg/m2/s]
!----------------------------------------------------------------------
! column level fluxes
!----------------------------------------------------------------------
! 95:104: mrlsl     ! mass of water of all phases in each soil layer [kg/m2]
! 105: mrsos        ! mass of water of all phases in the upper 0.1 meters of soil [kg/m2]
! 106: mrso         ! mass of water of all phases over all soil layers [kg/m2]
! 107: mrfso        ! mass of frozen water over all soil layers [kg/m2]
! 108: lwsnl        ! mass of liquid water of snow layers [kg/m2]
! 109: sm           ! surface snow melt [kg/m2/s]
! 110: tsn          ! snow internal temperature [K]
! 111: nsnow        ! number of snow events [-]
!----------------------------------------------------------------------
! 112: treeFrac     ! tree fraction [-]
! 113: shrubFrac    ! shrub fraction [-] 
! 114: grassFrac    ! grass fraction [-]
! 115: baresoilFrac ! bair soil fraction [-]
! 116: residualFrac ! residual fraction [-]
! 117: soilFrac     ! soil fraction [-]
! 118: urbanFrac    ! urban fraction [-]
! 119: wetlandFrac  ! wetland fraction [-]
! 120: iceFrac      ! ice fraction [-]
! 121: lakeFrac     ! lake & river fraction [-]
!----------------------------------------------------------------------

 subroutine flux_p2g(nump,numg,pgmap,wt,fp,fg)

      use precision
      implicit none

      integer ,intent(in) :: nump
      integer ,intent(in) :: numg
      integer ,intent(in) :: pgmap(nump)
      real(r8),intent(in) :: wt(nump)
      real(r8),intent(in) :: fp(nump)
      real(r8),intent(out):: fg(numg)

      real(r8) sumwt(numg)
      integer p, g

      fg   (:) = 0.
      sumwt(:) = 0.

      do p = 1, nump
         g = pgmap(p)
         fg(g) = fg(g) + wt(p)*fp(p)
         sumwt(g) = sumwt(g) + wt(p)
      end do

      do g = 1, numg
         if(sumwt(g).gt.0) then
            fg(g) = fg(g)/sumwt(g)
         else
            fg(g) = -9999.
         end if
      end do

 end subroutine flux_p2g

 subroutine flux_p2g_grass(nump,numg,pgmap,ivt,wt,fp,fg)

      use precision
      use paramodel, only: grasscateg, oceancateg
      implicit none

      integer ,intent(in) :: nump
      integer ,intent(in) :: numg
      integer ,intent(in) :: pgmap(nump)
      integer ,intent(in) :: ivt(nump)
      real(r8),intent(in) :: wt(nump)
      real(r8),intent(in) :: fp(nump)
      real(r8),intent(out):: fg(numg)

      real(r8) sumwt(numg)
      integer p, g

      fg   (:) = 0.
      sumwt(:) = 0.

      do p = 1, nump
         g = pgmap(p)
         fg(g) = fg(g) + wt(p)*fp(p)
         sumwt(g) = sumwt(g) + wt(p)
      end do

      do g = 1, numg
         if(sumwt(g).gt.0) then
            fg(g) = fg(g)/sumwt(g)
         else
            fg(g) = -9999.
         end if
      end do

      do p = 1, nump
         g = pgmap(p)
         if((ivt(p).eq.grasscateg .or. ivt(p).eq.oceancateg) .and. wt(p).gt.1.0E-2)then
            fg(g) = fp(p)
         endif 
      enddo

 end subroutine flux_p2g_grass

 subroutine flux_c2g(numc,numg,cgmap,wt,fc,fg)

      use precision
      implicit none

      integer ,intent(in) :: numc
      integer ,intent(in) :: numg
      integer ,intent(in) :: cgmap(numc)
      real(r8),intent(in) :: wt(numc)
      real(r8),intent(in) :: fc(numc)
      real(r8),intent(out):: fg(numg)

      real(r8) sumwt(numg)
      integer c, g

      fg   (:) = 0.
      sumwt(:) = 0.

      do c = 1, numc
         g = cgmap(c)
         fg(g) = fg(g) + wt(c)*fc(c)
         sumwt(g) = sumwt(g) + wt(c)
      end do

      do g = 1, numg
         if(sumwt(g).gt.0) then
            fg(g) = fg(g)/sumwt(g)
         else
            fg(g) = -9999.
         end if
      end do

 end subroutine flux_c2g

 subroutine flux_c2g_3d(numc,numg,nlev,cgmap,itypwat,wt,fc,fg)

      use precision
      implicit none

      integer ,intent(in) :: numc
      integer ,intent(in) :: numg
      integer ,intent(in) :: nlev
      integer ,intent(in) :: cgmap(numc)
      integer ,intent(in) :: itypwat(numc)
      real(r8),intent(in) :: wt(numc)
      real(r8),intent(in) :: fc(nlev,numc)
      real(r8),intent(out):: fg(nlev,numg)

      logical  iswater(numg) ! full grid is lake or ocean
      real(r8) sumwt(numg)
      integer c, g

      fg (:,:) = 0.
      sumwt(:) = 0.

      iswater(:) = .true.

      do c = 1, numc
         g = cgmap(c)
         if(itypwat(c).le.3) iswater(g) = .false.
      end do

      do c = 1, numc
         g = cgmap(c)
         if(itypwat(c).le.3) then
            fg(:,g) = fg(:,g) + wt(c)*fc(:,c)
            sumwt(g) = sumwt(g) + wt(c)
         end if
      end do

      do g = 1, numg
         if(sumwt(g).gt.0) then
            fg(:,g) = fg(:,g)/sumwt(g)
         else
            fg(:,g) = -9999.
         end if
      end do

 end subroutine flux_c2g_3d

 subroutine fluxave

      use precision
      use phycon_module
      use paramodel
      use spmd
      use spmd_decomp
      use colm_varMod

      implicit none

! local variables

      integer  :: c,p,g                               ! indices
      integer  :: ivt(numpatch)                       ! veg type
      real(r8) :: sumwt(numgrid), a                   ! sum of wt
      real(r8) rhoair,thm,th,thv,ur,displa,zldis,hu,ht,hq
      real(r8) z0m,z0h,z0q,us,vs,tm,qm,pbot,psrf
      real(r8) obu,temp1,temp2,temp12m,temp22m
      real(r8) um,thvstar,beta,zii,wc,wc2

![1-33] Grid averages by area-weight over grid patches
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%taux   ,fldv%taux   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%tauy   ,fldv%tauy   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%fsena  ,fldv%fsena  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%lfevpa ,fldv%lfevpa )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%fevpa  ,fldv%fevpa  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%fsenl  ,fldv%fsenl  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%fevpl  ,fldv%fevpl  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%etr    ,fldv%etr    )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%fseng  ,fldv%fseng  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%fevpg  ,fldv%fevpg  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%fgrnd  ,fldv%fgrnd  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%sabvsun,fldv%sabvsun)
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%sabvsha,fldv%sabvsha)
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%sabg   ,fldv%sabg   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%olrg   ,fldv%olrg   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%rnet   ,fldv%rnet   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%zerr   ,fldv%zerr   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%assim  ,fldv%assim  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%respc  ,fldv%respc  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%tlsun  ,fldv%tlsun  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%tlsha  ,fldv%tlsha  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%ldew   ,fldv%ldew   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%sigf   ,fldv%sigf   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%green  ,fldv%green  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%lai    ,fldv%lai    )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%sai    ,fldv%sai    )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%avsdr  ,fldv%avsdr  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%avsdf  ,fldv%avsdf  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%anidr  ,fldv%anidr  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%anidf  ,fldv%anidf  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%sols   ,fldv%sols   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%soll   ,fldv%soll   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%solsd  ,fldv%solsd  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%solld  ,fldv%solld  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%solrs  ,fldv%solrs  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%solrl  ,fldv%solrl  )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%solrsd ,fldv%solrsd )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%solrld ,fldv%solrld )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%emis   ,fldv%emis   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%z0ma   ,fldv%z0ma   )
      call flux_p2g(numpatch,numgrid,pgmap,wxy_patch,fldv_pft%qsubl  ,fldv%qsubl  )

![84-92] Meteorological forcing
      sumwt(:) = 0.0
      do p = 1, numpatch
         g = pgmap(p)
         sumwt(g) = sumwt(g) + wxy_patch(p)
      enddo

      c = 1
      do g = 1, numgrid
         do while(g.ne.cgmap(c))
            c = c+1
         enddo

         if(sumwt(g).gt.0.)then      !For land only defined
            fldv%us   (g) = forc(3,c)
            fldv%vs   (g) = forc(4,c)
            fldv%tm   (g) = forc(5,c)
            fldv%qm   (g) = forc(6,c)
            fldv%prc  (g) = forc(7,c)
            fldv%prl  (g) = forc(8,c)
            fldv%pbot (g) = forc(9,c)
            fldv%frl  (g) = forc(15,c)
            fldv%solar(g) = forc(11,c)+forc(12,c)+forc(13,c)+forc(14,c)
         else
            fldv%us   (g) = -9999.
            fldv%vs   (g) = -9999.
            fldv%tm   (g) = -9999.
            fldv%qm   (g) = -9999.
            fldv%prc  (g) = -9999.
            fldv%prl  (g) = -9999.
            fldv%pbot (g) = -9999.
            fldv%frl  (g) = -9999.
            fldv%solar(g) = -9999.
            write(6,*) 'impossible grid 4', sumwt(g)
            call abort
         endif
      enddo

      c = 1               ! column index
      do g = 1, numgrid   ! grid index
         do while(g.ne.cgmap(c))
            c = c+1
         enddo

         if(sumwt(g).gt.0.)then         !For land only defined
            z0m = fldv%z0ma(g)
            z0h = fldv%z0ma(g)
            z0q = fldv%z0ma(g)
            displa = 2./3.*z0m/0.07
       
            hu = max(forc(16,c),5.+displa)
            ht = max(forc(17,c),5.+displa)
            hq = max(forc(18,c),5.+displa)
            zldis = hu-displa
       
            us   = fldv%us  (g)
            vs   = fldv%vs  (g)
            tm   = fldv%tm  (g)
            qm   = fldv%qm  (g)
            pbot = fldv%pbot(g)

            psrf = forc  (10,c)
     
            rhoair = (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm)

            fldv%trad (g) = (fldv%olrg(g)/stefnc)**0.25 
            fldv%ustar(g) = sqrt(max(1.e-6,sqrt(fldv%taux(g)**2+fldv%tauy(g)**2))/rhoair)
            fldv%tstar(g) = -fldv%fsena(g)/(rhoair*fldv%ustar(g))/cpair
            fldv%qstar(g) = -fldv%fevpa(g)/(rhoair*fldv%ustar(g))
 
            thm = tm + 0.0098*ht
            th  = tm*(100000./psrf)**(rgas/cpair)
            thv = th*(1.+0.61*qm)       
 
            fldv%zol(g) = zldis*vonkar*grav&
                * (fldv%tstar(g)+0.61*th*fldv%qstar(g))&
                / (fldv%ustar(g)**2*thv)
 
            if(fldv%zol(g) .ge. 0.)then   !stable
               fldv%zol(g) = min(2.,max(fldv%zol(g),1.e-6))
            else                          !unstable
               fldv%zol(g) = max(-100.,min(fldv%zol(g),-1.e-6))
            endif

            beta = 1.
            zii = 1000.
            thvstar=fldv%tstar(g)+0.61*th*fldv%qstar(g)
            ur = sqrt(us*us+vs*vs)
            if(fldv%zol(g) .ge. 0.)then
               um = max(ur,0.1)
            else
               wc = (-grav*fldv%ustar(g)*thvstar*zii/thv)**(1./3.)
              wc2 = beta*beta*(wc*wc)
               um = max(0.1,sqrt(ur*ur+wc2))
            endif

            obu = zldis/fldv%zol(g)
            call moninobuk(hu,ht,hq,displa,z0m,z0h,z0q,&
                 obu,um,fldv%ustar(g),temp1,temp2,temp12m,temp22m,&
                 fldv%f10m(g),fldv%fm(g),fldv%fh(g),fldv%fq(g))
 
            fldv%rib(g) = fldv%zol(g)*vonkar**3*fldv%ustar(g)**2/(temp1*um**2)
            fldv%rib(g) = min(5.,fldv%rib(g)) 
         else
            fldv%trad (g) = -9999.
            fldv%ustar(g) = -9999.
            fldv%tstar(g) = -9999.
            fldv%qstar(g) = -9999.
            fldv%zol  (g) = -9999.
            fldv%rib  (g) = -9999.
            fldv%fm   (g) = -9999.
            fldv%fh   (g) = -9999.
            fldv%fq   (g) = -9999.
            write(6,*) 'impossible grid 2', sumwt(g)
            call abort
         endif
    
      enddo

![42-46] Modified by zhq 06/12/2009/ for matching the routine meteorological obs. If there is grass pft in a grid, fluxes are output as that calculated on grass pft; else fluxes are output as average on a grid.
#ifdef DGVM
      ivt = nint(fvar_pft(104,:))
#else
      ivt = nint(fcon_pft(1,:))
#endif

      call flux_p2g_grass(numpatch,numgrid,pgmap,ivt,wxy_patch,fldv_pft%tref,fldv%tref)
      call flux_p2g_grass(numpatch,numgrid,pgmap,ivt,wxy_patch,fldv_pft%qref,fldv%qref)
      call flux_p2g_grass(numpatch,numgrid,pgmap,ivt,wxy_patch,fldv_pft%u10m,fldv%u10m)
      call flux_p2g_grass(numpatch,numgrid,pgmap,ivt,wxy_patch,fldv_pft%v10m,fldv%v10m)
      call flux_p2g_grass(numpatch,numgrid,pgmap,ivt,wxy_patch,fldv_pft%f10m,fldv%f10m)

      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%xerr  ,fldv%xerr  )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%rsur  ,fldv%rsur  )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%rnof  ,fldv%rnof  )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%tg    ,fldv%tg    )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%scv   ,fldv%scv   )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%snowdp,fldv%snowdp)
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%fsno  ,fldv%fsno  )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%mrsos ,fldv%mrsos )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%mrso  ,fldv%mrso  )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%mrfso ,fldv%mrfso )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%lwsnl ,fldv%lwsnl )
      call flux_c2g(numcolumn,numgrid,cgmap,wxy_column,fldv_col%snm   ,fldv%snm   )

      call flux_c2g_3d(numcolumn,numgrid,nl_soil,cgmap,itypwat,wxy_column,fldv_col%tss  ,fldv%tss  )
      call flux_c2g_3d(numcolumn,numgrid,nl_soil,cgmap,itypwat,wxy_column,fldv_col%wliq ,fldv%wliq )
      call flux_c2g_3d(numcolumn,numgrid,nl_soil,cgmap,itypwat,wxy_column,fldv_col%wice ,fldv%wice )
      call flux_c2g_3d(numcolumn,numgrid,nl_soil,cgmap,itypwat,wxy_column,fldv_col%mrlsl,fldv%mrlsl)

    ! Special treatment for fmicr
      fldv%fmicr(:) = 0.

      do p = 1, numpatch
         c = pcmap(p)
         g = pgmap(p)
         if(itypwat(c).eq.0) then
            fldv%fmicr(g) = fldv%fmicr(g) + fldv_pft%fmicr(p)*wxy_column(c)
         end if
      end do

    ! Special treatment for tsn & nsnow
      sumwt(:) = 0.
      fldv%tsn(:) = 0.
      fldv%nsnow(:) = 0.

      do c = 1, numcolumn
         g = cgmap(c)
         if(fldv_col%tsn(c).gt.1.0E-6) then
            fldv%tsn(g) = fldv%tsn(g) + wxy_column(c)*fldv_col%scv(c)*fldv_col%tsn(c)
            sumwt(g) = sumwt(g) + wxy_column(c)*fldv_col%scv(c)
         end if
      enddo

      do g = 1, numgrid
         if(sumwt(g).gt.0.)then
            fldv%tsn(g)   = fldv%tsn(g)/sumwt(g)
            fldv%nsnow(g) = 1._r8
         else
            fldv%tsn(g)   = 0._r8
            fldv%nsnow(g) = 0._r8
         endif
      enddo

    ! Special treatment for land cover fraction
      fldv%treeFrac(:) = 0.
      fldv%shrubFrac(:) = 0.
      fldv%grassFrac(:) = 0.
      fldv%baresoilFrac(:) = 0.
      fldv%residualFrac(:) = 0.

      do p = 1, numpatch
         g = pgmap(p)

         if(ivt(p).ge.1 .and. ivt(p).le.8) then
            fldv%treeFrac(g) = fldv%treeFrac(g) + wxy_patch(p)
         else if(ivt(p).ge.9 .and. ivt(p).le.11) then
            fldv%shrubFrac(g) = fldv%shrubFrac(g) + wxy_patch(p)
         else if(ivt(p).ge.12 .and. ivt(p).le.14) then
            fldv%grassFrac(g) = fldv%grassFrac(g) + wxy_patch(p)
         else if(ivt(p).eq.17) then
            fldv%baresoilFrac(g) = fldv%baresoilFrac(g) + wxy_patch(p)
         else
            fldv%residualFrac(g) = fldv%residualFrac(g) + wxy_patch(p)
         end if
      end do

      do c = 1, numcolumn
         g = cgmap(c)

         if(itypwat(c).eq.0) then
            fldv%soilFrac(g) = wxy_column(c)
         else if(itypwat(c).eq.1) then
            fldv%urbanFrac(g) = wxy_column(c)
         else if(itypwat(c).eq.2) then
            fldv%wetlandFrac(g) = wxy_column(c)
         else if(itypwat(c).eq.3) then
            fldv%iceFrac(g) = wxy_column(c)
         else if(itypwat(c).eq.4) then
            fldv%lakeFrac(g) = wxy_column(c)
         end if
      end do

      do g = 1, numgrid
         a = fldv%treeFrac(g) + fldv%shrubFrac(g) + fldv%grassFrac(g) &
           + fldv%baresoilFrac(g) + fldv%residualFrac(g)

         if(abs(a-1._r8).gt.1.e-6_r8) then
            write(6,*) "Fraction error: ", fldv%treeFrac(g), fldv%shrubFrac(g),&
                        fldv%grassFrac(g), fldv%baresoilFrac(g), fldv%residualFrac(g), gxmap(g), gymap(g)
            call abort
         end if
      end do

 end subroutine fluxave
