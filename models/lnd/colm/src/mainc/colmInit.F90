#include <define.h>

subroutine colmInit()

   use precision
   use paramodel
   use spmd
   use spmd_decomp
   use colm_varMod
   use colm_cplMod
   use colm_varctl, only: nsrest
   use colm_ioMod, only: lurestart, fsbc, lusbc, readinidata, write_colm_restart
   use timemgr
   implicit none

   integer , parameter   :: numpftx = numpft+1

   real(r8), allocatable :: rockdep(:)
   real(r8), allocatable :: sand(:,:)
   real(r8), allocatable :: clay(:,:)
   real(r8), allocatable :: soc(:,:)
   integer,  allocatable :: isc(:)

   real(r8), allocatable :: avsdr(:)
   real(r8), allocatable :: avsdf(:)
   real(r8), allocatable :: anidr(:)
   real(r8), allocatable :: anidf(:)
   real(r8), allocatable ::  trad(:)
   real(r8), allocatable ::   scv(:)

   real(r8), allocatable :: avsdr_glob(:)
   real(r8), allocatable :: avsdf_glob(:)
   real(r8), allocatable :: anidr_glob(:)
   real(r8), allocatable :: anidf_glob(:)
   real(r8), allocatable ::  trad_glob(:)
   real(r8), allocatable ::   scv_glob(:)

   real(r8), allocatable :: avsdr_xy(:,:)
   real(r8), allocatable :: avsdf_xy(:,:)
   real(r8), allocatable :: anidr_xy(:,:)
   real(r8), allocatable :: anidf_xy(:,:)
   real(r8), allocatable ::  trad_xy(:,:)
   real(r8), allocatable ::   scv_xy(:,:)

   real(r8) zlnd    ! roughness length for soil [m]
   real(r8) zsno    ! roughness length for snow [m]
   real(r8) csoilc  ! drag coefficient for soil under canopy [-]
   real(r8) dewmx   ! maximum dew
   real(r8) wtfact  ! fraction of model area with high water table
   real(r8) capr    ! tuning factor to turn first layer T into surface T
   real(r8) cnfac   ! Crank Nicholson factor between 0 and 1 
   real(r8) ssi     ! irreducible water saturation of snow
   real(r8) wimp    ! water impremeable if porosity less than wimp
   real(r8) pondmx  ! ponding depth (mm)
   real(r8) smpmax  ! wilting point potential in mm
   real(r8) smpmin  ! restriction for min of soil poten. (mm)
   real(r8) trsmx0  ! max transpiration for moist soil+100% veg.  [mm/s]
   real(r8) tcrit   ! critical temp. to determine rain or snow

   real(r8) orb_coszen
   real(r8) calday                      ! Julian cal day (1.xx to 365.xx)
   integer  p1,p2                       ! pft index 
   integer  year                        ! current year of model run
   integer  jday                        ! current julian day of model run
   integer  msec                        ! current seconds of model run (0-86400)
   integer  n_pft                       ! number of pfts in a single column
   integer  i,j,c,p,g                   ! indices of column and pft

#ifdef DUST
   !Always read-in file: no matter whether there is a restart run.
   cvar%dustsource(:)     = dustsource_glob(ccmap)
   cvar%isltyp(:)         = isltyp_glob(ccmap)
   cvar%soil_top_cat(:,:) = soil_top_cat_glob(:,ccmap)
   cvar%mvegcov(:,:)      = mvegcov_glob(:,ccmap)
#endif

   if (lini .or. nsrest.gt.0) then
      call colm_srfvar_free()
      call readinidata()
      return
   end if

   allocate (rockdep     (numcolumn))
   allocate (sand(nl_soil,numcolumn))
   allocate (clay(nl_soil,numcolumn))
   allocate (soc (nl_soil,numcolumn))
   allocate (isc         (numcolumn))

   cvar%dlat       = dlat_glob(ccmap)
   cvar%dlon       = dlon_glob(ccmap)
   cvar%itypwat    = itypwat_glob(ccmap)
   cvar%wxy_column = wxy_column_glob(ccmap)
   pvar%wxy_patch  = wxy_patch_glob(ppmap)
   pvar%ivt        = ivt_glob(ppmap)

   rockdep(:)      = rockdep_glob(ccmap)
   sand(:,:)       = sand_glob(:,ccmap)
   clay(:,:)       = clay_glob(:,ccmap)
   soc(:,:)        = soc_glob(:,ccmap)
   isc(:)          = isc_glob(ccmap)

   call colm_srfvar_free()

! --------------------------------------------------------------------
! [3] 
! INITIALIZE TIME INVARIANT VARIABLES
! ----------------------------------------------------------------------
   write(6,*) 'initialize time invariant variables'

   p1 = 1
   p2 = 0
   do c = 1, numcolumn
#if(defined SINGLE)
      if(cvar%itypwat(c)==0)  n_pft = 1       ! natural vegetation + bare soil
#else
      if(cvar%itypwat(c)==0)  n_pft = numpftx ! natural vegetation + bare soil
#endif
      if(cvar%itypwat(c)==1)  n_pft = 1       ! urban and built-up
      if(cvar%itypwat(c)==2)  n_pft = 1       ! wetland
      if(cvar%itypwat(c)==3)  n_pft = 1       ! land ice
      if(cvar%itypwat(c)==4)  n_pft = 1       ! river or deep lake
      if(cvar%itypwat(c)==99)then             ! ocean
         write(6,*) 'ocean column',c
         stop 
      endif

      p2 = p2 + n_pft      

      CALL iniTimeConst(nl_soil,n_pft,isc(c),pvar%ivt(p1:p2),sand(1:,c),clay(1:,c),soc(1:,c)&
           ,rockdep(c),cvar%itypwat(c),cvar%z(1:,c),cvar%dz(1:,c)&
           ,cvar%albsat(1:,c),cvar%albdry(1:,c),cvar%csol(1:,c)&
           ,cvar%porsl(1:,c),cvar%phi0(1:,c),cvar%bsw(1:,c)&
           ,cvar%dkmg(1:,c),cvar%dksatu(1:,c),cvar%dkdry(1:,c),cvar%hksati(1:,c),cvar%lakedepth(c),cvar%dz_lake(1:,c)&
           ,pvar%z0m(p1:p2),pvar%displa(p1:p2),pvar%sqrtdi(p1:p2),pvar%effcon(p1:p2)&
           ,pvar%vmax25(p1:p2),pvar%slti(p1:p2),pvar%hlti(p1:p2)&
           ,pvar%shti(p1:p2),pvar%hhti(p1:p2),pvar%trda(p1:p2),pvar%trdm(p1:p2)&
           ,pvar%trop(p1:p2),pvar%gradm(p1:p2),pvar%binter(p1:p2)&
           ,pvar%extkn(p1:p2),pvar%chil(p1:p2),pvar%refl(1:,1:,p1:p2),pvar%refs(1:,1:,p1:p2)&
           ,pvar%tranl(1:,1:,p1:p2),pvar%trans(1:,1:,p1:p2),pvar%rootfr(1:,p1:p2)&
           ,zlnd,zsno,csoilc,dewmx,wtfact,capr,cnfac&
           ,ssi,wimp,pondmx,smpmax,smpmin,trsmx0,tcrit)

      p1 = p2 + 1
   enddo
!-----------------------------------------------------------------------
!IF DEFINED DGVM, INITIALIZE PFT's PARAMETER
!-----------------------------------------------------------------------

#ifdef DGVM
   CALL iniDGVMVar()
#endif

! ----------------------------------------------------------------------
! [4]
! INITIALIZE TIME-VARYING VARIABLES, as subgrid vectors of length [numpatch]
! initial run: create the time-varying variables based on :
!              i) observation (NOT CODING CURRENTLY), or
!             ii) some already-known information (NO CODING CURRENTLY), or
!            iii) arbitrarily 
! continuation run: time-varying data read in from restart file 
! ----------------------------------------------------------------------

!4.1 current time of model run
   year = startup_date(1)
   jday = startup_date(2)
   msec = startup_date(3)

   write(6,*) 'local time', year, jday, msec

   if(.not. greenwich)then
   ! convert local time to GMT
   ! off-line cases: input local time was assumed the time at the first grid point
   !                 if not this case, please make your change
   msec = msec - int(longxy(1,1)/15.*3600.)
   if(msec<0)then
      jday = jday -1
      msec = 86400 + msec
   endif
   if(jday<1)then
      year=year-1
      if((mod(year,4)==0 .AND. mod(year,100)/=0) .OR. mod(year,400)==0)then
          jday = 366
       else
          jday = 365
       endif
   endif
   endif

   write(6,*) 'Greenwich time', year, jday, msec

!4.2 cosine of solar zenith angle 
   calday = float(jday)+float(msec)/86400.
   do c = 1, numcolumn
      cvar%coszen(c) = orb_coszen(calday,cvar%dlon(c),cvar%dlat(c))
   enddo

#ifdef SOILINI
!4.3 READ in or GUSSES land state information
   read(lusoil) nl_soil_ini

   allocate (snow_d_grid(lon_points,lat_points))
   allocate (snow_d(numpatch))

   allocate (soil_z_grid(nl_soil_ini))
   allocate (soil_t_grid(nl_soil_ini,lon_points,lat_points))
   allocate (soil_w_grid(nl_soil_ini,lon_points,lat_points))

   allocate (soil_z(nl_soil_ini))
   allocate (soil_t(nl_soil_ini,numpatch))
   allocate (soil_w(nl_soil_ini,numpatch))

   read(lusoil) soil_z_grid  ! soil layer node depth (m)
   read(lusoil) soil_t_grid  ! soil layer temperature (K)
   read(lusoil) soil_w_grid  ! soil layer wetness (-)
   read(lusoil) snow_d_grid  ! snow depth (m)

   close (lusoil)
   soil_z(:) = soil_z_grid(:)

   do npatch = 1, numpatch
      i = ixy_patch(npatch)
      j = jxy_patch(npatch)

      snow_d(npatch) = snow_d_grid(i,j)

      do l = 1, nl_soil_ini
         soil_t(l,npatch) = soil_t_grid(l,i,j)
         soil_w(l,npatch) = soil_w_grid(l,i,j)
      enddo
   end do

   deallocate (snow_d_grid)
   deallocate (snow_d)

   deallocate (soil_z_grid)
   deallocate (soil_t_grid)
   deallocate (soil_w_grid)

   deallocate (soil_z)
   deallocate (soil_t)
   deallocate (soil_w)
#endif

!4.4 initialize time-varying variables, as subgrid vectors of length [numpatch]
   write(6,*) 'initialize time-varying variables'

   p1 = 1
   p2 = 0
   do c = 1, numcolumn
#ifdef SINGLE
      if(cvar%itypwat(c)==0) n_pft = 1       ! natural vegetation + bare soil
#else
      if(cvar%itypwat(c)==0) n_pft = numpftx ! natural vegetation + bare soil
#endif
      if(cvar%itypwat(c)==1) n_pft = 1       ! urban and built-up
      if(cvar%itypwat(c)==2) n_pft = 1       ! wetland
      if(cvar%itypwat(c)==3) n_pft = 1       ! land ice
      if(cvar%itypwat(c)==4) n_pft = 1       ! river or deep lake
      if(cvar%itypwat(c)==99)then                 ! ocean
         write(6,*) 'ocean column',i
         stop
      endif

      p2 = p2 + n_pft    

      CALL iniTimeVar(nl_soil,maxsnl,n_pft,cvar%itypwat(c),pvar%ivt(p1:p2)&
           ,cvar%porsl(1:,c),cvar%albsat(1:,c),cvar%albdry(1:,c),pvar%z0m(p1:p2),pvar%chil(p1:p2)&
           ,pvar%refl(1:,1:,p1:p2),pvar%refs(1:,1:,p1:p2),pvar%tranl(1:,1:,p1:p2),pvar%trans(1:,1:,p1:p2)&
           ,pvar%rootfr(1:,p1:p2),cvar%z(maxsnl+1:,c),cvar%dz(maxsnl+1:,c)&
           ,cvar%tss(maxsnl+1:,c),cvar%wliq(maxsnl+1:,c),cvar%wice(maxsnl+1:,c)&
           ,cvar%tg(c),pvar%tlsun(p1:p2),pvar%tlsha(p1:p2),pvar%ldew(p1:p2),cvar%sag(c),cvar%scv(c)&
           ,cvar%snowdp(c),cvar%t_lake(1:,c),cvar%lake_icefrac(1:,c),cvar%savedtke1(c),pvar%fveg(p1:p2),cvar%fsno(c),pvar%sigf(p1:p2),pvar%green(p1:p2),pvar%lai(p1:p2)&
           ,pvar%sai(p1:p2),cvar%coszen(c),zlnd&
           ,pvar%albg(1:,1:,p1:p2),pvar%albv(1:,1:,p1:p2),pvar%alb(1:,1:,p1:p2),pvar%ssun(1:,1:,p1:p2),pvar%ssha(1:,1:,p1:p2)&
           ,pvar%thermk(p1:p2),pvar%extkb(p1:p2),pvar%extkd(p1:p2)&
           ,pvar%trad(p1:p2),pvar%tref(p1:p2),pvar%qref(p1:p2),pvar%rst(p1:p2),pvar%emis(p1:p2),pvar%z0ma(p1:p2)&
           ,pvar%zol(p1:p2),pvar%rib(p1:p2)&
           ,pvar%ustar(p1:p2),pvar%qstar(p1:p2),pvar%tstar(p1:p2),pvar%fm(p1:p2),pvar%fh(p1:p2),pvar%fq(p1:p2)&
#ifdef SOILINI
           ,nl_soil_ini,soil_z,pvar%soil_t(1:,p1:p2),pvar%soil_w(1:,p1:p2),pvar%snow_d(p1:p2)&
#endif
#ifdef DGVM
           ,cvar%wxy_column(c),pvar%wxy_patch(p1:p2)&
           ,pvar%t10min(p1:p2),pvar%lai_ind(p1:p2)&
           ,pvar%dphen(p1:p2),pvar%leafon(p1:p2),pvar%leafof(p1:p2),pvar%firelength(p1:p2)&
#ifdef IAPDGVM
           ,pvar%afirefrac1(p1:p2),pvar%nfireg1(p1:p2),cvar%wliq6mon(c)&
#endif  
           ,pvar%litterag(p1:p2)  ,pvar%litterbg(p1:p2),pvar%cpool_fast(:,p1:p2),pvar%cpool_slow(:,p1:p2),pvar%k_fast_ave(p1:p2)&
           ,pvar%k_slow_ave(p1:p2),pvar%litter_decom_ave(p1:p2),pvar%fmicr(p1:p2)&
           ,pvar%ifpre(p1:p2),cvar%prec365(c),pvar%nind(p1:p2),pvar%lm_ind(p1:p2),pvar%sm_ind(p1:p2),pvar%hm_ind(p1:p2) &
           ,pvar%rm_ind(p1:p2),pvar%tmomin20(p1:p2),pvar%agdd0(p1:p2),pvar%agdd(p1:p2),pvar%agdd20(p1:p2)&
           ,pvar%t_mo(p1:p2),pvar%t_mo_sum(p1:p2)&
           ,pvar%t_mo_min(p1:p2),pvar%crownarea(p1:p2),pvar%htop(p1:p2),pvar%tsai(p1:p2),pvar%fpcgrid(p1:p2)&
           ,pvar%bm_inc(p1:p2),pvar%afmicr(p1:p2),pvar%annpsn(p1:p2),pvar%annpsnpot(p1:p2),pvar%tref10(p1:p2),pvar%tref_sum(p1:p2)&
           ,pvar%t10(1:,p1:p2),pvar%assimn10(p1:p2),pvar%assimn_sum(p1:p2),pvar%an10(1:,p1:p2),cvar%nday(c),cvar%nyr(c)&
           ,pvar%turnover_ind(p1:p2),pvar%fpc_inc(p1:p2),pvar%agddtw(p1:p2)&
           ,pvar%anngpp(p1:p2),pvar%annfrmf(p1:p2),pvar%annfrms(p1:p2),pvar%annfrmr(p1:p2),pvar%annfrg(p1:p2)&
#endif
#ifdef DyN
           ,pvar%afcton_leaf(p1:p2),pvar%afcton_root(p1:p2),pvar%afcton_sap(p1:p2) &
           ,pvar%litter_leaf(p1:p2),pvar%litter_root(p1:p2),pvar%litter_wood(p1:p2),pvar%litter_repr(p1:p2) &
           ,pvar%litter_leaf_n(p1:p2),pvar%litter_root_n(p1:p2),pvar%litter_wood_n(p1:p2),pvar%litter_repr_n(p1:p2) &
           ,pvar%lm_ind_n(p1:p2),pvar%sm_ind_n(p1:p2),pvar%hm_ind_n(p1:p2),pvar%rm_ind_n(p1:p2),pvar%an_up(p1:p2),pvar%an_stress(p1:p2) &
           ,cvar%soil_no3(c),cvar%soil_no2(c),cvar%soil_no(c),cvar%soil_n2o(c),cvar%soil_n2(c),cvar%soil_nh4(c) &
#endif
#if (defined FHNP) && (defined FTF)
!liruichao add
           ,cvar%frostdp(c),cvar%thawdp(c),cvar%frostdp0(c),cvar%D_temperature(c),cvar%N_time(c),cvar%frost_day(c),cvar%thaw_day(c)&
!end
#endif
           )

      p1 = p2 + 1

   enddo

 ! CLM time step and TUNABLE constants
   ftune(1)  = zlnd   
   ftune(2)  = zsno   
   ftune(3)  = csoilc 
   ftune(4)  = dewmx  
   ftune(5)  = wtfact 
   ftune(6)  = capr   
   ftune(7)  = cnfac  
   ftune(8)  = ssi    
   ftune(9)  = wimp   
   ftune(10) = pondmx 
   ftune(11) = smpmax 
   ftune(12) = smpmin 
   ftune(13) = trsmx0 
   ftune(14) = tcrit  

 ! the model variables for restart run 
   idate(1) = year
   idate(2) = jday
   idate(3) = msec
   idate_p  = idate   ! for COUP_CSM
   istep    = 0

 ! --------------------------------------------
 ! write out as a restart file
 ! --------------------------------------------
 ! CALL rstWrite()

   OPEN(unit=lurestart,file='restartfile',status='unknown',form='unformatted',action='write')

   CALL write_colm_restart()

   CLOSE(lurestart)

   write (6,*) ('successfully to initialize the land model')

 ! --------------------------------------------
 ! write out as a surface boundary condition file
 ! --------------------------------------------
 ! CALL sbcWrite()

 ! average subgrid albedos, srf temperature, etc. for atmospheric model

   allocate (avsdr(numgrid))
   allocate (avsdf(numgrid))
   allocate (anidr(numgrid))
   allocate (anidf(numgrid))
   allocate ( trad(numgrid))
   allocate (  scv(numgrid))

   allocate (avsdr_glob(numgrid_glob))
   allocate (avsdf_glob(numgrid_glob))
   allocate (anidr_glob(numgrid_glob))
   allocate (anidf_glob(numgrid_glob))
   allocate ( trad_glob(numgrid_glob))
   allocate (  scv_glob(numgrid_glob))

   allocate (avsdr_xy(lon_points,lat_points))
   allocate (avsdf_xy(lon_points,lat_points))
   allocate (anidr_xy(lon_points,lat_points))
   allocate (anidf_xy(lon_points,lat_points))
   allocate ( trad_xy(lon_points,lat_points))
   allocate (  scv_xy(lon_points,lat_points))

   call flux_p2g(numpatch ,numgrid,pgmap,pvar%wxy_patch ,pvar%alb(1,1,:),avsdr)
   call flux_p2g(numpatch ,numgrid,pgmap,pvar%wxy_patch ,pvar%alb(1,2,:),avsdf)
   call flux_p2g(numpatch ,numgrid,pgmap,pvar%wxy_patch ,pvar%alb(2,1,:),anidr)
   call flux_p2g(numpatch ,numgrid,pgmap,pvar%wxy_patch ,pvar%alb(2,2,:),anidf)
   call flux_p2g(numpatch ,numgrid,pgmap,pvar%wxy_patch ,pvar%trad      ,trad )
   call flux_c2g(numcolumn,numgrid,cgmap,cvar%wxy_column,cvar%scv       ,scv  )

   scv = scv/1000._r8

#ifdef COUP_CSM
   lnd_avsdr = avsdr
   lnd_avsdf = avsdf
   lnd_anidr = anidr
   lnd_anidf = anidf
   lnd_trad  = trad
   lnd_scv   = scv
#endif

#ifdef SPMD
   p_counts(:) = numgrid_proc(:)
   p_displs(0) = 0
   do i = 1, p_nprocs-1
      p_displs(i) = sum(numgrid_proc(0:i-1))
   end do

   call mpi_gatherv (avsdr,size(avsdr),mpi_real8,avsdr_glob,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
   call mpi_gatherv (avsdf,size(avsdf),mpi_real8,avsdf_glob,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
   call mpi_gatherv (anidr,size(anidr),mpi_real8,anidr_glob,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
   call mpi_gatherv (anidf,size(anidf),mpi_real8,anidf_glob,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
   call mpi_gatherv (trad ,size(trad ),mpi_real8, trad_glob,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
   call mpi_gatherv (scv  ,size(scv  ),mpi_real8,  scv_glob,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
   avsdr_glob = avsdr
   avsdf_glob = avsdf
   anidr_glob = anidr
   anidf_glob = anidf
    trad_glob = trad
     scv_glob = scv
#endif

   if (p_master) then
      open(lusbc,file="sbcfile",form="unformatted",status="unknown")

      avsdr_xy(:,:) = -9999.
      avsdf_xy(:,:) = -9999.
      anidr_xy(:,:) = -9999.
      anidf_xy(:,:) = -9999.
       trad_xy(:,:) = -9999.
        scv_xy(:,:) = -9999.

      do g = 1, numgrid_glob
         i = gxmap_glob(g)
         j = gymap_glob(g)
         avsdr_xy(i,j) = avsdr_glob(g)
         avsdf_xy(i,j) = avsdf_glob(g)
         anidr_xy(i,j) = anidr_glob(g)
         anidf_xy(i,j) = anidf_glob(g)
          trad_xy(i,j) =  trad_glob(g)
           scv_xy(i,j) =   scv_glob(g)
      end do

      write(lusbc) avsdr_xy
      write(lusbc) avsdf_xy
      write(lusbc) anidr_xy
      write(lusbc) anidf_xy
      write(lusbc)  trad_xy
      write(lusbc)   scv_xy

      close(lusbc)
   end if

   deallocate (avsdr)
   deallocate (avsdf)
   deallocate (anidr)
   deallocate (anidf)
   deallocate ( trad)
   deallocate (  scv)

   deallocate (avsdr_glob)
   deallocate (avsdf_glob)
   deallocate (anidr_glob)
   deallocate (anidf_glob)
   deallocate ( trad_glob)
   deallocate (  scv_glob)

   deallocate (avsdr_xy)
   deallocate (avsdf_xy)
   deallocate (anidr_xy)
   deallocate (anidf_xy)
   deallocate ( trad_xy)
   deallocate (  scv_xy)

   deallocate (rockdep)
   deallocate (sand   )
   deallocate (clay   )
   deallocate (soc    )
   deallocate (isc    )

end subroutine colmInit
