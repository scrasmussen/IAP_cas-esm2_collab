#include <define.h>

module colm_ioMod

 use precision
 use spmd
 use spmd_decomp
 use nchistMod
 use colm_varctl
 use colm_varMod
 use timemgr

 implicit none

 character(len=255) :: NLFilename = ""

 character(LEN=255) :: fsrf      = ""       ! file name of raw surface data
 character(LEN=255) :: flai      = ""       ! file name of time-varying vegetation data
 character(LEN=255) :: fmet      = ""       ! file name of meteorological data
 character(LEN=255) :: fco2      = ""       ! file name of global mean annual CO2 concentration data
 character(len=255) :: fini      = ""       ! file name of initial land surface condition 
 character(len=255) :: fsbc      = ""       ! file name of land surface boundary condition 
 character(len=255) :: frestart  = ""       ! file name of restart time-varying file
 character(len=255) :: foutdir   = "."      ! directory name of output files
#ifdef IAPDGVM
 character(len=255) :: fnig                 ! file name of ig for IAPDGVM
#endif

#ifdef DUST
 character(len=255) :: fdust     = ""       ! file name of land surface data for dust module
#endif

 integer            :: lusrf     = 100      ! logical unit of raw surface data
 integer            :: lulai     = 110      ! logical unit of vegetation data forcing
 integer            :: lumet     = 120      ! logical unit of meteorological forcing
 integer            :: lusbc     = 130      ! logical unit of surface boundary condition
 integer            :: lurestart = 140      ! logical unit of restart time-varying file

 integer            :: lugrid    = 180      ! logical unit of surface grid data
 character(len=255) :: fgrid     = ""       ! file name of surface grids info
#ifdef DUST
 integer            :: ludust    = 169      ! logical unit of surface data for dust module
#endif

 interface readnml
    module procedure readnml
 end interface

 interface readsrfdata
    module procedure readsrfdata
 end interface

 interface readgridata
    module procedure readgridata
 end interface

 interface readinidata
    module procedure readinidata
 end interface

 interface readc
    module procedure readc_1d_i4
    module procedure readc_1d_r8
    module procedure readc_2d_r8
    module procedure readc_3d_r8
 end interface

 interface readp
    module procedure readp_1d_i4
    module procedure readp_1d_r8
    module procedure readp_2d_r8
    module procedure readp_3d_r8
 end interface

 interface writec
    module procedure writec_1d_i4
    module procedure writec_1d_r8
    module procedure writec_2d_r8
    module procedure writec_3d_r8
 end interface

 interface writep
    module procedure writep_1d_i4
    module procedure writep_1d_r8
    module procedure writep_2d_r8
    module procedure writep_3d_r8
 end interface

#ifdef IAPDGVM
interface readigdata
      module procedure readigdata
   end interface

   interface getig
      module procedure getig
   end interface

   interface getrhm
      module procedure getrhm
   end interface
#endif

 interface writehistdata
    module procedure writehistdata
 end interface

CONTAINS

   subroutine readnml()

      use spmd
      use timemgr
      use colm_varMod
#ifdef RTM
      use colm_rtmVar
      use colm_rtmMod
#endif

#if (defined FHNP) && (defined NP)
      use pollution, only: fnpoint_rtm
#endif

#ifdef COUP_CSM
      use colm_cplMod, only: irad, nsrest
#ifdef CPL6
      use colm_csmMod, only: csm_dtime
#endif
#endif
      use landuse, only: fluc_emission
      use shr_sys_mod

      implicit none

! local variables:

      integer i, j

      namelist /clmexp/ caseid,           &
                        fsrf,             &
                        flai,             &
                        fmet,             &
                        fco2,             &
                        fini,             &
                        foutdir,          &
#ifdef IAPDGVM
                        fnig,             &
#endif

#ifdef RTM
                        frivinp_rtm,      &
                        rtm_nsteps,       &
#if (defined FHNP) && (defined NP)
                        fnpoint_rtm,      &     ! wangy
#endif

#endif
                        fluc_emission,    &
#ifdef DUST
                        fdust,            &
#endif
#ifdef COUP_CSM
                        fsbc,             &
                        irad,             &
                        csm_doflxave,     &
                        lnd_cflux_year,   &
#ifdef CPL6
                        nsrest,           &
                        co2_option,       &
#endif
#ifdef CPL7
                        co2_type,         &
                        co2_ppmv,         &
#endif
#endif
                        lhist_yearly,     &
                        lhist_monthly,    &
                        lhist_daily,      &
                        lhist_6hourly,    &
                        lhist_3hourly,    &
                        lhist_1hourly,    &
                        startup_date,     &
                        greenwich,        &
                        restart_freq,     &
                        lon_points,       &
                        lat_points,       &
                        dtime,            &
                        mstep

! routine:

      istep = 0    ! COUP_CSM : istep = 0; OFFLINE : istep will be advanced to istep=1.

#ifdef RTM
      rtm_nsteps = -9999
#endif

#ifdef COUP_CSM
      irad = -1
      csm_doflxave = .true.
#endif

      if (p_master) then
         if (len_trim(NLFilename)>0) then
            open(11,file=trim(NLFilename),form='formatted',status='old')
            read(11,nml=clmexp)
            close(11)
         else
            read(5,nml=clmexp)
         end if

         if (len_trim(fluc_emission).gt.0) luc_emission = .TRUE.
! ===========wangy==========0
#if (defined FHNP) && (defined NP)
         if (len_trim(fnpoint_rtm).gt.0) N_pollution = .TRUE.
#endif
! ===========wangy==========1
         if (len_trim(fsrf).eq.0) then
            write(6,*) 'fatal error: no surface data'
            call abort
         end if

#ifdef DUST
         if (len_trim(fdust).eq.0) then
            write(6,*) 'fatal error: no surface data for dust emission'
            call abort
         end if
#endif

         if (len_trim(fini).gt.0) lini = .true.
      end if

      if (p_master) then 
          write(6,clmexp)
          call shr_sys_flush(6)
          write(6,*) '-------CoLM namelist is ok------'
          call shr_sys_flush(6)
          write(6,*) 'fsrf =', trim(fsrf)
          call shr_sys_flush(6)
          write(6,*) 'flai =', trim(flai)
          call shr_sys_flush(6)
          write(6,*) 'fmet =', trim(fmet)
          call shr_sys_flush(6)
          write(6,*) 'fco2 =', trim(fco2)
          call shr_sys_flush(6)
          write(6,*) 'fini =', trim(fini)
          call shr_sys_flush(6)
#ifdef DUST
          write(6,*) 'fdust =', trim(fdust)
          call shr_sys_flush(6)
#endif
          write(6,*) 'lon =', lon_points
          call shr_sys_flush(6)
          write(6,*) 'lat =', lat_points
          call shr_sys_flush(6)
          write(6,*) 'dtime =', dtime
          call shr_sys_flush(6)
          write(6,*) 'mstep =', mstep
          call shr_sys_flush(6)
      end if 

#ifdef SPMD
      call mpi_bcast (dtime     ,1,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (mstep     ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lon_points,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lat_points,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lini      ,1,mpi_logical,0,p_comm,p_err)
#ifdef RTM
      call mpi_bcast (rtm_nsteps,1,mpi_integer,0,p_comm,p_err)
! ===========wangy==========0
#if (defined FHNP) && (defined NP)
      call mpi_bcast (N_pollution,1,mpi_logical,0,p_comm,p_err)
#endif
! ===========wangy==========1
#endif

#ifdef COUP_CSM
      call mpi_bcast (irad,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (nsrest,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (csm_doflxave,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lnd_cflux_year,1,mpi_integer,0,p_comm,p_err)
#ifdef CPL6
      call mpi_bcast (co2_option,len(co2_option),mpi_character,0,p_comm,p_err)
#endif
#ifdef CPL7
      call mpi_bcast (co2_type,len(co2_type),mpi_character,0,p_comm,p_err)
      call mpi_bcast (co2_ppmv,1,mpi_real8,0,p_comm,p_err)
#endif
#endif

      call mpi_bcast (lhist_yearly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_monthly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_daily,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_6hourly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_3hourly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_1hourly,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (lhist_steply,1,mpi_logical,0,p_comm,p_err)

      call mpi_bcast (greenwich,1,mpi_logical,0,p_comm,p_err)
      call mpi_bcast (startup_date,size(startup_date),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (restart_freq,len(restart_freq),mpi_character,0,p_comm,p_err)
      call mpi_bcast (history_freq,len(history_freq),mpi_character,0,p_comm,p_err)

      call mpi_bcast (luc_emission,1,mpi_logical,0,p_comm,p_err)
#endif

#ifdef COUP_CSM
      call ymdtod_to_yjs(start_ymd,start_tod,startup_date)
#endif
      call timeinit(startup_date,startup_date,istep)

#ifdef COUP_CSM
#ifdef CPL6
      csm_dtime = dtime
#endif

      if (irad < 0) irad = nint(-irad*3600./dtime)

      if (csm_doflxave .and. irad==1) then
         if (p_master) then
            write(6,*)'error: irad must be greater that one if', &
              ' flux averaging option is enabled'
         endif
         call abort
      endif
#endif

#ifdef RTM
      if (p_master) write(6,*) 'RTM river data = ',trim(frivinp_rtm)
#endif

   end subroutine readnml

   subroutine readsrfdata()

      use precision
      use phycon_module
      use paramodel, only: numpft
      use colm_varMod, only: lon_points, lat_points, numpatch_glob, numcolumn_glob
      use spmd, only: p_master
      implicit none

! ------------------------ local variables -----------------------------
! surface classification and soil information

      integer,  parameter   :: numpftx = numpft+1

      integer,  allocatable :: soic2d  (:,:)          ! soil color
      real(r8), allocatable :: rock2d  (:,:)          ! depth to bed rock
      real(r8), allocatable :: terrain (:,:)          ! 
      real(r8), allocatable :: glacier (:,:)          ! 
      real(r8), allocatable :: sand2d  (:,:,:)        ! percentage of sand
      real(r8), allocatable :: clay2d  (:,:,:)        ! percentage of clay
      real(r8), allocatable :: densoc  (:,:,:)        ! density of soil organic carbon [kg/m^3]
      integer,  allocatable :: surf2d  (:,:,:)        ! land cover type
      real(r8), allocatable :: fpatch2d(:,:,:)        ! subgrid weights

      integer,  allocatable :: landmaskx(:,:)

#ifdef DUST
      integer,  allocatable :: dustsource2d  (:,:)      ! index for potential dust source
      integer,  allocatable :: isltyp2d  (:,:)          ! dominant soil type
      real(r8), allocatable :: soil_top_cat2d  (:,:,:)  ! fraction for 12-categories soil
      real(r8), allocatable :: mvegcov2d (:,:,:)        ! read_in monthly vegetation cover [unit: %, i.e. 0-100]
#endif

#ifdef SOILINI
      integer, intent(in)   :: lusoil
      integer nl_soil_ini
      real(r8), allocatable :: snow_d_grid(:,:) 
      real(r8), allocatable :: snow_d(:)

      real(r8), allocatable :: soil_z_grid(:)
      real(r8), allocatable :: soil_t_grid(:,:,:)
      real(r8), allocatable :: soil_w_grid(:,:,:)
      real(r8), allocatable :: soil_z(:)
      real(r8), allocatable :: soil_t(:,:)
      real(r8), allocatable :: soil_w(:,:)
#endif

      integer   numpatch_lat(lat_points)    ! number of patches of grids at lon. strip
      integer   numcolumn_lat(lat_points)   ! number of columns of grids at lon. strip

      real(r8)  pi                          ! pie
      real(r8)  fcolumn_tmp                 ! weight of column relative to grid 
      integer   npatch,ncolumn,ncolumn0     ! indices
      integer   i,j,l,m,n,jm,p              ! indices			

#ifdef MYBUG
      real(r8), allocatable :: frac(:,:,:)
#endif

! ----------------------------------------------------------------------
! [1] READ IN LAND INFORMATION
! read time-invariant boundary data on [lon_points] x [lat_points] grid.

!    o first [lat_points] values: number of longitude points for each latitude.
!      this allows for variable longitudinal resolution for each latitude

! remaining data is for each grid cell:

!    o 1st : latitude at center of grid cell (degrees)
!    o 2th : longitude at center of grid cell (degrees)
!    o 3th : soil color (1 to 8) for use with soil albedos
!    o 4th : depth to bed rock
!    o 5th : soil texture, %sand, for thermal and hydraulic properties
!    o 6th : soil texture, %clay, for thermal and hydraulic properties
!    o 7th : surface type, for use as multiple subgrid point
!    o 8th : subgrid weight
!
! ------------------ note ------------------
! the model required the same longitudinal resolution for 
! all latitude strip. For the variable longitudinal resolution
! cases, please assign the surface types and soil character:
! 0 to ocean grids and -999 to the land grids 
! which are not included in the calcultion.
! ----------------------------------------------------------------------

      call colm_gridvar_alloc()

      allocate (soic2d  (lon_points,lat_points))
      allocate (rock2d  (lon_points,lat_points))
      allocate (terrain (lon_points,lat_points))
      allocate (glacier (lon_points,lat_points))
      allocate (sand2d  (nl_soil,lon_points,lat_points))
      allocate (clay2d  (nl_soil,lon_points,lat_points))
      allocate (densoc  (nl_soil,lon_points,lat_points))
      allocate (surf2d  (maxpatch,lon_points,lat_points))
      allocate (fpatch2d(maxpatch,lon_points,lat_points))

      allocate (landmaskx(lon_points,lat_points))

#ifdef DUST
      allocate (dustsource2d  (lon_points,lat_points))
      allocate (isltyp2d      (lon_points,lat_points))
      allocate (soil_top_cat2d(12,lon_points,lat_points))
      allocate (mvegcov2d     (12,lon_points,lat_points))
#endif

      if (p_master) then
         OPEN(unit=lusrf,file=trim(fsrf),status='old',form='unformatted',action='read')

         read(lusrf) numlon
         read(lusrf) lats
         read(lusrf) lonw
         read(lusrf) area
         read(lusrf) latixy
         read(lusrf) longxy
         read(lusrf) landfrac
         read(lusrf) landmask

         read(lusrf) soic2d
         read(lusrf) rock2d
         read(lusrf) terrain
         read(lusrf) glacier
         read(lusrf) sand2d
         read(lusrf) clay2d
         read(lusrf) densoc
         read(lusrf) surf2d
         read(lusrf) fpatch2d

         CLOSE (lusrf)

#ifdef DUST
         OPEN(unit=ludust,file=trim(fdust),status='old',form='unformatted',action='read')

         read(ludust) dustsource2d
         read(ludust) isltyp2d
         read(ludust) soil_top_cat2d
         read(ludust) mvegcov2d
 
         CLOSE (ludust)
#endif
 
      end if

#ifdef SPMD
      call mpi_bcast (numlon  ,size(numlon)  ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lats    ,size(lats)    ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (lonw    ,size(lonw)    ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (area    ,size(area)    ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (latixy  ,size(latixy)  ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (longxy  ,size(longxy)  ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (landfrac,size(landfrac),mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (landmask,size(landmask),mpi_integer,0,p_comm,p_err)

      call mpi_bcast (soic2d  ,size(soic2d)  ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (rock2d  ,size(rock2d)  ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (terrain ,size(terrain) ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (glacier ,size(glacier) ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (sand2d  ,size(sand2d)  ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (clay2d  ,size(clay2d)  ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (densoc  ,size(densoc)  ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (surf2d  ,size(surf2d)  ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (fpatch2d,size(fpatch2d),mpi_real8,  0,p_comm,p_err)

#ifdef DUST
      call mpi_bcast (dustsource2d  ,size(dustsource2d)  ,mpi_integer  ,0,p_comm,p_err)
      call mpi_bcast (isltyp2d      ,size(isltyp2d)      ,mpi_integer  ,0,p_comm,p_err)
      call mpi_bcast (soil_top_cat2d,size(soil_top_cat2d),mpi_real8    ,0,p_comm,p_err)
      call mpi_bcast (mvegcov2d     ,size(mvegcov2d)     ,mpi_real8    ,0,p_comm,p_err)
#endif

#endif

! ----------------------------------------------------------------------
! [2] MAPPING and ALLOCATE
! Build 1d subgrid patch <-> 2d grid mapping indices and weights
! 
! Build mapping indices and weights: [lon_points]x[lat_points] 2d grid <->
! <-> [numpatch] vector of subgrid patches. 
! The land surface model works by gathering all the land points on a
! [lon_points]x[lat_points] grid into a vector, and then expanded into 
! a vector of [numpatch] subgrid patches, allowing
! for up to [maxpatch] subgrid patches per land point. 
! [ixy], [jxy], [patch], and [land] are indices for the mapping: 
! [lon_points]x[lat_points] grid <-> [numpatch] vector of subgrid points. 
!
!-----------------------------------------------------------------------
! Find total number of patches [numpatch] allowing for multiple subgrid 
! patches in a grid cell.
! --------------------------------------------------------------------
      ncolumn = 0
      npatch  = 0
      numpatch_lat(:) = 0
      numcolumn_lat(:) = 0

      do j = 1, lat_points
         do i = 1, lon_points
            fcolumn_tmp = 0.
            do m = 1, maxpatch
               if(m/=oceancateg)then
#if(defined SINGLE)
                ! one pft on each column
                  if(fpatch2d(m,i,j)>1.0E-6)then
                     ncolumn = ncolumn +1
                     numcolumn_lat(j) = numcolumn_lat(j) +1
                     npatch = npatch +1
                     numpatch_lat(j) = numpatch_lat(j) +1
                  endif
#else
                  n = surf2d(m,i,j)
                  if(n<=numpftx)then
                     fcolumn_tmp = fcolumn_tmp + fpatch2d(m,i,j)
                     if(n==numpftx.and.fcolumn_tmp>1.0E-6)then
                        ncolumn = ncolumn + 1
                        numcolumn_lat(j) = numcolumn_lat(j) + 1
                        npatch  = npatch + numpftx
                        numpatch_lat(j) = numpatch_lat(j) + numpftx
                     endif    ! natural vegetation
                  else
                     if(fpatch2d(m,i,j)>1.0E-6)then
                        ncolumn = ncolumn +1
                        numcolumn_lat(j) = numcolumn_lat(j) +1
                        npatch = npatch +1
                        numpatch_lat(j) = numpatch_lat(j) +1
                     endif
                  endif
#endif
               endif
            enddo
         enddo
      enddo

      numpatch_glob = npatch
      numcolumn_glob = ncolumn

      if(numpatch_glob.ne.sum(numpatch_lat))then
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         call abort
      endif
      write(6,*) 'Total land patches = ', numpatch_glob

      if(numcolumn_glob.ne.sum(numcolumn_lat))then
         write(6,*) 'Total number of columns NOT as the summation of numcolumn_lat'
         call abort
      endif
      write(6,*) 'Total land columns = ', numcolumn_glob

! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch_glob] variables
! --------------------------------------------------------------------

      call colm_srfvar_alloc()

      call colm_subgridvar_alloc()

! --------------------------------------------------------------------
! Build 1d land vector and 1d patch vector mapping components
! Add a column vector BTW. patch and grid by zhq. 05/05/2009
! --------------------------------------------------------------------

! Determine land vector and patch vector mapping components

      wxy_patch_glob(:)  = 0.
      wxy_column_glob(:) = 0.
      npatch  = 0
      ncolumn = 0
      ncolumn0 = 0
      do j = 1, lat_points
         do i = 1, lon_points
            fcolumn_tmp = 0.
            do m = 1, maxpatch
               if(m/=oceancateg) then
#if(defined SINGLE)
                  if(fpatch2d(m,i,j)>1.0E-6) then
                     ncolumn = ncolumn + 1
                     wxy_column_glob(ncolumn) = fpatch2d(m,i,j)

                     npatch = npatch + 1
                     wxy_patch_glob(npatch) = fpatch2d(m,i,j)
                     ixy_patch_glob(npatch) = i               !patch longitude index
                     jxy_patch_glob(npatch) = j               !patch latitude index
                           ivt_glob(npatch) = surf2d(m,i,j)   !land cover type
                  endif
#else
                  n = surf2d(m,i,j)
                  if(n<=numpftx)then
                     fcolumn_tmp = fcolumn_tmp + fpatch2d(m,i,j)

                     if(n==numpftx.and.fcolumn_tmp>1.0E-6)then
                        ncolumn = ncolumn + 1
                        wxy_column_glob(ncolumn) = fcolumn_tmp

                        do jm = 1, numpftx
                           npatch = npatch + 1
                           wxy_patch_glob(npatch) = fpatch2d(jm,i,j)
                           ixy_patch_glob(npatch) = i               !patch longitude index
                           jxy_patch_glob(npatch) = j               !patch latitude index
                                 ivt_glob(npatch) = surf2d(jm,i,j)  !land cover type
                        enddo
                     endif
                  else
                     if(fpatch2d(m,i,j)>1.0E-6) then
                        ncolumn = ncolumn + 1
                        wxy_column_glob(ncolumn) = fpatch2d(m,i,j)

                        npatch = npatch + 1
                        wxy_patch_glob(npatch) = fpatch2d(m,i,j)
                        ixy_patch_glob(npatch) = i               !patch longitude index
                        jxy_patch_glob(npatch) = j               !patch latitude index
                              ivt_glob(npatch) = surf2d(m,i,j)   !land cover type
                     endif
                  endif
#endif
                  if(ncolumn/=ncolumn0) then
                      dlat_glob(ncolumn) = latixy(i,j)     !latitude in radians
                      dlon_glob(ncolumn) = longxy(i,j)     !longitude in radians
                ixy_column_glob(ncolumn) = i               !column longitude index
                jxy_column_glob(ncolumn) = j               !column latitude index
                       isc_glob(ncolumn) = soic2d(i,j)     !soil color index
                   rockdep_glob(ncolumn) = rock2d(i,j)     !depth to bed rock
                    sand_glob(:,ncolumn) = sand2d(:,i,j)   !percent of sand
                    clay_glob(:,ncolumn) = clay2d(:,i,j)   !percent of clay
                     soc_glob(:,ncolumn) = densoc(:,i,j)   !density of organic carbon

#ifdef DUST 
                     dustsource_glob(ncolumn)     = dustsource2d(i,j)     ! index for potential dust source
                     isltyp_glob(ncolumn)         = isltyp2d(i,j)         ! dominant soil type
                     soil_top_cat_glob(:,ncolumn) = soil_top_cat2d(:,i,j) ! fraction for 12-categories soil
                     mvegcov_glob(:,ncolumn)      = mvegcov2d(:,i,j)      ! read_in monthly vegetation cover [unit: %, i.e. 0-100]
#endif

!-----------------------------------------------------------------------
!land water type for PFT classification
!-----------------------------------------------------------------------
                     p = ivt_glob(npatch)
                     if(p<=17) itypwat_glob(ncolumn)=0  ! natural vegetation + soil	
                     if(p==21) itypwat_glob(ncolumn)=1  ! urban and built-up
                     if(p==19) itypwat_glob(ncolumn)=2  ! wetland
                     if(p==20) itypwat_glob(ncolumn)=3  ! land ice
                     if(p==18) itypwat_glob(ncolumn)=4  ! river or deep lake
                     if(p==22) itypwat_glob(ncolumn)=99 ! ocean
                  end if 

                  ncolumn0 = ncolumn

               end if
            end do
         end do
      end do

      if(numpatch_glob.ne.npatch)then
         write(6,*) 'the number of patches is not identical ', numpatch_glob, npatch
         call abort
      endif

      if(numcolumn_glob.ne.ncolumn)then
         write(6,*) 'the number of columnes is not identical ', numcolumn_glob, ncolumn
         call abort
      endif

    ! convert latitudes and longitudes from degrees to radians
      pi = 4.*atan(1.)
      dlat_glob(:) = dlat_glob(:)*pi/180. 
      dlon_glob(:) = dlon_glob(:)*pi/180. 

#ifdef MYBUG
      allocate(frac(lon_points,lat_points,5))
      frac(:,:,:) = -9999.

      do i = 1, numcolumn_glob
         if(itypwat_glob(i)<5) then 
            j = itypwat_glob(i)+1
            frac(ixy_column_glob(i),jxy_column_glob(i),j) = wxy_column_glob(i)
         end if
      end do

      do i = 1, 5 
         write(100), frac(:,:,i)
      end do

      deallocate(frac)
#endif

      print *, 'sum of landmask(checking on grid level)', sum(landmask)

      landmaskx(:,:) = 0
      l = ixy_column_glob(1) + (jxy_column_glob(1)-1)*lon_points
      landmaskx(ixy_column_glob(1),jxy_column_glob(1)) = 1

      do i = 1, numcolumn_glob
         j = ixy_column_glob(i) + (jxy_column_glob(i)-1)*lon_points

         if(i.gt.1 .and. j.ne.l) then
            landmaskx(ixy_column_glob(i),jxy_column_glob(i)) = 1
         end if

         l = j
      end do

      print *, 'sum of landmask(checking on column level)', sum(landmaskx)

      if(any(landmaskx.ne.landmask)) then
         stop 'landmask checking failed on column level'
      end if

      landmaskx(:,:) = 0
      l = ixy_patch_glob(1) + (jxy_patch_glob(1)-1)*lon_points
      landmaskx(ixy_patch_glob(1),jxy_patch_glob(1)) = 1

      do i = 1, numpatch_glob
         j = ixy_patch_glob(i) + (jxy_patch_glob(i)-1)*lon_points

         if(i.gt.1 .and. j.ne.l) then
            landmaskx(ixy_patch_glob(i),jxy_patch_glob(i)) = 1
         end if

         l = j
      end do

      print *, 'sum of landmask(checking on patch level)', sum(landmaskx)

      if(any(landmaskx.ne.landmask)) then
         stop 'landmask checking failed on patch level'
      end if

#ifdef MYBUG
      open(300,file='landmaskx',form='unformatted',status='unknown')
      write(300) landmaskx
      close(300)
#endif

      deallocate (soic2d  )
      deallocate (rock2d  )
      deallocate (terrain )
      deallocate (glacier )
      deallocate (sand2d  )
      deallocate (clay2d  )
      deallocate (densoc  )
      deallocate (surf2d  )
      deallocate (fpatch2d)

#ifdef DUST
      deallocate (dustsource2d  )
      deallocate (isltyp2d      )
      deallocate (soil_top_cat2d)
      deallocate (mvegcov2d     )
#endif

   end subroutine readsrfdata

   subroutine readgridata

      use spmd
      use colm_varMod
!*#ifdef COUP_CSM
!*      use colm_csmMod, only : csm_recvgrid
!*#endif

      implicit none

!*#ifdef COUP_CSM
!*      integer , pointer :: cam_numlon(:)      !cam number of longitudes
!*      real(r8), pointer :: cam_longxy(:,:)    !cam lon values
!*      real(r8), pointer :: cam_latixy(:,:)    !cam lat values
!*      real(r8), pointer :: cam_landfrac(:,:)  !cam fractional land
!*      integer , pointer :: cam_landmask(:,:)  !cam land mask
!*#endif

      integer i, j 
 
! routine:

    ! Read in surface grid info.
      if (p_master) then
         OPEN(lugrid,file=fgrid,form='unformatted',status='old')
         read(lugrid) numlon(:)
         read(lugrid) lats(:)
         read(lugrid) lonw(:,:)
         read(lugrid) area(:,:)
         read(lugrid) latixy(:,:)
         read(lugrid) longxy(:,:)
         read(lugrid) landfrac(:,:)
         read(lugrid) landmask(:,:)

       !*Following should be moved into <fvar> array and restart archieve.
       !*As they could change with DGVM enabled.
         read(lugrid) itypwat_glob(:)
         read(lugrid) ixy_column_glob(:) !longitude index for each patch point
         read(lugrid) jxy_column_glob(:) !latitude index for each patch point
         read(lugrid) wxy_column_glob(:) !subgrid weight for each patch point
         read(lugrid) ixy_patch_glob(:)  !longitude index for each patch point
         read(lugrid) jxy_patch_glob(:)  !latitude index for each patch point
         read(lugrid) wxy_patch_glob(:)  !subgrid weight for each patch point
       !**********************************************************************
         CLOSE(lugrid)
      end if

#ifdef SPMD
      call mpi_bcast (numlon         ,size(numlon)         ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (lats           ,size(lats)           ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (lonw           ,size(lonw)           ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (area           ,size(area)           ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (latixy         ,size(latixy)         ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (longxy         ,size(longxy)         ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (landfrac       ,size(landfrac)       ,mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (landmask       ,size(landmask)       ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (itypwat_glob   ,size(itypwat_glob)   ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (ixy_column_glob,size(ixy_column_glob),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (jxy_column_glob,size(jxy_column_glob),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (wxy_column_glob,size(wxy_column_glob),mpi_real8  ,0,p_comm,p_err)
      call mpi_bcast (ixy_patch_glob ,size(ixy_patch_glob) ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (jxy_patch_glob ,size(jxy_patch_glob) ,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (wxy_patch_glob ,size(wxy_patch_glob) ,mpi_real8  ,0,p_comm,p_err)
#endif

!*#ifdef COUP_CSM
!*      allocate (cam_numlon(lat_points))               !cam number of longitudes
!*      allocate (cam_longxy(lon_points,lat_points))    !cam lon values
!*      allocate (cam_latixy(lon_points,lat_points))    !cam lat values
!*      allocate (cam_landfrac(lon_points,lat_points))  !cam fractional land
!*      allocate (cam_landmask(lon_points,lat_points))  !cam land mask
!*
!*    ! Get grid and land mask back from flux coupler
!*
!*      call csm_recvgrid (cam_longxy, cam_latixy, cam_numlon, cam_landfrac, cam_landmask)
!*
!*    ! Determine land grid, land mask and land fraction
!*
!*    ! For CoLM making surface data by a offline fashion, consistance checking between 
!*    ! land grid info & cpl grid info is required here.
!*
!*      do j = 1, lat_points
!*         if(numlon(j).ne.cam_numlon(j)) then
!*            write(6,*) 'numlon(j).ne.cam_numlon(j), j=', j, numlon(j), cam_numlon(j)
!*            call abort
!*         end if
!*      end do
!*
!*      do j = 1,lat_points
!*         do i = 1,numlon(j)
!*            if(abs(longxy(i,j)-cam_longxy(i,j)).gt.1.0E-8) then
!*               write(6,*) 'longxy(i,j).ne.cam_longxy(i,j), i,j=', i,j, longxy(i,j), cam_longxy(i,j)
!*               call abort
!*            end if
!*            if(abs(latixy(i,j)-cam_latixy(i,j)).gt.1.0E-8) then
!*               write(6,*) 'latixy(i,j).ne.cam_latixy(i,j), i,j=', i,j, latixy(i,j), cam_latixy(i,j)
!*               call abort
!*            end if
!*            if(abs(landfrac(i,j)-cam_landfrac(i,j)).gt.1.0E-8) then
!*               write(6,*) 'landfrac(i,j).ne.cam_landfrac(i,j), i,j=', i,j, landfrac(i,j), cam_landfrac(i,j)
!*             ! call abort
!*            end if
!*            if(abs(landmask(i,j).ne.cam_landmask(i,j))) then
!*               write(6,*) 'landmask(i,j).ne.cam_landmask(i,j), i,j=', i,j, landmask(i,j), cam_landmask(i,j)
!*             ! call abort
!*            end if
!*
!*            longxy(i,j)   = cam_longxy(i,j)
!*            latixy(i,j)   = cam_latixy(i,j)
!*            landmask(i,j) = cam_landmask(i,j)
!*            landfrac(i,j) = cam_landfrac(i,j)
!*         end do
!*      end do
!*
!*      deallocate (cam_numlon  )
!*      deallocate (cam_longxy  )
!*      deallocate (cam_latixy  )
!*      deallocate (cam_landfrac)
!*      deallocate (cam_landmask)
!*#endif

   end subroutine readgridata
!======================================
! for IAPDGVM
!======================================
#if(defined IAPDGVM)
   subroutine readigdata

      use spmd
      use colm_varMod
      use netcdf
      implicit none

      integer nlon,nlat,fid,vid,vid_t

if(p_master) then
         write(6,*) 'Opening ig file...'

         call sanity(nf90_open(path=trim(fnig),mode=nf90_nowrite,ncid=fid))
         call sanity(nf90_inq_dimid(fid,'lsmlon',vid_t))
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nlon))
         call sanity(nf90_inq_dimid(fid,'lsmlat',vid_t))
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nlat))
         if(nlat.ne.lat_points .or. nlon.ne.lon_points) then
            write(6,*) 'Reading ig forcing data error', nlat, nlon
            call abort
         end if
         call sanity(nf90_inq_varid(fid,'ig',vid))
         call sanity(nf90_get_var(fid,vid,iglf,start=(/1,1,1/),count=(/lon_points,lat_points,365*8/)))
         call sanity(nf90_close(fid))
end if

#if(defined SPMD)
      call mpi_bcast(iglf,size(iglf),mpi_real8,0,p_comm,p_err)
#endif
   end subroutine readigdata

function getig(c,year,jday,msec)

      use colm_varMod
      use spmd_decomp, only: cgmap, gxmap, gymap

     implicit none

      integer, intent(in) :: c
      integer, intent(in) :: year
      integer, intent(in) :: jday
      integer, intent(in) :: msec

      real(r8) getig
      integer ind, lon, lat

      ind = (jday-1)*8 + msec/(3*3600) + 1

      ind = min(ind, 365*8)

      lon = gxmap(cgmap(c))
      lat = gymap(cgmap(c))
      
      getig = iglf(lon,lat,ind)
  end function getig


function getrhm(tm,pbot,qm)

      use colm_varMod

      implicit none

      real(r8), intent(in) :: tm
      real(r8), intent(in) :: pbot
      real(r8), intent(in) :: qm

      real(r8) getrhm
      real(r8) es, esdT, qs, qsdT

      call qsadv(tm,pbot,es,esdT,qs,qsdT)

      getrhm = min(qm/qs,1._r8)

   end function getrhm
subroutine sanity(ret)
      use netcdf
      implicit none

      integer, intent(in) :: ret

      if(ret /= nf90_noerr) then
         write(6,*) trim(nf90_strerror(ret))
         stop
      endif

   end subroutine sanity

#endif 
!end of IAPDGVM
!======================================


   subroutine readinidata

#ifdef RTM
      use RtmMod, only: restart_rtm
#endif
#ifdef COUP_CSM
      use colm_cplMod
#ifdef CPL6
      use colm_csmMod, only: csm_restart
#endif
#endif

      implicit none

! local variables:

      real(r8), pointer :: lndxy_avsdr(:,:) ! avsdr
      real(r8), pointer :: lndxy_avsdf(:,:) ! avsdf
      real(r8), pointer :: lndxy_anidr(:,:) ! anidr
      real(r8), pointer :: lndxy_anidf(:,:) ! anidf
      real(r8), pointer :: lndxy_scv(:,:)   ! snow water [mm]
      real(r8), pointer :: lndxy_trad(:,:)  ! radiative temp. [K]

#ifdef MYBUG
      real(r4), pointer :: frac(:,:,:)
#endif

      integer i, j, k, p, c, g

! routine:

#ifdef MYBUG
      write(6,*), p_iam, 'before read restart'
#endif

      if (p_master) then
         if(lini) frestart = trim(fini)
         if(nsrest.gt.0) call read_restart_pfile(frestart)

       ! Open for model time varying data (model state variables)
         OPEN(unit=lurestart,file=trim(frestart),form='unformatted',status='old',action='read')

         read(lurestart) istep, idate, idate_p

         if(lini) write(6,*) 'CoLM initialized time', idate
         if(nsrest.gt.0) write(6,*) 'CoLM restart time', idate
      end if

#ifdef SPMD
      call mpi_bcast (istep,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (idate,size(idate),mpi_integer,0,p_comm,p_err)
      call mpi_bcast (idate_p,size(idate_p),mpi_integer,0,p_comm,p_err)
#endif

      call timeinit(idate_p,idate,istep)

      if(startup_date(3) /= -1) then
#ifdef COUP_CSM
         if(nsrest == 0) then
            istep = 0  ! clear random istep in restart file
       !    call set_curr_date(startup_date)
       !    idate_p = startup_date
            call timeinit(startup_date,startup_date,istep)
         end if
#else
         istep = 0  ! clear random istep in restart file
       ! call set_curr_date(startup_date)
       ! idate_p = startup_date
         call timeinit(startup_date,startup_date,istep)
#endif
      end if

      if (p_master) read(lurestart) ftune ! clm tunable constants
#ifdef SPMD
      call mpi_bcast (ftune,size(ftune),mpi_real8,0,p_comm,p_err)
#endif

      call readc(lurestart,cvar%dlat)
      call readc(lurestart,cvar%dlon)
      call readc(lurestart,cvar%itypwat)
      call readc(lurestart,cvar%wxy_column)
      call readc(lurestart,cvar%albsat,2)
      call readc(lurestart,cvar%albdry,2)
      call readc(lurestart,cvar%csol,nl_soil)
      call readc(lurestart,cvar%porsl,nl_soil)
      call readc(lurestart,cvar%phi0,nl_soil)
      call readc(lurestart,cvar%bsw,nl_soil)
      call readc(lurestart,cvar%dkmg,nl_soil)
      call readc(lurestart,cvar%dksatu,nl_soil)
      call readc(lurestart,cvar%dkdry,nl_soil)
      call readc(lurestart,cvar%hksati,nl_soil)

      call readp(lurestart,pvar%ivt)
      call readp(lurestart,pvar%wxy_patch)
      call readp(lurestart,pvar%z0m)
      call readp(lurestart,pvar%displa)
      call readp(lurestart,pvar%sqrtdi)
      call readp(lurestart,pvar%effcon)
      call readp(lurestart,pvar%vmax25)
      call readp(lurestart,pvar%slti)
      call readp(lurestart,pvar%hlti)
      call readp(lurestart,pvar%shti)
      call readp(lurestart,pvar%hhti)
      call readp(lurestart,pvar%trda)
      call readp(lurestart,pvar%trdm)
      call readp(lurestart,pvar%trop)
      call readp(lurestart,pvar%gradm)
      call readp(lurestart,pvar%binter)
      call readp(lurestart,pvar%extkn)
      call readp(lurestart,pvar%chil)
      call readp(lurestart,pvar%refl,2,2)
      call readp(lurestart,pvar%refs,2,2)
      call readp(lurestart,pvar%tranl,2,2)
      call readp(lurestart,pvar%trans,2,2)
      call readp(lurestart,pvar%rootfr,nl_soil)
#ifdef DGVM
      call readp(lurestart,pvar%pftpara,npftpara)
      call readp(lurestart,pvar%vegclass)
      call readp(lurestart,pvar%summergreen)
      call readp(lurestart,pvar%raingreen)
      call readp(lurestart,pvar%sla)
      call readp(lurestart,pvar%lm_sapl)
      call readp(lurestart,pvar%sm_sapl)
      call readp(lurestart,pvar%hm_sapl)
      call readp(lurestart,pvar%rm_sapl)
#ifdef DyN
      call readp(lurestart,pvar%cton_soil)
      call readp(lurestart,pvar%cton_pro)
#endif
#endif

      call readc(lurestart,cvar%z,(nl_soil-maxsnl))
      call readc(lurestart,cvar%dz,(nl_soil-maxsnl))
      call readc(lurestart,cvar%tss,(nl_soil-maxsnl))
      call readc(lurestart,cvar%wliq,(nl_soil-maxsnl))
      call readc(lurestart,cvar%wice,(nl_soil-maxsnl))
      call readc(lurestart,cvar%tg)
      call readc(lurestart,cvar%sag)
      call readc(lurestart,cvar%scv)
      call readc(lurestart,cvar%snowdp)
      call readc(lurestart,cvar%fsno)
      call readc(lurestart,cvar%coszen)
      call readc(lurestart,cvar%lakedepth)
      call readc(lurestart,cvar%dz_lake,     nl_lake)
      call readc(lurestart,cvar%t_lake,      nl_lake)
      call readc(lurestart,cvar%lake_icefrac,nl_lake)
      call readc(lurestart,cvar%savedtke1)
#ifdef DGVM
      call readc(lurestart,cvar%nday)
      call readc(lurestart,cvar%nyr)
      call readc(lurestart,cvar%prec365)
#ifdef IAPDGVM
      call readc(lurestart,cvar%wliq6mon)
#endif
#ifdef DyN
      call readc(lurestart,cvar%soil_no3)
      call readc(lurestart,cvar%soil_no2)
      call readc(lurestart,cvar%soil_no)
      call readc(lurestart,cvar%soil_n2o)
      call readc(lurestart,cvar%soil_n2)
      call readc(lurestart,cvar%soil_nh4)
#endif
#endif


#if (defined FHNP) && (defined FTF)
!liruichao add     
      call readc(lurestart,cvar%frostdp)
      call readc(lurestart,cvar%thawdp)
      call readc(lurestart,cvar%frostdp0)
      call readc(lurestart,cvar%D_temperature)
      call readc(lurestart,cvar%N_time)
      call readc(lurestart,cvar%frost_day)
      call readc(lurestart,cvar%thaw_day)
!end
#endif


      call readp(lurestart,pvar%tlsun)
      call readp(lurestart,pvar%tlsha)
      call readp(lurestart,pvar%ldew)
      call readp(lurestart,pvar%fveg)
      call readp(lurestart,pvar%sigf)
      call readp(lurestart,pvar%green)
      call readp(lurestart,pvar%lai)
      call readp(lurestart,pvar%sai)
      call readp(lurestart,pvar%albg,2,2)
      call readp(lurestart,pvar%albv,2,2)
      call readp(lurestart,pvar%alb,2,2)
      call readp(lurestart,pvar%ssun,2,2)
      call readp(lurestart,pvar%ssha,2,2)
      call readp(lurestart,pvar%thermk)
      call readp(lurestart,pvar%extkb)
      call readp(lurestart,pvar%extkd)

 ! Additional variables required by reginal model
      call readp(lurestart,pvar%trad)
      call readp(lurestart,pvar%tref)
      call readp(lurestart,pvar%qref)
      call readp(lurestart,pvar%rst)
      call readp(lurestart,pvar%emis)
      call readp(lurestart,pvar%z0ma)
      call readp(lurestart,pvar%zol)
      call readp(lurestart,pvar%rib)
      call readp(lurestart,pvar%ustar)
      call readp(lurestart,pvar%qstar)
      call readp(lurestart,pvar%tstar)
      call readp(lurestart,pvar%fm)
      call readp(lurestart,pvar%fh)
      call readp(lurestart,pvar%fq)

#ifdef DGVM
      call readp(lurestart,pvar%t10min)
      call readp(lurestart,pvar%lai_ind)
      call readp(lurestart,pvar%dphen)
      call readp(lurestart,pvar%leafon)
      call readp(lurestart,pvar%leafof)
      call readp(lurestart,pvar%firelength)
#ifdef IAPDGVM
      call readp(lurestart,pvar%afirefrac1)
      call readp(lurestart,pvar%nfireg1)
#endif

      call readp(lurestart,pvar%litterag)
      call readp(lurestart,pvar%litterbg)
      call readp(lurestart,pvar%cpool_fast,nl_soil)
      call readp(lurestart,pvar%cpool_slow,nl_soil)
      call readp(lurestart,pvar%k_fast_ave)
      call readp(lurestart,pvar%k_slow_ave)
      call readp(lurestart,pvar%litter_decom_ave)
      call readp(lurestart,pvar%fmicr)
      call readp(lurestart,pvar%nind)
      call readp(lurestart,pvar%lm_ind)
      call readp(lurestart,pvar%sm_ind)
      call readp(lurestart,pvar%hm_ind)
      call readp(lurestart,pvar%rm_ind)
      call readp(lurestart,pvar%tmomin20)
      call readp(lurestart,pvar%agdd0)
      call readp(lurestart,pvar%agdd)
      call readp(lurestart,pvar%agdd20)
      call readp(lurestart,pvar%t_mo_min)
      call readp(lurestart,pvar%crownarea)
      call readp(lurestart,pvar%htop)
      call readp(lurestart,pvar%tsai)
      call readp(lurestart,pvar%fpcgrid)
      call readp(lurestart,pvar%bm_inc)
      call readp(lurestart,pvar%afmicr)
      call readp(lurestart,pvar%annpsn)
      call readp(lurestart,pvar%annpsnpot)
      call readp(lurestart,pvar%tref10)
      call readp(lurestart,pvar%tref_sum)
      call readp(lurestart,pvar%t10,10)
      call readp(lurestart,pvar%assimn10)
      call readp(lurestart,pvar%assimn_sum)
      call readp(lurestart,pvar%an10,10)
      call readp(lurestart,pvar%turnover_ind)
      call readp(lurestart,pvar%fpc_inc)
      call readp(lurestart,pvar%agddtw)
      call readp(lurestart,pvar%ifpre)
      call readp(lurestart,pvar%t_mo)
      call readp(lurestart,pvar%t_mo_sum)
#ifdef DyN
      call readp(lurestart,pvar%litter_leaf)
      call readp(lurestart,pvar%litter_wood)
      call readp(lurestart,pvar%litter_root)
      call readp(lurestart,pvar%litter_repr)
      call readp(lurestart,pvar%litter_leaf_n)
      call readp(lurestart,pvar%litter_wood_n)
      call readp(lurestart,pvar%litter_root_n)
      call readp(lurestart,pvar%litter_repr_n)
      call readp(lurestart,pvar%afcton_leaf)
      call readp(lurestart,pvar%afcton_sap)
      call readp(lurestart,pvar%afcton_root)
      call readp(lurestart,pvar%lm_ind_n)
      call readp(lurestart,pvar%sm_ind_n)
      call readp(lurestart,pvar%hm_ind_n)
      call readp(lurestart,pvar%rm_ind_n)
      call readp(lurestart,pvar%an_up
      call readp(lurestart,pvar%an_stress)
#endif 
#endif


#ifdef MYBUG
      write(6,*), p_iam, 'after read restart'
#endif

#ifdef RTM
      if (nsrest.gt.0) call restart_rtm(lurestart,'read')
#endif

#ifdef COUP_CSM
#ifdef CPL6
      if (nsrest.gt.0) call csm_restart(lurestart,'read') 
#endif
#endif

      if (p_master) CLOSE(lurestart)

#ifdef COUP_CSM
      if (nsrest.eq.0) then   ! initial run

         allocate (lndxy_avsdr(lon_points,lat_points))
         allocate (lndxy_avsdf(lon_points,lat_points))
         allocate (lndxy_anidr(lon_points,lat_points))
         allocate (lndxy_anidf(lon_points,lat_points))
         allocate (lndxy_trad (lon_points,lat_points))
         allocate (lndxy_scv  (lon_points,lat_points))

         if (p_master) then
            OPEN(unit=lusbc,file=fsbc,form='unformatted',status='old',action='read')
            read(lusbc) lndxy_avsdr(:,:) ! avsdr
            read(lusbc) lndxy_avsdf(:,:) ! avsdf
            read(lusbc) lndxy_anidr(:,:) ! anidr
            read(lusbc) lndxy_anidf(:,:) ! anidf
            read(lusbc) lndxy_trad
            read(lusbc) lndxy_scv
            CLOSE(lusbc)
         end if

#ifdef MYBUG
         write(6,*), p_iam, 'after sbc read'
#endif

#ifdef SPMD
         call mpi_bcast (lndxy_avsdr,size(lndxy_avsdr),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_avsdf,size(lndxy_avsdf),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_anidr,size(lndxy_anidr),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_anidf,size(lndxy_anidf),mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_trad ,size(lndxy_trad) ,mpi_real8,0,p_comm,p_err)
         call mpi_bcast (lndxy_scv  ,size(lndxy_scv)  ,mpi_real8,0,p_comm,p_err)
#endif

#ifdef MYBUG
         write(6,*), p_iam, 'after sbc bcast'
#endif

         do g = 1,numgrid
            i = gxmap(g)
            j = gymap(g)

            lnd_avsdr(g) = lndxy_avsdr(i,j)
            lnd_avsdf(g) = lndxy_avsdf(i,j)
            lnd_anidr(g) = lndxy_anidr(i,j)
            lnd_anidf(g) = lndxy_anidf(i,j)
            lnd_trad(g)  = lndxy_trad(i,j)
            lnd_scv(g)   = lndxy_scv(i,j)/1000._r8
         end do

#ifdef MYBUG
         write(6,*), 'min & max avsdr', minval(lnd_avsdr), maxval(lnd_avsdr)
         write(6,*), 'min & max avsdf', minval(lnd_avsdf), maxval(lnd_avsdf)
         write(6,*), 'min & max anidr', minval(lnd_anidr), maxval(lnd_anidr)
         write(6,*), 'min & max anidf', minval(lnd_anidf), maxval(lnd_anidf)
         write(6,*), 'min & max trad', minval(lnd_trad), maxval(lnd_trad)
         write(6,*), 'min & max scv', minval(lnd_scv), maxval(lnd_scv)
#endif

         deallocate (lndxy_avsdr)
         deallocate (lndxy_avsdf)
         deallocate (lndxy_anidr)
         deallocate (lndxy_anidf)
         deallocate (lndxy_trad)
         deallocate (lndxy_scv)
      end if
#endif

   end subroutine readinidata

#ifdef CPL7
   subroutine writehistdata(rstwr,nlend,rdate)
#else
   subroutine writehistdata
#endif

      use spmd, only: p_master
      use timemgr, only: idate
#ifdef RTM
      use RtmMod, only: restart_rtm
      use colm_rtmMod
#endif
#ifdef COUP_CSM
#ifdef CPL6
      use colm_csmMod, only: csm_restart
#endif
#endif
      implicit none

! arguments:
#ifdef CPL7
      logical,         intent(in) :: rstwr       ! true => write restart file this step
      logical,         intent(in) :: nlend       ! true => end of run on this step
      character(len=*),intent(in) :: rdate       ! restart file time stamp for name
#endif

! local variables:

      integer n, k
      character(len=255) :: cdate

! routine:

      frestart = "null"

      call nchist_update()

    ! this subroutine should be of module nchistMod.
      call lpwrite(foutdir,frestart,cdate)

      call nchist_write(idate,cdate)

#ifdef CPL7
      if (rstwr) then
         if (p_master) frestart = trim(foutdir)//'/'//trim(caseid)//'-colm-restart-'//trim(rdate)
         if (p_master) write(6,*) "writing restart: "//trim(frestart)
#else
      if (trim(frestart).ne."null") then
#endif
         if (p_master) OPEN(unit=lurestart,file=frestart,access='sequential', &
                            form='unformatted',status='unknown')

         call write_colm_restart

#ifdef RTM
         call restart_rtm(lurestart,'write')
#endif

#ifdef COUP_CSM
#ifdef CPL6
         if (p_master) call csm_restart(lurestart,'write')
#endif
#endif

         if (p_master) CLOSE(lurestart)

         call write_colm_sbc

         if (p_master) then
            call write_restart_nml(frestart)
            call write_restart_pfile(frestart)
         end if
      end if

    ! if((idate(1)*1000+idate(2)).gt.7360) lcycdump = .true.
      if (lcycdump) call cycdump

   end subroutine writehistdata

 ! Write out the model variables for restart run 
   subroutine write_colm_restart

      use spmd
      use spmd_decomp
      use timemgr, only: istep, idate, idate_p
      use colm_varMod
      implicit none

      if (p_master) write(lurestart) istep, idate, idate_p

      if (p_master) write(lurestart) ftune ! clm tunable constants

      call writec(lurestart,cvar%dlat)
      call writec(lurestart,cvar%dlon)
      call writec(lurestart,cvar%itypwat)
      call writec(lurestart,cvar%wxy_column)
      call writec(lurestart,cvar%albsat,2)
      call writec(lurestart,cvar%albdry,2)
      call writec(lurestart,cvar%csol,nl_soil)
      call writec(lurestart,cvar%porsl,nl_soil)
      call writec(lurestart,cvar%phi0,nl_soil)
      call writec(lurestart,cvar%bsw,nl_soil)
      call writec(lurestart,cvar%dkmg,nl_soil)
      call writec(lurestart,cvar%dksatu,nl_soil)
      call writec(lurestart,cvar%dkdry,nl_soil)
      call writec(lurestart,cvar%hksati,nl_soil)

      call writep(lurestart,pvar%ivt)
      call writep(lurestart,pvar%wxy_patch)
      call writep(lurestart,pvar%z0m)
      call writep(lurestart,pvar%displa)
      call writep(lurestart,pvar%sqrtdi)
      call writep(lurestart,pvar%effcon)
      call writep(lurestart,pvar%vmax25)
      call writep(lurestart,pvar%slti)
      call writep(lurestart,pvar%hlti)
      call writep(lurestart,pvar%shti)
      call writep(lurestart,pvar%hhti)
      call writep(lurestart,pvar%trda)
      call writep(lurestart,pvar%trdm)
      call writep(lurestart,pvar%trop)
      call writep(lurestart,pvar%gradm)
      call writep(lurestart,pvar%binter)
      call writep(lurestart,pvar%extkn)
      call writep(lurestart,pvar%chil)
      call writep(lurestart,pvar%refl,2,2)
      call writep(lurestart,pvar%refs,2,2)
      call writep(lurestart,pvar%tranl,2,2)
      call writep(lurestart,pvar%trans,2,2)
      call writep(lurestart,pvar%rootfr,nl_soil)
#ifdef DGVM
      call writep(lurestart,pvar%pftpara,npftpara)
      call writep(lurestart,pvar%vegclass)
      call writep(lurestart,pvar%summergreen)
      call writep(lurestart,pvar%raingreen)
      call writep(lurestart,pvar%sla)
      call writep(lurestart,pvar%lm_sapl)
      call writep(lurestart,pvar%sm_sapl)
      call writep(lurestart,pvar%hm_sapl)
      call writep(lurestart,pvar%rm_sapl)
#ifdef DyN
      call writep(lurestart,pvar%cton_soil)
      call writep(lurestart,pvar%cton_pro)
#endif
#endif

      call writec(lurestart,cvar%z,(nl_soil-maxsnl))
      call writec(lurestart,cvar%dz,(nl_soil-maxsnl))
      call writec(lurestart,cvar%tss,(nl_soil-maxsnl))
      call writec(lurestart,cvar%wliq,(nl_soil-maxsnl))
      call writec(lurestart,cvar%wice,(nl_soil-maxsnl))
      call writec(lurestart,cvar%tg)
      call writec(lurestart,cvar%sag)
      call writec(lurestart,cvar%scv)
      call writec(lurestart,cvar%snowdp)
      call writec(lurestart,cvar%fsno)
      call writec(lurestart,cvar%coszen)
      call writec(lurestart,cvar%lakedepth)
      call writec(lurestart,cvar%dz_lake,     nl_lake)
      call writec(lurestart,cvar%t_lake,      nl_lake)
      call writec(lurestart,cvar%lake_icefrac,nl_lake)
      call writec(lurestart,cvar%savedtke1)
#ifdef DGVM
      call writec(lurestart,cvar%nday)
      call writec(lurestart,cvar%nyr)
      call writec(lurestart,cvar%prec365)
#ifdef IAPDGVM
      call writec(lurestart,cvar%wliq6mon)
#endif

#ifdef DyN
      call writec(lurestart,cvar%soil_no3)
      call writec(lurestart,cvar%soil_no2)
      call writec(lurestart,cvar%soil_no)
      call writec(lurestart,cvar%soil_n2o)
      call writec(lurestart,cvar%soil_n2)
      call writec(lurestart,cvar%soil_nh4)
#endif
#endif

#if (defined FHNP) && (defined FTF)
!liruichao
      call writec(lurestart,cvar%frostdp)
      call writec(lurestart,cvar%thawdp)
      call writec(lurestart,cvar%frostdp0)
      call writec(lurestart,cvar%D_temperature)
      call writec(lurestart,cvar%N_time)
      call writec(lurestart,cvar%frost_day)
      call writec(lurestart,cvar%thaw_day)
!end
#endif
      call writep(lurestart,pvar%tlsun)
      call writep(lurestart,pvar%tlsha)
      call writep(lurestart,pvar%ldew)
      call writep(lurestart,pvar%fveg)
      call writep(lurestart,pvar%sigf)
      call writep(lurestart,pvar%green)
      call writep(lurestart,pvar%lai)
      call writep(lurestart,pvar%sai)
      call writep(lurestart,pvar%albg,2,2)
      call writep(lurestart,pvar%albv,2,2)
      call writep(lurestart,pvar%alb,2,2)
      call writep(lurestart,pvar%ssun,2,2)
      call writep(lurestart,pvar%ssha,2,2)
      call writep(lurestart,pvar%thermk)
      call writep(lurestart,pvar%extkb)
      call writep(lurestart,pvar%extkd)

 ! Additional variables required by reginal model
      call writep(lurestart,pvar%trad)
      call writep(lurestart,pvar%tref)
      call writep(lurestart,pvar%qref)
      call writep(lurestart,pvar%rst)
      call writep(lurestart,pvar%emis)
      call writep(lurestart,pvar%z0ma)
      call writep(lurestart,pvar%zol)
      call writep(lurestart,pvar%rib)
      call writep(lurestart,pvar%ustar)
      call writep(lurestart,pvar%qstar)
      call writep(lurestart,pvar%tstar)
      call writep(lurestart,pvar%fm)
      call writep(lurestart,pvar%fh)
      call writep(lurestart,pvar%fq)

#ifdef DGVM
      call writep(lurestart,pvar%t10min)
      call writep(lurestart,pvar%lai_ind)
      call writep(lurestart,pvar%dphen)
      call writep(lurestart,pvar%leafon)
      call writep(lurestart,pvar%leafof)
      call writep(lurestart,pvar%firelength)
#ifdef IAPDGVM
      call writep(lurestart,pvar%afirefrac1)
      call writep(lurestart,pvar%nfireg1)
#endif

      call writep(lurestart,pvar%litterag)
      call writep(lurestart,pvar%litterbg)
      call writep(lurestart,pvar%cpool_fast,nl_soil)
      call writep(lurestart,pvar%cpool_slow,nl_soil)
      call writep(lurestart,pvar%k_fast_ave)
      call writep(lurestart,pvar%k_slow_ave)
      call writep(lurestart,pvar%litter_decom_ave)
      call writep(lurestart,pvar%fmicr)
      call writep(lurestart,pvar%nind)
      call writep(lurestart,pvar%lm_ind)
      call writep(lurestart,pvar%sm_ind)
      call writep(lurestart,pvar%hm_ind)
      call writep(lurestart,pvar%rm_ind)
      call writep(lurestart,pvar%tmomin20)
      call writep(lurestart,pvar%agdd0)
      call writep(lurestart,pvar%agdd)
      call writep(lurestart,pvar%agdd20)
      call writep(lurestart,pvar%t_mo_min)
      call writep(lurestart,pvar%crownarea)
      call writep(lurestart,pvar%htop)
      call writep(lurestart,pvar%tsai)
      call writep(lurestart,pvar%fpcgrid)
      call writep(lurestart,pvar%bm_inc)
      call writep(lurestart,pvar%afmicr)
      call writep(lurestart,pvar%annpsn)
      call writep(lurestart,pvar%annpsnpot)
      call writep(lurestart,pvar%tref10)
      call writep(lurestart,pvar%tref_sum)
      call writep(lurestart,pvar%t10,10)
      call writep(lurestart,pvar%assimn10)
      call writep(lurestart,pvar%assimn_sum)
      call writep(lurestart,pvar%an10,10)
      call writep(lurestart,pvar%turnover_ind)
      call writep(lurestart,pvar%fpc_inc)
      call writep(lurestart,pvar%agddtw)
      call writep(lurestart,pvar%ifpre)
      call writep(lurestart,pvar%t_mo)
      call writep(lurestart,pvar%t_mo_sum)

#ifdef DyN
      call writep(lurestart,pvar%litter_leaf)
      call writep(lurestart,pvar%litter_wood)
      call writep(lurestart,pvar%litter_root)
      call writep(lurestart,pvar%litter_repr)
      call writep(lurestart,pvar%litter_leaf_n)
      call writep(lurestart,pvar%litter_wood_n)
      call writep(lurestart,pvar%litter_root_n)
      call writep(lurestart,pvar%litter_repr_n)
      call writep(lurestart,pvar%afcton_leaf)
      call writep(lurestart,pvar%afcton_sap)
      call writep(lurestart,pvar%afcton_root)
      call writep(lurestart,pvar%lm_ind_n)
      call writep(lurestart,pvar%sm_ind_n)
      call writep(lurestart,pvar%hm_ind_n)
      call writep(lurestart,pvar%rm_ind_n)
      call writep(lurestart,pvar%an_up
      call writep(lurestart,pvar%an_stress)
#endif 
#endif

   end subroutine write_colm_restart

   subroutine write_colm_sbc

      use spmd
      use spmd_decomp
      use colm_varMod

      implicit none

      real(r8), pointer :: fldv_glob(:,:)
      real(r8), pointer :: fldxy(:,:)

      integer i, j, g, L

      allocate (fldv_glob(numgrid_glob,6))
      allocate (fldxy(lon_points,lat_points))

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do
#endif

#ifdef SPMD
      call mpi_gatherv (fldv%avsdr,size(fldv%avsdr),mpi_real8,&
                        fldv_glob(:,1),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
      call mpi_gatherv (fldv%avsdf,size(fldv%avsdf),mpi_real8,&
                        fldv_glob(:,2),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
      call mpi_gatherv (fldv%anidr,size(fldv%anidr),mpi_real8,&
                        fldv_glob(:,3),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
      call mpi_gatherv (fldv%anidf,size(fldv%anidf),mpi_real8,&
                        fldv_glob(:,4),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
      call mpi_gatherv (fldv%trad ,size(fldv%trad ),mpi_real8,&
                        fldv_glob(:,5),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
      call mpi_gatherv (fldv%scv  ,size(fldv%scv  ),mpi_real8,&
                        fldv_glob(:,6),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      fldv_glob(:,1) = fldv%avsdr
      fldv_glob(:,2) = fldv%avsdf
      fldv_glob(:,3) = fldv%anidr
      fldv_glob(:,4) = fldv%anidf
      fldv_glob(:,5) = fldv%trad
      fldv_glob(:,6) = fldv%scv
#endif

      if (p_master) then
         open(lusbc,file=trim(frestart)//"-sbc",form="unformatted",status="unknown")

         do L = 1, 6
            fldxy(:,:) = -9999.

            do g = 1, numgrid_glob
               i = gxmap_glob(g)
               j = gymap_glob(g)
               fldxy(i,j) = fldv_glob(g,L)
            end do

            write(lusbc) fldxy
         end do

         close(lusbc)
      end if

      deallocate (fldxy)
      deallocate (fldv_glob)

   end subroutine write_colm_sbc

   subroutine write_restart_nml(frestart)

      use timemgr
      use colm_varMod
#ifdef RTM
      use colm_rtmVar
      use colm_rtmMod
#endif
#ifdef COUP_CSM
      use colm_cplMod, only: irad, nsrest
#endif

      implicit none

      character(len=255), intent(in) :: frestart

      integer :: lurestnml = 111

      open(unit=lurestnml,file="lnd.stdin.rst",form="formatted",status="unknown")

      write(lurestnml,10) "&clmexp"
      write(lurestnml,20) "caseid = ", "'"//trim(caseid)//"'"
      write(lurestnml,20) "fsrf = ", "'"//trim(fsrf)//"'"
#if(defined VEGDATA)
      write(lurestnml,20) "flai = ", "'"//trim(flai)//"'"
#endif
#ifdef OFFLINE
      write(lurestnml,20) "fmet = ", "'"//trim(fmet)//"'"
      write(lurestnml,20) "fco2 = ", "'"//trim(fco2)//"'"
#endif
      write(lurestnml,20) "foutdir = ", "'"//trim(foutdir)//"'"
#ifdef IAPDGVM
      write(lurestnml,20) "fnig = ", "'"//trim(fnig)//"'"
#endif
      write(lurestnml,20) "frestart = ", "'"//trim(frestart)//"'"
#ifdef RTM
      write(lurestnml,20) "frivinp_rtm = ", "'"//trim(frivinp_rtm)//"'"
#endif
#ifdef DUST
      write(lurestnml,20) "fdust = ", "'"//trim(fdust)//"'"
#endif
#ifdef COUP_CSM
      write(lurestnml,30) "fsbc = "//"''"
      write(lurestnml,30) "irad = ", irad
      write(lurestnml,30) "nsrest = ", 1
      write(lurestnml,60) "csm_doflxave = ", csm_doflxave
      write(lurestnml,40) "startup_date = ", idate(1), idate(2), idate(3)
      write(lurestnml,30) "lnd_cflux_year = ", lnd_cflux_year
#ifdef CPL6
      write(lurestnml,20) "co2_option = ", "'"//trim(co2_option)//"'"
#endif
#ifdef CPL7
      write(lurestnml,20) "co2_type = ", "'"//trim(co2_type)//"'"
      write(lurestnml,50) "co2_ppmv = ", co2_ppmv
#endif
#endif
      write(lurestnml,60) "lhist_yearly = ", lhist_yearly
      write(lurestnml,60) "lhist_monthly = ", lhist_monthly
      write(lurestnml,60) "lhist_daily = ", lhist_daily
      write(lurestnml,60) "lhist_3hourly = ", lhist_3hourly
      write(lurestnml,20) "restart_freq = ", "'"//trim(restart_freq)//"'"
      write(lurestnml,30) "lon_points = ", lon_points
      write(lurestnml,30) "lat_points = ", lat_points
      write(lurestnml,50) "dtime = ", dtime
      write(lurestnml,30) "mstep = ", mstep
      write(lurestnml,10) "/"

10 format(A)
20 format(A20,A)
30 format(A20,I10)
40 format(A20,I4,',',I3,',',I5)
50 format(A20,F10.2)
60 format(A20,L)

      close(lurestnml)

   end subroutine write_restart_nml

   subroutine write_restart_pfile(frestart)

      implicit none

      character(len=255), intent(in) :: frestart

      integer :: lupfile = 111

      open(unit=lupfile,file="rpointer.lnd",form="formatted",status="unknown")
      write(lupfile,'(a)') frestart
      close(lupfile)

      write(6,*) 'Successfully wrote local restart pointer file (lnd)'

   end subroutine write_restart_pfile

   subroutine read_restart_pfile(frestart)

      implicit none

      character(len=255), intent(out) :: frestart

      integer :: lupfile = 111
      integer i

      open(unit=lupfile,file="rpointer.lnd",form="formatted",status="old")
      read(lupfile,'(a255)') frestart
      close(lupfile)

      write(6,*) 'Reading restart data.....' // trim(frestart)
      write(6,'(72a1)') ("-",i=1,60)

   end subroutine read_restart_pfile

   subroutine cycdump

      use spmd
      use spmd_decomp
      use timemgr, only: istep, idate, idate_p
      use colm_varMod

      implicit none

      real(r8), pointer :: forc_loc(:)
      real(r8), pointer :: forc_glob(:,:)
      real(r8), pointer :: forc_buf(:)

      real(r8), pointer :: fvar_loc(:)
      real(r8), pointer :: fvar_glob(:,:)
      real(r8), pointer :: fvar_buf(:)

      integer i, L

      character(len=255) :: fcycdumpI

    ! forcing variables

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))
      end do
#endif

      allocate (forc_buf(numcolumn_glob))
      allocate (forc_glob(numcolumn_glob,20))

      do L = 1, 20
         if (L== 1) forc_loc => forc%pco2m
         if (L== 2) forc_loc => forc%po2m
         if (L== 3) forc_loc => forc%us
         if (L== 4) forc_loc => forc%vs
         if (L== 5) forc_loc => forc%tm
         if (L== 6) forc_loc => forc%qm
         if (L== 7) forc_loc => forc%prc
         if (L== 8) forc_loc => forc%prl
         if (L== 9) forc_loc => forc%rain
         if (L==10) forc_loc => forc%snow
         if (L==11) forc_loc => forc%pbot
         if (L==12) forc_loc => forc%psrf
         if (L==13) forc_loc => forc%sols
         if (L==14) forc_loc => forc%soll
         if (L==15) forc_loc => forc%solsd
         if (L==16) forc_loc => forc%solld
         if (L==17) forc_loc => forc%frl
         if (L==18) forc_loc => forc%hu
         if (L==19) forc_loc => forc%ht
         if (L==20) forc_loc => forc%hq
#ifdef SPMD
         call mpi_gatherv (forc_loc,size(forc_loc),mpi_real8,&
                           forc_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         forc_glob(:,L) = forc_loc(:)
#endif
      end do

      if (p_master) then
         write(fcycdumpI,'(I3.3)') mod(istep,ncycdump)
         fcycdumpI = trim(fcycdump)//"-"//trim(fcycdumpI)

         OPEN(unit=lucycdump,file=fcycdumpI,access='sequential', &
                             form='unformatted',status='unknown')
         write(lucycdump) istep, idate, idate_p

         do L = 1, nforc
            forc_buf(ccmap_glob) = forc_glob(:,L)
            write(lucycdump) forc_buf
         end do
      end if

      deallocate (forc_buf)
      deallocate (forc_glob)

!   ! column level variables

!     allocate (fvar_buf(numcolumn_glob))
!     allocate (fvar_loc(numcolumn))
!     allocate (fvar_glob(numcolumn_glob,nfvar_col))

#ifdef SPMD
!     p_counts(:) = numcolumn_proc(:)
!     p_displs(0) = 0
!     do i = 1, p_nprocs-1
!        p_displs(i) = sum(numcolumn_proc(0:i-1))
!     end do
#endif

!     do L = 1, nfvar_col
#ifdef SPMD
!        fvar_loc(:) = fvar_col(L,:)
!        call mpi_gatherv (fvar_loc,size(fvar_loc),mpi_real8,&
!                          fvar_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
!        fvar_glob(:,L) = fvar_col(L,:)
#endif
!     end do

!     if (p_master) then
!        do L = 1, nfvar_col
!           fvar_buf(ccmap_glob) = fvar_glob(:,L)
!           write(lucycdump) fvar_buf
!        end do
!     end if

!     deallocate (fvar_buf)
!     deallocate (fvar_loc)
!     deallocate (fvar_glob)

!   ! patch level variables

!     allocate (fvar_buf(numpatch_glob))
!     allocate (fvar_loc(numpatch))
!     allocate (fvar_glob(numpatch_glob,nfvar_pft))

#ifdef SPMD
!     p_counts(:) = numpatch_proc(:)
!     p_displs(0) = 0
!     do i = 1, p_nprocs-1
!        p_displs(i) = sum(numpatch_proc(0:i-1))
!     end do
#endif

!     do L = 1, nfvar_pft
#ifdef SPMD
!        fvar_loc(:) = fvar_pft(L,:)
!        call mpi_gatherv (fvar_loc,size(fvar_loc),mpi_real8,&
!                          fvar_glob(:,L),p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
!        fvar_glob(:,L) = fvar_pft(L,:)
#endif
!     end do

!     if (p_master) then
!        do L = 1, nfvar_pft
!           fvar_buf(ppmap_glob) = fvar_glob(:,L)
!           write(lucycdump) fvar_buf
!        end do

!        CLOSE(lucycdump)
!     end if

!     deallocate (fvar_buf)
!     deallocate (fvar_loc)
!     deallocate (fvar_glob)

   end subroutine cycdump

   subroutine readc_1d_i4(lures,var)

      integer, intent(in) :: lures
      integer, pointer, intent(inout) :: var(:)

      integer buf(numcolumn_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_integer,0,p_comm,p_err)
#endif

      var = buf(ccmap)

   end subroutine readc_1d_i4

   subroutine readc_1d_r8(lures,var)

      integer, intent(in) :: lures
      real(r8), pointer, intent(inout) :: var(:)

      real(r8) buf(numcolumn_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_real8,0,p_comm,p_err)
#endif

      write(6,*) 'ccmap1f1', p_iam, minval(ccmap), maxval(ccmap)

      if(associated(var)) then
         write(6,*) 'ccmap1f2', p_iam, size(var), size(ccmap)
      else
         write(6,*) 'ccmap1f3', p_iam, 'pointer is not associated'
      end if

      var = buf(ccmap)

   end subroutine readc_1d_r8

   subroutine readc_2d_r8(lures,var,dim1sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      real(r8), pointer, intent(inout) :: var(:,:)

      real(r8) buf(dim1sz,numcolumn_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_real8,0,p_comm,p_err)
#endif

      var = buf(:,ccmap)

   end subroutine readc_2d_r8

   subroutine readc_3d_r8(lures,var,dim1sz,dim2sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      integer, intent(in) :: dim2sz
      real(r8), pointer, intent(inout) :: var(:,:,:)

      real(r8) buf(dim1sz,dim2sz,numcolumn_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_real8,0,p_comm,p_err)
#endif

      var = buf(:,:,ccmap)

   end subroutine readc_3d_r8

   subroutine readp_1d_i4(lures,var)

      integer, intent(in) :: lures
      integer, pointer, intent(inout) :: var(:)

      integer buf(numpatch_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_integer,0,p_comm,p_err)
#endif

      var = buf(ppmap)

   end subroutine readp_1d_i4

   subroutine readp_1d_r8(lures,var)

      integer, intent(in) :: lures
      real(r8), pointer, intent(inout) :: var(:)

      real(r8) buf(numpatch_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_real8,0,p_comm,p_err)
#endif

      var = buf(ppmap)

   end subroutine readp_1d_r8

   subroutine readp_2d_r8(lures,var,dim1sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      real(r8), pointer, intent(inout) :: var(:,:)

      real(r8) buf(dim1sz,numpatch_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_real8,0,p_comm,p_err)
#endif

      var = buf(:,ppmap)

   end subroutine readp_2d_r8

   subroutine readp_3d_r8(lures,var,dim1sz,dim2sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      integer, intent(in) :: dim2sz
      real(r8), pointer, intent(inout) :: var(:,:,:)

      real(r8) buf(dim1sz,dim2sz,numpatch_glob)

      if (p_master) read(lures) buf

#ifdef SPMD
      call mpi_bcast (buf,size(buf),mpi_real8,0,p_comm,p_err)
#endif

      var = buf(:,:,ppmap)

   end subroutine readp_3d_r8

   subroutine writec_1d_i4(lures,var)

      integer, intent(in) :: lures
      integer, pointer, intent(in) :: var(:)

      integer buf(numcolumn_glob)
      integer glob(numcolumn_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))
      end do

      call mpi_gatherv (var,size(var),mpi_integer,&
                        buf,p_counts,p_displs,mpi_integer,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(ccmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writec_1d_i4

   subroutine writec_1d_r8(lures,var)

      integer, intent(in) :: lures
      real(r8), pointer, intent(in) :: var(:)

      real(r8) buf(numcolumn_glob)
      real(r8) glob(numcolumn_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))
      end do

      call mpi_gatherv (var,size(var),mpi_real8,&
                        buf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(ccmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writec_1d_r8

   subroutine writec_2d_r8(lures,var,dim1sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      real(r8), pointer, intent(in) :: var(:,:)

      real(r8) buf(dim1sz,numcolumn_glob)
      real(r8) glob(dim1sz,numcolumn_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)*dim1sz
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))*dim1sz
      end do

      call mpi_gatherv (var,size(var),mpi_real8,&
                        buf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(:,ccmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writec_2d_r8

   subroutine writec_3d_r8(lures,var,dim1sz,dim2sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      integer, intent(in) :: dim2sz
      real(r8), pointer, intent(in) :: var(:,:,:)

      real(r8) buf(dim1sz,dim2sz,numcolumn_glob)
      real(r8) glob(dim1sz,dim2sz,numcolumn_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numcolumn_proc(:)*dim1sz*dim2sz
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numcolumn_proc(0:i-1))*dim1sz*dim2sz
      end do

      call mpi_gatherv (var,size(var),mpi_real8,&
                        buf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(:,:,ccmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writec_3d_r8

   subroutine writep_1d_i4(lures,var)

      integer, intent(in) :: lures
      integer, pointer, intent(in) :: var(:)

      integer buf(numpatch_glob)
      integer glob(numpatch_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numpatch_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numpatch_proc(0:i-1))
      end do

      call mpi_gatherv (var,size(var),mpi_integer,&
                        buf,p_counts,p_displs,mpi_integer,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(ppmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writep_1d_i4

   subroutine writep_1d_r8(lures,var)

      integer, intent(in) :: lures
      real(r8), pointer, intent(in) :: var(:)

      real(r8) buf(numpatch_glob)
      real(r8) glob(numpatch_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numpatch_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numpatch_proc(0:i-1))
      end do

      call mpi_gatherv (var,size(var),mpi_real8,&
                        buf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(ppmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writep_1d_r8

   subroutine writep_2d_r8(lures,var,dim1sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      real(r8), pointer, intent(in) :: var(:,:)

      real(r8) buf(dim1sz,numpatch_glob)
      real(r8) glob(dim1sz,numpatch_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numpatch_proc(:)*dim1sz
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numpatch_proc(0:i-1))*dim1sz
      end do

      call mpi_gatherv (var,size(var),mpi_real8,&
                        buf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(:,ppmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writep_2d_r8

   subroutine writep_3d_r8(lures,var,dim1sz,dim2sz)

      integer, intent(in) :: lures
      integer, intent(in) :: dim1sz
      integer, intent(in) :: dim2sz
      real(r8), pointer, intent(in) :: var(:,:,:)

      real(r8) buf(dim1sz,dim2sz,numpatch_glob)
      real(r8) glob(dim1sz,dim2sz,numpatch_glob)
      integer i

#ifdef SPMD
      p_counts(:) = numpatch_proc(:)*dim1sz*dim2sz
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numpatch_proc(0:i-1))*dim1sz*dim2sz
      end do

      call mpi_gatherv (var,size(var),mpi_real8,&
                        buf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      glob = var
#endif

      if (p_master) then
         glob(:,:,ppmap_glob) = buf
         write(lures) glob
      end if

   end subroutine writep_3d_r8

END MODULE colm_ioMod
