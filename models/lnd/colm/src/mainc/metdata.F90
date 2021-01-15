!*******************************************************************************
! Initial author: Duoying Ji (2014/08/24)
!*******************************************************************************
!
! Assume that time step of forcing data could be divided by time step of model evenly.
! Default model grids lay from south to north, west to east.
!
! GSWP2     : 89.50N~59.50S, 179.50W~179.50E
! CRUNCEP   : 89.75N~89.75S, 179.75W~179.75E
! Princeton : 89.50S~89.50N, 179.50W~179.50E
!
!*******************************************************************************

#include <define.h>

MODULE metdata

#if(defined GSWP2 || defined PRINCETON || defined FLUXNET || defined CRUNCEP)

   use netcdf
   use spmd
   use precision
   use colm_varctl
   use timemgr, only: get_curr_date, get_curr_year, get_dates_range

   implicit none

   save
   private

   real(r8), parameter :: spv = -9999.        ! special value
   integer,  parameter :: max_ntpr = 200      ! maximum number of time records to pre-read

 ! grids mapping method: 
 ! mapping_mode=1 for coarse or equal land grids, set mapping_maxg according to grid resolutions
 ! mapping_mode=2 for fine land grids, mapping_maxg = 4
   integer,  parameter :: mapping_mode = 1 
   integer,  parameter :: mapping_maxg = 16   ! maximum number of mapping grids
 ! integer,  parameter :: mapping_mode = 2 
 ! integer,  parameter :: mapping_maxg = 4    ! maximum number of mapping grids

   character(len=256) :: fname_qair  = "null"
   character(len=256) :: fname_tair  = "null"
   character(len=256) :: fname_dlwrf = "null"
   character(len=256) :: fname_dswrf = "null"
   character(len=256) :: fname_pres  = "null"
   character(len=256) :: fname_wind  = "null"
   character(len=256) :: fname_windu = "null"
   character(len=256) :: fname_windv = "null"
   character(len=256) :: fname_snow  = "null"
   character(len=256) :: fname_rain  = "null"
   character(len=256) :: fname_rainc = "null"
   character(len=256) :: fname_co2   = "null"

   character(len=256) :: vname_qair  = "null"
   character(len=256) :: vname_tair  = "null"
   character(len=256) :: vname_dlwrf = "null"
   character(len=256) :: vname_dswrf = "null"
   character(len=256) :: vname_pres  = "null"
   character(len=256) :: vname_wind  = "null"
   character(len=256) :: vname_windu = "null"
   character(len=256) :: vname_windv = "null"
   character(len=256) :: vname_snow  = "null"
   character(len=256) :: vname_rain  = "null"
   character(len=256) :: vname_rainc = "null"
   character(len=256) :: vname_co2   = "null"

   character(len=256) :: metpath

   integer :: fid_qair  = -1
   integer :: fid_tair  = -1
   integer :: fid_dlwrf = -1
   integer :: fid_dswrf = -1
   integer :: fid_pres  = -1
   integer :: fid_wind  = -1
   integer :: fid_windu = -1
   integer :: fid_windv = -1
   integer :: fid_snow  = -1
   integer :: fid_rain  = -1
   integer :: fid_rainc = -1
   integer :: fid_co2   = -1

   integer :: vid_qair  = -1
   integer :: vid_tair  = -1
   integer :: vid_dlwrf = -1
   integer :: vid_dswrf = -1
   integer :: vid_pres  = -1
   integer :: vid_wind  = -1
   integer :: vid_windu = -1
   integer :: vid_windv = -1
   integer :: vid_snow  = -1
   integer :: vid_rain  = -1
   integer :: vid_rainc = -1
   integer :: vid_co2   = -1

   integer itime_pr  ! index of first record of pre-read data at time coordinate
   integer itpr      ! index of current record in pre-read data
   integer ntpr      ! number of records currently pre-read

!* Current forcing data calculated from pre-read data

   real(r4), pointer :: qair (:)
   real(r4), pointer :: tair (:)
   real(r4), pointer :: dlwrf(:)
   real(r4), pointer :: dswrf(:)
   real(r4), pointer :: pres (:)
   real(r4), pointer :: wind (:)
   real(r4), pointer :: windu(:)
   real(r4), pointer :: windv(:)
   real(r4), pointer :: snow (:)
   real(r4), pointer :: rain (:)
   real(r4), pointer :: rainc(:)
   real(r4), pointer :: co2  (:)

!* Pointers to maximum records of pre-read data buffer

   real(r4), pointer :: pr_qair (:,:)
   real(r4), pointer :: pr_tair (:,:)
   real(r4), pointer :: pr_dlwrf(:,:)
   real(r4), pointer :: pr_dswrf(:,:)
   real(r4), pointer :: pr_pres (:,:)
   real(r4), pointer :: pr_wind (:,:)
   real(r4), pointer :: pr_windu(:,:)
   real(r4), pointer :: pr_windv(:,:)
   real(r4), pointer :: pr_snow (:,:)
   real(r4), pointer :: pr_rain (:,:)
   real(r4), pointer :: pr_rainc(:,:)
   real(r4), pointer :: pr_co2  (:,:)

   real(r4), pointer :: lon    (:)
   real(r4), pointer :: lat    (:)
   real(r4), pointer :: nav_lon(:,:)
   real(r4), pointer :: nav_lat(:,:)
   real(r4), pointer :: time   (:)
   integer , pointer :: land   (:)

   integer,  pointer :: gridnum(:,:)
   integer,  pointer :: gridmap(:,:,:)
   real(r4), pointer :: gridwgt(:,:,:)

!IF enable_forcing_anomaly == .true.
   integer, parameter :: anomaly_begyr = 2006
   integer, parameter :: anomaly_endyr = 2300

   character(len=255) fscale_dswrf
   character(len=255) fscale_dlwrf
   character(len=255) fscale_prcp
   character(len=255) fanomaly_tair
   character(len=255) fanomaly_qair
   character(len=255) fanomaly_pres
   character(len=255) fanomaly_windu
   character(len=255) fanomaly_windv

   integer :: fid_scale_dswrf
   integer :: fid_scale_dlwrf
   integer :: fid_scale_prcp
   integer :: fid_anomaly_tair
   integer :: fid_anomaly_qair
   integer :: fid_anomaly_pres
   integer :: fid_anomaly_windu
   integer :: fid_anomaly_windv

   integer :: vid_scale_dswrf
   integer :: vid_scale_dlwrf
   integer :: vid_scale_prcp
   integer :: vid_anomaly_tair
   integer :: vid_anomaly_qair
   integer :: vid_anomaly_pres
   integer :: vid_anomaly_windu
   integer :: vid_anomaly_windv

   real(r4), allocatable :: dswrf_scale  (:,:)
   real(r4), allocatable :: dlwrf_scale  (:,:)
   real(r4), allocatable ::  prcp_scale  (:,:)
   real(r4), allocatable ::  tair_anomaly(:,:)
   real(r4), allocatable ::  qair_anomaly(:,:)
   real(r4), allocatable ::  pres_anomaly(:,:)
   real(r4), allocatable :: windu_anomaly(:,:)
   real(r4), allocatable :: windv_anomaly(:,:)

   integer :: itime_anomaly = -1

   integer :: nlon_anomaly
   integer :: nlat_anomaly

   integer,  allocatable :: xmap_anomaly(:,:)     ! x-grid mapping from land model to anomaly dataset
   integer,  allocatable :: ymap_anomaly(:,:)     ! y-grid mapping from land model to anomaly dataset
!ENDIF enable_forcing_anomaly == .true.

   integer , pointer :: forcmask(:,:)         ! mask land grids according to forcing data

   integer :: origin_year, origin_month, origin_day, origin_second
   integer :: nlon, nlat, nland, ntime        ! length of dimensions of forcing data
   integer :: itime        ! current time index on time coordinate, itime in [1~ntime]
   integer :: timestp      ! seconds of time step of forcing data

   integer :: za_t = 50    ! reference height of temperature [m]
   integer :: za_q = 50    ! reference height of humidity [m]
   integer :: za_u = 50    ! reference height of wind [m]

   logical :: annCO2_flag = .false.      ! annual global mean [CO2] dataset 
   integer :: annCO2_begyr
   integer :: annCO2_endyr
   real(r8), pointer :: annCO2(:)

   interface metdata_init
      module procedure metdata_init
   end interface

   interface metdata_read
      module procedure metdata_read
   end interface

   interface metdata_close
      module procedure metdata_close
   end interface

   public metdata_init, metdata_read, metdata_close, forcmask

CONTAINS

   SUBROUTINE gen_fnames(year,month,day,second)

      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day
      integer, intent(in) :: second

      character(len=256)  :: date

      integer fyear, fmonth, fday

      if (enable_forcing_looping) then
         forcing_looping_yr = forcing_looping_yr + 1

         if(forcing_looping_yr < forcing_looping_begyr .or. &
            forcing_looping_yr > forcing_looping_endyr) then
            forcing_looping_yr = forcing_looping_begyr
         end if
   
         fyear  = forcing_looping_yr
         fmonth = month
         fday   = day
      else
         fyear  = year
         fmonth = month
         fday   = day

       !*Fix for CRUNCEP like dataset
         if (month == 12 .and. day == 31 .and. (86400-second) < timestp) fyear = year+1
      end if

#if (defined GSWP2)

      write(date,"(I4.4,I2.2)") fyear,fmonth

      fname_qair  = trim(metpath)//"Qair_cru/Qair_cru"        //trim(date)//"_halfhour.nc"
      fname_tair  = trim(metpath)//"Tair_cru/Tair_cru"        //trim(date)//"_halfhour.nc"
      fname_dlwrf = trim(metpath)//"LWdown_srb/LWdown_srb"    //trim(date)//"_halfhour.nc"
      fname_dswrf = trim(metpath)//"SWdown_srb/SWdown_srb"    //trim(date)//"_halfhour.nc"
      fname_pres  = trim(metpath)//"PSurf_ecor/PSurf_ecor"    //trim(date)//"_halfhour.nc"
      fname_wind  = trim(metpath)//"Wind_ncep/Wind_ncep"      //trim(date)//"_halfhour.nc"
      fname_windu = "null"
      fname_windv = "null"
      fname_snow  = trim(metpath)//"Snowf_gswp/Snowf_gswp"    //trim(date)//"_halfhour.nc"
      fname_rain  = trim(metpath)//"Rainf_gswp/Rainf_gswp"    //trim(date)//"_halfhour.nc"
      fname_rainc = trim(metpath)//"Rainf_C_gswp/Rainf_C_gswp"//trim(date)//"_halfhour.nc"
      fname_co2   = "null"

      vname_qair  = "Qair"
      vname_tair  = "Tair"
      vname_dlwrf = "LWdown"
      vname_dswrf = "SWdown"
      vname_pres  = "PSurf"
      vname_wind  = "Wind"
      vname_windu = "null"
      vname_windv = "null"
      vname_snow  = "Snowf"
      vname_rain  = "Rainf"
      vname_rainc = "Rainf_C"
      vname_co2   = "null"

#elif (defined CRUNCEP)

      write(date,"(I4.4)") fyear

      fname_qair  = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_tair  = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_dlwrf = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_dswrf = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_pres  = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_wind  = "null"
      fname_windu = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_windv = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_snow  = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_rain  = trim(metpath)//"cruncep_"//trim(date)//".nc"
      fname_rainc = "null"
      fname_co2   = "null"

      vname_qair  = "Qair"
      vname_tair  = "Tair"
      vname_dlwrf = "LWdown"
      vname_dswrf = "SWdown"
      vname_pres  = "PSurf"
      vname_wind  = "null"
      vname_windu = "Wind_E"
      vname_windv = "Wind_N"
      vname_snow  = "Snowf"
      vname_rain  = "Rainf"
      vname_rainc = "null"
      vname_co2   = "null"

#elif (defined PRINCETON)

      write(date,"(I4.4)") fyear

      fname_qair  = trim(metpath)//"shum/shum_halfhour_"  //trim(date)//"-"//trim(date)//".nc"
      fname_tair  = trim(metpath)//"tas/tas_halfhour_"    //trim(date)//"-"//trim(date)//".nc"
      fname_dlwrf = trim(metpath)//"dlwrf/dlwrf_halfhour_"//trim(date)//"-"//trim(date)//".nc"
      fname_dswrf = trim(metpath)//"dswrf/dswrf_halfhour_"//trim(date)//"-"//trim(date)//".nc"
      fname_pres  = trim(metpath)//"pres/pres_halfhour_"  //trim(date)//"-"//trim(date)//".nc"
      fname_wind  = trim(metpath)//"wind/wind_halfhour_"  //trim(date)//"-"//trim(date)//".nc"
      fname_windu = "null"
      fname_windv = "null"
      fname_rain  = trim(metpath)//"prcp/prcp_halfhour_"  //trim(date)//"-"//trim(date)//".nc"
      fname_rainc = "null"
      fname_snow  = "null"
      fname_co2   = "null"

      vname_qair  = "shum"
      vname_tair  = "tas"
      vname_dlwrf = "dlwrf"
      vname_dswrf = "dswrf"
      vname_pres  = "pres"
      vname_wind  = "wind"
      vname_windu = "null"
      vname_windv = "null"
      vname_rain  = "prcp"
      vname_rainc = "null"
      vname_snow  = "null"
      vname_co2   = "null"

#elif (defined FLUXNET)

      fname_qair  = trim(metpath)
      fname_tair  = trim(metpath)
      fname_dlwrf = trim(metpath)
      fname_dswrf = trim(metpath)
      fname_pres  = trim(metpath)
      fname_wind  = trim(metpath)
      fname_windu = "null"
      fname_windv = "null"
      fname_snow  = "null" !trim(metpath)
      fname_rain  = trim(metpath)
      fname_rainc = "null"
      fname_co2   = "null" !trim(metpath)

      vname_qair  = "Qair"
      vname_tair  = "Tair"
      vname_dlwrf = "LWdown"
      vname_dswrf = "SWdown"
      vname_pres  = "PSurf"
      vname_wind  = "Wind"
      vname_windu = "null"
      vname_windv = "null"
      vname_snow  = "null" !"Snowf"
      vname_rain  = "Rainf"
      vname_rainc = "null"
      vname_co2   = "null" !"CO2air"

#endif

   END SUBROUTINE gen_fnames

   SUBROUTINE annCO2_init(fco2)

      character(len=255), intent(in) :: fco2 ! global mean annual CO2 concentration

      integer, pointer :: year(:)
      integer fid, vid, nyr

      if(fco2/="") then
         annCO2_flag = .true.
      else
         annCO2_flag = .false.
      end if

#if(defined SPMD)
      call mpi_bcast(annCO2_flag,1,mpi_logical,0,p_comm,p_err)
#endif

      if(.not. annCO2_flag) return

      if (p_master) then
         call sanity(nf90_open(path=trim(fco2),mode=nf90_nowrite,ncid=fid))

         call sanity(nf90_inq_dimid(fid,'year',vid))
         call sanity(nf90_inquire_dimension(fid,vid,len=nyr))

         allocate(year(nyr))

         call sanity(nf90_inq_varid(fid,'year',vid))
         call sanity(nf90_get_var(fid,vid,year))

         annCO2_begyr = year(1)
         annCO2_endyr = year(nyr)

         write(6,*) 'Annual CO2:', annCO2_begyr, annCO2_endyr

         deallocate(year)
      end if

#if(defined SPMD)
      call mpi_bcast(annCO2_begyr,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(annCO2_endyr,1,mpi_integer,0,p_comm,p_err)
#endif

      allocate(annCO2(annCO2_endyr-annCO2_begyr+1))

      if (p_master) then
         call sanity(nf90_inq_varid(fid,'CO2air',vid))
         call sanity(nf90_get_var(fid,vid,annCO2))     ! units of ppmv

         annCO2 = annCO2*1.E-6

         call sanity(nf90_close(fid))
      end if

#if(defined SPMD)
      call mpi_bcast(annCO2,size(annCO2),mpi_real8,0,p_comm,p_err)
#endif

   END SUBROUTINE annCO2_init

   SUBROUTINE annCO2_update

      use timemgr, only: idate_p, idate
      integer curryr

      if (annCO2_flag) then
         if (enable_forcing_looping) then
            curryr = forcing_looping_yr
         else
            curryr = get_curr_year()
         end if

         if (curryr<annCO2_begyr) then
            pco2 = annCO2(1)
         else if(curryr>annCO2_endyr) then
            pco2 = annCO2(annCO2_endyr-annCO2_begyr+1)
         else
            pco2 = annCO2(curryr-annCO2_begyr+1)
         end if

         if(p_master .and. idate_p(1).ne.idate(1)) &
            write(6,*) 'Annual mean CO2 concentration:', pco2
      end if

   END SUBROUTINE annCO2_update

   SUBROUTINE annCO2_exit

      if (annCO2_flag) then
         deallocate (annCO2)
      end if

   END SUBROUTINE annCO2_exit

   SUBROUTINE metdata_init(lon_points,lat_points,fmet,fco2)

      integer, intent(in) :: lon_points      ! longitude points of land model
      integer, intent(in) :: lat_points      ! latitude points of land model
      character(len=255), intent(in) :: fmet ! file name or path to meteorological data
      character(len=255), intent(in) :: fco2 ! global mean annual CO2 concentration

      itime    = 0
      ntime    = 0

      itime_pr = 0
      itpr     = 0
      ntpr     = 0

      nland    = 0

#if(defined GSWP2 || defined PRINCETON || defined FLUXNET)
      timestp = 1800
#endif

#if(defined CRUNCEP)
      timestp = 21600
#endif

#if(defined GSWP2 || defined PRINCETON || defined CRUNCEP)
      metpath = trim(adjustl(fmet))//"/"
#else
      metpath = trim(adjustl(fmet))
#endif

      allocate(gridnum(lon_points,lat_points))
      allocate(gridmap(mapping_maxg,lon_points,lat_points))
      allocate(gridwgt(mapping_maxg,lon_points,lat_points))

#ifdef AUTOMASK
      allocate(forcmask(lon_points,lat_points))
      forcmask(:,:) = 0
#endif

      call annCO2_init(fco2) 

   END SUBROUTINE metdata_init

   SUBROUTINE metdata_read(lon_points,lat_points,gmask,lonw,lats,longxy,latixy,&
                           tair2d,qair2d,pres2d,rainc2d,rainl2d,&
                           windu2d,windv2d,dswrf2d,dlwrf2d,tair_z,qair_z,wind_z)

      integer, intent(in)  :: lon_points
      integer, intent(in)  :: lat_points
      integer, pointer :: gmask(:,:)
      real(r8),pointer :: lonw (:)
      real(r8),pointer :: lats (:)
      real(r8),pointer :: longxy (:,:)
      real(r8),pointer :: latixy (:,:)
      real(r8),pointer :: tair2d (:,:)
      real(r8),pointer :: qair2d (:,:)
      real(r8),pointer :: pres2d (:,:)
      real(r8),pointer :: rainc2d(:,:)
      real(r8),pointer :: rainl2d(:,:)
      real(r8),pointer :: windu2d(:,:)
      real(r8),pointer :: windv2d(:,:)
      real(r8),pointer :: dswrf2d(:,:)
      real(r8),pointer :: dlwrf2d(:,:)
      real(r8),pointer :: tair_z (:,:)
      real(r8),pointer :: qair_z (:,:)
      real(r8),pointer :: wind_z (:,:)

      real(r8) tmpbuf1(mapping_maxg), tmpbuf2(mapping_maxg)

      integer  i, j, k

   !* -----------------------------------------------

      if (enable_forcing_anomaly) then
         if(itime_anomaly < 0) then
            call open_anomaly_files(lon_points,lat_points,longxy,latixy,gmask)
         end if

         call read_anomaly_files(lon_points,lat_points,gmask)
      end if

      if(itime < 1) then
         CALL open_files
         CALL grid_mapping(lon_points,lat_points,gmask,lonw,lats,longxy,latixy)
      endif

      CALL read_files

      CALL annCO2_update

      tair2d (1:lon_points,1:lat_points) = spv
      qair2d (1:lon_points,1:lat_points) = spv
      pres2d (1:lon_points,1:lat_points) = spv
      rainc2d(1:lon_points,1:lat_points) = spv
      rainl2d(1:lon_points,1:lat_points) = spv
      windu2d(1:lon_points,1:lat_points) = spv
      windv2d(1:lon_points,1:lat_points) = spv
      dswrf2d(1:lon_points,1:lat_points) = spv
      dlwrf2d(1:lon_points,1:lat_points) = spv
      tair_z (1:lon_points,1:lat_points) = spv
      qair_z (1:lon_points,1:lat_points) = spv
      wind_z (1:lon_points,1:lat_points) = spv

      do j = 1, lat_points
      do i = 1, lon_points
         k = gridnum(i,j)

         if (k == 0) cycle

         if (gmask(i,j) /= p_iam) stop 'gmask conflicts with gridmap'

         if (fid_rainc > 0) then
            tmpbuf1(1:k) = rain(gridmap(1:k,i,j))*gridwgt(1:k,i,j)
            tmpbuf2(1:k) = rainc(gridmap(1:k,i,j))*gridwgt(1:k,i,j)

            rainc2d(i,j) = sum(tmpbuf2(1:k))
            rainl2d(i,j) = max(0.,sum(tmpbuf1(1:k)-tmpbuf2(1:k)))
         else
            tmpbuf1(1:k) = rain(gridmap(1:k,i,j))*gridwgt(1:k,i,j)
            rainc2d(i,j) = 0.
            rainl2d(i,j) = sum(tmpbuf1(1:k))
         endif

         if (rainl2d(i,j) < 0) &
            write(6,*), 'rainl2d(i,j).lt.0', gridmap(1:k,i,j),gridwgt(1:k,i,j),rain(gridmap(1:k,i,j))

         if (fid_snow > 0) then
            tmpbuf1(1:k) = snow(gridmap(1:k,i,j))*gridwgt(1:k,i,j)
            rainl2d(i,j) = rainl2d(i,j) + sum(tmpbuf1(1:k))
         endif

         tair2d (i,j) = sum(tair(gridmap(1:k,i,j))*gridwgt(1:k,i,j))
         qair2d (i,j) = sum(qair(gridmap(1:k,i,j))*gridwgt(1:k,i,j))
         pres2d (i,j) = sum(pres(gridmap(1:k,i,j))*gridwgt(1:k,i,j))

         if (fid_wind > 0) then
            windu2d(i,j) = sum(wind(gridmap(1:k,i,j))*gridwgt(1:k,i,j))
            windv2d(i,j) = 0.
         end if

         if (fid_windu > 0 .and. fid_windv > 0) then
            windu2d(i,j) = sum(windu(gridmap(1:k,i,j))*gridwgt(1:k,i,j))
            windv2d(i,j) = sum(windv(gridmap(1:k,i,j))*gridwgt(1:k,i,j))
         end if

         dswrf2d(i,j) = sum(dswrf(gridmap(1:k,i,j))*gridwgt(1:k,i,j))
         dlwrf2d(i,j) = sum(dlwrf(gridmap(1:k,i,j))*gridwgt(1:k,i,j))

         if (enable_forcing_anomaly) then
            dswrf2d(i,j) = dswrf2d(i,j) * dswrf_scale  (i,j)
            dlwrf2d(i,j) = dlwrf2d(i,j) * dlwrf_scale  (i,j)
            rainc2d(i,j) = rainc2d(i,j) *  prcp_scale  (i,j)
            rainl2d(i,j) = rainl2d(i,j) *  prcp_scale  (i,j)
   
             tair2d(i,j) =  tair2d(i,j) +  tair_anomaly(i,j)
             qair2d(i,j) =  qair2d(i,j) +  qair_anomaly(i,j)
             pres2d(i,j) =  pres2d(i,j) +  pres_anomaly(i,j)
            windu2d(i,j) = windu2d(i,j) + windu_anomaly(i,j)
            windv2d(i,j) = windv2d(i,j) + windv_anomaly(i,j)

            if (qair2d(i,j) < 0.) qair2d(i,j) = 0.
         end if

         tair_z (i,j) = za_t
         qair_z (i,j) = za_q
         wind_z (i,j) = za_u

       ! further tuning here
         if (fid_co2 > 0) then
            pco2 = sum(co2(gridmap(1:k,i,j))*gridwgt(1:k,i,j))   ! units of ppmv
            pco2 = pco2*1.E-6
         endif
      end do 
      end do 

      if (itime == ntime) then
         if(itpr /= ntpr) then
            write(6,*), 'Error in pre-read data indexing', itpr, ntpr, itime, ntime
            call abort
         end if

         itime    = 0
         ntime    = 0

         itime_pr = 0
         itpr     = 0
         ntpr     = 0

         CALL close_files
      end if

   END SUBROUTINE metdata_read

   SUBROUTINE metdata_close

      call close_files

      call var_dealloc

      deallocate (gridnum)
      deallocate (gridmap)
      deallocate (gridwgt)

      call annCO2_exit

#ifdef AUTOMASK
      deallocate (forcmask)
#endif

      if (enable_forcing_anomaly) then
         call close_anomaly_files
      end if

   END SUBROUTINE metdata_close

   SUBROUTINE open_files

      integer year, month, day, second
      integer vid_t, ret, nval, HH, MM, SS

      character(len=255) time_units

      call get_curr_date(year,month,day,second)

      if(p_master) then
         call gen_fnames(year,month,day,second)

         call ncdf_open(fname_qair ,vname_qair ,fid_qair ,vid_qair )
         call ncdf_open(fname_tair ,vname_tair ,fid_tair ,vid_tair )
         call ncdf_open(fname_dlwrf,vname_dlwrf,fid_dlwrf,vid_dlwrf)
         call ncdf_open(fname_dswrf,vname_dswrf,fid_dswrf,vid_dswrf)
         call ncdf_open(fname_pres ,vname_pres ,fid_pres ,vid_pres )
         call ncdf_open(fname_wind ,vname_wind ,fid_wind ,vid_wind )
         call ncdf_open(fname_windu,vname_windu,fid_windu,vid_windu)
         call ncdf_open(fname_windv,vname_windv,fid_windv,vid_windv)
         call ncdf_open(fname_snow ,vname_snow ,fid_snow ,vid_snow )
         call ncdf_open(fname_rain ,vname_rain ,fid_rain ,vid_rain )
         call ncdf_open(fname_rainc,vname_rainc,fid_rainc,vid_rainc)
         call ncdf_open(fname_co2  ,vname_co2  ,fid_co2  ,vid_co2  )
      end if

#if(defined SPMD)
      call mpi_bcast(fid_qair ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_tair ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_dlwrf,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_dswrf,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_pres ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_wind ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_windu,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_windv,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_snow ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_rain ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_rainc,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(fid_co2  ,1,mpi_integer,0,p_comm,p_err)

      call mpi_bcast(nlon ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(nlat ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(nland,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(ntime,1,mpi_integer,0,p_comm,p_err)
#endif

   !* allocate memory for pre-read and forcing data grids lon&lat
      call var_alloc

   !* reference height for T,Q,U
      if(p_master) then
         ret = nf90_inq_varid(fid_tair,'za',vid_t)
         if (ret==nf90_NoErr) then
            call sanity(nf90_get_var(fid_tair,vid_t,za_t))
         end if

         ret = nf90_inq_varid(fid_qair,'za',vid_t)
         if (ret==nf90_NoErr) then
            call sanity(nf90_get_var(fid_qair,vid_t,za_q))
         end if

         ret = nf90_inq_varid(fid_wind,'za',vid_t)
         if (ret==nf90_NoErr) then
            call sanity(nf90_get_var(fid_wind,vid_t,za_u))
         end if
      end if

#if(defined SPMD)
      call mpi_bcast(za_t,1,mpi_real,0,p_comm,p_err)
      call mpi_bcast(za_q,1,mpi_real,0,p_comm,p_err)
      call mpi_bcast(za_u,1,mpi_real,0,p_comm,p_err)
#endif

   !* assume all files have same coordinates information
      if(p_master) then

         nval = 0

         ret = nf90_inq_varid(fid_tair,'origin_year',vid_t)
         if (ret==nf90_NoErr) then 
            nval = nval+1
            call sanity(nf90_get_var(fid_tair,vid_t,origin_year))
         end if

         ret = nf90_inq_varid(fid_tair,'origin_month',vid_t)
         if (ret==nf90_NoErr) then 
            nval = nval+1
            call sanity(nf90_get_var(fid_tair,vid_t,origin_month))
         end if

         ret = nf90_inq_varid(fid_tair,'origin_day',vid_t)
         if (ret==nf90_NoErr) then 
            nval = nval+1
            call sanity(nf90_get_var(fid_tair,vid_t,origin_day))
         end if

         ret = nf90_inq_varid(fid_tair,'origin_second',vid_t)
         if (ret==nf90_NoErr) then 
            nval = nval+1
            call sanity(nf90_get_var(fid_tair,vid_t,origin_second))
         end if

         call sanity(nf90_inq_varid(fid_tair,'time',vid_t))
         call sanity(nf90_get_var(fid_tair,vid_t,time))

         if (nval /= 4) then
            ret = nf90_get_att(fid_tair, vid_t, "units", time_units)
            if (ret==nf90_NoErr) then  !* units of format like "seconds since 2000-01-01 00:00:00"
               if (time_units( 1:7 ) == "seconds" .and. &
                   time_units(19:19) == "-" .and. time_units(22:22) == "-" .and. &
                   time_units(28:28) == ":" .and. time_units(31:31) == ":") then
                   read(time_units(15:18),"(I4.4)") origin_year 
                   read(time_units(20:21),"(I2.2)") origin_month 
                   read(time_units(23:24),"(I2.2)") origin_day 
                   read(time_units(26:27),"(I2.2)") HH 
                   read(time_units(29:30),"(I2.2)") MM
                   read(time_units(32:33),"(I2.2)") SS

                   origin_second = HH*3600 + MM*60 + SS
               else
                  STOP "NON-recognized time units information"
               end if
            else
               STOP 'NO units for time in forcing data'
            end if
         end if

         ret = nf90_inq_varid(fid_tair,'lon',vid_t)
         if (ret==nf90_NoErr) then 
            call sanity(nf90_get_var(fid_tair,vid_t,lon(:)))
         else
            call sanity(nf90_inq_varid(fid_tair,'nav_lon',vid_t))
            call sanity(nf90_get_var(fid_tair,vid_t,nav_lon(:,:)))
            lon(:) = nav_lon(:,1)
         end if

         ret = nf90_inq_varid(fid_tair,'lat',vid_t)
         if (ret==nf90_NoErr) then 
            call sanity(nf90_get_var(fid_tair,vid_t,lat(:)))
         else
            call sanity(nf90_inq_varid(fid_tair,'nav_lat',vid_t))
            call sanity(nf90_get_var(fid_tair,vid_t,nav_lat(:,:)))
            lat(:) = nav_lat(1,:)
         end if

         call sanity(nf90_inq_varid(fid_tair,'land',vid_t))
         call sanity(nf90_get_var(fid_tair,vid_t,land(:)))

         write (6,"(A, X, I4.4 ':' I2.2 ':' I2.2 ':' I6.6)"), &
                  'Initial date of forcing data: yy:mm:dd:ss', origin_year, origin_month, origin_day, origin_second
         write (6,"(A, X, I4.4 ':' I2.2 ':' I2.2 ':' I6.6)"), &
                  'Model date: yy:mm:dd:ss', year, month, day, second

      end if

#if(defined SPMD)
      call mpi_bcast(origin_year,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(origin_month,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(origin_day,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(origin_second,1,mpi_integer,0,p_comm,p_err)

      call mpi_bcast(lon,size(lon),mpi_real,0,p_comm,p_err)
      call mpi_bcast(lat,size(lat),mpi_real,0,p_comm,p_err)
      call mpi_bcast(land,size(land),mpi_integer,0,p_comm,p_err)
      call mpi_bcast(time,size(time),mpi_real,0,p_comm,p_err)
#endif

   END SUBROUTINE open_files

   SUBROUTINE read_files

      integer year, month, day, second, nseconds, delt, yearx
      integer tidx1, tidx2, itprx
      real(r4) frac

      call get_curr_date(year,month,day,second)

      if(enable_forcing_looping) then
       !*Fix for CRUNCEP like dataset
         if (month == 12 .and. day == 31 .and. (86400-second) < timestp) then
            yearx = origin_year-1
         else
            yearx = origin_year
         end if

         call get_dates_range(origin_year,origin_month,origin_day,origin_second, &
                              yearx,month,day,second,nseconds)
      else
         call get_dates_range(origin_year,origin_month,origin_day,origin_second, &
                              year,month,day,second,nseconds)
      end if

      delt = nseconds-time(1)

      if (delt >= 0) then
         itime = INT(delt/timestp) + 1
         frac  = MOD(delt,timestp)/FLOAT(timestp)
      else
         itime = 1
         frac  = MOD(-delt,timestp)/FLOAT(timestp)
      end if

      itpr = itime - itime_pr + 1

    ! read data when itpr==ntpr, no matter itpr(ntpr) equals zero or something else.
      if (itpr >= ntpr) then
         if (ntpr /= 0) then !* NOT READING A NEW FILE
            if(fid_qair .gt.0) pr_qair (:,0) = pr_qair (:,ntpr)
            if(fid_tair .gt.0) pr_tair (:,0) = pr_tair (:,ntpr)
            if(fid_dlwrf.gt.0) pr_dlwrf(:,0) = pr_dlwrf(:,ntpr)
            if(fid_dswrf.gt.0) pr_dswrf(:,0) = pr_dswrf(:,ntpr)
            if(fid_pres .gt.0) pr_pres (:,0) = pr_pres (:,ntpr)
            if(fid_wind .gt.0) pr_wind (:,0) = pr_wind (:,ntpr)
            if(fid_windu.gt.0) pr_windu(:,0) = pr_windu(:,ntpr)
            if(fid_windv.gt.0) pr_windv(:,0) = pr_windv(:,ntpr)
            if(fid_snow .gt.0) pr_snow (:,0) = pr_snow (:,ntpr)
            if(fid_rain .gt.0) pr_rain (:,0) = pr_rain (:,ntpr)
            if(fid_rainc.gt.0) pr_rainc(:,0) = pr_rainc(:,ntpr)
            if(fid_co2  .gt.0) pr_co2  (:,0) = pr_co2  (:,ntpr)

            itpr = 0
            itime_pr = itime + 1
         else
            itpr = 1
            itime_pr = itime
         end if

         ntpr  = min(ntime-itime_pr+1,max_ntpr)
         tidx1 = itime_pr
         tidx2 = itime_pr + ntpr - 1

         if (itime_pr <= ntime) then
            if (p_master) then
               call ncdf_read(fid_qair ,vid_qair ,tidx1,tidx2,pr_qair (:,1:))
               call ncdf_read(fid_tair ,vid_tair ,tidx1,tidx2,pr_tair (:,1:))
               call ncdf_read(fid_dlwrf,vid_dlwrf,tidx1,tidx2,pr_dlwrf(:,1:))
               call ncdf_read(fid_dswrf,vid_dswrf,tidx1,tidx2,pr_dswrf(:,1:))
               call ncdf_read(fid_pres ,vid_pres ,tidx1,tidx2,pr_pres (:,1:))
               call ncdf_read(fid_wind ,vid_wind ,tidx1,tidx2,pr_wind (:,1:))
               call ncdf_read(fid_windu,vid_windu,tidx1,tidx2,pr_windu(:,1:))
               call ncdf_read(fid_windv,vid_windv,tidx1,tidx2,pr_windv(:,1:))
               call ncdf_read(fid_snow ,vid_snow ,tidx1,tidx2,pr_snow (:,1:))
               call ncdf_read(fid_rain ,vid_rain ,tidx1,tidx2,pr_rain (:,1:))
               call ncdf_read(fid_rainc,vid_rainc,tidx1,tidx2,pr_rainc(:,1:))
               call ncdf_read(fid_co2  ,vid_co2  ,tidx1,tidx2,pr_co2  (:,1:))
            end if

#if(defined SPMD)
            if (fid_qair .gt.0) &
               call mpi_bcast(pr_qair ,size(pr_qair ),mpi_real,0,p_comm,p_err)
            if (fid_tair .gt.0) &
               call mpi_bcast(pr_tair ,size(pr_tair ),mpi_real,0,p_comm,p_err)
            if (fid_dlwrf.gt.0) &
               call mpi_bcast(pr_dlwrf,size(pr_dlwrf),mpi_real,0,p_comm,p_err)
            if (fid_dswrf.gt.0) &
               call mpi_bcast(pr_dswrf,size(pr_dswrf),mpi_real,0,p_comm,p_err)
            if (fid_pres .gt.0) &
               call mpi_bcast(pr_pres ,size(pr_pres ),mpi_real,0,p_comm,p_err)
            if (fid_wind .gt.0) &
               call mpi_bcast(pr_wind ,size(pr_wind ),mpi_real,0,p_comm,p_err)
            if (fid_windu.gt.0) &
               call mpi_bcast(pr_windu,size(pr_windu),mpi_real,0,p_comm,p_err)
            if (fid_windv.gt.0) &
               call mpi_bcast(pr_windv,size(pr_windv),mpi_real,0,p_comm,p_err)
            if (fid_snow .gt.0) &
               call mpi_bcast(pr_snow ,size(pr_snow ),mpi_real,0,p_comm,p_err)
            if (fid_rain .gt.0) &
               call mpi_bcast(pr_rain ,size(pr_rain ),mpi_real,0,p_comm,p_err)
            if (fid_rainc.gt.0) &
               call mpi_bcast(pr_rainc,size(pr_rainc),mpi_real,0,p_comm,p_err)
            if (fid_co2  .gt.0) &
               call mpi_bcast(pr_co2  ,size(pr_co2  ),mpi_real,0,p_comm,p_err)
#endif
         end if
      end if

      if(abs(frac) < 1.E-6) then  ! Do NO temporal interpolation
         itprx = itpr
      else if (delt > 0) then
         itprx = itpr + 1
         if (itprx > ntpr) then
            if(p_master) write(6,*) 'date', year, month, day, second
            if(p_master) write(6,*) 'itime, itpr, itprx, ntpr', itime, itpr, itprx, ntpr
            STOP 'itrpx > ntpr'
         end if
      else if (delt < 0) then
         itprx = itpr - 1
         if (itprx < 0) then
            if(p_master) write(6,*) 'date', year, month, day, second
            if(p_master) write(6,*) 'itime, itpr, itprx, ntpr', itime, itpr, itprx, ntpr
            STOP 'itprx <0'
         end if
      end if

      if(fid_qair .gt.0) qair (:) = pr_qair (:,itpr)*(1-frac) + pr_qair (:,itprx)*frac
      if(fid_tair .gt.0) tair (:) = pr_tair (:,itpr)*(1-frac) + pr_tair (:,itprx)*frac
      if(fid_dlwrf.gt.0) dlwrf(:) = pr_dlwrf(:,itpr)*(1-frac) + pr_dlwrf(:,itprx)*frac
      if(fid_dswrf.gt.0) dswrf(:) = pr_dswrf(:,itpr)*(1-frac) + pr_dswrf(:,itprx)*frac
      if(fid_pres .gt.0) pres (:) = pr_pres (:,itpr)*(1-frac) + pr_pres (:,itprx)*frac
      if(fid_wind .gt.0) wind (:) = pr_wind (:,itpr)*(1-frac) + pr_wind (:,itprx)*frac
      if(fid_windu.gt.0) windu(:) = pr_windu(:,itpr)*(1-frac) + pr_windu(:,itprx)*frac
      if(fid_windv.gt.0) windv(:) = pr_windv(:,itpr)*(1-frac) + pr_windv(:,itprx)*frac
      if(fid_snow .gt.0) snow (:) = pr_snow (:,itpr)
      if(fid_rain .gt.0) rain (:) = pr_rain (:,itpr)
      if(fid_rainc.gt.0) rainc(:) = pr_rainc(:,itpr)
      if(fid_co2  .gt.0) co2  (:) = pr_co2  (:,itpr)*(1-frac) + pr_co2  (:,itprx)*frac

   END SUBROUTINE read_files

   SUBROUTINE close_files
   
      if(p_master) then
         call ncdf_close(fid_qair ,vid_qair )
         call ncdf_close(fid_tair ,vid_tair )
         call ncdf_close(fid_dlwrf,vid_dlwrf)
         call ncdf_close(fid_dswrf,vid_dswrf)
         call ncdf_close(fid_pres ,vid_pres )
         call ncdf_close(fid_wind ,vid_wind )
         call ncdf_close(fid_windu,vid_windu)
         call ncdf_close(fid_windv,vid_windv)
         call ncdf_close(fid_snow ,vid_snow )
         call ncdf_close(fid_rain ,vid_rain )
         call ncdf_close(fid_rainc,vid_rainc)
         call ncdf_close(fid_co2  ,vid_co2  )
      end if

   END SUBROUTINE close_files

   SUBROUTINE ncdf_open(fname,vname,fid,vid)

      character(len=256), intent(in) :: fname
      character(len=256), intent(in) :: vname
      integer, intent(out) :: fid
      integer, intent(out) :: vid

      integer vid_t, nlon_t, nlat_t, nland_t, ret

   !* ---------------------------------

      if (fname /= "null") then
         write(6,*), 'Open file: '//trim(fname)

         call sanity(nf90_open(path=trim(fname),mode=nf90_nowrite,ncid=fid))
         call sanity(nf90_inq_varid(fid,trim(vname),vid))

         ret = nf90_inq_dimid(fid,'lon',vid_t)
         if (ret /= nf90_NoErr) then
            call sanity(nf90_inq_dimid(fid,'x',vid_t))
         end if
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nlon))

         ret = nf90_inq_dimid(fid,'lat',vid_t)
         if (ret /= nf90_NoErr) then
            call sanity(nf90_inq_dimid(fid,'y',vid_t))
         end if
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nlat))

         ret = nf90_inq_dimid(fid,'time',vid_t)
         if (ret /= nf90_NoErr) then
            call sanity(nf90_inq_dimid(fid,'tstep',vid_t))
         end if
         call sanity(nf90_inquire_dimension(fid,vid_t,len=ntime))

         call sanity(nf90_inq_dimid(fid,'land',vid_t))
         call sanity(nf90_inquire_dimension(fid,vid_t,len=nland))

         write(6,*), 'nlat :', nlat
         write(6,*), 'nlon :', nlon
         write(6,*), 'nland:', nland
         write(6,*), 'ntime:', ntime
      end if

   END SUBROUTINE ncdf_open

   SUBROUTINE ncdf_close(fid,vid)

      integer, intent(inout) :: fid
      integer, intent(inout) :: vid

   !* ---------------------------------

      if (fid > 0) then
         call sanity(nf90_close(fid))
         fid = -1
         vid = -1
      end if

   END SUBROUTINE ncdf_close

   SUBROUTINE ncdf_read(fid,vid,tidx1,tidx2,vbuf)

      integer,  intent(in) :: fid
      integer,  intent(in) :: vid
      integer,  intent(in) :: tidx1      ! begin time record to read
      integer,  intent(in) :: tidx2      ! end time record to read
      real(r4), intent(out):: vbuf(:,:)  ! pre-read buffer for time record between <tidx1:tidx2>

   !* ---------------------------------

      if (fid > 0) then
         call sanity(nf90_get_var(fid,vid,vbuf,start=(/1,tidx1/),count=(/nland,tidx2-tidx1+1/)))
      end if

   END SUBROUTINE ncdf_read

   SUBROUTINE grid_mapping(lon_points,lat_points,gmask,lonw,lats,longxy,latixy)
   !
   !--A simple grid mapping method --
   !
      integer,  intent(in) :: lon_points
      integer,  intent(in) :: lat_points
      integer,  pointer :: gmask(:,:)
      real(r8), pointer :: lonw(:)
      real(r8), pointer :: lats(:)
      real(r8), pointer :: longxy(:,:)
      real(r8), pointer :: latixy(:,:)

      integer,  allocatable :: xy2land(:,:)
      real(r8), allocatable :: longxyt(:,:)  ! transformation of longxy

      integer i, j, k, i1, j1, i2, j2

      integer x1, x2, y1, y2, ngrd

      real(r8) dx1, dx2, dy1, dy2, a, minlon, maxlon

      allocate(xy2land(nlon,nlat))

      xy2land(:,:) = 0

      gridnum(:,:) = 0
      gridmap(:,:,:) = 0
      gridwgt(:,:,:) = 0.

      do k = 1, nland
         i = mod(land(k),nlon)
         j = land(k)/nlon+1

         if(i.eq.0) then
            i = nlon
            j = j-1
         end if

         xy2land(i,j) = k  ! native 1d->2d mapping in forcing data
      end do

      IF(mapping_mode == 1) THEN

       ! Original longitude range of meteorology data is 180W~180E
       ! Transform forcing data grids according to land model grids
       ! Here one land grid contains several focing data grids.
         minlon = minval(lonw)
         maxlon = maxval(lonw)

         if(maxlon.gt.180.) then
            do i = 1, nlon
             ! if(lon(i).lt.minlon) lon(i) = lon(i)+360.
               if(lon(i).lt.minlon .and. lon(i).lt.0) lon(i) = lon(i)+360.
             ! if(lon(i).gt.maxlon) stop 'error in longitude transformation'
            end do
         end if

         do j1 = 1, lat_points
         do i1 = 1, lon_points
            if(gmask(i1,j1) /= p_iam) cycle

            ngrd = 0

            do j2 = 1, nlat
               if(lats(j1).lt.lat(j2) .and. lat(j2).lt.lats(j1+1)) then  !* land model grids in S->N order
                  do i2 = 1, nlon
                     if(lonw(i1).lt.lon(i2) .and. lon(i2).lt.lonw(i1+1)) then
                        if(xy2land(i2,j2).gt.0) then
                           ngrd = ngrd+1
                           gridmap(ngrd,i1,j1) = xy2land(i2,j2)
                        end if
                     end if
                  end do
               end if
            end do

            if(ngrd.eq.0)then
#ifdef AUTOMASK
               forcmask(i1,j1) = 1
#else
               write(6,*) "Fatal error in mapping(mode=1, ngrd=0):",lonw(i1),lonw(i1+1),lats(j1),lats(j1+1)
               call abort
#endif
            else if(ngrd.gt.mapping_maxg) then
               write(6,*) "Fatal error in mapping(mode=1, ngrd>mapping_maxg)",ngrd,mapping_maxg
               call abort
            end if

            if(ngrd.gt.0) then
               gridnum(i1,j1) = ngrd
               gridwgt(1:ngrd,i1,j1) = 1.0/ngrd
            end if

         end do
         end do

      ELSE IF(mapping_mode == 2) THEN

       ! Original longitude range of meteorology data is 180W~180E
       ! Transform land model grids according to forcing data grids
       ! Here one forcing data grid contains several land model grids

         allocate (longxyt(lon_points,lat_points))

       ! begin of transform
         minlon = minval(lonw)
         maxlon = maxval(lonw)

         if(maxlon.gt.180.) then
            do j = 1, lat_points
            do i = 1, lon_points
               if (longxy(i,j).gt.180.) then
                  longxyt(i,j) = longxy(i,j) - 360.
               else
                  longxyt(i,j) = longxy(i,j)
               end if
            end do
            end do
         else
            longxyt = longxy
         end if
       ! end of transform

         do j1 = 1, lat_points
         do i1 = 1, lon_points
            if(gmask(i1,j1) /= p_iam) cycle

            x1 = 0
            x2 = 0
            y1 = 0
            y2 = 0
   
            if(longxyt(i1,j1).lt.lon(1)) then
               x1 = 1 
               x2 = 1 
            else if(longxyt(i1,j1).ge.lon(nlon)) then
               x1 = nlon
               x2 = nlon
            else ! inner points
               do i2 = 1, nlon-1
                  if(longxyt(i1,j1).ge.lon(i2) .and. longxyt(i1,j1).lt.lon(i2+1)) then
                     x1 = i2
                     x2 = i2+1
                     exit
                  end if
               end do
            end if
   
            if(latixy(i1,j1).lt.lat(1)) then
               y1 = 1 
               y2 = 1 
            else if(latixy(i1,j1).ge.lat(nlat)) then
               y1 = nlat
               y2 = nlat
            else ! inner points
               do j2 = 1, nlat-1
                  if(latixy(i1,j1).ge.lat(j2) .and. latixy(i1,j1).lt.lat(j2+1)) then
                     y1 = j2
                     y2 = j2+1
                     exit
                  end if
               end do
            endif

            if(x1.gt.0 .and. y1.gt.0) then
               dx1 = longxyt(i1,j1)-lon(x1)
               dx2 = lon(x2)-longxyt(i1,j1)
               dy1 = latixy(i1,j1)-lat(y1)
               dy2 = lat(y2)-latixy(i1,j1)
   
               ngrd = 0
               a = (lon(x2)-lon(x1))*(lat(y2)-lat(y1))

               if (a.gt.1.e-6) then ! inner points
                  if(xy2land(x1,y1).gt.0) then
                     ngrd = ngrd+1
                     gridmap(ngrd,i1,j1) = xy2land(x1,y1)
                     gridwgt(ngrd,i1,j1) =  dx2*dy2/a
                  end if
      
                  if(xy2land(x2,y1).gt.0) then
                     ngrd = ngrd+1
                     gridmap(ngrd,i1,j1) = xy2land(x2,y1)
                     gridwgt(ngrd,i1,j1) = dx1*dy2/a
                  end if
      
                  if(xy2land(x2,y2).gt.0) then
                     ngrd = ngrd+1
                     gridmap(ngrd,i1,j1) = xy2land(x2,y2)
                     gridwgt(ngrd,i1,j1) = dx1*dy1/a
                  end if
      
                  if(xy2land(x1,y2).gt.0) then
                     ngrd = ngrd+1
                     gridmap(ngrd,i1,j1) = xy2land(x1,y2)
                     gridwgt(ngrd,i1,j1) = dx2*dy1/a
                  end if

                  a = sum(gridwgt(1:ngrd,i1,j1))

                  if(a.le.0 .or. a.gt.(1.0+1.e-6)) then
#ifdef AUTOMASK
                     forcmask(i1,j1) = 1
#else
                     write(6,*) 'None points near the land grid(a):', longxy(i1,j1), latixy(i1,j1)
                     call abort
#endif
                  end if

                  gridnum(i1,j1) = ngrd
                  gridwgt(1:ngrd,i1,j1) = gridwgt(1:ngrd,i1,j1)/a
               else   ! boundary points
                  if(xy2land(x1,y1).gt.0) then
                     gridnum(i1,j1) = 1
                     gridmap(1,i1,j1) = xy2land(x1,y1)
                     gridwgt(1,i1,j1) = 1.0
                  else if(xy2land(x1,y2).gt.0) then
                     gridnum(i1,j1) = 1
                     gridmap(1,i1,j1) = xy2land(x1,y2)
                     gridwgt(1,i1,j1) = 1.0
                  else if(xy2land(x2,y1).gt.0) then
                     gridnum(i1,j1) = 1
                     gridmap(1,i1,j1) = xy2land(x2,y1)
                     gridwgt(1,i1,j1) = 1.0
                  else if(xy2land(x2,y2).gt.0) then
                     gridnum(i1,j1) = 1
                     gridmap(1,i1,j1) = xy2land(x2,y2)
                     gridwgt(1,i1,j1) = 1.0
                  else
#ifdef AUTOMASK
                     forcmask(i1,j1) = 1
#else
                     write(6,*) 'None points near the land grid(b):', longxy(i1,j1), latixy(i1,j1)
                     call abort
#endif
                  end if
               end if
            else
#ifdef AUTOMASK
               forcmask(i1,j1) = 1
#else
               write(6,*) 'None points near the land grid(c):', longxy(i1,j1), latixy(i1,j1)
               call abort
#endif
            end if
         end do
         end do

         deallocate (longxyt)

      END IF

#ifdef AUTOMASK
      if(sum(forcmask).gt.0) then
         write(6,*) 'Total grids masked by (', p_iam, '):', sum(forcmask)
      end if
#endif

      deallocate (xy2land)

   END SUBROUTINE grid_mapping

   SUBROUTINE var_alloc

      if(.not. associated(lon     )) allocate(lon    (nlon     ))
      if(.not. associated(lat     )) allocate(lat    (nlat     ))
      if(.not. associated(land    )) allocate(land   (nland    ))
      if(.not. associated(time    )) allocate(time   (ntime    ))
      if(.not. associated(nav_lon )) allocate(nav_lon(nlon,nlat))
      if(.not. associated(nav_lat )) allocate(nav_lat(nlon,nlat))

      if(.not. associated(qair    )) allocate(qair   (nland))
      if(.not. associated(tair    )) allocate(tair   (nland))
      if(.not. associated(dlwrf   )) allocate(dlwrf  (nland))
      if(.not. associated(dswrf   )) allocate(dswrf  (nland))
      if(.not. associated(pres    )) allocate(pres   (nland))
      if(.not. associated(wind    )) allocate(wind   (nland))
      if(.not. associated(windu   )) allocate(windu  (nland))
      if(.not. associated(windv   )) allocate(windv  (nland))
      if(.not. associated(snow    )) allocate(snow   (nland))
      if(.not. associated(rain    )) allocate(rain   (nland))
      if(.not. associated(rainc   )) allocate(rainc  (nland))
      if(.not. associated(co2     )) allocate(co2    (nland))

   !* Allocate memory for pre-read
   !* the 0th record for saving last record from previous pre-reading

      if(.not. associated(pr_qair )) allocate(pr_qair (nland,0:max_ntpr))
      if(.not. associated(pr_tair )) allocate(pr_tair (nland,0:max_ntpr))
      if(.not. associated(pr_dlwrf)) allocate(pr_dlwrf(nland,0:max_ntpr))
      if(.not. associated(pr_dswrf)) allocate(pr_dswrf(nland,0:max_ntpr))
      if(.not. associated(pr_pres )) allocate(pr_pres (nland,0:max_ntpr))
      if(.not. associated(pr_wind )) allocate(pr_wind (nland,0:max_ntpr))
      if(.not. associated(pr_windu)) allocate(pr_windu(nland,0:max_ntpr))
      if(.not. associated(pr_windv)) allocate(pr_windv(nland,0:max_ntpr))
      if(.not. associated(pr_snow )) allocate(pr_snow (nland,0:max_ntpr))
      if(.not. associated(pr_rain )) allocate(pr_rain (nland,0:max_ntpr))
      if(.not. associated(pr_rainc)) allocate(pr_rainc(nland,0:max_ntpr))
      if(.not. associated(pr_co2  )) allocate(pr_co2  (nland,0:max_ntpr))

   ENDSUBROUTINE var_alloc

   SUBROUTINE var_dealloc

      if (associated(lon     )) deallocate(lon     )
      if (associated(lat     )) deallocate(lat     )
      if (associated(land    )) deallocate(land    )
      if (associated(time    )) deallocate(time    )
      if (associated(nav_lon )) deallocate(nav_lon )
      if (associated(nav_lat )) deallocate(nav_lat )

      if (associated(qair    )) deallocate(qair    )
      if (associated(tair    )) deallocate(tair    )
      if (associated(dlwrf   )) deallocate(dlwrf   )
      if (associated(dswrf   )) deallocate(dswrf   )
      if (associated(pres    )) deallocate(pres    )
      if (associated(wind    )) deallocate(wind    )
      if (associated(windu   )) deallocate(windu   )
      if (associated(windv   )) deallocate(windv   )
      if (associated(snow    )) deallocate(snow    )
      if (associated(rain    )) deallocate(rain    )
      if (associated(rainc   )) deallocate(rainc   )
      if (associated(co2     )) deallocate(co2     )

      if (associated(pr_qair )) deallocate(pr_qair )
      if (associated(pr_tair )) deallocate(pr_tair )
      if (associated(pr_dlwrf)) deallocate(pr_dlwrf)
      if (associated(pr_dswrf)) deallocate(pr_dswrf)
      if (associated(pr_pres )) deallocate(pr_pres )
      if (associated(pr_wind )) deallocate(pr_wind )
      if (associated(pr_windu)) deallocate(pr_windu)
      if (associated(pr_windv)) deallocate(pr_windv)
      if (associated(pr_snow )) deallocate(pr_snow )
      if (associated(pr_rain )) deallocate(pr_rain )
      if (associated(pr_rainc)) deallocate(pr_rainc)
      if (associated(pr_co2  )) deallocate(pr_co2  )

   END SUBROUTINE var_dealloc

   SUBROUTINE open_anomaly_files(lon_points,lat_points,longxy,latixy,gmask)

      integer, intent(in) :: lon_points      ! longitude points of land model
      integer, intent(in) :: lat_points      ! latitude points of land model
      real(r8),intent(in) :: longxy(lon_points,lat_points)
      real(r8),intent(in) :: latixy(lon_points,lat_points)
      integer, intent(in) :: gmask (lon_points,lat_points)

      real(r4), allocatable :: lon_anomaly(:)
      real(r4), allocatable :: lat_anomaly(:)
      real(r4), allocatable :: lonw_anomaly(:)
      real(r4), allocatable :: lats_anomaly(:)

      real(r8), allocatable :: longxy_cycle(:,:)

      real(r8) minlon, maxlon

      integer vid_t, ix, jx, i, j

      fscale_dswrf   = "/mnt/gfs/jidy/Data/RCN-Forcing/rsds.ccsm4.rcp45.2006-2300.nc"
      fscale_dlwrf   = "/mnt/gfs/jidy/Data/RCN-Forcing/rlds.ccsm4.rcp45.2006-2300.nc"
      fscale_prcp    = "/mnt/gfs/jidy/Data/RCN-Forcing/pr.ccsm4.rcp45.2006-2300.nc"
      fanomaly_tair  = "/mnt/gfs/jidy/Data/RCN-Forcing/tas.ccsm4.rcp45.2006-2300.nc"
      fanomaly_qair  = "/mnt/gfs/jidy/Data/RCN-Forcing/huss.ccsm4.rcp45.2006-2300.nc"
      fanomaly_pres  = "/mnt/gfs/jidy/Data/RCN-Forcing/ps.ccsm4.rcp45.2006-2300.nc"
      fanomaly_pres  = "/mnt/gfs/jidy/Data/RCN-Forcing/ps.ccsm4.rcp45.2006-2300.nc"
      fanomaly_windu = "/mnt/gfs/jidy/Data/RCN-Forcing/uas.ccsm4.rcp45.2006-2300.nc"
      fanomaly_windv = "/mnt/gfs/jidy/Data/RCN-Forcing/vas.ccsm4.rcp45.2006-2300.nc"

      fscale_dswrf   = "/mnt/gfs/jidy/Data/RCN-Forcing/rsds.ccsm4.rcp85.2006-2300.nc"
      fscale_dlwrf   = "/mnt/gfs/jidy/Data/RCN-Forcing/rlds.ccsm4.rcp85.2006-2300.nc"
      fscale_prcp    = "/mnt/gfs/jidy/Data/RCN-Forcing/pr.ccsm4.rcp85.2006-2300.nc"
      fanomaly_tair  = "/mnt/gfs/jidy/Data/RCN-Forcing/tas.ccsm4.rcp85.2006-2300.nc"
      fanomaly_qair  = "/mnt/gfs/jidy/Data/RCN-Forcing/huss.ccsm4.rcp85.2006-2300.nc"
      fanomaly_pres  = "/mnt/gfs/jidy/Data/RCN-Forcing/ps.ccsm4.rcp85.2006-2300.nc"
      fanomaly_pres  = "/mnt/gfs/jidy/Data/RCN-Forcing/ps.ccsm4.rcp85.2006-2300.nc"
      fanomaly_windu = "/mnt/gfs/jidy/Data/RCN-Forcing/uas.ccsm4.rcp85.2006-2300.nc"
      fanomaly_windv = "/mnt/gfs/jidy/Data/RCN-Forcing/vas.ccsm4.rcp85.2006-2300.nc"

      allocate (dswrf_scale  (lon_points,lat_points))
      allocate (dlwrf_scale  (lon_points,lat_points))
      allocate ( prcp_scale  (lon_points,lat_points))
      allocate ( tair_anomaly(lon_points,lat_points))
      allocate ( qair_anomaly(lon_points,lat_points))
      allocate ( pres_anomaly(lon_points,lat_points))
      allocate (windu_anomaly(lon_points,lat_points))
      allocate (windv_anomaly(lon_points,lat_points))

      allocate ( xmap_anomaly(lon_points,lat_points))
      allocate ( ymap_anomaly(lon_points,lat_points))

      allocate (longxy_cycle (lon_points,lat_points))

      xmap_anomaly = 0
      ymap_anomaly = 0

      if (p_master) then
         write(6,*), "Open anomaly file: "//trim(fscale_dswrf)
         call sanity(nf90_open(path=trim(fscale_dswrf),mode=nf90_nowrite,ncid=fid_scale_dswrf))
         call sanity(nf90_inq_varid(fid_scale_dswrf,"rsds",vid_scale_dswrf))
   
         write(6,*), "Open anomaly file: "//trim(fscale_dlwrf)
         call sanity(nf90_open(path=trim(fscale_dlwrf),mode=nf90_nowrite,ncid=fid_scale_dlwrf))
         call sanity(nf90_inq_varid(fid_scale_dlwrf,"rlds",vid_scale_dlwrf))
   
         write(6,*), "Open anomaly file: "//trim(fscale_prcp)
         call sanity(nf90_open(path=trim(fscale_prcp),mode=nf90_nowrite,ncid=fid_scale_prcp))
         call sanity(nf90_inq_varid(fid_scale_prcp,"pr",vid_scale_prcp))
   
         write(6,*), "Open anomaly file: "//trim(fanomaly_tair)
         call sanity(nf90_open(path=trim(fanomaly_tair),mode=nf90_nowrite,ncid=fid_anomaly_tair))
         call sanity(nf90_inq_varid(fid_anomaly_tair,"tas",vid_anomaly_tair))
   
         write(6,*), "Open anomaly file: "//trim(fanomaly_qair)
         call sanity(nf90_open(path=trim(fanomaly_qair),mode=nf90_nowrite,ncid=fid_anomaly_qair))
         call sanity(nf90_inq_varid(fid_anomaly_qair,"huss",vid_anomaly_qair))
   
         write(6,*), "Open anomaly file: "//trim(fanomaly_pres)
         call sanity(nf90_open(path=trim(fanomaly_pres),mode=nf90_nowrite,ncid=fid_anomaly_pres))
         call sanity(nf90_inq_varid(fid_anomaly_pres,"ps",vid_anomaly_pres))
   
         write(6,*), "Open anomaly file: "//trim(fanomaly_windu)
         call sanity(nf90_open(path=trim(fanomaly_windu),mode=nf90_nowrite,ncid=fid_anomaly_windu))
         call sanity(nf90_inq_varid(fid_anomaly_windu,"uas",vid_anomaly_windu))
   
         write(6,*), "Open anomaly file: "//trim(fanomaly_windv)
         call sanity(nf90_open(path=trim(fanomaly_windv),mode=nf90_nowrite,ncid=fid_anomaly_windv))
         call sanity(nf90_inq_varid(fid_anomaly_windv,"vas",vid_anomaly_windv))
   
         call sanity(nf90_inq_dimid(fid_anomaly_windv,"lon",vid_t))
         call sanity(nf90_inquire_dimension(fid_anomaly_windv,vid_t,len=nlon_anomaly))

         call sanity(nf90_inq_dimid(fid_anomaly_windv,"lat",vid_t))
         call sanity(nf90_inquire_dimension(fid_anomaly_windv,vid_t,len=nlat_anomaly))
      end if

#ifdef SPMD
      call mpi_bcast(nlon_anomaly,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast(nlat_anomaly,1,mpi_integer,0,p_comm,p_err)
#endif

      allocate (lon_anomaly (nlon_anomaly  ))
      allocate (lat_anomaly (nlat_anomaly  ))
      allocate (lonw_anomaly(nlon_anomaly+1))
      allocate (lats_anomaly(nlat_anomaly+1))
   
      if (p_master) then
         call sanity(nf90_inq_varid(fid_anomaly_windv,'lon',vid_t))
         call sanity(nf90_get_var(fid_anomaly_windv,vid_t,lon_anomaly))
         call sanity(nf90_inq_varid(fid_anomaly_windv,'lat',vid_t))
         call sanity(nf90_get_var(fid_anomaly_windv,vid_t,lat_anomaly))
      end if

#ifdef SPMD
      call mpi_bcast(lon_anomaly,size(lon_anomaly),mpi_real,0,p_comm,p_err)
      call mpi_bcast(lat_anomaly,size(lat_anomaly),mpi_real,0,p_comm,p_err)
#endif

    ! Simple grid mapping from land model to anomaly dataset.

      do j = 2, nlat_anomaly
         lats_anomaly(j) = (lat_anomaly(j-1) + lat_anomaly(j))/2.
      end do

      if(lat_anomaly(1) < 0.) then !* South => North
         lats_anomaly(1) = -90.
         lats_anomaly(nlat_anomaly+1) =  90.
      else                         !* North => South
         lats_anomaly(1) =  90.
         lats_anomaly(nlat_anomaly+1) = -90.
      end if

      do i = 2, nlon_anomaly
         lonw_anomaly(i) = (lon_anomaly(i-1) + lon_anomaly(i))/2.
      end do

      lonw_anomaly(1) = lon_anomaly(1) - (lonw_anomaly(2) - lon_anomaly(1))
      lonw_anomaly(nlon_anomaly+1) = lon_anomaly(nlon_anomaly) + (lon_anomaly(nlon_anomaly)-lonw_anomaly(nlon_anomaly))

    ! Assume longitude range is between 0~360

      if ((maxval(longxy) <= 180. .and. maxval(lon_anomaly) >  180.) .or.&
          (maxval(longxy) >  180. .and. maxval(lon_anomaly) <= 180.)) then
         STOP "Longitude range of anomaly dataset and model grid not consistent."
      end if

      minlon = minval(lonw_anomaly)
      maxlon = maxval(lonw_anomaly)

      do j = 1, lat_points
      do i = 1, lon_points
         if (longxy(i,j) < minlon) then
            longxy_cycle(i,j) = longxy(i,j) + 360
         else if (longxy(i,j) >= maxlon) then
            longxy_cycle(i,j) = longxy(i,j) - 360
         else
            longxy_cycle(i,j) = longxy(i,j)
         end if
      end do
      end do

      do j = 1, lat_points
      do i = 1, lon_points
         if (gmask(i,j) /= p_iam) cycle

         do jx = 1, nlat_anomaly
            if(lat_anomaly(1) < 0.) then  !* South => North
               if(lats_anomaly(jx) <= latixy(i,j) .and. latixy(i,j) <  lats_anomaly(jx+1)) then
                  ymap_anomaly(i,j) = jx
               end if
            else                          !* North => South
               if(lats_anomaly(jx) >  latixy(i,j) .and. latixy(i,j) >= lats_anomaly(jx+1)) then
                  ymap_anomaly(i,j) = jx
               end if
            end if
         end do

         do ix = 1, nlon_anomaly
            if(lonw_anomaly(ix) <= longxy_cycle(i,j) .and. longxy_cycle(i,j) < lonw_anomaly(ix+1)) then
               xmap_anomaly(i,j) = ix
            end if
         end do

         if (xmap_anomaly(i,j) == 0 .or. ymap_anomaly(i,j) == 0) then
            write (6,*) 'fatal error in anomaly grids mapping', i, j, latixy(i,j), longxy(i,j)
            STOP
         end if
      end do
      end do

      deallocate (lon_anomaly )
      deallocate (lat_anomaly )
      deallocate (lonw_anomaly)
      deallocate (lats_anomaly)

      deallocate (longxy_cycle)

   END SUBROUTINE open_anomaly_files

   SUBROUTINE read_anomaly_files(lon_points,lat_points,gmask)

      integer, intent(in) :: lon_points      ! longitude points of land model
      integer, intent(in) :: lat_points      ! latitude points of land model
      integer, intent(in) :: gmask(lon_points,lat_points)

      real(r4), allocatable :: dswrf_buf(:,:)
      real(r4), allocatable :: dlwrf_buf(:,:)
      real(r4), allocatable ::  prcp_buf(:,:)
      real(r4), allocatable ::  tair_buf(:,:)
      real(r4), allocatable ::  qair_buf(:,:)
      real(r4), allocatable ::  pres_buf(:,:)
      real(r4), allocatable :: windu_buf(:,:)
      real(r4), allocatable :: windv_buf(:,:)

      integer year, month, day, second
      integer itime_anomaly_last, i, j

      itime_anomaly_last = itime_anomaly

      call get_curr_date(year,month,day,second)

      itime_anomaly = (year-anomaly_begyr)*12 + month

      if (itime_anomaly >0 .and. itime_anomaly /= itime_anomaly_last) then

         allocate (dswrf_buf(nlon_anomaly,nlat_anomaly))
         allocate (dlwrf_buf(nlon_anomaly,nlat_anomaly))
         allocate ( prcp_buf(nlon_anomaly,nlat_anomaly))
         allocate ( tair_buf(nlon_anomaly,nlat_anomaly))
         allocate ( qair_buf(nlon_anomaly,nlat_anomaly))
         allocate ( pres_buf(nlon_anomaly,nlat_anomaly))
         allocate (windu_buf(nlon_anomaly,nlat_anomaly))
         allocate (windv_buf(nlon_anomaly,nlat_anomaly))

         if (p_master) then
            call sanity(nf90_get_var(fid_scale_dswrf  ,vid_scale_dswrf  ,dswrf_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
            call sanity(nf90_get_var(fid_scale_dlwrf  ,vid_scale_dlwrf  ,dlwrf_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
            call sanity(nf90_get_var(fid_scale_prcp   ,vid_scale_prcp   , prcp_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
            call sanity(nf90_get_var(fid_anomaly_tair ,vid_anomaly_tair , tair_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
            call sanity(nf90_get_var(fid_anomaly_qair ,vid_anomaly_qair , qair_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
            call sanity(nf90_get_var(fid_anomaly_pres ,vid_anomaly_pres , pres_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
            call sanity(nf90_get_var(fid_anomaly_windu,vid_anomaly_windu,windu_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
            call sanity(nf90_get_var(fid_anomaly_windv,vid_anomaly_windv,windv_buf,  &
                                     start=(/1,1,itime_anomaly/),count=(/nlon_anomaly,nlat_anomaly,1/)))
         end if

#if(defined SPMD)
         call mpi_bcast(dswrf_buf,size(dswrf_buf),mpi_real,0,p_comm,p_err)
         call mpi_bcast(dlwrf_buf,size(dlwrf_buf),mpi_real,0,p_comm,p_err)
         call mpi_bcast( prcp_buf,size( prcp_buf),mpi_real,0,p_comm,p_err)
         call mpi_bcast( tair_buf,size( tair_buf),mpi_real,0,p_comm,p_err)
         call mpi_bcast( qair_buf,size( qair_buf),mpi_real,0,p_comm,p_err)
         call mpi_bcast( pres_buf,size( pres_buf),mpi_real,0,p_comm,p_err)
         call mpi_bcast(windu_buf,size(windu_buf),mpi_real,0,p_comm,p_err)
         call mpi_bcast(windv_buf,size(windv_buf),mpi_real,0,p_comm,p_err)
#endif

         do j = 1, lat_points
         do i = 1, lon_points
            if (gmask(i,j) /= p_iam) cycle

            dswrf_scale  (i,j) = dswrf_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
            dlwrf_scale  (i,j) = dlwrf_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
             prcp_scale  (i,j) =  prcp_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
             tair_anomaly(i,j) =  tair_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
             qair_anomaly(i,j) =  qair_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
             pres_anomaly(i,j) =  pres_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
            windu_anomaly(i,j) = windu_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
            windv_anomaly(i,j) = windv_buf(xmap_anomaly(i,j),ymap_anomaly(i,j))
         end do
         end do
       
         deallocate (dswrf_buf)
         deallocate (dlwrf_buf)
         deallocate ( prcp_buf)
         deallocate ( tair_buf)
         deallocate ( qair_buf)
         deallocate ( pres_buf)
         deallocate (windu_buf)
         deallocate (windv_buf)

      end if

   END SUBROUTINE read_anomaly_files

   SUBROUTINE close_anomaly_files

      call sanity(nf90_close(fid_scale_dswrf  ))
      call sanity(nf90_close(fid_scale_dlwrf  ))
      call sanity(nf90_close(fid_scale_prcp   ))
      call sanity(nf90_close(fid_anomaly_tair ))
      call sanity(nf90_close(fid_anomaly_qair ))
      call sanity(nf90_close(fid_anomaly_pres ))
      call sanity(nf90_close(fid_anomaly_windu))
      call sanity(nf90_close(fid_anomaly_windv))

      deallocate (xmap_anomaly )
      deallocate (ymap_anomaly )

      deallocate (dswrf_scale  )
      deallocate (dlwrf_scale  )
      deallocate ( prcp_scale  )
      deallocate ( tair_anomaly)
      deallocate ( qair_anomaly)
      deallocate ( pres_anomaly)
      deallocate (windu_anomaly)
      deallocate (windv_anomaly)

   END SUBROUTINE close_anomaly_files

   SUBROUTINE sanity(ret)

      integer, intent(in) :: ret

      if (ret .ne. nf90_NoErr) then
         write(6,*) trim(nf90_strerror(ret)); stop
      end if
   END SUBROUTINE sanity

#endif

END MODULE metdata
