#include <define.h>

module nchistMod

 use precision
 use paramodel
 use colm_varctl
 use colm_varMod
 use spmd
 use spmd_decomp
 use netcdf
#ifdef RTM
 use colm_rtmVar
 use RunoffMod
#endif

 implicit none

 private

 real(r4), parameter :: missing_value = -9999.
 integer,  parameter :: nMaxHist = 6

 type HistVar
    character(len=31) :: vname
    character(len=31) :: units
    character(len=63) :: longname
    character(len=31) :: gridtype    ! LSM2D, LSM3D, LSM3D_LAKE, DGVM2D, DGVM3D, RTMLND2D, RTMOCN2D
    character(len=31) :: flag        ! AVG, INST
    real(r8), pointer :: ptr1d(:)    ! point to fldv%var(:)
    real(r8), pointer :: ptr2d(:,:)  ! point to fldv%var(:,:)
    real(r8), pointer :: buf1d(:)    ! local buffer for 1D calculation
    real(r8), pointer :: buf2d(:,:)  ! local buffer for 2D calculation
    integer           :: nac         ! number of accumulation
    integer           :: ncid        ! variable id of netCDF
    type(HistVar), pointer :: pNext
 end type HistVar

 type HistFile
    logical :: yearly                            ! year interval
    logical :: monthly                           ! month interval
    logical :: daily                             ! day interval
    logical :: sixhourly                         ! 6-hour interval
    logical :: threehourly                       ! 3-hour interval
    logical :: onehourly                         ! 1-hour interval
    logical :: is_valid                          ! valid history tape
    logical :: is_ready                          ! ready to write out
    logical :: is_newyear                        ! another model year coming
#ifdef RTM
    logical :: with_rtm                          ! with RTM coordinates & fluxes to write
#endif
#ifdef DGVM
    logical :: with_dgvm                         ! with DGVM coordinates & fluxes to write
#endif
    type(HistVar), pointer :: ptVar
    character(len=255)     :: fhistory           ! file name of history
    integer                :: fncid              ! file id of netCDF
 end type HistFile

 type(HistFile) :: histArray(nMaxHist)

 interface nchist_init
    module procedure nchist_init
 end interface

 interface nchist_update
    module procedure nchist_update
 end interface

 interface nchist_write
    module procedure nchist_write
 end interface

 interface nchist_exit
    module procedure nchist_exit
 end interface

 public nchist_init, nchist_update, nchist_write, nchist_exit
 public nMaxHist, HistFile, histArray

contains

   subroutine nchist_init

      use colm_varctl
      use colm_varMod, only: numgrid

      integer i

      do i = 1, nMaxHist
         histArray(i)%yearly      = .false.
         histArray(i)%monthly     = .false.
         histArray(i)%daily       = .false.
         histArray(i)%sixhourly   = .false.
         histArray(i)%threehourly = .false.
         histArray(i)%onehourly   = .false.
         histArray(i)%is_valid    = .false.
         histArray(i)%is_ready    = .false.
         histArray(i)%is_newyear  = .false.
         Nullify(histArray(i)%ptVar)
      end do

      i = 1

      if(lhist_1hourly) then
         histArray(i)%onehourly   = .true.
         histArray(i)%is_valid    = .true.
#ifdef RTM
         histArray(i)%with_rtm    = .false.
#endif
#ifdef DGVM
         histArray(i)%with_dgvm   = .false.
#endif
         call nchist_register_1hourly(histArray(i))
      end if

      i = i+1

      if(lhist_3hourly) then
         histArray(i)%threehourly = .true.
         histArray(i)%is_valid    = .true.
#ifdef RTM
         histArray(i)%with_rtm    = .false.
#endif
#ifdef DGVM
         histArray(i)%with_dgvm   = .false.
#endif
         call nchist_register_3hourly(histArray(i))
      end if

      i = i+1

      if(lhist_6hourly) then
         histArray(i)%sixhourly   = .true.
         histArray(i)%is_valid    = .true.
#ifdef RTM
         histArray(i)%with_rtm    = .false.
#endif
#ifdef DGVM
         histArray(i)%with_dgvm   = .false.
#endif
         call nchist_register_6hourly(histArray(i))
      end if

      i = i+1 

      if(lhist_daily) then
         histArray(i)%daily       = .true.
         histArray(i)%is_valid    = .true.
#ifdef RTM
         histArray(i)%with_rtm    = .false.
#endif
#ifdef DGVM
         histArray(i)%with_dgvm   = .true.
#endif
         call nchist_register_daily(histArray(i))
      end if

      i = i+1 

      if(lhist_monthly) then
         histArray(i)%monthly     = .true.
         histArray(i)%is_valid    = .true.
#ifdef RTM
         histArray(i)%with_rtm    = .true.
#endif
#ifdef DGVM
         histArray(i)%with_dgvm   = .true.
#endif
         call nchist_register_monthly(histArray(i))
      end if

      i = i+1

      if(lhist_yearly) then
         histArray(i)%yearly      = .true.
         histArray(i)%is_valid    = .true.
#ifdef RTM
         histArray(i)%with_rtm    = .true.
#endif
#ifdef DGVM
         histArray(i)%with_dgvm   = .true.
#endif
         call nchist_register_yearly(histArray(i))
      end if

   end subroutine nchist_init

   subroutine nchist_exit

      type(HistVar), pointer :: ptVar, ptVarNext
      integer i

      do i = 1, nMaxHist
         if(histArray(i)%is_valid) then
            histArray(i)%is_valid = .false.

            ptVar => histArray(i)%ptVar

            do while(associated(ptVar))
               if(associated(ptVar%buf1d)) deallocate(ptVar%buf1d)
               if(associated(ptVar%buf2d)) deallocate(ptVar%buf2d)

               ptVarNext => ptVar%pNext

               deallocate(ptVar)

               ptVar => ptVarNext
            end do
         end if
      end do

   end subroutine nchist_exit

   subroutine nchist_create(histF,idate,cdate)

      use colm_varMod
#ifdef RTM
      use colm_rtmVar, only: rtmlon, rtmlat
      use RtmMod, only: longxy_r, latixy_r
#endif
      implicit none

      type(HistFile), intent(inout) :: histF
      integer, intent(in) :: idate(3)
      character(len=255), intent(in) :: cdate

      real(r4), allocatable :: lon(:)
      real(r4), allocatable :: lat(:)
      real(r4), allocatable :: soilz(:)
      real(r4), allocatable :: soilzlk(:)
      real(r4), allocatable :: time(:)

      integer dim_lon, dim_lat, dim_soilz, dim_time, dim_soilzlk
      integer vid_lon, vid_lat, vid_soilz, vid_time, vid_soilzlk
      integer vid_year, vid_day, vid_second
      integer vid_longxy, vid_latixy, vid_area, vid_landfrac, vid_landmask

#ifdef DGVM
      integer dim_pft, vid_pft
      integer,allocatable :: pft(:)
#endif

#ifdef RTM
      real(r4), allocatable :: lonrof(:)
      real(r4), allocatable :: latrof(:)
      integer dim_lonrof, dim_latrof
      integer vid_lonrof, vid_latrof
#endif

      type(HistVar), pointer :: ptVar
      integer i, ntime

      ntime = 1

      allocate (lon(lon_points))
      allocate (lat(lat_points))
      allocate (soilz(nl_soil))
      allocate (soilzlk(nl_lake))
      allocate (time(ntime))

      lon = real(longxy(:,1))
      lat = real(latixy(1,:))

      if(nl_soil==10) then
      soilz = (/0.0071006, 0.0279250, 0.062258, 0.118865, 0.212193, &
                0.3660658, 0.6197585, 1.038027, 1.727635, 2.864607/)
      end if
      if(nl_soil==15) then
      soilz = (/0.0071006, 0.0279250, 0.062258, 0.118865, 0.212193, &
                0.3660658, 0.6197585, 1.038027, 1.727635, 2.864607, &
                4.7391570, 7.8297670, 12.92532, 21.32647, 35.17762/)
      end if

      soilzlk = (/0.1, 1., 2., 3., 4., 5., 7., 7., 10.45, 10.45/)  

      time = (/0/)

#ifdef RTM
      if(histF%with_rtm) then
         allocate (lonrof(rtmlon))
         allocate (latrof(rtmlat))

         lonrof = real(longxy_r(:,1))
         latrof = real(latixy_r(1,:))
      end if
#endif

#ifdef DGVM
      if(histF%with_dgvm) then
         allocate (pft(numpft_nat))
         pft(:) = (/(i,i=1,numpft_nat)/)
      end if
#endif

      call sanity(nf90_create(path=trim(histF%fhistory),cmode=nf90_clobber,ncid=histF%fncid))

      call sanity(nf90_def_dim(histF%fncid,'lon',lon_points,dim_lon))
      call sanity(nf90_def_dim(histF%fncid,'lat',lat_points,dim_lat))
      call sanity(nf90_def_dim(histF%fncid,'soilz',nl_soil,dim_soilz))
      call sanity(nf90_def_dim(histF%fncid,'soilzlk',nl_lake,dim_soilzlk))
      call sanity(nf90_def_dim(histF%fncid,'time',NF90_UNLIMITED,dim_time))
#ifdef RTM
      if(histF%with_rtm) then
         call sanity(nf90_def_dim(histF%fncid,'lonrof',rtmlon,dim_lonrof))
         call sanity(nf90_def_dim(histF%fncid,'latrof',rtmlat,dim_latrof))
      end if
#endif
#ifdef DGVM
      if(histF%with_dgvm) then 
         call sanity(nf90_def_dim(histF%fncid,'pft',numpft_nat,dim_pft))
      end if
#endif

      call sanity(nf90_def_var(histF%fncid,'origin_year',NF90_INT,varid=vid_year))
      call sanity(nf90_def_var(histF%fncid,'origin_day',NF90_INT,varid=vid_day))
      call sanity(nf90_def_var(histF%fncid,'origin_second',NF90_INT,varid=vid_second))

      call sanity(nf90_def_var(histF%fncid,'lon',NF90_FLOAT,dimids=(/dim_lon/),varid=vid_lon))
      call sanity(nf90_put_att(histF%fncid,vid_lon,"long_name","coordinate longitude"))
      call sanity(nf90_put_att(histF%fncid,vid_lon,"units","degrees_east"))

      call sanity(nf90_def_var(histF%fncid,'lat',NF90_FLOAT,dimids=(/dim_lat/),varid=vid_lat))
      call sanity(nf90_put_att(histF%fncid,vid_lat,"long_name","coordinate latitude"))
      call sanity(nf90_put_att(histF%fncid,vid_lat,"units","degrees_north"))

#ifdef RTM
      if(histF%with_rtm) then
         call sanity(nf90_def_var(histF%fncid,'lonrof',NF90_FLOAT,dimids=(/dim_lonrof/),varid=vid_lonrof))
         call sanity(nf90_put_att(histF%fncid,vid_lonrof,"long_name","coordinate longitude of RTM"))
         call sanity(nf90_put_att(histF%fncid,vid_lonrof,"units","degrees_east"))

         call sanity(nf90_def_var(histF%fncid,'latrof',NF90_FLOAT,dimids=(/dim_latrof/),varid=vid_latrof))
         call sanity(nf90_put_att(histF%fncid,vid_latrof,"long_name","coordinate latitude of RTM"))
         call sanity(nf90_put_att(histF%fncid,vid_latrof,"units","degrees_north"))
      end if
#endif

#ifdef DGVM
      if(histF%with_dgvm) then 
         call sanity(nf90_def_var(histF%fncid,'pft',NF90_INT,dimids=(/dim_pft/),varid=vid_pft))
         call sanity(nf90_put_att(histF%fncid,vid_pft,"long_name","PFT categories excluding 2 crops and 1 soil types"))
      end if
#endif

      call sanity(nf90_def_var(histF%fncid,'soilz',NF90_FLOAT,dimids=(/dim_soilz/),varid=vid_soilz))
      call sanity(nf90_put_att(histF%fncid,vid_soilz,"long_name","coordinate soil levels"))
      call sanity(nf90_put_att(histF%fncid,vid_soilz,"units","m"))
      call sanity(nf90_def_var(histF%fncid,'soilzlk',NF90_FLOAT,dimids=(/dim_soilzlk/),varid=vid_soilzlk))
      call sanity(nf90_put_att(histF%fncid,vid_soilzlk,"long_name","coordinate lake levels"))
      call sanity(nf90_put_att(histF%fncid,vid_soilzlk,"units","m"))

      call sanity(nf90_def_var(histF%fncid,'time',NF90_FLOAT,dimids=(/dim_time/),varid=vid_time))
      call sanity(nf90_put_att(histF%fncid,vid_time,"long_name","time"))
      call sanity(nf90_put_att(histF%fncid,vid_time,"units","days since "//trim(cdate)))

      call sanity(nf90_def_var(histF%fncid,"longxy",NF90_FLOAT,dimids=(/dim_lon,dim_lat/),varid=vid_longxy))
      call sanity(nf90_put_att(histF%fncid,vid_longxy,"long_name","longitude"))
      call sanity(nf90_put_att(histF%fncid,vid_longxy,"units","degrees_east"))

      call sanity(nf90_def_var(histF%fncid,"latixy",NF90_FLOAT,dimids=(/dim_lon,dim_lat/),varid=vid_latixy))
      call sanity(nf90_put_att(histF%fncid,vid_latixy,"long_name","latitude"))
      call sanity(nf90_put_att(histF%fncid,vid_latixy,"units","degrees_north"))

      call sanity(nf90_def_var(histF%fncid,"area",NF90_FLOAT,dimids=(/dim_lon,dim_lat/),varid=vid_area))
      call sanity(nf90_put_att(histF%fncid,vid_area,"long_name","grid cell areas"))
      call sanity(nf90_put_att(histF%fncid,vid_area,"units","km^2"))

      call sanity(nf90_def_var(histF%fncid,"landfrac",NF90_FLOAT,dimids=(/dim_lon,dim_lat/),varid=vid_landfrac))
      call sanity(nf90_put_att(histF%fncid,vid_landfrac,"long_name","land fraction"))
      call sanity(nf90_put_att(histF%fncid,vid_landfrac,"units","%"))

      call sanity(nf90_def_var(histF%fncid,"landmask",NF90_INT,dimids=(/dim_lon,dim_lat/),varid=vid_landmask))
      call sanity(nf90_put_att(histF%fncid,vid_landmask,"long_name","land/ocean mask"))

      ptVar => histF%ptVar

      do while(associated(ptVar))
         if(trim(ptVar%gridtype).eq."LSM2D" .or. trim(ptVar%gridtype).eq."DGVM2D") then
            call sanity(nf90_def_var(histF%fncid,trim(ptVar%vname),NF90_FLOAT,&
                                     dimids=(/dim_lon,dim_lat,dim_time/),varid=ptVar%ncid))
         else if(trim(ptVar%gridtype).eq."LSM3D")then
            call sanity(nf90_def_var(histF%fncid,trim(ptVar%vname),NF90_FLOAT,&
                                     dimids=(/dim_lon,dim_lat,dim_soilz,dim_time/),varid=ptVar%ncid))
         else if(trim(ptVar%gridtype).eq."LSM3D_LAKE")then
            call sanity(nf90_def_var(histF%fncid,trim(ptVar%vname),NF90_FLOAT,&
                                     dimids=(/dim_lon,dim_lat,dim_soilzlk,dim_time/),varid=ptVar%ncid))
         else if(trim(ptVar%gridtype).eq."DGVM3D")then
            call sanity(nf90_def_var(histF%fncid,trim(ptVar%vname),NF90_FLOAT,&
                                     dimids=(/dim_lon,dim_lat,dim_pft,dim_time/),varid=ptVar%ncid))
         else if(trim(ptVar%gridtype).eq."RTMLND2D" .or. trim(ptVar%gridtype).eq."RTMOCN2D")then
#ifdef RTM
            call sanity(nf90_def_var(histF%fncid,trim(ptVar%vname),NF90_FLOAT,&
                                     dimids=(/dim_lonrof,dim_latrof,dim_time/),varid=ptVar%ncid))
#endif
         end if

         call sanity(nf90_put_att(histF%fncid,ptVar%ncid,"long_name",trim(ptVar%longname)))
         call sanity(nf90_put_att(histF%fncid,ptVar%ncid,"units",trim(ptVar%units)))
         call sanity(nf90_put_att(histF%fncid,ptVar%ncid,"_FillValue",missing_value))

         ptVar => ptVar%pNext
      end do

      call sanity(nf90_enddef(histF%fncid))

      call sanity(nf90_put_var(histF%fncid,vid_lon,lon(:)))
      call sanity(nf90_put_var(histF%fncid,vid_lat,lat(:)))
#ifdef RTM
      if(histF%with_rtm) then
         call sanity(nf90_put_var(histF%fncid,vid_lonrof,lonrof(:)))
         call sanity(nf90_put_var(histF%fncid,vid_latrof,latrof(:)))
      end if
#endif
#ifdef DGVM
      if(histF%with_dgvm) then
         call sanity(nf90_put_var(histF%fncid,vid_pft,pft(:)))
      end if
#endif
      call sanity(nf90_put_var(histF%fncid,vid_time,time(:)))
      call sanity(nf90_put_var(histF%fncid,vid_soilz,soilz(:)))
      call sanity(nf90_put_var(histF%fncid,vid_soilzlk,soilzlk(:)))
      call sanity(nf90_put_var(histF%fncid,vid_year,(/idate(1)/)))
      call sanity(nf90_put_var(histF%fncid,vid_day,(/idate(2)/)))
      call sanity(nf90_put_var(histF%fncid,vid_second,(/idate(3)/)))
      call sanity(nf90_put_var(histF%fncid,vid_longxy,longxy(:,:)))
      call sanity(nf90_put_var(histF%fncid,vid_latixy,latixy(:,:)))
      call sanity(nf90_put_var(histF%fncid,vid_area,area(:,:)))
      call sanity(nf90_put_var(histF%fncid,vid_landfrac,landfrac(:,:)))
      call sanity(nf90_put_var(histF%fncid,vid_landmask,landmask(:,:)))

      deallocate (lon)
      deallocate (lat)
      deallocate (soilz)
      deallocate (soilzlk)
      deallocate (time)

#ifdef RTM
      if(histF%with_rtm) then
         deallocate (lonrof)
         deallocate (latrof)
      end if
#endif

#ifdef DGVM
      if(histF%with_dgvm) then
         deallocate (pft)
      end if
#endif
   end subroutine nchist_create

   subroutine nchist_close(histF)

      type(histFile), intent(inout) :: histF

      call sanity(nf90_close(histF%fncid))

   end subroutine nchist_close

   subroutine nchist_update

      type(HistVar), pointer :: ptVar
      integer i, beg_lnd, end_lnd, beg_ocn, end_ocn

      do i = 1, nMaxHist
         if(histArray(i)%is_valid) then
            ptVar => histArray(i)%ptVar

            do while(associated(ptVar))
               if(trim(ptVar%flag).eq."AVG") then
                  if(trim(ptVar%gridtype).eq."RTMLND2D") then
#ifdef RTM
                     call get_proc_rof_bounds (beg_lnd, end_lnd, beg_ocn, end_ocn)
                     ptVar%buf1d(beg_lnd:end_lnd) = ptVar%buf1d(beg_lnd:end_lnd) + ptVar%ptr1d(beg_lnd:end_lnd)
#endif
                  else if(trim(ptVar%gridtype).eq."RTMOCN2D") then
#ifdef RTM
                     call get_proc_rof_bounds (beg_lnd, end_lnd, beg_ocn, end_ocn)
                     ptVar%buf1d(beg_ocn:end_ocn) = ptVar%buf1d(beg_ocn:end_ocn) + ptVar%ptr1d(beg_ocn:end_ocn)
#endif
                  else if(trim(ptVar%gridtype).eq."LSM2D" .or. trim(ptVar%gridtype).eq."DGVM2D") then
                     ptVar%buf1d(:) = ptVar%buf1d(:) + ptVar%ptr1d(:)
                  else if(trim(ptVar%gridtype).eq."LSM3D" .or. trim(ptVar%gridtype).eq."LSM3D_LAKE" .or. trim(ptVar%gridtype).eq."DGVM3D") then
                     ptVar%buf2d(:,:) = ptVar%buf2d(:,:) + ptVar%ptr2d(:,:)
                  end if

                  ptVar%nac = ptVar%nac + 1
               end if

               ptVar => ptVar%pNext
            end do
         end if
      end do

   end subroutine nchist_update

   subroutine nchist_prewrite

      type(HistVar), pointer :: ptVar, ptVar_tsn, ptVar_nsnow
      integer i

      do i = 1, nMaxHist
         Nullify(ptVar_tsn)
         Nullify(ptVar_nsnow)

         if(histArray(i)%is_valid .and. histArray(i)%is_ready) then  ! %is_ready being set by lpwrite
            ptVar => histArray(i)%ptVar

            do while(associated(ptVar))
               if(trim(ptVar%flag).eq."AVG") then
                  if(trim(ptVar%gridtype).eq."LSM2D" .or. trim(ptVar%gridtype).eq."DGVM2D" .or. &
                     trim(ptVar%gridtype).eq."RTMLND2D" .or. trim(ptVar%gridtype).eq."RTMOCN2D") then
                     ptVar%buf1d(:) = ptVar%buf1d(:)/float(ptVar%nac)
                  else
                     ptVar%buf2d(:,:) = ptVar%buf2d(:,:)/float(ptVar%nac)
                  end if
               else if(trim(ptVar%flag).eq."INST") then
                  if(trim(ptVar%gridtype).eq."LSM2D" .or. trim(ptVar%gridtype).eq."DGVM2D" .or. &
                     trim(ptVar%gridtype).eq."RTMLND2D" .or. trim(ptVar%gridtype).eq."RTMOCN2D") then
                     ptVar%buf1d(:) = ptVar%ptr1d(:)
                  else
                     ptVar%buf2d(:,:) = ptVar%ptr2d(:,:)
                  end if
               end if

               if(trim(ptVar%vname).eq."tsn") ptVar_tsn => ptVar
               if(trim(ptVar%vname).eq."nsnow") ptVar_nsnow => ptVar

               ptVar => ptVar%pNext
            end do

          ! Special treatment for temperature of snow
            if(associated(ptVar_tsn) .and. associated(ptVar_nsnow)) then
               where(ptVar_nsnow%buf1d .gt. 0._r8) 
                  ptVar_tsn%buf1d = ptVar_tsn%buf1d / ptVar_nsnow%buf1d
               elsewhere
                  ptVar_tsn%buf1d = missing_value
               endwhere
            end if
         end if
      end do

   end subroutine nchist_prewrite

   subroutine nchist_write(idate,cdate)

      integer, intent(in) :: idate(3)
      character(len=255), intent(in) :: cdate

      type(HistVar), pointer :: ptVar
      integer i, nlev

      call nchist_prewrite

      do i = 1, nMaxHist
         if(histArray(i)%is_valid .and. histArray(i)%is_ready) then
            if(p_master) call nchist_create(histArray(i),idate,cdate)

            ptVar => histArray(i)%ptVar

            do while(associated(ptVar))
               if(trim(ptVar%gridtype).eq."LSM2D" .or. trim(ptVar%gridtype).eq."DGVM2D") then
                  nlev = 1
                  call nchist_write_field_LSM(histArray(i)%fncid,ptVar,nlev)
               else if(trim(ptVar%gridtype).eq."LSM3D") then
                  nlev = nl_soil
                  call nchist_write_field_LSM(histArray(i)%fncid,ptVar,nlev)
               else if(trim(ptVar%gridtype).eq."LSM3D_LAKE") then
                  nlev = nl_lake
                  call nchist_write_field_LSM(histArray(i)%fncid,ptVar,nlev)
               else if(trim(ptVar%gridtype).eq."DGVM3D") then
                  nlev = numpft_nat
                  call nchist_write_field_LSM(histArray(i)%fncid,ptVar,nlev)
               else if(trim(ptVar%gridtype).eq."RTMLND2D" .or. trim(ptVar%gridtype).eq."RTMOCN2D") then
                  nlev = 1
#ifdef RTM
                  call nchist_write_field_RTM(histArray(i)%fncid,ptVar,nlev)
#endif
               end if

               ptVar%nac = 0
               if(associated(ptVar%buf1d)) ptVar%buf1d = 0.
               if(associated(ptVar%buf2d)) ptVar%buf2d = 0.

               ptVar => ptVar%pNext
            end do

            if(p_master) call nchist_close(histArray(i))
         end if
      end do

   end subroutine nchist_write

   subroutine nchist_write_field_LSM(fncid,ptVar,nlev)

#ifdef AUTOMASK
      use metdata, only: forcmask
#endif

      integer, intent(in) :: fncid
      type(HistVar), pointer, intent(in) :: ptVar
      integer, intent(in) :: nlev

      real(r8) globbuf(nlev,numgrid_glob)
      real(r4) fldxy(lon_points,lat_points,nlev)

      integer i, j, g

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)*nlev
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))*nlev
      end do
#endif

      if(nlev.eq.1) then
#ifdef AUTOMASK
         do g = 1, numgrid
            i = gxmap(g)
            j = gymap(g)
            if(forcmask(i,j).eq.1) ptVar%buf1d(g) = missing_value
         end do
#endif
#ifdef SPMD
         call mpi_gatherv (ptVar%buf1d,size(ptVar%buf1d),mpi_real8,&
                           globbuf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         globbuf(1,:) = ptVar%buf1d
#endif
      else
#ifdef AUTOMASK
         do g = 1, numgrid
            i = gxmap(g)
            j = gymap(g)
            if(forcmask(i,j).eq.1) ptVar%buf2d(:,g) = missing_value
         end do
#endif
#ifdef SPMD
         call mpi_gatherv (ptVar%buf2d,size(ptVar%buf2d),mpi_real8,&
                           globbuf,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
         globbuf(:,:) = ptVar%buf2d
#endif
      end if

      fldxy(:,:,:) = missing_value
      do g = 1, numgrid_glob
         i = gxmap_glob(g)
         j = gymap_glob(g)
         fldxy(i,j,:) = globbuf(:,g)
      end do

      if(p_master) then
         if(nlev.eq.1) then
            call sanity(nf90_put_var(fncid,ptVar%ncid,fldxy(:,:,1), &
                                     start=(/1,1,1/),count=(/lon_points,lat_points,1/)))
         else
            call sanity(nf90_put_var(fncid,ptVar%ncid,fldxy(:,:,:), &
                                     start=(/1,1,1,1/),count=(/lon_points,lat_points,nlev,1/)))
         end if
      end if

   end subroutine nchist_write_field_LSM

#ifdef RTM
   subroutine nchist_write_field_RTM(fncid,ptVar,nlev)

      use RunoffMod
#ifdef SPMD
      use spmdGathScatMod, only : gather_data_to_master
#endif

      integer, intent(in) :: fncid
      type(HistVar), pointer, intent(in) :: ptVar
      integer, intent(in) :: nlev

      real(r8), pointer :: lnd_globdc(:)
      real(r8), pointer :: ocn_globdc(:)
      real(r4), pointer :: lnd_globxy(:,:)  ! RTM river flow [m3/s]
      real(r4), pointer :: ocn_globxy(:,:)  ! RTM river discharge into ocean [m3/s]

      integer num_lnd, num_ocn, i, j, g

      call get_proc_rof_global (num_lnd, num_ocn)

      if(trim(ptVar%gridtype).eq."RTMLND2D") then
         allocate (lnd_globdc(num_lnd))
         allocate (lnd_globxy(rtmlon,rtmlat))

#ifdef SPMD
         call gather_data_to_master (ptVar%buf1d, lnd_globdc, clmlevel='lndrof')
#else
         lnd_globdc(:) = ptVar%buf1d(:)
#endif

         lnd_globxy(:,:) = missing_value

         do g = 1,num_lnd
            i = runoff%lnd_ixy(g)
            j = runoff%lnd_jxy(g)
            lnd_globxy(i,j) = lnd_globdc(g)
         end do

         if (p_master) call sanity(nf90_put_var(fncid,ptVar%ncid,lnd_globxy, &
                                   start=(/1,1,1/),count=(/rtmlon,rtmlat,1/)))

         deallocate (lnd_globdc)
         deallocate (lnd_globxy)
      end if

      if(trim(ptVar%gridtype).eq."RTMOCN2D") then
         allocate (ocn_globdc(num_ocn))
         allocate (ocn_globxy(rtmlon,rtmlat))

#ifdef SPMD
         call gather_data_to_master (ptVar%buf1d, ocn_globdc, clmlevel='ocnrof')
#else
         ocn_globdc(:) = ptVar%buf1d(:)
#endif

         ocn_globxy(:,:) = missing_value

         do g = 1,num_ocn
            i = runoff%ocn_ixy(g)
            j = runoff%ocn_jxy(g)
            ocn_globxy(i,j) = ocn_globdc(g)
         end do

         if (p_master) call sanity(nf90_put_var(fncid,ptVar%ncid,ocn_globxy, &
                                   start=(/1,1,1/),count=(/rtmlon,rtmlat,1/)))

         deallocate (ocn_globdc)
         deallocate (ocn_globxy)
      end if

   end subroutine nchist_write_field_RTM
#endif

   subroutine nchist_register_field(histF,vname,units,gridtype,flag,longname,ptr1d,ptr2d)

      type(HistFile),   intent(inout) :: histF
      character(len=*), intent(in) :: vname
      character(len=*), intent(in) :: units
      character(len=*), intent(in) :: gridtype
      character(len=*), intent(in) :: flag
      character(len=*), intent(in) :: longname
      real(r8), pointer, optional, intent(in) :: ptr1d(:)
      real(r8), pointer, optional, intent(in) :: ptr2d(:,:)

      type(HistVar), pointer :: ptVar, ptVarPrev
      integer num_lnd, num_ocn, istat

      if(associated(histF%ptVar)) then
         ptVar => histF%ptVar

         do while (associated(ptVar))
            ptVarPrev => ptVar
            ptVar     => ptVar%pNext
         end do

         allocate(ptVarPrev%pNext)
         ptVar => ptVarPrev%pNext
      else
         allocate(histF%ptVar)
         ptVar => histF%ptVar
      end if

      ptVar%vname    = trim(vname)
      ptVar%units    = trim(units)
      ptVar%flag     = trim(flag)
      ptVar%gridtype = trim(gridtype)
      ptVar%longname = trim(longname)
      Nullify(ptVar%buf1d)
      Nullify(ptVar%buf2d)
      Nullify(ptVar%pNext)

#ifdef RTM
      call get_proc_rof_global (num_lnd, num_ocn)
#endif

      if(trim(gridtype).eq."LSM2D")then
         if(.not.present(ptr1d)) stop 'on register LSM2D'
         ptVar%ptr1d => ptr1d
         allocate(ptVar%buf1d(numgrid))
      else if(trim(gridtype).eq."RTMLND2D")then
         if(.not.present(ptr1d)) stop 'on register RTMLND2D'
         ptVar%ptr1d => ptr1d
         allocate(ptVar%buf1d(num_lnd))
      else if(trim(gridtype).eq."RTMOCN2D")then
         if(.not.present(ptr1d)) stop 'on register RTMOCN2D'
         ptVar%ptr1d => ptr1d
         allocate(ptVar%buf1d(num_ocn))
      else if(trim(gridtype).eq."DGVM2D")then
         if(.not.present(ptr1d)) stop 'on register DGVM2D'
         ptVar%ptr1d => ptr1d
         allocate(ptVar%buf1d(numgrid))
      else if(trim(gridtype).eq."LSM3D")then
         if(.not.present(ptr2d)) stop 'on register LSM3D'
         ptVar%ptr2d => ptr2d
         allocate(ptVar%buf2d(nl_soil,numgrid))
      else if(trim(gridtype).eq."LSM3D_LAKE")then
         if(.not.present(ptr2d)) stop 'on register LSM3D_LAKE'
         ptVar%ptr2d => ptr2d
         allocate(ptVar%buf2d(nl_lake,numgrid))
      else if(trim(gridtype).eq."DGVM3D")then
         if(.not.present(ptr2d)) stop 'on register DGVM3D'
         ptVar%ptr2d => ptr2d
         allocate(ptVar%buf2d(numpft_nat,numgrid))
      end if

      ptVar%nac = 0
      if(associated(ptVar%buf1d)) ptVar%buf1d(:)   = 0.
      if(associated(ptVar%buf2d)) ptVar%buf2d(:,:) = 0.

   end subroutine nchist_register_field

   subroutine nchist_register_yearly(histF)

      type(HistFile),intent(inout) :: histF

      call nchist_register_field(histF,vname="taux",units="kg/m/s^2",gridtype="LSM2D",flag="AVG",&
              longname="wind stress: E-W",ptr1d=fldv%taux)
      call nchist_register_field(histF,vname="tauy",units="kg/m/s^2",gridtype="LSM2D",flag="AVG",&
              longname="wind stress: N-S",ptr1d=fldv%tauy)
#ifdef DUST 
      call nchist_register_field(histF,vname="dustemis_bin_1",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 1",ptr1d=fldv%dustemis_bin_1)
      call nchist_register_field(histF,vname="dustemis_bin_2",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 2",ptr1d=fldv%dustemis_bin_2)
      call nchist_register_field(histF,vname="dustemis_bin_3",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 3",ptr1d=fldv%dustemis_bin_3)
      call nchist_register_field(histF,vname="dustemis_bin_4",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 4",ptr1d=fldv%dustemis_bin_4)

      call nchist_register_field(histF,vname="dustemis_total",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="total dust emission flux",ptr1d=fldv%dustemis_total)
#endif

      call nchist_register_field(histF,vname="fsena",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="sensible heat from canopy height to atmosphere",ptr1d=fldv%fsena)
      call nchist_register_field(histF,vname="lfevpa",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="latent heat flux from canopy height to atmosphere",ptr1d=fldv%lfevpa)
      call nchist_register_field(histF,vname="fevpa",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="evapotranspiration rate from canopy height to atmosphere",ptr1d=fldv%fevpa)
      call nchist_register_field(histF,vname="fsenl",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="sensible heat from leaves",ptr1d=fldv%fsenl)
      call nchist_register_field(histF,vname="fevpl",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="evaporation+transpiration rate from leaves",ptr1d=fldv%fevpl)
      call nchist_register_field(histF,vname="etr",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="transpiration rate",ptr1d=fldv%etr)
      call nchist_register_field(histF,vname="fseng",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="sensible heat flux from ground",ptr1d=fldv%fseng)
      call nchist_register_field(histF,vname="fevpg",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="evaporation rate from ground",ptr1d=fldv%fevpg)
      call nchist_register_field(histF,vname="fgrnd",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="ground heat flux",ptr1d=fldv%fgrnd)
      call nchist_register_field(histF,vname="sabvsun",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by sunlit canopy",ptr1d=fldv%sabvsun)
      call nchist_register_field(histF,vname="sabvsha",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by shaded canopy",ptr1d=fldv%sabvsha)
      call nchist_register_field(histF,vname="sabg",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by ground",ptr1d=fldv%sabg)
      call nchist_register_field(histF,vname="olrg",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="outgoing long-wave radiation from ground+canopy",ptr1d=fldv%olrg)
      call nchist_register_field(histF,vname="rnet",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="net radiation",ptr1d=fldv%rnet)
      call nchist_register_field(histF,vname="zerr",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="the error of energy balance",ptr1d=fldv%zerr)
      call nchist_register_field(histF,vname="assim",units="mol/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="canopy assimilation rate",ptr1d=fldv%assim)
      call nchist_register_field(histF,vname="respc",units="mol/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="respiration (plant+soil)",ptr1d=fldv%respc)
      call nchist_register_field(histF,vname="fmicr",units="mol/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="microbial respiration",ptr1d=fldv%fmicr)
      call nchist_register_field(histF,vname="tlsun",units="K",gridtype="LSM2D",flag="AVG",&
              longname="sunlit leaf temperature",ptr1d=fldv%tlsun)
      call nchist_register_field(histF,vname="tlsha",units="K",gridtype="LSM2D",flag="AVG",&
              longname="shaded leaf temperature",ptr1d=fldv%tlsha)
      call nchist_register_field(histF,vname="ldew",units="mm",gridtype="LSM2D",flag="AVG",&
              longname="depth of water on foliage",ptr1d=fldv%ldew)
      call nchist_register_field(histF,vname="sigf",units="%",gridtype="LSM2D",flag="AVG",&
              longname="fraction of veg cover, excluding snow-covered veg",ptr1d=fldv%sigf)
      call nchist_register_field(histF,vname="green",units="%",gridtype="LSM2D",flag="AVG",&
              longname="leaf greenness",ptr1d=fldv%green)
      call nchist_register_field(histF,vname="lai",units="m^2/m^2",gridtype="LSM2D",flag="AVG",&
              longname="leaf area index",ptr1d=fldv%lai)
      call nchist_register_field(histF,vname="sai",units="m^2/m^2",gridtype="LSM2D",flag="AVG",&
              longname="stem area index",ptr1d=fldv%sai)
      call nchist_register_field(histF,vname="avsdr",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="visible, direct averaged albedo",ptr1d=fldv%avsdr)
      call nchist_register_field(histF,vname="avsdf",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="visible, diffuse averaged albedo",ptr1d=fldv%avsdf)
      call nchist_register_field(histF,vname="anidr",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="near-infrared, direct averaged albedo",ptr1d=fldv%anidr)
      call nchist_register_field(histF,vname="anidf",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="near-infrared, diffuse averaged albedo",ptr1d=fldv%anidf)
      call nchist_register_field(histF,vname="emis",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="averaged bulk surface emissivity",ptr1d=fldv%emis)
      call nchist_register_field(histF,vname="z0ma",units="m",gridtype="LSM2D",flag="AVG",&
              longname="effective roughness",ptr1d=fldv%z0ma)
      call nchist_register_field(histF,vname="trad",units="K",gridtype="LSM2D",flag="AVG",&
              longname="radiative temperature of surface",ptr1d=fldv%trad)
      call nchist_register_field(histF,vname="ustar",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="u* in similarity theory",ptr1d=fldv%ustar)
      call nchist_register_field(histF,vname="tstar",units="kg/kg",gridtype="LSM2D",flag="AVG",&
              longname="t* in similarity theory",ptr1d=fldv%tstar)
      call nchist_register_field(histF,vname="qstar",units="kg/kg",gridtype="LSM2D",flag="AVG",&
              longname="q* in similarity theory",ptr1d=fldv%qstar)
      call nchist_register_field(histF,vname="zol",units="-",gridtype="LSM2D",flag="AVG",&
              longname="dimensionless height (z/L) used in Monin-Obukhov theory",ptr1d=fldv%zol)
      call nchist_register_field(histF,vname="rib",units="-",gridtype="LSM2D",flag="AVG",&
              longname="bulk Richardson number in surface layer",ptr1d=fldv%rib)
      call nchist_register_field(histF,vname="fm",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for momentum",ptr1d=fldv%fm)
      call nchist_register_field(histF,vname="fh",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for heat",ptr1d=fldv%fh)
      call nchist_register_field(histF,vname="fq",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for moisture",ptr1d=fldv%fq)
      call nchist_register_field(histF,vname="tref",units="K",gridtype="LSM2D",flag="AVG",&
              longname="2 m height air temperature",ptr1d=fldv%tref)
      call nchist_register_field(histF,vname="qref",units="kg/kg",gridtype="LSM2D",flag="AVG",&
              longname="2 m height air specific humidity",ptr1d=fldv%qref)
      call nchist_register_field(histF,vname="u10m",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="10m u-velocity",ptr1d=fldv%u10m)
      call nchist_register_field(histF,vname="v10m",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="10m v-velocity",ptr1d=fldv%v10m)
      call nchist_register_field(histF,vname="f10m",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for momentum at 10m",ptr1d=fldv%f10m)
      call nchist_register_field(histF,vname="xerr",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="the error of water banace",ptr1d=fldv%xerr)
      call nchist_register_field(histF,vname="rsur",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="surface runoff",ptr1d=fldv%rsur)
      call nchist_register_field(histF,vname="rnof",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="total runoff",ptr1d=fldv%rnof)

      call nchist_register_field(histF,vname="tss",units="K",gridtype="LSM3D",flag="AVG",&
              longname="soil temperature",ptr2d=fldv%tss)
      call nchist_register_field(histF,vname="wliq",units="kg/m^2",gridtype="LSM3D",flag="AVG",&
              longname="liquid water in soil layers",ptr2d=fldv%wliq)
      call nchist_register_field(histF,vname="wice",units="kg/m^2",gridtype="LSM3D",flag="AVG",&
              longname="ice lens in soil layers",ptr2d=fldv%wice)
      call nchist_register_field(histF,vname="mrlsl",units="kg/m^2",gridtype="LSM3D",flag="AVG",&
              longname="mass of water of all phases in each soil layer",ptr2d=fldv%mrlsl)

      call nchist_register_field(histF,vname="t_lake",units="K",gridtype="LSM3D_LAKE",flag="AVG",&
              longname="lake temperature",ptr2d=fldv%t_lake)
      call nchist_register_field(histF,vname="lake_icefrac",units="0-1",gridtype="LSM3D_LAKE",flag="AVG",&
              longname="mass fraction of lake layer that is frozen",ptr2d=fldv%lake_icefrac)

      call nchist_register_field(histF,vname="tg",units="K",gridtype="LSM2D",flag="AVG",&
              longname="ground surface temperature",ptr1d=fldv%tg)
      call nchist_register_field(histF,vname="scv",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="snow cover, water equivalent",ptr1d=fldv%scv)
      call nchist_register_field(histF,vname="snowdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="snow depth",ptr1d=fldv%snowdp)
      call nchist_register_field(histF,vname="fsno",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="fraction of snow cover on ground",ptr1d=fldv%fsno)
      call nchist_register_field(histF,vname="us",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="wind in eastward direction",ptr1d=fldv%us)
      call nchist_register_field(histF,vname="vs",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="wind in northward direction",ptr1d=fldv%vs)
      call nchist_register_field(histF,vname="tm",units="K",gridtype="LSM2D",flag="AVG",&
              longname="temperature at reference height",ptr1d=fldv%tm)
      call nchist_register_field(histF,vname="qm",units="K",gridtype="LSM2D",flag="AVG",&
              longname="specific humidity at reference height",ptr1d=fldv%qm)
      call nchist_register_field(histF,vname="prc",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="convective precipitation",ptr1d=fldv%prc)
      call nchist_register_field(histF,vname="prl",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="large scale precipitation",ptr1d=fldv%prl)
      call nchist_register_field(histF,vname="pbot",units="pa",gridtype="LSM2D",flag="AVG",&
              longname="atmospheric pressure at the surface",ptr1d=fldv%pbot)
      call nchist_register_field(histF,vname="frl",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="atmospheric infrared (longwave) radiation",ptr1d=fldv%frl)
      call nchist_register_field(histF,vname="solar",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="downward solar radiation at surface",ptr1d=fldv%solar)

      call nchist_register_field(histF,vname="qsubl",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="sublimation rate from snow pack",ptr1d=fldv%qsubl)
      call nchist_register_field(histF,vname="mrsos",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases in the upper 0.1 meters of soil",ptr1d=fldv%mrsos)
      call nchist_register_field(histF,vname="mrso",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases over all soil layers",ptr1d=fldv%mrso)
      call nchist_register_field(histF,vname="mrfso",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of frozen water over all soil layers",ptr1d=fldv%mrfso)
      call nchist_register_field(histF,vname="lwsnl",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of liquid water of snow layers",ptr1d=fldv%lwsnl)
      call nchist_register_field(histF,vname="snm",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="surface snow melt",ptr1d=fldv%snm)
      call nchist_register_field(histF,vname="tsn",units="K",gridtype="LSM2D",flag="AVG",&
              longname="snow internal temperature",ptr1d=fldv%tsn)
      call nchist_register_field(histF,vname="nsnow",units="-",gridtype="LSM2D",flag="AVG",&
              longname="number of snow events",ptr1d=fldv%nsnow)

      call nchist_register_field(histF,vname="treeFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="tree fraction",ptr1d=fldv%treeFrac)
      call nchist_register_field(histF,vname="shrubFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="shrub fraction",ptr1d=fldv%shrubFrac)
      call nchist_register_field(histF,vname="grassFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="grass fraction",ptr1d=fldv%grassFrac)
      call nchist_register_field(histF,vname="baresoilFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="bare soil fraction",ptr1d=fldv%baresoilFrac)
      call nchist_register_field(histF,vname="residualFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="residual land fraction",ptr1d=fldv%residualFrac)
      call nchist_register_field(histF,vname="soilFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="soil areal fraction on gridlevel",ptr1d=fldv%soilFrac)
      call nchist_register_field(histF,vname="urbanFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="urban areal fraction on gridlevel",ptr1d=fldv%urbanFrac)
      call nchist_register_field(histF,vname="wetlandFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="wetland areal fraction on gridlevel",ptr1d=fldv%wetlandFrac)
      call nchist_register_field(histF,vname="iceFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="ice areal fraction on gridlevel",ptr1d=fldv%icefrac)
      call nchist_register_field(histF,vname="lakeFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="lake areal fraction on gridlevel",ptr1d=fldv%lakefrac)

#ifdef RTM
    ! carefully check here!!!!
      call nchist_register_field(histF,vname="lndrof",units="m^3/s",gridtype="RTMLND2D",flag="AVG",&
              longname="river flow of land",ptr1d=runoff%lnd)
      call nchist_register_field(histF,vname="ocnrof",units="m^3/s",gridtype="RTMOCN2D",flag="AVG",&
              longname="river discharge into ocean",ptr1d=runoff%ocn)
#endif

#ifdef DGVM
    ! monthly variables
      call nchist_register_field(histF,vname="leafc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in leaves",ptr1d=fldv_dgvm%leafc)
      call nchist_register_field(histF,vname="woodc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in wood",ptr1d=fldv_dgvm%woodc)
      call nchist_register_field(histF,vname="rootc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in roots",ptr1d=fldv_dgvm%rootc)
      call nchist_register_field(histF,vname="vegc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in vegetation",ptr1d=fldv_dgvm%vegc)
      call nchist_register_field(histF,vname="litc_ag",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in above-ground litter",ptr1d=fldv_dgvm%litc_ag)
      call nchist_register_field(histF,vname="litc_bg",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in below-ground litter",ptr1d=fldv_dgvm%litc_bg)
      call nchist_register_field(histF,vname="litc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in litter pool",ptr1d=fldv_dgvm%litc)
      call nchist_register_field(histF,vname="soic_fast",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in fast soil pool",ptr1d=fldv_dgvm%soic_fast)
      call nchist_register_field(histF,vname="soic_slow",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in slow soil pool",ptr1d=fldv_dgvm%soic_slow)
      call nchist_register_field(histF,vname="soic",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in soil pool",ptr1d=fldv_dgvm%soic)
      call nchist_register_field(histF,vname="fveg2litter",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="total carbon mass flux from vegetation to litter",ptr1d=fldv_dgvm%fveg2litter)
      call nchist_register_field(histF,vname="flitter2soil",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="total carbon mass flux from litter to soil",ptr1d=fldv_dgvm%flitter2soil)
      call nchist_register_field(histF,vname="flitter2atmos",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="total carbon mass flux from litter to atmosphere",ptr1d=fldv_dgvm%flitter2atmos)
      call nchist_register_field(histF,vname="gpp",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to GPP on land",ptr1d=fldv_dgvm%gpp)
      call nchist_register_field(histF,vname="npp",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to NPP on land",ptr1d=fldv_dgvm%npp)
      call nchist_register_field(histF,vname="nep",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to NEP on land",ptr1d=fldv_dgvm%nep)
      call nchist_register_field(histF,vname="nbp",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to NBP on land",ptr1d=fldv_dgvm%nbp)
      call nchist_register_field(histF,vname="ra",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux into atmosphere due to autotrophic(plant) respiration on land",ptr1d=fldv_dgvm%ra)
      call nchist_register_field(histF,vname="rh",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux into atmosphere due to heterotrophic respiration on land",ptr1d=fldv_dgvm%rh)
      call nchist_register_field(histF,vname="ffirec",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux into atmosphere due to CO2 emission from fire",ptr1d=fldv_dgvm%ffirec)

      call nchist_register_field(histF,vname="pftFrac",units="percent",gridtype="DGVM3D",flag="AVG",&
              longname="faction of each PFT on gridlevel",ptr2d=fldv_dgvm%pftFrac)

    ! yearly variables
      call nchist_register_field(histF,vname="bare",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%bare)
      call nchist_register_field(histF,vname="afirec",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%afirec)
      call nchist_register_field(histF,vname="afiref",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%afiref)
      call nchist_register_field(histF,vname="avegc",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%avegc)
      call nchist_register_field(histF,vname="aestabc",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%aestabc)
      call nchist_register_field(histF,vname="anpp",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%anpp)
      call nchist_register_field(histF,vname="amrh",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%amrh)
      call nchist_register_field(histF,vname="alitcag",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%alitc_ag)
      call nchist_register_field(histF,vname="alitcbg",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%alitc_bg)
      call nchist_register_field(histF,vname="asoicf",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%asoic_fast)
      call nchist_register_field(histF,vname="asoics",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%asoic_slow)

      call nchist_register_field(histF,vname="fpc",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%fpcgrid)
      call nchist_register_field(histF,vname="npp_ind",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%npp_ind)
      call nchist_register_field(histF,vname="lm",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%lm_ind)
      call nchist_register_field(histF,vname="sm",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%sm_ind)
      call nchist_register_field(histF,vname="hm",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%hm_ind)
      call nchist_register_field(histF,vname="rm",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%rm_ind)
      call nchist_register_field(histF,vname="ca",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%crownarea)
      call nchist_register_field(histF,vname="htop",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%htop)
      call nchist_register_field(histF,vname="nind",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%nind)
      call nchist_register_field(histF,vname="laimx",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%lai_ind)
      call nchist_register_field(histF,vname="anngpp",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%gpp_ind)
      call nchist_register_field(histF,vname="annfrmf",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%frmf_ind)
      call nchist_register_field(histF,vname="annfrms",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%frms_ind)
      call nchist_register_field(histF,vname="annfrmr",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%frmr_ind)
      call nchist_register_field(histF,vname="annfrg",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%frg_ind)

#ifdef DyN
      call nchist_register_field(histF,vname="cnleaf",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%afcton_leaf)
      call nchist_register_field(histF,vname="cnsap",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%afcton_sap)
      call nchist_register_field(histF,vname="cnroot",units="-",gridtype="DGVM3D",flag="INST",&
              longname="-",ptr2d=fldv_dgvm%afcton_root)

      call nchist_register_field(histF,vname="an_up",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%an_up_total)
      call nchist_register_field(histF,vname="stress",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%an_stress_total)
      call nchist_register_field(histF,vname="avegn",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%avegn)
      call nchist_register_field(histF,vname="alitnag",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%alitn_ag)
      call nchist_register_field(histF,vname="alitnbg",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%alitn_bg)
      call nchist_register_field(histF,vname="asoin",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%asoin)
      call nchist_register_field(histF,vname="no3",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%soil_no3)
      call nchist_register_field(histF,vname="nh4",units="-",gridtype="DGVM2D",flag="INST",&
              longname="-",ptr1d=fldv_dgvm%soil_nh4)
#endif
#endif

#if (defined FHNP) && (defined GWUSE)
      !Note that GWUSE is still in the testing phase, don't define it.
      call nchist_register_field(histF,vname="GWINTAKE",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="goundwater intake (vegetated landunitsonly)",ptr1d=fldv%ggw ) !wanglh 2018/11
#endif

#if (defined FHNP) && (defined FTF)
!liruichao add
      call nchist_register_field(histF,vname="frostdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="frost depth",ptr1d=fldv%frostdp)
      call nchist_register_field(histF,vname="thawdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="thaw depth",ptr1d=fldv%thawdp)
      call nchist_register_field(histF,vname="frostdp0",units="m",gridtype="LSM2D",flag="AVG",&
              longname="second frost layer depth",ptr1d=fldv%frostdp0)
!end
#endif
! ==========wangy==========0
#if (defined FHNP) && (defined NP)
      call nchist_register_field(histF,vname="Q_N_point",units="kg/s",gridtype="RTMLND2D",&
              flag="AVG",longname="N_point flux:",ptr1d=runoff%solute)
      call nchist_register_field(histF,vname="S_N_point",units="kg",gridtype="RTMLND2D",&
              flag="AVG",longname="N_point storage:",ptr1d=runoff%sltvolr)
#endif
! ==========wangy==========1
   end subroutine nchist_register_yearly

   subroutine nchist_register_monthly(histF)

      type(HistFile),intent(inout) :: histF

      call nchist_register_field(histF,vname="taux",units="kg/m/s^2",gridtype="LSM2D",flag="AVG",&
              longname="wind stress: E-W",ptr1d=fldv%taux)
      call nchist_register_field(histF,vname="tauy",units="kg/m/s^2",gridtype="LSM2D",flag="AVG",&
              longname="wind stress: N-S",ptr1d=fldv%tauy)
#ifdef DUST 
      call nchist_register_field(histF,vname="dustemis_bin_1",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 1",ptr1d=fldv%dustemis_bin_1)
      call nchist_register_field(histF,vname="dustemis_bin_2",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 2",ptr1d=fldv%dustemis_bin_2)
      call nchist_register_field(histF,vname="dustemis_bin_3",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 3",ptr1d=fldv%dustemis_bin_3)
      call nchist_register_field(histF,vname="dustemis_bin_4",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 4",ptr1d=fldv%dustemis_bin_4)
      call nchist_register_field(histF,vname="dustemis_total",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="total dust emission flux",ptr1d=fldv%dustemis_total)
#endif

      call nchist_register_field(histF,vname="fsena",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="sensible heat from canopy height to atmosphere",ptr1d=fldv%fsena)
      call nchist_register_field(histF,vname="lfevpa",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="latent heat flux from canopy height to atmosphere",ptr1d=fldv%lfevpa)
      call nchist_register_field(histF,vname="fevpa",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="evapotranspiration from canopy to atmosphere",ptr1d=fldv%fevpa)
      call nchist_register_field(histF,vname="fsenl",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="sensible heat from leaves",ptr1d=fldv%fsenl)
      call nchist_register_field(histF,vname="fevpl",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="evaporation+transpiration from leaves",ptr1d=fldv%fevpl)
      call nchist_register_field(histF,vname="etr",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="transpiration rate",ptr1d=fldv%etr)
      call nchist_register_field(histF,vname="fseng",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="sensible heat flux from ground",ptr1d=fldv%fseng)
      call nchist_register_field(histF,vname="fevpg",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="evaporation heat flux from ground",ptr1d=fldv%fevpg)
      call nchist_register_field(histF,vname="fgrnd",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="ground heat flux",ptr1d=fldv%fgrnd)
      call nchist_register_field(histF,vname="sabvsun",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by sunlit canopy",ptr1d=fldv%sabvsun)
      call nchist_register_field(histF,vname="sabvsha",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by shaded canopy",ptr1d=fldv%sabvsha)
      call nchist_register_field(histF,vname="sabg",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by ground",ptr1d=fldv%sabg)
      call nchist_register_field(histF,vname="olrg",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="outgoing long-wave radiation from ground+canopy",ptr1d=fldv%olrg)
      call nchist_register_field(histF,vname="rnet",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="net radiation",ptr1d=fldv%rnet)
      call nchist_register_field(histF,vname="zerr",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="the error of energy balance",ptr1d=fldv%zerr)
      call nchist_register_field(histF,vname="assim",units="mol/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="canopy assimilation rate",ptr1d=fldv%assim)
      call nchist_register_field(histF,vname="respc",units="mol/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="respiration (plant+soil)",ptr1d=fldv%respc)
      call nchist_register_field(histF,vname="fmicr",units="mol/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="microbial respiration",ptr1d=fldv%fmicr)
      call nchist_register_field(histF,vname="tlsun",units="K",gridtype="LSM2D",flag="AVG",&
              longname="sunlit leaf temperature",ptr1d=fldv%tlsun)
      call nchist_register_field(histF,vname="tlsha",units="K",gridtype="LSM2D",flag="AVG",&
              longname="shaded leaf temperature",ptr1d=fldv%tlsha)
      call nchist_register_field(histF,vname="ldew",units="mm",gridtype="LSM2D",flag="AVG",&
              longname="depth of water on foliage",ptr1d=fldv%ldew)
      call nchist_register_field(histF,vname="sigf",units="%",gridtype="LSM2D",flag="AVG",&
              longname="fraction of veg cover, excluding snow-covered veg",ptr1d=fldv%sigf)
      call nchist_register_field(histF,vname="green",units="%",gridtype="LSM2D",flag="AVG",&
              longname="leaf greenness",ptr1d=fldv%green)
      call nchist_register_field(histF,vname="lai",units="m^2/m^2",gridtype="LSM2D",flag="AVG",&
              longname="leaf area index",ptr1d=fldv%lai)
      call nchist_register_field(histF,vname="sai",units="m^2/m^2",gridtype="LSM2D",flag="AVG",&
              longname="stem area index",ptr1d=fldv%sai)
      call nchist_register_field(histF,vname="avsdr",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="visible, direct averaged albedo",ptr1d=fldv%avsdr)
      call nchist_register_field(histF,vname="avsdf",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="visible, diffuse averaged albedo",ptr1d=fldv%avsdf)
      call nchist_register_field(histF,vname="anidr",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="near-infrared, direct averaged albedo",ptr1d=fldv%anidr)
      call nchist_register_field(histF,vname="anidf",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="near-infrared, diffuse averaged albedo",ptr1d=fldv%anidf)
      call nchist_register_field(histF,vname="emis",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="averaged bulk surface emissivity",ptr1d=fldv%emis)
      call nchist_register_field(histF,vname="z0ma",units="m",gridtype="LSM2D",flag="AVG",&
              longname="effective roughness",ptr1d=fldv%z0ma)
      call nchist_register_field(histF,vname="trad",units="K",gridtype="LSM2D",flag="AVG",&
              longname="radiative temperature of surface",ptr1d=fldv%trad)
      call nchist_register_field(histF,vname="ustar",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="u* in similarity theory",ptr1d=fldv%ustar)
      call nchist_register_field(histF,vname="tstar",units="kg/kg",gridtype="LSM2D",flag="AVG",&
              longname="t* in similarity theory",ptr1d=fldv%tstar)
      call nchist_register_field(histF,vname="qstar",units="kg/kg",gridtype="LSM2D",flag="AVG",&
              longname="q* in similarity theory",ptr1d=fldv%qstar)
      call nchist_register_field(histF,vname="zol",units="-",gridtype="LSM2D",flag="AVG",&
              longname="dimensionless height (z/L) used in Monin-Obukhov theory",ptr1d=fldv%zol)
      call nchist_register_field(histF,vname="rib",units="-",gridtype="LSM2D",flag="AVG",&
              longname="bulk Richardson number in surface layer",ptr1d=fldv%rib)
      call nchist_register_field(histF,vname="fm",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for momentum",ptr1d=fldv%fm)
      call nchist_register_field(histF,vname="fh",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for heat",ptr1d=fldv%fh)
      call nchist_register_field(histF,vname="fq",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for moisture",ptr1d=fldv%fq)
      call nchist_register_field(histF,vname="tref",units="K",gridtype="LSM2D",flag="AVG",&
              longname="2 m height air temperature",ptr1d=fldv%tref)
      call nchist_register_field(histF,vname="qref",units="kg/kg",gridtype="LSM2D",flag="AVG",&
              longname="2 m height air specific humidity",ptr1d=fldv%qref)
      call nchist_register_field(histF,vname="u10m",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="10m u-velocity",ptr1d=fldv%u10m)
      call nchist_register_field(histF,vname="v10m",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="10m v-velocity",ptr1d=fldv%v10m)
      call nchist_register_field(histF,vname="f10m",units="-",gridtype="LSM2D",flag="AVG",&
              longname="integral of profile function for momentum at 10m",ptr1d=fldv%f10m)
      call nchist_register_field(histF,vname="xerr",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="the error of water banace",ptr1d=fldv%xerr)
      call nchist_register_field(histF,vname="rsur",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="surface runoff",ptr1d=fldv%rsur)
      call nchist_register_field(histF,vname="rnof",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="total runoff",ptr1d=fldv%rnof)

      call nchist_register_field(histF,vname="tss",units="K",gridtype="LSM3D",flag="AVG",&
              longname="soil temperature",ptr2d=fldv%tss)
      call nchist_register_field(histF,vname="wliq",units="kg/m^2",gridtype="LSM3D",flag="AVG",&
              longname="liquid water in soil layers",ptr2d=fldv%wliq)
      call nchist_register_field(histF,vname="wice",units="kg/m^2",gridtype="LSM3D",flag="AVG",&
              longname="ice lens in soil layers",ptr2d=fldv%wice)
      call nchist_register_field(histF,vname="mrlsl",units="kg/m^2",gridtype="LSM3D",flag="AVG",&
              longname="mass of water of all phases in each soil layer",ptr2d=fldv%mrlsl)

      call nchist_register_field(histF,vname="t_lake",units="K",gridtype="LSM3D_LAKE",flag="AVG",&
              longname="lake temperature",ptr2d=fldv%t_lake)
      call nchist_register_field(histF,vname="lake_icefrac",units="0-1",gridtype="LSM3D_LAKE",flag="AVG",&
              longname="mass fraction of lake layer that is frozen",ptr2d=fldv%lake_icefrac)

      call nchist_register_field(histF,vname="tg",units="K",gridtype="LSM2D",flag="AVG",&
              longname="ground surface temperature",ptr1d=fldv%tg)
      call nchist_register_field(histF,vname="scv",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="snow cover, water equivalent",ptr1d=fldv%scv)
      call nchist_register_field(histF,vname="snowdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="snow depth",ptr1d=fldv%snowdp)
      call nchist_register_field(histF,vname="fsno",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="fraction of snow cover on ground",ptr1d=fldv%fsno)
      call nchist_register_field(histF,vname="us",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="wind in eastward direction",ptr1d=fldv%us)
      call nchist_register_field(histF,vname="vs",units="m/s",gridtype="LSM2D",flag="AVG",&
              longname="wind in northward direction",ptr1d=fldv%vs)
      call nchist_register_field(histF,vname="tm",units="K",gridtype="LSM2D",flag="AVG",&
              longname="temperature at reference height",ptr1d=fldv%tm)
      call nchist_register_field(histF,vname="qm",units="K",gridtype="LSM2D",flag="AVG",&
              longname="specific humidity at reference height",ptr1d=fldv%qm)
      call nchist_register_field(histF,vname="prc",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="convective precipitation",ptr1d=fldv%prc)
      call nchist_register_field(histF,vname="prl",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="large scale precipitation",ptr1d=fldv%prl)
      call nchist_register_field(histF,vname="pbot",units="pa",gridtype="LSM2D",flag="AVG",&
              longname="atmospheric pressure at the surface",ptr1d=fldv%pbot)
      call nchist_register_field(histF,vname="frl",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="atmospheric infrared (longwave) radiation",ptr1d=fldv%frl)
      call nchist_register_field(histF,vname="solar",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="downward solar radiation at surface",ptr1d=fldv%solar)

      call nchist_register_field(histF,vname="qsubl",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="sublimation rate from snow pack",ptr1d=fldv%qsubl)
      call nchist_register_field(histF,vname="mrsos",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases in the upper 0.1 meters of soil",ptr1d=fldv%mrsos)
      call nchist_register_field(histF,vname="mrso",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases over all soil layers",ptr1d=fldv%mrso)
      call nchist_register_field(histF,vname="mrfso",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of frozen water over all soil layers",ptr1d=fldv%mrfso)
      call nchist_register_field(histF,vname="lwsnl",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of liquid water of snow layers",ptr1d=fldv%lwsnl)
      call nchist_register_field(histF,vname="snm",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="surface snow melt",ptr1d=fldv%snm)
      call nchist_register_field(histF,vname="tsn",units="K",gridtype="LSM2D",flag="AVG",&
              longname="snow internal temperature",ptr1d=fldv%tsn)
      call nchist_register_field(histF,vname="nsnow",units="-",gridtype="LSM2D",flag="AVG",&
              longname="number of snow events",ptr1d=fldv%nsnow)

      call nchist_register_field(histF,vname="treeFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="tree fraction",ptr1d=fldv%treeFrac)
      call nchist_register_field(histF,vname="shrubFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="shrub fraction",ptr1d=fldv%shrubFrac)
      call nchist_register_field(histF,vname="grassFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="grass fraction",ptr1d=fldv%grassFrac)
      call nchist_register_field(histF,vname="baresoilFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="bare soil fraction",ptr1d=fldv%baresoilFrac)
      call nchist_register_field(histF,vname="residualFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="residual land fraction",ptr1d=fldv%residualFrac)
      call nchist_register_field(histF,vname="soilFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="soil areal fraction on gridlevel",ptr1d=fldv%soilFrac)
      call nchist_register_field(histF,vname="urbanFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="urban areal fraction on gridlevel",ptr1d=fldv%urbanFrac)
      call nchist_register_field(histF,vname="wetlandFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="wetland areal fraction on gridlevel",ptr1d=fldv%wetlandFrac)
      call nchist_register_field(histF,vname="iceFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="ice areal fraction on gridlevel",ptr1d=fldv%icefrac)
      call nchist_register_field(histF,vname="lakeFrac",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="lake areal fraction on gridlevel",ptr1d=fldv%lakefrac)

#ifdef RTM
    ! carefully check here!!!!
      call nchist_register_field(histF,vname="lndrof",units="m^3/s",gridtype="RTMLND2D",flag="AVG",&
              longname="river flow of land",ptr1d=runoff%lnd)
      call nchist_register_field(histF,vname="ocnrof",units="m^3/s",gridtype="RTMOCN2D",flag="AVG",&
              longname="river discharge into ocean",ptr1d=runoff%ocn)
#endif

#ifdef DGVM
    ! monthly variables
      call nchist_register_field(histF,vname="leafc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in leaves",ptr1d=fldv_dgvm%leafc)
      call nchist_register_field(histF,vname="woodc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in wood",ptr1d=fldv_dgvm%woodc)
      call nchist_register_field(histF,vname="rootc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in roots",ptr1d=fldv_dgvm%rootc)
      call nchist_register_field(histF,vname="vegc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in vegetation",ptr1d=fldv_dgvm%vegc)
      call nchist_register_field(histF,vname="litc_ag",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in above-ground litter",ptr1d=fldv_dgvm%litc_ag)
      call nchist_register_field(histF,vname="litc_bg",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in below-ground litter",ptr1d=fldv_dgvm%litc_bg)
      call nchist_register_field(histF,vname="litc",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in litter pool",ptr1d=fldv_dgvm%litc)
      call nchist_register_field(histF,vname="soic_fast",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in fast soil pool",ptr1d=fldv_dgvm%soic_fast)
      call nchist_register_field(histF,vname="soic_slow",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in slow soil pool",ptr1d=fldv_dgvm%soic_slow)
      call nchist_register_field(histF,vname="soic",units="kg/m^2",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass in soil pool",ptr1d=fldv_dgvm%soic)
      call nchist_register_field(histF,vname="fveg2litter",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="total carbon mass flux from vegetation to litter",ptr1d=fldv_dgvm%fveg2litter)
      call nchist_register_field(histF,vname="flitter2soil",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="total carbon mass flux from litter to soil",ptr1d=fldv_dgvm%flitter2soil)
      call nchist_register_field(histF,vname="flitter2atmos",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="total carbon mass flux from litter to atmosphere",ptr1d=fldv_dgvm%flitter2atmos)
      call nchist_register_field(histF,vname="gpp",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to GPP on land",ptr1d=fldv_dgvm%gpp)
      call nchist_register_field(histF,vname="npp",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to NPP on land",ptr1d=fldv_dgvm%npp)
      call nchist_register_field(histF,vname="nep",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to NEP on land",ptr1d=fldv_dgvm%nep)
      call nchist_register_field(histF,vname="nbp",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux out of atmosphere due to NBP on land",ptr1d=fldv_dgvm%nbp)
      call nchist_register_field(histF,vname="ra",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux into atmosphere due to autotrophic(plant) respiration on land",ptr1d=fldv_dgvm%ra)
      call nchist_register_field(histF,vname="rh",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux into atmosphere due to heterotrophic respiration on land",ptr1d=fldv_dgvm%rh)
      call nchist_register_field(histF,vname="ffirec",units="kg/m^2/s",gridtype="DGVM2D",flag="AVG",&
              longname="carbon mass flux into atmosphere due to CO2 emission from fire",ptr1d=fldv_dgvm%ffirec)

      call nchist_register_field(histF,vname="pftFrac",units="percent",gridtype="DGVM3D",flag="AVG",&
              longname="faction of each PFT on gridlevel",ptr2d=fldv_dgvm%pftFrac)
#endif

#if (defined FHNP) && (defined GWUSE)
      !Note that GWUSE is still in the testing phase, don't define it.
      call nchist_register_field(histF,vname="GWINTAKE",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="goundwater intake (vegetated landunitsonly)",ptr1d=fldv%ggw ) !wanglh 2018/11
#endif

#if (defined FHNP) && (defined FTF)
!liruichao add
      call nchist_register_field(histF,vname="frostdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="frost depth",ptr1d=fldv%frostdp)
      call nchist_register_field(histF,vname="thawdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="thaw depth",ptr1d=fldv%thawdp)
      call nchist_register_field(histF,vname="frostdp0",units="m",gridtype="LSM2D",flag="AVG",&
              longname="second frost layer depth",ptr1d=fldv%frostdp0)
      !end
#endif
! ==========wangy==========0
#if (defined FHNP) && (defined NP)
      call nchist_register_field(histF,vname="Q_N_point",units="kg/s",gridtype="RTMLND2D",&
              flag="AVG",longname="N_point flux:",ptr1d=runoff%solute)
      call nchist_register_field(histF,vname="S_N_point",units="kg",gridtype="RTMLND2D",&
              flag="AVG",longname="N_point storage:",ptr1d=runoff%sltvolr)
#endif
! ==========wangy==========1
   end subroutine nchist_register_monthly

   subroutine nchist_register_daily(histF)

      type(HistFile),intent(inout) :: histF

      call nchist_register_field(histF,vname="sabvsun",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by sunlit canopy",ptr1d=fldv%sabvsun)
      call nchist_register_field(histF,vname="sabvsha",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by shaded canopy",ptr1d=fldv%sabvsha)
      call nchist_register_field(histF,vname="sabg",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="solar absorbed by ground",ptr1d=fldv%sabg)
      call nchist_register_field(histF,vname="lai",units="m^2/m^2",gridtype="LSM2D",flag="AVG",&
              longname="leaf area index",ptr1d=fldv%lai)
      call nchist_register_field(histF,vname="sai",units="m^2/m^2",gridtype="LSM2D",flag="AVG",&
              longname="stem area index",ptr1d=fldv%sai)
      call nchist_register_field(histF,vname="solar",units="W/m^2",gridtype="LSM2D",flag="AVG",&
              longname="downward solar radiation at surface",ptr1d=fldv%solar)
      call nchist_register_field(histF,vname="snowdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="snow depth",ptr1d=fldv%snowdp)
      call nchist_register_field(histF,vname="sigf",units="%",gridtype="LSM2D",flag="AVG",&
              longname="fraction of veg cover, excluding snow-covered veg",ptr1d=fldv%sigf)

      call nchist_register_field(histF,vname="mrsos",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases in the upper 0.1 meters of soil",ptr1d=fldv%mrsos)
      call nchist_register_field(histF,vname="fsno",units="percent",gridtype="LSM2D",flag="AVG",&
              longname="fraction of snow cover on ground",ptr1d=fldv%fsno)
      call nchist_register_field(histF,vname="tg",units="K",gridtype="LSM2D",flag="AVG",&
              longname="ground surface temperature",ptr1d=fldv%tg)
      call nchist_register_field(histF,vname="tm",units="K",gridtype="LSM2D",flag="AVG",&
              longname="temperature at reference height",ptr1d=fldv%tm)
      call nchist_register_field(histF,vname="scv",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="snow cover, water equivalent",ptr1d=fldv%scv)
      call nchist_register_field(histF,vname="rnof",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="total runoff",ptr1d=fldv%rnof)

#ifdef DGVM
      call nchist_register_field(histF,vname="pftFrac",units="percent",gridtype="DGVM3D",flag="AVG",&
              longname="faction of each PFT on gridlevel",ptr2d=fldv_dgvm%pftFrac)
#endif


#if (defined FHNP) && (defined GWUSE)
      !Note that GWUSE is still in the testing phase, don't define it.
            call nchist_register_field(histF,vname="GWINTAKE",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="goundwater intake (vegetated landunitsonly)",ptr1d=fldv%ggw ) !wanglh 2018/11
#endif



#if (defined FHNP) && (defined FTF)
!liruichao add
      call nchist_register_field(histF,vname="frostdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="frost depth",ptr1d=fldv%frostdp)
      call nchist_register_field(histF,vname="thawdp",units="m",gridtype="LSM2D",flag="AVG",&
              longname="thaw depth",ptr1d=fldv%thawdp)
      call nchist_register_field(histF,vname="frostdp0",units="m",gridtype="LSM2D",flag="AVG",&
              longname="second frost layer depth",ptr1d=fldv%frostdp0)
      !end
#endif

#ifdef DUST 
      call nchist_register_field(histF,vname="dustemis_bin_1",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 1",ptr1d=fldv%dustemis_bin_1)
      call nchist_register_field(histF,vname="dustemis_bin_2",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 2",ptr1d=fldv%dustemis_bin_2)
      call nchist_register_field(histF,vname="dustemis_bin_3",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 3",ptr1d=fldv%dustemis_bin_3)
      call nchist_register_field(histF,vname="dustemis_bin_4",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="dust emission flux for bin 4",ptr1d=fldv%dustemis_bin_4)
      call nchist_register_field(histF,vname="dustemis_total",units="kg/m^2/s",gridtype="LSM2D",flag="AVG",&
              longname="total dust emission flux",ptr1d=fldv%dustemis_total)
#endif

   end subroutine nchist_register_daily

   subroutine nchist_register_6hourly(histF)

      type(HistFile),intent(inout) :: histF

      call nchist_register_field(histF,vname="mrsos",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases in the upper 0.1 meters of soil",ptr1d=fldv%mrsos)
      call nchist_register_field(histF,vname="tg",units="K",gridtype="LSM2D",flag="AVG",&
              longname="ground surface temperature",ptr1d=fldv%tg)
      call nchist_register_field(histF,vname="rnof",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="total runoff",ptr1d=fldv%rnof)

   end subroutine nchist_register_6hourly

   subroutine nchist_register_3hourly(histF)

      type(HistFile),intent(inout) :: histF

      call nchist_register_field(histF,vname="mrsos",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases in the upper 0.1 meters of soil",ptr1d=fldv%mrsos)
      call nchist_register_field(histF,vname="tg",units="K",gridtype="LSM2D",flag="AVG",&
              longname="ground surface temperature",ptr1d=fldv%tg)
      call nchist_register_field(histF,vname="rnof",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="total runoff",ptr1d=fldv%rnof)
#ifdef GWUSE
      !Note that GWUSE is still in the testing phase, don't define it.
      call nchist_register_field(histF,vname="GWINTAKE",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="goundwater intake (vegetated landunitsonly)",ptr1d=fldv%ggw ) !wanglh 2018/11
#endif
   end subroutine nchist_register_3hourly

   subroutine nchist_register_1hourly(histF)

      type(HistFile),intent(inout) :: histF

      call nchist_register_field(histF,vname="mrsos",units="kg/m^2",gridtype="LSM2D",flag="AVG",&
              longname="mass of water of all phases in the upper 0.1 meters of soil",ptr1d=fldv%mrsos)
      call nchist_register_field(histF,vname="tg",units="K",gridtype="LSM2D",flag="AVG",&
              longname="ground surface temperature",ptr1d=fldv%tg)
      call nchist_register_field(histF,vname="rnof",units="mm/s",gridtype="LSM2D",flag="AVG",&
              longname="total runoff",ptr1d=fldv%rnof)

   end subroutine nchist_register_1hourly

   subroutine sanity(ret)

      integer, intent(in) :: ret

      if(ret /= nf90_noerr) then
         write(6,*) trim(nf90_strerror(ret))
         stop
      endif

   end subroutine sanity

end module nchistMod
