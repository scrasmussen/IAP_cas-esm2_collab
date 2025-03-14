module da_define_structures

   !---------------------------------------------------------------------------
   ! Purpose: Collection of routines to define and allocate structures.
   !---------------------------------------------------------------------------

   use module_domain, only: vp_type, x_type

   use da_control, only : anal_type_randomcv, stdout, max_fgat_time, &
      vert_corr, global, vert_evalue,print_detail_be, maxsensor, &
      max_ob_levels, trace_use, num_ob_indexes, kms, kme, &
      vert_corr_1, vert_corr_2, vert_evalue_global, cv_options, &
      put_rand_seed, seed_array1, seed_array2, missing_r, &
      sound, synop, pilot, satem, geoamv, polaramv, airep, gpspw, gpsref, &
      metar, ships, ssmi_rv, ssmi_tb, ssmt1, ssmt2, qscat, profiler, buoy, bogus, &
      mtgirs, tamdar, tamdar_sfc, pseudo, radar, radiance, airsr, sonde_sfc, &
      trace_use_dull,comm, num_pseudo

   use da_tracing, only : da_trace_entry, da_trace_exit
   use da_tools_serial, only : da_array_print

   use da_reporting, only : da_error, message

   implicit none
   
   !--------------------------------------------------------------------------
   ! [2.0] Background field structure definition:
   !--------------------------------------------------------------------------

   type xbx_type
      character (len=4):: mminlu

      integer          :: fft_pad_i          ! Padding to get 2**p 3**q 5**r. (p>=1)
      integer          :: fft_pad_j          ! Padding to get 2**p 3**q 5**r.

      integer          :: pad_num            ! Splitted fft_pad_i on this processor.
      integer          :: pad_inc            ! Pad increment (split over v2y).
      integer, pointer :: pad_loc(:)         ! pad location on this processor.
      integer, pointer :: pad_pos(:)         ! pad position beyond ide for this processor.

      integer          :: fft_ix             ! x-direction FFT number, in 2**p 3**q 5**r.
      integer          :: fft_jy             ! y-direction FFT number, in 2**p 3**q 5**r.

      integer, pointer :: fft_factors_x(:)   ! FFT factors in x direction.
      integer, pointer :: fft_factors_y(:)   ! FFT factors in y direction.

      real, pointer    :: trig_functs_x(:)   ! Trig functions in x direction.
      real, pointer    :: trig_functs_y(:)   ! Trig functions in y direction.

      real             :: psac_mean          ! Mean pressure.
      real, pointer    :: latc_mean(:)       ! Mean latitude.

      real, pointer    :: fft_coeffs(:,:)    ! FFT Coefficients

      real             :: fft_adjoint_factor ! FFT Adjoint factor
      ! spectral transform related variables
      integer          :: inc                ! Vector array increment 
      integer          :: ni
      integer          :: nj
      integer          :: nk
      integer          :: max_wavenumber
      integer          :: lenr
      integer          :: lensav
      integer          :: lenwrk
      integer          :: alp_size
      real, pointer       :: wsave(:)          ! Primes for FFT.
      real, pointer       :: lon(:)            ! Longitude (radians).
      real, pointer       :: sinlon(:)         ! sine(longitude).
      real, pointer       :: coslon(:)         ! cosine(longitude).
      real, pointer       :: lat(:)            ! Latitude (radians, from south).
      real, pointer       :: sinlat(:)         ! sine(latitude).
      real, pointer       :: coslat(:)         ! cosine(latitude).
      real, pointer       :: int_wgts(:)       ! Legendre integration weights.
      real, pointer       :: alp(:)            ! Associated Legendre Polynomial.
   end type xbx_type

   !--------------------------------------------------------------------------
   ! [3.0] Innovation vector structure definition:
   !--------------------------------------------------------------------------

   ! [3.1] Generic sub-structures used in iv_type:

   type field_type
      real                   :: inv             ! Innovation vector
      integer                :: qc              ! Observation QC
      real                   :: error           ! Observational error
      real                   :: sens            ! Sensitivity vector
      real                   :: imp             ! Impact vector
   end type field_type

   type model_loc_type
      type (field_type)       :: slp            ! Pressure in Pa
      ! type (field_type)       :: psfc           ! Pressure in Pa
      ! Remove the following in future (needed now for obs i/o only):
      type (field_type)       :: pw             ! Toatl precipitable water cm

      real                    :: x
      real                    :: y
      integer                 :: i
      integer                 :: j
      real                    :: dx
      real                    :: dxm
      real                    :: dy
      real                    :: dym
      logical                 :: proc_domain
      ! obs_global_index is the original index of this obs in the serial 
      ! code.  It is used to reassemble obs in serial-code-order to replicate 
      ! summation order for bitwise-exact testing of distributed-memory 
      ! parallel configurations.  
      integer                 :: obs_global_index
   end type model_loc_type

   type each_level_type
      real                    :: height         ! Height in m
      integer                 :: height_qc      ! Height QC
      real                    :: zk             ! k-coordinates
      type (field_type)       :: u              ! Wind x-component in m/s
      type (field_type)       :: v              ! Wind y-component in m/s
      type (field_type)       :: p              ! Pressure in Pa
      type (field_type)       :: t              ! Temperature in K
      type (field_type)       :: q              ! Mixing ratio (kg/kg).
      type (field_type)       :: rh             ! Relative humidity (%).
      type (field_type)       :: td             ! dew-point in K
      type (field_type)       :: Speed          ! Wind speed m/s
   end type each_level_type

   type radar_each_level_type
      real                   :: height         ! Height in m
      integer                :: height_qc      ! Height QC
      real                   :: zk             ! MM5 k-coordinates
      type (field_type)      :: rv
      type (field_type)      :: rf
   end type radar_each_level_type

   type info_type
      character (len = 40)   :: name          ! Station name
      character (len = 12)   :: platform      ! Instrument platform
      character (len = 40)   :: id            ! 5 digit station identifer
      character (len = 19)   :: date_char     ! CCYY-MM-DD_HH:MM:SS date
      integer                :: levels        ! number of levels
      real                   :: lat           ! Latitude in degree
      real                   :: lon           ! Longitude in degree
      real                   :: elv           ! Elevation in m
      real                   :: pstar         ! Surface pressure
      real                   :: dhr           ! obs time minus analysis time in hour
   end type info_type

   type infa_type
      integer                             :: max_lev
      integer                             :: nlocal
      integer                             :: ntotal
      integer                             :: plocal(0:max_fgat_time)
      integer                             :: ptotal(0:max_fgat_time)
      integer                             :: n1
      integer                             :: n2
      character (len = 40) , allocatable  :: name(:)       ! Station name
      character (len = 12), allocatable   :: platform(:)   ! Instrument platform
      character (len = 40), allocatable   :: id(:)         ! 5 digit station identifer
      character (len = 19), allocatable   :: date_char(:)  ! CCYY-MM-DD_HH:MM:SS date
      integer, allocatable                :: levels(:)     ! number of levels
      real, allocatable                   :: lat(:,:)      ! Latitude in degree
      real, allocatable                   :: lon(:,:)      ! Longitude in degree
      real, allocatable                   :: elv(:)        ! Elevation in m
      real, allocatable                   :: pstar(:)      ! Surface pressure
      type (field_type), allocatable :: slp(:)         ! Pressure in Pa
      ! type (field_type)       :: psfc(:)           ! Pressure in Pa
      ! Remove the following in future (needed now for obs i/o only):
      type (field_type), allocatable :: pw(:)          ! Total precipitable water cm

      real, allocatable       :: x  (:,:)
      real, allocatable       :: y  (:,:)
      integer, allocatable    :: i  (:,:)
      integer, allocatable    :: j  (:,:)
      integer, allocatable    :: k  (:,:)
      real, allocatable       :: dx (:,:)
      real, allocatable       :: dxm(:,:)
      real, allocatable       :: dy (:,:)
      real, allocatable       :: dym(:,:)
      real, allocatable       :: dz (:,:)
      real, allocatable       :: dzm(:,:)
      real, allocatable       :: zk(:,:)
      logical, allocatable    :: proc_domain(:,:)
      logical, allocatable    :: thinned(:,:)
      ! obs_global_index is the original index of this obs in the serial 
      ! code.  It is used to reassemble obs in serial-code-order to replicate 
      ! summation order for bitwise-exact testing of distributed-memory 
      ! parallel configurations.  
      integer, allocatable                 :: obs_global_index(:)
   end type infa_type

   type stn_loc_type
      real                    :: lon                  ! radar site loc
      real                    :: lat                  ! radar site loc
      real                    :: elv                  ! radar site loc
      real                    :: x                    ! radar site loc
      real                    :: y                    ! radar site loc
      real                    :: zk                   ! radar site loc
   end type stn_loc_type
 
   type radar_type
      type (stn_loc_type)     :: stn_loc

      real, pointer           :: model_p(:)
      real, pointer           :: model_rho(:)
      real, pointer           :: model_qrn(:)
      real                    :: model_ps

      real                  , pointer :: height   (:) ! Height in m
      integer               , pointer :: height_qc(:) ! Height QC

      type (field_type)     , pointer :: rv       (:) ! Radial Velocity
      type (field_type)     , pointer :: rf       (:) ! Reflectivity
   end type radar_type

   type multi_level_type
      type (info_type)                        :: info
      type (model_loc_type)                   :: loc
      type (each_level_type)                  :: each(max_ob_levels)
   end type multi_level_type

   type multi_level_type_BUFR
      type (info_type)                        :: info
      type (model_loc_type)                   :: loc
      type (each_level_type), pointer         :: each(:)
   end type multi_level_type_BUFR

   type radar_stn_type
      character (len = 5)    :: platform      ! Data type
      character (len = 12)   :: name          ! Station name
      character (len = 19)   :: date_char     ! CCYY-MM-DD_HH:MM:SS date
      integer                :: numobs        ! number of Obs
      integer                :: levels        ! number of levels
      real                   :: lat           ! Latitude in degree
      real                   :: lon           ! Longitude in degree
      real                   :: elv           ! Elevation in m
   end type radar_stn_type

   type radar_multi_level_type
      type (radar_stn_type)                   :: stn
      type (info_type)                        :: info
      type (model_loc_type)                   :: loc
      type (radar_each_level_type)            :: each(max_ob_levels)
   end type radar_multi_level_type

   ! [3.2] Innovation vector structure:

   type airep_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! pressure
      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
      type (field_type)     , pointer :: t        (:) ! temperature.
   end type airep_type

   type pilot_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! pressure
      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
   end type pilot_type

   type bogus_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! pressure.
      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
      type (field_type)     , pointer :: t        (:) ! temperature.
      type (field_type)     , pointer :: q        (:) ! q.
      type (field_type)               :: slp          ! sea level pressure.
   end type bogus_type

   type satem_type
      real                            :: ref_p        ! Reference pressure
      real                  , pointer :: p        (:) ! Multi-level pressure

      type (field_type)     , pointer :: thickness(:)     ! Thickness.
      type (field_type)     , pointer :: org_thickness(:) ! To store original Thickness info.
   end type satem_type

   type geoamv_type
      real                  , pointer :: p        (:) ! Height in Pa
      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
   end type geoamv_type

   type polaramv_type
      real                  , pointer :: p        (:) ! Height in Pa
      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
   end type polaramv_type

   type gpsref_type
      real             , pointer :: h  (:)      ! Multi-level height
      type (field_type), pointer :: ref(:)      ! GPS Refractivity
      type (field_type), pointer :: p  (:)      ! Retrieved P from Ref.
      type (field_type), pointer :: t  (:)      ! Retrieved T from Ref.
      type (field_type), pointer :: q  (:)      ! From NCEP analysis.
   end type gpsref_type

   type synop_type
      real                    :: h              ! Height in m
      type (field_type)       :: u              ! u-wind.
      type (field_type)       :: v              ! v-wind.
      type (field_type)       :: t              ! temperature.
      type (field_type)       :: p              ! pressure.
      type (field_type)       :: q              ! q.
   end type synop_type

   type sound_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! pressure.

      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
      type (field_type)     , pointer :: t        (:) ! temperature.
      type (field_type)     , pointer :: q        (:) ! q.
   end type sound_type
     
   type mtgirs_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! pressure.

      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
      type (field_type)     , pointer :: t        (:) ! temperature.
      type (field_type)     , pointer :: q        (:) ! q.
   end type mtgirs_type

   type tamdar_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! pressure.

      type (field_type)     , pointer :: u        (:) ! u-wind.
      type (field_type)     , pointer :: v        (:) ! v-wind.
      type (field_type)     , pointer :: t        (:) ! temperature.
      type (field_type)     , pointer :: q        (:) ! q.
   end type tamdar_type

   type airsr_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! pressure.
      type (field_type)     , pointer :: t        (:) ! temperature.
      type (field_type)     , pointer :: q        (:) ! q.
   end type airsr_type

   type gpspw_type
      type (field_type)       :: tpw  ! Toatl precipitable water cm from GPS
  end type gpspw_type

   type ssmi_rv_type
      type (field_type)       :: Speed          ! Wind speed in m/s
      type (field_type)       :: tpw            ! Toatl precipitable water cm
   end type ssmi_rv_type

   type ssmi_tb_type

      type (field_type)       :: tb19v          ! Brightness T (k) 19V
      type (field_type)       :: tb19h          ! Brightness T (k) 19H
      type (field_type)       :: tb22v          ! Brightness T (k) 22V
      type (field_type)       :: tb37v          ! Brightness T (k) 37V
      type (field_type)       :: tb37h          ! Brightness T (k) 37H
      type (field_type)       :: tb85v          ! Brightness T (k) 85V
      type (field_type)       :: tb85h          ! Brightness T (k) 85H
   end type ssmi_tb_type
   
   type ssmt1_type   
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! Pressure in Pa.
      type (field_type)     , pointer :: t        (:) ! temperature.
   end type ssmt1_type

   type ssmt2_type
      real                  , pointer :: h        (:) ! Height in m
      real                  , pointer :: p        (:) ! Pressure in Pa.
      type (field_type)     , pointer :: rh       (:) ! Relative humidity.
   end type ssmt2_type

   type pseudo_type
      type (field_type)       :: u              ! u-wind.
      type (field_type)       :: v              ! v-wind.
      type (field_type)       :: t              ! Temperature.
      type (field_type)       :: p              ! Pressure.
      type (field_type)       :: q              ! Specific Humidity.
   end type pseudo_type

   type qscat_type
      real                    :: h              ! Height in m
      type (field_type)       :: u              ! u-wind.
      type (field_type)       :: v              ! v-wind.
   end type qscat_type

   type varbc_info_type
      integer              :: platform_id, satellite_id, sensor_id
      integer              :: npredmax
      integer              :: gammapred
      integer              :: nchanl
      integer, pointer     :: nbgerr(:) 
      real,    pointer     :: pred(:,:)
      real,    pointer     :: pred_mean(:)
      real,    pointer     :: pred_std(:)
   end type varbc_info_type
   
   type varbc_type
      integer              :: nobs
      integer              :: npred 
      integer              :: ichanl
      integer, pointer     :: pred_use(:)
      integer, pointer     :: ipred(:)
      integer, pointer     :: index(:)
      real,    pointer     :: param(:)
      real,    pointer     :: bgerr(:) 
      real,    pointer     :: vtox(:,:)
   end type varbc_type
   
   type cv_index_type
      integer              :: ts
      integer              :: nclouds
      integer              :: ncv
      integer, pointer     :: cc(:)
      real, pointer        :: vtox(:,:)
   end type cv_index_type

   type instid_type
      ! Instrument triplet, follow the convension of RTTOV
      integer              :: platform_id, satellite_id, sensor_id
      integer              :: rad_monitoring ! 0 (monitor_off): assimilating
                                             !    (default in Registry.wrfvar),
                                             ! 1 (monitor_on):  monitoring
                                             ! monitor_on and monitor_off defined in da_control.f90
      character(len=20)    :: rttovid_string
      integer              :: num_rad, nchan, nlevels
      integer              :: nchannels, nfrequencies,nbtout
      integer              :: num_rad_glo
      integer, pointer     :: ichan(:)
      real,    pointer     :: tb_inv(:,:)
      integer, pointer     :: tb_qc(:,:)
      real,    pointer     :: tb_error(:,:)
      real,    pointer     :: tb_xb(:,:) 
      real,    pointer     :: tb_sens(:,:)
      real,    pointer     :: tb_imp(:,:)
      real,    pointer     :: rad_xb(:,:)
      real,    pointer     :: rad_obs(:,:)
      real,    pointer     :: rad_ovc(:,:,:)
      integer, pointer     :: scanpos(:)
      integer, pointer     :: scanline(:)
      integer, pointer     :: cloud_flag(:,:)
      integer, pointer     :: rain_flag(:)
      real,    pointer     :: satzen(:) 
      real,    pointer     :: satazi(:) 
      real,    pointer     :: solzen(:) 
      real,    pointer     :: solazi(:) 
      real,    pointer     :: t(:,:)
      real,    pointer     :: q(:,:)
      real,    pointer     :: mr(:,:)
      real,    pointer     :: tm(:,:)
      real,    pointer     :: qm(:,:)
      real,    pointer     :: qrn(:,:)
      real,    pointer     :: qcw(:,:)
      real,    pointer     :: qci(:,:)
      real,    pointer     :: qsn(:,:)
      real,    pointer     :: qgr(:,:)
      real,    pointer     :: qhl(:,:)
      real,    pointer     :: pm(:,:)
      real,    pointer     :: rcw(:,:) ! cloud water effectiv radius
      real,    pointer     :: rci(:,:) ! cloud ice effective radius
      real,    pointer     :: rrn(:,:) ! rain effective radius
      real,    pointer     :: rsn(:,:) ! snow effective radius
      real,    pointer     :: rgr(:,:) ! graupel effective radius
      real,    pointer     :: rhl(:,:) ! hail effective radius
      real,    pointer     :: pf(:,:)  ! full level pressure for CRTM
      real,    pointer     :: emiss(:,:)
      real,    pointer     :: u10(:)
      real,    pointer     :: v10(:)
      real,    pointer     :: t2m(:)
      real,    pointer     :: q2m(:)
      real,    pointer     :: mr2m(:)
      real,    pointer     :: psfc(:)
      real,    pointer     :: ps(:)
      real,    pointer     :: ts(:)
      real,    pointer     :: smois(:)
      real,    pointer     :: tslb(:)
      real,    pointer     :: snowh(:)
      integer, pointer     :: isflg(:)
      integer, pointer     :: ifgat(:)
      integer, pointer     :: landsea_mask(:)
      real,    pointer     :: elevation(:)
      real,    pointer     :: soiltyp(:)
      real,    pointer     :: vegtyp(:)
      real,    pointer     :: vegfra(:)
      real,    pointer     :: clwp(:)
      real,    pointer     :: ps_jacobian(:,:) ! only RTTOV
      real,    pointer     :: ts_jacobian(:,:) ! only over water CRTM
      real,    pointer     :: windspeed_jacobian(:,:) ! only MV and over water CRTM
      real,    pointer     :: emiss_jacobian(:,:)
      real,    pointer     :: gamma_jacobian(:,:)
      real,    pointer     :: t_jacobian(:,:,:)
      real,    pointer     :: q_jacobian(:,:,:)
      real,    pointer     :: water_jacobian(:,:,:) ! water content jacobian
      real,    pointer     :: ice_jacobian(:,:,:)
      real,    pointer     :: rain_jacobian(:,:,:)
      real,    pointer     :: snow_jacobian(:,:,:)
      real,    pointer     :: graupel_jacobian(:,:,:)
      real,    pointer     :: hail_jacobian(:,:,:)
      real,    pointer     :: water_r_jacobian(:,:,:) ! effective radius jacobian
      real,    pointer     :: ice_r_jacobian(:,:,:)
      real,    pointer     :: rain_r_jacobian(:,:,:)
      real,    pointer     :: snow_r_jacobian(:,:,:)
      real,    pointer     :: graupel_r_jacobian(:,:,:)
      real,    pointer     :: hail_r_jacobian(:,:,:)
      real,    pointer     :: water_coverage(:)
      real,    pointer     :: land_coverage(:)
      real,    pointer     :: ice_coverage(:)
      real,    pointer     :: snow_coverage(:)

      type (varbc_info_type)        :: varbc_info
      type (varbc_type),pointer     :: varbc(:)
      type (cv_index_type), pointer :: cv_index(:)
      type (infa_type)              :: info
   end type instid_type

   type iv_type
      integer :: nstats(num_ob_indexes)

      integer :: time

      integer :: num_inst, total_rad_pixel, total_rad_channel

      real    :: synop_ef_u, synop_ef_v, synop_ef_t, synop_ef_p, synop_ef_q
      real    :: metar_ef_u, metar_ef_v, metar_ef_t, metar_ef_p, metar_ef_q
      real    :: ships_ef_u, ships_ef_v, ships_ef_t, ships_ef_p, ships_ef_q
      real    :: geoamv_ef_u, geoamv_ef_v
      real    :: polaramv_ef_u, polaramv_ef_v
      real    :: gpspw_ef_tpw
      real    :: sound_ef_u, sound_ef_v, sound_ef_t, sound_ef_q
      real    :: mtgirs_ef_u, mtgirs_ef_v, mtgirs_ef_t, mtgirs_ef_q
      real    :: tamdar_ef_u, tamdar_ef_v, tamdar_ef_t, tamdar_ef_q
      real    :: tamdar_sfc_ef_u, tamdar_sfc_ef_v, tamdar_sfc_ef_t, tamdar_sfc_ef_p, tamdar_sfc_ef_q
      real    :: airep_ef_u, airep_ef_v, airep_ef_t
      real    :: pilot_ef_u, pilot_ef_v
      real    :: ssmir_ef_speed, ssmir_ef_tpw
      real    :: satem_ef_thickness, ssmt1_ef_t, ssmt2_ef_rh
      real    :: gpsref_ef_ref, gpsref_ef_p, gpsref_ef_t, gpsref_ef_q
      real    :: qscat_ef_u, qscat_ef_v
      real    :: profiler_ef_u, profiler_ef_v
      real    :: buoy_ef_u, buoy_ef_v, buoy_ef_t, buoy_ef_p, buoy_ef_q
      real    :: radar_ef_rv, radar_ef_rf
      real    :: bogus_ef_u, bogus_ef_v, bogus_ef_t, bogus_ef_p, bogus_ef_q, bogus_ef_slp
      real    :: airsr_ef_t,  airsr_ef_q

      type (infa_type) :: info(num_ob_indexes)

      type (airsr_type)    , pointer :: airsr(:)
      type (sound_type)    , pointer :: sound(:)
      type (synop_type)    , pointer :: sonde_sfc(:)
      type (airep_type)    , pointer :: airep(:)
      type (pilot_type)    , pointer :: pilot(:)
      type (satem_type)    , pointer :: satem(:)
      type (geoamv_type)   , pointer :: geoamv(:)
      type (polaramv_type) , pointer :: polaramv(:)
      type (synop_type)    , pointer :: synop(:)
      type (synop_type)    , pointer :: metar(:)
      type (synop_type)    , pointer :: ships(:)
      type (gpspw_type)    , pointer :: gpspw(:)
      type (gpsref_type)   , pointer :: gpsref(:)
      type (ssmi_tb_type)  , pointer :: ssmi_tb(:)
      type (ssmi_rv_type)  , pointer :: ssmi_rv(:)
      type (ssmt1_type)    , pointer :: ssmt1(:)
      type (ssmt2_type)    , pointer :: ssmt2(:)
      type (pseudo_type)   , pointer :: pseudo(:)
      type (qscat_type)    , pointer :: qscat(:)
      type (synop_type)    , pointer :: buoy(:)
      type (pilot_type)    , pointer :: profiler(:)
      type (bogus_type)    , pointer :: bogus(:)
      type (radar_type)    , pointer :: radar(:)
      type (instid_type)   , pointer :: instid(:)
      type (mtgirs_type)   , pointer :: mtgirs(:)
      type (tamdar_type)   , pointer :: tamdar(:)
      type (synop_type)    , pointer :: tamdar_sfc(:)

      real :: missing
      real :: ptop
   end type iv_type

   type number_type
      integer                    :: bad
      integer                    :: miss
      integer                    :: use
   end type number_type

   type bad_info_type
      type (number_type)         :: num
      integer                    :: nn(100000)
      integer                    :: kk(100000)
   end type bad_info_type

   type bad_data_type
      type (bad_info_type)       :: u
      type (bad_info_type)       :: v
      type (bad_info_type)       :: t
      type (bad_info_type)       :: p
      type (bad_info_type)       :: q
      type (bad_info_type)       :: tpw
      type (bad_info_type)       :: Speed
      type (bad_info_type)       :: gpsref
      type (bad_info_type)       :: thickness
      type (bad_info_type)       :: rh
      type (bad_info_type)       :: rv
      type (bad_info_type)       :: rf
      type (bad_info_type)       :: slp
      type (bad_info_type)       :: rad
   end type bad_data_type

   type count_obs_number_type
      integer                                 :: num_used
      integer                                 :: num_outside_iyjx
      integer                                 :: num_max_err_chk
      integer                                 :: num_missing
   end type count_obs_number_type

   !--------------------------------------------------------------------------
   ! [3.0] Observation/residual structure definition:
   !--------------------------------------------------------------------------

   type residual_synop_type
      real :: u                                 ! u-wind.
      real :: v                                 ! v-wind.
      real :: t                                 ! temperature.
      real :: p                                 ! pressure.
      real :: q                                 ! q.
   end type residual_synop_type

   type residual_qscat_type
      real :: u                                 ! u-wind.
      real :: v                                 ! v-wind.
   end type residual_qscat_type

   type residual_geoamv_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
   end type residual_geoamv_type

   type residual_polaramv_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
   end type residual_polaramv_type

   type residual_gpspw_type
      real :: tpw                               ! Total precipitable water.
   end type residual_gpspw_type

   type residual_sound_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
      real, pointer :: t(:)                     ! temperature.
      real, pointer :: q(:)                     ! specific humidity.
   end type residual_sound_type
     
   type residual_mtgirs_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
      real, pointer :: t(:)                     ! temperature.
      real, pointer :: q(:)                     ! specific humidity.
   end type residual_mtgirs_type

   type residual_tamdar_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
      real, pointer :: t(:)                     ! temperature.
      real, pointer :: q(:)                     ! specific humidity.
   end type residual_tamdar_type

   type residual_airsr_type
      real, pointer :: t(:)                     ! temperature.
      real, pointer :: q(:)                     ! specific humidity.
   end type residual_airsr_type

   type residual_airep_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
      real, pointer :: t(:)                     ! temperature.
   end type residual_airep_type

   type residual_pilot_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
   end type residual_pilot_type

   type residual_bogus_type
      real, pointer :: u(:)                     ! u-wind.
      real, pointer :: v(:)                     ! v-wind.
      real, pointer :: t(:)                     ! temperature.
      real, pointer :: q(:)                     ! specific humidity.
      real          :: slp                      ! sea-level pressure.
   end type residual_bogus_type

   type residual_satem_type
      real, pointer :: thickness(:)             ! Thickness.
   end type residual_satem_type

   type residual_gpsref_type
      real, pointer :: ref(:)         ! GPS Refractivity
      real, pointer :: p  (:)         ! GPS Retrived p from Refractivity
      real, pointer :: t  (:)         ! GPS Retrived t from Refractivity
      real, pointer :: q  (:)         ! q from NCEP used by CDAAC in retrieval
   end type residual_gpsref_type

   type residual_ssmi_rv_type
      real                    :: tpw      ! Toatl precipitable water cm
      real                    :: Speed    ! Wind speed m/s
   end type residual_ssmi_rv_type

   type residual_ssmi_tb_type
      real                    :: tb19v          ! Brightness T (k) 19V
      real                    :: tb19h          ! Brightness T (k) 19H
      real                    :: tb22v          ! Brightness T (k) 22V
      real                    :: tb37v          ! Brightness T (k) 37V
      real                    :: tb37h          ! Brightness T (k) 37H
      real                    :: tb85v          ! Brightness T (k) 85V
      real                    :: tb85h          ! Brightness T (k) 85H
   end type residual_ssmi_tb_type
   
   type residual_ssmt1_type
      real, pointer :: t(:)                       ! temperature.
   end type residual_ssmt1_type
   
   type residual_ssmt2_type
      real, pointer :: rh(:)                      ! Relative Humidity.
   end type residual_ssmt2_type

   type residual_pseudo_type
      real :: u                                   ! u-wind.
      real :: v                                   ! v-wind.
      real :: t                                   ! temperature.
      real :: p                                   ! pressure.
      real :: q                                   ! specific humidity.
   end type residual_pseudo_type

   type residual_radar_type
      real, pointer :: rv(:)                    ! rv
      real, pointer :: rf(:)                    ! rf
   end type residual_radar_type

   type residual_instid_type
      integer                          :: num_rad
      integer                          :: nchan
      integer, pointer                 :: ichan (:)
      real, pointer                    :: tb(:,:)
   end type residual_instid_type

   type y_type
      integer :: nlocal(num_ob_indexes)
      integer :: ntotal(num_ob_indexes)

      integer :: num_inst

      type (residual_synop_type),    pointer :: synop(:)
      type (residual_synop_type),    pointer :: metar(:) ! Same as synop type
      type (residual_synop_type),    pointer :: ships(:) ! Same as synop type
      type (residual_geoamv_type),   pointer :: geoamv(:)
      type (residual_polaramv_type), pointer :: polaramv(:)
      type (residual_gpspw_type),    pointer :: gpspw (:)
      type (residual_gpsref_type),   pointer :: gpsref(:)
      type (residual_sound_type),    pointer :: sound(:)
      type (residual_mtgirs_type),   pointer :: mtgirs(:)
      type (residual_tamdar_type),   pointer :: tamdar(:)
      type (residual_synop_type),    pointer :: tamdar_sfc(:)
      type (residual_airsr_type),    pointer :: airsr(:)
      type (residual_bogus_type),    pointer :: bogus(:)
      type (residual_synop_type),    pointer :: sonde_sfc(:) ! Same as synop type
      type (residual_airep_type),    pointer :: airep(:)
      type (residual_pilot_type),    pointer :: pilot(:)
      type (residual_satem_type),    pointer :: satem(:)
      type (residual_ssmi_tb_type),  pointer :: ssmi_tb(:)
      type (residual_ssmi_rv_type),  pointer :: ssmi_rv(:)
      type (residual_ssmt1_type),    pointer :: ssmt1(:)
      type (residual_ssmt2_type),    pointer :: ssmt2(:)
      type (residual_pseudo_type),   pointer :: pseudo(:)
      type (residual_qscat_type),    pointer :: qscat(:)
      type (residual_synop_type),    pointer :: buoy(:) ! Same as synop type
      type (residual_pilot_type),    pointer :: profiler(:) ! Same as pilot type
      type (residual_radar_type),    pointer :: radar(:)
      type (residual_instid_type),   pointer :: instid(:)
   end type y_type

   !--------------------------------------------------------------------------
   ! [4.0] Control variable structure:
   !--------------------------------------------------------------------------

   ! Max/Min type:

   type maxmin_type
      real                       :: value
      integer                    :: n, l
   end type maxmin_type

   !--------------------------------------------------------------------------
   ! [5.0] Control variable structure:
   !--------------------------------------------------------------------------
   
   type jo_type_rad
      integer, pointer :: num_ichan(:)
      real, pointer    :: jo_ichan(:)
   end type jo_type_rad

   type jo_type
      real                :: total
      real                :: synop_u, synop_v, synop_t, synop_p, synop_q
      real                :: metar_u, metar_v, metar_t, metar_p, metar_q
      real                :: ships_u, ships_v, ships_t, ships_p, ships_q
      real                :: geoamv_u, geoamv_v
      real                :: polaramv_u, polaramv_v
      real                :: gpspw_tpw, satem_thickness, gpsref_ref
      real                :: sound_u, sound_v, sound_t, sound_q
      real                :: sonde_sfc_u, sonde_sfc_v, sonde_sfc_t, &
                             sonde_sfc_p, sonde_sfc_q
      real                :: mtgirs_u, mtgirs_v, mtgirs_t, mtgirs_q
      real                :: tamdar_u, tamdar_v, tamdar_t, tamdar_q
      real                :: tamdar_sfc_u, tamdar_sfc_v, tamdar_sfc_t, &
                             tamdar_sfc_p, tamdar_sfc_q
      real                :: airep_u, airep_v, airep_t
      real                :: pilot_u, pilot_v
      real                :: ssmir_speed, ssmir_tpw
      real                :: ssmi_tb19v, ssmi_tb19h, ssmi_tb22v, ssmi_tb37v, &
                             ssmi_tb37h, ssmi_tb85v, ssmi_tb85h
      real                :: ssmt1_t, ssmt2_rh
      real                :: pseudo_u, pseudo_v, pseudo_t, pseudo_p, pseudo_q
      real                :: qscat_u, qscat_v
      real                :: profiler_u, profiler_v
      real                :: buoy_u, buoy_v, buoy_t, buoy_p, buoy_q
      real                :: radar_rv, radar_rf
      real                :: bogus_u, bogus_v, bogus_t, bogus_q, bogus_slp
      real                :: airsr_t, airsr_q
      type(jo_type_rad), pointer       :: rad(:)
   end type jo_type

   type j_type
      real             :: total
      real             :: jb
      real             :: jc
      real             :: je
      real             :: jp
      real             :: js
      type (jo_type)   :: jo
   end type j_type

   type cv_type
      integer :: size        ! Total size of control variable.
      integer :: size_jb     ! Size of CV array for Jb term.
      integer :: size_je     ! Size of CV array for Je term.
      integer :: size_jp     ! Size of CV array for Jp term.
      integer :: size_js     ! Size of CV array for Js term.
      integer :: size1c      ! Complex size of CV array of 1st variable error.
      integer :: size2c      ! Complex size of CV array of 2nd variable error.
      integer :: size3c      ! Complex size of CV array of 3rd variable error.
      integer :: size4c      ! Complex size of CV array of 4th variable error.
      integer :: size5c      ! Complex size of CV array of 5th variable error.
      integer :: size_alphac ! Size of alpha control variable (complex).
      integer :: size1       ! Size of CV array of 1st variable error.
      integer :: size2       ! Size of CV array of 2nd variable error.
      integer :: size3       ! Size of CV array of 3rd variable error.
      integer :: size4       ! Size of CV array of 4th variable error.
      integer :: size5       ! Size of CV array of 5th variable error.
   end type cv_type

   type be_subtype
      integer           :: mz          ! Vertical truncation of errors.
      integer           :: max_wave    ! Global only - horizontal spectral truncation.
      character*5       :: name        ! Variable name.
      real*8, pointer     :: rf_alpha(:) ! RF scale length.
      real*8, pointer     :: val(:,:)    ! Local Standard dev./sqrt(eigenvalue).
      real*8, pointer     :: evec(:,:,:) ! Local Vertical eigenvectors.
      real*8, pointer     :: val_g(:)    ! Global Standard dev./sqrt(eigenvalue).
      real*8, pointer     :: evec_g(:,:) ! Global Vertical eigenvectors.
      real*8, pointer     :: power(:,:)  ! Power spectrum
   end type be_subtype

   type be_type
      integer           :: ne
      integer           :: max_wave           ! Smallest spectral mode (global).
      integer           :: mix
      integer           :: mjy
      type (be_subtype) :: v1
      type (be_subtype) :: v2
      type (be_subtype) :: v3
      type (be_subtype) :: v4
      type (be_subtype) :: v5
      type (be_subtype) :: alpha
      real*8, pointer     :: pb_vert_reg(:,:,:)

      ! Control variable space errors:
      type (cv_type)    :: cv

      real, pointer     :: reg_chi(:,:)
      real, pointer     :: reg_t  (:,:,:)
      real, pointer     :: reg_ps (:,:)

!-----For cv option 3:
      INTEGER          :: ndeg,nta
      REAL             :: swidth
      REAL, POINTER    :: be(:)
      REAL, POINTER    :: rate(:)
      REAL, POINTER    :: table(:,:)
      REAL, POINTER    :: agvz(:,:,:,:)
      REAL, POINTER    :: bvz(:,:,:)
      REAL, POINTER    :: wgvz(:,:,:)
      REAL, POINTER    :: slix(:,:,:,:)
      REAL, POINTER    :: slipx(:,:)
      REAL, POINTER    :: sljy(:,:,:,:)
      REAL, POINTER    :: sljpy(:,:)
      REAL, POINTER    :: vz(:,:,:,:)
      REAL, POINTER    :: corz(:,:,:,:)
      REAL, POINTER    :: corp(:,:)
   end type be_type

   ! Analysis_Stats maximum-minumum structure.

   type maxmin_field_type
      real                         :: value
      integer                      :: i, j
   end type maxmin_field_type

   ! vp_type is defined in the Registry
   ! x_type  is defined in the Registry
   ! The framework allocates the (local-grid) xa structure.
   ! The framework allocates the (local-grid) xb structure.
   ! The framework (de)allocates the vv structure.
   ! The framework (de)allocates the vp structure.

contains

#include "da_allocate_background_errors.inc"
#include "da_allocate_observations.inc"
#include "da_allocate_y.inc"
#include "da_deallocate_background_errors.inc"
#include "da_deallocate_observations.inc"
#include "da_deallocate_y.inc"
#include "da_zero_x.inc"
#include "da_zero_y.inc"
#include "da_zero_vp_type.inc"
#include "da_initialize_cv.inc"
#include "da_random_seed.inc"
#include "da_gauss_noise.inc"

end module da_define_structures

