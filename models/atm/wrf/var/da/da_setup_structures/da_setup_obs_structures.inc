subroutine da_setup_obs_structures( grid, ob, iv, j_cost)

   !---------------------------------------------------------------------------
   ! Purpose: Allocate and read in components of observation structure.
   !---------------------------------------------------------------------------

   implicit none
   
   type (domain),     intent(inout) :: grid   ! Model data
   type (y_type),     intent(out)   :: ob     ! Observation structure.
   type (iv_type),    intent(out)   :: iv     ! O-B structure.
   type (j_type),     intent(inout) :: j_cost ! Cost function.

   integer :: i,j
   character (len=14)   :: ob_name

   if (trace_use) call da_trace_entry("da_setup_obs_structures")

! Default values for fmt_info, fmt_srfc, and fmt_each(YRG, 02/25/2009):

   fmt_info = '(a12,1x,a19,1x,a40,1x,i6,3(f12.3,11x),6x,a40)'
   fmt_srfc = '(7(:,f12.3,i4,f7.2))'
   fmt_each = '(3(f12.3,i4,f7.2),11x,3(f12.3,i4,f7.2),11x,1(f12.3,i4,f7.2))'

   call da_message((/'Set up observations (ob)'/))

   ! Adjust obs switches

   if (use_synopobs .OR. use_shipsobs .OR. use_metarobs .OR. use_pilotobs .OR. &
      use_profilerobs .OR. use_buoyobs .OR. use_soundobs .OR. use_mtgirsobs .OR. use_bogusobs .OR. &
      use_radarobs .OR. use_radar_rv .OR. use_radar_rf .OR. use_satemobs .OR. &
      use_geoamvobs .OR. use_polaramvobs .OR. use_airepobs .OR. use_tamdarobs .OR. &
      use_gpspwobs .OR. use_gpsztdobs .OR. use_gpsrefobs .OR. use_ssmiretrievalobs .OR. &
      use_ssmitbobs .OR. use_ssmt1obs .OR. use_ssmt2obs .OR. use_qscatobs .or. &
      use_airsretobs ) then
 
      use_obsgts = .true.
   else
      use_obsgts = .false.
   end if

   if (use_hirs2obs .OR. use_hirs3obs .OR. use_msuobs .OR. use_amsuaobs .OR. &
      use_amsubobs .OR. use_airsobs .OR. use_eos_amsuaobs .OR. &
      use_hsbobs .OR. use_kma1dvar .OR. use_filtered_rad .or. &
      use_ssmisobs .or. use_hirs4obs .or. use_mhsobs .or. use_pseudo_rad .or. &
      use_simulated_rad ) then
      use_rad = .true.
   else
      use_rad = .false.
   end if

   ! test_dm_exact can be set to .true. to force DM_PARALLEL runs 
   ! to produce results that are bitwise-identical regardless of the number of 
   ! MPI tasks used.  This is useful for validating that multi-processor runs 
   ! are not a source of bugs.  Runtime will be much longer.  This option is 
   ! automatically overridden to .false. for serial or 1-MPI-task runs.  

   if (test_dm_exact) then
      if (num_procs == 1) then
         test_dm_exact = .false.
         write(unit=stdout,fmt='(A)') &
            ' test_dm_exact overridden to .false. for serial or 1-MPI-task run'
      end if
      ! only implmenented for Sound and Synop, so switch other types off
      use_shipsobs         = .false.
      use_metarobs         = .false.
      use_bogusobs         = .false.
      use_pilotobs         = .false.
      use_airepobs         = .false.
      use_geoamvobs        = .false.
      use_polaramvobs      = .false.
      use_buoyobs          = .false.
      use_profilerobs      = .false.
      use_satemobs         = .false.
      use_gpspwobs         = .false.
      use_gpsztdobs        = .false.
      use_gpsrefobs        = .false.
      use_ssmiretrievalobs = .false.
      use_ssmitbobs        = .false.
      use_ssmt1obs         = .false.
      use_ssmt2obs         = .false.
      use_qscatobs         = .false.
      use_hirs2obs         = .false.
      use_hirs3obs         = .false.
      use_hirs4obs         = .false.
      use_mhsobs           = .false.
      use_msuobs           = .false.
      use_amsuaobs         = .false.
      use_amsubobs         = .false.
      use_airsobs          = .false.
      use_eos_amsuaobs     = .false.
      use_hsbobs           = .false.
      use_obsgts           = .false.
      use_rad              = .false.
      use_airsretobs       = .false.
   end if
    
   if (use_pseudo_rad .or. num_pseudo > 0) then
      call da_message((/"Single OBS Test:: Turn off all the OBS switches ***"/))
      use_synopobs         = .false.
      use_shipsobs         = .false.
      use_metarobs         = .false.
      use_soundobs         = .false.
      use_mtgirsobs        = .false.
      use_tamdarobs        = .false.
      use_bogusobs         = .false.
      use_pilotobs         = .false.
      use_airepobs         = .false.
      use_geoamvobs        = .false.
      use_polaramvobs      = .false.
      use_buoyobs          = .false.
      use_profilerobs      = .false.
      use_satemobs         = .false.
      use_gpspwobs         = .false.
      use_gpsztdobs        = .false.
      use_gpsrefobs        = .false.
      use_ssmiretrievalobs = .false.
      use_ssmitbobs        = .false.
      use_ssmt1obs         = .false.
      use_ssmt2obs         = .false.
      use_qscatobs         = .false.
      use_hirs2obs         = .false.
      use_hirs3obs         = .false.
      use_hirs4obs         = .false.
      use_mhsobs           = .false.
      use_msuobs           = .false.
      use_amsuaobs         = .false.
      use_amsubobs         = .false.
      use_airsobs          = .false.
      use_eos_amsuaobs     = .false.
      use_hsbobs           = .false.
      use_obsgts           = .true.
      use_rad              = .false.
      use_airsretobs       = .false.
      ob_format = ob_format_ascii
   end if

   if ( use_pseudo_rad ) then
      use_rad = .true.
      thinning         = .false.
      qc_rad           = .false.
      rtminit_nsensor  = 1
      rtminit_platform = pseudo_rad_platid
      rtminit_satid    = pseudo_rad_satid
      rtminit_sensor   = pseudo_rad_senid
   end if

   if (sfc_assi_options < 1 .OR. sfc_assi_options > 2) then
      write(unit=message(1),fmt='(A,I3)') &
         'Invalid sfc_assi_option = ', sfc_assi_options
      call da_error(__FILE__,__LINE__,message(1:1))
   end if

   !---------------------------------------------------------------------------      
   ! [1.0] Setup and read in fields from first guess:
   !----------------------------------------------------------------------------     

   iv%missing = missing
   ! iv%ptop    = grid%xb%ptop

   iv%total_rad_pixel   = 0
   iv%total_rad_channel = 0

   iv%info(:)%nlocal = 0
   iv%info(:)%ntotal = 0
   do i=1,num_ob_indexes
      iv%info(i)%plocal(:) = 0
      iv%info(i)%ptotal(:) = 0
   end do
   iv%num_inst  = 0 

   ob%nlocal(:) = 0
   ob%ntotal(:) = 0
   ob%num_inst  = 0

   iv%info(:)%max_lev = 1

   ! get FGAT time slots

   allocate ( time_slots(0:num_fgat_time) )
   call da_get_time_slots(num_fgat_time,time_window_min,analysis_date, &
                          time_window_max, time_slots)

   if (use_obsgts) then
      ! Conventional obs can be in BUFR or ascii format
      if (ob_format == ob_format_bufr) then
         call da_message((/'Using BUFR format observation input'/))
         call da_setup_obs_structures_bufr (grid, ob, iv)
      else if (ob_format == ob_format_ascii) then
         call da_message((/'Using ASCII format observation input'/))
         call da_setup_obs_structures_ascii (ob, iv, grid)
      else if (ob_format == ob_format_madis) then
         call da_message((/'Using MADIS format observation input'/))
         call da_setup_obs_structures_madis (ob, iv) 
      end if
   end if


   if (use_rad) then
#if defined(CRTM) || defined(RTTOV)
      ! Radiance files can only be in BUFR
      call da_message((/'Using NCEP BUFR radiance 1b input'/))
      call da_setup_radiance_structures(grid, ob, iv)
#else
      call da_error(__FILE__,__LINE__,(/"options to use radiance data are turned on in wrfvar4 namelist, the code needs to be compiled with CRTM or RTTOV library"/))
#endif
   end if

   ! Summarize observations 

   write(unit=stdout, fmt='(a)')  'Observation summary'
   do i=1,num_fgat_time
      write(unit=stdout, fmt='(3x,a,i2)') 'ob time ', i
      do j=1,num_ob_indexes
         if(j.eq.27)cycle
         if (use_gpsztdobs .and.obs_names(j) == 'gpspw         ' ) then
            ob_name = 'gpsztd        '
         else
            ob_name = obs_names(j)
         endif
         write(unit=stdout, fmt='(6x,a,i10,a,i8,a)') &
            ob_name, iv%info(j)%ptotal(i) - iv%info(j)%ptotal(i-1), ' global,', &
            iv%info(j)%plocal(i) - iv%info(j)%plocal(i-1), ' local'
      end do
   end do
   write(unit=stdout, fmt='(a)') ' '
  
   ! Get horizontal interpolation weights.

   call da_setup_obs_interp_wts (iv) 

   if (trace_use) call da_trace_exit("da_setup_obs_structures")    

end subroutine da_setup_obs_structures


