module da_wrfvar_top

   !-----------------------------------------------------------------------
   ! Purpose: 
   !-----------------------------------------------------------------------

   use module_configure, only : grid_config_rec_type,model_config_rec, &
      model_to_grid_config_rec, get_config_as_buffer,set_config_as_buffer, &
      initial_config
   use module_domain, only : domain,alloc_and_configure_domain, head_grid, &
      program_name
   use module_driver_constants, only : max_comms
   use module_symbols_util, only : wrfu_finalize, wrfu_initialize, &
      wrfu_cal_gregorian

#if defined(RTTOV) || defined(CRTM)
   use module_radiance, only : satinfo
#endif
#ifdef RTTOV
   use module_radiance, only : coefs, sensor_descriptor
#endif

   use module_state_description, only : num_moist, num_a_moist, num_g_moist, &
      num_dfi_moist, num_scalar, num_a_scalar, num_g_scalar, num_dfi_scalar, &
      num_fdda3d, num_fdda2d, num_ozmixm, num_aerosolc
   use module_tiles, only : set_tiles

#ifdef DM_PARALLEL
   use module_dm, only : local_communicator, local_communicator_x, &
      local_communicator_y, ntasks_x, ntasks_y, data_order_xyz, mytask, &
      ntasks, data_order_xy,wrf_dm_initialize
   use module_comm_dm, only : halo_radar_xa_w_sub,halo_ssmi_xa_sub, &
      halo_sfc_xa_sub, halo_xa_sub, halo_psichi_uv_adj_sub, halo_bal_eqn_adj_sub, &
      halo_psichi_uv_sub, halo_init_sub, halo_psichi_uv_adj_sub, halo_2d_work_sub
#endif

   ! too many namelist options to list
   use da_control
   use da_define_structures, only : y_type, j_type, iv_type, be_type, &
      xbx_type,da_deallocate_background_errors,da_initialize_cv, &
      da_zero_vp_type,da_allocate_y,da_deallocate_observations, &
      da_deallocate_y, da_zero_x
   use da_minimisation, only : da_get_innov_vector,da_minimise_cg, &
      da_minimise_lz, da_write_diagnostics, da_calculate_residual, &
      da_calculate_grady
   use da_obs, only : da_transform_xtoy_adj 
   use da_obs_io, only : da_write_filtered_obs, da_write_obs, da_final_write_obs , &
                         da_write_obs_etkf, da_write_modified_filtered_obs  !cys_change
   use da_par_util, only : da_system,da_copy_tile_dims,da_copy_dims
   use da_physics, only : da_uvprho_to_w_lin
#if defined (CRTM) || defined (RTTOV)
   use da_radiance1, only : num_tovs_before, tovs_recv_pe,tovs_copy_count, &
      tovs_send_pe,tovs_send_count,tovs_recv_start, num_tovs_after, &
      tovs_send_start, da_oi_stats_rad, da_write_oa_rad_ascii, da_setup_satcv
   use da_varbc, only : da_varbc_init,da_varbc_update
#endif
   use da_reporting, only : message, da_warning, da_error, da_message
   use da_setup_structures, only : da_setup_obs_structures, &
      da_setup_background_errors,da_setup_flow_predictors, &
      da_setup_cv, da_scale_background_errors
   use da_test, only : da_check
   use da_tools_serial, only : da_get_unit, da_free_unit
   use da_tracing, only : da_trace_entry, da_trace_exit, da_trace
   use da_transfer_model, only : da_transfer_xatoanalysis,da_setup_firstguess, &
       da_transfer_wrftltoxa_adj
   use da_vtox_transforms, only : da_transform_vtox, da_transform_xtoxa, &
      da_transform_xtoxa_adj
   use da_wrfvar_io, only : da_med_initialdata_input, da_med_initialdata_output
   use da_tools, only : da_set_randomcv

#ifdef CRTM
   use module_radiance, only : crtm_destroy
   use da_crtm, only : channelinfo, sensor_descriptor
#endif

   use da_airep, only : da_oi_stats_airep
   use da_airsr , only : da_oi_stats_airsr
   use da_bogus, only : da_oi_stats_bogus
   use da_buoy , only : da_oi_stats_buoy
   use da_geoamv, only : da_oi_stats_geoamv
   use da_gpspw, only : da_oi_stats_gpspw
   use da_gpsref, only : da_oi_stats_gpsref
   use da_metar, only : da_oi_stats_metar
   use da_pilot, only : da_oi_stats_pilot
   use da_polaramv, only : da_oi_stats_polaramv
   use da_profiler, only : da_oi_stats_profiler
   use da_qscat, only : da_oi_stats_qscat
   use da_mtgirs, only : da_oi_stats_mtgirs
   use da_radar, only : da_oi_stats_radar, da_write_oa_radar_ascii
   use da_satem, only : da_oi_stats_satem
   use da_ships, only : da_oi_stats_ships
   use da_sound, only : da_oi_stats_sound, da_oi_stats_sonde_sfc
   use da_ssmi, only : da_oi_stats_ssmt1, da_oi_stats_ssmt2, da_oi_stats_ssmi_tb, da_oi_stats_ssmi_rv
   use da_synop, only : da_oi_stats_synop  

   use da_wrf_interfaces

   implicit none

   integer :: loop, levels_to_process

   type (domain) , pointer :: keep_grid, grid_ptr, null_domain
   type (grid_config_rec_type), save :: config_flags
   integer                 :: number_at_same_level
   integer                 :: time_step_begin_restart

   integer :: domain_id , fid , oid , idum1 , idum2

#ifdef DM_PARALLEL
   integer                 :: nbytes
   integer, parameter      :: configbuflen = 4* CONFIG_BUF_LEN
   integer                 :: configbuf( configbuflen )
#endif

   character (len=80)      :: rstname

contains

#include "da_wrfvar_init1.inc"
#include "da_wrfvar_init2.inc"
#include "da_wrfvar_run.inc"
#include "da_wrfvar_interface.inc"
#include "da_wrfvar_finalize.inc"
#include "da_solve.inc"

end module da_wrfvar_top
