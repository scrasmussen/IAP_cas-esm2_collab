module da_profiler

   use module_domain, only : domain
   
   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      check_max_iv_print, check_max_iv_unit, v_interp_p, v_interp_h, &
      check_max_iv, missing_data, max_error_uv, max_error_t, rootproc, &
      profiler, max_error_p,max_error_q, fails_error_max, &
      max_stheight_diff, anal_type_verify, kms,kme,kts,kte, trace_use_dull,&
      ob_vars, qcstat_conv_unit
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, x_type, number_type, bad_data_type
   use da_interpolation, only : da_interp_lin_3d, da_to_zk, &
      da_interp_lin_3d_adj
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk,da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   ! The "stats_profiler_type" is ONLY used locally in da_profiler:

   type residual_profiler1_type
      real          :: u                        ! u-wind.
      real          :: v                        ! v-wind.
   end type residual_profiler1_type

   type maxmin_profiler_stats_type
      type (maxmin_type)         :: u, v
   end type maxmin_profiler_stats_type

   type stats_profiler_type
      type (maxmin_profiler_stats_type)  :: maximum, minimum
      type (residual_profiler1_type)     :: average, rms_err
   end type stats_profiler_type

contains

#include "da_ao_stats_profiler.inc"
#include "da_jo_and_grady_profiler.inc"
#include "da_residual_profiler.inc"
#include "da_oi_stats_profiler.inc"
#include "da_print_stats_profiler.inc"
#include "da_transform_xtoy_profiler.inc"
#include "da_transform_xtoy_profiler_adj.inc"
#include "da_check_max_iv_profiler.inc"
#include "da_get_innov_vector_profiler.inc"
#include "da_calculate_grady_profiler.inc"


end module da_profiler

