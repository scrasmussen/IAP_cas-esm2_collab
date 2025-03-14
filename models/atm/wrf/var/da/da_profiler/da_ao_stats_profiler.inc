subroutine da_ao_stats_profiler (stats_unit, iv, re)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,        intent(in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent(inout) :: iv            ! iv
   type (y_type),  intent(in)    :: re            ! A - O

   type (stats_profiler_type)       :: stats
   integer                          :: nu, nv
   integer                          :: n, k
   
   if (trace_use_dull) call da_trace_entry("da_ao_stats_profiler")

   nu = 0
   nv = 0

   stats%maximum%u = maxmin_type (missing_r, 0, 0)
   stats%maximum%v = maxmin_type (missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_profiler1_type(0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(profiler)%nlocal
      if (iv%info(profiler)%proc_domain(1,n)) then
         do k=1, iv%info(profiler)%levels(n)
            call da_stats_calculate (n, k, iv%profiler(n)%u(k)%qc, & 
               re%profiler(n)%u(k), nu, & 
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, k, iv%profiler(n)%v(k)%qc, & 
               re%profiler(n)%v(k), nv, & 
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
         end do
      end if    ! end if (iv%info(profiler)%proc_domain(1,n))
   end do 

   ! do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   iv%nstats(profiler) = nu + nv
   
   call da_proc_stats_combine(stats%average%u, stats%rms_err%u, &
      stats%minimum%u%value, stats%maximum%u%value, &
      stats%minimum%u%n, stats%maximum%u%n, &
      stats%minimum%u%l, stats%maximum%u%l)
   call da_proc_stats_combine(stats%average%v, stats%rms_err%v, &
      stats%minimum%v%value, stats%maximum%v%value, &
      stats%minimum%v%n, stats%maximum%v%n, &
      stats%minimum%v%l, stats%maximum%v%l)
   
   if (rootproc) then 
      if (nu /= 0 .or. nv /= 0) then 
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for profiler'
         call da_print_stats_profiler(stats_unit, nu, nv, stats)
      end if 
   end if
   
   if (trace_use_dull) call da_trace_exit("da_ao_stats_profiler")

end subroutine da_ao_stats_profiler


