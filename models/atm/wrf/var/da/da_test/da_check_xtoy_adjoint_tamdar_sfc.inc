subroutine da_check_xtoy_adjoint_tamdar_sfc(iv, y, adjtest_lhs, pertile_lhs)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)     :: iv            ! obs. inc. vector (o-b).
   type (y_type) , intent(inout)  :: y             ! y = h (xa)
   real          , intent(inout)  :: adjtest_lhs, pertile_lhs   

   integer :: n             ! Loop counter.

   if (trace_use_dull) call da_trace_entry("da_check_xtoy_adjoint_tamdar_sfc")

   do n=1, iv%info(tamdar)%nlocal
      if (iv%info(tamdar)%zk(1,n) < 1.0) then
         y%tamdar_sfc(n)%u = 0.0
         y%tamdar_sfc(n)%v = 0.0
         y%tamdar_sfc(n)%t = 0.0
         y%tamdar_sfc(n)%p = 0.0
         y%tamdar_sfc(n)%q = 0.0

         cycle
      end if

      y%tamdar_sfc(n)%u = y%tamdar_sfc(n)%u/typical_u_rms
      y%tamdar_sfc(n)%v = y%tamdar_sfc(n)%v/typical_v_rms
      y%tamdar_sfc(n)%t = y%tamdar_sfc(n)%t/typical_t_rms
      y%tamdar_sfc(n)%p = y%tamdar_sfc(n)%p/typical_p_rms
      y%tamdar_sfc(n)%q = y%tamdar_sfc(n)%q/typical_q_rms

      if (iv%info(tamdar)%proc_domain(1,n)) then
         adjtest_lhs = adjtest_lhs  &
                     + (y%tamdar_sfc(n)%u)**2 &
                     + (y%tamdar_sfc(n)%v)**2 &
                     + (y%tamdar_sfc(n)%t)**2 &
                     + (y%tamdar_sfc(n)%p)**2 &
                     + (y%tamdar_sfc(n)%q)**2
      end if

      pertile_lhs = pertile_lhs &
                  + (y%tamdar_sfc(n)%u)**2 &
                  + (y%tamdar_sfc(n)%v)**2 &
                  + (y%tamdar_sfc(n)%t)**2 &
                  + (y%tamdar_sfc(n)%p)**2 &
                  + (y%tamdar_sfc(n)%q)**2

      y%tamdar_sfc(n)%u = y%tamdar_sfc(n)%u/typical_u_rms
      y%tamdar_sfc(n)%v = y%tamdar_sfc(n)%v/typical_v_rms
      y%tamdar_sfc(n)%t = y%tamdar_sfc(n)%t/typical_t_rms
      y%tamdar_sfc(n)%p = y%tamdar_sfc(n)%p/typical_p_rms
      y%tamdar_sfc(n)%q = y%tamdar_sfc(n)%q/typical_q_rms
   end do

   if (trace_use_dull) call da_trace_exit("da_check_xtoy_adjoint_tamdar_sfc")

end subroutine da_check_xtoy_adjoint_tamdar_sfc


