subroutine da_transform_xtoy_tamdar_sfc (grid, iv, y)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   type (domain),     intent(inout) :: grid
   type (iv_type),    intent(in)    :: iv       ! Innovation vector (O-B).
   type (y_type),     intent(inout) :: y        ! y = h (grid%xa) (linear)

   integer :: n        ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_psfc(:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_tamdar_sfc")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
      allocate (model_v(1,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
      allocate (model_t(1,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
      allocate (model_q(1,iv%info(tamdar)%n1:iv%info(tamdar)%n2))
      allocate (model_psfc(iv%info(tamdar)%n1:iv%info(tamdar)%n2))
      ! [1.2] Interpolate horizontally:
#ifdef A2C
      call da_interp_lin_3d (grid%xa%u, iv%info(tamdar_sfc), model_u,'u')
      call da_interp_lin_3d (grid%xa%v, iv%info(tamdar_sfc), model_v,'v')
#else
      call da_interp_lin_3d (grid%xa%u, iv%info(tamdar_sfc), model_u)
      call da_interp_lin_3d (grid%xa%v, iv%info(tamdar_sfc), model_v)
#endif
      call da_interp_lin_3d (grid%xa%t, iv%info(tamdar_sfc), model_t)
      call da_interp_lin_3d (grid%xa%q, iv%info(tamdar_sfc), model_q)

      call da_interp_lin_2d (grid%xa%psfc, iv%info(tamdar_sfc), 1, model_psfc)

      do n=iv%info(tamdar)%n1,iv%info(tamdar)%n2
         y%tamdar_sfc(n)%u = model_u(1,n)
         y%tamdar_sfc(n)%v = model_v(1,n)
         y%tamdar_sfc(n)%t = model_t(1,n)
         y%tamdar_sfc(n)%q = model_q(1,n)
         y%tamdar_sfc(n)%p = model_psfc(n)
      end do
      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc(grid, iv, tamdar_sfc, iv%tamdar_sfc(:), y%tamdar_sfc(:))
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_tamdar_sfc")

end subroutine da_transform_xtoy_tamdar_sfc


