subroutine da_transform_xtoy_metar_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid          ! first guess state.
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(inout) :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: n        ! Loop counter.

   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_q(:,:)
   real, allocatable :: model_psfc(:)

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_metar_adj")

   if (sfc_assi_options == sfc_assi_options_1) then
      allocate (model_u(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_v(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_t(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_q(1,iv%info(metar)%n1:iv%info(metar)%n2))
      allocate (model_psfc(iv%info(metar)%n1:iv%info(metar)%n2))
      ! [1.2] Interpolate horizontally:
      do n=iv%info(metar)%n1, iv%info(metar)%n2
         model_u(1,n)  = jo_grad_y%metar(n)%u
         model_v(1,n)  = jo_grad_y%metar(n)%v
         model_t(1,n)  = jo_grad_y%metar(n)%t
         model_q(1,n)  = jo_grad_y%metar(n)%q
         model_psfc(n) = jo_grad_y%metar(n)%p
      end do
#ifdef A2C
      call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(metar), model_u,'u')
      call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(metar), model_v,'v')
#else
      call da_interp_lin_3d_adj (jo_grad_x%u, iv%info(metar), model_u)
      call da_interp_lin_3d_adj (jo_grad_x%v, iv%info(metar), model_v)
#endif
      call da_interp_lin_3d_adj (jo_grad_x%t, iv%info(metar), model_t)
      call da_interp_lin_3d_adj (jo_grad_x%q, iv%info(metar), model_q)
      call da_interp_lin_2d_adj (jo_grad_x%psfc, iv%info(metar), 1, model_psfc)

      deallocate (model_u)
      deallocate (model_v)
      deallocate (model_t)
      deallocate (model_q)
      deallocate (model_psfc)
   else if (sfc_assi_options == sfc_assi_options_2) then
      call da_transform_xtopsfc_adj(grid,iv,metar,iv%metar(:), jo_grad_y%metar(:),jo_grad_x)
   end if

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_metar_adj")

end subroutine da_transform_xtoy_metar_adj


