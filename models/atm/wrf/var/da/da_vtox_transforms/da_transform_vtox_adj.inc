subroutine da_transform_vtox_adj(grid, cv_size, xbx, be, ep, vp, vv, cv)

   !--------------------------------------------------------------------------
   ! Purpose: Control variable transform Adjoint
   !--------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   integer,        intent(in)    :: cv_size ! Size of cv array.
   type(xbx_type), intent(in)    :: xbx  ! For header & non-grid arrays.
   type(be_type),  intent(in)    :: be   ! background errors.
   type(ep_type),  intent(in)    :: ep   ! ensemble perturbation structure.
   type(vp_type),  intent(out)   :: vp   ! grdipt/level cv (local).
   type(vp_type),  intent(out)   :: vv   ! grdipt/eof cv (local).
   real,           intent(inout) :: cv(1:cv_size) ! control variables.

   if (trace_use) call da_trace_entry("da_transform_vtox_adj")

   if (.not. use_background_errors) return

   if (  cv_options == 3 ) then

      call da_transform_bal_adj( vp, be, grid)

      call da_apply_be_adj( be, cv, vp, grid)

   else
  
   !-------------------------------------------------------------------------
   ! [3.0]: Perform x = u_p vp transform::
   !-------------------------------------------------------------------------

   call da_zero_vp_type (vp)
   call da_transform_vptox_adj(grid, vp, be, ep)

   !-------------------------------------------------------------------------
   ! [2.0]: Perform vp = u_v vv transform:
   !-------------------------------------------------------------------------

   call da_zero_vp_type (vv)

   if (vert_corr == 2) then      
      call da_vertical_transform(grid, 'u_adj', be, &
         grid%xb % vertical_inner_product, vv, vp)
   else
      vv % v1(its:ite,jts:jte,kts:kte) = vp % v1(its:ite,jts:jte,kts:kte)
      vv % v2(its:ite,jts:jte,kts:kte) = vp % v2(its:ite,jts:jte,kts:kte)
      vv % v3(its:ite,jts:jte,kts:kte) = vp % v3(its:ite,jts:jte,kts:kte)
      vv % v4(its:ite,jts:jte,kts:kte) = vp % v4(its:ite,jts:jte,kts:kte)
      vv % v5(its:ite,jts:jte,kts:kte) = vp % v5(its:ite,jts:jte,kts:kte)

   !  Uv for alpha is a null transform:
      if (be % ne > 0) then
         vv % alpha(its:ite,jts:jte,kts:kte,1:be%ne) = &
         vp % alpha(its:ite,jts:jte,kts:kte,1:be%ne) 
      end if
   end if
   
   end if

   !-------------------------------------------------------------------------
   ! [1.0]: perform vv = u_h cv transform:
   !-------------------------------------------------------------------------

   if (global) then
      call da_transform_vtovv_global_adj(cv_size, xbx, be, cv, vv)
   else if ( (fg_format == fg_format_wrf_arw_regional .or. &
              fg_format == fg_format_wrf_nmm_regional) .and. &
              (.not. cv_options == 3) ) then
      call da_transform_vtovv_adj(grid, cv_size, be, cv, vv)
   end if

   if (trace_use) call da_trace_exit("da_transform_vtox_adj")

end subroutine da_transform_vtox_adj


