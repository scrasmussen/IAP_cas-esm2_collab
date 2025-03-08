!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE SMOOTH_MODULE
!
! This module provides routines for smoothing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module smooth_module

   use parallel_module


   contains

 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: one_two_one
   !
   ! Purpose: Apply the 1-2-1 smoother from the MM5 program TERRAIN 
   !   (found in smth121.F) to array.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine one_two_one(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                          start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array 
  
      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch
  
      allocate(scratch(start_x+1:end_x-1, start_y:end_y, start_z:end_z))
  
      do ipass=1,npass

         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 0.50*array(ix,iy,iz)+0.25*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do
   
         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 0.50*scratch(ix,iy,iz)+0.25*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
             end do
          end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do
  
      deallocate(scratch)
 
   end subroutine one_two_one 
 
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: smth_desmth
   !
   ! Purpose: Apply the smoother-desmoother from the MM5 program TERRAIN 
   !   (found in smther.F) to array.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine smth_desmth(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                          start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval)
 
      implicit none
  
      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array 
  
      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch
  
      allocate(scratch(start_x+1:end_x-1, start_y:end_y, start_z:end_z))
  
      do ipass=1,npass

         !
         ! Smoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 0.5*array(ix,iy,iz) + 0.25*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do
   
         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 0.5*scratch(ix,iy,iz) + 0.25*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)
   
         !
         ! Desmoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 1.52*array(ix,iy,iz) - 0.26*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do
   
         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 1.52*scratch(ix,iy,iz) - 0.26*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do
  
      deallocate(scratch)
 
   end subroutine smth_desmth


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Name: smth_desmth_special
   !
   ! Purpose: Apply the smoother-desmoother from the MM5 program TERRAIN 
   !   (found in smther.F) to array; however, any grid points that were not
   !   originally negative but which have been smoothed to a negative value
   !   will be restored to their original values.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine smth_desmth_special(array, start_dom_x, end_dom_x, start_dom_y, end_dom_y, &
                                  start_x, end_x, start_y, end_y, start_z, end_z, npass, msgval)

      implicit none

      ! Arguments
      integer, intent(in) :: start_dom_x, start_dom_y, start_x, start_y, start_z
      integer, intent(in) :: end_dom_x, end_dom_y, end_x, end_y, end_z
      integer, intent(in) :: npass
      real, intent(in) :: msgval
      real, dimension(start_x:end_x, start_y:end_y, start_z:end_z), intent(inout) :: array

      ! Local variables
      integer :: ix, iy, iz, ipass
      real, pointer, dimension(:,:,:) :: scratch, orig_array

      allocate(scratch(start_x+1:end_x-1, start_y:end_y, start_z:end_z))
      allocate(orig_array(start_x:end_x, start_y:end_y, start_z:end_z))

      orig_array = array

      do ipass=1,npass

         !
         ! Smoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 0.5*array(ix,iy,iz) + 0.25*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do

         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 0.5*scratch(ix,iy,iz) + 0.25*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

         !
         ! Desmoothing pass
         !
         do iy=start_y,end_y
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  scratch(ix,iy,iz) = 1.52*array(ix,iy,iz) - 0.26*(array(ix-1,iy,iz)+array(ix+1,iy,iz))
               end do
            end do
         end do

         do iy=start_y+1,end_y-1
            do ix=start_x+1,end_x-1
               do iz=start_z,end_z
                  array(ix,iy,iz) = 1.52*scratch(ix,iy,iz) - 0.26*(scratch(ix,iy-1,iz)+scratch(ix,iy+1,iz))
               end do
            end do
         end do

         call exchange_halo_r(array, &
                              start_x, end_x, start_y, end_y, start_z, end_z, &
                              start_dom_x, end_dom_x, start_dom_y, end_dom_y, start_z, end_z)

      end do

      ! Remove artificially negative values
      do iy=start_y,end_y
         do ix=start_x,end_x
            do iz=start_z,end_z
               if (array(ix,iy,iz) < 0. .and. orig_array(ix,iy,iz) >= 0.) then
                  array(ix,iy,iz) = orig_array(ix,iy,iz)
               end if
            end do
         end do
      end do

      deallocate(scratch)
      deallocate(orig_array)

   end subroutine smth_desmth_special

   end module smooth_module
