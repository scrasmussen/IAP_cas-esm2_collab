#include <define.h>

module colm_rtmMod

#ifdef RTM

   use spmd
   use colm_rtmVar
   implicit none

   interface colm_rtm_init
      module procedure colm_rtm_init
   end interface

   interface colm_rtm_drv
      module procedure colm_rtm_drv
   end interface

   interface colm_rtm_exit
      module procedure colm_rtm_exit
   end interface

CONTAINS

   subroutine colm_rtm_init

      use colm_varctl, only: rtm_dtime
      use colm_varMod, only: numgrid
      use timemgr    , only: get_step_size
      use RtmMod
      use RunoffMod
      implicit none

      integer dtime

    ! Allocate memory to store CoLM flux

      allocate (fevpa(numgrid))
      allocate (rnof (numgrid))
      allocate (prc  (numgrid))
      allocate (prl  (numgrid))

    ! If rtm_nsteps was not entered in the namelist, 
    ! give it the following default value

      if (rtm_nsteps .eq. -9999) then
         dtime = get_step_size()
         rtm_nsteps = rtm_dtime/dtime
      end if

      if (p_master) then
         if (rtm_nsteps > 1) then
            write(6,*)'river runoff calculation performed only every ',rtm_nsteps,' nsteps'
         else
            write(6,*)'river runoff calculation performed every time step'
         endif
      endif

    ! Initialize RTM river routing grid and mask

      call Rtmgridini

    ! Initialize river routing model

      call Rtmlandini

      call Rtmfluxini

   end subroutine colm_rtm_init

   subroutine colm_rtm_drv

      use colm_varMod, only : fldv, numgrid
      use RtmMod

      implicit none

      integer g

      do g = 1, numgrid
         fevpa(g) = fldv%fevpa(g)  ! [mm/s]
          rnof(g) = fldv%rnof(g)   ! [mm/s]
           prc(g) = fldv%prc(g)    ! [mm/s]
           prl(g) = fldv%prl(g)    ! [mm/s]
      end do

      CALL Rtmriverflux()

   end subroutine colm_rtm_drv

   subroutine colm_rtm_exit

      deallocate (fevpa)
      deallocate ( rnof)
      deallocate (  prc)
      deallocate (  prl)

   end subroutine colm_rtm_exit

#endif

end module colm_rtmMod
