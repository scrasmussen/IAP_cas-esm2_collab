subroutine da_deallocate_y(y)

   !---------------------------------------------------------------------------
   ! Purpose: Deallocate arrays used in y and residual obs structures.
   !
   ! Method:  Deallocate component in turn.
   !---------------------------------------------------------------------------

   implicit none
   
   type (y_type), intent(inout)          :: y      ! residual type structure.
   integer                               :: n,i  ! Loop counter.


   if (trace_use) call da_trace_entry("da_deallocate_y")

   !---------------------------------------------------------------------------
   ! [1.0] Deallocate:
   !---------------------------------------------------------------------------

   if (y % nlocal(synop) > 0) deallocate (y % synop)

   if (y % nlocal(ships) > 0) deallocate (y % ships)

   if (y % nlocal(metar) > 0) deallocate (y % metar)


   if (y % nlocal(sound) > 0) then
      do n = 1, y % nlocal(sound)
         deallocate (y % sound(n)%u)
         deallocate (y % sound(n)%v)
         deallocate (y % sound(n)%t)
         deallocate (y % sound(n)%q)
      end do

      deallocate (y % sound)

      deallocate (y % sonde_sfc)
   end if
      
   if (y % nlocal(mtgirs) > 0) then
      do n = 1, y % nlocal(mtgirs)
         deallocate (y % mtgirs(n)%u)
         deallocate (y % mtgirs(n)%v)
         deallocate (y % mtgirs(n)%t)
         deallocate (y % mtgirs(n)%q)
      end do

      deallocate (y % mtgirs)

   end if

   if (y % nlocal(tamdar) > 0) then
      do n = 1, y % nlocal(tamdar)
         deallocate (y % tamdar(n)%u)
         deallocate (y % tamdar(n)%v)
         deallocate (y % tamdar(n)%t)
         deallocate (y % tamdar(n)%q)
      end do

      deallocate (y % tamdar)
      deallocate (y % tamdar_sfc)
   end if


   if (y % nlocal(pilot) > 0) then
      do n = 1, y % nlocal(pilot)
         deallocate (y % pilot(n)%u)
         deallocate (y % pilot(n)%v)
      end do
      deallocate (y % pilot)
   end if

   if (y % nlocal(bogus) > 0) then
      do n = 1, y % nlocal(bogus)
         deallocate (y % bogus(n)%u)
         deallocate (y % bogus(n)%v)
         deallocate (y % bogus(n)%t)
         deallocate (y % bogus(n)%q)
      end do
      deallocate (y % bogus)
   end if

    if (y % nlocal(radar) > 0) then
       do n = 1, y % nlocal(radar)
          deallocate (y % radar(n)%rv)
          deallocate (y % radar(n)%rf)
       end do
       deallocate (y % radar)
    end if


   if (y % nlocal(airep) > 0) then
      do n = 1, y % nlocal(airep)
         deallocate (y % airep(n)%u)
         deallocate (y % airep(n)%v)
         deallocate (y % airep(n)%t)
      end do
      deallocate (y % airep)
   end if

   if (y % nlocal(geoamv) > 0) then
      do n=1, y % nlocal(geoamv)
         deallocate (y % geoamv(n) % u)
         deallocate (y % geoamv(n) % v)
      end do
      deallocate (y % geoamv)
   end if

   if (y % nlocal(polaramv) > 0) then
      do n=1, y % nlocal(polaramv)
         deallocate (y % polaramv(n) % u)
         deallocate (y % polaramv(n) % v)
      end do
      deallocate (y % polaramv)
   end if

   if (y % nlocal(gpspw) > 0) deallocate (y % gpspw)

   if (y % nlocal(gpsref) > 0) then
      do n = 1, y % nlocal(gpsref)
         deallocate (y % gpsref(n)%ref)
         deallocate (y % gpsref(n)%  p)
         deallocate (y % gpsref(n)%  t)
         deallocate (y % gpsref(n)%  q)
      end do
      deallocate (y % gpsref)
   end if

   if (y % nlocal(satem) > 0) then
      do n = 1, y % nlocal(satem)
         deallocate (y % satem(n) % thickness)
      end do
      deallocate (y % satem)
   end if

   if (y % nlocal(ssmi_tb) > 0) deallocate (y % ssmi_tb)
   if (y % nlocal(ssmi_rv) > 0) deallocate (y % ssmi_rv)
   if (y % nlocal(pseudo)  > 0) deallocate (y % pseudo)

   if (y % nlocal(ssmt1) > 0) then
      do n = 1, y % nlocal(ssmt1)
         deallocate (y % ssmt1(n) % t)
      end do
      deallocate (y % ssmt1)
   end if

   if (y % nlocal(ssmt2) > 0) then
      do n = 1, y % nlocal(ssmt2)
         deallocate (y % ssmt2(n) % rh)
      end do
      deallocate (y % ssmt2)
   end if

   if (y % nlocal(qscat) > 0) deallocate (y % qscat)

   if (y % nlocal(profiler) > 0) then
      do n = 1, y % nlocal(profiler)
         deallocate (y % profiler(n)%u)
         deallocate (y % profiler(n)%v)
      end do
      deallocate (y % profiler)
   end if

   if (y % nlocal(buoy)  > 0) deallocate (y % buoy)

   !  radiance:
 
   if (y % num_inst > 0) then
      do i = 1,  y % num_inst
        if (y % instid(i) % num_rad < 1) cycle
        ! deallocate (y % instid(i) % ichan)
        deallocate ( y % instid(i) % tb )
      end do
      deallocate (y % instid)
   end if
   if (y % nlocal(airsr) > 0) then
      do n = 1, y % nlocal(airsr)
         deallocate (y % airsr(n)%t)
         deallocate (y % airsr(n)%q)
      end do
      deallocate (y % airsr)
   end if

   if (trace_use) call da_trace_exit("da_deallocate_y")

end subroutine da_deallocate_y


