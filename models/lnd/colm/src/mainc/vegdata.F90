#include <define.h>

module vegdata

#ifdef VEGDATA

   use precision
   use spmd
   use paramodel, only : numpft
   use timemgr, only : get_curr_date, get_ndays_of_month
   use colm_varMod, only : lon_points, lat_points, numpatch, mlai, msai, mhtop, cvar
   use colm_ioMod, only : flai, lulai

   implicit none

   interface readvegdata
      module procedure readvegdata
   end interface

   interface interp_vegdata
      module procedure interp_vegdata
   end interface

contains 

   subroutine readvegdata

      use spmd_decomp, only : pcmap, pgmap, gxmap, gymap

      real(r8), allocatable :: laixy(:,:,:,:)
      real(r8), allocatable :: saixy(:,:,:,:)
      real(r8), allocatable :: htopxy(:,:,:,:)

      integer c, p, g, x, y, c0, p1, p2, p3

      allocate(laixy (lon_points,lat_points,numpft+1,12))   ! including bare soil
      allocate(saixy (lon_points,lat_points,numpft+1,12))   ! including bare soil
      allocate(htopxy(lon_points,lat_points,numpft+1,12))   ! including bare soil

      if (p_master) then
         open(lulai,file=trim(flai),form='unformatted',status='old')
         read(lulai) laixy
         read(lulai) saixy
         read(lulai) htopxy
         close(lulai)
      end if

#ifdef SPMD
      call mpi_bcast(laixy,size(laixy),mpi_real8,0,p_comm,p_err)
      call mpi_bcast(saixy,size(saixy),mpi_real8,0,p_comm,p_err)
      call mpi_bcast(htopxy,size(htopxy),mpi_real8,0,p_comm,p_err)
#endif

      mlai = 0.
      msai = 0.
      mhtop = 0.

      c0 = 0

      do p = 1, numpatch
         c = pcmap(p)
         g = pgmap(p)
         x = gxmap(g)
         y = gymap(g)

         if(c.ne.c0 .and. cvar%itypwat(c).eq.0) then
            p1 = p
            p2 = p+numpft
            do p3 = p1, p2
               mlai(:,p3) = laixy(x,y,p3-p1+1,:)
               msai(:,p3) = saixy(x,y,p3-p1+1,:)
               mhtop(:,p3) = htopxy(x,y,p3-p1+1,:)
            end do
         end if

         c0 = c
      end do

      deallocate(laixy)
      deallocate(saixy)
      deallocate(htopxy)

   end subroutine readvegdata

#ifdef DUST
   !Vegetation cover for dust emission added by Chenglai Wu: 05/25/2014
   subroutine interp_vegdata(ctype,p1,p2,lai_r,sai_r,htop_r,mvegcov,vegcov_r)
#else
   subroutine interp_vegdata(ctype,p1,p2,lai_r,sai_r,htop_r)
#endif

      integer,  intent(in)  :: ctype
      integer,  intent(in)  :: p1
      integer,  intent(in)  :: p2
      real(r8), intent(out) :: lai_r(numpft+1)
      real(r8), intent(out) :: sai_r(numpft+1)
      real(r8), intent(out) :: htop_r(numpft+1)
#ifdef DUST
      real(r8), intent(in)  :: mvegcov(1:12)
      real(r8), intent(out) :: vegcov_r
#endif
 
      integer p, m1, m2, yr, mon, day, sec, ndays, ndays1, ndays2
      real(r8) r

      if(ctype.ne.0) then ! non-vegetation
         lai_r(:) = 0.
         sai_r(:) = 0.
         htop_r(:) = 0.
#ifdef DUST
         vegcov_r =0. !added by Chenglai Wu: 05/25/2014
#endif
         return
      end if

      if((p2-p1).ne.numpft) stop 'Error on checking numpft in interp_vegdata'

      call get_curr_date(yr,mon,day,sec)

      ndays = get_ndays_of_month(mon,yr)

      if(day.le.ndays/2.) then
         m1 = mon-1
         m2 = mon
         if(m1.lt.1) m1 = 12
      else
         m1 = mon
         m2 = mon+1
         if(m2.gt.12) m2 = 1
      end if

      ndays1 = get_ndays_of_month(m1,yr)
      ndays2 = get_ndays_of_month(m2,yr)

      if(m1.eq.mon) then
         r = (day-ndays1/2.)*2/(ndays1+ndays2)
      else
         r = (day+ndays1/2.)*2/(ndays1+ndays2)
      end if

      do p = p1, p2
         lai_r(p-p1+1) = mlai(m1,p)*(1-r)+mlai(m2,p)*r
         sai_r(p-p1+1) = msai(m1,p)*(1-r)+msai(m2,p)*r
         htop_r(p-p1+1) = mhtop(m1,p)*(1-r)+mhtop(m2,p)*r
      end do
#ifdef DUST
      vegcov_r         = mvegcov (m1)*(1.-r)+mvegcov(m2)*r
#endif

   end subroutine interp_vegdata

#endif

end module vegdata

