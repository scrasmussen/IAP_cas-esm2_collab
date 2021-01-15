#include <define.h>
#if (defined FHNP) && (defined NP)
module pollution

   use precision
   use netcdf
   use spmd
   use colm_varctl, only: N_pollution
   use colm_rtmVar, only: rtmlat, rtmlon 

   implicit none

   character(len=255) :: fnpoint_rtm ="/public/xzh/wangy/N1960-2013_CLM_RTM0.5.nc"

   integer :: nlon, nlat, nyear

   real(r8), pointer :: lon(:)
   real(r8), pointer :: lat(:)
   real(r8), pointer :: year(:)
   real(r8), pointer :: N_point(:,:,:) 

   interface pollution_init
      module procedure pollution_init
   end interface

   interface solute_flux
      module procedure solute_flux
   end interface

   interface pollution_exit
      module procedure pollution_exit
   end interface

   public N_pollution, N_point
   public pollution_init, pollution_exit, solute_flux

contains

   subroutine pollution_init

      integer fid, vid
     
     if (p_master) then
         write(6,*)'Read in N pollution file name: '//trim(fnpoint_rtm)

         call sanity(nf90_open(path=trim(fnpoint_rtm),mode=nf90_nowrite,ncid=fid))

         call sanity(nf90_inq_dimid(fid,'lon',vid))
         call sanity(nf90_inquire_dimension(fid,vid,len=nlon))
         call sanity(nf90_inq_dimid(fid,'lat',vid))
         call sanity(nf90_inquire_dimension(fid,vid,len=nlat))
         call sanity(nf90_inq_dimid(fid,'year',vid))
         call sanity(nf90_inquire_dimension(fid,vid,len=nyear))
      end if

#ifdef SPMD
      call mpi_bcast (nlon ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (nlat ,1,mpi_integer,0,p_comm,p_err)
      call mpi_bcast (nyear,1,mpi_integer,0,p_comm,p_err)
#endif
      
      allocate (lon(nlon))
      allocate (lat(nlat))
      allocate (year(nyear))
      allocate (N_point(nlon,nlat,nyear))

      if (p_master) then
         call sanity(nf90_inq_varid(fid,'lon',vid))
         call sanity(nf90_get_var(fid,vid,lon(:)))

         call sanity(nf90_inq_varid(fid,'lat',vid))
         call sanity(nf90_get_var(fid,vid,lat(:)))

         call sanity(nf90_inq_varid(fid,'year',vid))
         call sanity(nf90_get_var(fid,vid,year(:)))

         call sanity(nf90_inq_varid(fid,'N_point',vid))
         call sanity(nf90_get_var(fid,vid,N_point(:,:,:)))

         call sanity(nf90_close(fid))

         if (nlon.ne.rtmlon .or. nlat.ne.rtmlat) then
            write(6,*) 'fatal error on dimensions of N pollution file', nlon, nlat
            call sanity(nf90_close(fid))
            call abort
         end if
      end if

#ifdef SPMD
      call mpi_bcast (lon,size(lon),mpi_real8,0,p_comm,p_err)
      call mpi_bcast (lat,size(lat),mpi_real8,0,p_comm,p_err)
      call mpi_bcast (year,size(year),mpi_real8,0,p_comm,p_err)
      call mpi_bcast (N_point,size(N_point),mpi_real8,0,p_comm,p_err)
#endif

   end subroutine

   subroutine solute_flux(yr,sltflux)

      integer, intent(in) :: yr
        real(r8), intent(out) :: sltflux(rtmlon,rtmlat)

      integer i, itime

      sltflux(:,:) = 0._r8

      if (N_pollution) then
         itime = -9999

         if (yr.lt.year(1)) then
            itime = -1
         else if (yr.ge.year(nyear)) then
            itime = nyear
         else
            do i = 1, nyear-1
               if(yr.ge.year(i) .and. yr.lt.year(i+1)) then
                  itime = i
                  exit
               end if 
            end do
         end if

         if (itime.lt.0) then
            write(6,*) 'no available N pollution data', nlon, nlat
            call abort
         end if

         sltflux(:,:) = N_point(:,:,itime)
      end if

   end subroutine

   subroutine pollution_exit

      if (.not. N_pollution) return

      deallocate (lon)
      deallocate (lat)
      deallocate (year)
      deallocate (N_point)

   end subroutine

   subroutine sanity(ret)

      integer, intent(in) :: ret

      if (ret .ne. nf90_noerr) then
         write(6,*) trim(nf90_strerror(ret)); stop
      end if
   end subroutine


end module pollution
#endif
