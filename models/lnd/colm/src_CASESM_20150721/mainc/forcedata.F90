!*******************************************************************************
! Initial author: Duoying Ji (2014/08/25)
!*******************************************************************************

#include <define.h>

MODULE forcedata

   use precision
   use phycon_module, only: rgas, grav
   use colm_varctl, only: pco2
   use spmd_decomp, only : gmask, cgmap
   use colm_varMod, only : lats, lonw, latixy, longxy
#if(defined GSWP2 || defined PRINCETON || defined FLUXNET || defined CRUNCEP)
   use metdata, only: metdata_init, metdata_read, metdata_close, forcmask
#endif

   implicit none

   integer            :: lumet              ! logical unit number of meteorological forcing
   character(LEN=255) :: fmet  = ""         ! file name meteorological data
   character(LEN=255) :: fco2  = ""         ! file name global mean annual CO2 concentration data

   interface open_forcedata
      module procedure open_forcedata
   end interface

   interface read_forcedata
      module procedure read_forcedata
   end interface

   interface close_forcedata
      module procedure close_forcedata
   end interface

CONTAINS

   subroutine open_forcedata

      use spmd, only: p_master
      use colm_varMod, only: lon_points, lat_points

      implicit none

       ! Open for meteorological forcing data
#if(defined GSWP2 || defined PRINCETON || defined FLUXNET || defined CRUNCEP)
      CALL metdata_init(lon_points,lat_points,fmet,fco2)
#else
      if (p_master) then
         OPEN(unit=lumet,file=fmet,form='formatted',status='old',action='read')
      end if
#endif

   end subroutine open_forcedata


   subroutine close_forcedata

      use spmd, only : p_master

      implicit none

#if(defined GSWP2 || defined PRINCETON || defined FLUXNET || defined CRUNCEP)
      CALL metdata_close
#else
      if (p_master) then
         CLOSE(lumet)
      end if
#endif

   end subroutine close_forcedata


   subroutine read_forcedata

      use spmd
      use spmd_decomp, only: pgmap, gxmap, gymap
      use timemgr    , only: idate
      use colm_varMod
      implicit none

! arguments:

! local variables:

      integer  :: year             ! current year of model run
      integer  :: jday             ! current julian day of model run 
      integer  :: msec             ! current seconds of model run (0 - 86400)
      real(r8) :: dlat             ! latitude in radians
      real(r8) :: dlon             ! longitude in radians
      real(r8) :: coszen           ! cosine of solar zenith angle
      real(r8) :: calday           ! Julian cal day (1.xx to 365.xx)
      real(r8) :: sunang, cloud, difrat, vnrat
      real(r8) :: orb_coszen

      integer  :: i, j, k, m
      real(r8), pointer :: lonw1d(:)

! routine:

    ! Read in the atmospheric forcing
#if(defined GSWP2 || defined PRINCETON || defined FLUXNET || defined CRUNCEP)
    ! Take cautions against lonw & lats
      lonw1d => lonw(:,1)
      call metdata_read(lon_points,lat_points,gmask,lonw1d,lats,&
                        longxy,latixy,tair,qair,pres,rainc,rainl,&
                        windu,windv,dswrf,dlwrf,tair_z,qair_z,wind_z)
#else
      if (p_master) then
         call getmet(lumet,lon_points,lat_points,tair,qair,pres, &
                     rainc,rainl,windu,windv,dswrf,dlwrf,tair_z,qair_z,wind_z)
      end if

#if(defined SPMD)
      call mpi_bcast(tair,  size(tair),  mpi_real8,0,p_comm,p_err)
      call mpi_bcast(qair,  size(qair),  mpi_real8,0,p_comm,p_err)
      call mpi_bcast(pres,  size(pres),  mpi_real8,0,p_comm,p_err)
      call mpi_bcast(rainc, size(rainc), mpi_real8,0,p_comm,p_err)
      call mpi_bcast(rainl, size(rainl), mpi_real8,0,p_comm,p_err)
      call mpi_bcast(windu, size(windu), mpi_real8,0,p_comm,p_err)
      call mpi_bcast(windv, size(windv), mpi_real8,0,p_comm,p_err)
      call mpi_bcast(dswrf, size(dswrf), mpi_real8,0,p_comm,p_err)
      call mpi_bcast(dlwrf, size(dlwrf), mpi_real8,0,p_comm,p_err)
      call mpi_bcast(tair_z,size(tair_z),mpi_real8,0,p_comm,p_err)
      call mpi_bcast(qair_z,size(qair_z),mpi_real8,0,p_comm,p_err)
      call mpi_bcast(wind_z,size(wind_z),mpi_real8,0,p_comm,p_err)
#endif

#endif

   !* Mapping atmospheric fields to force clm: [lon_points]x[lat_points] grid
   !*     -> [numpatch] vector of subgrid points

      year = idate(1)
      jday = idate(2)
      msec = idate(3)

      calday = float(jday) + float(msec)/86400.

      do k = 1, numcolumn                     !clm vector index

         i = gxmap(cgmap(k))
         j = gymap(cgmap(k))

#ifdef AUTOMASK
         if(forcmask(i,j).eq.1) cycle
#endif

       !---------------------------------------------------------------
       ! as the downward solar is in full band, an empirical expression
       ! will be used to divide fractions of band and incident
       ! (visible, near-infrad, dirct, diffuse)
       ! Julian calday (1.xx to 365.xx)

         dlat = fcon_col(1,k)
         dlon = fcon_col(2,k)

         sunang = orb_coszen(calday,dlon,dlat)

         cloud = (1160.*sunang-dswrf(i,j))/(963.*sunang)
         cloud = max(cloud,0.)
         cloud = min(cloud,1.)
         cloud = max(0.58,cloud)

         difrat = 0.0604/(sunang-0.0223)+0.0683
         if(difrat.lt.0.) difrat = 0.
         if(difrat.gt.1.) difrat = 1.

         difrat = difrat+(1.0-difrat)*cloud
         vnrat = (580.-cloud*464.)/((580.-cloud*499.)+(580.-cloud*464.))
       !---------------------------------------------------------------

         forc%pco2m(k) = pres(i,j)*pco2
         forc%po2m (k) = pres(i,j)*0.209
         forc%us   (k) = windu(i,j)
         forc%vs   (k) = windv(i,j)
         forc%tm   (k) = tair(i,j)
         forc%qm   (k) = qair(i,j)
         forc%prc  (k) = rainc(i,j)
         forc%prl  (k) = rainl(i,j)
         forc%pbot (k) = pres(i,j)
         forc%psrf (k) = pres(i,j)
         forc%sols (k) = (1.0-difrat)*vnrat*dswrf(i,j)
         forc%soll (k) = (1.0-difrat)*(1.0-vnrat)*dswrf(i,j)
         forc%solsd(k) = difrat*vnrat*dswrf(i,j)
         forc%solld(k) = difrat*(1.0-vnrat)*dswrf(i,j)
         forc%frl  (k) = dlwrf(i,j)
         forc%hu   (k) = wind_z(i,j)
         forc%ht   (k) = tair_z(i,j)
         forc%hq   (k) = qair_z(i,j)

#ifndef AUTOMASK
         if(rainc(i,j).lt.0 .or. rainl(i,j).lt.0) then
            write(6,*) 'forcedata.F90', dlat,dlon,rainc(i,j),rainl(i,j), &
                        gmask(i,j),tair(i,j),qair(i,j),pres(i,j),dswrf(i,j),dlwrf(i,j)
         end if
#endif

      end do

   end subroutine read_forcedata

end module forcedata
