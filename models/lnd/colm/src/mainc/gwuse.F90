!*******************************************************************************
! Initial author: Duoying Ji (2014/08/25)
! wanglh (2019/01/10)
!*******************************************************************************

#include <define.h>
#if (defined FHNP) && (defined GWUSE)
MODULE gwdata

   use precision
   use netcdf
   use phycon_module
   use colm_varctl
   use spmd_decomp
   use colm_varMod, only : lats, lonw, latixy, longxy, lon_points, lat_points
   use ncdio 


   implicit none

   character(LEN=255) :: fgw  ='/public/xzh/wanglh/mydata/global_work/final_Global_gw_use_f14_f14_new3.nc'        ! file name meteorological data
   integer :: nlon, nlat, ntime


   real(r8), pointer :: lon_gw(:)
   real(r8), pointer :: lat_gw(:)
   real(r8), pointer :: ind_hum_frac(:,:)
   real(r8), pointer :: gw_intake(:,:,:,:) ! 
   real(r8), pointer :: gw_uptake(:,:,:)
   real(r8), pointer :: noirr_frac(:)
   interface open_gwdata
      module procedure open_gwdata
   end interface

   interface read_gwdata
      module procedure read_gwdata
   end interface

   interface get_cgwdata
      module procedure get_cgwdata
   end interface

   interface close_gwdata
      module procedure close_gwdata
   end interface

CONTAINS

   subroutine open_gwdata

      use spmd, only: p_master
      use colm_varMod, only: lon_points, lat_points
      
      implicit none
     integer fid, vid
     integer :: nlon, nlat
     character(len=20):: name1,name2
  

      if (p_master) then
         print *, 'open groundwater use  file: '//trim(fgw)

         call sanity(nf_open(trim(fgw),nf_nowrite,fid))

         call sanity(nf_inq_dimid(fid,'lon',vid))
         call sanity(nf_inq_dim(fid,vid,name1,nlon))
         call sanity(nf_inq_dimid(fid,'lat',vid))
         call sanity(nf_inq_dim(fid,vid,name2,nlat))

         write(6,*) 'nlon=',nlon
                  write(6,*) 'nlat=',nlat
         if (nlon.ne.lon_points .or. nlat.ne.lat_points) then
            write(6,*) 'fatal error on dimensions of groundwater use  file', nlon, nlat, lon_points,lat_points
            call sanity(nf_close(fid))
            call abort
         end if 
         call sanity(nf_close(fid))
      end if
  
   end subroutine open_gwdata


   subroutine close_gwdata

      use spmd, only : p_master

      implicit none

     

      deallocate (lon_gw)
      deallocate (lat_gw)
      deallocate (ind_hum_frac)
      deallocate (gw_intake)

   end subroutine close_gwdata


   subroutine read_gwdata

      use spmd
      use spmd_decomp, only: pgmap, gxmap, gymap,ggmap,gxmap_glob, gymap_glob,cgmap
      use timemgr, only : get_step_size, get_curr_date!wanglh
      !use colm_varMod, only : latixy, longxy, numcolumn
      use colm_varMod
      implicit none

! arguments:

! local variables:
   real(r8), pointer ::  lon_gw(:)
   real(r8), pointer ::  lat_gw(:)
   real(r8), pointer :: ind_hum_frac(:,:)
   real(r8), pointer :: gw_intake(:,:,:,:) 
   real(r8), pointer :: gw_uptake(:,:,:)
   real(r8), pointer :: noirr_frac(:)
   real(r8), pointer :: cgwintake(:)
   real(r8), pointer :: cgw(:)
   real(r8), pointer :: llatixy(:,:)
   real(r8), pointer :: llongxy(:,:)
      integer  :: step_size,yrp1,monp1,dayp1,secp1
      integer  :: i, j, k, m, numlat1,numlon1
      real(r8) :: er=1.e-1_r8
      real(r8) :: ggw_temp
      integer :: myrank,g
      integer fid1, vid1,vid2,vid3,vid4, nlon, nlat,begg,endg
! routine:
      nlon=256
      nlat=128
      allocate (lon_gw(nlon))
      allocate (lat_gw(nlat))
      allocate (ind_hum_frac(nlon,nlat))
      allocate (gw_intake(nlon,nlat,12,46))
!gw_uptake=>fldv%gw_uptake
!noirr_frac=>fldv%noirr_frac
llatixy=>latixy
llongxy=>longxy 
 ggw_temp=0._r8
      if (p_master) then
          write(6,*)'now p_master begin get the grounwater data '
          write(6,*),fgw
         call sanity(nf_open(trim(fgw),nf_nowrite,fid1))
          write(6,*)'now p_master begin get lon '
         call sanity(nf_inq_varid(fid1,'lon',vid1))
         call sanity(nf_get_var_double(fid1,vid1,lon_gw))
         !write(6,*),lon_gw
         write(6,*)'now p_master begin get lat '
         call sanity(nf_inq_varid(fid1,'lat',vid2))
         call sanity(nf_get_var_double(fid1,vid2,lat_gw))
         !write(6,*),lat_gw
         write(6,*)'now p_master begin get frac ' 
         call sanity(nf_inq_varid(fid1,'frac',vid3))
         call sanity(nf_get_var_double(fid1,vid3,ind_hum_frac))
         write(6,*)'now p_master begin get data '
         call sanity(nf_inq_varid(fid1,'datag',vid4))
         call sanity(nf_get_var_double(fid1,vid4,gw_intake))
         write(6,*)'now p_master end read data '
         call sanity(nf_close(fid1))
      end if


      if(p_master) write(6,*)'now p_master begin bcast data '
               
     ! call mpi_barrier(p_comm,p_err)

      !call mpi_bcast (real(lon_gw,8),size(lon_gw),mpi_real8,0,p_comm,p_err)
      !if(p_master) write(6,*)'1 '
      !call mpi_bcast (real(lat_gw,8),size(lat_gw),mpi_real8,0,p_comm,p_err)
      !if(p_master) write(6,*)'2 '
#ifdef SPMD
      call mpi_bcast (lon_gw,size(lon_gw),mpi_real8,0,p_comm,p_err)
      if(p_master) write(6,*)'1 '
      call mpi_bcast (lat_gw,size(lat_gw),mpi_real8,0,p_comm,p_err)
      if(p_master) write(6,*)'2 '
      call mpi_bcast (ind_hum_frac,size(ind_hum_frac),mpi_real8,0,p_comm,p_err)
      if(p_master) write(6,*)'3 ' 
      call mpi_bcast (gw_intake,size(gw_intake),mpi_real8,0,p_comm,p_err)
      !call mpi_barrier(p_comm,p_err)
#endif
         if(p_master)       write(6,*)'now p_master end bcast data '



      

   !* Mapping atmospheric fields to force clm: [lon_points]x[lat_points] grid
   !*     -> [numpatch] vector of subgrid points
      if (p_master) then
          write(6,*)'now p_master begin print size of lat_gw ',size(lat_gw)
          write(6,*)'now p_master begin print size of lon_gw ',size(lon_gw)
          !write(6,*)'now p_master begin print  lat_gw ',lat_gw
          !write(6,*)'now p_master begin print  lon_gw ',lon_gw
          !write(6,*),size(lat_gw)
          !write(6,*),size(lon_gw)
          write(6,*)'now p_master begin print size of ind_hum_frac ',size(ind_hum_frac)
          write(6,*)'now p_master begin print size of gw_intake ',size(gw_intake)
          write(6,*),ind_hum_frac(1,1)
          write(6,*),gw_intake(186,92,1,1)
         !  write(6,*),size(gw_uptake)
                      write(6,*),size(llongxy)
                                 write(6,*),size(llatixy)
            ! write(6,*),gw_intake(256,128,12,46)                   !  write(6,*),size(noirr_frac)
      end if
!wanlh map to grid
        

#ifdef SPMD
               call mpi_comm_rank(p_comm,myrank,p_err)
#endif

      do k = 1, numgrid                  !clm vector index

               

         i = gxmap(k)
         j = gymap(k)

                       
                    
                        if((gw_intake(i,j,1,1).ge.0_r8).and.(gw_intake(i,j,1,1).le.5.e-5_r8))then            


                              fldv%gw_uptake(k,:,:)=gw_intake(i,j,:,:)
                           
                              fldv%noirr_frac(k)=ind_hum_frac(i,j)
                  !      !      fldv%ggw(k)=gw_intake(numlon1,numlat1,1,1)
                   !        write(6,*)'myrank=',myrank,' k=',k
                   !        write(6,*) 'glat=',llatixy(i,j),'glon=',llongxy(i,j),&
                   !           'flat=',lat_gw(j),'flon=',lon_gw(i),&
                   !              'gwuptake=',fldv%gw_uptake(k,1,15),'noirr_frac=',fldv%noirr_frac(k)
                                 
                        else
                              fldv%gw_uptake(k,:,:)=0.
                              fldv%noirr_frac(k)=0.
                         end if


      end do
   

!if(myrank.eq.10) write(6,*)'fldv%ggw(3) is',fldv%ggw(3),'fldv%gw_uptake(k,1,1)',fldv%gw_uptake(3,1,1)




    !  deallocate (lon_gw)
    !  deallocate (lat_gw)
    !  deallocate (ind_hum_frac)
    !  deallocate (gw_intake)




if(p_master) write(6,*)'now end read gwdata,I am coming '
   end subroutine read_gwdata

! grid to column 
   subroutine get_cgwdata

      use spmd
      use spmd_decomp, only: pgmap, gxmap, gymap,ggmap,cgmap,gxmap_glob, gymap_glob
      use timemgr, only : get_step_size, get_curr_date!wanglh
      !use colm_varMod, only : latixy, longxy, numcolumn
      use colm_varMod
      implicit none

! arguments:

! local variables:
   real(r8), pointer :: ind_hum_frac(:,:)
   real(r8), pointer :: gw_intake(:,:,:,:) 
   real(r8), pointer :: gw_uptake(:,:,:)
   real(r8), pointer :: noirr_frac(:)
   real(r8), pointer :: cgwintake(:)
   real(r8), pointer :: cgw(:)
   real(r8), pointer :: llatixy(:,:)
   real(r8), pointer :: llongxy(:,:)
      integer  :: yrp1,monp1,dayp1,secp1
      integer  :: i, j, k, m, numlat1,numlon1
      real(r8) :: er=1.e-2_r8
      integer :: myrank,g,g1
      integer fid, vid,countt
! routine:
countt=0

gw_uptake=>fldv%gw_uptake
!cgwintake=>fldv_col%cgwintake
llatixy=>latixy
llongxy=>longxy 


!wanglh map to column 
#ifdef SPMD
  call mpi_comm_rank(p_comm,myrank,p_err)
#endif
   call get_curr_date(yrp1, monp1, dayp1, secp1)
        do k = 1, numcolumn   
          fldv_col%cgwintake(k) = 0._r8
          
           g = cgmap(k)
                        g1=ggmap(g)

                        i = gxmap_glob(g1)
                        j = gymap_glob(g1)

            if ((yrp1 >= 1965 .and. yrp1 <= 2010) .and. (cvar%itypwat(k)<=1)) then
     !      ! if ((yrp1 >= 1965 .and. yrp1 <= 2010) ) then
                 fldv_col%cgwintake(k) = fldv%gw_uptake(g,monp1,yrp1-1964)
      !               if((abs(llatixy(i,j)-44.6456).lt.0.01).and.(abs(llongxy(i,j)-285.4687).lt.0.01)) write(6,*) 'CGW :yrp1',yrp1,'fldv%gw_uptake is :',fldv%gw_uptake(g,monp1,yrp1-1964),'fldv_col%cgwintake is :',fldv_col%cgwintake(k)
                     
             end if

                  if ((yrp1 >= 1948 .and. yrp1 < 1965) .and.(cvar%itypwat(k)<=1))then
       !        ! if ((yrp1 >= 1948 .and. yrp1 < 1965) )then
       !             !if(p_master) write(6,*)'yrp1 is :',yrp1
                    
                     fldv_col%cgwintake(k) = fldv%gw_uptake(g,monp1,1)
        !            if((abs(llatixy(i,j)-44.6456).lt.0.01).and.(abs(llongxy(i,j)-285.4687).lt.0.01)) write(6,*) 'CGW :yrp1',yrp1,'fldv%gw_uptake is :',fldv%gw_uptake(g,monp1,1),'fldv_col%cgwintake is :',fldv_col%cgwintake(k)
                    
                     !fldv_col%cgwintake(k) = 0.0000025
                 
                  end if


          

if (cvar%itypwat(k)<=1) then
countt=countt+1
end if






         end do

!write(6,*)'myrank is ',myrank,'countt is',countt


        !       write(6,*)'now p_master begin bcast data '


   !   call mpi_bcast (cgwintake,size(cgwintake),mpi_real8,0,p_comm,p_err)

!
         !      write(6,*)'now p_master end bcast data '


   end subroutine get_cgwdata


   subroutine sanity(ret)

      integer, intent(in) :: ret

      if (ret .ne. nf_noerr) then
         write(6,*) trim(nf_strerror(ret)); stop
      end if
   end subroutine


end module gwdata
#endif 
