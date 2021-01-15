!----------------------------------------------------------------------------
!    This program is designed to use PIO output GEATM results.
!    Juanxiong He, 2013-06-30
!----------------------------------------------------------------------------

  module inout3
  use pio, only : pio_write_darray, io_desc_t, file_desc_t, var_desc_t, iosystem_desc_t, &
      pio_int, pio_real, pio_double, pio_noerr, &
      PIO_iotype_binary, pio_iotype_netcdf, &
      pio_char, pio_write, pio_clobber,pio_noclobber,&
      pio_initdecomp,  pio_openfile, pio_closefile, pio_createfile,pio_freedecomp, &
      pio_finalize,  PIO_enddef,  PIO_def_dim,  PIO_def_var, PIO_put_att, PIO_put_var
  use geatm_vartype, only : gc_name,gc_unit
  type(iosystem_desc_t), pointer, private :: pio_subsystem
  contains        

  subroutine output_atm(iyear, imonth, iday, ihour, ne)

   use geatm_vartype, only: sx,ex,sy,ey,nx,ny,nz,iaer,igas,isize,nzz,ntt
   use naqpms_varlist,only: ip2mem,ip3mem,ip4mem,ip5mem, &
       duo3, duno2, ana, aso4, anh4, ano3, &
       acl, cph, ope, aer, gas, duso2, &
       jo1d, jno2, uvb, uvbs, uva, vis,&
       ssa, aod, cldopd
   use met_fields,only: u,v,t,w,rh1,plev,h,rhsfc,u10,v10,ust,t2,psfc
   use naqpms_gridinfo, only: heiz,terrain,dz
   use work_vars,only: tropp 

   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, km, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   real :: xlat(ny(ne)),xlon(nx(ne)),sigma(nzz),b(isize),c(iaer),fill_value
   real, dimension(:,:),allocatable :: chem2d
   real, dimension(:,:,:),allocatable :: chem3d
   integer, pointer :: ldof2d(:),ldof3d(:)
   integer         :: m,mlen, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48) 
   fillvalue = -99999.9
   
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo

   allocate(varid_x2a(18))
   
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))
   
   fname='geatm.atm.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
  
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m)) 
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)   
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   dim2d(1:2)=dim3d(1:2)
   rcode = pio_def_dim(File,'lev',nzz,dim3d(3))
   
   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east') 
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'dz',PIO_REAL,dim3d,varid_x2a(4))   
   rcode = pio_def_var(File,'u',PIO_REAL,dim3d,varid_x2a(5))   
   rcode = pio_def_var(File,'v',PIO_REAL,dim3d,varid_x2a(6))   
   rcode = pio_def_var(File,'w',PIO_REAL,dim3d,varid_x2a(7))
   rcode = pio_def_var(File,'t',PIO_REAL,dim3d,varid_x2a(8))   
   rcode = pio_def_var(File,'rh1',PIO_REAL,dim3d,varid_x2a(9))
   rcode = pio_def_var(File,'plev',PIO_REAL,dim3d,varid_x2a(10))   
   rcode = pio_def_var(File,'heiz',PIO_REAL,dim3d,varid_x2a(11))
   rcode = pio_def_var(File,'tropp',PIO_REAL,dim2d,varid_x2a(12))
   rcode = pio_def_var(File,'rhsfc',PIO_REAL,dim2d,varid_x2a(13))
   rcode = pio_def_var(File,'u10',PIO_REAL,dim2d,varid_x2a(14))
   rcode = pio_def_var(File,'v10',PIO_REAL,dim2d,varid_x2a(15))
   rcode = pio_def_var(File,'ust',PIO_REAL,dim2d,varid_x2a(16))
   rcode = pio_def_var(File,'t2',PIO_REAL,dim2d,varid_x2a(17))
   rcode = pio_def_var(File,'psfc',PIO_REAL,dim2d,varid_x2a(18))

   mm=18
   do i=1,mm
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   
   ! dz
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=dz(km) 
          end do
          end do
   end do  
   call pio_write_darray(File, varid_x2a(4), iodesc3d,chem3d, rcode)
   ! u
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=u(km)                        
          end do
          end do
   end do     
   call pio_write_darray(File, varid_x2a(5), iodesc3d, chem3d, rcode)
   ! v
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=v(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(6), iodesc3d, chem3d, rcode)
   ! w
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=w(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(7), iodesc3d, chem3d, rcode)
   ! t
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=t(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(8), iodesc3d, chem3d, rcode)
   ! rh
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=rh1(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(9), iodesc3d, chem3d, rcode)
   ! p
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=plev(km)                        
          end do
          end do
   end do        
   call pio_write_darray(File, varid_x2a(10), iodesc3d, chem3d, rcode)
   ! z
   do k=1,nz(ne)
        i0=ip3mem(k,ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=h(km) 
          end do
          end do
   end do  
   call pio_write_darray(File, varid_x2a(11), iodesc3d, chem3d, rcode)
   ! tropp   
        i0=ip2mem(ne)
          do j= sy(ne),ey(ne)
          do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=tropp(km)                        
          end do
          end do   
   call pio_write_darray(File, varid_x2a(12), iodesc2d, chem2d, rcode)
   ! rhsfc
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=rhsfc(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(13), iodesc2d, chem2d, rcode)
   ! u10
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=u10(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(14), iodesc2d, chem2d, rcode)
   ! v10
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=v10(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(15), iodesc2d, chem2d, rcode)
   ! ust
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=ust(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(16), iodesc2d, chem2d, rcode)
   ! t2
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=t2(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(17), iodesc2d, chem2d, rcode)
   ! psfc
        i0=ip2mem(ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=psfc(km)
        end do
        end do
   call pio_write_darray(File, varid_x2a(18), iodesc2d, chem2d, rcode)



   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     

   deallocate(ldof2d)
   deallocate(ldof3d)
   deallocate(chem2d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_atm
!--------------------------------------------------------------------------------
  subroutine output_aerosol(iyear, imonth, iday, ihour, it1, ne)
   use geatm_vartype, only: sx,ex,sy,ey,nx,ny,nz,iaer,igas,isize,nzz,ntt
   use naqpms_varlist,only: ip2mem,ip3mem,ip4mem,ip5mem,aerom,iedgas,ip4mem_aer
   use met_fields,only: u,v,t,w,rh1,plev,h,rhsfc,u10,v10,ust,t2,psfc
   use naqpms_gridinfo, only: heiz,terrain,dz
   use work_vars,only: tropp 

   integer :: iyear, imonth, iday, ihour, it1,  ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, km, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   real :: xlat(ny(ne)),xlon(nx(ne)),sigma(nzz),b(isize),c(iaer),fill_value
   real, dimension(:,:),allocatable :: chem2d
   real, dimension(:,:,:),allocatable :: chem3d
   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer, pointer :: ldof2d(:),ldof3d(:),ldof5d(:)
   integer         :: m,mlen, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48) 
   fillvalue = -99999.9
   
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo

   do i=1,isize
   b(i)=i
   enddo
   do i=1,iaer
   c(i)=i
   enddo
      
   allocate(varid_x2a(10))
   
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))
   
   fname='geatm.aerosol.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
  
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m)) 
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)   
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   dim2d(1:2)=dim3d(1:2)
   rcode = pio_def_dim(File,'lev',nzz,dim3d(3))
   dim5d(1:3)=dim3d(1:3)
   
   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east') 
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))

   rcode = pio_def_var(File,'ppm',PIO_REAL,dim3d,varid_x2a(4))
   rcode = pio_put_att(File,varid_x2a(4),"_fillvalue",fillvalue)
   rcode = pio_def_var(File,'bc',PIO_REAL,dim3d,varid_x2a(5))
   rcode = pio_put_att(File,varid_x2a(5),"_fillvalue",fillvalue)
   rcode = pio_def_var(File,'oc',PIO_REAL,dim3d,varid_x2a(6))
   rcode = pio_put_att(File,varid_x2a(6),"_fillvalue",fillvalue)
   rcode = pio_def_var(File,'so4',PIO_REAL,dim3d,varid_x2a(7))
   rcode = pio_put_att(File,varid_x2a(7),"_fillvalue",fillvalue)
   rcode = pio_def_var(File,'nh4',PIO_REAL,dim3d,varid_x2a(8))
   rcode = pio_put_att(File,varid_x2a(8),"_fillvalue",fillvalue)
   rcode = pio_def_var(File,'no3',PIO_REAL,dim3d,varid_x2a(9))
   rcode = pio_put_att(File,varid_x2a(9),"_fillvalue",fillvalue)
   rcode = pio_def_var(File,'pm10',PIO_REAL,dim3d,varid_x2a(10))
   rcode = pio_put_att(File,varid_x2a(10),"_fillvalue",fillvalue)

   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   
   do ia=1,6
   do is=1,1
   do k=1,nzz
     i0=ip4mem_aer(k,is,ia,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
           km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
           chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aerom(km)
        end do
        end do
   enddo
   call pio_write_darray(File, varid_x2a(3+ia), iodesc3d, chem3d, rcode)   
   enddo
   enddo

   do ia=1,1
   do is=2,2
   do k=1,nzz
     i0=ip4mem_aer(k,is,ia,ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
           km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
           chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aerom(km)
        end do
        end do
   enddo
   call pio_write_darray(File, varid_x2a(10), iodesc3d, chem3d, rcode)  
   enddo
   enddo


   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     

   deallocate(ldof2d)
   deallocate(ldof3d)
   deallocate(chem2d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_aerosol
!--------------------------------------------------------------------------------
  subroutine output_gas(iyear, imonth, iday, ihour, it1, ne)
   use geatm_vartype, only: sx,ex,sy,ey,nx,ny,nz,iaer,igas,&
       isize,nzz,ntt,PrintGas
   use naqpms_varlist,only: ip2mem,ip3mem,ip4mem,ip5mem,iedgas,gas
   use met_fields,only: u,v,t,w,rh1,plev,h,rhsfc,u10,v10,ust,t2,psfc
   use naqpms_gridinfo, only: heiz,terrain,dz
   use work_vars,only: tropp

   integer :: iyear, imonth, iday, ihour, it1,  ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc3d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, km, dim3d(3)
   real :: xlat(ny(ne)),xlon(nx(ne)),sigma(nzz),fill_value
   real, dimension(:,:,:),allocatable :: chem3d
   integer, pointer :: ldof3d(:)
   integer         :: m,mlen, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
   character*3 cgas(iedgas)
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48) 
   do i=1,iedgas
   cgas(i)=char(i/100+48)//char((i-i/100*100)/10+48)//char(i-i/10*10+48)
   end do
   fillvalue = -99999.9
   
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nz(ne)
   sigma(i)=i
   enddo

   mlen=0
     do i=1,iedgas
     if(PrintGas(i)==1)then
     mlen=mlen+1
     endif
     enddo

   allocate(varid_x2a(3+mlen))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))
   
   fname='geatm.gas.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
  
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m))
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   rcode = pio_def_dim(File,'lev',nz(ne),dim3d(3))
   
   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east') 
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))


   do ig=1,iedgas
     if(PrintGas(ig)==1)then
     if(ig.eq.1) rcode = pio_def_var(File,'h2so4',PIO_REAL,dim3d,varid_x2a(3+ig))
     if(ig.eq.2) rcode = pio_def_var(File,'h2o2',PIO_REAL,dim3d,varid_x2a(3+ig))
     if(ig.eq.3) rcode = pio_def_var(File,'so2',PIO_REAL,dim3d,varid_x2a(3+ig))
     if(ig.eq.4) rcode = pio_def_var(File,'dms',PIO_REAL,dim3d,varid_x2a(3+ig))
     if(ig.eq.5) rcode = pio_def_var(File,'nh3',PIO_REAL,dim3d,varid_x2a(3+ig))
     endif
   end do
  
   do i=1,3+iedgas
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   
   ! gas
   mlen=0
   do ig=1,iedgas        
     if(PrintGas(ig)==1)then        
     mlen=mlen+1    
     do k=1,nz(ne)
        i0=ip4mem(k,ig,ne)   
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
           km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
           chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=gas(km)
        end do
        end do        
     end do
     call pio_write_darray(File, varid_x2a(3+mlen), iodesc3d, chem3d, rcode)
     endif
   end do

   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     

   deallocate(ldof3d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_gas
!--------------------------------------------------------------------

  subroutine output_dust(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only: &
      sx,ex,sy,ey,nx,ny,nz,iaer,igas,isize,nzz,ntt
   use naqpms_varlist, only: ip2mem,ip3mem,ip5mem,aer      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i0, ig, is, ia, dim2d(2), dim3d(3)
   integer :: mem5d, mem3d, mem2d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nz(ne)), fill_value
   real, dimension(:,:,:),allocatable :: chem3d
   integer, pointer :: ldof2d(:),ldof3d(:)
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname

   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)

   fillvalue = -99999.9
   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo
      
   allocate(varid_x2a(7))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))   
   
   fname='geatm.dust.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m)) 
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   dim2d(1:2)=dim3d(1:2)
   rcode = pio_def_dim(File,'lev',nz(ne),dim3d(3))  

   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'dust01',PIO_REAL,dim3d,varid_x2a(4))
   rcode = pio_def_var(File,'dust02',PIO_REAL,dim3d,varid_x2a(5))
   rcode = pio_def_var(File,'dust03',PIO_REAL,dim3d,varid_x2a(6))
   rcode = pio_def_var(File,'dust04',PIO_REAL,dim3d,varid_x2a(7))
   mm=7
   
   do i=1,mm
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
          
   
      do k=1,nzz        
        mem5d=ip5mem(k,1,2, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)  
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(4), iodesc3d, chem3d, rcode)

      do k=1,nzz
        mem5d=ip5mem(k,2,2, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(5), iodesc3d, chem3d, rcode)

      do k=1,nzz
        mem5d=ip5mem(k,3,2, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(6), iodesc3d, chem3d, rcode)

      do k=1,nzz
        mem5d=ip5mem(k,4,2, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(7), iodesc3d, chem3d, rcode)

   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     

   deallocate(ldof2d)
   deallocate(ldof3d)
   deallocate(chem3d)   
   deallocate(varid_x2a)
   end subroutine output_dust
!--------------------------------------------------------------------------------

subroutine output_dry(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only: sx,ex,sy,ey,nx,ny,nz,iaer,igas,isize,nzz,ntt
   use naqpms_varlist, only: ip2mem,ip2memgas,ip2memaer,ip3memaer,naerbin,naersp,iedgas,&
                             dryvelaer,DryVelGas,dryveldust
      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, dim2d(2), dim3d(3)
   integer :: mem3d, mem2d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), fill_value
   real, dimension(:,:),allocatable :: chem2d
   integer,pointer :: ldof2d(:)
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
   character*3 cdnum
   
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
      
   allocate(varid_x2a(17+iedgas))  
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
      
   fname='geatm.dry.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
      
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d) 
 
   rcode = pio_def_dim(File,'lon',nx(ne),dim2d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim2d(2))

   rcode = pio_def_var(File,'lon',PIO_real,dim2d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim2d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
        
   do k=1,iedgas
   write(cdnum(1:3),'(i3.3)')k
   rcode = pio_def_var(File,'dryvelgas'//cdnum,PIO_REAL,dim2d,varid_x2a(2+k))
   rcode = pio_put_att(File,varid_x2a(2+k),"_fillvalue",fillvalue)
   end do
   do k=1,7
   write(cdnum(1:3),'(i3.3)')k
   rcode = pio_def_var(File,'dryvelaer'//cdnum,PIO_REAL,dim2d,varid_x2a(iedgas+2+k))
   end do
   k=0
   do ia=1,iaer
   do is=1,isize
     k=k+1
     write(cdnum(1:3),'(i3.3)')is
     if(ia.eq.1) rcode = pio_def_var(File,'dryvelseasalt'//cdnum,PIO_REAL,dim2d,varid_x2a(iedgas+9+k))
     if(ia.eq.2) rcode = pio_def_var(File,'dryveldust'//cdnum,PIO_REAL,dim2d,varid_x2a(iedgas+9+k))
   enddo
   enddo
   rcode = pio_enddef(File)
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
          
   do k=1,iedgas        
   mem2d=ip2memGas(k,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=DryVelGas(km)    
          end do
          end do
          call pio_write_darray(File, varid_x2a(2+k), iodesc2d, chem2d, rcode)    
   end do  

   k=0
   do ia=1,6
   do is=1,1
     mem2d=ip2memaer(is,ia,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=DryVelaer(km)   
          end do
          end do
     k=k+1
     write(cdnum(1:3),'(i3.3)')k
     call pio_write_darray(File, varid_x2a(iedgas+2+k), iodesc2d, chem2d, rcode)
   enddo
   enddo

   do ia=1,1
   do is=2,2
     mem2d=ip2memaer(is,ia,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=dryvelaer(km) 
          end do
          end do
     call pio_write_darray(File, varid_x2a(iedgas+9), iodesc2d, chem2d, rcode)
   enddo
   enddo
        
  k=0
  do ia=1,iaer
  do is=1,isize
    mem2d = ip3memaer(is,ia,ne)
    do i= sx(ne),ex(ne)
    do j= sy(ne),ey(ne)
       km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
       chem2d(i-sx(ne)+1,j-sy(ne)+1)=dryveldust(km) 
    end do
    end do
     k=k+1
     call pio_write_darray(File, varid_x2a(iedgas+9+k), iodesc2d, chem2d, rcode)
  enddo
  enddo

   call pio_freedecomp(File,iodesc2d)
   call pio_closefile(file)     
   
   deallocate(ldof2d)
   deallocate(chem2d)
   deallocate(varid_x2a)
   end subroutine output_dry
   
!-------------------------------------------------------------------------------        

subroutine output_wet(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only: sx,ex,sy,ey,nx,ny,nz,iaer,igas,isize,nzz,ntt
   use naqpms_varlist, only: ip2mem,ip3mem,ip2memgas,WETDEP,WETDEP2,iedgas,&
                             ip2memaer,aerom_wdepfld,aerom_wdepfld2,naersp,naerbin
   use met_fields,only: RAINCON, RAINNON
      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, dim2d(2), dim3d(3), dim4d(4), dim5d(5)
   integer :: mem3d, mem2d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(igas), fill_value
   real, dimension(:,:),allocatable :: chem2d
   real, dimension(:,:,:),allocatable :: chem3d
   integer,pointer :: ldof3d(:),ldof2d(:)
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   character*10 date
   character*30 fname
   character*3 cdnum
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
      
   allocate(varid_x2a(18+2*iedgas))
   allocate(chem2d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1))
      
   fname='geatm.wet.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)   
 
   rcode = pio_def_dim(File,'lon',nx(ne),dim2d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim2d(2))
   rcode = pio_def_var(File,'lon',PIO_real,dim2d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim2d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'RAINNON',PIO_REAL,dim2d,varid_x2a(3))
   rcode = pio_put_att(File,varid_x2a(3),"_fillvalue",fillvalue)
   rcode = pio_def_var(File,'RAINCON',PIO_REAL,dim2d,varid_x2a(4))
   rcode = pio_put_att(File,varid_x2a(4),"_fillvalue",fillvalue)
   do i=1,iedgas
   write(cdnum(1:3),'(i3.3)')i
   rcode = pio_def_var(File,'wetdep'//cdnum,PIO_REAL,dim2d,varid_x2a(4+i))       
   rcode = pio_def_var(File,'wetdep2'//cdnum,PIO_REAL,dim2d,varid_x2a(iedgas+4+i))      
   rcode = pio_put_att(File,varid_x2a(4+i),"_fillvalue",fillvalue)
   rcode = pio_put_att(File,varid_x2a(iedgas+4+i),"_fillvalue",fillvalue)
   end do
   do i=1,7
   write(cdnum(1:3),'(i3.3)')i
   rcode = pio_def_var(File,'aerom_wdepfld'//cdnum,PIO_REAL,dim2d,varid_x2a(2*iedgas+4+i))
   rcode = pio_def_var(File,'aerom_wdepfld2'//cdnum,PIO_REAL,dim2d,varid_x2a(2*iedgas+11+i))
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)

   ! rainnon
          i0=ip2mem(ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=rainnon(km)                           
          end do
          end do
   call pio_write_darray(File, varid_x2a(3), iodesc2d, chem2d, rcode)
   ! raincon
          i0=ip2mem(ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=i0+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=raincon(km)                           
          end do
          end do
   call pio_write_darray(File, varid_x2a(4), iodesc2d, chem2d, rcode)   

   do k=1,iedgas        
   ! WEPDEP
   mem2d=ip2memGas(k,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=WETDEP(km)    
          end do
          end do
          call pio_write_darray(File, varid_x2a(4+k), iodesc2d, chem2d, rcode)    
   ! WEPDEP2
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=wetdep2(km)   
          end do
          end do
         call pio_write_darray(File, varid_x2a(iedgas+4+k), iodesc2d, chem2d, rcode)
   end do  

   k=0
   do ia=1,6
   do is=1,1
     mem2d=ip2memaer(is,ia,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=aerom_wdepfld(km)   
          end do
          end do
     k=k+1
     call pio_write_darray(File, varid_x2a(2*iedgas+4+k), iodesc2d, chem2d, rcode)
   enddo
   enddo

   do ia=1,1
   do is=2,2
     mem2d=ip2memaer(is,ia,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=aerom_wdepfld(km) 
          end do
          end do
     call pio_write_darray(File, varid_x2a(2*iedgas+11), iodesc2d, chem2d, rcode)
   enddo
   enddo

   k=0
   do ia=1,6
   do is=1,1
     mem2d=ip2memaer(is,ia,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=aerom_wdepfld2(km)   
          end do
          end do
     k=k+1
     call pio_write_darray(File, varid_x2a(2*iedgas+11+k), iodesc2d, chem2d, rcode)
   enddo
   enddo

   do ia=1,1
   do is=2,2
     mem2d=ip2memaer(is,ia,ne)
          do i= sx(ne),ex(ne)
          do j= sy(ne),ey(ne)
            km=mem2d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem2d(i-sx(ne)+1,j-sy(ne)+1)=aerom_wdepfld2(km) 
          end do
          end do
     call pio_write_darray(File, varid_x2a(2*iedgas+18), iodesc2d, chem2d, rcode)
   enddo
   enddo
 
   call pio_freedecomp(File,iodesc2d)
   call pio_closefile(file)     
   
   deallocate(ldof2d)
   deallocate(chem2d)
   deallocate(varid_x2a)
   end subroutine output_wet
   
!-------------------------------------------------------------------------------

subroutine output_sea(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only: sx,ex,sy,ey,nx,ny,nz,iaer,igas,&
       isize,nzz,ntt
   use naqpms_varlist, only: ip2mem,ip3mem,ip5mem,aer

   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, i0, dim2d(2), dim3d(3)
   integer :: mem5d,mem2d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), fill_value
   real, dimension(:,:,:),allocatable :: chem3d
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   integer,pointer :: ldof2d(:),ldof3d(:)
   character*10 date
   character*30 fname
   
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nz(ne)
   sigma(i)=i
   enddo
   
   allocate(varid_x2a(7))
   allocate(chem3d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne)))      
   
   fname='geatm.sea.'//date//'.nc'
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)
   allocate(ldof2d(m)) 
   m=0
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof2d(m)=(j-1)*nx(ne)+i
   end do
   end do

   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*nz(ne)
   allocate(ldof3d(m))
   m=0
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof3d(m)=(j-1)*nx(ne)+i+(k-1)*nx(ne)*ny(ne)
   end do
   end do
   end do

   dim2d=(/nx(ne),ny(ne)/)
   dim3d=(/nx(ne),ny(ne),nz(ne)/)
   call pio_initdecomp(pio_subsystem, pio_real, dim2d, ldof2d, iodesc2d)    
   call pio_initdecomp(pio_subsystem, pio_real, dim3d, ldof3d, iodesc3d)    
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim3d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim3d(2))
   dim2d(1:2)=dim3d(1:2)
   rcode = pio_def_dim(File,'lev',nz(ne),dim3d(3))  

   rcode = pio_def_var(File,'lon',PIO_real,dim3d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim3d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim3d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'sea01',PIO_REAL,dim3d,varid_x2a(4))     
   rcode = pio_def_var(File,'sea02',PIO_REAL,dim3d,varid_x2a(5))     
   rcode = pio_def_var(File,'sea03',PIO_REAL,dim3d,varid_x2a(6))     
   rcode = pio_def_var(File,'sea04',PIO_REAL,dim3d,varid_x2a(7))     
   mm=7
        
   do i=1,mm
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)

      do k=1,nzz        
        mem5d=ip5mem(k,1,1, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)  
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(4), iodesc3d, chem3d, rcode)

      do k=1,nzz
        mem5d=ip5mem(k,2,1, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(5), iodesc3d, chem3d, rcode)

      do k=1,nzz
        mem5d=ip5mem(k,3,1, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(6), iodesc3d, chem3d, rcode)

      do k=1,nzz
        mem5d=ip5mem(k,4,1, ne)
        do j= sy(ne),ey(ne)
        do i= sx(ne),ex(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem3d(i-sx(ne)+1,j-sy(ne)+1,k)=aer(km)
        end do
        end do
      end do
     call pio_write_darray(File, varid_x2a(7), iodesc3d, chem3d, rcode)
        
   call pio_freedecomp(File,iodesc2d)
   call pio_freedecomp(File,iodesc3d)
   call pio_closefile(file)     
   
   deallocate(ldof2d)
   deallocate(ldof3d)
   deallocate(chem3d)
   deallocate(varid_x2a)
   end subroutine output_sea
   
!-------------------------------------------------------------------------------

subroutine output_mark(iyear, imonth, iday, ihour, ne, restart)
   use geatm_vartype, only:  sx,ex,sy,ey,nx,ny,nz,&
       iaer,igas,isize,nzz,ntt,idmset,ismmax
   use naqpms_varlist, only:ip2mem,ip3mem,sourcemark
      
   integer :: iyear, imonth, iday, ihour, ne
   logical :: restart
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, ig, is, ia, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   integer :: mem5d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(ismmax), b(idmset)
   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer         :: m, mm        ! numbr of variables
   integer         :: rcode        ! return error code
   integer,pointer :: ldof5d(:)
   character*10 date
   character*30 fname

   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)   
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo
   do i=1,ismMax
   a(i)=i
   enddo
   do i=1,idmSet
   b(i)=i
   enddo
   
   allocate(varid_x2a(6))
      
   allocate(chem5d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne),ismmax,idmset))
   
   if(restart) then
   fname='geatm.restd.'//date//'.nc'
   else
   fname='geatm.mark.'//date//'.nc'
   endif
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)   
   
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*ismmax*nz(ne)*idmset
   allocate(ldof5d(m))
   m=0
   do ia=1,idmset
   do is=1,ismmax
   do k=1,nz(ne)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof5d(m)=(j-1)*nx(ne)+(k-1)*nx(ne)*ny(ne)+(is-1)*nx(ne)*ny(ne)*nz(ne)+&
             (ia-1)*nx(ne)*ny(ne)*nz(ne)*ismmax
   end do
   end do
   end do
   end do
   end do

   dim5d=(/nx(ne),ny(ne),nz(ne),ismmax,idmset/)
   call pio_initdecomp(pio_subsystem, pio_real, dim5d, ldof5d, iodesc5d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim5d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim5d(2))
   rcode = pio_def_dim(File,'lev',nzz,dim5d(3))  
   rcode = pio_def_dim(File,'ismmax',ismmax,dim5d(4))
   rcode = pio_def_dim(File,'idmset',idmset,dim5d(5))   

   rcode = pio_def_var(File,'lon',PIO_real,dim5d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim5d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim5d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'ismmax',PIO_REAL,dim5d(4:4),varid_x2a(4))
   rcode = pio_def_var(File,'idmset',PIO_REAL,dim5d(5:5),varid_x2a(5))   
   rcode = pio_def_var(File,'sourcemark',PIO_REAL,dim5d,varid_x2a(6))
        
   do i=1,6
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   rcode = pio_put_var(File, varid_x2a(4), a)
   rcode = pio_put_var(File, varid_x2a(5), b)   
          
   ! seacomp
   mem5d=1        ! number of memory for 5D variables  
   do ia=1,idmset
   do is=1,ismmax  
   do k=1,nzz        
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=mem5d+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem5d(i-sx(ne)+1,j-sy(ne)+1,k,is,ia)=sourcemark(km)                           
        end do
        end do
        mem5d=mem5d+(nx(ne)+2)*(ny(ne)+2)
   end do
   end do
   end do
   call pio_write_darray(File, varid_x2a(6), iodesc5d, chem5d, rcode)
        
   call pio_freedecomp(File,iodesc5d)
   call pio_closefile(file)     
      
   deallocate(ldof5d)
   deallocate(chem5d)
   deallocate(varid_x2a)
   end subroutine output_mark
!-------------------------------------------------------------------------------

subroutine output_term(iyear, imonth, iday, ihour, ne)
   use geatm_vartype, only:  sx,ex,sy,ey,nx,ny,nz,&
       iaer,igas,isize,nzz,ntt,iprocess,PrintTermGas

   use naqpms_varlist, only: ip2mem,ip3mem,ip5mem,IGOPos,&
       ipGasTermBal,GasTermBal
      
   integer :: iyear, imonth, iday, ihour, ne
   
   type(file_desc_t) :: file
   type(io_desc_t) :: iodesc2d,iodesc3d,iodesc4d,iodesc5d
   type(var_desc_t), pointer :: varid_x2a(:)
   integer :: i, j, k, i05, igo, is, ia, mlen, dim2d(2), dim3d(3), dim4d(4),dim5d(5)
   integer :: mem5d
   real :: xlat(ny(ne)), xlon(nx(ne)), sigma(nzz), a(iprocess), b(igas), fill_value
   real, dimension(:,:,:,:,:),allocatable :: chem5d
   integer         :: m,mm        ! numbr of variables
   integer         :: rcode        ! return error code
   integer,pointer :: ldof5d(:)
   character*10 date
   character*30 fname
  
   date(1:4)=char(iyear/1000+48)//char(iyear/100-iyear/1000*10+48)//char(iyear/10-iyear/100*10+48)//&
             char(iyear-iyear/10*10+48)
   date(5:6)=char(imonth/10+48)//char(imonth-imonth/10*10+48)
   date(7:8)=char(iday/10+48)//char(iday-iday/10*10+48)
   date(9:10)=char(ihour/10+48)//char(ihour-ihour/10*10+48)
   fillvalue = -99999.9

   do i=1,ny(ne)
   xlat(i)=(i-1)*1.0-89.5
   enddo
   do i=1,nx(ne)
   xlon(i)=(i-1)*1.0+0.5
   enddo
   do i=1,nzz
   sigma(i)=i
   enddo
   do i=1,iprocess
   a(i)=i
   enddo

   mlen=0
   do i=1,igas
   if(PrintTermGas(i).eq.1)then
   mlen=mlen+1
   b(i)=mlen
   end if
   enddo
   if (mlen.eq.0) then
   go to 9999
   end if
   
   allocate(varid_x2a(6))
   
   allocate(chem5d(ex(ne)-sx(ne)+1,ey(ne)-sy(ne)+1,nz(ne),iprocess,mlen))
   
   fname='geatm.term.'//date//'.nc'   
   rcode = pio_createfile(pio_subsystem, file, pio_iotype_netcdf, trim(adjustl(fname)), PIO_CLOBBER)   
       
   m=(ex(ne)-sx(ne)+1)*(ey(ne)-sy(ne)+1)*iprocess*nz(ne)*m
   allocate(ldof5d(m))
   m=0
   do ia=1,mlen
   do is=1,iprocess
   do k=1,nzz
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
   m=m+1
   ldof5d(m)=(j-1)*nx(ne)+(k-1)*nx(ne)*ny(ne)+(is-1)*nx(ne)*ny(ne)*nz(ne)+&
             (ia-1)*nx(ne)*ny(ne)*nz(ne)*isize
   end do
   end do
   end do
   end do
   end do

   dim5d=(/nx(ne),ny(ne),nz(ne),iprocess,mlen/)
   call pio_initdecomp(pio_subsystem, pio_real, dim5d, ldof5d, iodesc5d)
   
   rcode = pio_def_dim(File,'lon',nx(ne),dim5d(1))
   rcode = pio_def_dim(File,'lat',ny(ne),dim5d(2))
   rcode = pio_def_dim(File,'lev',nzz,dim5d(3))  
   rcode = pio_def_dim(File,'iprocess',iprocess,dim5d(4))
   rcode = pio_def_dim(File,'igas',mlen,dim5d(5))   

   rcode = pio_def_var(File,'lon',PIO_real,dim5d(1:1),varid_x2a(1))
   rcode = pio_put_att (File, varid_x2a(1), 'long_name', 'longitude')
   rcode = pio_put_att (File, varid_x2a(1), 'units', 'degrees_east')
   rcode = pio_def_var(File,'lat',PIO_REAL,dim5d(2:2),varid_x2a(2))
   rcode = pio_put_att (File, varid_x2a(2), 'long_name', 'latitude')
   rcode = pio_put_att (File, varid_x2a(2), 'units', 'degrees_north')
   rcode = pio_def_var(File,'lev',PIO_REAL,dim5d(3:3),varid_x2a(3))
   rcode = pio_def_var(File,'iprocess',PIO_REAL,dim5d(4:4),varid_x2a(4))
   rcode = pio_def_var(File,'igas',PIO_REAL,dim5d(5:5),varid_x2a(5))   
   rcode = pio_def_var(File,'GasTermBal',PIO_REAL,dim5d,varid_x2a(6))
        
   do i=1,6
   rcode = pio_put_att(File,varid_x2a(i),"_fillvalue",fillvalue)
   end do
   rcode = pio_enddef(File) 
   
   rcode = pio_put_var(File, varid_x2a(1), xlon)
   rcode = pio_put_var(File, varid_x2a(2), xlat)
   rcode = pio_put_var(File, varid_x2a(3), sigma)
   rcode = pio_put_var(File, varid_x2a(4), a)
   rcode = pio_put_var(File, varid_x2a(5), b(1:mlen))   
          
   ! GasTermBal
   do ia=1,igas
   if(PrintTermGas(ia)==1)then
   igo = IGOPos(ia) 
   do is=1,iprocess
   do k=1,nzz
        i05=ipGasTermBal(k,is,igo,ne)
        do i= sx(ne),ex(ne)
        do j= sy(ne),ey(ne)
            km=i05+(ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
            chem5d(i-sx(ne)+1,j-sy(ne)+1,k,is,ia)=GasTermBal(km)                           
        end do
        end do   
   end do
   end do
   endif
   end do
   call pio_write_darray(File, varid_x2a(6), iodesc5d, chem5d, rcode)
        
   call pio_freedecomp(File,iodesc5d)
   call pio_closefile(file)     
  
   deallocate(ldof5d)    
   deallocate(chem5d)
   deallocate(varid_x2a)
   9999 continue
   end subroutine output_term

!--------------------------------------------------------------------------------------   
  subroutine init_gea_pio()
    use seq_io_mod,       only: seq_io_getiosys, seq_io_getiotype

    pio_subsystem => seq_io_getiosys('GEA')

  end subroutine init_gea_pio

  end module inout3
