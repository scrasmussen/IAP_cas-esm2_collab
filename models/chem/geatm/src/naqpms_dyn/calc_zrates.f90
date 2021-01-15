
! read meteorological fields 
! and calculate vertical mixing coefficient, cloud parameters

subroutine calc_zrates( myid,iyear,imonth,iday,ihour,iminute &
                       ,nest,nzz,nx,ny,nz,sx,ex,sy,ey,ne )

use naqpms_varlist, only: ip2mem,ip3mem
use met_fields, only : u,v,roair3d,entrn3d,dilut3d,dsdt3d
use naqpms_gridinfo, only: dx,dy,dz,mpfac,pzps 


implicit none 

logical,parameter :: lrain_1h=.false.

real,parameter :: deg2rad=0.01745329
real,parameter :: gamma=0.286

real,parameter :: kvmin=0.5

real,parameter :: denh2o=1.0e6 ! g/m3
real,parameter :: fspfa=0.86 ! scattering phase function asymmetry factor


integer :: myid
integer :: ne,nest
integer :: nx(5),ny(5),nz(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)


integer :: iyear,imonth,iday,ihour,iminute


integer :: i,j,k,i0,i02,i03,ixy,i03_kp1

character :: cdnum*1
character :: date*10,fname*100

integer :: irec
integer :: funit
logical :: lexist

!=========

real :: windw

real :: deltax, deltay, deltaz,dvgc,dvgcx,dvgcy,rho_z,areayz,areaxz,m2fac

real :: dvgss,deltas

real,allocatable,dimension(:,:,:) :: dat3d 

real,allocatable,dimension(:,:) :: den3d,u3d,v3d,uwd_x,vwd_y,rho_x,rho_y

real,allocatable,dimension(:,:) :: dx3d,dy3d,dz3d,areayz_x,areaxz_y 

real,allocatable,dimension(:,:) :: mpfac2d,fac2d_x,fac2d_y

real,allocatable,dimension(:,:) :: pzps3d,pzps_x,pzps_y

real,allocatable,dimension(:,:) :: dvgc3d,dvgss3d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


allocate( den3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
         ,u3d  (sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
         ,v3d  (sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
         ,dx3d (sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
         ,dy3d (sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
         ,dz3d (sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) &
         ,pzps3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1) )

allocate(mpfac2d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1))



allocate( uwd_x(sx(ne)-1:ex(ne),sy(ne):ey(ne)) &
         ,rho_x(sx(ne)-1:ex(ne),sy(ne):ey(ne)) &
         ,areayz_x(sx(ne)-1:ex(ne),sy(ne):ey(ne)) &
         ,pzps_x(sx(ne)-1:ex(ne),sy(ne):ey(ne)) )


allocate( fac2d_x(sx(ne)-1:ex(ne),sy(ne):ey(ne)) ) 


allocate( vwd_y(sx(ne):ex(ne),sy(ne)-1:ey(ne)) &
         ,rho_y(sx(ne):ex(ne),sy(ne)-1:ey(ne)) &
         ,areaxz_y(sx(ne):ex(ne),sy(ne)-1:ey(ne)) &
         ,pzps_y(sx(ne):ex(ne),sy(ne)-1:ey(ne)) )

allocate( fac2d_y(sx(ne):ex(ne),sy(ne)-1:ey(ne))  )

allocate( dvgc3d(sx(ne):ex(ne),sy(ne):ey(ne)) )
allocate( dvgss3d(sx(ne):ex(ne),sy(ne):ey(ne)) )

!allocate( dat3d(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz)  )

allocate( dat3d(sx(ne):ex(ne),sy(ne):ey(ne),nzz)  )

!stop



 do j=sy(ne),ey(ne)
 do i=sx(ne),ex(ne)
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do k=1,nzz
     i03=ip3mem(k,ne)
     entrn3d(i03+ixy)=0.0 ! phpt-windw
     dilut3d(i03+ixy)=0.0 ! phpt=0.0
     dsdt3d(i03+ixy)=0.0
   enddo
 enddo
 enddo



!===================================
!=================================
! (2) calculate cloud parameters

i02=ip2mem(ne)


 do j=sy(ne)-1,ey(ne)+1
 do i=sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   mpfac2d(i,j)=mpfac(i02+ixy)
 enddo
 enddo

 do j=sy(ne),ey(ne)
   do i=sx(ne)-1,ex(ne)
     fac2d_x(i,j)=0.5*(mpfac2d(i,j)+mpfac2d(i+1,j))
   enddo
 enddo


 do j=sy(ne)-1,ey(ne)
   do i=sx(ne),ex(ne)
     fac2d_y(i,j)=0.5*(mpfac2d(i,j)+mpfac2d(i,j+1))
   enddo
 enddo


dvgc3d=0.0
dvgss3d=0.0

do k=1,nzz-1

  i03=ip3mem(k,ne)

  i03_kp1=ip3mem(k+1,ne)

  do j=sy(ne)-1,ey(ne)+1
  do i=sx(ne)-1,ex(ne)+1
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    den3d(i,j)=roair3d(i03+ixy)
    u3d(i,j)=u(i03+ixy)
    v3d(i,j)=v(i03+ixy)
    dx3d(i,j)=dx(i03+ixy)
    dy3d(i,j)=dy(i03+ixy)
    dz3d(i,j)=dz(i03+ixy)
    pzps3d(i,j)=pzps(i03+ixy)
  enddo
  enddo

  do j=sy(ne),ey(ne)
    do i=sx(ne)-1,ex(ne)
      rho_x(i,j)=0.5*(den3d(i,j)+den3d(i+1,j))
      uwd_x(i,j)=0.5*(u3d(i,j)+u3d(i+1,j))
      areayz_x(i,j)=0.25*(dy3d(i,j)+dy3d(i+1,j))*(dz3d(i,j)+dz3d(i+1,j))
      pzps_x(i,j)=0.5*(pzps3d(i,j)+pzps3d(i+1,j))
    enddo
  enddo

  do i=sx(ne),ex(ne)
    do j=sy(ne)-1,ey(ne)
      rho_y(i,j)=0.5*(den3d(i,j)+den3d(i,j+1))
      vwd_y(i,j)=0.5*(v3d(i,j)+v3d(i,j+1))
      areaxz_y(i,j)=0.25*(dx3d(i,j)+dx3d(i,j+1))*(dz3d(i,j)+dz3d(i,j+1))
      pzps_y(i,j)=0.5*(pzps3d(i,j)+pzps3d(i,j+1))
    enddo
  enddo



  do j=sy(ne),ey(ne)
  do i=sx(ne),ex(ne)

    ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1

    rho_z = 0.5*(roair3d(i03+ixy)+roair3d(i03_kp1))

    deltaz=dz(i03+ixy)
    deltay=dy(i03+ixy)
    deltax=dx(i03+ixy)
 
    deltas=dz(i03+ixy)/pzps3d(i,j)


    areayz=deltay*deltaz
    areaxz=deltax*deltaz
    m2fac=mpfac2d(i,j)**2.0


!for camx ppm advection 
    dvgcx =  uwd_x(i,j)*rho_x(i,j)*areayz_x(i,j)/fac2d_x(i,j) &
            -uwd_x(i-1,j)*rho_x(i-1,j)*areayz_x(i-1,j)/fac2d_x(i-1,j) 

    dvgcx = dvgcx*m2fac/areayz/deltax
  
    dvgcy = vwd_y(i,j)*rho_y(i,j)*areaxz_y(i,j)/fac2d_y(i,j) &
           -vwd_y(i,j-1)*rho_y(i,j-1)*areaxz_y(i,j-1)/fac2d_y(i,j-1)

    dvgcy = dvgcy*m2fac/areaxz/deltay

    dvgc3d(i,j)  = dvgc3d(i,j) + (dvgcx+dvgcy)*deltaz

    windw = -dvgc3d(i,j)/rho_z

    entrn3d(i03+ixy)=dilut3d(i03+ixy)-windw

!    dat3d(i,j,k)=entrn3d(i03+ixy)


! for walcek advection

    dvgcx = uwd_x(i,j)*rho_x(i,j)*pzps_x(i,j)/fac2d_x(i,j) &
           -uwd_x(i-1,j)*rho_x(i-1,j)*pzps_x(i-1,j)/fac2d_x(i-1,j)
    dvgcx = dvgcx*mpfac2d(i,j)**2.0/deltax

    dvgcy = vwd_y(i,j)*rho_y(i,j)*pzps_y(i,j)/fac2d_y(i,j) &
           -vwd_y(i,j-1)*rho_y(i,j-1)*pzps_y(i,j-1)/fac2d_y(i,j-1)
    dvgcy = dvgcy*mpfac2d(i,j)**2.0/deltay

    dvgss3d(i,j) = dvgss3d(i,j) + (dvgcx + dvgcy)*deltas  ! integrating

!cycle

    dsdt3d(i03+ixy)=-dvgss3d(i,j)/rho_z/pzps3d(i,j)

    dat3d(i,j,k)=dsdt3d(i03+ixy)

  enddo
  enddo

enddo

!print*,'dat3d=',dat3d(:,:,1)

deallocate(u3d,v3d,dx3d,dy3d,dz3d)

deallocate(uwd_x,vwd_y,rho_x,rho_y,dat3d)
deallocate(areayz_x,areaxz_y)
deallocate(mpfac2d,fac2d_x,fac2d_y)
deallocate(pzps3d,pzps_x,pzps_y)

deallocate(dvgc3d,dvgss3d)

!stop 'zrates ok'

end subroutine calc_zrates

