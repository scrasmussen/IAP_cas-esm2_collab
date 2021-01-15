
subroutine cal_trop_pbltop &
 & ( myid &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d,ktop,tropp )

use naqpms_varlist, only : ip3mem,ip2mem 
use naqpms_gridinfo, only : heiz,terrain
use met_fields, only : Plev,t,PBL_HGT,NPBL
implicit none

integer :: myid

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: mem2d
real,dimension(mem2d) :: ktop,tropp


integer :: i,j,k
integer :: iwb,ieb,jsb,jeb

integer :: i0,i03,i02,ixy,i03_t

real    :: plimu,pliml,plimlex,tperr
logical :: dofill

real,allocatable,dimension(:,:,:) :: ppp,ttn
real,allocatable,dimension(:,:) :: tp



  ! find the layer of top conditons from global model 

  ! (1) trop top

  iwb=sx(ne)-1;ieb=ex(ne)+1
  jsb=sy(ne)-1;jeb=ey(ne)+1

  allocate(ppp(iwb:ieb,jsb:jeb,nzz),ttn(iwb:ieb,jsb:jeb,nzz),tp(iwb:ieb,jsb:jeb))

  do j = sy(ne),ey(ne)
  do i = sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i03=ip3mem(k,ne)
       ppp(i,j,k)=Plev(i03+ixy)
       ttn(i,j,k)=t(i03+ixy)
     enddo  !K
  enddo  !i
  enddo !j   

  plimu    = 45000.
  pliml    = 7500.
  plimlex  = 10000.
  dofill   = .TRUE.

  CALL tropo( myid,ttn, sx(ne),ex(ne),sy(ne),ey(ne),nzz &
             ,ppp, plimu, pliml, plimlex,dofill, tp, tperr)

  ! (2) pbl top

    
  do j = sy(ne),ey(ne)
  do i = sx(ne),ex(ne)

      i0 = ip2mem(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      tropp(i0+ixy)=tp(i,j)/100. !PA-->HPA
      NPBL (i0+ixy) = 1

      do k=1,nzz

         i03 = ip3mem(k,ne)
         i02=ip2mem(ne)
         if(k==nzz) then
          i03_t = ip3mem(k,ne)
         else 
          i03_t = ip3mem(k+1,ne)  
         endif     
         if(Plev(i03+ixy)>=tropp(i0+ixy).and.Plev(i03_t+ixy)<tropp(i0+ixy)) then
           ktop(i02+ixy)=float(k)   
         endif 

         IF(heiz(i03+ixy)-terrain(i02+ixy).LE. PBL_HGT(i0+ixy)) NPBL (i0+ixy) = K ! get the  layer of pbl eight

      enddo !k

        
      if(i==53.and.j==30.and.ne==1) then
!         print*,'ktop=',ktop(i02+ixy),tropp(i0+ixy)
      endif

  enddo !i
  enddo !j

  deallocate(ttn,ppp,tp)


end subroutine cal_trop_pbltop


