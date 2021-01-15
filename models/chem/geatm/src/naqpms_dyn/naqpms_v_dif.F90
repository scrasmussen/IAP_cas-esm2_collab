subroutine naqpms_v_dif &
  ( myid &
   ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
   ,igasCBM &
   ,iwb,ieb,jsb,jeb &
   ,dzz &
   ,rkv,ttn,ppp,atm,kktop &
   ,igas,iaer,isize,nseacom,ndustcom &
   ,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use naqpms_gridinfo

implicit none

integer :: myid

real :: dt

integer :: igasCBM

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom


integer :: ifsm(5)

integer :: ifsmt

integer :: idmSet,ismMax

integer :: igMark(idmSet)

real,allocatable,dimension(:,:,:,:) :: smconv
real,allocatable,dimension(:,:,:) :: concmark

integer :: letdoit,idm,ism,i04sm,i04aer

integer :: i,j,k,is

integer :: ixy,iapm,i03

integer :: ig,i04,i05,i05c,iduc,ia
integer :: jsb1,jeb1,iwb1,ieb1


integer :: iwb,ieb,jsb,jeb
real,dimension(iwb:ieb,jsb:jeb,nzz) :: ppp,ttn,conc,atm,rkv,dzz
real,dimension(iwb:ieb,jsb:jeb    ) :: kktop

real :: SEA(iwb:ieb,jsb:jeb,nzz,nseacom),DUST(iwb:ieb,jsb:jeb,nzz,ndustcom)

integer :: ispflag

DO ig=1,iedgas

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i04=ip4mem(k,ig,ne)
       conc(i,j,k)=gas(i04+ixy)
       if(ifsm(ne).eq.1) concmark(i,j,k)=gas(i04+ixy)
     enddo !k
   enddo!i
   enddo !j
   jsb1=sy(ne);jeb1=ey(ne)
   iwb1=sx(ne);ieb1=ex(ne)
   ispflag=0
   call  diffus(myid,ispflag,iwb1,ieb1,jsb1,jeb1,nzz,dt,rkv,dzz,ttn,ppp,conc,atm,kktop)
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
        i04=ip4mem(k,ig,ne)
        i03=ip3mem(k,ne)
        gas(i04+ixy)=conc(i,j,k)
     enddo !k
   enddo !i
   enddo !j

   if(ifsm(ne)==1)then
     letdoit=0
     do idm=1,idmSet
      if(igMark(idm)==ig) letdoit=idm
     enddo

     if(letdoit>0)then
       iwb=sx(ne)-1;ieb=ex(ne)+1
       jsb=sy(ne)-1;jeb=ey(ne)+1
       allocate(smconv(ismMax,iwb:ieb,jsb:jeb,nzz))
       do j=sy(ne),ey(ne)
       do i=sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz
         do ism=1,ismMax
           i04sm=ipSMmem(k,ism,letdoit,ne)
           smconv(ism,i,j,k)=SourceMark(i04sm+ixy)
         enddo !is
         enddo !k
       enddo!i
       enddo !j
       jsb1=sy(ne);jeb1=ey(ne)
       iwb1=sx(ne);ieb1=ex(ne)
       call diffus_mark(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dt,rkv,dzz,ttn,ppp,concmark,atm,ismMax,smconv,kktop)

       do j=sy(ne),ey(ne)
       do i=sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz
         do ism=1,ismMax
           i04sm=ipSMmem(k,ism,letdoit,ne)
           SourceMark(i04sm+ixy)=smconv(ism,i,j,k)
         enddo !is
         enddo !k
       enddo!i
       enddo !j

       if(allocated(smconv)) deallocate(smconv)
   endif !letdoit
  endif !ifsm(ne)
  !!!!!!!!!!!!!! 
 
ENDDO !ig


if(laerv2) then

 do ia=1,naersp
 do is=1,naerbin

   if(ia.gt.1.and.is.gt.1) cycle ! skip zero aersol tracer

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
        i04aer=ip4mem_aer(k,is,ia,ne)
        conc(i,j,k)=aerom(i04aer+ixy)
     enddo
   enddo
   enddo

   jsb1=sy(ne);jeb1=ey(ne)
   iwb1=sx(ne);ieb1=ex(ne)
   ispflag=9999
   call diffus(myid,ispflag,iwb1,ieb1,jsb1,jeb1,nzz,dt,rkv,dzz,ttn,ppp,conc,atm,kktop)

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
       i04aer=ip4mem_aer(k,is,ia,ne)
       aerom(i04aer+ixy)=conc(i,j,k)
     enddo !k
   enddo !i
   enddo !j

 enddo
 enddo
endif



 do IA=1,IAER
 do IS=1,ISIZE

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)

      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

      do k=1,nzz
        i05=ip5mem(k,is,ia,ne)
        conc(i,j,k)=aer(i05+ixy)

        IF(IA==1) THEN
          do iduc = 1,nseacom
            i05c = ip5memcs (k,is,iduc,ne)
            sea(i,j,k,iduc) = SEACOMP(i05c+ixy)
          enddo
        ELSE IF(IA==2) THEN
          do iduc = 1, ndustcom
            i05c = ip5memc (k,is,iduc,ne)
            dust(i,j,k,iduc) = DUSTCOMP(i05c+ixy)
          enddo ! iduc 
        ENDIF
      enddo !k

   enddo!i
   enddo !j

   jsb1=sy(ne);jeb1=ey(ne)
   iwb1=sx(ne);ieb1=ex(ne)

   if(ia==1) then
      call diffus_ds(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dt,rkv,dzz,ttn,ppp,conc,sea,atm,kktop,nseacom,is,ia)
   else if(ia==2)then
      call diffus_ds(myid,ig,iwb1,ieb1,jsb1,jeb1,nzz,dt,rkv,dzz,ttn,ppp,conc,dust,atm,kktop,ndustcom,is,ia)
   endif

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)

     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

     do k=1,nzz
       i05=ip5mem(k,is,ia,ne)
       aer(i05+ixy)=conc(i,j,k)

       IF(IA==1) THEN
        do iduc = 1,nseacom
         i05c = ip5memcs (k,is,iduc,ne)
         SEACOMP(i05c+ixy) = sea(i,j,k,iduc)
        enddo
       ELSE IF(IA==2) THEN
        do iduc = 1, ndustcom
         i05c = ip5memc (k,is,iduc,ne)
         DUSTCOMP(i05c+ixy) = dust(i,j,k,iduc)
        enddo ! iduc
       ENDIF
      enddo !k

   enddo!i
   enddo !j

 ENDDO ! IS
 ENDDO ! IA


end subroutine naqpms_v_dif

