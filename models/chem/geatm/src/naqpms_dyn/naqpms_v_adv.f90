
subroutine naqpms_v_adv &
 & ( myid &
 &  ,imasskeep &
 &  ,ne,dt,nstep &
 &  ,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ktop &
 &  ,mem2d &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,ifsm,idmSet,ismMax,igMark )


use naqpms_varlist
use naqpms_gridinfo
use met_fields
implicit none

integer :: myid

real :: dt
integer :: nstep
real :: dtt0

logical :: lapm
integer :: imasskeep

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)


integer :: ii
integer :: i,j,k,is

integer :: iwb,ieb,jsb,jeb
real,allocatable,dimension(:,:) :: kktop
real,allocatable,dimension(:,:,:) :: uu,vv,ww,conc,ddx,ddz

integer :: mem2d,mem3d

real,dimension(mem3d) :: kpmass_m2


real,dimension(mem2d) :: ktop

integer :: ixy,i03,i02,iapm


integer :: ig,i04,ia,iduc,i05,i05c,i0,idm,ism,i04sm,i04aer

real,allocatable,dimension(:,:,:,:) :: smconv
real,allocatable,dimension(:,:,:) :: concmark

integer :: letdoit


!print*,'w_ws=',w

!print*,'mem2d 3d=',mem2d,mem3d

!if(myid.eq.0) print*,'vadv : nstep=',nstep ! print ok(8.)

!print*,'igas,iaer,isize,nseacom,ndustcom=',igas,iaer,isize,nseacom,ndustcom

!print*,'nx,ny,nzz,nest,sy,ey,sx,ex=',nx,ny,nzz,nest,sy,ey,sx,ex

!stop

iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1

allocate( kktop(iwb:ieb,jsb:jeb) )
allocate( uu(iwb:ieb,jsb:jeb,nzz) )
allocate( vv(iwb:ieb,jsb:jeb,nzz) )
allocate( ww(iwb:ieb,jsb:jeb,nzz) )
allocate( conc(iwb:ieb,jsb:jeb,nzz) )
allocate( ddx(iwb:ieb,jsb:jeb,nzz) )
allocate( ddz(iwb:ieb,jsb:jeb,nzz) )

if(ifsm(ne)==1) allocate( concmark(iwb:ieb,jsb:jeb,nzz) )

!print*,'vdav dim',iwb,ieb,jsb,jeb

!print*,'ktop=',ktop



do j=sy(ne)-1,ey(ne)+1
do i=sx(ne)-1,ex(ne)+1
     !print*,'k-j-i',j,i
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     do k=1,nzz
         i03=ip3mem(k,ne)
         uu(i,j,k) = u (i03+ixy)
         vv(i,j,k) = v (i03+ixy)
         ww(i,j,k) = w (i03+ixy)
         ddx(i,j,k)= dx(i03+ixy)
         ddz(i,j,k)= dz(i03+ixy)
         !if(k.eq.1) print*,uu(i,j,k),vv(i,j,k),ww(i,j,k),ddx(i,j,k),ddz(i,j,k)
     enddo
     i02=ip2mem(ne)
     kktop(i,j) = ktop(i02+ixy)
     !print*,'ktop',kktop(i,j)
enddo
enddo



!loop_spc :  do ig=1,igas
loop_spc :  do ig=1,iedgas

   if(imasskeep==1)then
       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i04=ip4mem(k,ig,ne)
            kpmass_m2(i03+ixy)=gas(i04+ixy)
         enddo
       enddo
       enddo
   endif

   loop_tstep : do ii=1,nstep

      dtt0=dt/float(nstep)

      do j=sy(ne)-1,ey(ne)+1
      do i=sx(ne)-1,ex(ne)+1

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
          i03=ip3mem(k,ne)
          i04=ip4mem(k,ig,ne)
          uu(i,j,k) = u (i03+ixy)
          vv(i,j,k) = v (i03+ixy)
          ww(i,j,k) = w (i03+ixy)
          ddx(i,j,k)= dx(i03+ixy)
          ddz(i,j,k)= dz(i03+ixy)
          conc(i,j,k)=gas(i04+ixy)
          if(ifsm(ne)==1) concmark(i,j,k)=gas(i04+ixy)
        enddo  !k
        i0=ip2mem(ne)
        kktop(i,j) = ktop(i0+ixy)
      enddo   !j
      enddo    !i

      CALL ADV_VERT(MYID,conc,ww,uu,vv,ddz,ddx,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
         i04=ip4mem(k,ig,ne)
         gas(i04+ixy) = conc(i,j,k)
        enddo !k
      enddo  !i
      enddo    !j      

  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  ! to close vertical source mark

  if(ifsm(ne)==1)then

    letdoit=0

    do idm=1,idmSet
       if(igMark(idm)==ig)letdoit=idm
    enddo

    if(letdoit>0)then

      iwb=sx(ne)-1;ieb=ex(ne)+1
      jsb=sy(ne)-1;jeb=ey(ne)+1

      allocate(smconv(ismMax,iwb:ieb,jsb:jeb,nzz))

      do j=sy(ne)-1,ey(ne)+1
      do i=sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
        do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
          smconv(ism,i,j,k)=SourceMark(i04sm+ixy)
        enddo !is
        enddo !k
      enddo!i
      enddo !j

      CALL ADV_VERT_MARK(MYID,concmark,ww,uu,vv,ddz,ddx,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig,ismMax,smconv,kktop)

      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz
         do ism=1,ismMax
           i04sm=ipSMmem(k,ism,letdoit,ne)
           SourceMark(i04sm+ixy)=smconv(ism,i,j,k)
         enddo !is
         enddo !k
      enddo !i
      enddo !j

      if(allocated(smconv)) deallocate(smconv)

    endif !letdoit

  endif !ifsm(ne)

!!!!!!!!!!!!!!!!!!!!!!!!
  enddo loop_tstep !ii  

enddo loop_spc !ig


if(allocated(concmark)) deallocate(concmark)

!===================


if(laerv2) then
do ia=1,naersp
do is=1,naerbin
   if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i04aer=ip4mem_aer(k,is,ia,ne)
            kpmass_m2(i03+ixy)=aerom(i04aer+ixy)
         enddo
      enddo
      enddo
   endif

   do ii=1,nstep

     dtt0=dt/float(nstep)

     do j=sy(ne)-1,ey(ne)+1
     do i=sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         i04=ip4mem(k,ig,ne)
         i04aer=ip4mem_aer(k,is,ia,ne)
         uu(i,j,k) = u (i03+ixy)
         vv(i,j,k) = v (i03+ixy)
         ww(i,j,k) = w (i03+ixy)
         ddx(i,j,k)= dx(i03+ixy)
         ddz(i,j,k)= dz(i03+ixy)
         conc(i,j,k)=aerom(i04aer+ixy)
         !concmark(i,j,k)=AER(i05+ixy)
       enddo  !k
       i0=ip2mem(ne)
       kktop(i,j) = ktop(i0+ixy)
     enddo   !j
     enddo    !i

     CALL ADV_VERT(MYID,conc,ww,uu,vv,ddz,ddx,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
       i04aer=ip4mem_aer(k,is,ia,ne)
       aerom(i04aer+ixy) = conc(i,j,k)
       enddo !k
     enddo  !i
     enddo    !j

   enddo

enddo
enddo
endif



!--- FOR DUST AND SEA SALT

 DO IA=1,IAER
 DO IS=1,ISIZE


   if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i04=ip4mem(k,ig,ne)
            I05=IP5MEM(K,IS,IA,NE)
            kpmass_m2(i03+ixy)=AER(i05+ixy)
         enddo
      enddo
      enddo
   endif

   do ii=1,nstep

     dtt0=dt/float(nstep)

     do j=sy(ne)-1,ey(ne)+1
     do i=sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         i04=ip4mem(k,ig,ne)
         I05=IP5MEM(K,IS,IA,NE)
         uu(i,j,k) = u (i03+ixy)
         vv(i,j,k) = v (i03+ixy)
         ww(i,j,k) = w (i03+ixy)
         ddx(i,j,k)= dx(i03+ixy)
         ddz(i,j,k)= dz(i03+ixy)
         conc(i,j,k)=AER(i05+ixy)
         !concmark(i,j,k)=AER(i05+ixy)
       enddo  !k
       i0=ip2mem(ne)
       kktop(i,j) = ktop(i0+ixy)
     enddo   !j
     enddo    !i

     CALL ADV_VERT(MYID,conc, ww,uu,vv,ddz,ddx,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
       i04=ip4mem(k,ig,ne)
       i05=ip5mem(k,is,ia,ne)
       aer(i05+ixy) = conc(i,j,k)
       enddo !k
     enddo  !i
     enddo    !j      


    IF(ia==1) THEN ! for Sea salt

     do iduc = 1, nseacom

      do j=sy(ne)-1,ey(ne)+1
      do i=sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
          i03=ip3mem(k,ne)
          i05c = ip5memcs (k,is,iduc,ne)
          uu(i,j,k) = u (i03+ixy)
          vv(i,j,k) = v (i03+ixy)
          ww(i,j,k) = w (i03+ixy)
          ddx(i,j,k)= dx(i03+ixy)
          ddz(i,j,k)= dz(i03+ixy)
          conc(i,j,k)= SEACOMP(i05c+ixy)
        enddo
        i0=ip2mem(ne)
        kktop(i,j) = ktop(i0+ixy)
      enddo  !i
      enddo    !j

      CALL ADV_VERT (MYID,conc,ww,uu,vv,ddz,ddx,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz
          i05c = ip5memcs (k,is,iduc,ne)
          SEACOMP(i05c+ixy) = AMAX1(conc(i,j,k), 1.E-20)
        enddo !k
      enddo !i
      enddo    !j   

    enddo ! iduc

    ELSE IF (ia==2) then
      do iduc = 1, ndustcom
        do j=sy(ne)-1,ey(ne)+1
        do i=sx(ne)-1,ex(ne)+1
          ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
          do k=1,nzz
            i03=ip3mem(k,ne)
            i05c = ip5memc (k,is,iduc,ne)
            uu(i,j,k) = u (i03+ixy)
            vv(i,j,k) = v (i03+ixy)
            ww(i,j,k) = w (i03+ixy)
            ddx(i,j,k)= dx(i03+ixy)
            ddz(i,j,k)= dz(i03+ixy)
            conc(i,j,k)= DUSTCOMP(i05c+ixy)
          enddo
          i0=ip2mem(ne)
          kktop(i,j) = ktop(i0+ixy)
        enddo   !j
        enddo    !i

        CALL ADV_VERT (MYID,conc,ww,uu,vv,ddz,ddx,sx(ne),ex(ne),sy(ne),ey(ne),nzz,dtt0,ig)

        do j=sy(ne),ey(ne)
        do i=sx(ne),ex(ne)
          ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
          do k=1,nzz
            i05c = ip5memc (k,is,iduc,ne)
            DUSTCOMP(i05c+ixy) = AMAX1(conc(i,j,k), 1.E-20)
          enddo !k
        enddo  !i
        enddo    !j   

      enddo ! iduc


    ENDIF    ! IA IF

    enddo !ii nstep

  ENDDO !IS
  ENDDO ! IA


deallocate(uu,vv,ww,ddz,ddx,conc,kktop)

end subroutine naqpms_v_adv





