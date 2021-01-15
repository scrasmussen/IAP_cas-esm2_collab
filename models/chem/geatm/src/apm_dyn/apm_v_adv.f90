
subroutine apm_v_adv &
 & ( myid &
 &  ,lapm,imasskeep &
 &  ,ne,dt,nstep &
 &  ,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dx,dy,dz &
 &  ,u,v,w &
 &  ,ktop &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real :: dt
integer :: nstep
real :: dtt0

logical :: lapm
integer :: imasskeep

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: ii
integer :: i,j,k,is

integer :: iwb,ieb,jsb,jeb
real,allocatable,dimension(:,:) :: kktop
real,allocatable,dimension(:,:,:) :: uu,vv,ww,conc,ddx,ddz

integer :: mem2d,mem3d

real,dimension(mem3d) :: dx,dy,u,v,kpmass_m2

real,dimension(mem3d) :: dz,w

real,dimension(mem2d) :: ktop

integer :: ixy,i03,i02,iapm

integer :: ip2mem(nest)
integer :: ip3mem(nzz,nest)



!print*,'w_ws=',w

!print*,'nstep=',nstep ! print ok(8.)

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


!print*,'vdav dim',iwb,ieb,jsb,jeb

!print*,'ktop=',ktop

!return

!IF(lapm) THEN ! apm flag

   do j=sy(ne)-1,ey(ne)+1
   do i=sx(ne)-1,ex(ne)+1
       !print*,'k-j-i',j,i
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         iapm=ip_sulf(k,is,ne)
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




 !> shun : apm sulfate vadv
  if(lfor_sulf) then
   !print*,'sulf vadv'
   loop_so4_vadv : do is=1,NSO4

    if(imasskeep==1) then
       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)
          ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
          do k=1,nzz
            i03=ip3mem(k,ne)
            iapm=ip_sulf(k,is,ne)
            kpmass_m2(i03+ixy)=apm_sulf(iapm+ixy)
!            wk(i03+ixy)=apm_sulf(iapm+ixy)
          enddo
       enddo
       enddo
    endif

    do ii=1,nstep
     dtt0=dt/float(nstep)
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)

       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         iapm=ip_sulf(k,is,ne)
         conc(i,j,k)=apm_sulf(iapm+ixy)
       enddo
     enddo
     enddo

!if(sy(ne).le.30.and.30.le.ey(ne).and.sx(ne).le.30.and.30.le.ex(ne).and..false.) then
!print*,'inputs'
!print*,'ww=',ww(30,30,:)
!print*,'uu=',uu(30,30,:)
!print*,'vv=',vv(30,30,:)
!print*,'ddx=',ddx(30,30,:)
!print*,'ddz=',ddz(30,30,:)


!print*,'0000'
!print*,conc(30,30,:)
!endif


     CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                  & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                  & ,nzz,dtt0,is )

!if(sy(ne).le.30.and.30.le.ey(ne).and.sx(ne).le.30.and.30.le.ex(ne).and..false.) then
!print*,'1111'
!print*,conc(30,30,:)

!stop 'kkconc'
!endif

     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         iapm=ip_sulf(k,is,ne)
         apm_sulf(iapm+ixy)=conc(i,j,k)
       enddo
     enddo
     enddo
    enddo !! ii

   enddo loop_so4_vadv
  endif
 !< shun : end of apm sulfate vadv


!return


 !> shun : apm seasalt vadv
  if(lfor_salt) then
   !print*,'salt vadv'
   loop_salt_vadv : do is=1,NSEA

    if(imasskeep==1) then
       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)
          ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
          do k=1,nzz-1
            i03=ip3mem(k,ne)
            iapm=ip_salt(k,is,ne)
            kpmass_m2(i03+ixy)=apm_salt(iapm+ixy)
          enddo
       enddo
       enddo
    endif

    do ii=1,nstep
     dtt0=dt/float(nstep)
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         iapm=ip_salt(k,is,ne)
         conc(i,j,k)=apm_salt(iapm+ixy)
       enddo
     enddo
     enddo
     CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                  & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                  & ,nzz,dtt0,is )
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         iapm=ip_salt(k,is,ne)
         apm_salt(iapm+ixy)=conc(i,j,k)
       enddo
     enddo
     enddo
    enddo !! ii

   enddo loop_salt_vadv
  endif
 !< shun : end of apm seasalt vadv

 !> shun : apm dust vadv
  if(lfor_dust) then
   !print*,'dust vadv'
   loop_dust_vadv : do is=1,NDSTB

    if(imasskeep==1) then
       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)
          ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
          do k=1,nzz-1
            i03=ip3mem(k,ne)
            iapm=ip_dust(k,is,ne)
            kpmass_m2(i03+ixy)=apm_dust(iapm+ixy)
          enddo
       enddo
       enddo
    endif

    do ii=1,nstep
     dtt0=dt/float(nstep)
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         iapm=ip_dust(k,is,ne)
         conc(i,j,k)=apm_dust(iapm+ixy)
       enddo
     enddo
     enddo
     CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                  & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                  & ,nzz,dtt0,is )
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         iapm=ip_dust(k,is,ne)
         apm_dust(iapm+ixy)=conc(i,j,k)
       enddo
     enddo
     enddo
    enddo !! ii

   enddo loop_dust_vadv
  endif
 !< shun : end of apm dust vadv

 !> shun : apm bcoc vadv
  if(lfor_bcoc) then
   !print*,'bcoc vadv'
   loop_bcoc_vadv : do is=1,NBCOCT

    if(imasskeep==1) then
       do j = sy(ne),ey(ne)
       do i = sx(ne),ex(ne)
          ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
          do k=1,nzz-1
            i03=ip3mem(k,ne)
            iapm=ip_bcoc(k,is,ne)
            kpmass_m2(i03+ixy)=apm_bcoc(iapm+ixy)
          enddo
       enddo
       enddo
    endif

    do ii=1,nstep
     dtt0=dt/float(nstep)
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         iapm=ip_bcoc(k,is,ne)
         conc(i,j,k)=apm_bcoc(iapm+ixy)
       enddo
     enddo
     enddo
     CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                  & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                  & ,nzz,dtt0,is )
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         iapm=ip_bcoc(k,is,ne)
         apm_bcoc(iapm+ixy)=conc(i,j,k)
       enddo
     enddo
     enddo
    enddo !! ii

   enddo loop_bcoc_vadv
  endif
 !< shun : end of apm bcoc vadv


!=======================================================
!=======================================================

! apm coated species

IF(lcoated_dyn) THEN
!-> sulfate coated on seasalt
   !print*,'coated vadv'
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=msltsulf(i03+ixy)
         endif
      enddo
   enddo
   enddo

   do ii=1,nstep
    dtt0=dt/float(nstep)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        conc(i,j,k)=msltsulf(i03+ixy)
      enddo
    enddo
    enddo
    CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                 & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                 & ,nzz,dtt0,is )
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        msltsulf(i03+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo !! ii

!-> sulfate coated on dust
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mdstsulf(i03+ixy)
         endif
      enddo
   enddo
   enddo

   do ii=1,nstep
    dtt0=dt/float(nstep)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        conc(i,j,k)=mdstsulf(i03+ixy)
      enddo
    enddo
    enddo
    CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                 & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                 & ,nzz,dtt0,is )
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        mdstsulf(i03+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo !! ii

!-> sulfate coated on BC
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mbcsulf(i03+ixy)
         endif
      enddo
   enddo
   enddo

   do ii=1,nstep
    dtt0=dt/float(nstep)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        conc(i,j,k)=mbcsulf(i03+ixy)
      enddo
    enddo
    enddo
    CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                 & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                 & ,nzz,dtt0,is )
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        mbcsulf(i03+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo !! ii


!-> sulfate coated on POC
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mocsulf(i03+ixy)
         endif
      enddo
   enddo
   enddo

   do ii=1,nstep
    dtt0=dt/float(nstep)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        conc(i,j,k)=mocsulf(i03+ixy)
      enddo
    enddo
    enddo
    CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                 & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                 & ,nzz,dtt0,is )
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        mocsulf(i03+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo !! ii

ENDIF ! lcoated_dyn


! sulferic acid vapor 
 if(lfor_h2so4) then
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=h2so4_gas(i03+ixy)
         endif
      enddo
   enddo
   enddo

   do ii=1,nstep
    dtt0=dt/float(nstep)
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        i03=ip3mem(k,ne)
        conc(i,j,k)=h2so4_gas(i03+ixy)
      enddo
    enddo
    enddo
    CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                 & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                 & ,nzz,dtt0,is )
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz
        h2so4_gas(i03+ixy)=conc(i,j,k)
      enddo
    enddo
    enddo
   enddo !! ii
 endif
!!!!!!!!!!


!stop 'over vadvsub'

!ENDIF ! apm flag

   loop_binbc : do is=1,nbincb

    do ii=1,nstep
     dtt0=dt/float(nstep)
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         iapm=ip_cbbin(k,is,ne)
         conc(i,j,k)=apm_binbc(iapm+ixy)
       enddo
     enddo
     enddo
     CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                  & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                  & ,nzz,dtt0,is )
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         iapm=ip_cbbin(k,is,ne)
         apm_binbc(iapm+ixy)=conc(i,j,k)
       enddo
     enddo
     enddo
    enddo !! ii

   enddo loop_binbc

   loop_binoc : do is=1,nbincb

    do ii=1,nstep
     dtt0=dt/float(nstep)
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         i03=ip3mem(k,ne)
         iapm=ip_cbbin(k,is,ne)
         conc(i,j,k)=apm_binoc(iapm+ixy)
       enddo
     enddo
     enddo
     CALL ADV_VERT( MYID,conc,ww,uu,vv,ddz,ddx &
                  & ,sx(ne),ex(ne),sy(ne),ey(ne) &
                  & ,nzz,dtt0,is )
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz
         iapm=ip_cbbin(k,is,ne)
         apm_binoc(iapm+ixy)=conc(i,j,k)
       enddo
     enddo
     enddo
    enddo !! ii

   enddo loop_binoc




end subroutine apm_v_adv





