
subroutine apm_shun_check &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip3mem,mem3d &
 &  ,prgname )

use apm_varlist
implicit none
include 'apm_parm.inc'

character*25 :: prgname

integer :: myid

real :: dt

logical :: lapm

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is
!integer :: k,is

integer :: mem3d

integer :: ixy,i02,i03,iapm

integer :: ip3mem(nzz,nest)

!real,parameter :: upper=1.0e4,lower=0

real :: tmp


!print*,sy(ne)-1,ey(ne)+1,sx(ne)-1,ex(ne)+1
!print*,ne

!stop

   loop_so4_hdif : do is=1,NSO4

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          iapm=ip_sulf(k,is,ne)
          tmp=apm_sulf(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            print*
            print*,'var : ','so4'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'is=',is
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_so4_hdif



   loop_salt_hdif : do is=1,NSEA

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          iapm=ip_salt(k,is,ne)
          tmp=apm_salt(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            print*
            print*,'var : ','salt'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'is=',is
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_salt_hdif




   loop_dust_hdif : do is=1,NDSTB

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          iapm=ip_dust(k,is,ne)
          tmp=apm_dust(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            print*
            print*,'var : ','dust'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'is=',is
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_dust_hdif





   loop_bcoc_hdif : do is=1,NBCOCT

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          tmp=apm_bcoc(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            print*
            print*,'var : ','bcoc'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'is=',is
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_bcoc_hdif


!==========================================================================
!==========================================================================

! apm coated species

!-> sulfate on seasalt
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=msltsulf(i03+ixy)
          if(tmp.lt.lower.or.tmp.gt.upper) then
            print*
            print*,'var : ','sltsulf'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
          endif
       enddo
    enddo
    enddo

!-> sulfate on dust
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=mdstsulf(i03+ixy)
         if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            print*
            print*,'var : ','dstsulf'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
          endif
       enddo
    enddo
    enddo

!-> sulfate on BC
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=mbcsulf(i03+ixy)
         if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            print*
            print*,'var : ','bcsulf'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
         endif
       enddo
    enddo
    enddo

!-> sulfate on POC
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=mocsulf(i03+ixy)
         if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            print*
            print*,'var : ','ocsulf'
            print*,'erro after ',trim(prgname)
            print*,'location:'
            print*,'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            print*,'k =',k
            print*,'j =',j
            print*,'i =',i
            print*,'value=',tmp
            stop 'shun_check stop'
         endif
       enddo
    enddo
    enddo

!ENDIF ! apm flag


end subroutine apm_shun_check





