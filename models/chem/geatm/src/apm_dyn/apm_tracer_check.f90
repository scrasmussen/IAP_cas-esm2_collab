
subroutine apm_shun_check &
 & ( myid &
 &  ,lapm &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,ip3mem,mem3d &
 &  ,prgname )

use apm_varlist
use cputime
implicit none
include 'apm_parm.inc'

character(len=*) :: prgname

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

character :: cmyid*3
integer   :: funit

logical   :: lerr


!print*,sy(ne)-1,ey(ne)+1,sx(ne)-1,ex(ne)+1
!print*,ne

!stop

   lerr=.false.

   write(cmyid,'(i3.3)') myid


   loop_so4_hdif : do is=1,NSO4

!     do j = sy(ne)-1,ey(ne)+1
!     do i = sx(ne)-1,ex(ne)+1
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)

        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          iapm=ip_sulf(k,is,ne)
          tmp=apm_sulf(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.apm_sulf.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','so4'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_so4_hdif



   loop_salt_hdif : do is=1,NSEA

!     do j = sy(ne)-1,ey(ne)+1
!     do i = sx(ne)-1,ex(ne)+1
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          iapm=ip_salt(k,is,ne)
          tmp=apm_salt(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.apm_salt.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','salt'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_salt_hdif




   loop_dust_hdif : do is=1,NDSTB

!     do j = sy(ne)-1,ey(ne)+1
!     do i = sx(ne)-1,ex(ne)+1
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          iapm=ip_dust(k,is,ne)
          tmp=apm_dust(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.apm_dust.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','dust'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_dust_hdif





   loop_bcoc_hdif : do is=1,NBCOCT

!     do j = sy(ne)-1,ey(ne)+1
!     do i = sx(ne)-1,ex(ne)+1
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          tmp=apm_bcoc(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.apm_bcoc.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','bcoc'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_bcoc_hdif


!==========================================================================
!==========================================================================

! apm coated species

!-> sulfate on seasalt
!    do j = sy(ne)-1,ey(ne)+1
!    do i = sx(ne)-1,ex(ne)+1
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=msltsulf(i03+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.msltsulf.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','sltsulf'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
       enddo
    enddo
    enddo

!-> sulfate on dust
!    do j = sy(ne)-1,ey(ne)+1
!    do i = sx(ne)-1,ex(ne)+1
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=mdstsulf(i03+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.mdstsulf.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','mdstsulf'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
       enddo
    enddo
    enddo

!-> sulfate on BC
!    do j = sy(ne)-1,ey(ne)+1
!    do i = sx(ne)-1,ex(ne)+1
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=mbcsulf(i03+ixy)
         if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.bcsulf.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','bcsulf'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
         endif
       enddo
    enddo
    enddo

!-> sulfate on POC
!    do j = sy(ne)-1,ey(ne)+1
!    do i = sx(ne)-1,ex(ne)+1
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         tmp=mocsulf(i03+ixy)
         if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.ocsulf.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','mocsulf'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
         endif
       enddo
    enddo
    enddo

   loop_binbc : do is=1,nbincb

!     do j = sy(ne)-1,ey(ne)+1
!     do i = sx(ne)-1,ex(ne)+1
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_cbbin(k,is,ne)
          tmp=apm_binbc(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.apm_binbc.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','binbc'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_binbc


   loop_binoc : do is=1,nbincb

!     do j = sy(ne)-1,ey(ne)+1
!     do i = sx(ne)-1,ex(ne)+1
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_cbbin(k,is,ne)
          tmp=apm_binoc(iapm+ixy)
          if(.not.(tmp.gt.lower.and.tmp.le.upper)) then
            call get_funitnaqpms(funit)
            open(funit,file='out/rsl.'//cmyid//'.apm_binbc.'//trim(prgname) &
                      ,position='append')
            write(funit,*)
            write(funit,*) my_year,my_month,my_day,my_hour,my_minute
            write(funit,*) 'var : ','binbc'
            write(funit,*) 'erro after ',trim(prgname)
            write(funit,*) 'location:'
            write(funit,*) 'domain:',sx(ne),ex(ne),sy(ne),ey(ne)
            write(funit,*) 'is=',is
            write(funit,*) 'k =',k
            write(funit,*) 'j =',j
            write(funit,*) 'i =',i
            write(funit,*) 'value=',tmp
            close(funit)
            lerr=.true.
!            stop 'shun_check stop'
          endif
        enddo
     enddo
     enddo

   enddo loop_binoc


   if(lerr) then
     stop 'apm_check stop'
   endif


!ENDIF ! apm flag


end subroutine apm_shun_check





