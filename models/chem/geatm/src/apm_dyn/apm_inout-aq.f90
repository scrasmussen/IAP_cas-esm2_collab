
!!!!
subroutine openfile_aqf( myid,ne,iyear,imonth,iday,ihour)

 integer      :: myid
 integer      :: ne
 character*1  :: cdnum
 character*20 :: fname
 character*10 :: date

 write(date(1:4),'(i4)')iyear
 write(date(5:6),'(i2.2)')imonth
 write(date(7:8),'(i2.2)')iday
 write(date(9:10),'(i2.2)')ihour

 print*,'wr myid=',myid

 if ( myid .lt. 0 ) then
     close( 11 )
     return
 endif

 write(cdnum(1:1),'(i1)')ne
 write (fname(1:3),'(i3.3)') myid

 call system("mkdir -p out/aq_tmp")
 call system("mkdir -p out/aq_tmp/"//date(1:8))

 !call get_funit(apmfunit)

 close(445+ne)
 open( 445+ne, file='out/aq_tmp/'//date(1:8)// &
       '/aq_'//cdnum//'.'//fname(1:3)//'.'//date,&
         form='unformatted',status='UNKNOWN' )

 !print*,'apm out file opened'

 return

end subroutine openfile_aqf
!!!!



subroutine wr_aq_var( myid,lapm,nest,nzz,sx,ex,sy,ey,ne  &
                     ,ip2mem,ip3mem,mem2d,mem3d &
                     ,clw_3d )

 use apm_varlist
 use aqchem_varlist
 implicit none
 include 'apm_parm.inc'

 logical :: lapm
 integer :: myid
 integer :: ne,nest
 integer :: nx(5),ny(5),nzz
 integer :: sy(5),ey(5),sx(5),ex(5)

 integer :: i,j,k,is
 integer :: mem2d,mem3d

 integer :: ixy,i02,i03,iapm
 integer :: ip3mem(nzz,nest),ip2mem(nest)

 integer :: kk

 real    :: clw_3d(mem3d)

 print*,'write aqchem data :',myid,sx(ne),ex(ne),sy(ne),ey(ne)
 write(445+ne) myid,sx(ne),ex(ne),sy(ne),ey(ne)
 !print*, 'dims=',myid,sx(ne),ex(ne),sy(ne),ey(ne),nzz
 !print*,'so4 bins=',NSO4

!goto 34567

 ! print*,'write wind speed, cloud and rain water content'
 do k=1,nzz
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      iapm=ip3mem(k,ne)
      write(445+ne)  &
                   ,clw_3d(iapm+ixy) &
                   ,so4_aqchem(iapm+ixy) &
                   ,clw_ph(iapm+ixy) &
                   ,so4_ccn2(iapm+ixy) &
                   ,salt_ccn2(iapm+ixy) &
                   ,dust_ccn2(iapm+ixy) &
                   ,bc_ccn2(iapm+ixy) &
                   ,oc_ccn2(iapm+ixy)
   enddo
   enddo
 enddo

 loop_fccn : do is = 1,NSO4+4
  do k=1,nzz
    do j=sy(ne),ey(ne)
    do i=sx(ne),ex(ne)
       ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       iapm=ip_fccn (k,is, ne)
       write(445+ne)  &
                    ,frac_ccn(iapm+ixy)
    enddo
    enddo
  enddo
 enddo loop_fccn

 close(445+ne)

end subroutine wr_aq_var

