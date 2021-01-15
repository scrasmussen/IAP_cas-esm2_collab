
subroutine getvalue_apm( myid,a,sx,ex,sy,ey,ix,iy,value)
  integer :: myid, sx, ex, sy, ey, ix, iy
  real    :: a(sx-1:ex+1,sy-1:ey+1)
  real    :: value
  value=a(ix,iy)
  return
end subroutine getvalue_apm


subroutine putvalue_apm( myid,a,sx,ex,sy,ey,ix,iy,value)
  integer :: myid, sx, ex, sy, ey,ix,iy
  real    :: a(sx-1:ex+1,sy-1:ey+1)
  real    :: value
  a(ix,iy)=value
  return
end


!!!!
subroutine openfile_apm( myid,ne,iyear,imonth,iday,ihour)

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

 call system("mkdir -p out/apm_tmp")
 call system("mkdir -p out/apm_tmp/"//date(1:8))

 !call get_funit(apmfunit)

 close(345+ne)
 open( 345+ne, file='out/apm_tmp/'//date(1:8)// &
       '/apm_'//cdnum//'.'//fname(1:3)//'.'//date,&
         form='unformatted',status='UNKNOWN' )

 !print*,'apm out file opened'

 return

end subroutine openfile_apm
!!!!


subroutine openfile_apm_5min( myid,ne,iyear,imonth,iday,ihour,timestr )

 integer      :: myid
 integer      :: ne
 character*1  :: cdnum
 character*31 :: fname
 character*10 :: date

 character*19 :: timestr

 character*3  :: cid

!apm_d01_000_2011-03-05_00:25:00

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
 write (cid,'(i3.3)') myid

 call system("mkdir -p out/apm_tmp")
 call system("mkdir -p out/apm_tmp/"//date(1:8))

 !call get_funit(apmfunit)

 fname='apm_d0'//cdnum//'_'//cid//'_'//timestr

 print*,'apm output file : ',fname

!stop

 close(345+ne)
 open( 345+ne, file='out/apm_tmp/'//date(1:8)// &
       '/'//trim(fname),&
         form='unformatted',status='UNKNOWN' )

 !print*,'apm out file opened'

 return

end subroutine openfile_apm_5min



subroutine wr_apm_var( myid,lapm &
                      ,iyear,imonth,iday,ihour,iminute &
                      ,nest,nzz,sx,ex,sy,ey,ne  &
                      ,ip2mem,ip3mem,mem2d,mem3d &
                      ,u,v,rh1,clw,rnw,RAINCON,RAINNON ) 

 use apm_varlist
 use aqchem_varlist,  only : so4_aqchem,clw_ph
 implicit none
 include 'apm_parm.inc'

 logical :: lapm
 integer :: myid
 integer :: ne,nest
 integer :: nx(5),ny(5),nzz
 integer :: sy(5),ey(5),sx(5),ex(5)

 integer :: i,j,k,is
 integer :: mem2d,mem3d

 integer :: ixy,i02,i03,iapm,i0
 integer :: ip3mem(nzz,nest),ip2mem(nest)

 integer :: kk

 real    :: u(mem3d),v(mem3d),clw(mem3d),rnw(mem3d),rh1(mem3d)
 real    :: RAINCON(mem2d),RAINNON(mem2d)

 integer :: iyear,imonth,iday,ihour,iminute
 character*1  :: cdnum
 character*20 :: fname
 character*10 :: date
 character*3  :: cmyid
 integer      :: funit

 write(date(1:4),'(i4)')iyear
 write(date(5:6),'(i2.2)')imonth
 write(date(7:8),'(i2.2)')iday
 write(date(9:10),'(i2.2)')ihour

 write(cdnum(1:1),'(i1)')ne

! print*,'wr myid=',myid

 if ( myid .lt. 0 ) then
     close( 11 )
     return
 endif

!apm_1.007.2007032022

 call system("mkdir -p out/apm_tmp")
 call system("mkdir -p out/apm_tmp/"//date(1:8))

 write (cmyid,'(i3.3)') myid

! write (fname(1:3),'(i3.3)') myid

 fname='apm_'//cdnum//'.'//cmyid//'.'//date

! print*,'output_fname : ',trim(fname)

 call get_funit(funit)
 open(funit,file='out/apm_tmp/'//date(1:8)//'/'//trim(fname) &
           ,form='unformatted',status='UNKNOWN' )


 if(myid.eq.0) print*,'write apm data :',myid,sx(ne),ex(ne),sy(ne),ey(ne)

 write(funit) myid,sx(ne),ex(ne),sy(ne),ey(ne)
 !print*, 'dims=',myid,sx(ne),ex(ne),sy(ne),ey(ne),nzz
 !print*,'so4 bins=',NSO4

!goto 34567

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,rh1(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,clw(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,rnw(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,so4_aqchem(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,clw_ph(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,cn10nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,npf3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,h2so4_gas(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,apm_cacid(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,apm_pacid(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,rgf_sulf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,rgf_salt(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,rgf_dust(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,rgf_bc(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,rgf_oc(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,apm_xq3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

!cn3nm
  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,cn3nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,spgf_3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,bcagt(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,ocagt(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,bcgetsp_rate(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call write2d_v2(myid,ocgetsp_rate(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  i0=ip2mem(ne)
  call write2d_v2(myid,RAINCON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,RAINNON(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,vsblt2d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,wrfvsb(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,real(iffog(i0)),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,real(ifhaze(i0)),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,real(ifwrffog(i0)),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AOD_390nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AOD_500nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AOD_530nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AOD_550nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AOD_700nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AOD_1010nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AAOD_390nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AAOD_500nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AAOD_530nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AAOD_550nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AAOD_700nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  i0=ip2mem(ne)
  call write2d_v2(myid,AAOD_1010nm(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)

  loop_sulf : do is=1,NSO4
  do k=1,nzz
    i0=ip_sulf(k,is,ne)
    call write2d_v2(myid,apm_sulf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo loop_sulf

  loop_salt : do is=1,NSEA
  do k=1,nzz
    i0=ip_salt(k,is,ne)
    call write2d_v2(myid,apm_salt(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo loop_salt

  loop_dust : do is=1,NDSTB
  do k=1,nzz
    i0=ip_dust(k,is,ne)
    call write2d_v2(myid,apm_dust(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo loop_dust

  loop_bcoc : do is=1,NBCOCT
  do k=1,nzz
    i0=ip_bcoc(k,is,ne)
    call write2d_v2(myid,apm_bcoc(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo loop_bcoc

  do k=1,nzz
    i0=ip3mem(k,ne)
    call write2d_v2(myid,msltsulf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
    i0=ip3mem(k,ne)
    call write2d_v2(myid,mdstsulf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
    i0=ip3mem(k,ne)
    call write2d_v2(myid,mbcsulf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
    i0=ip3mem(k,ne)
    call write2d_v2(myid,mocsulf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
    i0=ip3mem(k,ne)
    call write2d_v2(myid,bulk_msp(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo

  do k=1,nzz
    i0=ip3mem(k,ne)
    call write2d_v2(myid,tot_sulf(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo


  do is=1,5
    do k=1,nzz
      i0=ip_type(k,is,ne)
      call write2d_v2(myid,number_ccn1(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
    enddo
  enddo

  do is=1,5
    do k=1,nzz
      i0=ip_type(k,is,ne)
      call write2d_v2(myid,number_ccn2(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
    enddo
  enddo

  do is=1,5
    do k=1,nzz
      i0=ip_type(k,is,ne)
      call write2d_v2(myid,number_ccn3(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
    enddo
  enddo

  do is=1,5
    do k=1,nzz
      i0=ip_type(k,is,ne)
      call write2d_v2(myid,ztn3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
    enddo
  enddo

  do is=1,5
    do k=1,nzz
      i0=ip_type(k,is,ne)
      call write2d_v2(myid,ycs3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
    enddo
  enddo

  do is=1,5
    do k=1,nzz
      i0=ip_type(k,is,ne)
      call write2d_v2(myid,surf3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
    enddo
  enddo

  do is=1,5
    do k=1,nzz
      i0=ip_type(k,is,ne)
      call write2d_v2(myid,rgfdry3d(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
    enddo
  enddo

  do is=1,5
    i0=ip_2dtype(is, ne)
    call write2d_v2(myid,TAOD(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo


  do is=1,nbincb
  do k=1,nzz
    i0=ip_cbbin(k,is,ne)
    call write2d_v2(myid,apm_binbc(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit) 
  enddo
  enddo

  do is=1,nbincb
  do k=1,nzz
    i0=ip_cbbin(k,is,ne)
    call write2d_v2(myid,apm_binoc(i0),sx(ne),ex(ne),sy(ne),ey(ne),k,ne,funit)
  enddo
  enddo


 close(funit)

!
! 201406042albany
! clean up deposition mass

i02 = ip2mem(ne)
do j = sy(ne),ey(ne)
do i = sx(ne),ex(ne)

  ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
  ddep_apm_spsulf(i02+ixy)=0.0
  ddep_apm_ppsulf(i02+ixy)=0.0
  ddep_apm_sulf(i02+ixy)=0.0

enddo
enddo



end subroutine wr_apm_var








subroutine openfile_apm_bph( myid,ne,iyear,imonth,iday,ihour,timestr )

 integer      :: myid
 integer      :: ne
 character*1  :: cdnum
 character*31 :: fname
 character*10 :: date

 character*19 :: timestr

 character*3  :: cid

!apm_d01_000_2011-03-05_00:25:00

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
 write (cid,'(i3.3)') myid

 call system("mkdir -p out/apm_tmp.bph")
 call system("mkdir -p out/apm_tmp.bph/"//date(1:8))

 !call get_funit(apmfunit)


 fname='apm_d0'//cdnum//'_'//cid//'_'//timestr

 print*,'apm output file : ',fname

!stop

 close(345+ne)
 open( 345+ne, file='out/apm_tmp.bph/'//date(1:8)// &
       '/'//trim(fname),&
         form='unformatted',status='UNKNOWN' )

 !print*,'apm out file opened'

 return

end subroutine openfile_apm_bph
