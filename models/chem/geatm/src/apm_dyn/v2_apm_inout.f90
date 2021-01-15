
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

 integer :: ixy,i02,i03,iapm
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

 ! print*,'write wind speed, cloud and rain water content'
 do k=1,nzz
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      iapm=ip3mem(k,ne)
      write(funit) &
                   ,u(iapm+ixy),v(iapm+ixy) &
                   ,rh1(iapm+ixy) &
                   ,clw(iapm+ixy),rnw(iapm+ixy) &
                   ,so4_aqchem(iapm+ixy) &
                   ,clw_ph(iapm+ixy) &
!                   ,number_ccn1(iapm+ixy),number_ccn2(iapm+ixy) &
!                   ,number_ccn3(iapm+ixy) &
                   ,cn10nm(iapm+ixy) &
                   ,npf3d(iapm+ixy)  &
                   ,h2so4_gas(iapm+ixy) &
                   ,apm_cacid(iapm+ixy),apm_pacid(iapm+ixy)  &
                   ,rgf_sulf(iapm+ixy),rgf_salt(iapm+ixy) &
                   ,rgf_dust(iapm+ixy),rgf_bc(iapm+ixy) &
                   ,rgf_oc(iapm+ixy) &
                   ,apm_xq3d(iapm+ixy) &
                   ,cn3nm(iapm+ixy) &
                   ,spgf_3d(iapm+ixy) &
                   ,bcagt(iapm+ixy) &
                   ,ocagt(iapm+ixy) &
                   ,bcgetsp_rate(iapm+ixy) &
                   ,ocgetsp_rate(iapm+ixy)
   enddo
   enddo
 enddo

 ! print*,'write cloud and rain water content'
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      iapm=ip2mem(ne)
      write(funit) &
                   ,RAINCON(iapm+ixy),RAINNON(iapm+ixy) &
                   ,vsblt2d(iapm+ixy),wrfvsb(iapm+ixy)  &
                   ,real(iffog(iapm+ixy)),real(ifhaze(iapm+ixy)) &
                   ,real(ifwrffog(iapm+ixy)) &
                   ,AOD_390nm(iapm+ixy)   &
                   ,AOD_500nm(iapm+ixy)   &
                   ,AOD_530nm(iapm+ixy)   &
                   ,AOD_550nm(iapm+ixy)   &
                   ,AOD_700nm(iapm+ixy)   &
                   ,AOD_1010nm(iapm+ixy)  &
                   ,AAOD_390nm(iapm+ixy)  &
                   ,AAOD_500nm(iapm+ixy)  &
                   ,AAOD_530nm(iapm+ixy)  &
                   ,AAOD_550nm(iapm+ixy)  &
                   ,AAOD_700nm(iapm+ixy)  &
                   ,AAOD_1010nm(iapm+ixy)  

!      if(i.eq.33.and.j.eq.33) then      
!        print*,'vis_33&33=',vsblt2d(iapm+ixy)
!      endif
   enddo
   enddo


!print*,'write sulfate'
 loop_sulf : do is=1,NSO4
   do k=1,nzz
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
        ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_sulf(k,is,ne)
        !print*,k,j,i,apm_sulf(iapm+ixy)
        write(funit) apm_sulf(iapm+ixy)
        !kk=kk+1
        !print*,k,j,i
        !print*,kk
     enddo
     enddo
   enddo 
 enddo loop_sulf

!print*,'write salt'
 loop_salt : do is=1,NSEA
   do k=1,nzz
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
        ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_salt(k,is,ne)
        write(funit) apm_salt(iapm+ixy)
     enddo
     enddo
   enddo
 enddo loop_salt

!print*,'write dust'
 loop_dust : do is=1,NDSTB
   do k=1,nzz
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
        ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_dust(k,is,ne)
        write(funit) apm_dust(iapm+ixy)
     enddo
     enddo
   enddo
 enddo loop_dust

!print*,'write bcoc'
 loop_bcoc : do is=1,NBCOCT
   do k=1,nzz
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
        ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_bcoc(k,is,ne)
        write(funit) apm_bcoc(iapm+ixy)
     enddo
     enddo
   enddo
 enddo loop_bcoc

! print*,'write emit data'

! do k=1,nzz
!   do j=sy(ne),ey(ne)
!   do i=sx(ne),ex(ne)
!      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!      iapm=ip3mem(k,ne)
!      !if(k.eq.1) print*,sulf_emit(iapm+ixy),bc_emit(iapm+ixy),oc_emit(iapm+ixy)
!      write(345+ne) k,j,i &
!                   ,sulf_emit(iapm+ixy),bc_emit(iapm+ixy),oc_emit(iapm+ixy)
!   enddo
!   enddo
! enddo


!34567 continue

 !print*,'write salt emit'
! do is=1,NSEA
!   do k=1,nzz
!     do j=sy(ne),ey(ne)
!     do i=sx(ne),ex(ne)
!        ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!        iapm=ip_salt(k,is,ne)
!        write(345+ne) is,k,j,i,salt_emit(iapm+ixy)
!     enddo
!     enddo
!   enddo
! enddo

 !print*,'write dust emit'
! do k=1,nzz
!   do j=sy(ne),ey(ne)
!   do i=sx(ne),ex(ne)
!      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!      iapm=ip3mem(k,ne)
!      write(345+ne) k,j,i,dust_emit(iapm+ixy)
!   enddo
!   enddo
! enddo

 !print*,'write H2SO4 #/cm3'
! do k=1,nzz
!   do j=sy(ne),ey(ne)
!   do i=sx(ne),ex(ne)
!      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
!      iapm=ip3mem(k,ne)
!      write(345+ne) k,j,i,apm_cacid(iapm+ixy),apm_pacid(iapm+ixy),apm_xq3d(iapm+ixy)
!   enddo
!   enddo
! enddo

 !print*,'write coated sulfate'
 do k=1,nzz
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      iapm=ip3mem(k,ne)
      write(funit) msltsulf(iapm+ixy),mdstsulf(iapm+ixy) &
                   ,mbcsulf(iapm+ixy),mocsulf(iapm+ixy),bulk_msp(iapm+ixy) &
                   ,tot_sulf(iapm+ixy)
   enddo
   enddo
 enddo

 loop_type : do is=1,5
   do k=1,nzz
     do j=sy(ne),ey(ne)
     do i=sx(ne),ex(ne)
        ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        iapm=ip_type(k,is,ne)
        write(funit) number_ccn1(iapm+ixy),number_ccn2(iapm+ixy) &
                     ,number_ccn3(iapm+ixy),ztn3d(iapm+ixy) &
                     ,ycs3d(iapm+ixy),surf3d(iapm+ixy) &
                     ,rgfdry3d(iapm+ixy)
     enddo
     enddo
   enddo
 enddo loop_type


 do is=1,5
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     iapm=ip_2dtype(is, ne)
     write(funit) TAOD(iapm+ixy)
   enddo
   enddo
 enddo


 do is=1,nbincb
  do k=1,nzz
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
     ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     iapm=ip_cbbin(k,is,ne)
     write(funit) apm_binbc(iapm+ixy),apm_binoc(iapm+ixy)
   enddo
   enddo
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
