
subroutine naqpms_rd_emfl &
 & ( myid &
 &  ,naqpms_dir &
 &  ,iyear,imonth,iday,ihour &
 &  ,lnaqpms_ems,lfrac_so2_emit &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,NEMIT,NLAY_EM,IPIG,emsfrc &
 &  ,igas,iaer,isize,nseacom,ndustcom )

use naqpms_varlist
use naqpms_gridinfo

implicit none

integer :: myid

integer :: iyear,imonth,iday,ihour

character    :: naqpms_dir*500

logical :: lapm,lnaqpms_ems,lfrac_so2_emit

real :: dt

integer :: ne,nest
integer :: irecg(5)
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: i,j,k,is
integer :: ictg

integer :: iemittype

integer :: i0,i02,i03,ig,img,i02emt
integer :: ixy

integer  ::  NEMIT,NLAY_EM
integer,dimension(NEMIT) :: IPIG

real :: emsfrc(NEMIT,NLAY_EM)

character :: ctime*100,fname*100
character :: cdnum*1
character :: date*10

integer :: funit,irec
logical :: lexist

real,dimension(nzz) :: ratioem1,ratioem2,ratioem3,ratioem4,ratioem5,ratioem6


!print*,'dim=',sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne)

 write(date(1:4),'(i4)')iyear
 write(date(5:6),'(i2.2)')imonth
 write(date(7:8),'(i2.2)')iday
 write(date(9:10),'(i2.2)')ihour

 write(cdnum(1:1),'(i1)') ne

 ctime=date(5:6)

 call get_funit(funit)

 fname='emit/data.emit/emitgrid_'//trim(ctime)//'.d'//cdnum

 if(myid.eq.0) then
!   print*,'read '//trim(fname)
!   print*
 endif

  inquire(file=fname,exist=lexist)
  if(.not.lexist) then
    print*,trim(fname)//' NOT exist'
    stop
  endif


 open(funit,file=trim(fname),form='unformatted',access='direct' &
           ,recl=nx(ne)*ny(ne),status='old')

 irec=1

 do img=1,NEMIT!  NEMIT

      ig= IPIG(img)

      if(ig.eq.57) then ! skip dms
        irec=irec+NLAY_EM
        cycle
      endif

      if(ig>0) then
       !i0=ip2memGas(ig,ne)
       !call read2d(myid,EmtaGas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
       !        nx(ne),ny(ne),irec,funit)
        do ictg=1,NLAY_EM

          i02emt=ip_emit2d(ig,ictg,ne)
          call read2d(myid,emit2d(i02emt),sx(ne),ex(ne),sy(ne),ey(ne) &
                     ,nx(ne),ny(ne),irec,funit)

          if(lnaqpms_ems) then
             do j=sy(ne),ey(ne)
             do i=sx(ne),ex(ne)
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
               emit2d(i02emt+ixy)=emit2d(i02emt+ixy)*emsfrc(img,ictg)
             enddo
             enddo
          endif

        enddo
      else
         print*,'ig for emit err'
!          irecg(ne)=irecg(ne)+3       
      endif

 enddo !img


 close(funit)



end subroutine naqpms_rd_emfl





