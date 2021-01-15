
subroutine apm_rst( myid,lapm,ne,nx,ny,nzz,nest,sy,ey,sx,ex,ip3mem &
                   ,naqpms_dir,iyear,imonth,iday,ihour)
use apm_varlist !, only : apm_sulf,apm_salt,apm_dust,apm_bcoc
implicit none
include 'apm_parm.inc'
integer :: myid
logical :: lapm
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)
integer :: i,j,k,is
integer :: iapm,ixy,iccn,i03
integer :: ip3mem(nzz,nest)

integer :: imode

character :: naqpms_dir*500,rstfile*500
character :: rst_time*19,cdd*3
integer   :: iyear,imonth,iday,ihour
integer   :: iminute=0,isecond=0
integer   :: funit
integer   :: id
integer   :: irec
!equivalence(id,ne)

real      :: var2d(nx(ne),ny(ne)) 

character(len=*),parameter :: mbar='-',bbar='_',colon=':'

real :: so4tmp
real,parameter :: kapa(1:5)=(/0.90,1.28,0.0,0.0,0.1/)

real,allocatable,dimension(:,:,:) :: ttbc,ttoc,ttdust,ttsalt

!apm_sulf=1.0d0
!return
!print*,'kk_init',myid,ne
!print*,'kk_dim ',sx(ne),ex(ne),sy(ne),ey(ne)
!apm_d02-2011-10-14_00:00:00

!print*,'kk apm rst'

!print*,iyear,imonth,iday,ihour

allocate(ttbc(sx(ne):ex(ne),sy(ne):ey(ne),nzz))
allocate(ttoc(sx(ne):ex(ne),sy(ne):ey(ne),nzz))
allocate(ttdust(sx(ne):ex(ne),sy(ne):ey(ne),nzz))
allocate(ttsalt(sx(ne):ex(ne),sy(ne):ey(ne),nzz))

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! if restart apm, do not distribute sulfate in crude method
   do k=1,nzz
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         linit_box(iapm+ixy)=.false.
      enddo
      enddo
    enddo
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




write(rst_time(1:4),'(i4.4)') iyear
rst_time(5:5)= mbar
write(rst_time(6:7),'(i2.2)') imonth
rst_time(8:8)=mbar
write(rst_time(9:10),'(i2.2)') iday
rst_time(11:11)=bbar
write(rst_time(12:13),'(i2.2)') ihour
rst_time(14:14)=colon
write(rst_time(15:16),'(i2.2)') iminute
rst_time(17:17) = colon
write(rst_time(18:19),'(i2.2)') isecond

cdd(1:1)='d'
write(cdd(2:3),'(i2.2)') ne

if(myid.eq.0) print*,'retsart_time : ',rst_time

!stop
!print*,cdd
!print*,trim(naqpms_dir)
!print*,'k1','/rst_apm.init/'
!print*,'k2','apm_'//cdd//'-'//rst_time

rstfile = trim(naqpms_dir)//'/rst_apm.init/'//'apm_'//cdd//'.'//&
          &rst_time(1:11)//rst_time(12:13)//rst_time(15:16)//rst_time(18:19)//&
          &'.dat'

!print*,'restart file :',trim(rstfile)

!stop



call get_funit(funit)
open( funit,file=trim(rstfile),form='unformatted',access='direct' &
     ,recl=nx(ne)*ny(ne) )
! read tracers only

    id=ne

    irec=0
    ! write u_ws
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((u_ws(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    ! write v_ws
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((v_ws(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    ! write rh
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((rh(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    ! write clw
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((clw(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    ! write rnw
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((rnw(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    ! write raincon
    irec=irec+1
    !read(funit,rec=irec) ((raincon(i,j),i=1,nx(id)),j=1,ny(id))
    ! write rainnon
    irec=irec+1
    !read(funit,rec=irec) ((rainnon(i,j),i=1,nx(id)),j=1,ny(id))


    irec=irec+1
    !write(funit,rec=irec) ((vsb(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((wrfvsb(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((iffog(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((ifhaze(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((ifwrffog(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aod390(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aod500(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aod530(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aod550(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aod700(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aod1010(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aaod390(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aaod500(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aaod530(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aaod550(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aaod700(i,j),i=1,nx(id)),j=1,ny(id))

    irec=irec+1
    !write(funit,rec=irec) ((aaod1010(i,j),i=1,nx(id)),j=1,ny(id))

    do is=1,5 !ntyp
       irec=irec+1
       !write(funit,rec=irec) ((taod(i,j,is),i=1,nx(id)),j=1,ny(id))
    enddo

    do is=1,5 !ntyp
    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((ycs3d(i,j,k,is),i=1,nx(id)),j=1,ny(id))
    enddo
    enddo

    do is=1,5 !ntyp
    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((surf3d(i,j,k,is),i=1,nx(id)),j=1,ny(id))
    enddo
    enddo

    do is=1,5 !ntyp
    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((rgfdry3d(i,j,k,is),i=1,nx(id)),j=1,ny(id))
    enddo
    enddo



    ! write sulfate
    do is=1,NSO4
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( ( var2d(i,j),i=1,nx(id) ), j=1,ny(id) )
      iapm=ip_sulf(k,is,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         apm_sulf(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    enddo
    ! write sea salt
    do is=1,NSEA
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip_salt(k,is,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         apm_salt(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    enddo
    ! write dust
    do is=1,NDSTB
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip_dust(k,is,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         apm_dust(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    enddo
    ! write bcoc
    do is=1,NBCOCT
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip_bcoc(k,is,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         apm_bcoc(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    enddo

    ! write XQ (ions pairs : #/cm3)
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((apm_xq3d(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo

    ! write coated sulfate
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         msltsulf(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         mdstsulf(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         mbcsulf(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         mocsulf(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo

! write bulk sp sulfate mass
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) (( spsulf(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) (( ppsulf(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) (( ttsulf(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    ! write total sulfate mass
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((tot_sulf(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo

!    do k=1,nzz
!      irec=irec+1
!      !write(funit,rec=irec) ((ttdust(i,j,k),i=1,nx(id)),j=1,ny(id))
!    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         ttdust(i,j,k)=var2d(i,j)
      enddo
      enddo
    enddo



!    do k=1,nzz
!      irec=irec+1
!      !write(funit,rec=irec) ((ttsalt(i,j,k),i=1,nx(id)),j=1,ny(id))
!    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         ttsalt(i,j,k)=var2d(i,j)
      enddo
      enddo
    enddo




    ! write ccn number
    do is=1,5
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((number_ccn1(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    enddo

    do is=1,5
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((number_ccn2(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    enddo

    do is=1,5
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((number_ccn3(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    enddo

    do is=1,5
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((ztn3d(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    enddo

    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((real(cn3nm(i,j,k)),i=1,nx(id)),j=1,ny(id))
    enddo

   ! write cn10nm
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((cn10nm(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    ! write npf3d
    do k=1,nzz
      irec=irec+1
      !read(funit,rec=irec) ((npf3d(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo

    ! write h2so4 gas
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         h2so4_gas(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo

    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((real(apm_cacid(i,j,k)),i=1,nx(id)),j=1,ny(id))
    enddo
    !
    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((real(apm_pacid(i,j,k)),i=1,nx(id)),j=1,ny(id))
    enddo



    ! write particle radius growth factor
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         rgf_sulf(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         rgf_salt(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         rgf_dust(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         rgf_bc(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         rgf_oc(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo

    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((so4aq(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((clwph(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo
    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((spgf(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo


    do is=1,nbincb
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip_cbbin(k,is,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         apm_binbc(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    enddo


    do is=1,nbincb
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ((var2d(i,j),i=1,nx(id)),j=1,ny(id))
      iapm=ip_cbbin(k,is,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         apm_binoc(iapm+ixy)=var2d(i,j)
      enddo
      enddo
    enddo
    enddo

!    do k=1,nzz
!      irec=irec+1
!      !write(funit,rec=irec) ((ttbc(i,j,k),i=1,nx(id)),j=1,ny(id))
!    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         ttbc(i,j,k)=var2d(i,j)
      enddo
      enddo
    enddo



!    do k=1,nzz
!      irec=irec+1
!      !write(funit,rec=irec) ((ttoc(i,j,k),i=1,nx(id)),j=1,ny(id))
!    enddo
    do k=1,nzz
      irec=irec+1
      read(funit,rec=irec) ( (var2d(i,j),i=1,nx(id)),j=1,ny(id) )
      iapm=ip3mem(k,ne)
      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
         ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         ttoc(i,j,k)=var2d(i,j)
      enddo
      enddo
    enddo



    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((bcagt(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo

    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((ocagt(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo


    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((bcsprt(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo

    do k=1,nzz
      irec=irec+1
      !write(funit,rec=irec) ((ocsprt(i,j,k),i=1,nx(id)),j=1,ny(id))
    enddo


if(myid.eq.0)    print*,'restart record number =',irec


   do k=1,nzz
   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
      ixy=(ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

!      ttccn2=0.0
!      do is=1,5 ! ntyp
!        iapm=ip_type(k,is,ne)
!        ttccn2=ttccn2+number_ccn2(iapm+ixy)
!      enddo
      i03=ip3mem(k,ne)
      so4tmp=msltsulf(i03+ixy)*kapa(1)+ttsalt(i,j,k)*kapa(2) &
            +mdstsulf(i03+ixy)*kapa(1)+ttdust(i,j,k)*kapa(3) &
            +mbcsulf(i03+ixy)*kapa(1)+ttbc(i,j,k)*kapa(4) &
            +mocsulf(i03+ixy)*kapa(1)+ttoc(i,j,k)*kapa(5)
      do is=1,NSO4 
        iapm=ip_sulf(k,is,ne)
        so4tmp=so4tmp+apm_sulf(iapm+ixy)*kapa(1)
      enddo

      if(so4tmp.gt.0) then
 
        do is=1,NSO4 ! size
          iccn=ip_fccn(k,is,ne)
          frac_ccn(iccn+ixy)=apm_sulf(iapm+ixy)*kapa(1)/so4tmp
        enddo

        is=NSO4+1 ! bc
        iccn=ip_fccn(k,is,ne)
        frac_ccn(iccn+ixy)=(mbcsulf(i03+ixy)*kapa(1)+ttbc(i,j,k)*kapa(4))/so4tmp

        is=NSO4+2 ! oc
        iccn=ip_fccn(k,is,ne)
        frac_ccn(iccn+ixy)=(mocsulf(i03+ixy)*kapa(1)+ttoc(i,j,k)*kapa(5))/so4tmp

        is=NSO4+3 ! dust
        iccn=ip_fccn(k,is,ne)
        frac_ccn(iccn+ixy)=(mdstsulf(i03+ixy)*kapa(1)+ttdust(i,j,k)*kapa(3))/so4tmp

        is=NSO4+1 ! salt
        iccn=ip_fccn(k,is,ne)
        frac_ccn(iccn+ixy)=(msltsulf(i03+ixy)*kapa(1)+ttsalt(i,j,k)*kapa(2))/so4tmp

 
      else
        do is=1,NSO4+4 ! size
          iccn=ip_fccn(k,is,ne)
          frac_ccn(iccn+ixy)=0.0
        enddo
      endif
   enddo
   enddo
   enddo

    close(funit)

!stop

deallocate(ttbc,ttoc,ttdust,ttsalt)

return

end subroutine apm_rst




