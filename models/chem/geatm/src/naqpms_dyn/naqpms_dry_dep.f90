
subroutine naqpms_dry_dep &
 & ( myid &
 &  ,dt &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,igasCBM,igas,iaer,isize,nseacom,ndustcom &
 &  ,MSIZDIS,MSIZDID &
 &  ,iyear,imonth,iday,ihour,iminute )

use naqpms_varlist
use met_fields, only : Plev,t,u,v,QVAPOR
use met_fields, only : T2,SWDOWN,PBL_HGT,tskwrf
use met_fields, only : rmol,cldfrc3d
use met_fields, only : clw,rnw
use met_fields, only : wrflai
use naqpms_gridinfo, only : dz,heiz,HGT1,LAND_USE,LATITCRS,LONGICRS

use smpsulf_var, only: idx4dep_smpsulf

implicit none


!integer,parameter :: itzon=0

integer :: myid

real :: dt

integer :: imonth2

integer :: iyear,imonth,iday,ihour,iminute

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igasCBM,igas,iaer,isize,nseacom,ndustcom
!integer :: idry

integer :: i,j,k,is,ia,iduc,ig

integer :: mem2d,mem3d,mem2dgas,mem3daer

integer :: ixy,i02,i03,i04,I03AER,i02Gas,i05c,i05,i02aer,i04aer
integer :: i03_z1

real :: dryvel=0.0 ! need to be modified

real,dimension(isize+1) :: MSIZDIS,MSIZDID


!real,dimension(mem2d)    :: DUSTDRY,DUSTDRYSO4,DUSTDRYNO3 &
!                           ,DUSTDRYFEII,DUSTDRYFEIII
!real,dimension(mem3daer) :: DryVeldust

integer :: iwb,ieb,jsb,jeb

real   ,allocatable,dimension(:,:)   :: lai2d

real   ,allocatable,dimension(:,:)   :: tsurf,xlat,xlon,SWDOWN1,pblhgt !,vdep
integer,allocatable,dimension(:,:)   :: land
real   ,allocatable,dimension(:,:)   :: vdep,gas_vdep
real   ,allocatable,dimension(:,:,:) :: ppp,ttn,uu,vv,height1,pwc,cwc,water
real   ,allocatable,dimension(:,:)   :: diam2d,rhop2d

real   ,allocatable,dimension(:,:,:) :: fcloud

real   ,allocatable,dimension(:,:)   :: tmp2d

!real   ,allocatable,dimension(:,:)   :: zmol


real :: convfac

real :: rgf_tmp
real :: diam

real :: growfac

!real,dimension(mem3d) :: wk,wks
!real,dimension(mem2d) :: area0,area,ratio_area

real*8,allocatable,dimension(:) :: numcc,tmp,radius
real*8 :: totalbc,totaloc

! 1:FF 2:BB

real    :: vdepm
integer :: imlo,jmlo


!+++ sulferic acid vapor ++++!
logical :: lrddrt,lrdsno     !
integer,allocatable,dimension(:,:) :: icddrt,icdsno     !
!++++++++++++++++++++++++++++!


! shun@20161227
integer :: ig_ddep


!IF(lapm) THEN ! apm flag
!print*,'dims in sub-ddep',sx(ne),ex(ne),sy(ne),ey(ne),nzz

iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1


allocate( tmp2d(sx(ne):ex(ne),sy(ne):ey(ne)) )

allocate( diam2d(iwb:ieb,jsb:jeb) )
allocate( rhop2d(iwb:ieb,jsb:jeb) )

allocate( vdep(iwb:ieb,jsb:jeb) )

allocate( gas_vdep(iwb:ieb,jsb:jeb) )


IF(idry.gt.1) THEN ! calculate deposition velocity online

allocate( lai2d(iwb:ieb,jsb:jeb) )

allocate( tsurf(iwb:ieb,jsb:jeb) )
allocate( xlat(iwb:ieb,jsb:jeb) )
allocate( xlon(iwb:ieb,jsb:jeb) )
allocate( SWDOWN1(iwb:ieb,jsb:jeb) )
allocate( pblhgt(iwb:ieb,jsb:jeb) )

allocate( land(iwb:ieb,jsb:jeb) )

allocate( ppp(iwb:ieb,jsb:jeb,nzz) )
allocate( ttn(iwb:ieb,jsb:jeb,nzz) )
allocate( uu(iwb:ieb,jsb:jeb,nzz) )
allocate( vv(iwb:ieb,jsb:jeb,nzz) )
allocate( height1(iwb:ieb,jsb:jeb,nzz) )
!allocate( QQ(iwb:ieb,jsb:jeb,nzz) )
allocate( water(iwb:ieb,jsb:jeb,nzz) )


allocate( pwc(iwb:ieb,jsb:jeb,nzz) )
allocate( cwc(iwb:ieb,jsb:jeb,nzz) )
allocate( fcloud(iwb:ieb,jsb:jeb,nzz))

allocate( icddrt(iwb:ieb,jsb:jeb) )
allocate( icdsno(iwb:ieb,jsb:jeb) )

!allocate( zmol(iwb:ieb,jsb:jeb) )


! get local variables


!print*,ddep_flag

i02=ip2mem(ne)

do j=sy(ne)-1,ey(ne)+1
do i=sx(ne)-1,ex(ne)+1
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   
!   if(rmol(i02+ixy).ne.0) then
!     zmol(i,j)=1.0/rmol(i02+ixy)
!   else
!!     print*,'zmol abnormal:',i,j
!     zmol(i,j)=1.0e20 
!   endif


   do k=1,nzz
      i03=ip3mem(k,ne)
      i02=ip2mem(ne)
      ppp(i,j,k) = Plev(i03+ixy)
      ttn(i,j,k) = t(i03+ixy)
!      QQ(i,j,k)  = QVAPOR(i03+ixy)
      uu(i,j,k)  = u(i03+ixy)
      vv(i,j,k)  = v(i03+ixy)

      convfac=44.9*(273./t(i03+ixy))*(Plev(i03+ixy)/100/1013.)*29.0

      cwc(i,j,k)=clw(i03+ixy)*convfac
      pwc(i,j,k)=rnw(i03+ixy)*convfac

      water(i,j,k)  = QVAPOR(i03+ixy)*1.E06*29./18.!ppm

      fcloud(i,j,k) = cldfrc3d(i03+ixy)


      if(i.eq.83.and.j.eq.63.and.k.eq.1) then
!         print*,hlayer(i03+ixy)*1000,HGT1(i02+ixy)
      endif

      if(k==1) then
         height1(i,j,k) = 2*(heiz(i03+ixy)-HGT1(i02+ixy))
      endif
      if(k/=1) then
         height1(i,j,k) = 2*(heiz(i03+ixy)-HGT1(i02+ixy))-height1(i,j,k-1)
      endif
   enddo    !k     

   if(trim(ddep_flag).eq.'W89') then
     CALL getland(LAND_USE(i02+ixy),land(i,j))
   elseif(trim(ddep_flag).eq.'Z03') then
     CALL getland_usgs2zhang03 (LAND_USE(i02+ixy),land(i,j))
   endif

   tsurf(i,j)= tskwrf(i02+ixy)
   xlat(i,j) = LATITCRS(i02+ixy)
   xlon(i,j) = LONGICRS(i02+ixy)
   SWDOWN1(i,j)=SWDOWN(i02+ixy)
   pblhgt(i,j) = PBL_HGT(i02+ixy)

   if(lrd_lai) then
   lai2d(i,j)=wrflai(i02+ixy)
   endif

   if(i.eq.83.and.j.eq.63) then
!     print*,'83&63 dz=',height1(83,63,1)
   endif

enddo
enddo



  lrddrt=.false.;icddrt=0  !if droughtness index
  lrdsno=.false.;icdsno=0  !if snow index

!!!! ***    GAS SPECIES ****
!  do ig=1,igasCBM
do ig=1,iedgas

!print*,'ig=',ig

   ig_ddep=ig
   if(lgaschemsmp) then 
     ig_ddep=idx4dep_smpsulf(ig) !ig
   endif

   CALL drydep_gas(myid,ddep_flag,fcloud, &
               imonth2,tsurf,xlat,xlon,pwc,cwc,height1,ppp,uu,vv,&
               SWDOWN1,land,water, &
               lrd_lai,lai2d, &
               ttn,pblhgt,lrddrt,icddrt,lrdsno,icdsno,vdep,&
               itzon,iyear,imonth,iday,ihour,iminute, &
               sx(ne),ex(ne),sy(ne),ey(ne),nzz,ig_ddep )

   i02Gas = ip2memGas(ig,ne)
   do j=sy(ne)-1,ey(ne)+1
   do  i=sx(ne)-1,ex(ne)+1
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      DryVelGas(i02Gas+ixy) = vdep(i,j)
   enddo
   enddo

  do j=sy(ne),ey(ne)
   do  i=sx(ne),ex(ne)
      tmp2d(i,j) = vdep(i,j)
   enddo
   enddo



!print*,'vdep=',tmp2d
!stop


  enddo !ig


  i03_z1=ip3mem(1,ne)
!  growfac=hgfac(i03_z1+ixy)
  growfac=1.0

  if(laerv2) then
    do ia=1,naersp
    do is=1,naerbin 

      if(ia.gt.1.and.is.gt.1) cycle ! skip zero aersol tracer

      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        i03_z1=ip3mem(1,ne)
        growfac=hgfac(i03_z1+ixy)
        diam2d(i,j)=1.0E-6*diamlgm(is)*growfac

if(i.eq.21.and.j.eq.2) then
  print*,'kk hgfac=',growfac,diamlgm(is),diam2d(i,j)
endif
      enddo
      enddo
!      ig=77

      CALL drydep_aer(myid,ddep_flag,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv,&
                   land,ttn,lrd_lai,lai2d, &
                   pblhgt,diam2d,vdep,sx(ne),ex(ne),sy(ne),ey(ne),nzz,ig,aer_rhop(ia))
      
      i04aer = ip4mem_aer(1,is,ia,ne)
      i02aer = ip2memaer(is,ia,ne)

      do j=sy(ne)-1,ey(ne)+1
      do  i=sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        dryvelaer(i02aer+ixy) = vdep(i,j)
      enddo
      enddo

  do j=sy(ne),ey(ne)
   do  i=sx(ne),ex(ne)
      tmp2d(i,j) = vdep(i,j)
   enddo
   enddo

!print*,'vdep=',tmp2d
!stop


      i03=ip3mem(1,ne)
      call dry_dep_gas_zhu( myid, aerom(i04aer),dryvelaer(i02aer),dz(i03)  &
                         ,sx(ne), ex(ne), sy(ne), ey(ne), dt )

    enddo
    enddo
  endif

if(laerv1) then
  do ig = igasCBM + 1, igas

   IF(ig==igasCBM+2) THEN ! PM10
      diam = 1.0E-6*sqrt(2.5*10.0)*growfac
   ELSE
      diam = 1.0E-6*sqrt(0.1*2.5)*growfac  ! THINK OTHER AEROSOLS AS PM25
   ENDIF

      do j=sy(ne),ey(ne)
      do i=sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        i03_z1=ip3mem(1,ne)
        growfac=hgfac(i03_z1+ixy)
        diam2d(i,j)=1.0E-6*diam*growfac
      enddo
      enddo

   CALL drydep_aer_old(myid,ddep_flag,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv,&
                   land,ttn,lrd_lai,lai2d, &
                   pblhgt,diam2d,vdep,sx(ne),ex(ne),sy(ne),ey(ne),nzz,ig)


   i02Gas = ip2memGas(ig,ne)

   do j=sy(ne)-1,ey(ne)+1
    do  i=sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     DryVelGas(i02Gas+ixy) = vdep(i,j)
    enddo
   enddo
!   print*,'vdeppp=',vdep
!   stop

  enddo !ig

endif

!stop

  DO IA=1,IAER
  DO IS=1,ISIZE
    IF(IA==1) THEN
      DIAM = 1.0E-6*SQRT(MSIZDIS(IS+1)*MSIZDIS(IS))
    ELSE IF(IA==2) THEN
      DIAM = 1.0E-6*SQRT(MSIZDID(IS+1)*MSIZDID(IS))
    ENDIF
    CALL drydep_dust_salt(myid,ddep_flag &
                ,imonth2,tsurf,xlat,xlon,height1,ppp,uu,vv,lrd_lai,lai2d &
                ,land,ttn,pblhgt,diam,vdep,sx(ne),ex(ne),sy(ne),ey(ne),nzz,IA )
    I03AER = IP3MEMAER(IS,IA,NE)
    do j=sy(ne)-1,ey(ne)+1
    do  i=sx(ne)-1,ex(ne)+1
     ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
     DryVeldust(I03AER+ixy) = vdep(i,j)*0.45 ! ???
    enddo
    enddo
  ENDDO ! IS
  ENDDO ! IA

deallocate( ppp,ttn,land,tsurf,xlat,xlon,pwc,cwc,height1,uu,vv &
           ,water,SWDOWN1,vdep,pblhgt,fcloud,icddrt,icdsno)

ENDIF ! whether to calculate deposition velocity


!=======================================================================

! do ig=1,igas !Gas
do ig=1,iedgas !Gas

    i04 = ip4mem(1,ig,ne)
    i02Gas = ip2memGas(ig,ne)
    i03=ip3mem(1,ne)
    call dry_dep_gas_zhu( myid, gas(i04),DryVelGas(i02gas),dz(i03)  &
                         ,sx(ne), ex(ne), sy(ne), ey(ne), dt )
 enddo !ig


 DO IA=1,IAER
 DO IS=1,ISIZE
   I03AER = IP3MEMAER(IS,IA,NE)
   I05=IP5MEM(1,IS,IA,NE)
   i03=ip3mem(1,ne)
   I02 = IP2MEM(NE)


   IF(IA==2) THEN
      CALL DUSTDRYDEP( MYID, DUSTDRY(I02), AER(I05),DryVeldust(I03AER) & 
                      ,dz(i03),sx(ne), ex(ne),sy(ne), ey(ne), dt ) 
                      ! to calculate the DUST DRY DEPOSITION  UG/M2/HOUR
      do iduc = 1, ndustcom
        i05c = ip5memc (1,is,iduc,ne)

        IF(iduc==7) then
          CALL DUSTDRYDEP( MYID, DUSTDRYSO4(I02),DUSTCOMP(I05C) &
                          ,DryVeldust(I03AER) &
                          ,dz(i03),sx(ne),ex(ne),sy(ne),ey(ne), dt)
                 !  to calculate the DUST DSO4 DRY DEPOSITION  UG/M2/HOUR
         ELSE IF(iduc == 8) then
          CALL DUSTDRYDEP( MYID, DUSTDRYNO3(I02),DUSTCOMP(I05C) &
                          ,DryVeldust(I03AER) &
                          ,dz(i03), sx(ne),ex(ne), sy(ne),ey(ne), dt) 
                 !  to calculate the DUST DNO3 DRY DEPOSITION  UG/M2/HOUR
         ELSE IF (iduc==9) then
          CALL DUSTDRYDEP( MYID, DUSTDRYFeII(I02),DUSTCOMP(I05C) &
                          ,DryVeldust(I03AER) &
                          ,dz(i03),sx(ne),ex(ne), sy(ne),ey(ne), dt) 
                 !  to calculate the DUST DFeII DRY DEPOSITION  UG/M2/HOUR
         ELSE IF(iduc==10) then
          CALL DUSTDRYDEP( MYID, DUSTDRYFeIII(I02),DUSTCOMP(I05C) &
                          ,DryVeldust(I03AER) &
                          ,dz(i03),sx(ne),ex(ne),sy(ne),ey(ne),dt ) 
                 !  to calculate the DUST DFeIII DRY DEPOSITION  UG/M2/HOUR
         ENDIF
        enddo

      ENDIF

      call dry_dep_gas_zhu( myid, AER(i05),DryVeldust(I03AER),dz(i03)  &
                           ,sx(ne), ex(ne), sy(ne), ey(ne), dt)

      IF(IA==1)  THEN !SEA SALT 
       do iduc = 1,nseacom
        i05c = ip5memcs (1,is,iduc,ne)
        CALL  dry_dep_gas_zhu( myid, SEACOMP(i05c),DryVeldust(I03AER) &
                              ,dz(i03), sx(ne), ex(ne), sy(ne), ey(ne), dt)
       enddo  ! iduc
      ELSE IF(IA==2) THEN ! DUST
       do iduc = 1,ndustcom
        i05c = ip5memc (1,is,iduc,ne)
        CALL  dry_dep_gas_zhu( myid, DUSTCOMP(i05c),DryVeldust(I03AER) &
                              ,dz(i03), sx(ne), ex(ne), sy(ne), ey(ne), dt)
       enddo  ! iduc       
      ENDIF

    ENDDO ! IS
    ENDDO ! IA


end subroutine naqpms_dry_dep



