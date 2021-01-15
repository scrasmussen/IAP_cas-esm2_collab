
subroutine naqpms_cld_convect &
 & ( myid &
 &  ,ne &
 &  ,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d,ktop &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,ifsm,ifsmt,idmSet,ismMax,igMark )

use naqpms_varlist
use met_fields
!use naqpms_gridinfo
implicit none

real    :: dtt,dt0,ttime
integer :: nstep

integer :: myid

integer :: IFLAG

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: ifsm(5)

integer :: ifsmt

integer :: idmSet,ismMax

integer :: igMark(idmSet)


integer :: igas,iaer,isize,nseacom,ndustcom

integer :: iconv
integer :: i,j,k,is

integer :: mem2d,mem3d

real :: ktop(mem2d)

integer :: nn
real :: dttnn
real :: deltc1,deltc2,deltc3,deltc4

real,dimension(mem2d) :: CBMF1

real    :: FCBMF

integer :: ixy,i02,i03,i03_1,i0,iapm


integer :: ig,i04,ia,iduc,i05,i05c

integer :: ism,idm,i04sm

integer :: letdoit

real,allocatable,dimension(:,:)  :: smconv1

real,allocatable,dimension(:) :: c00,smthis,smother


real :: T1(nzz),Q1(nzz),QS1(nzz),U1(nzz),V1(nzz),TRA1(nzz,1),P1(nzz) &
     & ,PH1(nzz+1),FTRA1(nzz,1),FTRA1D(nzz,1),FTRA1U(nzz,1),FTRA1O(nzz,1) &
     & ,FTRA1E(nzz,1)


   i0=ip2mem(ne)

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
!do j=30,30
!do i=30,30

      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
   do ig=1,iedgas
!   do ig=1,1
     dt0=600.   !600 seconds,cannot alter unless you have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=1.0E-20
     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)= AMAX1(gas(i04+ixy) , 1.E-20)
      FTRA1(k,1)=1.0E-20
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
     enddo !k=nzz
     i03_1=ip3mem(nzz,ne)
     PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)

!stop 'onegas'

     DO iconv=1,1000 !1000 no meaning ,just for safty
990  CONTINUE
     dtt=dt0/nstep
     ttime=ttime+dtt
     FCBMF=CBMF1(i0+ixy)
     IFLAG=0

!stop 'kk1'
if(.false.) then 
print*,'T1=',T1
print*,'Q1=',Q1
print*,'QS1=',QS1
print*,'U1=',U1
print*,'V1=',V1
print*,'TRA1=',TRA1
print*,'P1=',P1
print*,'PH1=',PH1
print*,'nzz=',nzz,nzz-1,1,dtt
print*,'FTRA1=',FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E
endif
!stop
     CALL  CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG)

!stop 'kk2'

     IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 992 ! No moist convection
     IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2
      ttime=ttime-dtt
      GOTO 990
     ENDIF
     CBMF1(i0+ixy)=FCBMF


     IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

 !!!!!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  ! to close vertical source mark
    if(ifsm(ne)==1)then
      allocate(smconv1(ismMax,nzz), c00(nzz), smthis(ismMax),smother(ismMax))

      letdoit=0
      do idm=1,idmSet
        if(igMark(idm)==ig)letdoit=idm
      enddo
      if(letdoit>0)then

        do ism=1,ismMax
        do k=1,nzz
         i04sm=ipSMmem(k,ism,letdoit,ne)
         i04=ip4mem(k,ig,ne)
         smconv1 (ism,k)= SourceMark(i04sm +ixy)
         c00(k)=gas(i04+ixy)
         c00(k)=amax1(c00(k), 1.E-20)
        enddo !k
        enddo !ism
        DO nn=1,40 !!no meaning, just reduce the errors due to the order of up and down
         dttnn=dtt/40.
         i02=ip2mem(ne)
         do k=1,ktop(i02+ixy)
           deltc1=dttnn*FTRA1D(k,1)
           deltc2=dttnn*FTRA1U(k,1)
           deltc3=dttnn*FTRA1O(k,1)
           deltc4=dttnn*FTRA1E(k,1)
    !-----------downdraft
           IF (deltc1>0.0) THEN
            do is=1,ismMax
             smthis(is)=smconv1(is,k)
             if(k==nzz)  smother(is)=smconv1(is,k)
             if(k<nzz)   smother(is)=smconv1(is,k+1)
            enddo
            if(k==ktop(i02+ixy)) then
             call GetBounChange(deltc1,c00(k),smthis,2,ismMax)
            else if(k>ktop(i02+ixy)) then
             smthis(2)=1.
            else if(k<ktop(i02+ixy))then
             call SMmixing(c00(k),smthis,deltc1,smother,ismMax)
            endif !k

            do is=1,ismMax
             smconv1(is,k)=smthis(is)
            enddo
           ENDIF
           c00(k)=c00(k)+deltc1
           c00(k)=amax1(c00(k),1.E-20)
    !---------updraft
           IF(deltc2>0.0)THEN
            do is=1,ismMax
             smthis(is)=smconv1(is,k)
             if(k==1) smother(is)=smconv1(is,1)
             if(k>1) smother(is)=smconv1(is,k-1)
            enddo

            if(k>ktop(i02+ixy)) then
             smthis(2)=1.
            else if(k<=ktop(i02+ixy)) then
             call SMmixing(c00(k),smthis,deltc2,smother,ismMax)
            endif!k

            do is=1,ismMax
              smconv1(is,k)=smthis(is)
            enddo
           ENDIF
           c00(k)=c00(k)+deltc2
           c00(k)=c00(k)+deltc3
           c00(k)=c00(k)+deltc4
           c00(k)=amax1(c00(k),1.E-20)
         enddo !k=ktop
        ENDDO !nn

        do ism=1,ismMax
        do k=1,nzz
         i04sm=ipSMmem(k,ism,letdoit,ne)
         SourceMark(i04sm +ixy)=smconv1(ism,k)
        enddo !k
        enddo !ism

      endif !letdoit
      deallocate(smconv1,c00,smthis,smother)

     endif !ifsm
 !!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!


     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      gas(i04+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
     enddo!k

    ENDIF !ttime

     IF(nint(ttime)>dt0*7.) GOTO 991  ! 1 hour lifetime

     ENDDO !iconv

 991  CONTINUE
   enddo !igas   
 992 CONTINUE

   enddo ! j
   enddo ! i

  !deallocate(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E)



!stop 'naqpms_cld_gas'

! ---- gas and sec aerosols

! ---- FOR SEA SALT AND DUST EMISSIONS ----
 !allocate(T1(nzz),Q1(nzz),QS1(nzz),U1(nzz),V1(nzz),TRA1(nzz,1),P1(nzz),&
 !          PH1(nzz+1),FTRA1(nzz,1),FTRA1D(nzz,1),FTRA1U(nzz,1),FTRA1O(nzz,1),FTRA1E(nzz,1))
           !FTRA1 is net, FTRA1D is downdraft, FTRA1U is updraft, 
           !FTRA1O is the export form the cell, FTRA1E is enchange of
           !environment 

  ! in 1-3-5-7-9-11 itt if block 

    i0=ip2mem(ne)

   do j=sy(ne),ey(ne)
   do i=sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1


  DO IA=1,IAER
  DO IS=1,ISIZE

     dt0=600.   !600 seconds,cannot alter unless yo have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=0.0
   do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      I05=IP5MEM(K,IS,IA,NE)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)=AMAX1(AER(i05+ixy), 1.E-20)
      FTRA1(k,1)=0.0
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
    enddo !k=nzz

      i03_1=ip3mem(nzz,ne)
      PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)

   DO iconv=1,1000 !1000 no meaning ,just for safty
993  CONTINUE
    dtt=dt0/nstep
    ttime=ttime+dtt
    FCBMF=CBMF1(i0+ixy)
    IFLAG=0
    CALL CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG)
    IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 995
    IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2
      ttime=ttime-dtt
      GOTO 993
    ENDIF
        CBMF1(i0+ixy)=FCBMF

    IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      I05=IP5MEM(K,IS,IA,NE)
      AER(i05+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
     enddo!k
    ENDIF !ttime

     IF(nint(ttime)>dt0*7.) GOTO 994

   ENDDO !iconv
 994  CONTINUE


   IF(IA==1) THEN ! FOR SEA SALT
     DO iduc = 1, nseacom
     dt0=600.   !600 seconds,cannot alter unless yo have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=0.0
   do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      i05c = ip5memcs(k,is,iduc,ne)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)=AMAX1(SEACOMP(i05c+ixy), 1.E-20)
      FTRA1(k,1)=0.0
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
    enddo !k=nzz
      i03_1=ip3mem(nzz,ne)
      PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)
   DO iconv=1,1000 !1000 no meaning ,just for safty
996  CONTINUE
    dtt=dt0/nstep
    ttime=ttime+dtt
    FCBMF=CBMF1(i0+ixy)
    IFLAG=0
    CALL CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG)
    IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 995
    IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2
      ttime=ttime-dtt
      GOTO 996
    ENDIF
        CBMF1(i0+ixy)=FCBMF


   IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i05c = ip5memcs(k,is,iduc,ne)
      SEACOMP(i05c+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
     enddo!k
    ENDIF !ttime

     IF(nint(ttime)>dt0*7.) GOTO 997

   ENDDO !iconv
 997  CONTINUE

     ENDDO ! iduc

    ELSE IF (IA==2) THEN ! FOR DUST


     DO iduc = 1, ndustcom
     dt0=600.   !600 seconds,cannot alter unless yo have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=0.0
   do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i03=ip3mem(k,ne)
      i05c = ip5memc (k,is,iduc,ne)
      if(k/=1) i03_1=ip3mem(k-1,ne)
      U1(k)=u(i03+ixy)
      V1(k)=v(i03+ixy)
      T1(k)=t(i03+ixy)
      Q1(k)=QVAPOR(i03+ixy)/(QVAPOR(i03+ixy)+1.) !specific huminity
      QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
      P1(k)=Plev(i03+ixy)
      TRA1(k,1)=AMAX1(DUSTCOMP(i05c+ixy), 1.E-20)
      FTRA1(k,1)=0.0
      if(k==1) PH1(k)=PSFC(i0+ixy)/100.
      if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
    enddo !k=nzz
      i03_1=ip3mem(nzz,ne)
      PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)
   DO iconv=1,1000 !1000 no meaning ,just for safty
892  CONTINUE
    dtt=dt0/nstep
    ttime=ttime+dtt
    FCBMF=CBMF1(i0+ixy)
    IFLAG=0
    CALL CONVECT(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG)
    IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 995
    IF(IFLAG==4.and.nstep<=16) THEN
      nstep=nstep*2
      ttime=ttime-dtt
      GOTO 892
    ENDIF
        CBMF1(i0+ixy)=FCBMF


    IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN !after initial (2hr)time

     do k=1,nzz
      i04=ip4mem(k,ig,ne)
      i05c = ip5memc (k,is,iduc,ne)
      DUSTCOMP(i05c+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
     enddo!k
    ENDIF !ttime

     IF(nint(ttime)>dt0*7.) GOTO 998

   ENDDO !iconv
 998  CONTINUE
     ENDDO ! iduc

    ENDIF ! IAIF

    ENDDO ! IS
    ENDDO ! IA  

 995 CONTINUE

   enddo !i
  enddo !j
    !deallocate(T1,Q1,QS1,U1,V1,TRA1,P1,PH1,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E)

end subroutine naqpms_cld_convect





