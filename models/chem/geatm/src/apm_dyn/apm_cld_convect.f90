
subroutine apm_cld_convect &
 & ( myid &
 &  ,lapm &
 &  ,ne &
 &  ,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,u,v,t,QVAPOR,rh1,Plev,PSFC &
 &  ,ip2mem,mem2d &
 &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'

real    :: dtt,dt0,ttime
integer :: nstep

integer :: myid

logical :: lapm
integer :: IFLAG

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: iconv
integer :: i,j,k,is

integer :: mem2d,mem3d

real,dimension(mem3d) :: u,v,t,QVAPOR,rh1,Plev

real,dimension(mem2d) :: PSFC,CBMF1

real    :: FCBMF

integer :: ixy,i03,i03_1,i0,iapm

integer :: ip2mem(nest)
integer :: ip3mem(nzz,nest)

real :: T1(nzz),Q1(nzz),QS1(nzz),U1(nzz),V1(nzz),TRA1(nzz,1),P1(nzz) &
     & ,PH1(nzz+1),FTRA1(nzz,1),FTRA1D(nzz,1),FTRA1U(nzz,1),FTRA1O(nzz,1) &
     & ,FTRA1E(nzz,1)



!IF(lapm) THEN ! apm flag

 i0=ip2mem(ne)

 !loop_j : do j=sy(ne)-1,ey(ne)+1
 !loop_i : do i=sx(ne)-1,ex(ne)+1

  loop_j : do j=sy(ne),ey(ne)
  loop_i : do i=sx(ne),ex(ne)


   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

   !> meteorological inputs
   do k=1,nzz
     i03=ip3mem(k,ne)
     U1(k) = u (i03+ixy)
     V1(k) = v (i03+ixy)
     T1(k) = t (i03+ixy)
     Q1(k) = QVAPOR(i03+ixy)/(1.0+QVAPOR(i03+ixy))
     QS1(k)=Q1(k)/rh1(i03+ixy)*100.  !sturation specific huminity
     P1(k)=Plev(i03+ixy)
     if(k/=1) i03_1=ip3mem(k-1,ne)
     if(k==1) PH1(k)=PSFC(i0+ixy)/100.
     if(k>1)  PH1(k)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(k-1)
   enddo
   i03_1=ip3mem(nzz,ne)
   PH1(nzz+1)=Plev(i03_1+ixy)+Plev(i03_1+ixy)-PH1(nzz)
   !


   !> sulfate
  if(lfor_sulf) then
   !if(i.eq.1.and.j.eq.1) print*,'sulf cloud convect'
   loop_so4 : do is=1,NSO4

      dt0=600.   !600 seconds,cannot alter unless you have enough reasons
      nstep=1
      ttime=0
      CBMF1(i0+ixy)=1.0E-20
      do k=1,nzz
        iapm=ip_sulf(k,is,ne)
        TRA1(k,1) = AMAX1(apm_sulf(iapm+ixy) , 1.E-20)
        FTRA1(k,1)=1.0E-20
      enddo

!print*,'kktr',TRA1(:,1)

!stop

      !
      DO iconv=1,1000 !1000 no meaning ,just for safty
990     continue
        dtt=dt0/nstep
        ttime=ttime+dtt
        FCBMF=CBMF1(i0+ixy)
        IFLAG=0

!pause

!print*,'T1',T1
!print*,'Q1',Q1
!PRINT*,'QS1',QS1
!PRINT*,'U1',U1
!PRINT*,'V1',V1
!PRINT*,'TRA1',TRA1
!PRINT*,'P1',P1
!PRINT*,'PH1',PH1

!PAUSE

        CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                    & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )

!stop

        ! No moist convection occured at this grid point, goto 992
        ! to conduct calculation at next grid point
        IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1001
        IF(IFLAG==4.and.nstep<=16) THEN
          nstep=nstep*2
          ttime=ttime-dtt
          GOTO 990
        ENDIF
        CBMF1(i0+ixy)=FCBMF
        IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
          do k=1,nzz
            iapm=ip_sulf(k,is,ne)
            apm_sulf(iapm+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
            apm_sulf(iapm+ixy)=amax1(apm_sulf(iapm+ixy),0.0)
          enddo!k
        ENDIF !ttime
        IF(nint(ttime)>dt0*7.) GOTO 991  ! 1 hour lifetime
      ENDDO !iconv
991   continue
   enddo loop_so4
  endif
   !< end of sulfate

1001 continue

!stop
!cycle

   !> seasalt
  if(lfor_salt) then
   !if(i.eq.1.and.j.eq.1) print*,'salt cloud convect'
   loop_salt : do is=1,NSEA
      dt0=600.   !600 seconds,cannot alter unless you have enough reasons
      nstep=1
      ttime=0
      CBMF1(i0+ixy)=1.0E-20
      do k=1,nzz
        iapm=ip_salt(k,is,ne)
        TRA1(k,1)= AMAX1(apm_salt(iapm+ixy) , 1.E-20)
        FTRA1(k,1)=1.0E-20
      enddo
      !
      DO iconv=1,1000 !1000 no meaning ,just for safty
890     continue
        dtt=dt0/nstep
        ttime=ttime+dtt
        FCBMF=CBMF1(i0+ixy)
        IFLAG=0
        CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                    & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
        ! No moist convection occured at this grid point, goto 992
        ! to conduct calculation at next grid point
        IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1002
        IF(IFLAG==4.and.nstep<=16) THEN
          nstep=nstep*2
          ttime=ttime-dtt
          GOTO 890
        ENDIF
        CBMF1(i0+ixy)=FCBMF
        IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
          do k=1,nzz
            iapm=ip_salt(k,is,ne)
            apm_salt(iapm+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
            apm_salt(iapm+ixy)=amax1(apm_salt(iapm+ixy),0.0)
          enddo!k
        ENDIF !ttime
        IF(nint(ttime)>dt0*7.) GOTO 891  ! 1 hour lifetime
      ENDDO !iconv
891   continue
   enddo loop_salt
  endif
   !< end of seasalt

1002 continue


   !> dust
  if(lfor_dust) then
   !if(i.eq.1.and.j.eq.1) print*,'dust cloud convect'
   loop_dust : do is=1,NDSTB
      dt0=600.   !600 seconds,cannot alter unless you have enough reasons
      nstep=1
      ttime=0
      CBMF1(i0+ixy)=1.0E-20
      do k=1,nzz
        iapm=ip_dust(k,is,ne)
        TRA1(k,1)= AMAX1(apm_dust(iapm+ixy) , 1.E-20)
        FTRA1(k,1)=1.0E-20
      enddo
      !
      DO iconv=1,1000 !1000 no meaning ,just for safty
790     continue
        dtt=dt0/nstep
        ttime=ttime+dtt
        FCBMF=CBMF1(i0+ixy)
        IFLAG=0
        CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                    & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
        ! No moist convection occured at this grid point, goto 992
        ! to conduct calculation at next grid point
        IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1003
        IF(IFLAG==4.and.nstep<=16) THEN
          nstep=nstep*2
          ttime=ttime-dtt
          GOTO 790
        ENDIF
        CBMF1(i0+ixy)=FCBMF
        IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
          do k=1,nzz
            iapm=ip_dust(k,is,ne)
            apm_dust(iapm+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
            apm_dust(iapm+ixy)=amax1(apm_dust(iapm+ixy),0.0)
          enddo!k
        ENDIF !ttime
        IF(nint(ttime)>dt0*7.) GOTO 791  ! 1 hour lifetime
      ENDDO !iconv
791   continue
   enddo loop_dust
  endif
   !< end of dust

1003 continue

   !> bcoc
  if(lfor_bcoc) then
   !if(i.eq.1.and.j.eq.1) print*,'bcoc cloud convect'
   loop_bcoc : do is=1,NBCOCT
      dt0=600.   !600 seconds,cannot alter unless you have enough reasons
      nstep=1
      ttime=0
      CBMF1(i0+ixy)=1.0E-20
      do k=1,nzz
        iapm=ip_bcoc(k,is,ne)
        TRA1(k,1)= AMAX1(apm_bcoc(iapm+ixy) , 1.E-20)
        FTRA1(k,1)=1.0E-20
      enddo
      !
      DO iconv=1,1000 !1000 no meaning ,just for safty
690     continue
        dtt=dt0/nstep
        ttime=ttime+dtt
        FCBMF=CBMF1(i0+ixy)
        IFLAG=0
        CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                    & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
        ! No moist convection occured at this grid point, goto 992
        ! to conduct calculation at next grid point
        IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1004
        IF(IFLAG==4.and.nstep<=16) THEN
          nstep=nstep*2
          ttime=ttime-dtt
          GOTO 690
        ENDIF
        CBMF1(i0+ixy)=FCBMF
        IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
          do k=1,nzz
            iapm=ip_bcoc(k,is,ne)
            apm_bcoc(iapm+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
            apm_bcoc(iapm+ixy)=amax1(apm_bcoc(iapm+ixy),0.0)
          enddo!k
        ENDIF !ttime
        IF(nint(ttime)>dt0*7.) GOTO 691  ! 1 hour lifetime
      ENDDO !iconv
691   continue
   enddo loop_bcoc
  endif
   !< end of bcoc

1004 continue

!========================================================================
!========================================================================

! apm coated species

 IF(lcoated_dyn) THEN
     !if(i.eq.1.and.j.eq.1) print*,'coating cloud convect'
!-> sulfate on seasalt
     dt0=600.   !600 seconds,cannot alter unless you have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=1.0E-20
     do k=1,nzz
       i03=ip3mem(k,ne)
       TRA1(k,1)= AMAX1(msltsulf(i03+ixy) , 1.E-20)
       FTRA1(k,1)=1.0E-20
     enddo
     !
     DO iconv=1,1000 !1000 no meaning ,just for safty
590    continue
       dtt=dt0/nstep
       ttime=ttime+dtt
       FCBMF=CBMF1(i0+ixy)
       IFLAG=0
       CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                   & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
       ! No moist convection occured at this grid point, goto 992
       ! to conduct calculation at next grid point
       IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1005
       IF(IFLAG==4.and.nstep<=16) THEN
         nstep=nstep*2
         ttime=ttime-dtt
         GOTO 590
       ENDIF
       CBMF1(i0+ixy)=FCBMF
       IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
         do k=1,nzz
           i03=ip3mem(k,ne)
           msltsulf(i03+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
           msltsulf(i03+ixy)=amax1(msltsulf(i03+ixy),0.0)
         enddo!k
       ENDIF !ttime
       IF(nint(ttime)>dt0*7.) GOTO 591  ! 1 hour lifetime
     ENDDO !iconv
591  continue
1005 continue
!-> sulfate on dust
     dt0=600.   !600 seconds,cannot alter unless you have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=1.0E-20
     do k=1,nzz
       i03=ip3mem(k,ne)
       TRA1(k,1)= AMAX1(mdstsulf(i03+ixy) , 1.E-20)
       FTRA1(k,1)=1.0E-20
     enddo
     !
     DO iconv=1,1000 !1000 no meaning ,just for safty
490    continue
       dtt=dt0/nstep
       ttime=ttime+dtt
       FCBMF=CBMF1(i0+ixy)
       IFLAG=0
       CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                   & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
       ! No moist convection occured at this grid point, goto 992
       ! to conduct calculation at next grid point
       IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1006
       IF(IFLAG==4.and.nstep<=16) THEN
         nstep=nstep*2
         ttime=ttime-dtt
         GOTO 490
       ENDIF
       CBMF1(i0+ixy)=FCBMF
       IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
         do k=1,nzz
           i03=ip3mem(k,ne)
           mdstsulf(i03+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
           mdstsulf(i03+ixy)=amax1(mdstsulf(i03+ixy),0.0)
         enddo!k
       ENDIF !ttime
       IF(nint(ttime)>dt0*7.) GOTO 491  ! 1 hour lifetime
     ENDDO !iconv
491  continue
1006 continue
!-> sulfate on BC
     dt0=600.   !600 seconds,cannot alter unless you have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=1.0E-20
     do k=1,nzz
       i03=ip3mem(k,ne)
       TRA1(k,1)= AMAX1(mbcsulf(i03+ixy) , 1.E-20)
       FTRA1(k,1)=1.0E-20
     enddo
     !
     DO iconv=1,1000 !1000 no meaning ,just for safty
390    continue
       dtt=dt0/nstep
       ttime=ttime+dtt
       FCBMF=CBMF1(i0+ixy)
       IFLAG=0
       CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                   & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
       ! No moist convection occured at this grid point, goto 992
       ! to conduct calculation at next grid point
       IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1007
       IF(IFLAG==4.and.nstep<=16) THEN
         nstep=nstep*2
         ttime=ttime-dtt
         GOTO 390
       ENDIF
       CBMF1(i0+ixy)=FCBMF
       IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
         do k=1,nzz
           i03=ip3mem(k,ne)
           mbcsulf(i03+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
           mbcsulf(i03+ixy)=amax1(mbcsulf(i03+ixy),0.0)
         enddo!k
       ENDIF !ttime
       IF(nint(ttime)>dt0*7.) GOTO 391  ! 1 hour lifetime
     ENDDO !iconv
391  continue
1007 continue

!-> sulfate on POC
     dt0=600.   !600 seconds,cannot alter unless you have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=1.0E-20
     do k=1,nzz
       i03=ip3mem(k,ne)
       TRA1(k,1)= AMAX1(mocsulf(i03+ixy) , 1.E-20)
       FTRA1(k,1)=1.0E-20
     enddo
     !
     DO iconv=1,1000 !1000 no meaning ,just for safty
290    continue
       dtt=dt0/nstep
       ttime=ttime+dtt
       FCBMF=CBMF1(i0+ixy)
       IFLAG=0
       CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                   & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
       ! No moist convection occured at this grid point, goto 992
       ! to conduct calculation at next grid point
       IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1008
       IF(IFLAG==4.and.nstep<=16) THEN
         nstep=nstep*2
         ttime=ttime-dtt
         GOTO 290
       ENDIF
       CBMF1(i0+ixy)=FCBMF
       IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
         do k=1,nzz
           i03=ip3mem(k,ne)
           mocsulf(i03+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
           mocsulf(i03+ixy)=amax1(mocsulf(i03+ixy),0.0)
         enddo!k
       ENDIF !ttime
       IF(nint(ttime)>dt0*7.) GOTO 291  ! 1 hour lifetime
     ENDDO !iconv
291  continue
1008 continue 
  ENDIF ! lcoated_dyn

!-> sulferic acid vapor
     dt0=600.   !600 seconds,cannot alter unless you have enough reasons
     nstep=1
     ttime=0
     CBMF1(i0+ixy)=1.0E-20
     do k=1,nzz
       i03=ip3mem(k,ne)
       TRA1(k,1)= AMAX1(h2so4_gas(i03+ixy) , 1.E-20)
       FTRA1(k,1)=1.0E-20
     enddo
     !
     DO iconv=1,1000 !1000 no meaning ,just for safty
190    continue
       dtt=dt0/nstep
       ttime=ttime+dtt
       FCBMF=CBMF1(i0+ixy)
       IFLAG=0
       CALL  CONVECT( T1,Q1,QS1,U1,V1,TRA1,P1,PH1,nzz,nzz-1,1,dtt &
                   & ,FTRA1,FTRA1D,FTRA1U,FTRA1O,FTRA1E,FCBMF,IFLAG )
       ! No moist convection occured at this grid point, goto 992
       ! to conduct calculation at next grid point
       IF (IFLAG .ne. 1 .and.IFLAG .ne. 4) goto 1009
       IF(IFLAG==4.and.nstep<=16) THEN
         nstep=nstep*2
         ttime=ttime-dtt
         GOTO 190
       ENDIF
       CBMF1(i0+ixy)=FCBMF
       IF(ttime<=7.*dt0.and.ttime>=dt0*6.) THEN ! conv cld lifetime : 1hr
         do k=1,nzz
           i03=ip3mem(k,ne)
           h2so4_gas(i03+ixy)=TRA1(k,1)+FTRA1(k,1)*dtt
           h2so4_gas(i03+ixy)=amax1(h2so4_gas(i03+ixy),0.0)
         enddo!k
       ENDIF !ttime
       IF(nint(ttime)>dt0*7.) GOTO 191  ! 1 hour lifetime
     ENDDO !iconv
191 continue
1009 continue


992 CONTINUE
 enddo loop_i
 enddo loop_j


!ENDIF ! apm flag


end subroutine apm_cld_convect





