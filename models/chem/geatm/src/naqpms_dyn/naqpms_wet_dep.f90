

subroutine naqpms_wet_dep &
 & ( myid &
 &  ,itt_in &
 &  ,dt &
 &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,CLW,RNW,temp,Plev &
 &  ,RAINCON,RAINNON &
 &  ,mem2d &
 &  ,mem3d &
 &  ,igasCBM,igas,iaer,isize,nseacom,ndustcom &
 &  ,mem2dgas ) !&
! &  ,DUSTWET,DUSTWETSO4,DUSTWETNO3,DUSTWETFEII,DUSTWETFEIII )

use wdepflags, only: laerincld
use naqpms_varlist
!use met_fields
use naqpms_gridinfo, only : dx,dy,dz 

use smpsulf_var, only: idx4dep_smpsulf

implicit none

integer :: myid

integer :: itt_in,itt
real    :: dt


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igasCBM,igas,iaer,isize,nseacom,ndustcom

integer :: i,j,k,is,ia,ig,iduc

integer :: mem2d,mem3d,mem2dgas

real,dimension(mem2d) :: RAINCON,RAINNON

!real,dimension(mem2d) :: DUSTWET,DUSTWETSO4,DUSTWETNO3 &
!                        ,DUSTWETFEII,DUSTWETFEIII

real,dimension(mem3d) :: CLW,RNW,temp,Plev

integer :: ixy,i02,i03,i04,i05,i05c,i02gas,i04aer

integer :: i02aer


!===============================================================
! local vars
logical :: LPREC
real    :: TMASS
integer :: iwb,ieb,jsb,jeb 
integer :: KBOTC,KTOPC
real,dimension(nzz) :: clwc,rnwc,twet,pwet,RR,VOLRAT,pwr_c,cwc_c
real,dimension(nzz) :: con,con_s

real :: wcav

real,dimension(nzz) :: diam1d,density1d

real :: cwph_1d ! 2013-03-06

real :: rgf_tmp
real :: diam

real,allocatable,dimension(:,:,:) :: DEPFLD,DEPFLD2

!=================================================================

real,allocatable,dimension(:,:) :: aer_depfld,aer_depfld2

!=================================================================


! shun@20161227
integer :: ig_wdep


itt=itt_in

! allocate memory for work arrays


iwb=sx(ne)-1;ieb=ex(ne)+1
jsb=sy(ne)-1;jeb=ey(ne)+1
ALLOCATE(DEPFLD(iwb:ieb,jsb:jeb,IGAS),DEPFLD2(iwb:ieb,jsb:jeb,IGAS))

allocate(aer_depfld(iwb:ieb,jsb:jeb),aer_depfld2(iwb:ieb,jsb:jeb) )


!>-----------------------------------------------------------------

loop_j : do j = sy(ne),ey(ne)
loop_i : do i = sx(ne),ex(ne)

  ixy = (ex(ne)-sx(ne)+3)*(j-sy(ne)+1)+i-sx(ne)+1
  i02 = ip2mem(ne)

  DO k = 1, nzz
      i03 = ip3mem(k,ne)
      clwc(k) = CLW(i03+ixy)
      rnwc(k) = RNW(i03+ixy)
      TWET(k) = temp(i03+ixy)
      PWET(k) = Plev(i03+ixy)
  ENDDO ! k 

!cycle

! CCC  GET THE LAYER containin precipitation bottom/top
  CALL GETCLOUDDEPTH ( MYID,TWET,PWET,clwc,rnwc,CWC_C, PWR_C,KBOTC,KTOPC,NZZ &
                      ,RR,VOLRAT,RAINCON(i02+ixy),RAINNON(i02+ixy),LPREC)

!cycle

  ! gas
  IF(LPREC) THEN
!         DO IG = 1, IGASCBM
          do ig=1,iedgas

            i02gas = ip2memGas(ig,ne)

            ig_wdep=ig
            if(lgaschemsmp) then
              ig_wdep=idx4dep_smpsulf(ig) !18 !SO2
            endif

            IF(ITT.EQ.1) THEN
               DEPFLD(I,J,IG) =  0.0
               DEPFLD2(I,J,IG) = 0.0
            ELSE
               DEPFLD(I,J,IG) =  WETDEP  (i02gas+ixy)
               DEPFLD2(I,J,IG) = WETDEP2 (i02gas+ixy)
            ENDIF


            DO K = KTOPC, KBOTC,-1
              i03 = ip3mem(k,ne)
              i04 = ip4mem(k,ig,ne)

              con(k) = GAS(I04+IXY)

              CALL WETDEP_GAS ( MYID, KBOTC, KTOPC, DT, DX(I03+IXY),DY(I03+IXY),&
                 DZ(I03+IXY),TWET(K), PWET(K),CWC_C(K), PWR_C(K), 0., 0., VOLRAT(K), &
                 RR(K),CPH(I03+IXY),con(k), TMASS, DEPFLD(I,J,IG), DEPFLD2(I,J,IG), &
                 wcav,I, J, K, ig_wdep )

              if(ig.eq.18) then
                 gscav_so2(i03+ixy)=wcav
              endif

              GAS(I04+IXY) = con(k)
            ENDDO ! K

            WETDEP  (i02gas+ixy) =  DEPFLD(I,J,IG)
            WETDEP2 (i02gas+ixy) =  DEPFLD2(I,J,IG)
         ENDDO ! IG

  ! aerosol
if(laerv1) then
         DO IG = IGASCBM + 1, IGAS
            i02gas = ip2memGas(ig,ne)
            IF(ITT.EQ.1) THEN
              DEPFLD(I,J,IG) =  0.0
              DEPFLD2(I,J,IG) = 0.0
            ELSE
               DEPFLD(I,J,IG) =  WETDEP  (i02gas+ixy)
               DEPFLD2(I,J,IG) = WETDEP2 (i02gas+ixy)
            ENDIF


            DO K = KTOPC, KBOTC,-1
             i03 = ip3mem(k,ne)
             i04 = ip4mem(k,ig,ne)
             con(k) = GAS(I04+IXY)
             ! particle size
             IF(ig==igasCBM+2) THEN ! PM10
                diam = 1.0E-6*sqrt(2.5*10.0)*hgfac(i03+ixy)
             ELSE
                diam = 1.0E-6*sqrt(0.1*2.5)*hgfac(i03+ixy)  ! THINK OTHER AEROSOLS AS PM25
             ENDIF

             CALL WETDEP_AER_old (MYID, KBOTC, KTOPC, DT, DX(I03+IXY), DY(I03+IXY),&
                 DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k) ,VOLRAT(K), &
                 RR(K),con(k),diam, TMASS, DEPFLD(I,J,IG), DEPFLD2(I,J,IG), &
                 wcav ,I, J, K, IG   )

             
              if(ig.eq.75) then
                 ascav_so4(i03+ixy)=wcav
              endif


             GAS(I04+IXY) = con(k)

            ENDDO ! K

            WETDEP  (i02gas+ixy) =  DEPFLD(I,J,IG)
            WETDEP2 (i02gas+ixy) =  DEPFLD2(I,J,IG)

         ENDDO ! IG
endif

if(laerv2) then
         do ia=1,naersp
!         if(ia.eq.13) cycle ! skip hydrophobic bc
         if(ia.eq.13) then
            laerincld=.false.
         else
            laerincld=.true.
         endif         

         do is=1,naerbin

            if(ia.gt.1.and.is.gt.1) cycle ! skip zero aersol tracer

            i02aer=ip2memaer(is,ia,ne) 
     
            if(itt.eq.1) then
              aer_depfld(i,j)  = 0.0
              aer_depfld2(i,j) = 0.0
            else
              aer_depfld(i,j)  = aerom_wdepfld(i02aer+ixy)
              aer_depfld2(i,j) = aerom_wdepfld2(i02aer+ixy)
            endif


            do k=KTOPC, KBOTC,-1
              i03 = ip3mem(k,ne)
              i04aer = ip4mem_aer(k,is,ia,ne)
              con(k)=aerom(i04aer+ixy)
              diam=1.0E-6*diamlgm(is)*hgfac(i03+ixy)
              IG=77
              CALL WETDEP_AER (MYID, KBOTC, KTOPC, DT, DX(I03+IXY),DY(I03+IXY),&
                 DZ(I03+IXY),TWET(K), PWET(K),clwc(k),rnwc(k) ,VOLRAT(K), &
                 RR(K),con(k),diam, TMASS, aer_depfld(i,j), aer_depfld2(i,j), &
                 wcav ,I, J, K, aer_rhop(ia) )
              aerom(i04aer+ixy)=con(k)
            enddo

            aerom_wdepfld(i02aer+ixy)=aer_depfld(i,j) 
            aerom_wdepfld2(i02aer+ixy)=aer_depfld2(i,j)

         enddo
         enddo

endif


  ELSE
         DO IG = 1, IGAS
             i02gas = ip2memGas(ig,ne)
!             WETDEP  (i02gas+ixy) =  -1.E20
!             WETDEP2 (i02gas+ixy) =  -1.E20
         ENDDO

  ENDIF ! LPRC

enddo loop_i
enddo loop_j


!return

!$$$$$$$$$$$$$$$$$$$$$$$$$
! dust and sea salt

 DO K=1,11   ! WET DEP

     i03 = ip3mem(k,ne )

     DO IA=1,IAER
     DO IS=1,ISIZE

        I05=IP5MEM(K,IS,IA,NE)
        I02 = IP2MEM(NE)

        IF( IA == 2 ) THEN
          !!! TO CALCULATE THE DUST WET DEPOSITION UG/M2/HR
          CALL DUSTWETDEP(MYID, DUSTWET(I02),AER(I05),RAINCON(I02),RAINNON(I02),&
                          SX(NE),EX(NE),SY(NE),EY(NE),DT,K,dz(i03))
          do iduc = 1 , ndustcom
            i05c = ip5memc (k,is,iduc,ne)
            IF (iduc ==7 ) then
              CALL DUSTWETDEP( MYID,DUSTWETSO4(I02),DUSTCOMP(I05C) &
                              ,RAINCON(I02),RAINNON(I02) &
                              ,SX(NE),EX(NE),SY(NE),EY(NE),DT,K,dz(i03))
            ELSE IF(iduc == 8) then
              CALL DUSTWETDEP( MYID,DUSTWETNO3(I02),DUSTCOMP(I05C) &
                              ,RAINCON(I02),RAINNON(I02) &
                              ,SX(NE),EX(NE),SY(NE),EY(NE),DT,K,dz(i03) )
            ELSE  IF(iduc==9) then
              CALL DUSTWETDEP( MYID,DUSTWETFeII(I02),DUSTCOMP(I05C) &
                              ,RAINCON(I02),RAINNON(I02) &
                              ,SX(NE),EX(NE),SY(NE),EY(NE),DT,K,dz(i03) )
            ELSE IF(iduc==10) then
              CALL DUSTWETDEP( MYID,DUSTWETFeIII(I02),DUSTCOMP(I05C) &
                              ,RAINCON(I02),RAINNON(I02) &
                              ,SX(NE),EX(NE),SY(NE),EY(NE),DT,K,dz(i03))
            ENDIF
          enddo

        ENDIF ! IA IF

        CALL  WET_DEP_DUST( MYID,AER(I05),RAINCON(I02),RAINNON(I02) &
                           ,SX(NE),EX(NE),SY(NE),EY(NE),DT )

        IF(IA==1)  THEN ! SEA SALT
          do iduc = 1 ,nseacom
            i05c = ip5memcs (k,is,iduc,ne)
            CALL WET_DEP_DUST( MYID,SEACOMP(I05C),RAINCON(I02),RAINNON(I02) &
                              ,SX(NE),EX(NE),SY(NE),EY(NE),DT )
          enddo

        ELSE IF(IA==2) THEN
          do iduc = 1 ,ndustcom
            i05c = ip5memc (k,is,iduc,ne)
            CALL WET_DEP_DUST( MYID,DUSTCOMP(I05C),RAINCON(I02),RAINNON(I02) &
                              ,SX(NE),EX(NE),SY(NE),EY(NE),DT )
          enddo
        ENDIF

     ENDDO        !ISIZE
     ENDDO        !IAER

  ENDDO !K


  deallocate(aer_depfld,aer_depfld2,depfld,depfld2)

end subroutine naqpms_wet_dep


