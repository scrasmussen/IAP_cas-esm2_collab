           SUBROUTINE GETCLOUDDEPTH ( myid,T, PRESS, cwc, pwr, CWC_C, PWR_C,KBOTC, KTOPC, NZZ,&
                       RR,VOLRAT,RAINN,RAINC,LPREC, ix, iy)
!c-----Scan column for layers containing precipitation bottom/top
           INTEGER :: MYID,KBOTC , KTOPC ! THE BOTTOM AND TOP LAYER OF containing precipitation
           REAL    :: cwc(NZZ) ! CLOUD WATER CONTENT  (g/m3)
           REAL    :: pwr(NZZ) !  rain water content (g/m3)

           REAL    :: pws(nzz),pwg(nzz)
           REAL    :: pws_c(nzz),pwg_c(nzz)


           REAL    :: T(NZZ), PRESS(NZZ) ! IN K and HPA
           REAL    :: convfac ! conversion factor: umol/m3 = ppm * convfac   
           LOGICAL :: LPREC
           REAL    :: RR(NZZ),VOLRAT(NZZ)
           REAL    :: CWC_C(NZZ),PWR_C(NZZ) ! g/m3
           REAL    :: CWMIN,rhoh2o
           REAL    :: RAINN,RAINC ! WRF OUTPUT NON-CONVECTION and CONVECTION  RAIN (mm/hr)

           integer :: ix,iy
           logical :: lwrf_prec

           real    :: pp(nzz),dz(nzz),pcirsg(nzz),prsg(nzz)
           real    :: wrfrain_cm,wrfrain_1d,ppcol,pccol

           DATA CWMIN/ 0.05 /
           data rhoh2o /1.e6/     ! water density (g/m3) 

           ! to be modified
           pws=0.0
           pwg=0.0

           wrfrain_cm=RAINN+RAINC

           lwrf_prec=.false.
           if(wrfrain_cm.gt.0.0) then
             lwrf_prec=.true.
           endif


           LPREC = .false.

           if(lwrf_prec) then  ! if precipitation exist from wrf
            
           LPREC = .TRUE.
           ktopc = 0
           kbotc = 0
           
           do k = 1, nzz
!--   cpnvert kg/kg to g/m3
           convfac = 44.9 * (273./t(k))*(press(k)/1013.)         
           cwc_c (k) = cwc (k) * 29. * convfac
           pwr_c (k) = pwr (k)* 29. * convfac           
           pws_c(k)=pws(k)* 29. * convfac
           pwg_c(k)=pwg(k)* 29. * convfac
           enddo
 
           do k = 1, nzz
            IF( pwr_c(k).ge.cwmin.or.pws_c(k).ge.cwmin.or. &
              & pwg_c(k).ge.cwmin ) then
             kbotc = k
             goto 25
            ENDIF
           enddo
           
           LPREC = .FALSE.
           GOTO 20 
  25       CONTINUE
           IF (kbotc.eq.nzz) THEN
           LPREC = .FALSE.
           goto 20       
           ENDIF

           NCNT = 1
           DO K = KBOTC+1, NZZ
            IF( pwr_c(k).LT.cwmin.and.pws_c(k).LT.cwmin.and. &
                pwg_c(k).LT.cwmin ) then
             KTOPC = K - 1
             GOTO 26
            ENDIF
            NCNT = NCNT + 1
           ENDDO     
           KTOPC = NZZ
  26       CONTINUE
           IF (KBOTC.gt.1 .and. ncnt.eq.1) THEN
            LPREC = .FALSE. 
            goto 20    
           ENDIF
      
  20       CONTINUE ! no precipitation water

           endif ! if precipitation exist from wrf


!          IF(LPREC.AND.KBOTC.EQ.1.AND.KTOPC.GE.5) PRINT*,KBOTC,KTOPC
          do k=1,nzz
            volrat(k) = 0.
            rr(k) = 0.
            prsg(k)  =(pwr_c(k)+pws_c(k)+pwg_c(k))*dz(k)
            pcirsg(k)=(pwr_c(k)+pws_c(k)+pwg_c(k)+cwc_c(k))*dz(k) !+ cwc_i(k)
          enddo


          if(LPREC) then
            ppcol=0.0
            pccol=0.0
            do k=kbotc,ktopc
              ppcol=ppcol+prsg(k)
              pccol=pccol+pcirsg(k)
            enddo
            do k = kbotc,ktopc
!              rr(k) = (volrat(k)/1.0e-7)**1.27   ! rainfall rate (mm/hr)
              if(ppcol.gt.0) then
                volrat(k) = prsg(k)/dz(k)/rhoh2o  ! drop volume/air volume
                rr(k) = wrfrain_cm*10.0*prsg(k)/ppcol
              elseif(pccol.gt.0)
                rr(k) = wrfrain_cm*10.0*pcirsg(k)/pccol
                volrat(k)=pcirsg(k)/dz(k)/rhoh2o
              else
                ! reconstruct rain water content ???
                rr(k) = 0.0
                volrat(k)=0.0
              endif
            enddo
          endif
 
  30       CONTINUE         

            
           RETURN
           ENDSUBROUTINE
