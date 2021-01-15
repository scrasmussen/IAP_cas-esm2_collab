        SUBROUTINE CALCLF2(MYID,CLFLO,CLFMI,CLFHI,FCLD,RAINCV,RAINNCV,COPD,&
                       HE,TER,LON,LAT,TBEG_DD,TBEG_HH,TBEG_MM,TBEG_SS,&
                       I,J,K )
         INTEGER   :: MYID,SX,EX,SY,EY,NE
         REAL      :: CLFLO,CLFMI,CLFHI,CLF,FCLD,CTYPE,COPD !COPD:OPITICAL DEPTH
         REAL      :: H,TER,HE,LON,LAT,ITOP  ! ITOP IS THE TYPE OF CLOUD 
         REAL      :: cos_sa,tmid_sec1 !cosine of zenith angle 
         REAL      :: TBEG_DD,TBEG_HH,TBEG_MM,TNEG_SS
         REAL      :: HTOP(7) ! HEIGHT OF CLOUD M
         REAL      :: TR,AC,RAINCV,RAINNCV
         DATA HTOP / 2000., 4000., 6000., 4000., 6000., 6000., 6000./

         FCLD = 1.0
         tmid_sec1  = ((tbeg_dd*24.+tbeg_hh)*60.+tbeg_mm)*60.+tbeg_ss 
   
         CALL  ZenithAngle(tmid_sec1,LON,LAT,cos_sa)
 
         COS_SA = MAX(0.5,COS_SA)
         
         CLF = MAX(CLFLO,CLFMI,CLFHI)
         
         IF(CLFLO.GE.0.5.AND.CLFMI.GE.0.5.AND.CLFHI.GE.0.5) THEN
           CTYPE = 7.      
           ITOP  = 7
         ELSE IF(CLFLO.GE.0.5.AND.CLFMI.GE.0.5) THEN
           CTYPE = 6.
           ITOP  = 4
         ELSE IF(CLFLO.GE.0.5.AND.CLFHI.GE.0.5) THEN     
           CTYPE = 5.
           ITOP  = 5
         ELSE IF(CLFMI.GE.0.5.AND.CLFHI.GE.0.5) THEN
           CTYPE = 4.      
           ITOP  = 6
         ELSE IF(CLF == CLFLO.AND.CLF.GE.0.1) THEN
           CTYPE = 3.      
           ITOP  = 1
         ELSE IF(CLF == CLFMI.AND.CLF.GE.0.1) THEN 
           CTYPE = 2.
           ITOP  = 2 
         ELSE IF(CLF == CLFHI.AND.CLF.GE.0.1) THEN
           CTYPE = 1.
           ITOP  = 3 
         ELSE 
           CTYPE = 0.0     
           ITOP  = 3 
         ENDIF

!         IF(CTYPE == 7 ) COPD = 82.
!         IF(CTYPE == 6 ) COPD = 80.
!         IF(CTYPE == 5 ) COPD = 52.
!         IF(CTYPE == 4 ) COPD = 32.
!         IF(CTYPE == 3 ) COPD = 50.
!         IF(CTYPE == 2 ) COPD = 30.
!         IF(CTYPE == 1 ) COPD = 2.

         IF(COPD>=5.) THEN
          TR = (5.- 1./EXP(COPD))/(4.+3.*COPD*(1.-0.86)) ! THE ENEGERGY TRANSMISSION COFFICIENT
          H = HE - TER
          IF(H.GT.HTOP(ITOP)) THEN         
           AC = 1. +  (1.-TR)*COS_SA
          ELSE  
           AC = 1.6*TR*COS_SA 
          ENDIF
           FCLD = 1.+ CLF*(AC -1. )    
         ELSE
           FCLD = 1.0
         ENDIF       
        RETURN
        END        

        subroutine ZenithAngle(tmid_sec1,lon,lat,cos_sa)
          real :: tmid_sec1,rlon,rlat,cos_sa,lon,lat
          real :: tlocal,tdec,codec,sidec,tloc,thou
          data deg2rad/0.017453293/ 
          rlon = lon*deg2rad
          rlat = lat*deg2rad
          tlocal=tmid_sec1
          tdec=0.4092797*sin(1.992385E-7*tlocal)
          sidec=sin(tdec)
          codec=cos(tdec)
          tloc=7.272205E-5*tlocal
          thou=cos(rlon+tloc)
          cos_sa=sin(rlat)*sidec+cos(rlat)*codec*thou 
        return
        end
