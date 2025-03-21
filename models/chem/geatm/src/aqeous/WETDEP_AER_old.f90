      SUBROUTINE WETDEP_AER_old (MYID,KBOT,KTOP,DELTAT,DELTAX,DELTAY,DEPTH,&
            TEMPK,PRESS,cwc,pwr,VOLRAT, RR, CONC,PSIZE,TMASS,&
             DEPFLD,DEPFLD2,wcav,I,J,K,IG)

!   THIS PROGRAM IS FROM CAMX FOR WET DEPTION OF GAS  
!     WETDEP modifies vertical concentration profiles for a given grid
!     via 
!     precipitation processes.  This subroutine has been completely
!     rewritten
!     for CAMx v4. 
!     Input arguments:

!        deltat              time step (s)
!        deltax              cell size in x-direction (m)
!        deltay              cell size in y-direction (m)
!        mapscl              map scale factor
!        depth               cell depth (m)
!        tempk               temperature field (K)
!        press               pressure field (mb)
!        cwc                 cloud water content (g/m3)
!        pwr                 rain water content (g/m3)
!        pws                 snow water content (g/m3)
!        pwg                 graupel water content (g/m3)
!        cph                 cloud water pH
!        conc                concentration field (ppbv)
!        ktop                the top layer of precipatatiom 
!        kbot                thr bottom layer of precipation 
!        volrat              drop volume/air volume
!        rr                  rainfall rate (mm/hr)
!        tmass               amounts umol in rainfall dropets
!        PSIZE               

!     Output arguments: 
!        conc                concentration field (umol/m3, ug/m3) 
!        depfld              2-D array of wet deposited mass (g/ha)
!        depfld2             and surface liquid concentrations (g/l)
        
         IMPLICIT NONE


         real :: wcav
         logical :: laerincld=.true.

         integer :: nspcs,igrid,j,jycl,i1,i2,i,ixcl,kbot,ktop, &
             k,ncnt,kzcl,l,isemptyf,isemptyc,kwtr,isec,isempty,&
             iaero,ll,knh3,khno3,kso2,ko3,ig,MYID
         real :: deltax,tempk,press,cwc,pwr,pws,pwg,cph,depth
         real conc,depfld,depfld2
         real c0,pp,rr,volrat,tmass
         real delr
         real rd,rhoh2o,deltat,deltay,densfac,dtout,cellvol,rainvol,&
              rhoair, delc,delm,cmin,hlaw,gscav,ascav,c00,totc,&
              totw,ceq,qtf,qtc,vtf,vtc,psizec,ascavf,rhop,psize,  &
              ascavc,qt,vt,cwat,pwat,rconst,cwmin,tamin
         real convfac   !          conversion factor: umol/m3 = ppm *convfac

         real cwc_c,pwr_c
         logical lcloud,ltop,lgraupl

         data rd /287./         ! Dry air gas constant (J/K/kg)
         data rhoh2o /1.e6/     ! water density (g/m3)
         data rconst /8.206e-2/ ! gas constant (l.atm/mol.K)
         data cwmin /0.05/
         DATA TAMIN / 243./

!-----Entry point

!--   cpnvert kg/kg to g/m3 20151014
           convfac = 44.9 * (273./TEMPK)*(PRESS/1013.)
           cwc_c  = cwc  * 29. * convfac
           pwr_c  = pwr  * 29. * convfac
           tmass = 0.0



          IF(K.EQ.KTOP)   ltop = .true.
          lcloud  = .false.
          lgraupl = .false.
          if (cwc_c.ge.cwmin) lcloud = .true.
!          if (pwr_c.ge.cwmin) lgraupl = .true.
          cellvol = deltax*deltay*depth
          rainvol = volrat*cellvol
          rhoair = 100.*press/(rd*tempk)


      IF(ig == 75 ) rhop = 1.0E06 ! Primary PM25
      IF(ig == 76 ) rhop = 1.0E06 ! Primary PM10
      IF(ig == 77 ) rhop = 2.0E06 ! BC
      IF(ig == 78 ) rhop = 1.0E06 ! Primary OC
      IF(ig == 79 ) rhop = 0.0    ! H+
      IF(ig == 80 ) rhop = 2.0E06 ! Na+
      IF(ig == 81 ) rhop = 1.5E06 ! NH4+
      IF(ig == 82 ) rhop = 2.0E06 ! CL-
      IF(ig == 83 ) rhop = 1.5E06 ! SO42-
      IF(ig == 84 ) rhop = 1.5E06 ! HSO4-
      IF(ig == 85 ) rhop = 1.5E06 ! NO3-
      IF(ig == 86 ) rhop = 2.0E06 ! NACL
      IF(ig == 87 ) rhop = 1.5E06 ! NA2SO4
      IF(ig == 88 ) rhop = 1.5E06 ! NANO3
      IF(ig == 89 ) rhop = 1.5E06 ! NH42SO4
      IF(ig == 90 ) rhop = 1.5E06 ! NH4NO3
      IF(ig == 91 ) rhop = 1.5E06 ! NH4CL
      IF(ig == 92 ) rhop = 1.5E06 ! H2SO4
      IF(ig == 93 ) rhop = 1.5E06 ! NH4HSO4
      IF(ig == 94)  rhop = 1.5E06 ! NAHSO4
      IF(ig == 95)  rhop = 1.5E06 ! (NH4)4H(SO4)2(S)
      IF(ig == 96)  rhop = 1.0E06 ! SOA1
      IF(ig == 97)  rhop = 1.5E06 ! SOA2
      IF(ig == 98)  rhop = 1.5E06 ! SOA3
      IF(ig == 99)  rhop = 1.5E06 ! SOA4
      IF(ig == 100) rhop = 1.5E06 ! SOA5
      IF(ig == 101) rhop = 1.5E06 ! SOA6
      IF(ig == 102) rhop = 1.5E06 ! AH2O 


          call scavrat( .true., lcloud, lgraupl, tamin, rr, tempk, 0.,&
                       depth, rhoair, 0., 0., 0., 0., psize, rhop,&
                       gscav, ascav,laerincld)


          wcav=ascav

          delr = 1. - exp(-ascav*deltat)    

          delc = 0.
          delm = 0.
         
          if (ltop) tmass = 0.
          cmin = 1.E-20
          conc = amax1(cmin, conc)
          delc = conc * delr
          delc = amin1(delc, conc-cmin)
          
          conc = conc - delc
          delm = delc * cellvol
          tmass = tmass + delm
 
          ltop = .false.

          dtout = 60. ! minutes out frequency
!-----If rain evaporates before reaching the ground, return all mass
!     back to layer KBOT
          IF(K.EQ.KBOT) THEN
           if (kbot.gt.1) then
             cellvol = deltax*deltay*depth
             conc = conc + tmass/cellvol
             tmass = 0.0
            else
            depfld  = depfld  + 1.e-2*tmass/(deltax*deltay) !
            depfld2 = depfld2 + 1.e-9*(tmass/rainvol)*deltat/(60.*dtout)

           endif
          ENDIF


         RETURN
         ENDSUBROUTINE
    
