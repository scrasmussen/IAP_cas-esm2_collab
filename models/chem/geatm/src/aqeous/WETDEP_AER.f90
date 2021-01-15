         SUBROUTINE WETDEP_AER (MYID,KBOT,KTOP,DELTAT,DELTAX,DELTAY,DEPTH,&
            TEMPK,PRESS,cwc,pwr,VOLRAT, RR, CONC,PSIZE,TMASS,&
             DEPFLD,DEPFLD2,wcav,I,J,K,denaer)

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

         use wdepflags, only: laerincld
        
         IMPLICIT NONE
  
         real ::denaer

         real :: wcav

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

           
          ltop = .false. ! juanxiong he
          IF(K.EQ.KTOP)   ltop = .true.
          lcloud  = .false.
          lgraupl = .false.
          if (cwc_c.ge.cwmin) lcloud = .true.
!          if (pwr_c.ge.cwmin) lgraupl = .true.
          cellvol = deltax*deltay*depth
          rainvol = volrat*cellvol
          rhoair = 100.*press/(rd*tempk)

          rhop = denaer


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
    
