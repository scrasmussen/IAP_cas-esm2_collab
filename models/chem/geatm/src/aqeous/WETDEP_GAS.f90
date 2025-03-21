         SUBROUTINE WETDEP_GAS (MYID,KBOT, KTOP,DELTAT,DELTAX,DELTAY,DEPTH,&
            TEMPK,PRESS,CWC, PWR, PWS, PWG,VOLRAT, RR,CPH, CONC,TMASS,&
             DEPFLD,DEPFLD2,wcav,I,J,K,IG)

!   THIS PROGRAM IS FROM CAMX FOR WET DEPTION OF GAS  
!     WETDEP modifies vertical concentration profiles for a given grid via 
!     precipitation processes.  This subroutine has been completely rewritten
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


!     Output arguments: 
!        conc                concentration field (umol/m3, ug/m3) 
!        depfld              2-D array of wet deposited mass (mol/ha,g/ha)
!        depfld2             and surface liquid concentrations (mol/l,g/l)
      
         IMPLICIT NONE

         real :: wcav
         logical :: laerincld=.false.  ! not used in gas wdep 
 
         integer :: nspcs,igrid,j,jycl,i1,i2,i,ixcl,kbot,ktop, &
             k,ncnt,kzcl,l,isemptyf,isemptyc,kwtr,isec,isempty,&
             iaero,ll,knh3,khno3,kso2,ko3,ig,MYID
         real :: deltax,tempk,press,cwc,pwr,pws,pwg,cph,depth
         real conc,depfld,depfld2
         real*8 fluxes(11)        
         real c0,pp,rr,volrat,tmass
         real delr
         real rd,rhoh2o,deltat,deltay,densfac,dtout,cellvol,rainvol,&
              rhoair, delc,delm,cmin,hlaw,gscav,ascav,c00,totc,&
              totw,ceq,qtf,qtc,vtf,vtc,psizec,ascavf,roprta,psize,  &
              ascavc,qt,vt,cwat,pwat,rconst,cwmin,tamin
         real convfac   !          conversion factor: umol/m3 = ppm * convfac
         REAL  henry0(74),tfact(74),f0(74),rscale(74), diffrat(74)

         real cwc_c,pwr_c        

         logical lcloud,ltop,lgraupl

         data rd /287./         ! Dry air gas constant (J/K/kg)
         data rhoh2o /1.e6/     ! water density (g/m3)
         data rconst /8.206e-2/ ! gas constant (l.atm/mol.K)
         data cwmin /0.05/
         DATA TAMIN / 243./
         data henry0&
          /1.00e+10,  2.10e+05, 1.00e+05, 5.76e+01, 1.90e-03, 1.00e-02,&
           0.0,       3.20e+04, 5.90e+01, 0.0     , 1.10e-02, 0.0,&
           0.0,       0.0     , 0.0,      7.40e+04, 1.00e-10, 1.22e+00,&
           0.0,       1.00e-03, 0.0,      0.0,      6.30e+03, 0.0,&
           0.0,       0.0,      2.20e+02, 6.30e+03, 0.0,      0.0,&
           0.0,       3.60e+00, 1.00e-03, 0.0,      2.70e+03, 1.00e-02,&
           5.00e-03,  5.00e-03, 1.20e+00, 1.40e+00, 2.70e+03, 0.0,&
           0.0,       2.70e+03, 9.40e+03, 0.0,      0.0,      0.0,&
           0.0,       0.0,      1.00e-03, 1.00e-02, 6.30e+03, 6.30e+03,&
           6.30e+03,  6.30e+03, 0.0,      0.0,      0.0,      0.0,&
           0.0,       0.0,      0.0,      0.0,      0.0,      0.0,&
           1.00e+10, 2.7e+03,  2.7e+03,  2.7e+03, 2.7e+03,   2.7e+03,&
           2.7e+03,   2.7e+03/
         data tfact&
          /0.0,    -8707., 0.0,    -4100.,  -1480.,  -2516.,&
           0.0,    -8706., -4781., 0.0,     -2415.,  0.0,&
           0.0,    0.0,    0.0,    -6643.,  0.0,     -3156.,&
           0.0,    0.0,    0.0,     0.0,    -6492.,  0.0,&
           0.0,    0.0,    -4932., -6492.,  0.0,     0.0,&
           0.0,    -5910., 0.0,    0.0,     -6492.,  0.0,&
           0.0,    0.0,    0.0,    0.0,     -6492.,  0.0,&
           0.0,    -6492., -8706., 0.0,     0.0,     0.0,&
           0.0,    0.0,    0.0,    0.0,     -6492.,  -6492.,&
         -6492., -6492., 0.0,    0.0,      0.0,      0.0,&
           0.0,    0.0,    0.0,    0.0,      0.0,      0.0,&
           0.0,    -6492., -6492., -6492.,  -6492.,  -6492.,&
         -6492., -6492./
         data rscale&
         /0.0,    0.0,    0.0,    0.0,      1.0,      1.0,&
          1.0,    0.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0,    1.0,    1.0,      1.0,      1.0,&
          1.0,    1.0 /
        data diffrat&
         /1.00,   1.87,   1.42,   0.97,     1.29,    1.29,&
          0.0,    2.45,   1.62,   0.0,      1.63,    0.0,&
          0.0,    0.0,    0.0,    1.37,     1.25,    1.89,&
          0.0,    2.00,   0.0,    0.0 ,     1.29,    0.0,&
          0.0,    0.0,    1.60,   1.56,     0.0,     0.0,&
          0.0,    2.59,   2.00,   0.0,      2.00,    1.25,&
          1.80,   1.80,   2.26,   2.43,     2.45,    0.0,&
          0.0,    2.47,   2.72,   0.0,      0.0,     0.0,&
          0.0,    0.0,    2.00,   1.94,     1.97,    1.97,&
          1.97,   1.97,   0.0,    0.0,      0.0,     0.0,&
          0.0,    0.0,    0.0,    0.0,      0.0,     0.0,&
          1.00,   2.50,   2.50,   2.50,      2.50,   2.50,&
          2.50,   2.50/

!-----Entry point

!--   cpnvert kg/kg to g/m3 20150728
           convfac = 44.9 * (273./TEMPK)*(PRESS/1013.)
           cwc_c  = cwc  * 29. * convfac
           pwr_c  = pwr  * 29. * convfac
           tmass = 0.0

          ltop= .false. ! juanxiong he
          IF(K.EQ.KTOP)   ltop = .true.
          lcloud  = .false.
          lgraupl = .false.    
          if (cwc_c.ge.cwmin) lcloud = .true.
!          if (pwg.ge.cwmin) lgraupl = .true.
          cellvol = deltax*deltay*depth
          rainvol = volrat*cellvol
          rhoair = 100.*press/(rd*tempk)

!-----Calculate scavenging for soluble gas species
          knh3  =  4
          khno3 =  2
          kso2  = 18
          ko3   = 11

          
          if( henry0(ig).LT.1.e-6 ) then
          tmass   = 0.0
          depfld  = 0.0
          depfld2 = 0.0             
          goto 40      
          endif

          IF(cph.eq.-1.e20) cph = 5.6

          call henryfnc(ig,henry0(ig),tfact(ig),tempk, &
                cph,knh3,khno3,kso2,hlaw)

          hlaw = hlaw*rconst*tempk
          cwat = cwc
          pwat = pwr + pws + pwg

          if (tempk.lt.273. .and. rscale(ig).gt.0.) then
            cwat = amax1(0.,cwc* &
                         (tempk - tamin)/(273. - tamin))
            pwat = pwr
          endif

          delc = 0.
          delm = 0.
          convfac = 44.9 * (273./tempk)*(press/1013.)
!    UNIT CHANGE FROM PPB TO UMOL/M3 
          conc = conc * convfac * 1.E-3

          if (ltop) tmass = 0.

          c0 = tmass/rainvol          
          cmin = 1.E-20*convfac
          conc = amax1(cmin,conc)


          call scavrat(.false.,lcloud,lgraupl,tamin,rr,&
                         tempk,cwat,depth,rhoair,&
                         conc,hlaw,diffrat, &
                         rscale,0.,0.,gscav,ascav,laerincld) 


          wcav=gscav

          c00  = c0*volrat
          totc = conc + c00
          totw = cwat + pwat
          ceq  = totc/(1. + hlaw*totw/rhoh2o)
          ceq  = totc - ceq
          delc = (ceq - c00)*(1. - exp(-gscav*deltat))


          if (delc.gt.0.) delc = amin1(delc,conc-cmin)
          conc = conc - delc
          delm = delc*cellvol
          tmass= tmass + delm


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

!    UNIT CHANGE FROM UMOL/M3 to PPB
          conc = conc / convfac /1.E-3
  40   continue
        
         RETURN
         ENDSUBROUTINE
