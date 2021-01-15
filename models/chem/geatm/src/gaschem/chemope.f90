      SUBROUTINE CHEMOPE(MYID, OPE, DT, I, J, K, NE)
!C  TO CALCULATE THE NET OPE (OZONE PRODUCTION EFFICIENCY)FOR O3      
!C     OPE = (P(O3)-P(O3))/P(NOZ) 
   
      include 'chm1.inc'
      include 'gas1.inc'
      integer myid,i,j,k,ne
      real    DT,OPE,PO3,LO3,PNOZ

!ozone production and loss,and loss due to nox,radicals from radical-radical reactions
!P(o3)=k4[NO][HO2]+k5[NO][CH3O2]+k6[RO2][NO]
!L(o3)=k3[O1D][H2O]+k8[O3][HO2]+k7[O3][OH]+(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
      
      PO3  = 60.*dt*(rk_com(33)*cnn(kno)*cnn(kho2)& !k4[NO][HO2] inmolec/cc
            +rk_com(57)*cnn(kch3o2)*cnn(kno)& !k5[NO][CH3O2])
            +rk_urb(17)*cnn(kto2)*cnn(kno)&   !k[TO2][NO]
            +rk_com(58)*cnn(kethp)*cnn(kno)&  !K[ETHP][NO]
            +rk_urb(28)*cnn(kro2)*cnn(kno)&   !K[RO2][NO]
            +rk_com(71)*cnn(kc2o3)*cnn(kno)&  !K[c2o3][no]
            +rk_urb(29)*cnn(kano2)*cnn(kno)&  !k[ano2][no]
            +rk_urb(30)*cnn(knap)*cnn(kno)&   !k[nap][no]
            +rk_urb(31)*cnn(kxo2)*cnn(kno)&   !k[xo2][no]
            +rk_bio(8)*cnn(kisopp)*cnn(kno)&  !k[isopp][no]
            +rk_bio(9)*cnn(kisopn)*cnn(kno)&  !k[isopn][no]
            +rk_bio(10)*cnn(kisopo2)*cnn(kno)&!k[isopo2][no]
            )/cair_mlc*1.e+9

 
      LO3 = 60.*dt*(rk_com(12)*cnn(ko1d)*h2o& ! k3[O1D][H2O]
            +rk_com(21)*cnn(ko3)*cnn(kho2)&  ! k8[O3][HO2]
            +rk_com(20)*cnn(ko3)*cnn(koh)&   ! k7[O3][OH]
            +(rk_com(18)*cnn(kno)*cnn(ko3)*& !k10[NO][O3]
             rk_com(24)*cnn(kno2)*cnn(koh)/& !*k12[NO2][OH]
            (rk_com(1)*cnn(kno2)+rk_com(24)*cnn(kno2)*cnn(koh))& !k11[NO2]+k12[NO2][OH]
            )&      !(k10[NO][O3]*k12[NO2][OH]/{k11[NO2]+k12[NO2][OH]})
            )/cair_mlc*1.e+9

! CCCCCCCCCCC   NOZ PRODUCTION
! OH + NO2       = HNO3    
! NO3 + NO2      = N2O5      
! C2O3 + NO2     = PAN           
! CRO + NO2      = ONIT


      PNOZ =  60.*dt*(rk_com(24)*cnn(koh)*cnn(kno2) & ! OH +NO2
                  + rk_com(39)*cnn(kno3)*cnn(kno2) & ! NO3+NO2
                  + rk_com(69)*cnn(kc2o3)*cnn(kno2)& ! C2O3+NO2
                  + rk_urb(20)*cnn(kcro)*cnn(kno2) &  ! CRO +NO2 
                  ) / cair_mlc*1.e+9

      IF ( (PO3-LO3) .GE. 1.E-20 .AND. PNOZ.GE.1.E-20) THEN
       OPE = (PO3 - LO3) / PNOZ
      ELSE
       OPE = -1.E20
      ENDIF

      RETURN 
      ENDSUBROUTINE
