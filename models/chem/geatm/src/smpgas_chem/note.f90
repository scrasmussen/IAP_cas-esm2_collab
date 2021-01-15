

h2o      = WaterVapor(RH, cair_mlc, te, pr_atm)


      function WaterVapor(rh, cair_mlc, te, pr_atm)

      t_steam = 373.15                  ! steam temperature  [K]
      pr_std   = 1.0                    ! standard pressure  [atm]

      a      = 1.0 - t_steam/te
      arg    = (((-.1299*a -.6445)*a -1.976)*a +13.3185)*a
      pr_h2o = pr_std*exp(arg)                          ! [atm]
      WaterVapor = RH*(pr_h2o/pr_atm)*cair_mlc/100.     ! [molec/cc]

      return
      end



      RK_2HO2    = 2.3e-13 * exp(600./te) +     ! ho2 + ho2 --> h2o2
     &             1.7e-33 * exp(1000./te)*cair_mlc                           
      rk_com(31) = RK_2HO2

      RK_2HO2_H2O= RK_2HO2*1.4e-21*exp(2200./te)! ho2 + ho2 + h2o --> h2o2                                                          
      rk_com(32) = RK_2HO2_H2O






      cair_mlc = avogad*pr_atm/(82.056*te) ! air conc [molec/cc]



      rk0 = 3.0e-31
      rnn = 3.3
      rki = 1.5e-12
      rmm = 0.0
      rk_com(45) = Troe(cair_mlc,te,rk0,rnn,rki,rmm)

end



      Function Troe(cair_mlc,te,rk0,rnn,rki,rmm)
      rk0 = rk0*cair_mlc*(te/300.)**(-rnn)
      rki = rki*(te/300.)**(-rmm)
      expo= 1./(1. + (ALOG10(rk0/rki))**2)
      troe  = (rk0*rki/(rk0+rki))*.6**expo
      return
      end
