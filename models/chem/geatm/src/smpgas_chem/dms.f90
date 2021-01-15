      rk_tot_num =       te * exp(-234./te) +
     &             8.46e-10 * exp(7230./te) +
     &             2.68e-10 * exp(7810./te)
      rk_tot_den = 1.04e+11 * te + 88.1 * exp(7460./te)
      rk_tot     = rk_tot_num/rk_tot_den
c
      rk_mar(1)   = 9.60e-12 * exp(-234./te) ! ch3sch3 + oh --> ch3sch2
      Babs       = rk_mar(1)/rk_tot
      Badd       = 1. - Babs
      rk_mar(2)   = 1.40e-13 * exp(500./te)  ! ch3sch3 + no3 --> 
      rk_mar(3)   = 1.26e-11 * exp(409./te)  ! ch3sch3 + o3p -->


      rk_mar(4)   = Badd*rk_tot




      p_mar(idms)= 0.0
      d_mar(idms)= r_mar(1)+r_mar(2)+r_mar(3)+r_mar(4)




      r_mar(1)    = rk_mar(1)   *s(idms)*s(ioh)
      r_mar(2)    = rk_mar(2)   *s(idms)*s(ino3)
      r_mar(3)    = rk_mar(3)   *s(idms)*s(io3p)
      r_mar(4)    = rk_mar(4)   *s(idms)*s(ioh)







