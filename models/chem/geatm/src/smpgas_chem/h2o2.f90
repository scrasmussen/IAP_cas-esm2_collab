c h2o2 + hv --> 2 oh
      alpha0=2.540534e-9-4.e-11*(z-0.75)**2.
      alpha8=-3.e-5+sqrt(alpha0)
      beta0=.539284-.16*(z-1.75)**2.
      beta8=-1+sqrt(beta0)
      rk_photo(jphoto_h2o2)   = alpha8*exp(beta8/cos_sza)







r_com(30) = rk_com(30)*s(ioh)*s(ih2o2) 


rk_com(30) = ARR(2.9e-12, -160.)


      Function ARR(AA,BB)
      include 'chm.inc' ! for te
      ARR = AA*exp(BB/te)
      return
      end

