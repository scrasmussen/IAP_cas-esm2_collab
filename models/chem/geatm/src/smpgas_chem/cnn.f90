
ppb = 1.e+9

cppb(l) = (cnn(l)/cair_mlc)*ppb


      data pi           /3.141592654/
      data avogad       /6.02217e+23/
      data deg2rad      /0.017453293/


cair_mlc = avogad*pr_atm/(82.056*te)


pr_atm =  PA2ATM*Plev(i03+ixy)*100.  ! the pressure but as the atm




