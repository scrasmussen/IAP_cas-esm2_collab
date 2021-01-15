           subroutine calfconv(fconv,xmol,pbl)
             real :: fconv,xmol,pbl 
             const1=(0.4**(-2./3.))/(0.1*7.2) 
             const2=(-pbl/(xmol+1.e-6))
             if(const2>0.) then
              const3=const2**(-1./3.)
             else
              const3=-abs(const2)**(-1./3.)
             endif
             if(xmol<0.0) then
              fconv=1./(1.+const1*const3)
             else
              fconv=0.0
             endif

           return
           end
