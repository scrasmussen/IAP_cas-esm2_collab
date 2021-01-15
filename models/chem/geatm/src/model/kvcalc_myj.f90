      subroutine kv_tkemyj(nz,tv,uu,vv,zh,zf,tke,el,kvmin,pblc,kv)
!
!-----Determines Kv at a given layer interface using the TKE methodology  
!     employed in the MYJ (Mellor-Yamada-Janjic) scheme in WRF
!     Assumes that Kv follows diffusivity for heat (Kh)

!     TKE and mixing lengths are read directly from WRF outputs 

      implicit none
!
!     Inputs
      integer nz
      real tv(nz),uu(nz),vv(nz),zh(nz),zf(nz),tke(nz),el(nz),kvmin,pblc

!     Output
      real kv(nz)

!     Other variables
      integer k
      real q(nz)
      real g,dz,dtemp,du,dv,duv,dth

      real a1,a2x,b1,b2,c1
      real btg,aeqh,aeqm,adnm,adnh,ubry3
      real bdnm,bdnh,bubm,bubh
      real bshm,bshh,besh
      real eloq2,sh
      real gm(nz),gh(nz)
      real aden,bden,cden,rden
      real grav,beta,vk


      data a1 /0.659888514560862645/
      data a2x /0.6574209922667784586/
      data b1 /11.87799326209552761/
      data b2 /7.226971804046074028/
      data c1 /0.000830955950095854396/

      data  grav/9.81/, beta/273./
      data  vk/0.4/

!
!     diffusivities are computed as kv = q * el * f(stability)
!     where q = sqrt(2*tke) and el = mixing lengths
!
!
!---  Compute some coefficients
!
      btg = grav/beta
      aeqh = 9.*a1*a2x*a2x*b1*btg*btg+  &
             9.*a1*a2x*a2x*(12.*a1+3.*b2)*btg*btg
      aeqm = 3.*a1*a2x*b1*(3.*a2x+3.*b2*c1+18.*a1*c1-b2)*btg &
             +18.*a1*a1*a2x*(b2-3.*a2x)*btg
      adnm = 18.*a1*a1*a2x*(b2-3.*a2x)*btg
      adnh =  9.*a1*a2x*a2x*(12.*a1+3*b2)*btg*btg

      ubry3 = 3.*(1.+1.e-7)*(18.*(aeqh/aeqm)*a1*a1*a2x*b2*c1*btg &
              + 9.*a1*a2x*a2x*b2*btg*btg)/((aeqh/aeqm)*adnm+adnh)

      bdnm = 6.*a1*a1
      bdnh = 3.*a2x*(7*a1+b2)*btg

!
!---  Solve for q
!
      do k=1,nz
        q(k) = sqrt(2.0 * tke(k))
      enddo

!
!-----Solve the stability determinant, sh
!
      do k=1,nz-1
        dz = (zh(k+1) - zh(k))
        dth = tv(k+1) - tv(k)
        du = (uu(k+1) - uu(k))
        dv = (vv(k+1) - vv(k))
        duv = (du*du + dv*dv)/(dz*dz)

        gm(k) = max(duv,aeqh/aeqm*1.e-9)
        gh(k) = dth/dz
        if (abs(gh(k)) .le. 1.e-9) gh(k) = 1.e-9

        eloq2 = el(k)*el(k)/(q(k)**2.)
        bshm = 18.*a1*a1*a2x*c1
        bshh =  9.*a1*a2x*a2x*btg
        besh = bshm*gm(k) + bshh*gh(k)

        aden = (adnm*gm(k) + adnh*gh(k)) * gh(k)
        bden = bdnm * gm(k) + bdnh* gh(k)
        cden = 1.0
!       Place a cap on the denominator
        rden = 1.0/(max(0.1,aden*eloq2*eloq2 + bden*eloq2 + cden))

        sh = (besh*eloq2 + a2x) * rden

!
!----   Compute vertical diffusivities
!
        kv(k) = q(k) * el(k) * sh
        kv(k) = amax1(kv(k),kvmin)
      enddo
      kv(nz) = 0.
      
      return
      end

