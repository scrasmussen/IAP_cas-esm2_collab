      subroutine kv_ob70(nz,kvmin,pbl,zz,tv,uu,vv,rkv)
!
!-----Determines Kv profile using the profile methodology of O'Brien (1970). 
!
      implicit none
!
      integer nz,k
      real kvmin,pbl,zz(nz),tv(nz),uu(nz),vv(nz),rkv(nz)
      real dkdzs,obk,rktop,vonk,xinf,g
      real dzm,dtdz,tbar,dudz,dvdz,dwdz,ri,sbl,ksbl,xl,sh,cc1,cc2,ccc
      data vonk,xinf,g,rktop /0.4,100.,9.8,1.0/
!
!-----Determine Kv at top of surface layer based on Louis (1979) relationships
!
      dzm  = zz(2)/2.
      dtdz = (tv(2) - tv(1))/dzm
      tbar = (tv(2) + tv(1))/2.
      dudz = (uu(2) - uu(1))/dzm
      dvdz = (vv(2) - vv(1))/dzm
      dwdz = dudz*dudz + dvdz*dvdz
      ri = (g/tbar)*dtdz/(dwdz + 1.e-20)

      pbl = max(pbl,zz(1))
      sbl = max(zz(1),min(100.,vonk*pbl/10.))
      xl = vonk*sbl/(1. + vonk*sbl/xinf)

      if (ri.ge.0.) then
        sh = 1./(1. + 4.7*ri)**2
      else
        cc1 = (zz(1) + dzm/zz(1))**0.333
        cc1 = (cc1 - 1.)**1.5
        cc2 = 1./sqrt(zz(1)*dzm**3)
        ccc = 49.8*xl*xl*cc1*cc2
        sh = 1. - 9.4*ri/(1. + ccc*sqrt(abs(ri)))
      endif
      rkv(1) = max(kvmin,sh*xl*xl*sqrt(dwdz))
      dkdzs  = rkv(1)/zz(1)
      ksbl   = sbl*dkdzs
!
!-----Determine Kv through PBL using profile method
!
      do k = 2,nz-1
        if (zz(k).lt.sbl) then
          rkv(k) = max(zz(k)*dkdzs,kvmin)
        elseif (zz(k).lt.pbl) then
          obk = rktop + ((pbl - zz(k))**2/(pbl - sbl)**2)* &
                         (ksbl - rktop + (zz(k) - sbl)*    &
                         (dkdzs + (2.*(ksbl - rktop)/    &
                         (pbl - sbl))))
          rkv(k) = max(obk,kvmin)
        else
          rkv(k) = kvmin
        endif
      enddo
      rkv(nz) = 0. 
!
      return
      end
