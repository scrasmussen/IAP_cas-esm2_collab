
!subroutine carbon_bin_info

implicit none
integer,parameter :: nbincb=28,nmode=2
real*8 ,parameter :: radmed=3.0d-8   !m
real*8 ,parameter :: sigmacb=1.8
real*8 ,parameter :: third=1.0d0/3.0d0
real*8 ,parameter :: cbrad1= 3.0d-8,cbrad2=7.5d-8
real*8 ,parameter :: cbsigma1=1.8,cbsigma2=1.8
real*8 ,parameter :: dencb=1800 ! kg/m3
real*8 ,parameter :: gen2pi=1.0/sqrt(2.d0*3.1415926d0)
real*8  :: radcbm(nbincb),vcbm3(nbincb),coagparcb(nbincb,nbincb,nbincb)
real*8  :: radedge(nbincb+1),deltdp(nbincb)
real*8  :: nfrac1(nbincb),nfrac2(nbincb)
real*8  :: cbemt2bin(nbincb,nmode)

real*8  :: vrat,rmin,vij,narb,nleft,nrght,ntmp,xx1,xx2,xxx
real*8  :: ttfrac1,ttfrac2,ttmass1,ttmass2

real*8  :: xtemp

integer :: ibin,ii,jj,kk,nn

rmin = 5.0d-9

vrat = 2.0d0

narb=1000.0d0

radcbm(1)=rmin
do ibin=2,nbincb
 radcbm(ibin)=radcbm(ibin-1)*vrat**third
enddo

radedge(1)=rmin/sqrt(vrat**third)
do ibin=1,nbincb
radedge(ibin+1)=radcbm(ibin)*sqrt(vrat**third)
deltdp(ibin)=radedge(ibin+1)-radedge(ibin)
enddo


vcbm3=4.0d0/3.0d0*3.1416*radcbm**3.0d0


do ibin=1,nbincb
 write(*,'(i3,3(2x,f10.7))') ibin,radcbm(ibin)*1.0e6,radedge(ibin)*1.0e6,radedge(ibin+1)*1.0e6
enddo


ttfrac1=0.0
ttfrac2=0.0
ttmass1=0.0
ttmass2=0.0
do ibin=1,nbincb
xx1=log(radcbm(ibin)/cbrad1)**2.0/(2.0*log(cbsigma1)**2.0)
xx2=log(radcbm(ibin)/cbrad2)**2.0/(2.0*log(cbsigma2)**2.0)
nfrac1(ibin)=deltdp(ibin)/radcbm(ibin)*exp(-xx1)/(sqrt(2.0*3.1415926)*log(cbsigma1))
nfrac2(ibin)=deltdp(ibin)/radcbm(ibin)*exp(-xx2)/(sqrt(2.0*3.1415926)*log(cbsigma2))
ttfrac1=ttfrac1+nfrac1(ibin)
ttfrac2=ttfrac2+nfrac2(ibin)
enddo

nfrac1=nfrac1/ttfrac1
nfrac2=nfrac2/ttfrac2

print*,sum(nfrac1),sum(nfrac2)

do ibin=1,nbincb

  xtemp=narb*nfrac1(ibin)/ttfrac1
  cbemt2bin(ibin,1)=xtemp
  ttmass1=ttmass1+xtemp*dencb*vcbm3(ibin)

  xtemp=narb*nfrac2(ibin)/ttfrac2
  cbemt2bin(ibin,2)=xtemp 
  ttmass2=ttmass2+xtemp*dencb*vcbm3(ibin)

enddo

do ibin=1,nbincb
  cbemt2bin(ibin,1)=cbemt2bin(ibin,1)/ttmass1
  cbemt2bin(ibin,2)=cbemt2bin(ibin,2)/ttmass2
enddo

do ibin=1,nbincb
  write(*,'(i2,2x,f12.8,2(2x,e12.6))') ibin,radcbm(ibin)*1.0e6,cbemt2bin(ibin,1),cbemt2bin(ibin,2)
enddo


nn=0

do ii=1,nbincb
  do jj=1,ii
    vij=vcbm3(ii)+vcbm3(jj)
    if(vij.ge.vcbm3(nbincb)) then
      coagparcb(ii,jj,nbincb)=1.0d0
      coagparcb(jj,ii,nbincb)=1.0d0
      print*,'supmax'
    else
      if(ii.le.nbincb-1) then
         do kk=ii,nbincb-1
           if(vij.ge.vcbm3(kk).and.vij.lt.vcbm3(kk+1)) then
             coagparcb(ii,jj,kk)=(vcbm3(kk+1)-vij)/(vcbm3(kk+1)-vcbm3(kk))* &
                                 vcbm3(kk)/vij
             coagparcb(ii,jj,kk+1)=1.0d0-coagparcb(ii,jj,kk)
             coagparcb(jj,ii,kk)=coagparcb(ii,jj,kk)
             coagparcb(jj,ii,kk+1)=coagparcb(ii,jj,kk+1)
             nn=nn+1
             print863,nn,ii,jj,kk,kk+1,coagparcb(ii,jj,kk),coagparcb(ii,jj,kk+1)
             if(nn.eq.2) then
               !print*,vcbm3(kk),vcbm3(kk+1),vij
             endif
           endif
         enddo
      endif
    endif
  enddo
enddo

863 format(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,f5.3,2x,f5.3)




end

!end subroutine carbon_bin_info
