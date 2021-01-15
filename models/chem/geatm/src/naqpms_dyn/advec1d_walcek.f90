subroutine advec1d_walcek( sx,ex,Q00,QN,U,DEN0,DEN1,DXX,DD0,DT &
                          ,p1d,t1d,rhoa,mwgt,unitlab)

use adv1d_comv

character(len=*) :: unitlab

integer                      :: sx,ex
integer,dimension(sx-1:ex+1) :: IMXMN
real,   dimension(sx-1:ex+1) :: Q0,Q00
real,   dimension(sx:ex)     :: QN,DXX
real,   dimension(sx-1:ex)   :: U,DD0,FLUX
real,   dimension(sx:ex)     :: DEN0,DEN1,VCMAX,VCMIN

real                         :: mwgt
real,   dimension(sx-1:ex+1) :: p1d,t1d,rhoa,trfac

real,   dimension(sx-1:ex+1) :: jcbcel

DATA ZR0,LSHP/0.,0./

! identify local max and min, specigy mxing ratio limit at new time

! Q0 is advected in mass mixing ratio (ppbm)

jcbcel=jcbin


do i=sx-1,ex+1
!if(1==2) then
  if(trim(unitlab).eq.'ppbv') then
    trfac(i)=mwgt/( 0.08206*t1d(i)/(p1d(i)/1000.0) )/rhoa(i)
    Q0(i)=Q00(i)*trfac(i)
  elseif(trim(unitlab).eq.'ugom3') then
    trfac(i)=1.0/rhoa(i)
    Q0(i)=Q00(i)*trfac(i) 
  elseif(trim(unitlab).eq.'ppbm') then
    trfac(i)=1.0
    Q0(i)=Q00(i)
  else
    print*,'unitlab : ',trim(unitlab)
    stop 'unitlab err in advec1d_walcek.f90'
  endif
!else
!  trfac(i)=1.0
!  Q0(i)=Q00(i)
!endif
enddo

Q0=Q0*jcbcel


!print*,'subQ0=',Q0



DO 5 I=sx,ex
   IMXMN(i)=0
   if(Q0(I).GE.AMAX1(Q0(I-1),Q0(I+1)).or.  &
      Q0(I).LE.AMIN1(Q0(I-1),Q0(I+1))    ) IMXMN(I)=1
   CK1=Q0(I)
   IF(U(I ).lt.ZR0)  CK1=Q0(I+1)
   IF(U(I-1).ge.ZR0) CK1=Q0(I-1)
   DEN1(I) = amax1(DEN1(I) , 1.E-20)
   DXX(I) = amax1(DXX(I), 1.E-20)
   VCMAX(I)=aMAX1(Q0(I),CK1)
 5 VCMIN(I)=amin1(q0(i),ck1)


! update mixing ratios and limit Flux going up where u>0
 IF(U(sx-1).GE.ZR0) FLUX(sx-1) = Q0(sx-1)*U(sx-1)*DT*DD0(sx-1)
 DO 10  I=sx,ex
   IF(U(I).lt.ZR0) GOTO 10
   IF(U(I-1).lt.ZR0) THEN
     FLUX(I)=Q0(I)*U(I)*DT*DD0(I)  ! outflow only cell
   ELSE
     X1=DT*U(I)/DXX(I)
     CF1=CFCALC22(X1,Q0(I),Q0(I-1),Q0(I+1),IMXMN(I-1),IMXMN(I+1),LSHP)
     FLUXI=DD0(I)*U(I)*DT*CF1/DXX(I)
     QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
           (Q0(I)*DEN0(I) - FLUXI + FLUX(I-1)/DXX(I))/DEN1(I) ))
     QN(I)=AMAX1(QN(I),0.0)
     FLUX(I)= DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I)) + FLUX(I-1)
   ENDIF
 10 CONTINUE


! update mixing ratios and limit Flux going down where u<0

 IF(U(ex).lt.ZR0) FLUX(ex)=Q0(ex+1)*U(ex)*DT*DD0(ex)
 DO 20 I=ex,sx,-1
    IF(U(I-1).GE.ZR0)THEN
      if(U(I).lt.ZR0) then
           QN(I) =    &     ! inflow-only cell
                (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I)
           QN(I)=AMAX1(QN(I),0.0)
      endif
    ELSE
      X1=DT*ABS(U(I-1))/DXX(I)
      CF1=CFCALC22(X1,Q0(I),Q0(I+1),Q0(I-1),IMXMN(I+1),IMXMN(I-1),LSHP)
      IF(U(I).ge.ZR0) CF1=Q0(I)   ! outflow only cell
      FLUXex1=DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
      QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXex1)/DEN1(I)))
      QN(I)=AMAX1(QN(I),0.0)
      FLUX(I-1)= DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)
   ENDIF
 20 CONTINUE

   do i=sx,ex
     QN(i)=QN(i)/jcbcel(i)/trfac(i)
   enddo

   RETURN


end subroutine advec1d_walcek



 FUNCTION CFCALC22(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
  ! function to calculate mixing ratio in the fuild moving into neighnor
  ! cell durint DT; X1=Courant No. set LSHP=1 for full sharpening

  ALFA=1.

  IF(IMXDN1.gt.0) ALFA=1.7763-0.5*X1
  IF(IMXUP2.gt.0) ALFA=1.2578+0.58502*X1
  IF(LSHP==1)     ALFA=5.

  IF(X1.lt.0.5)THEN
     ALFA=MIN(ALFA,.98*6./(1.-2.*X1))
     CF=VCUP1*(1.+ALFA/6.*(2.*X1-1.)) + &
        VCDW1*ALFA/12.*(4-5.*X1) + VCUP2*ALFA/12.*(X1-2.)
  ELSE
     X1=MAX(X1,1.E-20)
     CF=VCUP1*(1.-ALFA*(1./X1 + 2.*X1 -3.)/6.) + &
        VCDW1*ALFA*(1./X1-X1)/12. + VCUP2*ALFA*(1./X1+5.*X1-6.)/12.
  ENDIF

! limit outflow mixing ratio to reasonable mixing ratio

  CFCALC22=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))

  RETURN

 END


