
  subroutine adv_hori( myid, c, u, v, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt)
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, u, v, deltx, delty
  real, dimension(sx-1:ex+1) :: Q0,QN,u1D,DXX,DD0,DEN0,DEN1 
  real, dimension(sy-1:ey+1) :: Q0y,QNy,u1Dy,DXXy,DD0y,DEN0y,DEN1y 
  integer i, j 
  
  data d0/1./

  dt0 = 300.
  do j = sy,ey
  do i = sx,ex
   dt0 = min(dt0, deltx(i,j)/abs(u(i,j)),delty(i,j)/abs(v(i,j)))
  enddo
  enddo

  istep=dt/dt0+1
  dtt=dt/istep
  
do 10 itt=1,istep
 
  do j = sy,ey        ! do I advection
      do i = sx-1,ex+1
           Q0(i)=c(i,j)
           u1D(i)=u(i,j) 
           DXX(i)=deltx(i,j)
           DEN0(i)=d0
!           DEN1(i)=1.-dtt/DXX(i)*(d0*u(i,j)-d0*u(i-1,j)) ! juanxiong he
           DD0(i)=DEN0(i)
      enddo
      do i = sx,ex+1
           DEN1(i)=1.-dtt/DXX(i)*(d0*u(i,j)-d0*u(i-1,j)) !juanxiong he
      enddo
      call advec1d(sx,ex,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,dtt)
      do i = sx,ex
          c(i,j)=amax1(QN(i),1.E-20)
      enddo
  enddo
   do i = sx,ex   ! do j direction
        do j = sy-1,ey+1
           Q0y(j)=c(i,j)
           u1Dy(j)=v(i,j)
           DXXy(j)=delty(i,j)
           DEN0y(j)=1.-dtt/deltx(i,j)*(d0*u(i,j)-d0*u(i-1,j))
           DD0y(j)=d0
         enddo
      if(sy==1)then
           DEN1y(0)=DEN0y(0)
        do j = sy,ey+1
           DEN1y(j)=DEN0y(j)-dtt/delty(i,j)*(d0*v(i,j)-d0*v(i,j-1))
        enddo
       else
      !  do j = sy-1,ey+1 
        do j = sy,ey+1 ! juanxiong he
           DEN1y(j)=DEN0y(j)-dtt/delty(i,j)*(d0*v(i,j)-d0*v(i,j-1))
        enddo
       endif
      call advec1d(sy,ey,Q0y,QNy,u1Dy,DEN0y,DEN1y,DXXy,DD0y,dtt)
      do j = sy,ey
          c(i,j)=amax1(QNy(j),1.0E-20)
      enddo
  enddo

 10 continue
 return
  end


subroutine advec1d(sx,ex,Q0,QN,U,DEN0,DEN1,DXX,DD0,DT)
integer sx,ex
real, dimension(sx-1:ex+1) :: Q0,QN,U,DEN0,DEN1,DXX,DD0
real, dimension(sx-1:ex+1) :: FLUX,VCMAX,VCMIN
integer,dimension(sx-1:ex+1) :: IMXMN
DATA ZR0,LSHP/0.,0./

! identify local max and min, specigy mxing ratio limit at new time

DO 5 I=sx,ex
   IMXMN(i)=0
   if(Q0(I).GE.AMAX1(Q0(I-1),Q0(I+1)).or.  &
      Q0(I).LE.AMIN1(Q0(I-1),Q0(I+1))    ) IMXMN(I)=1
   CK1=Q0(I)
   IF(U(I ).lt.ZR0)CK1=Q0(I+1)
   IF(U(I-1).ge.ZR0)CK1=Q0(I-1)
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
        CF1=CFCALC(X1,Q0(I),Q0(I-1),Q0(I+1),IMXMN(I-1),IMXMN(I+1),LSHP)
        FLUXI=DD0(I)*U(I)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I) - FLUXI + FLUX(I-1)/DXX(I))/DEN1(I) ))
        FLUX(I)= DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I)) + FLUX(I-1)
   ENDIF   
 10 CONTINUE 

! update mixing ratios and limit Flux going down where u<0

  IF(U(ex).lt.ZR0)FLUX(ex)=Q0(ex+1)*U(ex)*DT*DD0(ex)
  DO 20 I=ex,sx,-1
  IF(U(I-1).GE.ZR0)THEN   
         if(U(I).lt.ZR0) QN(I) =    &     ! inflow-only cell
         (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I)  
   ELSE
     X1=DT*ABS(U(I-1))/DXX(I)
     CF1=CFCALC(X1,Q0(I),Q0(I+1),Q0(I-1),IMXMN(I+1),IMXMN(I-1),LSHP)
     IF(U(I).ge.ZR0)CF1=Q0(I)   ! outflow only cell
        FLUXex1=DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXex1)/DEN1(I)))
        FLUX(I-1)= DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)
   ENDIF 
 20 continue
   RETURN
   END
     
  FUNCTION CFCALC(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
  ! function to calculate mixing ratio in the fuild moving into neighnor
  ! cell durint DT; X1=Courant No. set LSHP=1 for full sharpening
  ALFA=1.
  IF(IMXDN1.gt.0)ALFA=1.7763-0.5*X1
  IF(IMXUP2.gt.0)ALFA=1.2578+0.58502*X1
  IF(LSHP==1)ALFA=5.
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
    CFCALC=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))
    RETURN
    END 
  

 subroutine  GetMassRatio(myid,c1,RatioMass,sx, ex, sy, ey ,nx,ny)
  integer myid, sx, ex, sy, ey,nx,ny
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1,RatioMass
  integer i, j

  do j = sy,ey
  do i = sx,ex
     RatioMass(i,j) = (c1(i,j)-1000.)/1000.
     RatioMass(i,j) = amin1(0.5,RatioMass(i,j))
     RatioMass(i,j) = amax1(-0.5,RatioMass(i,j))
  enddo
  enddo

  return
  end

  subroutine  balance(myid,c,c2,ratio,sx, ex, sy, ey ,nx,ny)
  integer myid, sx, ex, sy, ey,nx,ny
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,ratio,c2
  integer i, j
  do j = sy,ey
  do i = sx,ex
     c(i,j)=c(i,j)-c2(i,j)*ratio(i,j)
     c(i,j)=amax1(c(i,j),1.E-20)
  enddo
  enddo

  return
  end


