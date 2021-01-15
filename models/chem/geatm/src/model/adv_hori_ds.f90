
  subroutine adv_hori_ds( myid, c, dcom1,dcom2,dcom3,dcom4,dcom5,dcom6,dcom7,dcom8,&
                          dcom9,dcom10,scom1,scom2,scom3,scom4,scom5,scom6,scom7,scom8,&
                           ndcom,nscom, u, v, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt)
! cccc THIS ROUTINE IS TO CALCULATE THE ADVECTION OF DUST AND SEA SALT AEROSOLS AND AEROSOLS ccccc
!       COMPOSITIONS (ASO4, ANO3,FEII,FEIII) ASSUMING DUST (OR. SEA  SALT) INTERNAL MIXTURE cccccc

  integer myid, sx, ex, sy, ey, it,ndcom,nscom ! compositions numbers
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c, u, v, deltx, delty
  real, dimension(sx-1:ex+1) :: Q0,QN,u1D,DXX,DD0,DEN0,DEN1 
  real, dimension(sy-1:ey+1) :: Q0y,QNy,u1Dy,DXXy,DD0y,DEN0y,DEN1y 
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: dcom1,dcom2,dcom3,dcom4,dcom5,dcom6,dcom7,dcom8,dcom9,dcom10 !  dust  particles compositions
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: scom1,scom2,scom3,scom4,scom5,scom6,scom7,scom8 ! sea salt compositon
  real,dimension(sx-1:ex+1,ndcom) :: Dcom,Dcomy
  real,dimension(sx-1:ex+1,nscom) :: Scom,Scomy
!  real,dimension(sy-1:ey+1) :: DSSO4y, DSNO3y, DSFEIIy, DSFEIIIy
 
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
           DEN1(i)=1.-dtt/DXX(i)*(d0*u(i,j)-d0*u(i-1,j))
           DD0(i)=DEN0(i)
           Dcom(i,1) = dcom1(i,j)
           Dcom(i,2) = dcom2(i,j)
           Dcom(i,3) = dcom3(i,j)
           Dcom(i,4) = dcom4(i,j)
           Dcom(i,5) = dcom5(i,j)
           Dcom(i,6) = dcom6(i,j)
           Dcom(i,7) = dcom7(i,j)
           Dcom(i,8) = dcom8(i,j)
           Dcom(i,9) = dcom9(i,j)
           Dcom(i,10) =dcom10(i,j)       
           Scom(i,1) = scom1(i,j)
           Scom(i,2) = scom2(i,j)
           Scom(i,3) = scom3(i,j)
           Scom(i,4) = scom4(i,j)
           Scom(i,5) = scom5(i,j)
           Scom(i,6) = scom6(i,j)
           Scom(i,7) = scom7(i,j)
           Scom(i,8) = scom8(i,j)
      enddo
      call advec1d_ds(sx,ex,Q0,QN,Dcom,Scom,ndcom,nscom,u1D,DEN0,DEN1,DXX,DD0,dtt)
      do i = sx,ex
          c(i,j)=amax1(QN(i),1.E-20)
         
            dcom1(i,j) = Dcom(i,1)
            dcom2(i,j) = Dcom(i,2)
            dcom3(i,j) = Dcom(i,3)
            dcom4(i,j) = Dcom(i,4)
            dcom5(i,j) = Dcom(i,5)
            dcom6(i,j) = Dcom(i,6)
            dcom7(i,j) = Dcom(i,7)
            dcom8(i,j) = Dcom(i,8)
            dcom9(i,j) = Dcom(i,9)
            dcom10(i,j)= Dcom(i,10)
            scom1(i,j) = Scom(i,1)
            scom2(i,j) = Scom(i,2)
            scom3(i,j) = Scom(i,3)
            scom4(i,j) = Scom(i,4)
            scom5(i,j) = Scom(i,5)
            scom6(i,j) = Scom(i,6)
            scom7(i,j) = Scom(i,7)
            scom8(i,j) = Scom(i,8)
      enddo
 
  enddo

   do i = sx,ex   ! do j direction
        do j = sy-1,ey+1
           Q0y(j)=c(i,j)
           u1Dy(j)=v(i,j)
           DXXy(j)=delty(i,j)
           DEN0y(j)=1.-dtt/deltx(i,j)*(d0*u(i,j)-d0*u(i-1,j))
           DD0y(j)=d0
           Dcomy(j,1) = dcom1(i,j)
           Dcomy(j,2) = dcom2(i,j)
           Dcomy(j,3) = dcom3(i,j)
           Dcomy(j,4) = dcom4(i,j)
           Dcomy(j,5) = dcom5(i,j)
           Dcomy(j,6) = dcom6(i,j)
           Dcomy(j,7) = dcom7(i,j)
           Dcomy(j,8) = dcom8(i,j)
           Dcomy(j,9) = dcom9(i,j)
           Dcomy(j,10) =dcom10(i,j)
           Scomy(j,1) = scom1(i,j)
           Scomy(j,2) = scom2(i,j)
           Scomy(j,3) = scom3(i,j)
           Scomy(j,4) = scom4(i,j)
           Scomy(j,5) = scom5(i,j)
           Scomy(j,6) = scom6(i,j)
           Scomy(j,7) = scom7(i,j)
           Scomy(j,8) = scom8(i,j)
         enddo
      if(sy==1)then
           DEN1y(0)=DEN0y(0)
        do j = sy,ey+1
           DEN1y(j)=DEN0y(j)-dtt/delty(i,j)*(d0*v(i,j)-d0*v(i,j-1))
        enddo
       else
        do j = sy-1,ey+1
           DEN1y(j)=DEN0y(j)-dtt/delty(i,j)*(d0*v(i,j)-d0*v(i,j-1))
        enddo
       endif
      call advec1d(sy,ey,Q0y,QNy,Dcomy,Scomy,ndcom,nscom,u1Dy,DEN0y,DEN1y,DXXy,DD0y,dtt)

      do j = sy,ey
          c(i,j)=amax1(QNy(j),1.E-20)

            dcom1(i,j) = Dcomy(j,1)
            dcom2(i,j) = Dcomy(j,2)
            dcom3(i,j) = Dcomy(j,3)
            dcom4(i,j) = Dcomy(j,4)
            dcom5(i,j) = Dcomy(j,5)
            dcom6(i,j) = Dcomy(j,6)
            dcom7(i,j) = Dcomy(j,7)
            dcom8(i,j) = Dcomy(j,8)
            dcom9(i,j) = Dcomy(j,9)
            dcom10(i,j)= Dcomy(j,10)
            scom1(i,j) = Scomy(j,1)
            scom2(i,j) = Scomy(j,2)
            scom3(i,j) = Scomy(j,3)
            scom4(i,j) = Scomy(j,4)
            scom5(i,j) = Scomy(j,5)
            scom6(i,j) = Scomy(j,6)
            scom7(i,j) = Scomy(j,7)
            scom8(i,j) = Scomy(j,8)

      enddo

  enddo
   


 10 continue
 return
  end


subroutine advec1d_ds(sx,ex,Q0,QN,Dcom,Scom,ndcom ,nscom,U,DEN0,DEN1,DXX,DD0,DT)
integer sx,ex,ndcom,nscom
real, dimension(sx-1:ex+1) :: Q0,QN,U,DEN0,DEN1,DXX,DD0
real, dimension(sx-1:ex+1,ndcom) :: Dcom
real, dimension(sx-1:ex+1,nscom) :: Scom
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
        CF1=CFCALC_dust(X1,Q0(I),Q0(I-1),Q0(I+1),IMXMN(I-1),IMXMN(I+1),LSHP)
        FLUXI=DD0(I)*U(I)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I) - FLUXI + FLUX(I-1)/DXX(I))/DEN1(I) ))
        FLUX(I)= DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I)) + FLUX(I-1)

        IF(FLUX(I).GT. 0. ) THEN

        DO IDUC = 1, NDCOM
         Q0(I) = AMAX1(Q0(I),1.E-20)
         DCOM(I,IDUC) = DCOM(I,IDUC) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DCOM(I-1,IDUC)/Q0(I-1) -&
                   FLUX(I)/DXX(I)/DEN1(I) * DCOM(I,IDUC)/Q0(I)  )
        ENDDO

        DO IDUC = 1, NSCOM
         Q0(I) = AMAX1(Q0(I),1.E-20)
         SCOM(I,IDUC) = SCOM(I,IDUC) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * SCOM(I-1,IDUC)/Q0(I-1) -&
                   FLUX(I)/DXX(I)/DEN1(I) * SCOM(I,IDUC)/Q0(I)  )
        ENDDO


        ELSE IF( FLUX(I).LT. 0. ) THEN
        DO IDUC = 1, NDCOM
         Q0(I) = AMAX1(Q0(I),1.E-20)
        DCOM(I,IDUC) = DCOM(I,IDUC) + &
                  ( FLUX(I-1)/DXX(I)/DEN1(I) * DCOM(I-1,IDUC)/Q0(I-1) -&
                   FLUX(I)/DXX(I)/DEN1(I) * DCOM(I+1,IDUC)/Q0(I+1))
        ENDDO
        DO IDUC = 1, NSCOM
         Q0(I) = AMAX1(Q0(I),1.E-20)
         SCOM(I,IDUC) = SCOM(I,IDUC) + &
                  ( FLUX(I-1)/DXX(I)/DEN1(I) * SCOM(I-1,IDUC)/Q0(I-1) -&
                   FLUX(I)/DXX(I)/DEN1(I) * SCOM(I+1,IDUC)/Q0(I+1))
        ENDDO

        ENDIF ! FLUX               


        DO IDUC = 1, NDCOM
         DCOM(I,IDUC) = AMAX1 (  DCOM(I,IDUC), 1.E-20)
        ENDDO 
        DO IDUC = 1, NSCOM
         SCOM(I,IDUC) = AMAX1 (  SCOM(I,IDUC), 1.E-20)
        ENDDO

   ENDIF   
 10 CONTINUE 

! update mixing ratios and limit Flux going down where u<0

  IF(U(ex).lt.ZR0)FLUX(ex)=Q0(ex+1)*U(ex)*DT*DD0(ex)
  DO 20 I=ex,sx,-1
  IF(U(I-1).GE.ZR0)THEN   
         if(U(I).lt.ZR0) then
              QN(I) =    &     ! inflow-only cell
          (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I) 
        
           IF(FLUX(I).GT. 0. ) THEN           
         
            DO IDUC = 1, NDCOM
             Q0(I) = AMAX1(Q0(I),1.E-20)

            DCOM(I,IDUC) = DCOM(I,IDUC) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DCOM(I-1,IDUC)/Q0(I-1) - &
                     FLUX(I)/DXX(I)/DEN1(I) * DCOM(I,IDUC)/Q0(I))
            ENDDO

            DO IDUC = 1, NSCOM
             Q0(I) = AMAX1(Q0(I),1.E-20)

             SCOM(I,IDUC) = SCOM(I,IDUC) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * SCOM(I-1,IDUC)/Q0(I-1) - &
                     FLUX(I)/DXX(I)/DEN1(I) * SCOM(I,IDUC)/Q0(I))
            ENDDO

           ELSE IF( FLUX(I).LT. 0. ) THEN

            DO IDUC = 1, NDCOM
              Q0(I) = AMAX1(Q0(I),1.E-20)
            DCOM(I,IDUC) = DCOM(I,IDUC) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * DCOM(I-1,IDUC)/Q0(I-1) -&
                     FLUX(I)/DXX(I)/DEN1(I) * DCOM(I+1,IDUC)/Q0(I+1))
            ENDDO
            DO IDUC = 1, NSCOM
              Q0(I) = AMAX1(Q0(I),1.E-20)
              SCOM(I,IDUC) = SCOM(I,IDUC) + &
                    ( FLUX(I-1)/DXX(I)/DEN1(I) * SCOM(I-1,IDUC)/Q0(I-1) -&
                     FLUX(I)/DXX(I)/DEN1(I) * SCOM(I+1,IDUC)/Q0(I+1))
            ENDDO

        ENDIF ! FLUX   
       
        DO IDUC = 1, NDCOM
         DCOM(I,IDUC) = AMAX1 (  DCOM(I,IDUC), 1.E-20)
        ENDDO
       DO IDUC = 1, NSCOM
         SCOM(I,IDUC) = AMAX1 (  SCOM(I,IDUC), 1.E-20)
        ENDDO

         endif 
 
   ELSE
     X1=DT*ABS(U(I-1))/DXX(I)
     CF1=CFCALC_dust(X1,Q0(I),Q0(I+1),Q0(I-1),IMXMN(I+1),IMXMN(I-1),LSHP)
     IF(U(I).ge.ZR0)CF1=Q0(I)   ! outflow only cell
        FLUXex1=DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXex1)/DEN1(I)))
        FLUX(I-1)= DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)

      IF( FLUX(I).LT. 0. ) THEN        
        DO IDUC = 1,NDCOM
           Q0(I) = AMAX1(Q0(I),1.E-20)
            DCOM(I,IDUC) = DCOM(I,IDUC) -  &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DCOM(I,IDUC)/Q0(I) -&
                     FLUX(I)/DXX(I)/DEN1(I) * DCOM(I+1,IDUC)/Q0(I+1))
        ENDDO
        DO IDUC = 1,NSCOM
           Q0(I) = AMAX1(Q0(I),1.E-20)
           SCOM(I,IDUC) = SCOM(I,IDUC) -  &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * SCOM(I,IDUC)/Q0(I) -&
                     FLUX(I)/DXX(I)/DEN1(I) * SCOM(I+1,IDUC)/Q0(I+1))
        ENDDO

       ELSE IF (  FLUX(I).GT. 0. ) THEN
         DO IDUC = 1, NDCOM
            Q0(I) = AMAX1(Q0(I),1.E-20)
            DCOM(I,IDUC) = DCOM(I,IDUC) -  &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * DCOM(I,IDUC)/Q0(I) -&
                     FLUX(I)/DXX(I)/DEN1(I) * DCOM(I,IDUC)/Q0(I))
         ENDDO

         DO IDUC = 1, NSCOM
            Q0(I) = AMAX1(Q0(I),1.E-20)
            SCOM(I,IDUC) = SCOM(I,IDUC) -  &
                    ( ABS(FLUX(I-1)/DXX(I)/DEN1(I)) * SCOM(I,IDUC)/Q0(I) -&
                     FLUX(I)/DXX(I)/DEN1(I) * SCOM(I,IDUC)/Q0(I))
         ENDDO

       ENDIF      ! FLUX

        DO IDUC = 1, NDCOM
         DCOM(I,IDUC) = AMAX1 (  DCOM(I,IDUC), 1.E-20)
        ENDDO
        DO IDUC = 1, NSCOM
         SCOM(I,IDUC) = AMAX1 (  SCOM(I,IDUC), 1.E-20)
        ENDDO


   ENDIF 
 20 continue
   RETURN
   END
     
  FUNCTION CFCALC_dust(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
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
    CFCALC_dust=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))
    RETURN
    END 
  



