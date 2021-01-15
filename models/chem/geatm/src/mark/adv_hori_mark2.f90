
  subroutine adv_hori_mark( myid, c, u, v, deltx, delty, &
                       sx, ex, sy, ey ,ne,nx,ny,dt,k,ktop,&
                       ISMMAX,SM)
  integer myid, sx, ex, sy, ey, it,nx,ISMMAX
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c,c1, u, v, deltx, delty
  real, dimension(sx-1:ex+1) :: Q0,QN,u1D,DXX,DD0,DEN0,DEN1 
  real, dimension(sy-1:ey+1) :: Q0y,QNy,u1Dy,DXXy,DD0y,DEN0y,DEN1y 
  integer i, j, k
  real,dimension(sx-1:ex+1,sy-1:ey+1) :: ktop
  REAL,ALLOCATABLE,DIMENSION(:,:)     ::  TMPsmconv
  REAL,DIMENSION(ISMMAX,SX-1:EX+1,SY-1:EY+1) :: SM
  INTEGER                   :: ITYPE
  
  data d0/1./

  dt0 = 300.
  do j = sy,ey
  do i = sx,ex
   dt0 = min(dt0, deltx(i,j)/abs(u(i,j)),delty(i,j)/abs(v(i,j)))
  enddo
  enddo

  istep=dt/dt0+1
  dtt=dt/istep

  do j=sy-1,ey+1
   do i=sx-1,ex+1
    c1(i,j)= c(i,j)
   enddo
  enddo   
  
  
do 10 itt=1,istep
 
  do j = sy,ey        ! do I advection
      do i = sx-1,ex+1
           Q0(i)=c1(i,j)
           u1D(i)=u(i,j) 
           DXX(i)=deltx(i,j)
           DEN0(i)=d0
           DEN1(i)=1.-dtt/DXX(i)*(d0*u(i,j)-d0*u(i-1,j))
           DD0(i)=DEN0(i)
      enddo
       KKKTOP   =  KTOP(I,J)

       ITYPE = 1

      ALLOCATE(TMPsmconv(ISMMAX,SX-1:EX+1)) 
      DO  is=1,ismMax
       DO I=SX-1,EX+1
        TMPsmconv(is,I) = SM(IS,I,J)
!         IF(I==84.AND.J==43.AND.K==1.AND.IS==3) PRINT*, SM(3,I,J) , 'inner1'
       ENDDO 
      ENDDO 
      
      
      call advec1d_mark2(sx,ex,Q0,QN,u1D,DEN0,DEN1,DXX,DD0,dtt,ISMMAX,&
                 TMPsmconv,KKKTOP,NE,NX,NY,ITYPE,K)

      DO  is=1,ismMax
       DO I=SX-1,EX+1
          SM(IS,I,J) = TMPsmconv(is,I)
!          IF(I==84.AND.J==43.AND.K==1.AND.IS==3) PRINT*, SM(3,I,J) , 'inner2'
       ENDDO   
      ENDDO 
                 
      do i = sx,ex
          c1(i,j)=amax1(QN(i),0.)
      enddo
      DEALLOCATE(TMPsmconv)
  enddo
  
   do i = sx,ex   ! do j direction
        do j = sy-1,ey+1
           Q0y(j)=c1(i,j)
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
        do j = sy-1,ey+1
           DEN1y(j)=DEN0y(j)-dtt/delty(i,j)*(d0*v(i,j)-d0*v(i,j-1))
        enddo
       endif

       ITYPE = 2
       ALLOCATE(TMPsmconv(ISMMAX,SY-1:EY+1))
         
      DO  is=1,ismMax
        DO J=SY-1,EY+1
          TMPsmconv(is,J) = SM(IS,I,J)
        ENDDO  
      ENDDO  
      call advec1d_mark2(sy,ey,Q0y,QNy,u1Dy,DEN0y,DEN1y,DXXy,DD0y,dtt,ISMMAX,&
              TMPsmconv,KKKTOP,NE,NX,NY,ITYPE,K)
      
      DO  is=1,ismMax        
       DO J=SY-1,EY+1
        SM(IS,I,J) = TMPsmconv(is,J)
       ENDDO 
      ENDDO 

      DEALLOCATE(TMPsmconv)
      do j = sy,ey
          c1(i,j)=amax1(QNy(j),0.)
      enddo
  enddo

 10 continue
 return
  end


subroutine advec1d_mark2(sx,ex,Q0,QN,U,DEN0,DEN1,DXX,DD0,DT,ISMMAX,&
              TMPsmconv,KKKTOP,NE,NX,NY,ITYPE,K)
integer sx,ex
real, dimension(sx-1:ex+1)   :: Q0,QN,U,DEN0,DEN1,DXX,DD0
real, dimension(sx-1:ex+1)   :: FLUX,VCMAX,VCMIN
integer,dimension(sx-1:ex+1) :: IMXMN
INTEGER                      :: ITYPE,IMONTHERBOUNDATY
INTEGER                      :: ISMMAX,KKKTOP,NE,NX,NY
REAL,DIMENSION(ISMMAX,SX-1:EX+1) :: TMPsmconv
REAL,DIMENSION(ISMMAX):: smthis,smother
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
        CF1=CFCALC_mark2(X1,Q0(I),Q0(I-1),Q0(I+1),IMXMN(I-1),IMXMN(I+1),LSHP)
        FLUXI=DD0(I)*U(I)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I) - FLUXI + FLUX(I-1)/DXX(I))/DEN1(I) ))
!!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF(FLUX(I-1)>=0.0) THEN
          DO is=1,ismMax
            smthis(is) = TMPsmconv(is,I)
            smother(is)= TMPsmconv(is,I-1)
!            IF(I==84.AND.IJ==43.AND.K==1) PRINT*,IS,smthis(IS),'inner inner 01'
          ENDDO

                    
!         deltc1 = (Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
          deltc1 = FLUX(I-1)/DXX(I)/DEN1(I)
          deltc1 = MAX(deltc1,0.0)
          Q0(I) = MAX (Q0(I), 1.E-20)

!         IF(I==84.AND.IJ==43.AND.K==1) PRINT*,Q0(I),deltc1,smother(3)
      

         IMONTHERBOUNDARY = 0.0 
         IF(NE==1) THEN
          IF(ITYPE==1) THEN
            IF(I==1) IMONTHERBOUNDARY=1
          ELSE           
            IF(I==1) IMONTHERBOUNDARY=1
          ENDIF      
         ENDIF
          IF(IMONTHERBOUNDARY==1) THEN
             CALL GetBounChange(deltc1,Q0(I),smthis,1,ismMax)
          ELSE
             CALL SMmixing(Q0(I),smthis,deltc1,smother,ismMax)     
          ENDIF  ! IMONTHERBOUNDARY 

          DO is=1,ismMax
            TMPsmconv(is,I)=smthis(is)
          ENDDO 
!          IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 1'                  
        ENDIF !FLUX(I-1) 

        QQ0 = Q0(I) + FLUX(I-1)/DXX(I)/DEN1(I) 
        QQ0 = MAX(QQ0, 1.E-20)

        IF(FLUXI<=0.) THEN
          DO is=1,ismMax    
           smthis(is) = TMPsmconv(is,I)
           smother(is)= TMPsmconv(is,I+1) 
          ENDDO 

!          deltc1 = (Q0(I)*DEN0(I)/2. - FLUXI)/DEN1(I)-Q0(I)/2.
          deltc1 = - FLUXI/DEN1(I)
          
         IMONTHERBOUNDARY = 0.0
         IF(NE==1) THEN
           IF(ITYPE==1) THEN       
             IF(I==NX) IMONTHERBOUNDARY=1      
           ELSE
             IF(I==NY) IMONTHERBOUNDARY=1      
           ENDIF  
         ENDIF  

         IF(IMONTHERBOUNDARY==1) THEN
           CALL GetBounChange(deltc1,QQ0,smthis,1,ismMax)
         ELSE
           CALL SMmixing(QQ0,smthis,deltc1,smother,ismMax)
         ENDIF

         DO is=1,ismMax
           TMPsmconv(is,I)=smthis(is)
         ENDDO
  
       ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!         IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 2'
        FLUX(I)= DXX(I)*(Q0(I)*DEN0(I)-QN(I)*DEN1(I)) + FLUX(I-1)
   ENDIF
  
 10 CONTINUE 

! update mixing ratios and limit Flux going down where u<0

  IF(U(ex).lt.ZR0)FLUX(ex)=Q0(ex+1)*U(ex)*DT*DD0(ex)
  DO 20 I=ex,sx,-1
  IF(U(I-1).GE.ZR0)THEN   
         if(U(I).lt.ZR0) QN(I) =    &     ! inflow-only cell
         (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUX(I-1)/DXX(I))/DEN1(I)  
!!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(FLUX(I)<=0) THEN

        DO is=1,ismMax
         smthis(is) = TMPsmconv(is,I)
         smother(is)= TMPsmconv(is,I+1)
        ENDDO 

!        deltc1=(Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)-Q0(I)/2.
         deltc1 = -FLUX(I)/DXX(I)/DEN1(I)
        deltc1 = MAX(deltc1,0.0)
        Q0(I) = MAX(Q0(I), 1.E-20)
              
        IMONTHERBOUNDARY = 0.0
        IF(NE==1) THEN
          IF(ITYPE==1) THEN       
           IF(I==NX) IMONTHERBOUNDARY=1
          ELSE 
           IF(I==NY) IMONTHERBOUNDARY=1   
          ENDIF
        ENDIF

        IF(IMONTHERBOUNDARY==1) THEN   
          CALL GetBounChange(deltc1,Q0(I),smthis,1,ismMax)       
        ELSE  
          CALL SMmixing(Q0(I),smthis,deltc1,smother,ismMax)
        ENDIF

        DO is=1,ismMax
         TMPsmconv(is,I)=smthis(is)
        ENDDO
         IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 3'

      ENDIF  ! FLUX(I)
        
        QQ0 = (Q0(I)*DEN0(I)- FLUX(I)/DXX(I))/DEN1(I) 
        QQ0 =MAX(1.E-20, QQ0)       
 
      IF(FLUX(I-1)>=0) THEN
        DO is=1,ismMax      
          smthis(is) = TMPsmconv(is,I)
          smother(is)= TMPsmconv(is,I-1)
        ENDDO  

!       deltc1=(Q0(I)*DEN0(I)/2.+FLUX(I-1)/DXX(I))/DEN1(I)-Q0(I)/2.
        deltc1 = FLUX(I-1)/DXX(I)/DEN1(I)
       deltc1 = MAX(deltc1,0.0)

       
       IMONTHERBOUNDARY=0.0
       IF(NE==1) THEN
        IF(ITYPE==1) THEN
         IF(I==1) IMONTHERBOUNDARY=1
        ELSE
         IF(I==1) IMONTHERBOUNDARY=1
        ENDIF 
       ENDIF 

       IF(IMONTHERBOUNDARY==1) THEN
         CALL GetBounChange(deltc1,QQ0,smthis,1,ismMax)
       ELSE
         CALL SMmixing(QQ0,smthis,deltc1,smother,ismMax)
       ENDIF

       DO is=1,ismMax
        TMPsmconv(is,I)=smthis(is)   
       ENDDO
!        IF(I==84.AND.IJ==43.AND.K==1) PRINT*,smthis(3),'inner inner 4'

     ENDIF ! FLUX(I-1)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
   ELSE
     X1=DT*ABS(U(I-1))/DXX(I)
     CF1=CFCALC_mark2(X1,Q0(I),Q0(I+1),Q0(I-1),IMXMN(I+1),IMXMN(I-1),LSHP)
     IF(U(I).ge.ZR0)CF1=Q0(I)   ! outflow only cell
        FLUXex1=DD0(I-1)*U(I-1)*DT*CF1/DXX(I)
        QN(I)=AMAX1( VCMIN(I) , AMIN1( VCMAX(I),  &
             (Q0(I)*DEN0(I)-FLUX(I)/DXX(I)+FLUXex1)/DEN1(I)))
!!!!!!!!!!!!!!!!!!!!!!!!SOURCE MARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(FLUXex1>=0) THEN
      DO is=1,ismMax      
        smthis(is) = TMPsmconv(is,I)
        smother(is)= TMPsmconv(is,I-1)
      ENDDO

!      deltc1 = (Q0(I)*DEN0(I)/2.+FLUXex1)/DEN1(I) - Q0(I)/2.
       deltc1 = FLUXex1/DEN1(I)
       deltc1 = MAX(deltc1,0.0)
       Q0(I) = MAX( Q0(I), 1.E-20)

      
      IMONTHERBOUNDARY=0.0
      IF(NE==1) THEN
       IF(ITYPE==1) THEN       
        IF(I==1) IMONTHERBOUNDARY=1      
       ELSE 
        IF(I==1) IMONTHERBOUNDARY=1       
       ENDIF
      ENDIF

      IF(IMONTHERBOUNDARY==1) THEN
        CALL GetBounChange(deltc1,Q0(I),smthis,1,ismMax)
      ELSE
        CALL SMmixing(Q0(I),smthis,deltc1,smother,ismMax)
      ENDIF

      DO is=1,ismMax
        TMPsmconv(is,I)=smthis(is)
      ENDDO
       
     ENDIF  !FLUX(I-1)

     QQ0 = Q0(I)+FLUXex1/DEN1(I) 
     QQ0 = MAX( QQ0, 1.E-20)
     
     IF(FLUX(I)<=0) THEN
       DO is=1,ismMax       
        smthis(is) = TMPsmconv(is,I)
        smother(is)= TMPsmconv(is,I+1)
       ENDDO 
       
!       deltc1 = (Q0(I)*DEN0(I)/2.-FLUX(I)/DXX(I))/DEN1(I)- Q0(I)/2.
        deltc1 = -FLUX(I)/DXX(I)/DEN1(I) 
        deltc1 = MAX(deltc1,0.0)       
       IMONTHERBOUNDARY=0.0
       IF(NE==1) THEN
         IF(ITYPE==1) THEN
          IF(I==NX) IMONTHERBOUNDARY=1       
         ELSE 
          IF(I==NY) IMONTHERBOUNDARY=1  
         ENDIF
       ENDIF

       IF(IMONTHERBOUNDARY==1) THEN
         CALL GetBounChange(deltc1,QQ0,smthis,1,ismMax)
       ELSE
         CALL SMmixing(QQ0,smthis,deltc1,smother,ismMax)
       ENDIF
      
       DO is=1,ismMax
         TMPsmconv(is,I)=smthis(is)
       ENDDO
       
     ENDIF ! FLUX(I)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
        FLUX(I-1)= DXX(I)*(QN(I)*DEN1(I)-Q0(I)*DEN0(I)) + FLUX(I)
   ENDIF 
 20 continue
   RETURN
   END
     
  FUNCTION CFCALC_mark2(X1,VCUP1,VCUP2,VCDW1,IMXUP2,IMXDN1,LSHP)
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
    CFCALC_mark2=MIN(MAX(CF,MIN(VCUP1,VCDW1)),MAX(VCUP1,VCDW1))
    RETURN
    END 
  



