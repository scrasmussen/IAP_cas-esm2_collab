       subroutine dif_vert_conv_mark( myid,c1m,cm,cp1m,c0m,kz,kz1,kzb,deltz,deltz1,&
                                deltzt,deltzb,pbl,&
                                xmol,heiz,heiz1,heizb,heizt,TERRAIN,&
                              sx, ex, sy, ey,dt,ig,k,nzz,ismMax,sm,smb,smp,sm0)
      integer,parameter :: ismMaxSet=100
      integer myid, sx, ex, sy, ey, it,ig
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: kz,pbl,xmol,heiz,kz1,heiz1,heizb,heizt,kzb,TERRAIN
      !kz1 and heiz1 is the kz and height of l layers,kzb and heizb is the i-1/2 layers,heizt i+1/2 layer
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1m,cm,cp1m,c0m !c0m is the mixing ratios in the 1st layer
      real, dimension(sx-1:ex+1,sy-1:ey+1) :: deltz,deltz1,deltzt,deltzb !delta1 is depth of 1st layer      
      real :: M2u,M2di,M2dit,fconv,kzz,kzzb 
      real :: const1,const2,karman,alpha1
      integer i, j, k     
      real :: smthis(ismMaxSet) 
      real :: smother(ismMaxSet) 
      real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1):: sm,smb,smp,sm0

      real, allocatable, dimension(:,:) ::c1,c,cp1,c0
      allocate(c1(sx-1:ex+1,sy-1:ey+1))
      allocate(c(sx-1:ex+1,sy-1:ey+1))
      allocate(cp1(sx-1:ex+1,sy-1:ey+1))
      allocate(c0(sx-1:ex+1,sy-1:ey+1))

      c1=c1m
      c=cm
      cp1=cpm
      c0=c0m

! this scheme is using ACM2 scheme(non-local scheme)
      do j = sy,ey
          do i = sx,ex
             karman=0.4
             alpha1=7.2
             const1=(karman**(-2./3.))/(0.1*alpha1)
              if(xmol(i,j)<0) then
             const2=(-pbl(i,j)/(xmol(i,j)+1.e-6))**(-1./3.)
             fconv=1./(1.+const1*const2)
              else 
             fconv=0.
              endif 
              if((heiz(i,j)-TERRAIN(i,j))>pbl(i,j)) fconv=0.0 !confirm the highest layer in pbl
              fconv=min(fconv,0.3)
               fconv =0.      
             M2u=fconv*kz1(i,j)/(deltz1(i,j)*(pbl(i,j)-heiz1(i,j)+TERRAIN(i,j)))
             M2di=M2u*(pbl(i,j)-heizb(i,j)+TERRAIN(i,j))/deltz(i,j)
             M2dit=M2u*(pbl(i,j)-heiz(i,j)+TERRAIN(i,j))/deltzt(i,j)
             kzz=kz(i,j)*(1.-fconv)
             kzzb=kzb(i,j)*(1.-fconv)
             dt0=deltz(i,j)*deltz(i,j)/(10.*kz(i,j))
             iii= dt/dt0+1
             dtt=dt/iii
             do it=1,iii
                 c00=c(i,j)
                 deltam01=M2u*c0(i,j)
                 deltam02=-M2di*c(i,j)
                 deltam03= M2dit*cp1(i,j)*deltzt(i,j)/deltz(i,j)
                 deltam04a= (kzz*(cp1(i,j)-c(i,j))/deltz(i,j))/deltz(i,j)
                 deltam05a= -(kzzb*(c(i,j)-c1(i,j))/deltzb(i,j))/deltz(i,j)
                 deltam04b = (kzz*(cp1(i,j))/deltz(i,j))/deltz(i,j)
                 deltam05b = -(kzzb*(-c1(i,j))/deltz(i,j))/deltz(i,j)
               if(k==1) then                 
                  c(i,j)=c(i,j)+dtt*(deltam04+deltam05+M2dit*(cp1(i,j)-c(i,j)))
               else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j)) then
                  c(i,j)=c(i,j)+dtt*(deltam04+deltam05+deltam01+deltam02)
               else      
                  c(i,j)=c(i,j)+dtt*(deltam01+deltam02+deltam03+deltam04+deltam05)
               endif
                 
!----------------to estimate the contribution from 1st layers--------------
              if(k==1) then                 
                deltc1=0.0
              else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j)) then
                deltc1=dtt*(deltam01+deltam02)
!                 deltc1=0.0
              else      
                deltc1=dtt*(deltam01+deltam02)
              endif
               IF (deltc1>0.0) then
                do is=1,ismMax
                 smthis(is)=sm(is,i,j)
                 smother(is)=sm0(is,i,j)
                enddo
                if(k>nzz-1) then
                 smthis(2)=1.0
                else
                  call SMmixing(c00,smthis,deltc1,smother,ismMax)
                endif
                                                      
                 do  is=1,ismMax
                   sm(is,i,j)=smthis(is)
                   sm(is,i,j)=max(smthis(is),0.)
                   sm(is,i,j)=min(smthis(is),1.)
                 enddo
                 
                ENDIF
             c00=c00+deltc1
!!----------------to estimate the diffusion contribution from next lower layers--------------
             if(k==1) then
                 deltc2a=0.0
                 deltc2b=0.0
             else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j))  then
                 deltc2a=dtt*deltam05a
                 deltc2b=dtt*deltam05b
             else
                 deltc2a=dtt*deltam05a
                 deltc2b=dtt*deltam05b
             endif    

               IF(deltc2a>0.0) THEN
                 do is=1,ismMax
                   smthis(is)=sm(is,i,j)
                   smother(is)=smb(is,i,j)
                 enddo
                 if(k>nzz-1) then
                   smthis(2)=1.0
                 else
                   call SMmixing(c00,smthis,deltc2a,smother,ismMax)
                 endif

                 do  is=1,ismMax
                    sm(is,i,j)=smthis(is)
                    sm(is,i,j)=max(smthis(is),0.)
                    sm(is,i,j)=min(smthis(is),1.)
                 enddo
               ENDIF        
              c00=c00+deltc2a
!!!----------------to estimate the convective contribution from next higher layers--------------

             if(k==1) then
                deltc3=dtt*deltam03
             else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j))   then
                deltc3=0.0
             else
                deltc3=dtt*deltam03
             endif

             IF(deltc3>0.0) THEN
              do is=1,ismMax
                 smthis(is)=sm(is,i,j)
                 smother(is)=smp(is,i,j)
              enddo
              if(k==nzz-1)then
                 call GetBounChange(deltc3,c00,smthis,2,ismMax)
              else if(k>nzz-1) then
                 smthis(2)=1.0
              else
                 call SMmixing(c00,smthis,deltc3,smother,ismMax)
              endif 

              do  is=1,ismMax
                sm(is,i,j)=smthis(is)
                sm(is,i,j)=max(smthis(is),0.)
                sm(is,i,j)=min(smthis(is),1.)
              enddo
             ENDIF               
              c00=c00+deltc3
!!!----------------to estimate the diffusion  contribution from next higher layers--------------             

             if(k==1) then
                deltc4a=dtt*deltam04a
                deltc4b=dtt*deltam04b
             else if((heizt(i,j)-TERRAIN(i,j))>pbl(i,j).and.(heiz(i,j)-TERRAIN(i,j))<pbl(i,j))    then
                deltc4a=dtt*deltam04a
                deltc4b=dtt*deltam04b                
             else
                deltc4a=dtt*deltam04a
                deltc4b=dtt*deltam04b                
             endif
                                                       
             IF(deltc4a>0.0) THEN
              do is=1,ismMax
                  smthis(is)=sm(is,i,j)
                  smother(is)=smp(is,i,j)
              enddo
              if(k==nzz-1)then
                call GetBounChange(deltc4a,c00,smthis,2,ismMax)      
              else if(k>nzz-1) then
                smthis(2)=1.0      
              else
                call SMmixing(c00,smthis,deltc4a,smother,ismMax)
              endif
               
              do  is=1,ismMax
                sm(is,i,j)=smthis(is)
                sm(is,i,j)=max(smthis(is),0.)
                sm(is,i,j)=min(smthis(is),1.)
              enddo
             ENDIF
             c00=c00+deltc4a
!--------------------------             

             c(i,j)=max(0.,c(i,j))
 
          enddo !it 
         enddo !i
        enddo  !j
          deallocate(c1,c,cp1,c0)
        return
       end                       
