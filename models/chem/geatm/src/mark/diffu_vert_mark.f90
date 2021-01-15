            subroutine diffu_vert_mark(myid,&
              c1tmp,c_1tmp,ctmp,cp1tmp,&              
              kv1,kv_1,kv,kvp1,h1,h_1,h,hp1,terrain,&
              dep1,dep_1,dep,depp1,&
              fcon,pbl,kpbl,Gc_MOLWT,ig,igas,&
              ismMax,sm1, sm,smb,smp,ktop,&
              sx,ex,sy,ey,dt,k)
              
            integer, parameter :: ismMaxSet=500
            real smthis(ismMaxSet)
            real smother(ismMaxSet)
            integer  :: myid, sx, ex, sy, ey,it
            real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1) :: sm1,sm,smb,smp
            
            
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: c1tmp,c_1tmp,ctmp,cp1tmp
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: kv1,kv_1,kv,kvp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: h1,h_1,h,hp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: dep1,dep_1,dep,depp1
            real, dimension(sx-1:ex+1,sy-1:ey+1) :: fcon,pbl,terrain
            integer, dimension(sx-1:ex+1,sy-1:ey+1) :: kpbl,ktop
            real :: M2u,M2di,M2dii,cm1tmp,cm_1tmp,cmtmp,cmp1tmp
            real :: dt,dtt,dt0,delta1,delta2,delta3,delta4,delta5
            real :: GC_MOLWT(igas)
           do j = sy,ey
           do i = sx,ex
           
            M2u= fcon(i,j)*kv1(i,j)/dep1(i,j)/(pbl(i,j)-h1(i,j)+terrain(i,j)-dep1(i,j)/2.)
            M2di=M2u*(pbl(i,j)-h(i,j)+terrain(i,j)+dep(i,j)/2.)/dep(i,j)
            M2dii=M2u*(pbl(i,j)-hp1(i,j)+terrain(i,j)+depp1(i,j)/2.)/depp1(i,j)
           

            dtt=dt
!------convert ppbv --> mass mixing ratio------

             cm1tmp=c1tmp(i,j)*GC_MOLWT(ig)/29.
             cm_1tmp=c_1tmp(i,j)*GC_MOLWT(ig)/29.
             cmtmp=ctmp(i,j)*GC_MOLWT(ig)/29.
             cmp1tmp=cp1tmp(i,j)*GC_MOLWT(ig)/29.
!-----------------------convert  over---------------------                  
!              if(i==46.and.j==32.and.k==2.and.ig==17.and.ne==1)print*,cm_1,cm,cmp1,'222' 
           delta1=0.0
           delta2=0.0
           delta3=0.0
           delta4=0.0
           delta5=0.0
           IF(kpbl(i,j)==1) goto 999
           if(k==1) then
            delta1=-cm1tmp*M2u*(pbl(i,j)-h1(i,j)+terrain(i,j)-dep1(i,j)/2.)/dep1(i,j)
            delta2=0.0
            delta3= M2dii*cmp1tmp*depp1(i,j)/dep(i,j)
           else if(k==kpbl(i,j)) then
            delta1= M2u*cm1tmp
            delta2=-M2di*cmtmp
            delta3= 0.0      
           else if(k<kpbl(i,j)) then 
            delta1= M2u*cm1tmp
            delta2=-M2di*cmtmp
            delta3= M2dii*cmp1tmp*depp1(i,j)/dep(i,j)
           endif 
999         continue
           if(k<=kpbl(i,j)) then
           else
            delta1=0.0;delta2=0.0;delta3=0.0 
           endif
             delta4=0.0; delta5=0.0

             c00=cmtmp
           if(delta3>0.0) then
             do is=1,ismMax
               smthis(is)=sm(is,i,j)
               smother(is)=smp(is,i,j)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta3'
             enddo
                 
             if(k==ktop(i,j))then                   
               call GetBounChange(delta3*dtt,c00,smthis,2,ismMax)
             else if(k>ktop(i,j)) then
              do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0 
               endif
              enddo
             else
               call SMmixing(c00,smthis,delta3*dtt,smother,ismMax)
             endif

             do is=1,ismMax
               sm(is,i,j)=smthis(is)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),'delta3'    
             enddo 
           endif !delta3

            c00=c00+delta3*dtt 

            if(delta2>0.0) then
            do is=1,ismMax
              smthis(is)=sm(is,i,j)
              smother(is)=smb(is,i,j)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta2'
            enddo       
  
            if(k>ktop(i,j)) then              
               do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0
               endif
               enddo
            else
              call SMmixing(c00,smthis,delta2*dtt,smother,ismMax)      
            endif    

            do is=1,ismMax                               
              sm(is,i,j)=smthis(is)                           
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),'delta2'
            enddo  
           endif !delta2 
             c00=c00+delta2*dtt
  
           
           if(delta1>0.0) then     
            do is=1,ismMax                                              
             smthis(is)=sm(is,i,j)                              
             smother(is)=sm1(is,i,j)                                  
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta1'
            enddo      

            if(k>ktop(i,j)) then                                 
               do is = 1 ,ismMax
               if(is==2) then
                 smthis(is)=1.0
               else
                 smthis(is) =0.0
               endif
               enddo
            else  
             call SMmixing(c00,smthis,delta1*dtt,smother,ismMax)    
            endif          
                    
            do is=1,ismMax
              sm(is,i,j)=smthis(is)
               IF(i==46.and.j==55.and.k==1.and.is==3) PRINT*, smthis(is),smother(is),'delta1'
            enddo          

           endif !delta1
             c00=c00+delta1*dtt

                
           enddo !i
           enddo !j
            return 
            end


