! Zifa 2006/10/05 revised
! For Source Mark
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  subroutine dif_hori_mark( myid, c0, kh, deltx, delty, &
                       sx, ex, sy, ey ,nx,ny,dt,   &
                      ismMax,sm, ne,k,ktop)
  integer,parameter :: ismMaxSet=500
  integer myid, sx, ex, sy, ey, it
  real, dimension(sx-1:ex+1,sy-1:ey+1) :: c0, kh, deltx, delty,ktop
  real, dimension(ismMax,sx-1:ex+1,sy-1:ey+1) :: sm

  integer i, j,k
  integer sx1,ex1,sy1,ey1
  real, allocatable, dimension(:,:) :: c
  real smthis(ismMaxSet)
  real smother(ismMaxSet)

  allocate(c(sx-1:ex+1,sy-1:ey+1))

  c=c0

  sx1=sx
  ex1=ex
  sy1=sy
  ey1=ey

    do j = sy1,ey1
      do i = sx1,ex1
         IF(k>ktop(i,j)) then
          do is =1 , ismMax
          if(is==2) then
               sm(is,i,j)=1.0
          else
               sm(is,i,j)=0.0
          endif
          enddo
         ENDIF
      enddo
    enddo


  
  do j = sy1,ey1
    do i = sx1,ex1

      deltc1=dt*kh(i,j)*(c(i+1,j)-c(i,j))/(deltx(i,j)*deltx(i,j))
       
  
      if(deltc1 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i+1,j)
         enddo
          imotherboundary =0  
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(i==ex1) then 
            call GetBounChange(deltc1,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc1,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       

      deltc2=dt*kh(i,j)*(c(i-1,j)-c(i,j))/(deltx(i,j)*deltx(i,j))
     if(deltc2 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax    
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i-1,j)
         enddo 
          imotherboundary =0
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1 
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(i==sx1) then
            call GetBounChange(deltc2,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc2,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      c(i,j)=c(i,j)+deltc1+deltc2
    enddo
  enddo

  do j = sy1,ey1
    do i = sx1,ex1

     deltc1=dt*kh(i,j)*(c(i,j+1)-c(i,j))/(delty(i,j)*delty(i,j))

     if(deltc1 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i,j+1)
         enddo
          imotherboundary =0
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(j==ey1) then
            call GetBounChange(deltc1,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc1,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     deltc2=dt*kh(i,j)*(c(i,j-1)-c(i,j))/(delty(i,j)*delty(i,j))

     if(deltc2 > 0 ) then !!!!!!!!!!!!!!!!!!
      c00=max(c(i,j),1.E-20)
         do is=1,ismMax
           smthis(is)=sm(is,i,j)
           smother(is)=sm(is,i,j-1)
         enddo
          imotherboundary =0
          ! to judge if it is mother domain
            if(ne==1)then
              if(i==1)imotherboundary =1
              if(i==nx)imotherboundary =1
              if(j==1)imotherboundary =1
              if(j==ny)imotherboundary =1
            endif

          if(imotherboundary==1)then
!--------------------------lijie modify for sm----------------------
            IF(j==sy1) then
            call GetBounChange(deltc2,c00,smthis,1,ismMax)
            ENDIF
!-------------------------finish------------------------------------
           else
            call SMmixing(c00,smthis,deltc2,smother,ismMax)
          endif
          do is=1,ismMax
           sm(is,i,j)=smthis(is)
          enddo
      endif  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      c(i,j)=c(i,j)+deltc1+deltc2
    enddo
  enddo

 deallocate(c)


  return
  end

