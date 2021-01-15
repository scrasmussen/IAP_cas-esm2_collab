      subroutine calpv(myid,pv,u,u0,u00,v,v0,v00,theta,theta0,theta00,&
                 xlat,deltx,delty,sx,ex,sy,ey,nx,ny,k,it,ne)
        integer :: myid,sx,ex,sy,ey,nx,ny,k,it
        real,dimension(sx-1:ex+1,sy-1:ey+1) :: pv,u,u0,u00,v,v0,v00,theta,&
                                       theta0, theta00,deltx,delty
        real,dimension(sx-1:ex+1,sy-1:ey+1) :: vorp,dthedp,dthedx,dthedy, &
                                       xlat,dudp,dvdp,f

        real,dimension(20) :: p
        DATA p/1001,1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,&
               500.,450.,400.,350.,300.,250.,200.,150.,100. /
        integer :: i,j   
        real :: g,scale
        parameter (g=9.81,scale=-1.E6)
        integer :: k0,k2


        if(k==1) then
          k0=k
          k2=k+1
         elseif(k==20) then
          k0=k-1
          k2=k
         else
          k0=k-1
          k2=k+1
         endif
      

        !*************compute theta& vor & f*****************
loop11:   do  i=sx,ex
loop12:     do  j=sy,ey
!          theta(i,j,k0)=t(i,j,k0)*(1000./p(k0))**0.286
!          theta(i,j,k2)=t(i,j,k2)*(1000./p(k1))**0.286

         IF(i/=ex.and.j/=ey) THEN
          i0=i
          j0=j
          i1=i+1
          j1=j+1
          f(i,j)=2.*7.29E-05*sin(xlat(i,j)*3.14/180.)
          vorp(i,j)= (-(u(i0,j1)-u(i0,j0)+u(i1,j1)-u(i1,j0)) &
                        +(v(i1,j0)-v(i0,j0)+v(i1,j1)-v(i0,j1)))&
                       /(deltx(i,j)*2.)+f(i,j)

         ELSE
          f(i,j)=2.*7.29E-05*sin(xlat(i,j)*3.14/180.)
          vorp(i,j)= (-(u(i0,j0)-u(i0,j0-1)) &
                        +(v(i0,j0)-v(i0-1,j0)))&
                       /deltx(i,j)+f(i,j)

         ENDIF

        enddo  loop12
       enddo   loop11


!     print*,vorp(46,32,k) , f(46,32),k,u(46,33,k)
     

     !**************compute dthedx & dthe/dy
loop21:   do  i=sx,ex
loop22:     do  j=sy,ey
           if(i/=sx.and.i/=ex) then
            i0=i-1
            i1=i+1
            dthedx(i,j)=(theta(i1,j)-theta(i0,j))/(2.*deltx(i,j))
           else if(i==sx) then
            i0=i
            i1=i+1
            dthedx(i,j)=(theta(i1,j)-theta(i0,j))/(deltx(i,j))
           else if(i==ex) then
            i0=i-1
            i1=i
            dthedx(i,j)=(theta(i1,j)-theta(i0,j))/(deltx(i,j))
           endif

           IF(j/=sy.and.j/=ey) then
            j0=j-1
            j1=j+1
            dthedy(i,j)=(theta(i,j1)-theta(i,j0))/(2.*deltx(i,j))
           ELSE IF(j==sy) then
            j0=j
            j1=j+1
            dthedy(i,j)=(theta(i,j1)-theta(i,j0))/(deltx(i,j))
           ELSE IF(j==ey) then
            j0=j-1
            j1=j
            dthedy(i,j)=(theta(i,j1)-theta(i,j0))/(deltx(i,j))
           ENDIF
         enddo loop22
       enddo loop21

!*****************compute dthe/dp******
loop31:  do  i=sx,ex
loop32:     do  j=sy,ey
             dthedp(i,j)=(theta00(i,j)-theta0(i,j))/(p(k2)-p(k0))
            enddo loop32
         enddo loop31
!           print*,dthedp(46,32,k),k

!************compute du/dp&dv/dp********
loop41:  do  i=sx,ex
loop42:     do  j=sy,sy
         if(i/=ex.and.j/=ey) then
           i0=i
           i1=i+1
           j0=j
           j1=j+1
           dudp(i,j)=0.25*(u00(i1,j1)-u0(i1,j1)+ &
                             u00(i1,j0)-u0(i1,j0)+ &
                             u00(i0,j1)-u0(i0,j1)+ &
                             u00(i0,j0)-u0(i0,j0)) &
                             /(p(k2)-p(k0))
           dvdp(i,j)=0.25*(v00(i1,j1)-v0(i1,j1)+ &
                             v00(i1,j0)-v0(i1,j0)+ &
                             v00(i0,j1)-v0(i0,j1)+ &
                             v00(i0,j0)-v0(i0,j0)) &
                             /(p(k2)-p(k0))
         else
           dudp(i,j)=(u00(i,j)-u0(i,j))/(p(k2)-p(k0))
           dvdp(i,j)=(v00(i,j)-v0(i,j))/(p(k2)-p(k0))
         endif
         enddo loop42
      enddo loop41

!***************compute pv**************
loop51:  do   i=sx,ex
loop52:    do  j=sy,ey
               pv(i,j)=g*scale/100.*                     &
                   (vorp(i,j)*dthedp(i,j)+dudp(i,j)*dthedy(i,j)    &
                    -dvdp(i,j)*dthedx(i,j))
           enddo loop52
         enddo loop51
     
       
        return
       end 






