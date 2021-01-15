            subroutine caltheta(myid,theta,t,sx,ex,sy,ey,nx,ny,k)
             integer :: myid,sx,sy,ex,ey,nx,ny,k
             real,dimension(sx-1:ex+1,sy-1:ey+1) :: theta,t
             integer :: i,j
             real,dimension(20) ::p
             DATA p/1001.,1000.,950.,900.,850.,800.,750.,700.,650.,     &
                  600.,550.,500.,450.,400.,350.,300.,250.,200.,150.,100. /

               do j=sy,ey
                 do i=sx,ex
                   theta(i,j)=t(i,j)*(1000./p(k))**0.286
!                   print*,theta(i,j),i,j,sx,t(i,j)
                 enddo
               enddo

              return

            end 





