
subroutine apm_4bdy(myid,lapm,ne,nx,ny,nzz,nest,sy,ey,sx,ex,ip3mem)
use apm_varlist !, only : apm_sulf,apm_salt,apm_dust,apm_bcoc
implicit none
include 'apm_parm.inc'
integer :: myid
logical ::lapm
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)
integer :: i,j,k,is
integer :: iapm,ixy
real :: value
integer :: ip3mem(nzz,nest)

!return

loop_k : do k=1,nzz-1

  if(lfor_sulf) then
   print*,'set bdy for sulf'
   loop_so4 : do is=1,NSO4
      iapm=ip_sulf(k,is,ne)
      if(lapm_nogdt) then
        call setwindbound( myid,apm_sulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                          ,nx(ne),ny(ne) ) 
      else
        value=0.0
        call set_apm_sbdy(myid,apm_sulf(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_nbdy(myid,apm_sulf(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_wbdy(myid,apm_sulf(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_ebdy(myid,apm_sulf(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
      endif
   enddo loop_so4
  endif

  if(lfor_salt) then
   print*,'set bdy for sea salt'
   do is=1,NSEA
      iapm=ip_salt(k,is,ne)
      if(lapm_nogdt) then
        call setwindbound( myid,apm_salt(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                          ,nx(ne),ny(ne) )
      else
        value=0.0
        call set_apm_sbdy(myid,apm_salt(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_nbdy(myid,apm_salt(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_wbdy(myid,apm_salt(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_ebdy(myid,apm_salt(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
      endif
   enddo
  endif

  if(lfor_dust) then
   print*,'set bdy for dust'
   do is=1,NDSTB
      iapm=ip_dust(k,is,ne)
      if(lapm_nogdt) then
        call setwindbound( myid,apm_dust(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                          ,nx(ne),ny(ne) )
      else
        value=0.0
        call set_apm_sbdy(myid,apm_dust(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_nbdy(myid,apm_dust(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_wbdy(myid,apm_dust(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_ebdy(myid,apm_dust(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
      endif
   enddo
  endif

  if(lfor_bcoc) then
   print*,'set bdy for bcoc'
   do is=1,NBCOCT
      iapm=ip_bcoc(k,is,ne)
      if(lapm_nogdt) then
        call setwindbound( myid,apm_bcoc(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                          ,nx(ne),ny(ne) )
      else
        value=0.0
        call set_apm_sbdy(myid,apm_bcoc(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_nbdy(myid,apm_bcoc(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_wbdy(myid,apm_bcoc(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
        call set_apm_ebdy(myid,apm_bcoc(iapm) &
                         ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
      endif
   enddo
  endif

  if(lcoated_dyn) then
   print*,'set bdy for coated sulfate'
   iapm=ip3mem(k,ne)
   if(lapm_nogdt) then
     call setwindbound( myid,msltsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                       ,nx(ne),ny(ne) )
     call setwindbound( myid,mdstsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                       ,nx(ne),ny(ne) )
     call setwindbound( myid,mbcsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                       ,nx(ne),ny(ne) )
     call setwindbound( myid,mocsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                       ,nx(ne),ny(ne) )

   else

   value=0.0

   call set_apm_sbdy(myid,msltsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_nbdy(myid,msltsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_wbdy(myid,msltsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_ebdy(myid,msltsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value) 

   call set_apm_sbdy(myid,mdstsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_nbdy(myid,mdstsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_wbdy(myid,mdstsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_ebdy(myid,mdstsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)

   call set_apm_sbdy(myid,mbcsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_nbdy(myid,mbcsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_wbdy(myid,mbcsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_ebdy(myid,mbcsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)

   call set_apm_sbdy(myid,mocsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_nbdy(myid,mocsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_wbdy(myid,mocsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)
   call set_apm_ebdy(myid,mocsulf(iapm) &
                    ,sx(ne),ex(ne),sy(ne),ey(ne),nx,ny,value)

   endif
  endif

enddo loop_k 
 
return

end subroutine apm_4bdy



subroutine apm_tbdy(myid,lapm,ne,nx,ny,nzz,nest,sy,ey,sx,ex,ip3mem)
use apm_varlist !, only : apm_sulf,apm_salt,apm_dust,apm_bcoc
implicit none
include 'apm_parm.inc'
integer :: myid
logical :: lapm
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)
integer :: i,j,k,is
integer :: iapm,ixy
real :: value
integer :: ip3mem(nzz,nest)

k=nzz

value=0.0

if(lfor_sulf) then
loop_so4 : do is=1,NSO4
 exit loop_so4
 iapm=ip_sulf(k,is,ne)
 call set_apm_tbdy(myid,apm_sulf(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
enddo loop_so4
endif

if(lfor_salt) then
do is=1,NSEA
 iapm=ip_salt(k,is,ne)
 call set_apm_tbdy(myid,apm_salt(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
enddo
endif

if(lfor_dust) then
do is=1,NDSTB
 iapm=ip_dust(k,is,ne)
 call set_apm_tbdy(myid,apm_dust(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
enddo
endif

if(lfor_bcoc) then
do is=1,NBCOCT
 iapm=ip_bcoc(k,is,ne)
 call set_apm_tbdy(myid,apm_bcoc(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
enddo
endif


if(lcoated_dyn) then
 iapm=ip3mem(k,ne)
 call set_apm_tbdy(myid,msltsulf(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
 call set_apm_tbdy(myid,mdstsulf(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
 call set_apm_tbdy(myid,mbcsulf(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
 call set_apm_tbdy(myid,mocsulf(iapm),sy(ne),ey(ne),sx(ne),ex(ne),value)
endif

end subroutine apm_tbdy





subroutine set_apm_sbdy(myid,c,sx,ex,sy,ey,nx,ny,d)
  implicit none
  integer myid, sx, ex, sy, ey, nx, ny
  real*8, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real :: d
  integer i, j
  if(sy==1)then         ! south boundary
    do i=sx,ex
      c(i,sy-1) = d
    enddo
   endif
  return
end subroutine set_apm_sbdy


subroutine set_apm_nbdy(myid,c,sx,ex,sy,ey,nx,ny,d)
  implicit none
  integer myid, sx, ex, sy, ey, nx, ny
  real*8, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real :: d
  integer i, j
  if(ey==ny)then         ! north boundary
    do i=sx,ex
      c(i,ey+1) = d
    enddo
  endif
end subroutine set_apm_nbdy


subroutine set_apm_wbdy(myid,c,sx,ex,sy,ey,nx,ny,d)
  implicit none
  integer myid, sx, ex, sy, ey, nx, ny
  real*8, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real :: d
  integer i, j
  if(sx==1)then         ! west boundary
    do j=sy,ey
       c(sx-1,j) = d
    enddo
  endif
end subroutine set_apm_wbdy


subroutine set_apm_ebdy(myid,c,sx,ex,sy,ey,nx,ny,d)
  implicit none
  integer myid, sx, ex, sy, ey, nx, ny
  real*8, dimension(sx-1:ex+1,sy-1:ey+1) :: c
  real :: d
  integer i, j
  if(ex==nx)then         ! east boundary
    do j=sy,ey
      c(ex+1,j) = d
    enddo
  endif
end subroutine set_apm_ebdy


subroutine set_apm_tbdy(myid,c,sy,ey,sx,ex,d)
 implicit none
 integer :: myid
 integer :: sy,ey,sx,ex
 integer :: i,j,k,is
 integer :: iapm,ixy
 real :: d
 real*8,dimension(sx-1:ex+1,sy-1:ey+1) :: c

 do j=sy,ey
 do i=sx,ex
    c(i,j)=d
 enddo
 enddo

 return
end subroutine set_apm_tbdy

