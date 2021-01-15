
subroutine print_test(prgname,varname,ne,nest,nzz,sy,ey,sx,ex,ip3mem)

use apm_varlist
implicit none
include 'apm_parm.inc'
character :: prgname*50
character :: varname*50
integer :: is,i,j,k
integer :: ixy,iapm
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)
integer :: ip3mem(nzz,nest)

 
print*,'start test '//trim(prgname)//'-'//trim(varname)

IF(trim(varname).eq.'sulf') THEN
 do is=1,NSO4
   print*,'size=',is
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip_sulf(k,is,ne)
         print*,k,j,i,apm_sulf(iapm+ixy)
      enddo
   enddo
   enddo
 enddo 
ELSEIF(trim(varname).eq.'salt') THEN
 do is=1,NSEA
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip_salt(k,is,ne)
         print*,k,j,i,apm_salt(iapm+ixy)
      enddo
   enddo
   enddo
 enddo
ELSEIF(trim(varname).eq.'dust') THEN
 do is=1,NDSTB
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip_dust(k,is,ne)
         print*,k,j,i,apm_dust(iapm+ixy)
      enddo
   enddo
   enddo
 enddo
ELSEIF(trim(varname).eq.'bcoc') THEN
 do is=1,NBCOCT
   do j = sy(ne),ey(ne)
   do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip_bcoc(k,is,ne)
         print*,k,j,i,apm_bcoc(iapm+ixy)
      enddo
   enddo
   enddo
 enddo
ELSEIF(trim(varname).eq.'sltsulf') THEN
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip3mem(k,ne)
         print*,k,j,i,msltsulf(iapm+ixy)
      enddo
 enddo
 enddo
ELSEIF(trim(varname).eq.'dstsulf') THEN
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip3mem(k,ne)
         print*,k,j,i,mdstsulf(iapm+ixy)
      enddo
 enddo
 enddo
ELSEIF(trim(varname).eq.'bcsulf') THEN
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip3mem(k,ne)
         print*,k,j,i,mbcsulf(iapm+ixy)
      enddo
 enddo
 enddo
ELSEIF(trim(varname).eq.'ocsulf') THEN
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do k=1,nzz-1
         iapm=ip3mem(k,ne)
         print*,k,j,i,mocsulf(iapm+ixy)
      enddo
 enddo
 enddo
ENDIF

print*,'stop test '//trim(prgname)//'-'//trim(varname)
stop



end subroutine print_test
