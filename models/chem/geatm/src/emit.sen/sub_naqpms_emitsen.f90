
subroutine naqpms_snstvy_emit &
 & ( myid &
 &  ,ictg,ig &
 &  ,NLAY_EM &
 &  ,NEMIT   &
 &  ,img     &
 &  ,IPIG    &
 &  ,ne,nx,ny,nzz,nest,sx,ex,sy,ey &
 &  ,igas &
 &  ,ip_emit2d,mem_emt2d,emit2d &
 &  ,emsfrc )

implicit none

integer :: myid

integer :: iemittype
integer :: i,j,k,ig,i02Gas

integer :: ixy

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: NEMIT,NLAY_EM
integer :: IPIG(NEMIT),img
real    :: emsfrc(NEMIT,NLAY_EM)

integer :: igas
integer,dimension(igas,nest) :: ip2memGas

integer :: mem_emt2d
real,dimension(mem_emt2d) :: emit2d



 i02emt = ip_emit2d(ig,ictg,ne)

 do j=sy(ne),ey(ne)
 do i=sx(ne),ex(ne)

   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

   emit2d(i02emt+ixy)=emit2d(i02emt+ixy)*emsfrc(img,ictg)

 enddo
 enddo

!stop

end subroutine naqpms_snstvy_emit


