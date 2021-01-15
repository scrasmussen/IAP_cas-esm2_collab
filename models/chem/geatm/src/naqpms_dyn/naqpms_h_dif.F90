
subroutine naqpms_h_dif &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,ktop &
 &  ,mem3d &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,ifsm,idmSet,ismMax,igMark)

use naqpms_varlist
use naqpms_gridinfo
use met_fields
implicit none

integer :: myid

real :: dt


integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)

integer :: i,j,k,is
!integer :: k,is

integer :: mem3d

integer :: mem2d
real,dimension(mem2d) :: ktop


integer :: ixy,i02,i03,i04,i04aer

integer :: ig,iduc,ia,i05,i05c


real,dimension(mem3d) :: wk

integer :: letdoit,idm,ism,i04sm

real,allocatable,dimension(:,:,:) :: sm

! do ig=1,igas    ! for gas phase
do ig=1,iedgas

  do k=1,nzz-1
    i03=ip3mem(k,ne)
    i04=ip4mem(k,ig,ne)
    i02=ip2mem(ne)
    call  dif_hori( myid, gas(i04), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,ktop(i02),k,ig)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  if(ifsm(ne)==1)then

    letdoit=0

    do idm=1,idmSet
       if(igMark(idm)==ig) letdoit=idm
    enddo

    if(letdoit>0)then

      allocate(sm(ismMax,sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1))

      do k=1,nzz-1

        i03=ip3mem(k,ne)
        i04=ip4mem(k,ig,ne)
        do ism=1,ismMax
           i04sm=ipSMmem(k,ism,letdoit,ne)
           do i=sx(ne)-1,ex(ne)+1
           do j=sy(ne)-1,ey(ne)+1
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
               sm(ism,i,j)= SourceMark(i04sm+ixy)
           enddo
           enddo
        enddo

        i02=ip2mem(ne)

        call dif_hori_mark( myid, gas(i04), kh(i03), dx(i03),dy(i03), &
                  sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt,&
                  ismMax,sm, ne,k,ktop(i02))

        do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
          do i=sx(ne)-1,ex(ne)+1
          do j=sy(ne)-1,ey(ne)+1
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
               SourceMark(i04sm+ixy)=sm(ism,i,j)
          enddo
          enddo
        enddo

      enddo
      if(allocated(sm)) deallocate(sm)
    endif

   endif
   !!!!!!!!!!!!!!!!!!!!!!!!

 enddo    !ig



if(laerv2) then
 do ia=1,naersp
 do is=1,naerbin

   if(ia.gt.1.and.is.gt.1) cycle ! skip zero aersol tracer

   do k=1,nzz-1
     i03=ip3mem(k,ne)
     i04aer=ip4mem_aer(k,is,ia,ne)
     call  dif_hori( myid, aerom(i04aer), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt)
   enddo
 enddo
 enddo
endif



 do ia=1,iaer
 do is=1,isize

   do k=1,nzz-1
     i03=ip3mem(k,ne)
     i05=ip5mem(k,is,ia,ne)
     call  dif_hori( myid, aer(i05), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt)
   enddo

   IF(ia==1) THEN ! For sea salt
     do iduc = 1, nseacom
     do k = 1, nzz-1
       i03=ip3mem(k,ne)
       i05c = ip5memcs (k,is,iduc,ne)
       call dif_hori( myid, SEACOMP(i05c), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt)
     enddo
     enddo ! iduc
   ELSE IF(IA==2) THEN ! for dust
     do iduc = 1, ndustcom
     do k = 1, nzz-1
       i03=ip3mem(k,ne)
       i05c = ip5memc (k,is,iduc,ne)
       call  dif_hori( myid, DUSTCOMP(i05c), kh(i03), dx(i03), dy(i03), &
            sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt)
     enddo
     enddo ! iduc
   ENDIF ! IA IF

  enddo        !isize
  enddo        !iaer


end subroutine naqpms_h_dif





