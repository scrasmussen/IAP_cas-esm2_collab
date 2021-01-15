
subroutine naqpms_h_adv &
 & ( myid &
 &  ,imasskeep &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem3d,RatioMass &
 &  ,mem2d,ktop &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use naqpms_gridinfo, only : dx,dy
use met_fields, only : u,v
implicit none

integer :: myid

real :: dt

logical :: lapm
integer :: imasskeep

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is

integer :: mem3d

real,dimension(mem3d) :: RatioMass,kpmass_m2

integer :: mem2d
real,dimension(mem2d) :: ktop

integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)

integer :: ixy,i03,iapm



integer :: igas,iaer,isize,nseacom,ndustcom


real,dimension(mem3d) :: wk

integer :: ig,i04,i02,ia,iduc,i05,i05c,i04aer
integer :: idm,ism,i04sm

integer :: letdoit

real,allocatable,dimension(:,:,:) :: sm

logical :: lflag 
real :: uws,vws,deltx,delty

real :: rrr

if(1==2) then
 do j = sy(ne)-1,ey(ne)+1
 do i = sx(ne)-1,ex(ne)+1

if((j.eq.sy(ne)-1.or.j.eq.ey(ne)+1).and.(i.eq.sx(ne)-1.or.i.eq.ex(ne)+1)) cycle

    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    do k=1,nzz
        i03=ip3mem(k,ne)
        uws=u(i03+ixy)
        vws=v(i03+ixy)
        deltx=dx(i03+ixy)
        delty=dy(i03+ixy)

rrr=3**(ne-1)

        lflag=(abs(uws).lt.1.0e3).and.(abs(vws).lt.1.0e3).and.&
             &(deltx.eq.27000/rrr).and.(delty.eq.27000/rrr)
        if(.not.lflag) then
          print*,'erro-hadv',sx(ne),ex(ne),sy(ne),ey(ne)
          print*,'erro-hadv11',i,j,k
          print*,'erro-hadv22',uws, vws,deltx,delty,27000/rrr 
        endif
    enddo
 enddo
 enddo
endif


!return

!loop_gas_spc : do ig=1,igas    ! for gas phase
loop_gas_spc : do ig=1,iedgas

    if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i04=ip4mem(k,ig,ne)
            kpmass_m2(i03+ixy)=gas(i04+ixy)
         enddo
      enddo
      enddo
    endif

    do k=1,nzz-1
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)
       call  adv_hori(myid,gas(i04),u(i03),v(i03),dx(i03),dy(i03),&
          sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne),dt)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!

    if(imasskeep==1)then
     do k=1,nzz-1   !mass conservation 
       i03=ip3mem(k,ne)
       i04=ip4mem(k,ig,ne)

       call  balance(myid,gas(i04),kpmass_m2(i03),RatioMass(i03), &
                 sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
     enddo
    endif


   !===========================================
   ! for Source Mark
    if(ifsm(ne).eq.1) then
      letdoit=0
      do idm=1,idmSet
        if(igMark(idm)==ig)letdoit=idm
      enddo
      if(letdoit>0)then
       allocate(sm(ismMax,sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1))
       do k=1,nzz-1
        i03=ip3mem(k,ne)
        i04=ip4mem(k,ig,ne)
        i02=ip2mem(ne)
        do ism=1,ismMax
          i04sm=ipSMmem(k,ism,letdoit,ne)
          do i=sx(ne)-1,ex(ne)+1
          do j=sy(ne)-1,ey(ne)+1
               ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
               sm(ism,i,j)= SourceMark(i04sm+ixy)
          enddo
          enddo
        enddo
        call adv_hori_mark(myid,gas(i04),u(i03),v(i03),&
            dx(i03),dy(i03),sx(ne), ex(ne), sy(ne),&
            ey(ne),ne,nx(ne),ny(ne),dt,k,ktop(i02),ismMax,sm)

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
    endif !ifsm
   !=================================

enddo loop_gas_spc    !ig
!========================


if(laerv2) then
   do ia=1,naersp
   do is=1,naerbin

      if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i04aer=ip4mem_aer(k,is,ia,ne)
            kpmass_m2(i03+ixy)=aerom(i04aer+ixy)
         enddo
      enddo
      enddo
      endif

      do k=1,nzz-1
       i04aer=ip4mem_aer(k,is,ia,ne)
       i03=ip3mem(k,ne)
         CALL ADV_HORI(MYID, aerom(i04aer),u(i03),v(i03),dx(i03),dy(i03) &
                      ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt)
      enddo ! k


       if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
          i03=ip3mem(k,ne)
          i04aer=ip4mem_aer(k,is,ia,ne)
        call  balance(myid,aerom(i04aer),kpmass_m2(i03),RatioMass(i03), &
                 sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
       endif

   enddo        !isize
   enddo

endif


   ! for DUST AND SEA SALT  aerosols
    do ia=1,iaer
    do is=1,isize

      if(imasskeep==1)then
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         do k=1,nzz-1
            i03=ip3mem(k,ne)
            i05=ip5mem(k,is,ia,ne)
            kpmass_m2(i03+ixy)=aer(i05+ixy)
         enddo
      enddo
      enddo
      endif

      do k=1,nzz-1
       i05=ip5mem(k,is,ia,ne)
       i03=ip3mem(k,ne)
         CALL ADV_HORI(MYID, AER(I05),u(i03),v(i03),dx(i03),dy(i03) &
                      ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt)
      enddo ! k

      do k = 1, nzz - 1

        IF(ia==1) THEN ! sea salt
         do iduc = 1, nseacom
            i05c = ip5memcs (k,is,iduc,ne)
            i03=ip3mem(k,ne)
          call ADV_HORI(MYID, SEACOMP(i05c),u(i03),v(i03),dx(i03),dy(i03) &
                       ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
         enddo
        ELSE  IF(ia==2) THEN ! for dust
         do iduc = 1, ndustcom
          i05c = ip5memc (k,is,iduc,ne)
          i03=ip3mem(k,ne)
          call ADV_HORI(MYID, DUSTCOMP(i05c),u(i03),v(i03),dx(i03),dy(i03) &
                       ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
         enddo ! iduc   

        ENDIF ! ia

       enddo ! k

       if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
          i03=ip3mem(k,ne)
          i05=ip5mem(k,is,ia,ne)
        call  balance(myid,aer(i05),kpmass_m2(i03),RatioMass(i03), &
                 sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
       endif

     enddo        !isize
    enddo




end subroutine naqpms_h_adv





