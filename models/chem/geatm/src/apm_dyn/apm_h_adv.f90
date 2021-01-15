
subroutine apm_h_adv &
 & ( myid &
 &  ,lapm,imasskeep &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,dx,dy &
 &  ,u,v &
 &  ,ip3mem,mem3d,RatioMass )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid

real :: dt

logical :: lapm
integer :: imasskeep

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: i,j,k,is

integer :: mem3d

real,dimension(mem3d) :: dx,dy,u,v,RatioMass,kpmass_m2

integer :: ixy,i03,iapm

integer :: ip3mem(nzz,nest)

real,dimension(mem3d) :: wk

!return

!IF(lapm) THEN ! apm flag

 !> shun : apm sulfate hadv
  if(lfor_sulf) then
   !print*,'sulf hadv'
   loop_so4_hadv : do is=1,NSO4
!exit loop_so4_hadv
     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_sulf(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_sulf(iapm+ixy)
          endif
          wk(i03+ixy)=apm_sulf(iapm+ixy)
        enddo
     enddo
     enddo
!goto 123
     do k=1,nzz-1
       i03=ip3mem(k,ne)
       CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo
!goto 123
     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

!123 continue

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_sulf(k,is,ne)
          apm_sulf(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_so4_hadv
  endif
 !< shun : end of apm sulfate hadv




 !> shun : apm seasalt hadv
  if(lfor_salt) then
   !print*,'salt hadv'
   loop_salt_hadv : do is=1,NSEA

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_salt(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_salt(iapm+ixy)
          endif
          wk(i03+ixy)=apm_salt(iapm+ixy)
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_salt(k,is,ne)
          apm_salt(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_salt_hadv
  endif
 !< shun : end of apm seasalt hadv

 !> shun : apm dust hadv
  if(lfor_dust) then
   !print*,'dust hadv'
   loop_dust_hadv : do is=1,NDSTB

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_dust(k,is,ne)
          wk(i03+ixy)=apm_dust(iapm+ixy)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_dust(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_dust(k,is,ne)
          apm_dust(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_dust_hadv
  endif

 !< shun : end of apm dust hadv

 !> shun : apm bcoc hadv
 if(lfor_bcoc) then
   !print*,'bcoc hadv'
   loop_bcoc_hadv : do is=1,NBCOCT

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          wk(i03+ixy)=apm_bcoc(iapm+ixy)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_bcoc(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

     do j = sy(ne)-1,ey(ne)+1
     do i = sx(ne)-1,ex(ne)+1
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          apm_bcoc(iapm+ixy)=wk(i03+ixy)
        enddo
     enddo
     enddo

   enddo loop_bcoc_hadv
  endif
 !< shun : end of apm bcoc hadv


!===============================================================
!===============================================================

!return


! apm coated species
IF(lcoated_dyn) THEN
    !print*,'coated hadv'
!-> sulfate coated on seasalt
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=msltsulf(i03+ixy)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=msltsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         msltsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo

!-> sulfate coated on dust

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=mdstsulf(i03+ixy)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mdstsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         mdstsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo


!-> sulfate coated on BC

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=mbcsulf(i03+ixy)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mbcsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         mbcsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo

!-> sulfate coated on POC

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=mocsulf(i03+ixy)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mocsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         mocsulf(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo
ENDIF ! lcoated_dyn


! sulferic acid vapor
  if(lfor_h2so4) then
    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         wk(i03+ixy)=h2so4_gas(i03+ixy)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=h2so4_gas(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,wk(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,wk(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif

    do j = sy(ne)-1,ey(ne)+1
    do i = sx(ne)-1,ex(ne)+1
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         h2so4_gas(i03+ixy)=wk(i03+ixy)
       enddo
    enddo
    enddo 
  endif

!ENDIF ! apm flag


end subroutine apm_h_adv





