
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


!IF(lapm) THEN ! apm flag

 !> shun : apm sulfate hadv
  if(lfor_sulf) then
   !print*,'sulf hadv'
   loop_so4_hadv : do is=1,NSO4
!exit loop_so4_hadv
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_sulf(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_sulf(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       iapm=ip_sulf(k,is,ne)
       CALL ADV_HORI( MYID,apm_sulf(iapm),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation
         i03=ip3mem(k,ne) 
         iapm=ip_sulf(k,is,ne)
         call balance(myid,apm_sulf(iapm),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

   enddo loop_so4_hadv
  endif
 !< shun : end of apm sulfate hadv




 !> shun : apm seasalt hadv
  if(lfor_salt) then
   !print*,'salt hadv'
   loop_salt_hadv : do is=1,NSEA

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_salt(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_salt(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       iapm=ip_salt(k,is,ne)
       CALL ADV_HORI( MYID,apm_salt(iapm),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         iapm=ip_salt(k,is,ne)
         call balance(myid,apm_salt(iapm),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

   enddo loop_salt_hadv
  endif
 !< shun : end of apm seasalt hadv

 !> shun : apm dust hadv
  if(lfor_dust) then
   !print*,'dust hadv'
   loop_dust_hadv : do is=1,NDSTB

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_dust(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_dust(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       iapm=ip_dust(k,is,ne)
       CALL ADV_HORI( MYID,apm_dust(iapm),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         iapm=ip_dust(k,is,ne)
         call balance(myid,apm_dust(iapm),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif


   enddo loop_dust_hadv
  endif

 !< shun : end of apm dust hadv

 !> shun : apm bcoc hadv
 if(lfor_bcoc) then
   !print*,'bcoc hadv'
   loop_bcoc_hadv : do is=1,NBCOCT

     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_bcoc(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_bcoc(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       iapm=ip_bcoc(k,is,ne)
       CALL ADV_HORI( MYID,apm_bcoc(iapm),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         iapm=ip_bcoc(k,is,ne)
         call balance(myid,apm_bcoc(iapm),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif


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
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=msltsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,msltsulf(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,msltsulf(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif


!-> sulfate coated on dust

    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mdstsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,mdstsulf(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,mdstsulf(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif


!-> sulfate coated on BC

    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mbcsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,mbcsulf(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,mbcsulf(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif


!-> sulfate coated on POC

    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=mocsulf(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,mocsulf(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,mocsulf(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif

ENDIF ! lcoated_dyn


! sulferic acid vapor
  if(lfor_h2so4) then
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
       ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
       do k=1,nzz-1
         i03=ip3mem(k,ne)
         if(imasskeep==1) then
           kpmass_m2(i03+ixy)=h2so4_gas(i03+ixy)
         endif
       enddo
    enddo
    enddo

    do k=1,nzz-1
      i03=ip3mem(k,ne)
      CALL ADV_HORI( MYID,h2so4_gas(i03),u(i03),v(i03),dx(i03),dy(i03) &
           & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
    enddo

    if(imasskeep==1)then
      do k=1,nzz-1   !mass conservation 
        i03=ip3mem(k,ne)
        call balance(myid,h2so4_gas(i03),kpmass_m2(i03),RatioMass(i03) &
             &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
      enddo
    endif

  endif

!ENDIF ! apm flag

   loop_binbc : do is=1,nbincb
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_cbbin(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_binbc(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       iapm=ip_cbbin(k,is,ne)
       CALL ADV_HORI( MYID,apm_binbc(iapm),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         iapm=ip_cbbin(k,is,ne)
         call balance(myid,apm_binbc(iapm),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

   enddo loop_binbc


   loop_binoc : do is=1,nbincb
     do j = sy(ne),ey(ne)
     do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        do k=1,nzz-1
          i03=ip3mem(k,ne)
          iapm=ip_cbbin(k,is,ne)
          if(imasskeep==1) then
            kpmass_m2(i03+ixy)=apm_binoc(iapm+ixy)
          endif
        enddo
     enddo
     enddo

     do k=1,nzz-1
       i03=ip3mem(k,ne)
       iapm=ip_cbbin(k,is,ne)
       CALL ADV_HORI( MYID,apm_binoc(iapm),u(i03),v(i03),dx(i03),dy(i03) &
            & ,sx(ne),ex(ne),sy(ne),ey(ne),nx(ne),ny(ne),dt )
     enddo

     if(imasskeep==1)then
       do k=1,nzz-1   !mass conservation 
         i03=ip3mem(k,ne)
         iapm=ip_cbbin(k,is,ne)
         call balance(myid,apm_binoc(iapm),kpmass_m2(i03),RatioMass(i03) &
              &  ,sx(ne), ex(ne), sy(ne), ey(ne),nx(ne),ny(ne))
       enddo
     endif

   enddo loop_binoc


end subroutine apm_h_adv





