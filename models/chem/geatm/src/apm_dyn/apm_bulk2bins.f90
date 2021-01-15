

!>------------------------------------------
! allocate memory for coated sulfate in bins

 allocate(con_s(nzz))

! allocate(apm_wdeps_salt(salt2dmem))
! allocate(apm_wdep2s_salt(salt2dmem))
! allocate(apm_wdflds_salt(salt2dmem))
! allocate(apm_wdfld2s_salt(salt2dmem))

! allocate(apm_wdeps_dust(dust2dmem))
! allocate(apm_wdep2s_dust(dust2dmem))
! allocate(apm_wdflds_dust(dust2dmem))
! allocate(apm_wdfld2s_dust(dust2dmem))

! allocate(apm_wdeps_bcoc(bcoc2dmem))
! allocate(apm_wdep2s_bcoc(bcoc2dmem))
! allocate(apm_wdflds_bcoc(bcoc2dmem))
! allocate(apm_wdfld2s_bcoc(bcoc2dmem))

!=============================================


! IF(lapm) THEN
   ! distribute bulk seasalt coated sulfate to bins
   if(.not.allocated(apm_msltsulf)) then
      allocate(apm_msltsulf(saltmem))
   else
      stop 'allocate(apm_msltsulf) erro'
   endif
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius1))) then
      allocate(numcc(NSEA))
      allocate(tmp(NSEA))
      allocate(radius1(NSEA))
   else
      stop 'allocate(numcc) erro'
   endif

   DO k=1,nzz
    if(k.eq.1) print*,'salt so4 bulk2bins'
    i03=ip3mem(k,ne)
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NSEA
        iapm=ip_salt(k,is,ne)
        radius1(is)=rd_salt(is)
        numcc(is)= apm_salt(iapm+ixy)/ &
                 & radius1(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius1(is)**3.0d0)
        tmp(is)=numcc(is)*radius1(is)**2
      enddo
      do is=1,NSEA
        iapm=ip_salt(k,is,ne)
        if(sum(tmp(:)).gt.0) then
          apm_msltsulf(iapm+ixy)=msltsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
          apm_msltsulf(iapm+ixy)=0
        endif
      enddo
    enddo
    enddo
   ENDDO
   if(allocated(numcc).and.allocated(tmp).and.allocated(radius1)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius1)
   else
      stop 'deallocate(numcc) erro'
   endif
   !!!!!

   ! distribute bulk dust coated sulfate to bins
   if(.not.allocated(apm_mdstsulf)) then
      allocate(apm_mdstsulf(dustmem))
   else
      stop 'allocate(apm_mdstsulf) erro'
   endif
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius1))) then
      allocate(numcc(NDSTB))
      allocate(tmp(NDSTB))
      allocate(radius1(NDSTB))
   else
      stop 'allocate(numcc) erro'
   endif

   DO k=1,nzz
    if(k.eq.1) print*,'dust so4 bulk2bins'
    i03=ip3mem(k,ne)
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NDSTB
        iapm=ip_dust(k,is,ne)
        radius1(is)=rd_dust(is)
        numcc(is)= apm_dust(iapm+ixy)/ &
                 & radius1(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius1(is)**3.0d0)
        tmp(is)=numcc(is)*radius1(is)**2
      enddo
      do is=1,NDSTB
        iapm=ip_dust(k,is,ne)
        if(sum(tmp(:)).gt.0) then
          apm_mdstsulf(iapm+ixy)=mdstsulf(i03+ixy)*tmp(is)/sum(tmp(:))
        else
          apm_mdstsulf(iapm+ixy)=0
        endif
      enddo
    enddo
    enddo
   ENDDO
   if(allocated(numcc).and.allocated(tmp).and.allocated(radius1)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius1)
   else
      stop 'deallocate(numcc) erro'
   endif
   !!!!!

   ! distribute bulk BC coated sulfate to bins
   if(.not.allocated(apm_mbcocsulf)) then
      allocate(apm_mbcocsulf(bcocmem))
   else
      stop 'allocate(apm_mbcocsulf) erro'
   endif
   if(.not.(allocated(numcc).and.allocated(tmp).and.allocated(radius1))) then
      allocate(numcc(NBCOCT))
      allocate(tmp(NBCOCT))
      allocate(radius1(NBCOCT))
   else
      stop 'allocate(numcc) erro'
   endif
   DO k=1,nzz
    if(k.eq.1) print*,'bcoc so4 bulk2bins'
    i03=ip3mem(k,ne)
    do j = sy(ne),ey(ne)
    do i = sx(ne),ex(ne)
      ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        radius1(is)=ref1d_bcoc(ridx(is))
        numcc(is)= apm_bcoc(iapm+ixy)/ &
                 & radius1(is)**3.0d0
                 !& (densalt*1.0e+12*4.0d0/3.0d0*apm_pi*radius1(is)**3.0d0)
        tmp(is)=numcc(is)*radius1(is)**2
      enddo
      totalbc=tmp(1)+tmp(5)+tmp(2)+tmp(6)
      totaloc=tmp(3)+tmp(7)+tmp(4)+tmp(8)
      do is=1,NBCOCT
        iapm=ip_bcoc(k,is,ne)
        if(flagbcoc(is).eq.'bc') then
          if(totalbc.gt.0) then
            apm_mbcocsulf(iapm+ixy)=mbcsulf(i03+ixy)*tmp(is)/totalbc
          else
            apm_mbcocsulf(iapm+ixy)=0
          endif
        elseif(flagbcoc(is).eq.'oc') then
          if(totaloc.gt.0) then
            apm_mbcocsulf(iapm+ixy)=mocsulf(i03+ixy)*tmp(is)/totaloc
          else
            apm_mbcocsulf(iapm+ixy)=0
          endif
        else
          stop 'erro'
        endif
      enddo
    enddo
    enddo
   ENDDO
   if(allocated(numcc).and.allocated(tmp).and.allocated(radius1)) then
      deallocate(numcc)
      deallocate(tmp)
      deallocate(radius1)
   else
      stop 'deallocate(numcc) erro'
   endif
   !!!!!

! ENDIF 


