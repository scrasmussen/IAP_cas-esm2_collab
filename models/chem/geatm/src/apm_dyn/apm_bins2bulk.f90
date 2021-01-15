

 if(allocated(con_s)) deallocate(con_s)

! if(allocated(apm_wdeps_salt))   deallocate(apm_wdeps_salt)
! if(allocated(apm_wdep2s_salt))  deallocate(apm_wdep2s_salt)
! if(allocated(apm_wdflds_salt))  deallocate(apm_wdflds_salt)
! if(allocated(apm_wdfld2s_salt)) deallocate(apm_wdfld2s_salt)

! if(allocated(apm_wdeps_dust))   deallocate(apm_wdeps_dust)
! if(allocated(apm_wdep2s_dust))  deallocate(apm_wdep2s_dust)
! if(allocated(apm_wdflds_dust))  deallocate(apm_wdflds_dust)
! if(allocated(apm_wdfld2s_dust)) deallocate(apm_wdfld2s_dust)

! if(allocated(apm_wdeps_bcoc))   deallocate(apm_wdeps_bcoc)
! if(allocated(apm_wdep2s_bcoc))  deallocate(apm_wdep2s_bcoc)
! if(allocated(apm_wdflds_bcoc))  deallocate(apm_wdflds_bcoc)
! if(allocated(apm_wdfld2s_bcoc)) deallocate(apm_wdfld2s_bcoc)


 !==============!
 ! bins to bulk !
 !==============!

 DO k=1,nzz

 i03=ip3mem(k,ne)

 !-> salt
 if(k.eq.1) print*,'salt so4 bins2bulk'
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    msltsulf(i03+ixy)=0.0d0
    do is=1,NSEA
      iapm=ip_salt(k,is,ne)
      msltsulf(i03+ixy)=msltsulf(i03+ixy)+apm_msltsulf(iapm+ixy)
    enddo
 enddo
 enddo

 !-> dust
 if(k.eq.1) print*,'salt so4 bins2bulk'
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    mdstsulf(i03+ixy)=0.0d0
    do is=1,NDSTB
      iapm=ip_dust(k,is,ne)
      mdstsulf(i03+ixy)=mdstsulf(i03+ixy)+apm_mdstsulf(iapm+ixy)
    enddo
 enddo
 enddo

 !-> BCOC 
 if(k.eq.1) print*,'salt so4 bins2bulk'
 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)
    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    mbcsulf(i03+ixy)=0.0d0
    mocsulf(i03+ixy)=0.0d0
    do is=1,NBCOCT
      iapm=ip_bcoc(k,is,ne)
      if(flagbcoc(is).eq.'bc') then
        mbcsulf(i03+ixy)=mbcsulf(i03+ixy)+apm_mbcocsulf(iapm+ixy)
      elseif(flagbcoc(is).eq.'oc') then
        mocsulf(i03+ixy)=mocsulf(i03+ixy)+apm_mbcocsulf(iapm+ixy)
      else
        stop 'erro'
      endif
    enddo
 enddo
 enddo

 ENDDO ! k


 if(allocated(apm_msltsulf)) then
    deallocate(apm_msltsulf)
 else 
    stop 'deallocate(apm_msltsulf) erro'
 endif

 if(allocated(apm_mdstsulf)) then
    deallocate(apm_mdstsulf)
 else
    stop 'deallocate(apm_mdstsulf) erro'
 endif


 if(allocated(apm_mbcocsulf)) then
    deallocate(apm_mbcocsulf)
 else
    stop 'deallocate(apm_mbcocsulf) erro'
 endif


! coated deposited mass

 do j = sy(ne),ey(ne)
 do i = sx(ne),ex(ne)

    ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
    i02 = ip2mem(ne)

    wdep_sltsulf(i02+ixy)=0
    wdep2_sltsulf(i02+ixy)=0
    wdep_dstsulf(i02+ixy)=0
    wdep2_dstsulf(i02+ixy)=0
    wdep_bcsulf(i02+ixy)=0
    wdep2_bcsulf(i02+ixy)=0
    wdep_ocsulf(i02+ixy)=0
    wdep2_ocsulf(i02+ixy)=0

    do is=1,NSEA
      iapm=ip2_salt(is,ne)
      wdep_sltsulf(i02+ixy)=wdep_sltsulf(i02+ixy)+apm_wdeps_salt(iapm+ixy)
      wdep2_sltsulf(i02+ixy)=wdep2_sltsulf(i02+ixy)+apm_wdep2s_salt(iapm+ixy)
    enddo 

    do is=1,NDSTB
      iapm=ip2_dust(is,ne)
      wdep_dstsulf(i02+ixy)=wdep_dstsulf(i02+ixy)+apm_wdeps_dust(iapm+ixy)
      wdep2_dstsulf(i02+ixy)=wdep2_dstsulf(i02+ixy)+apm_wdep2s_dust(iapm+ixy)
    enddo

    do is=1,NBCOCT
      iapm=ip2_bcoc(is,ne)
      if(flagbcoc(is).eq.'bc') then
        wdep_bcsulf(i02+ixy)=wdep_bcsulf(i02+ixy)+apm_wdeps_bcoc(iapm+ixy)
        wdep2_bcsulf(i02+ixy)=wdep2_bcsulf(i02+ixy)+apm_wdep2s_bcoc(iapm+ixy)
      elseif(flagbcoc(is).eq.'oc') then
        wdep_ocsulf(i02+ixy)=wdep_ocsulf(i02+ixy)+apm_wdeps_bcoc(iapm+ixy)
        wdep2_ocsulf(i02+ixy)=wdep2_ocsulf(i02+ixy)+apm_wdep2s_bcoc(iapm+ixy)
      else
        stop 'erro'
      endif
    enddo

 enddo
 enddo



!ENDIF ! apm flag






