

SUBROUTINE apm_eff_rad 

!
use apm_varlist
implicit none
include 'apm_parm.inc'
integer :: ibin
! in


!IF(lapm) THEN ! apm flag

 loop_j : DO j = sy(ne),ey(ne)
 loop_i : DO i = sx(ne),ex(ne)

   i02 = ip2mem(ne)
   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

   loop_k : DO k=1,nzz

     i03 = ip3mem(k,ne)
     RH = rh1(i03+ixy)

     MNIT = 0
     MNH4 = 0
     MMSA = 0

     MSULFLV=0;MBCLV=0;MOCLV=0;MDSTLV=0;MSALTLV=0

     MBCS=mbcsulf(i03+ixy)
     MOCS=mocsulf(i03+ixy)
     MDSTS=mdstsulf(i03+ixy)
     MSALTS=msltsulf(i03+ixy)


     loop1_sulf : do is=1,NSO4
       iapm=ip_sulf(k,is,ne)
       rd1d_sulf(is)=rd_sulf(iapm+ixy)
       XM1D(is)=apm_sulf(iapm+ixy)
     enddo loop1_sulf

     loop1_salt : do is=1,NSEA
       iapm=ip_salt(k,is,ne)
       rd1d_salt(is)=rd_salt(iapm+ixy)
       XM1D(is+NSO4)=apm_salt(iapm+ixy)
     enddo loop_salt
 
     loop1_dust : do is=1,NSEA
       iapm=ip_dust(k,is,ne)
       rd1d_dust(is)=rd_dust(iapm+ixy)
       XMDST(is)=apm_dust(iapm+ixy)
     enddo loop1_dust

     loop1_bcoc : do is=1,NBCOCT
       iapm=ip_bcoc(k,is,ne)
       MBCOC8(is)=apm_bcoc(iapm+ixy)
     enddo loop1_bcoc


     ! box model
     call eff_rad &
       & ( rd1d_sulf,rd1d_salt,rd1d_dust &
       &  ,RH,grfac &
       &  ,XM1D,XMDST,MBCOC8 &
       &  ,MNIT,MNH4,MMSA &
       &  ,MSULFLV,MBCLV,MOCLV,MDSTLV,MSALTLV &
       &  ,MBCS,MOCS,MDSTS,MSALTS &
       &  ,r1d_sulf,r1d_salt,r1d_dust ) !&

     loop2_sulf : do is=1,NSO4
       iapm=ip_sulf(k,is,ne)
       r_sulf(iapm+ixy)=r1d_sulf(is)
     enddo loop2_sulf

     loop2_salt : do is=1,NSEA
       iapm=ip_salt(k,is,ne)
       r_salt(iapm+ixy)=r1d_salt(is)
     enddo loop2_dust

     loop2_dust : do is=1,NDSTB
       iapm=ip_dust(k,is,ne)
       r_dust(iapm+ixy)=r1d_dust(is)
     enddo loop2_dust



   enddo loop_k

 enddo loop_i
 enddo loop_k


!ENDIF ! apm flag


return


END SUBROUTINE apm_eff_rad
