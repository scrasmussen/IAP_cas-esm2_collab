

subroutine get_sp_distribution(myid,i,j,k,ne,ixy,msulf_sp)
use apm_varlist
implicit none
include 'apm_parm.inc'
integer :: myid
integer :: is,i,j,k
integer :: ne
integer :: iapm,ixy
real    :: msulf_sp
real    :: wfac


wfac=0.5

do is=1,NSO4
 iapm=ip_sulf(k,is,ne)
 ! kg/m3 -> ug/m3
 apm_sulf(iapm+ixy)=msulf_sp*kg2ug*(sulf_unit_em2(is,1)*wfac+(1.0-wfac)*sulf_unit_em2(is,2))
enddo

end subroutine get_sp_distribution



subroutine get_coated_sulfate( myid,i,j,k,ne,ixy &
                              ,msulf_pp,MSALTS,MDSTS,MBCS,MOCS)
use apm_varlist
implicit none
include 'apm_parm.inc'
integer :: myid
integer :: ne,nest

real    :: msulf_pp ! kg/m3
real*8  :: MSALTS,MDSTS,MBCS,MOCS ! kg/m3
real*8  :: radius
real*8  :: totarea,area_salt(NSEA),area_dust(NDSTB),area_bcoc(NBCOCT)
integer :: is,i,j,k
integer :: iapm,ixy

integer,parameter :: ridx(1:8) = (/1,2,1,2,1,2,1,2/)
character*2,parameter :: flagbcoc(1:8) = (/'bc','bc','oc','oc' &
                                          ,'bc','bc','oc','oc'/)

totarea=0.0d0

do is=1,NSEA
 radius=rd_salt(is)
 iapm=ip_salt(k,is,ne)
 area_salt(is)=apm_salt(iapm+ixy)/ &
        & (densalt*1.0e+12/3.0d0*radius)
enddo

do is=1,NDSTB
 radius=rd_dust(is)
 iapm=ip_dust(k,is,ne)
 area_dust(is)=apm_dust(iapm+ixy)/ &
        & (dendust*1.0e+12/3.0d0*radius)
enddo

do is=1,NBCOCT
 radius=ref1d_bcoc(ridx(is))
 iapm=ip_bcoc(k,is,ne)
 area_bcoc(is)=apm_bcoc(iapm+ixy)/ &
        & (denbcoc*1.0e+12/3.0d0*radius)
! print*,'kkbocarea',i,j,k
! print*,is,radius
! print*,apm_bcoc(iapm+ixy),area_bcoc(is)
enddo

!stop

totarea=sum(area_salt(:))+sum(area_dust(:))+sum(area_bcoc(:))

if(totarea.gt.0) then
 MSALTS=msulf_pp*sum(area_salt(:))/totarea
 MDSTS =msulf_pp*sum(area_dust(:))/totarea
 MBCS  =msulf_pp*(area_bcoc(1)+area_bcoc(2)+area_bcoc(5)+area_bcoc(6))/totarea
 MOCS  =msulf_pp*(area_bcoc(3)+area_bcoc(4)+area_bcoc(7)+area_bcoc(8))/totarea
else
 MSALTS=0
 MDSTS =0
 MBCS  =0
 MOCS  =0
endif

end subroutine get_coated_sulfate






