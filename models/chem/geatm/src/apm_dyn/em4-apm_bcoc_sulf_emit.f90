
subroutine apm_bcoc_sulf_emit &
 & ( myid &
 &  ,lapm &
 &  ,iemittype &
 &  ,ne,nx,ny,nzz,nest,sx,ex,sy,ey &
 &  ,igas &
 &  ,ip3mem &
 &  ,ip4mem &
 &  ,ip2memGas,mem2dgas &
 &  ,ratioemit,ratioemitP,ratioemitL,ratioemitB &
 &  ,EmtaGas,EmtpGas,EmttGas,EmtbGas )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid
logical :: lapm

integer :: iemittype
integer :: i,j,k,ig,ig1,ig2,i02Gas_1,i02Gas_2

integer :: ixy,i03

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: ip3mem(nzz,nest)
integer :: ip4mem(nzz,igas,nest)

integer :: igas
integer,dimension(igas,nest) :: ip2memGas

integer :: mem2dgas
real,dimension(mem2dgas) :: EmtaGas,EmtpGas,EmttGas,EmtbGas

real :: ratioemit(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas)
real :: ratioemitP(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas)
real :: ratioemitL(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas)
real :: ratioemitB(sx(ne)-1:ex(ne)+1,sy(ne)-1:ey(ne)+1,nzz,igas)

integer :: flag(4)

flag=0

! 1:anthoropogenic 2:power plant 3:biomass burning 4:biogenic

if    (iemittype.eq.1) then
 flag(1)=1
elseif(iemittype.eq.2) then
 flag(2)=1
elseif(iemittype.eq.3) then
 flag(3)=1
elseif(iemittype.eq.4) then
 flag(4)=1
else
 stop 'iemittype input erro'
endif


!print*
!print*,'ne emit type:',ne,iemittype
!print*,'flag=',flag

!print*,'dim:',sx(ne),ex(ne),sy(ne),ey(ne)

!print*,ip2memGas(1,ne)
!print*,'sub_ratioemit:',ratioemit(30,30,:,1)
!print*,'sub_EmtaGas:',EmtaGas(ip2memGas(1,ne):ip2memGas(1,ne)+500)



do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)
  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
  do k=1,nzz

    i03=ip3mem(k,ne)

!if(k.eq.1.and.i.eq.30.and.j.eq.30) then
!print*,'last call xy3030 bc '
!print*,k,j,i
!print*,oc_emit(i03+ixy)
!endif

    ig1=92 ! H2SO4(AQ)
    i02Gas_1= ip2memGas(ig1,ne)
    sulf_emit(i03+ixy) = sulf_emit(i03+ixy) + &
                       & EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)*flag(1) + &
                       & EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)*flag(2) + &
                       & EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)*flag(3) + &
                       & EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)*flag(4)

    !sulf_emit(i03+ixy) = sulf_emit(i03+ixy)+1*flag(1)+10*flag(2)+100*flag(3)+1000*flag(4)




if(ne.eq.1.and.k.eq.1.and.i.eq.30.and.j.eq.30.and..false.) then
!print*,'xy3030 so4'
print*,k,j,i
print*,'sulf emit'
print*,sulf_emit(i03+ixy)
print*,ig1,i02Gas_1+ixy
print*,EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)
print*,EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)
print*,EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)
print*,EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)
print*,EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)+ &
     & EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)+ &
     & EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)+ &
     & EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)
endif
!cycle 


    ig1=77 ! BC
    i02Gas_1= ip2memGas(ig1,ne)
    bc_emit(i03+ixy)   = bc_emit(i03+ixy) + &
                       & EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)*flag(1) + &
                       & EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)*flag(2) + &
                       & EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)*flag(3) + &
                       & EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)*flag(4)


    bc_emit_bb(i03+ixy)= bc_emit_bb(i03+ixy) + &
                       & EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)*flag(3)

    bc_emit_ff(i03+ixy)= bc_emit_ff(i03+ixy) + &
                       & EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)*flag(1) + &
                       & EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)*flag(2) + &
                       & EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)*flag(4)



if(ne.eq.1.and.k.eq.1.and.i.eq.30.and.j.eq.30.and..false.) then
!print*,'xy3030 bc'
print*,k,j,i
print*,'bc_emit:',bc_emit(i03+ixy)
!print*,ig1,i02Gas_1
print*,i02Gas_1+ixy
print*,EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)
print*,EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)
print*,EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)
print*,EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)
endif
!cycle

    ig1=78 ! OC
    i02Gas_1= ip2memGas(ig1,ne)
    oc_emit(i03+ixy)   = oc_emit(i03+ixy) + &
                       & EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)*flag(1) + &
                       & EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)*flag(2) + &
                       & EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)*flag(3) + &
                       & EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)*flag(4)

    oc_emit_bb(i03+ixy)= oc_emit_bb(i03+ixy) + &
                       & EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)*flag(3)

    oc_emit_ff(i03+ixy)= oc_emit_ff(i03+ixy) + &
                       & EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)*flag(1) + &
                       & EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)*flag(2) + &
                       & EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)*flag(4)



if(ne.eq.1.and.k.eq.1.and.i.eq.30.and.j.eq.30.and..false.) then
!print*,'xy3030 oc'
print*,k,j,i
print*,'oc_emit',oc_emit(i03+ixy)
!print*,ig1,i02Gas_1
print*,i02Gas_1+ixy
print*,EmtaGas(i02Gas_1+ixy)*ratioemit (i,j,k,ig1)
print*,EmtpGas(i02Gas_1+ixy)*ratioemitP(i,j,k,ig1)
print*,EmttGas(i02Gas_1+ixy)*ratioemitL(i,j,k,ig1)
print*,EmtbGas(i02Gas_1+ixy)*ratioemitB(i,j,k,ig1)
endif

  if(k.eq.1) then
    !print*,'kkem',sulf_emit(i03+ixy),bc_emit(i03+ixy),oc_emit(i03+ixy)
  endif

  enddo
enddo
enddo

!sulf_emit,bcoc_emit
end subroutine apm_bcoc_sulf_emit


