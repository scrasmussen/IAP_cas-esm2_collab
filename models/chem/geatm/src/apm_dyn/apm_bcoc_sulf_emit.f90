
subroutine apm_bcoc_sulf_emit &
 & ( myid &
 &  ,lapm &
 &  ,ne,nx,ny,nzz,nest,sx,ex,sy,ey &
 &  ,igas &
 &  ,ip3mem &
 &  ,ip4mem &
 &  ,ip2memGas,mem2dgas &
 &  ,NLAY_EM,emt2d_zfrc,cfmode,ip_emit2d,mem_emt2d,emit2d )

!use naqpms_varlist, only:
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

integer :: NLAY_EM
integer :: ip_emit2d(igas,NLAY_EM,nest)
integer :: mem_emt2d
real    :: emt2d_zfrc(nzz,NLAY_EM)

real,dimension(mem_emt2d)           :: emit2d
character(len=*),dimension(NLAY_EM) :: cfmode


integer :: mem2dgas

integer :: flag(10)

integer :: i02emt,ictg



!print*,cfmode
!print*,'ne emit type:',ne,iemittype
!print*,'flag=',flag

!print*,'dim:',sx(ne),ex(ne),sy(ne),ey(ne)




do j=sy(ne),ey(ne)
do i=sx(ne),ex(ne)
  ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
  do k=1,nzz

    i03=ip3mem(k,ne)


!------------------------
    ig1=92 ! H2SO4(AQ)

    sulf_emit(i03+ixy) = 0.0
    do ictg=1,NLAY_EM
       i02emt=ip_emit2d(ig1,ictg,ne)
       sulf_emit(i03+ixy) = sulf_emit(i03+ixy) + &
                          & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
    enddo


!------------------------
    ig1=77 ! BC

    bc_emit(i03+ixy)   = 0.0
    do ictg=1,NLAY_EM
      i02emt=ip_emit2d(ig1,ictg,ne)
      bc_emit(i03+ixy) = bc_emit(i03+ixy) + &
                       & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
    enddo


    bc_emit_bb(i03+ixy)=0.0
    bc_emit_ff(i03+ixy)=0.0
    do ictg=1,NLAY_EM
      i02emt=ip_emit2d(ig1,ictg,ne)
      if(trim(cfmode(ictg)).eq.'bb') then
        bc_emit_bb(i03+ixy) = bc_emit_bb(i03+ixy) + &
                            & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
      elseif(trim(cfmode(ictg)).eq.'ff') then
        bc_emit_ff(i03+ixy) = bc_emit_ff(i03+ixy) + &
                            & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
      else
        stop 'cfmode err'
      endif
    enddo

!---------------------------
    ig1=78 ! OC

    oc_emit(i03+ixy) = 0.0
    do ictg=1,NLAY_EM
      i02emt=ip_emit2d(ig1,ictg,ne)
      oc_emit(i03+ixy) = oc_emit(i03+ixy) + &
                       & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
    enddo

    oc_emit_bb(i03+ixy)=0.0
    oc_emit_ff(i03+ixy)=0.0
    do ictg=1,NLAY_EM
      i02emt=ip_emit2d(ig1,ictg,ne)
      if(trim(cfmode(ictg)).eq.'bb') then
        oc_emit_bb(i03+ixy) = oc_emit_bb(i03+ixy) + &
                            & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
      elseif(trim(cfmode(ictg)).eq.'ff') then
        oc_emit_ff(i03+ixy) = oc_emit_ff(i03+ixy) + &
                            & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
      else
        stop 'cfmode err'
      endif
    enddo


!------------------------
    ig1=75 ! ppm2.5

    ppmfine_emit(i03+ixy)=0.0
    do ictg=1,NLAY_EM
      i02emt=ip_emit2d(ig1,ictg,ne)
      ppmfine_emit(i03+ixy) = ppmfine_emit(i03+ixy) + &
                            & emit2d(i02emt+ixy)*emt2d_zfrc(k,ictg)
    enddo

  enddo
enddo
enddo

end subroutine apm_bcoc_sulf_emit


