
subroutine naqpms_put_emit &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,GC_MOLWT &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,NLAY_EM &
 &  ,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use naqpms_gridinfo, only : dz

implicit none

real,parameter :: wgt_h2so4=98,wgt_nh3=17,wgt_hno3=63,wgt_hcl=36.5
real,parameter :: wgt_so4=96,wgt_nh4=18,wgt_no3=62,wgt_na=23,wgt_cl=35.5

integer :: myid

real :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: NLAY_EM

integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)


integer :: letdoit,MapS
real :: DeltSpeMark,OrgCon

integer :: i,j,k,is

integer :: mem3d
integer :: mem2d

real,dimension(igas) :: GC_MOLWT


integer :: ixy,i02,i03,i04
integer :: i02Gas
integer :: i0emt

integer :: ig,iduc,ia,i05,i05c

integer :: idm,ism,i04sm

integer :: iemittype


real :: contem0(idmSet)
real :: TmpSM(ismMax)

integer :: IHgtLev,iHgtLMax,ISrcDefined

real :: deltc,flag


loop_emtyp :  do iemittype=1,NLAY_EM


  do ig=1,igas   ! for gas

    i02Gas= ip2memGas(ig,ne)

    do k=1,nzz-1     

      i03=ip3mem(k,ne)
      i04=ip4mem(k,ig,ne)

      !!!!!!!!!!!!!!!!
      ! for Source Mark
      if(ifsm(ne)==1) then
        letdoit=0
        do idm=1,idmSet
          if(igMark(idm)==ig)letdoit=idm
        enddo
        if(letdoit>0)then
          i02=ip2mem(ne)
          do j = sy(ne),ey(ne)
          do i = sx(ne),ex(ne)
             ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1 
             tmpMarkCon(i02+ixy)=gas(i04+ixy)
          enddo
          enddo
        endif
      endif !ifsm
      !!!!!!!!!!!!!!!!

      i0emt=ip_emit2d(ig,iemittype,ne)

      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
        deltc=emit2d(i0emt+ixy)*emt2d_zfrc(k,iemittype)*dt/dz(i03+ixy)
        flag=1.0
        if(deltc.lt.0) then
          deltc=0.0
          print*,'emit-erro',i,j,k,is,emit2d(i0emt+ixy)*emt2d_zfrc(k,iemittype)
        endif
        gas(i04+ixy)=gas(i04+ixy)+deltc*flag
      enddo
      enddo


     !!!!!!!!!!!!!!!!
     ! for Source Mark
     if(ifsm(ne)==1)then 
     if(letdoit>0)then
      i02=ip2mem(ne)
      i04=ip4mem(k,ig,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1 
        OrgCon = tmpMarkCon(i02+ixy)
        DeltSpeMark=gas(i04+ixy)-OrgCon
        if(DeltSpeMark > 0 )then 
             do ism=1,ismMax
               i04sm=ipSMmem(k,ism,letdoit,ne)
               TmpSM(ism)=SourceMark(i04sm+ixy)
             enddo

             MapS=int(MapSource(i02+ixy))

             IHgtLev=iemittype-1
             call GetSMarkChange(DeltSpeMark,OrgCon,TmpSM,MapS, &
                  ISrcDefined,IHgtLev,ismMax)

             do ism=1,ismMax
               i04sm=ipSMmem(k,ism,letdoit,ne)
               SourceMark(i04sm+ixy)=TmpSM(ism)
             enddo

        endif
     enddo
     enddo
    endif
    endif !ifsm
    !!!!!!!!!!!!!!!!
   enddo ! k
 enddo ! ig

enddo loop_emtyp ! iemittype Zifa 2007/03/03

end subroutine naqpms_put_emit
