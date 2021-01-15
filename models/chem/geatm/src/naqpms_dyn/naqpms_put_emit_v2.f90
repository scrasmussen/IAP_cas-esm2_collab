
subroutine naqpms_put_emit_v2 &
 & ( myid &
 &  ,ne,dt,nx,ny,nzz,nest,sy,ey,sx,ex &
 &  ,mem2d &
 &  ,mem3d &
 &  ,GC_MOLWT &
 &  ,igas,iaer,isize,nseacom,ndustcom &
 &  ,igasCBM &
 &  ,NLAY_EM &
 &  ,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use naqpms_gridinfo, only: dz
use met_fields, only: t,Plev

use smpsulf_var, only: idx4dep_smpsulf

implicit none

real,parameter :: wgt_h2so4=98,wgt_nh3=17,wgt_hno3=63,wgt_hcl=36.5
real,parameter :: wgt_so4=96,wgt_nh4=18,wgt_no3=62,wgt_na=23,wgt_cl=35.5

integer :: myid

real :: dt

integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: igas,iaer,isize,nseacom,ndustcom

integer :: igasCBM

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

real :: trfac

real :: fac2modes(2)
real :: fac2modes_ppm(2),fac2modes_bc(2),fac2modes_oc(2)

real :: frcbin

integer :: ig_emt


fac2modes(1)=1.0
fac2modes(2)=1.0-fac2modes(1)


loop_emtyp :  do iemittype=1,NLAY_EM

  do ia=1,naersp
  if(ia.gt.4.and.ia.le.12) cycle ! skip secondaery aerosols
  do is=1,naerbin ! == 2 (two modes)

    do k=1,nzz-1

      i04=ip4mem_aer(k,is,ia,ne)

      i03=ip3mem(k,ne)

      if(ia.eq.1) then
        if(is.eq.1) ig=75
        if(is.eq.2) ig=76
        frcbin=1.0
      elseif(ia.eq.2) then
        ig=77
        frcbin=fac2modes(is)
      elseif(ia.eq.3) then
        ig=78
        frcbin=fac2modes(is)
      elseif(ia.eq.4) then
        ig=92
        frcbin=fac2modes(is)!*0.0
      elseif(ia.eq.13) then  ! hydrophobic bc
        ig=77
        frcbin=fac2modes(is)*1.0
      elseif(ia.eq.14) then  ! hydrophilic bc
        ig=77
        frcbin=fac2modes(is)*0.0
      endif


      i0emt=ip_emit2d(ig,iemittype,ne) ! for test

      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

        trfac=1.0

        deltc=frcbin*emit2d(i0emt+ixy)*emt2d_zfrc(k,iemittype)*dt/dz(i03+ixy)
        flag=1.0
        if(deltc.lt.0) then
          deltc=0.0
          print*,'emit-erro',i,j,k,is,emit2d(i0emt+ixy)*emt2d_zfrc(k,iemittype)
        endif
        aerom(i04+ixy)=aerom(i04+ixy)*trfac+deltc*flag
        aerom(i04+ixy)=aerom(i04+ixy)/trfac
      enddo
      enddo

    enddo

  enddo
  enddo

!  do ig=1,igas   ! for gas
   do ig=1,iedgas

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

      ig_emt=ig
      if(lgaschemsmp) then
        ig_emt=idx4dep_smpsulf(ig)  !18 ! SO2
      endif

      i0emt=ip_emit2d(ig_emt,iemittype,ne)

      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
        ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

        if(ig.le.iedgas) then ! ppb->ug/m3
          trfac=GC_MOLWT(ig_emt)/( (0.08206*T(i03+ixy))/(Plev(i03+ixy)/1000.) )
        else
          trfac=1.0
        endif

        deltc=emit2d(i0emt+ixy)*emt2d_zfrc(k,iemittype)*dt/dz(i03+ixy)
        flag=1.0
        if(deltc.lt.0) then
          deltc=0.0
          print*,'emit-erro',i,j,k,is,emit2d(i0emt+ixy)*emt2d_zfrc(k,iemittype)
        endif
        gas(i04+ixy)=gas(i04+ixy)*trfac+deltc*flag
        gas(i04+ixy)=gas(i04+ixy)/trfac
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

end subroutine naqpms_put_emit_v2
