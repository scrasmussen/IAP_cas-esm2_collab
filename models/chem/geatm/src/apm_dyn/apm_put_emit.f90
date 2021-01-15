
subroutine put_apm_emit &
 & ( myid &
 &  ,lapm &
 &  ,dt,ne &
 &  ,nx,ny,nzz,nest,sx,ex,sy,ey &
 &  ,dz &
 &  ,ip3mem,mem3d )

use apm_varlist
implicit none
include 'apm_parm.inc'
integer :: myid
logical :: lapm
real    :: dt
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)
 
integer :: i,j,k,is

integer :: ixy,i03,iapm

integer :: ip3mem(nzz,nest)

integer :: mem3d
real    :: dz(mem3d)

real    :: dz1d

real*8  :: m1p,radius

integer,parameter :: ridx(1:8) = (/1,2,1,2,1,2,1,2/)

! NOTE : sulf_emit, bcoc_emit are bulk type
!      : salt_emit has been scaled to each bins

!real,parameter :: const=1 ! dyn ok
real,parameter :: const=1

real,parameter :: frac95=1.00
real,parameter :: frac05=1-frac95

real :: emission_rate

real :: mass_frac


integer :: ssss



!=========================
!-------------------------
! add primary sulfate to apm sulfate particle
if(lem_sulf) then
!print*,'put sulf emit'
 loop_so4 : do is=1,NSO4
   do k=1,nzz-1
      i03=ip3mem(k,ne)
      iapm=ip_sulf(k,is,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)

         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

         if(ls2sulf) then ! add nucleation and accumulation mode
           emission_rate = sulf_emit(i03+ixy)*( frac95*sulf_unit_em2(is,2) + &
                         & frac05*sulf_unit_em2(is,1) )
         else             ! add nucleation mode
           emission_rate = sulf_emit(i03+ixy)*( frac05*sulf_unit_em2(is,1) )
         endif

!20140601@albany
 emission_rate=0.0

         dz1d=dz(i03+ixy)
         apm_sulf(iapm+ixy) = apm_sulf(iapm+ixy)+emission_rate*dt/dz1d

      enddo
      enddo
   enddo
 enddo loop_so4
endif
!-------------------------
!=========================


!=========================
!-------------------------

if(lem_salt) then
!print*,'put salt emit'
 loop_salt : do is=1,NSEA
   do k=1,nzz-1
      i03=ip3mem(k,ne)
      iapm=ip_salt(k,is,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         dz1d=dz(i03+ixy)

         apm_salt(iapm+ixy) = apm_salt(iapm+ixy)+ &
                            & dt*salt_emit(iapm+ixy)/dz1d
      enddo
      enddo
   enddo
 enddo loop_salt
endif
!-------------------------
!=========================

!return

!========================== 
!--------------------------
if(lem_dust) then
!print*,'put dust emit'
 loop_dust : do is=1,NDSTB
   do k=1,nzz-1
      i03=ip3mem(k,ne)
      iapm=ip_dust(k,is,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)
         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         dz1d=dz(i03+ixy)
         apm_dust(iapm+ixy) = apm_dust(iapm+ixy)+ &
                            & dust_emit(i03+ixy)*fdstbin(is)*dt/dz1d
      enddo
      enddo
   enddo
 enddo loop_dust
endif
!-------------------------
!=========================

!return

!===========================
!---------------------------
if(lem_bcoc) then

 if(lppmass_ch) then
    mass_frac = frc_mass_ch
    !print*,'shun_kk lppmass_ch=',lppmass_ch
    !print*,'frc_mass_ch=',frc_mass_ch
    !stop 'shun_kk'
 else
    mass_frac = 1.0
 endif


 !print*,'put bcoc emit'
 loop_bcoc : do is=1,NBCOCT ! 8
   !cycle
   do k=1,nzz-1
      i03=ip3mem(k,ne)
      iapm=ip_bcoc(k,is,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)

         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1
         dz1d=dz(i03+ixy)

!bcoc_unit_em2 : #/kg(emit)

! modify code below

!print*,'shunkkbbcoc',sum(bcoc_unit_em2(:,1))*dt/dz1d,sum(bcoc_unit_em2(:,2))*dt/dz1d
!print*,sum(bcoc_unit_em2(:,1)),sum(bcoc_unit_em2(:,2))
!stop


! apm_bcoc : ug/m3
         radius=ref1d_bcoc(ridx(is))*1.0e2       ! m->cm
         m1p=denbcoc*apm_Pi*4.0/3.0*radius**3.0*1.0e6  ! g->ug


         !apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+k
         !cycle

         if(is.eq.1) then     ! FF(mode1) QW BC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*bc_emit_ff(i03+ixy)*fbc_qw*dt/dz1d
         elseif(is.eq.2) then ! BB(mode2) QW BC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*bc_emit_bb(i03+ixy)*fbc_qw*dt/dz1d
         elseif(is.eq.3) then ! FF(mode1) QW OC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*oc_emit_ff(i03+ixy)*foc_qw*dt/dz1d
         elseif(is.eq.4) then ! BB(mode2) QW OC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*oc_emit_bb(i03+ixy)*foc_qw*dt/dz1d
         elseif(is.eq.5) then ! FF(mode1) SW BC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*bc_emit_ff(i03+ixy)*fbc_sw*dt/dz1d
         elseif(is.eq.6) then ! BB(mode2) SW BC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*bc_emit_bb(i03+ixy)*fbc_sw*dt/dz1d
         elseif(is.eq.7) then ! FF(mode1) SW OC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*oc_emit_ff(i03+ixy)*foc_sw*dt/dz1d
         elseif(is.eq.8) then ! BB(mode2) SW OC
           apm_bcoc(iapm+ixy) = apm_bcoc(iapm+ixy)+ &
            & mass_frac*oc_emit_bb(i03+ixy)*foc_sw*dt/dz1d
         endif

      enddo
      enddo
   enddo
 enddo loop_bcoc
endif
!---------------------------
!===========================




 loop_binbc : do is=1,nbincb 
   !cycle
   do k=1,nzz-1
      i03=ip3mem(k,ne)
      iapm=ip_cbbin(k,is,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)

         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

         dz1d=dz(i03+ixy)

         apm_binbc(iapm+ixy)=apm_binbc(iapm+ixy) &
                          & +dt*bc_emit_ff(i03+ixy)*bc_emit2bin_mfrc(is,1)/dz1d &
                          & +dt*bc_emit_bb(i03+ixy)*bc_emit2bin_mfrc(is,2)/dz1d

      enddo
      enddo
   enddo
 enddo loop_binbc


 loop_binoc : do is=1,nbincb
   !cycle
   do k=1,nzz-1
      i03=ip3mem(k,ne)
      iapm=ip_cbbin(k,is,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)

         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

         dz1d=dz(i03+ixy)

         apm_binoc(iapm+ixy)=apm_binoc(iapm+ixy) &
                          & +dt*oc_emit_ff(i03+ixy)*oc_emit2bin_mfrc(is,1)/dz1d &
                          & +dt*oc_emit_bb(i03+ixy)*oc_emit2bin_mfrc(is,2)/dz1d

      enddo
      enddo
   enddo
 enddo loop_binoc


if(1==2) then
do is=1,nbincb
write(*,'(i2,3(2x,f15.12))') is,rd_binoc(is),oc_emit2bin_mfrc(is,1),oc_emit2bin_mfrc(is,2)
write(*,'(i2,3(2x,f15.12))') is,rd_binbc(is),bc_emit2bin_mfrc(is,1),bc_emit2bin_mfrc(is,2)
enddo

print*,sum(bc_emit2bin_mfrc(:,1))
print*,sum(bc_emit2bin_mfrc(:,2))
print*,sum(oc_emit2bin_mfrc(:,1))
print*,sum(oc_emit2bin_mfrc(:,2))

stop

endif
















!===========================
!---------------------------
if(ls2bcoc) then ! add primary sulfate to bcoc coated sulfate mass
   do k=1,nzz-1
      i03=ip3mem(k,ne)
      do j = sy(ne),ey(ne)
      do i = sx(ne),ex(ne)

         ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

         dz1d=dz(i03+ixy)

         emission_rate=sulf_emit(i03+ixy)*frac95*fs2bc
         mbcsulf(i03+ixy)=mbcsulf(i03+ixy)+emission_rate*dt/dz1d

         emission_rate=sulf_emit(i03+ixy)*frac95*fs2oc
         mocsulf(i03+ixy)=mocsulf(i03+ixy)+emission_rate*dt/dz1d

      enddo
      enddo
   enddo
endif
!---------------------------
!===========================



return
end subroutine put_apm_emit





