

subroutine apm_isor_so4 & 
                        & ( myid &
                        &  ,lapm &
                        &  ,ne,nx,ny,nzz,nest,sy,ey,sx,ex &
                        &  ,ip3mem,mem3d &
                        &  ,ip4mem,mem4d &
                        &  ,igas,gas,GC_MOLWT &
                        &  ,flag )

use apm_varlist
implicit none
include 'apm_parm.inc'

integer :: myid
logical :: lapm
integer :: igas
integer :: ne,nest
integer :: nx(5),ny(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

character :: flag*2

integer :: mem3d,mem4d

!real,dimension(mem3d) :: acid_gas

real,dimension(mem4d) :: gas ! ppb,ug/m3

real, dimension(102)  :: GC_MOLWT

integer :: ip3mem(nzz,nest),ip4mem(nzz,igas,nest)

integer :: i03,i04,ixy

integer :: i,j,k,ig,idx,is

integer :: isulf

real              :: msulf,mwght
integer,parameter :: nsulf = 8
integer,parameter :: sulf_o_gas(1:nsulf) = (/ 1, 1, 1, 1, 1, 1, 1, 1/)
integer,parameter :: sulf_index(1:nsulf) = (/83,84,87,89,92,93,94,95/)


if(flag.eq.'11') then
 print*,'cal before gas-aerosol'
elseif(flag.eq.'ch') then
 print*,'cal after gas-aerosol'
endif

 loop_j : DO j = sy(ne),ey(ne)
 loop_i : DO i = sx(ne),ex(ne)

   ixy = (ex(ne)-sx(ne)+3)*(j -sy(ne)+1)+i-sx(ne)+1

   loop_k : DO k=1,nzz-1

     i03 = ip3mem(k,ne)


     msulf=0.0d0
     do isulf=1,nsulf
      idx   = sulf_index(isulf)
      mwght = 96.0 ! NH4
      ! kg/m3
      msulf = msulf + gas(ip4mem(k,idx,ne)+ixy)* &
                    & sulf_o_gas(isulf)*mwght/GC_MOLWT(idx)
     enddo

     if(flag.eq.'11') then
       ch_ratio(i03+ixy)=msulf
       acid_gas_1(i03+ixy)=gas(ip4mem(k,1,ne)+ixy)
       if(k.eq.1) then
         !print*,ch_ratio(i03+ixy),msulf,j,i
       endif
     elseif(flag.eq.'ch') then
       acid_gas_2(i03+ixy)=gas(ip4mem(k,1,ne)+ixy)
       if(ch_ratio(i03+ixy).ne.0) then
         ch_ratio(i03+ixy)=msulf/ch_ratio(i03+ixy)
       else
         ch_ratio(i03+ixy)=1
       endif
       if(k.eq.1) then
         !print*,ch_ratio(i03+ixy),msulf,j,i
       endif
       if(k.eq.1) then
         !print*,'i&j',i,j
         !print*,'sulferic acid',acid_gas_1(i03+ixy),acid_gas_2(i03+ixy)
       endif
     else
       stop 'flag erro in apm_isor_so4.f90 '
     endif

   enddo loop_k

 enddo loop_i
 enddo loop_j

!stop

end 






