
module smpsulf_var

integer,parameter :: nmoxdt=4,nm12=12,igassmp=5

integer :: idx_smpsulf_h2so4,idx_smpsulf_h2o2,idx_smpsulf_so2,idx_smpsulf_dms,idx_smpsulf_nh3

integer :: idx_oxdt_oh,idx_oxdt_ho2,idx_oxdt_no3,idx_oxdt_o3

integer,allocatable,dimension(:,:,:,:) :: ip4mem_oxdt

integer,allocatable,dimension(:,:,:) :: ip4mem_ox3d

real,allocatable,dimension(:) :: mmean_oxdt,oxdt3d,tmp4d

integer,allocatable,dimension(:) :: idx4dep_smpsulf,idx_smpsulf,idx4oxdt

contains

subroutine allo_smpchem_var(lgaschemsmp,nx,ny,nzz,nest,sx,ex,sy,ey,mem_per_block)

logical :: lgaschemsmp
integer :: mem4d_oxdt

integer :: mem4d_tmp

integer :: nest
integer :: nx(5),ny(5),nzz
integer :: sx(5), ex(5), sy(5), ey(5)

integer :: mem_per_block(5)

integer :: isp

allocate(idx4oxdt(nmoxdt))

do isp=1,nmoxdt
 idx4oxdt(isp)=isp
enddo

idx_oxdt_oh=1
idx_oxdt_ho2=2
idx_oxdt_no3=3
idx_oxdt_o3=4

allocate(idx_smpsulf(igassmp))

do isp=1,igassmp
idx_smpsulf(isp)=isp
enddo

idx_smpsulf_h2so4=1
idx_smpsulf_h2o2=2
idx_smpsulf_so2=3
idx_smpsulf_dms=4
idx_smpsulf_nh3=5

! map index to cbmz gaschem spc

allocate(idx4dep_smpsulf(igassmp))


idx4dep_smpsulf(idx_smpsulf_h2so4)=1
idx4dep_smpsulf(idx_smpsulf_h2o2)=16
idx4dep_smpsulf(idx_smpsulf_so2)=18
idx4dep_smpsulf(idx_smpsulf_dms)=57
idx4dep_smpsulf(idx_smpsulf_nh3)=4


if(lgaschemsmp) then
  allocate(ip4mem_oxdt(nzz,nmoxdt,nm12,nest))
  ii=1
  do ne=1,nest
  do imon=1,nm12
  do ig=1,nmoxdt
  do k=1,nzz
   ip4mem_oxdt(k,ig,imon,ne)=ii
   ii=ii+mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo
  mem4d_oxdt=0        ! number of memory for 4d oxidants
  do ne=1,nest
  do imon=1,nm12
  do ig=1,nmoxdt
  do k=1,nzz
    mem4d_oxdt=mem4d_oxdt+mem_per_block(ne)
  enddo
  enddo
  enddo
  enddo
  allocate(mmean_oxdt(mem4d_oxdt)); mmean_oxdt=0.0


  allocate(ip4mem_ox3d(nzz,nmoxdt,nest))
  ii=1
  do ne=1,nest
  do ig=1,nmoxdt
  do k=1,nzz
   ip4mem_ox3d(k,ig,ne)=ii
   ii=ii+mem_per_block(ne)
  enddo
  enddo
  enddo
  mem4d_oxdt=0        ! number of memory for 3d oxidants
  do ne=1,nest
  do ig=1,nmoxdt
  do k=1,nzz
    mem4d_oxdt=mem4d_oxdt+mem_per_block(ne)
  enddo
  enddo
  enddo
  allocate(oxdt3d(mem4d_oxdt)); oxdt3d=0.0



  mem4d_tmp=0
  do ne=1,nest
  do k=1,nzz
    mem4d_tmp=mem4d_tmp+mem_per_block(ne)
  enddo
  enddo
  allocate(tmp4d(mem4d_tmp)); tmp4d=0.0


endif



end subroutine allo_smpchem_var

end module smpsulf_var 

