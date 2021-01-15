!-----------------------------------------------------------------------
! Update: TianyiWang,JinrongJiang,2016.01.16, variables about Processor partitioning
!                                             to nameliset
!-----------------------------------------------------------------------
module dyn_mod
#include <def-undef.h>
use precision_mod
use param_mod
!     ------------------------------------------------------------------
!     U V T S H0 W RHO
!     ------------------------------------------------------------------

!TianyiWang,JinrongJiang,2016.01.16
      real(r8),dimension(:,:),allocatable::ub,vb,ubp,vbp,h0p
      real(r8),dimension(:,:,:),allocatable::up,vp
      real(r8),dimension(:,:,:),allocatable::ws
      real(r8),dimension(:,:),allocatable::h0l,h0f,h0bl,h0bf
      real(r8),dimension(:,:,:),allocatable::utl,utf,vtl,vtf

     real(r8),dimension(:,:,:),allocatable::viso,uiso,wiso
#ifdef omipout
     real(r8),dimension(:,:,:),allocatable::psur,psur1,ritsave,rict1
#endif

!TianyiWang,JinrongJiang,2016.01.16

!#ifdef SPMD
!      real(r8),dimension(imt_global,jmt_global,km)::u_io,v_io
!      real(r8),dimension(imt_global,jmt_global)::h0_io
!mohr
      real(r8),allocatable,dimension(:,:):: buffer
      real(r8),allocatable,dimension(:,:)::h0
      real(r8),allocatable,dimension(:,:,:)::u,v
!
#ifdef COUP
!TianyiWang,JinrongJiang,2016.01.16
      real(r4),dimension(:,:),allocatable::t_cpl_io,s_cpl_io,u_cpl_io,v_cpl_io,dhdx_io,dhdy_io,q_io
!TianyiWang,JinrongJiang,2016.01.16
#endif
!#endif
!
!
!     ------------------------------------------------------------------
!     Pressure gradient
!     ------------------------------------------------------------------
      real(r8),dimension(:,:,:),allocatable::gg,dlu,dlv
      real(r8),dimension(:,:),allocatable::dlub,dlvb
end module dyn_mod

