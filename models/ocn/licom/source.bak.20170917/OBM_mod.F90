! CVS: $Id: carbon_mod.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  module carbon_mod
!========================
! carbon_mod
!---------------------------------------------------------------------
!
! author: Zhao Liang@lapc 2004/02/29
!
!----------------------------------------------------------------------
! IMT,JMT,KM      - is defined in module param
! NPTRA           - number of tracer
! sge             - air-sea exchange coefficient of carbon
! ssfc            - carbon flux
! kmmix           - the number of mixing layer (9 temperature layer)
! ze              - the depth of mixing layer (100 meter in this study)
! pco2a           - atmospheric partial pressure of carbon dioxide
! dpco2a          - pertubation of atmospheric partial pressure
! pco2o           - oceanic partial pressure of carbon dioxide
! dpco2o          - atmospheric partial pressure minus oceanic one
!----------------------------------------------------------------------
! r_cp,r_np,r_cap - read from file in ctrl.F90
! r_o2p,r_fep           - read from file in ctrl.F90
! delta_a         -
! kappa_b         - calculated by the variablity of different variables
! r_rain          - the parameter named rain ratio (the ratio of calcite to POC)
! d_cal           - the eholdting length scale of calcite flux
!----------------------------------------------------------------------
! a0_b            - the production of POP
! a1_b            - new production
! a2_b            - flux of POP below the euphotic zone
! b0_b            - the source of LDOC
! b1_b            - the production of LDCO in the euphotic zone
! b2_b            - the LDOC remineralization
! c_b             - the source minus sink of CaCO3
!---------------------------------------------------------------------
! ldocga          - the global average concentration of ldoc
 !---------------------------------------------------------------------
! ar0             - a proportional factor( bio-production efficiency)
! hf              - a half-saturation constant
! zdc             - the decading coeffient by the depth
!--------------------------------------------------------------------------

#include <def-undef.h>
!
      USE param_mod
!
!----------------------------------------------------------------------
      IMPLICIT NONE
#if(defined carbonC) || (defined carbonAbio)
      integer,parameter::nptra=1
      real,dimension(:,:),allocatable::dpco2a
      real,dimension(:,:),allocatable::nta_pt
#endif      
!
#if (defined carbonBio)
!lyc,2013.05
      integer,parameter::nptra=6 !+oxygen +iron
      integer,parameter::kmmix=10
      real::ldocga=4.2
      real,parameter::POP_k=-0.858
!      real::ze=110.0    ! unit is meter
! add biological model variable 
      real::r_cp,r_np,tau_b,r_o2p,r_fec0
      real,parameter::r_rain=0.07
      real,parameter::d_cal=3500.0
      real*8::delldoc
      real::sigma_b
      real,dimension(:),allocatable::r_cap,delta_a
      real,dimension(:),allocatable::kappa_b
      real,dimension(:,:),allocatable::nta_pt
! biological source, suffix 'b' mean biological variables       
      real,dimension(:,:,:),allocatable::a0_b,a1_b,a2_b,b0_b,b1_b,b2_b,c_b
      real,parameter::hf=0.02
      real,parameter::ar0=0.5/30/86400
      real,parameter::zdc=23
!lyc,2013.05
!---some constant for iron cycle------------------------------------
     logical:: DUST_DATA,FE_FLUX_DATA
     real,parameter::hf_Fe=0.3E-4 !umol/l
     real,parameter::fepos=0.6  !only 60% of Fe_salv can be treated as pofe
     real,parameter::Fe_s=250 !m for the Fe scavenge
     real,parameter::Fe_a =-0.9  !for the single Martin power-law curve
     real,parameter::Fe_bioava=0.02 !the  bioavailable fe from dust
     real,dimension(:),allocatable::Fe_scav_prof,ReFe !ReFe the remineralization of POFe(including 0.6*Fe_scav)
     real,dimension(:,:,:),allocatable::dust_in
     real,dimension(:,:,:),allocatable::Fe_scav
     real,dimension(:,:,:),allocatable::P_Fe
     real,dimension(:,:,:),allocatable::Fe_source
     real,dimension(:,:,:),allocatable::TFe_scav !the column integrate
     real,dimension(:,:),allocatable::r_fec,r_fep
!     real,dimension(imt,jmt,km)::Fe_free
#ifdef SPMD
     real,dimension(:,:,:),allocatable::flux_Fe_io
     real,dimension(:,:,:),allocatable::Fe_scav_io
!     real,dimension(imt,jmt_global,km)::Fe_free_io
#endif
#endif
!lyc
#ifdef cfc
      integer,parameter::nptra=1
      integer,parameter::LM=2
      real::t0,s0
#endif
!
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
      integer::nsg1
      real,parameter::unitc=1025.0*12.011*1.0E-21
      real,dimension(:,:),allocatable::pco2a,pco2o,dpco2o
#endif
!
#ifdef carbonC14
      integer,parameter::nptra=1
      integer::ny
#endif
!
! global variables for carbonC, carbonC14 and carbonBio
      real,dimension(:,:,:,:),allocatable::pt
      real,dimension(:,:,:,:),allocatable::ptb
      real,dimension(:,:,:,:),allocatable::ptf
      real,dimension(:,:),allocatable::sge,ssfc
      integer::yearR,monthR
      integer::NSTARTC
      integer::isp
!
#ifdef SPMD
      real,dimension(:,:,:,:),allocatable::pt_io
      real,dimension(:,:,:,:),allocatable::ptb_io
#endif      
 real,dimension(:,:,:),allocatable::wstmon
 real,dimension(:,:,:),allocatable::wst_io
!
!add the output of a0 and caco3
real,dimension(:,:,:),allocatable::tocaco3,toa0
#ifdef SPMD
      real,dimension(:,:,:),allocatable::tocaco3_io,toa0_io
#endif   

integer:: NPP,IDTP
!
  end module carbon_mod


! CVS: $Id: cforce_mod.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  module cforce_mod
!========================
! cforce_mod
!---------------------------------------------------------------------
!
! author: Zhao Liang@lapc 2004/02/29
!
!---------------------------------------------------------------------
! IMT,JMT,KM     - is defined in module param
! ny,nsg1        - the year number of force data
! boml,bomm,bomh - atmospheric bomb carbon(C14)(low,mediate,high)
! csg            - observed atmospheric partial pressure of carbon dioxide (uatm)
! catm           - atmospheric carbon of every day
! w22np          - wind speed of every day intepolated from winds
! po4force       - po4 forcing data of every day intepolated from po4obs
! taforce        - ta forcing data of every day intepolated from taobs
! winds          - monthly mean observed wind speed
! po4obs         - monthly mean observed concentration of PO4 (WOA2001)
! taobs          - global horizontal mean observed TA (GEOSECS)
!----------------------------------------------------------------------
#include <def-undef.h> 
!
      USE param_mod
#ifdef carbonBio      
      USE carbon_mod
#endif      
!---------------------------------------------------------------------
      IMPLICIT NONE
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
      real,parameter::pco2dry0=280.77
      real::pco2dry
      real,allocatable,dimension(:)::csgn,csg,csgx
#endif
!
#ifdef carbonBio
! biological model forcing data
      real,dimension(:,:,:,:),allocatable::po4obs
      real,dimension(:,:),allocatable::taobs
      real,dimension(:,:,:),allocatable::po4force
      real,dimension(:),allocatable::taforce
!lyc 2013,07
      real,dimension(:,:,:),allocatable::fe_flux,dust_flux !umol/cm2/s->10m*umol/l/s->10m*umol/kg/s/1.025 
      real,dimension(:,:),allocatable::fe_f,dust_f  !iron flux in a day
#endif
!
#ifdef carbonC14
      integer,allocatable,dimension(:)::kyear
      real,allocatable,dimension(:)::boml,bomm,bomh
      real,dimension(3)::catm
!cm080418
      real,parameter::rdca=1.2097e-4/3.1536e7  !ratio of carbonC14's decay
!cm080418
#endif
!lyc
#ifdef cfc
     real,dimension(2)::cfcatm
     real::sol,xsol
     integer,parameter::ny=54
     real,dimension(ny)::bomn,boms
#endif
!
! global variables for carbonC, carbonC14 and carbonBio
      integer::yearData,yearStart
      real,dimension(:,:),allocatable::w22np
      real,dimension(:,:,:),allocatable::winds
!cm090330---------------------------------
      real,dimension(:,:,:),allocatable::pressure
      real,dimension(:,:),allocatable::pressureday
      real,dimension(:,:),allocatable::iceday
!cm090330---------------------------------
!      
#ifdef SPMD
      real,dimension(:,:,:),allocatable::winds_global
!cm090330---------------------------------
      real,dimension(:,:,:),allocatable::pressure_global
!cm090330---------------------------------
#ifdef carbonBio
      real,dimension(:,:,:,:),allocatable::po4obs_global
#endif
#endif      
!
  end module cforce_mod
! CVS: $Id: coutput_mod.F90,v 2.1 2004/06/10 07:45:17 cvsroot Exp $
  module coutput_mod
!========================
! coutput_mod
!---------------------------------------------------------------------
!
! author: Zhao Liang@lapc 2004/02/29
!
!---------------------------------------------------------------------
! IMT,JMT,KM     - is defined in module param
! tpco2o         - the sum of oceanic partial pressure of carbon dioxide
! tdpco2o        - the sum of atmospheric partial pressure minus oceanic one
! totup          - the sum of flux of carbon dioxide
! jpopmon        - a0_b
! prodmon        - a1_b
! fpopmon        - a2_b
! pldocmon       - b1_b
! remimon        - b2_b
! caco3mon       - c_b
! ccmon          - monthly mean carbon
! ssfcmon        - the monthly mean fluxes of cfc
!----------------------------------------------------------------------
#include <def-undef.h>
!
      USE param_mod
!
!----------------------------------------------------------------------
      IMPLICIT NONE
#ifdef carbonC

#endif
!
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
      real,dimension(:,:),allocatable::tpco2o,tdpco2o,totup
#endif
!
#ifdef carbonBio
      real,dimension(:,:,:),allocatable::po4mon,ldocmon,tamon 
      real,dimension(:,:,:),allocatable::prodmon,fpopmon,pldocmon,remimon
      real,dimension(:,:,:),allocatable::jpopmon,caco3mon
      real,dimension(:,:,:),allocatable::o2mon,femon  !lyc,2013.05
!for flux
      real,dimension(:,:),allocatable::o2up,feup
#ifdef COUP
      real*8,dimension(:,:),allocatable::uptake
#endif
#endif
!
#ifdef carbonC14
   real,dimension(:,:),allocatable::ssfcmon
#ifdef SPMD
   real,dimension(:,:),allocatable::ssfcmon_io
#endif
#endif
!lyc
#ifdef cfc
   real,dimension(:,:),allocatable::ssfcmon
#ifdef SPMD
   real,dimension(:,:),allocatable::ssfcmon_io
#endif
#endif
!lyc
!
! global variables for carbonC, carbonC14 and carbonBio
      integer::IO_out,yearStore
      real,dimension(:,:,:),allocatable::ccmon
!      
#ifdef SPMD
      real,dimension(:,:,:),allocatable::ccmon_io
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
      real,dimension(:,:),allocatable::tpco2o_io,tdpco2o_io,totup_io
#endif
#ifdef carbonBio
     real,dimension(:,:,:),allocatable::po4mon_io,ldocmon_io,tamon_io
     real,dimension(:,:,:),allocatable::prodmon_io,fpopmon_io,pldocmon_io,remimon_io
     real,dimension(:,:,:),allocatable::jpopmon_io,caco3mon_io
     real,dimension(:,:,:),allocatable::o2mon_io,femon_io !lyc 2013.06
!for flux,lyc 2013.06
     real,dimension(:,:),allocatable::O2up_io,Feup_io
#endif
#endif      
!
! lyc  variables for netcdf
!   integer::ttVarID,ssVarID,ccVarID
   integer::ncFileID
   integer::ccVarID,tcVarID,po4VarID,ldocVarID,taVarID,o2VarID,feVarID
   integer::prodVarID,fpopVarID,pldocVarID,remiVarID,jpopVarID,caco3VarID
   real,dimension(:,:),allocatable::cct

  end module coutput_mod
