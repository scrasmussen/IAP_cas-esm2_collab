      subroutine alloc_init
#include<def-undef.h>
      use  param_mod
      use  pmix_mod
      use dyn_mod
      use pconst_mod
      use output_mod
       use work_mod
       use forc_mod
       use tracer_mod
       use sw_mod
       use isopyc_mod
      use carbon_mod
      use cforce_mod
      use coutput_mod
      use grids_pt_mod
       implicit none
      allocate(ricdt(imt,jmt,kmm1))
!      real(r8),dimension(imt,jmt,kmm1):: ridu,ridt,s2u,s2t
      allocate (ridt(imt,jmt,kmm1))
      allocate (s2u(imt,jmt,kmm1))
      allocate (s2t(imt,jmt,kmm1))
#if (defined SOLAR)
      allocate( pen(kmm1))
#endif

#if(defined carbonC) || (defined carbonAbio)
      allocate(dpco2a(imt,jmt),nta_pt(imt,jmt))
#endif

#if (defined carbonBio)
      allocate(r_cap(km),delta_a(km),kappa_b(km))
      allocate(nta_pt(imt,jmt))
      allocate(a0_b(imt,jmt,km),a1_b(imt,jmt,km),a2_b(imt,jmt,km),b0_b(imt,jmt,km),&
               b1_b(imt,jmt,km),b2_b(imt,jmt,km),c_b(imt,jmt,km))
      allocate(Fe_scav_prof(km),ReFe(km))
      allocate(dust_in(imt,jmt,km),Fe_scav(imt,jmt,km))
      allocate(P_Fe(imt,jmt,kmp1))
      allocate(Fe_source(imt,jmt,km),TFe_scav(imt,jmt,km))
      allocate(r_fec(imt,jmt),r_fep(imt,jmt))
#ifdef SPMD
     allocate(flux_Fe_io(imt_global,jmt_global,km),Fe_scav_io(imt_global,jmt_global,km))
#endif      
#endif
#ifdef SPMD 
      allocate(DYR_IN_global(jmt_global),OUY_global(jmt_global),OTX_global(jmt_global),&
               OUX_global(jmt_global),SOTX_global(jmt_global),SOUX_global(jmt_global),&
               FF_global(jmt_global),CV1_global(jmt_global),CV2_global(jmt_global), &
               SNLAT_global(jmt_global),SINT_global(jmt_global), &
               SINU_global(jmt_global),DYT_global(jmt_global),DYR_global(jmt_global), &
               DXDYU_global(jmt_global),DXDYT_global(jmt_global), &
               R1A_global(jmt_global),R1B_global(jmt_global),R2A_global(jmt_global), &
               R2B_global(jmt_global),R1C_global(jmt_global), &
               R1D_global(jmt_global),R2C_global(jmt_global),R2D_global(jmt_global), &
               EBEA_global(jmt_global),EBEB_global(jmt_global), &
               EBLA_global(jmt_global),EBLB_global(jmt_global),EPEA_global(jmt_global),&
               EPEB_global(jmt_global),EPLA_global(jmt_global),EPLB_global(jmt_global))
      allocate(FF1_global(jmt_global))
#if (defined TSPAS)
      allocate(dtdy_global(jmt_global),dtdx_global(jmt_global),RAA_global(jmt_global),RBB_global(jmt_global))
#endif

#endif
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
      allocate(pco2a(imt,jmt),pco2o(imt,jmt),dpco2o(imt,jmt))
#endif
      
      allocate(pt(imt,jmt,km,nptra),ptb(imt,jmt,0:km,nptra),ptf(imt,jmt,km,nptra),sge(imt,jmt),ssfc(imt,jmt))
      
#ifdef SPMD
      allocate(pt_io(imt_global,jmt_global,km,nptra),ptb_io(imt_global,jmt_global,km+1,nptra))
#endif

      allocate(wstmon(imt_global,jmt_global,km),wst_io(imt_global,jmt_global,kmp1))
      allocate(tocaco3(imt,jmt,km),toa0(imt,jmt,km))
#ifdef SPMD
      allocate(tocaco3_io(imt_global,jmt_global,km),toa0_io(imt_global,jmt_global,km))
#endif

#ifdef carbonBio
      allocate(po4obs(imt,jmt,kmmix,12),taobs(km,12),po4force(imt,jmt,kmmix),taforce(km),&
              fe_flux(imt,jmt,12),dust_flux(imt,jmt,12),fe_f(imt,jmt),dust_f(imt,jmt))
#endif 
      allocate(w22np(imt,jmt),winds(imt,jmt,12),pressure(imt,jmt,12),pressureday(imt,jmt),iceday(imt,jmt))
#ifdef SPMD
      allocate(winds_global(imt_global,jmt_global,12))
      allocate(pressure_global(imt_global,jmt_global,12))
#ifdef carbonBio
      allocate(po4obs_global(imt_global,jmt_global,kmmix,12))
#endif
#endif
      
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
      allocate(tpco2o(imt,jmt),tdpco2o(imt,jmt),totup(imt,jmt))
#endif

#ifdef carbonBio
      allocate(po4mon(imt,jmt,km),ldocmon(imt,jmt,km),tamon(imt,jmt,km))
      allocate(prodmon(imt,jmt,km),fpopmon(imt,jmt,km),pldocmon(imt,jmt,km),remimon(imt,jmt,km))
      allocate(jpopmon(imt,jmt,km),caco3mon(imt,jmt,km))
      allocate(o2mon(imt,jmt,km),femon(imt,jmt,km))
      allocate(o2up(imt,jmt),feup(imt,jmt))
#ifdef COUP
      allocate(uptake(imt,jmt))
#endif
#endif

#ifdef carbonC14
      allocate(ssfcmon(imt,jmt))
#ifdef SPMD
      allocate(ssfcmon_io(imt_global,jmt_global))
#endif
#endif

#ifdef cfc
      allocate(ssfcmon(imt,jmt))
#ifdef SPMD
      allocate(ssfcmon_io(imt_global,jmt_global))
#endif
#endif
      allocate(ccmon(imt,jmt,km))

#ifdef SPMD
      allocate(ccmon_io(imt_global,jmt_global,km))
#if (defined carbonC) || (defined carbonBio) || (defined carbonAbio)
      allocate(tpco2o_io(imt_global,jmt_global),tdpco2o_io(imt_global,jmt_global),totup_io(imt_global,jmt_global))
#endif
#ifdef carbonBio
      allocate(po4mon_io(imt_global,jmt_global,km),ldocmon_io(imt_global,jmt_global,km),&
               tamon_io(imt_global,jmt_global,km))
      allocate(prodmon_io(imt_global,jmt_global,km),fpopmon_io(imt_global,jmt_global,km),&
               pldocmon_io(imt_global,jmt_global,km),remimon_io(imt_global,jmt_global,km))
      allocate(jpopmon_io(imt_global,jmt_global,km),caco3mon_io(imt_global,jmt_global,km))
      allocate(o2mon_io(imt_global,jmt_global,km),femon_io(imt_global,jmt_global,km))
!for flux,lyc 2013.06
      allocate(O2up_io(imt_global,jmt_global),Feup_io(imt_global,jmt_global))
#endif
#endif
      allocate(cct(imt,jmt_global))

      allocate(zb(km),ztop(km),zm(km),sdxdy(km),sdxdy0(km))
      allocate(kmt(imt,jmt))

#if (defined SOLARCHLORO)
      allocate(pen_chl(imt,jmt,km))   !  be different from
#endif

      allocate(ub(imt,jmt),vb(imt,jmt),ubp(imt,jmt),vbp(imt,jmt),h0p(imt,jmt))
      allocate(up(imt,jmt,km),vp(imt,jmt,km))
      allocate(ws(imt,jmt,kmp1))
      allocate(h0l(imt,jmt),h0f(imt,jmt),h0bl(imt,jmt),h0bf(imt,jmt))
      allocate(utl(imt,jmt,km),utf(imt,jmt,km),vtl(imt,jmt,km),vtf(imt,jmt,km))
!      allocate(vetiso(imt,jmt,km),vntiso(imt,jmt,km),vbtiso(imt,jmt,km))
#ifdef COUP
      allocate(t_cpl_io(imt_global,jmt_global),s_cpl_io(imt_global,jmt_global), &
               u_cpl_io(imt_global,jmt_global),v_cpl_io(imt_global,jmt_global), &
               dhdx_io(imt_global,jmt_global),dhdy_io(imt_global,jmt_global),q_io(imt_global,jmt_global))
#endif
      allocate( j_global(jmt))
      allocate(i_global(imt))
!#endif
      allocate( vit(imt,jmt,km),viv(imt,jmt,km))
      allocate(vit_global_surface(imt_global,jmt_global),viv_global_surface(imt_global,jmt_global))
!      allocate(vit_ori(imt,jmt))
!jjr      real(r8),dimension(imt_global,jmt,km):: vit_1d,viv_1d
      allocate(vit_1d(imt_global,jmt,km),viv_1d(imt_global,jmt,km))
      allocate( ahv_back(jmt_global))
      allocate(basin(imt_global,jmt_global))
      allocate(itnu(imt,jmt),na(imt,jmt))
#if (defined NETCDF) || (defined ALL)
      allocate( lon(imt_global))
      allocate( lat(jmt_global))
      allocate( lev(km))
      allocate( lev1(km+1))
!      allocate( lon_local(imt))
!      allocate( lat_local(jmt))
#endif
!lhl090729
      allocate(s_lon(s_imt))
      allocate( s_lat(s_jmt))
!lhl090729
      allocate(zkt(km),dzp(km),odzp(km),odzt(km))
      allocate(zkp(kmp1))


      allocate(DYR_IN(jmt))
      allocate(OUY(jmt),OTX(jmt),OUX(jmt),SOTX(jmt),SOUX(jmt))
      allocate(FF(jmt),CV1(jmt),CV2(jmt),SNLAT(jmt),SINT(jmt))
      allocate(SINU(jmt),DYT(jmt),DYR(jmt),DXDYU(jmt),DXDYT(jmt))
      allocate(R1A(jmt),R1B(jmt),R2A(jmt),R2B(jmt),R1C(jmt))
      allocate(R1D(jmt),R2C(jmt),R2D(jmt),EBEA(jmt),EBEB(jmt))
      allocate(EBLA(jmt),EBLB(jmt),EPEA(jmt),EPEB(jmt),EPLA(jmt),EPLB(jmt))
!lhl060506
      allocate(FF1(jmt),RRD1(jmt),RRD2(jmt))
!lhl060506

!XC
#if (defined TSPAS)
      allocate(dtdy(jmt),dtdx(jmt),RAA(jmt),RBB(jmt))
#endif
!XC



#if ( defined SMAG)
      allocate(CXT(jmt),CXU(jmt),CYT(jmt),CYU(jmt) &
                   ,R1E(jmt),R1F(jmt),R2E(jmt),R2F(jmt) &
                   ,R3E(jmt),R3F(jmt),R4E(jmt),R4F(jmt))
#endif

!Yu

      allocate(ohbt(imt,jmt),ohbu(imt,jmt),dzph(imt,jmt),hbx(imt,jmt),hby(imt,jmt))
      allocate(COSU(jmt),COST(jmt))
      allocate(CF1(imt),CF2(imt),SF1(imt),SF2(imt))
!
!
!     -----------------------------------------------------------
!     Reference T S & coefficients for calculation of d(density)
!     -----------------------------------------------------------
!YU
      allocate( TO(KM),SO(KM),PO(KM))
!
      allocate(C(km,9))
      allocate(akmu(imt,jmt,km),akmt(imt,jmt,km))
       allocate(akt(imt,jmt,km,NTRA))

      allocate(AM(jmt),AH(jmt))
      allocate(am3(imt,jmt,km),ah3(imt,jmt,km),amx(imt,jmt,km),amy(imt,jmt,km))
#if ( defined TIDEMIX )
!     -----------------------------------------------------------
!     Tidal mixig
!     -----------------------------------------------------------
       allocate(fz_tide(imt,jmt,km),ak_tide(imt,jmt,km))
#endif

      allocate(z0mon(imt,jmt),himon(imt,jmt),hdmon(imt,jmt))
!
!linpf091126
      allocate(lthfmon(imt,jmt),sshfmon(imt,jmt),lwvmon(imt,jmt),swvmon(imt,jmt))
      allocate(sumon(imt,jmt),svmon(imt,jmt))
!linpf091126
!
      allocate(wsmon(imt,jmt,km),tsmon(imt,jmt,km),&
         ssmon(imt,jmt,km),usmon(imt,jmt,km),vsmon(imt,jmt,km))
      allocate(icmon(imt,jmt,2))
      allocate(netmon(imt,jmt,NTRA))
!
!lhl1204
      allocate(mldmon(imt,jmt))
      allocate(akmmon(imt,jmt,km),aktmon(imt,jmt,km),aksmon(imt,jmt,km))
!lhl1204
!lhl20131025
#if ( defined TIDEMIX )
      allocate(aktidemon(imt,jmt,km))
#endif
!lhl20131025
!      allocate(vetisomon(imt,jmt,km),vntisomon(imt,jmt,km),vbtisomon(imt,jmt,km))
!
!wangty modify
!      allocate(mth(jmt_global,3,ntra),mth_adv(jmt_global,3,ntra))
!      allocate(mth_dif(jmt_global,3,ntra),mth_adv_iso(jmt_global,3,ntra))
      allocate(mth(jmt_global,2,ntra),mth_adv(jmt_global,2,ntra))
      allocate(mth_dif(jmt_global,2,ntra),mth_adv_iso(jmt_global,2,ntra))

!
      allocate(trendmon(imt,jmt,km,ntra),axmon(imt,jmt,km,ntra))
      allocate(aymon(imt,jmt,km,ntra),azmon(imt,jmt,km,ntra))
      allocate(dxmon(imt,jmt,km,ntra),dymon(imt,jmt,km,ntra))
      allocate(dzmon(imt,jmt,km,ntra))
      allocate(penmon(imt,jmt,km))
!
      allocate(ddymon(imt,jmt,km,ntra))

#ifdef ISO
      allocate(axmon_iso(imt,jmt,km,ntra),aymon_iso(imt,jmt,km,ntra))
      allocate(azmon_iso(imt,jmt,km,ntra),dxmon_iso(imt,jmt,km,ntra))
      allocate(dymon_iso(imt,jmt,km,ntra),dzmon_iso(imt,jmt,km,ntra))
      allocate(aaymon_iso(imt,jmt,km,ntra),ddymon_iso(imt,jmt,km,ntra))
#endif
!
#if (defined SMAG_OUT)
      allocate(a3mon(imt,jmt,km))
#endif
      allocate(atb(imt,jmt,0:km,NTRA))
      allocate(net(imt,jmt,NTRA))
!mohr
      allocate(pdensity(imt,jmt,km))
      allocate(amld(imt,jmt))
!lhl1204
!
      allocate(trend(imt,jmt,km,NTRA),ax(imt,jmt,km,ntra))
      allocate(ay(imt,jmt,km,ntra),az(imt,jmt,km,ntra))
      allocate(dx(imt,jmt,km,ntra),dy(imt,jmt,km,ntra),dz(imt,jmt,km,ntra))
      allocate(penetrate(imt,jmt,km))
!
      allocate(ddy(imt,jmt,km,NTRA))
!
#ifdef ISO
      allocate(aay_iso(imt,jmt,km,ntra),ddy_iso(imt,jmt,km,ntra))
      allocate(ax_iso(imt,jmt,km,ntra),ay_iso(imt,jmt,km,ntra))
      allocate(az_iso(imt,jmt,km,ntra),dx_iso(imt,jmt,km,ntra))
      allocate(dy_iso(imt,jmt,km,ntra),dz_iso(imt,jmt,km,ntra))
#endif
!
!
!#endif
!
!     ------------------------------------------------------------------
!     Sea Ice
!     ------------------------------------------------------------------
!      real(r8),dimension(imt,jmt):: ITICE,ALEAD,TLEAD, HI
      allocate(TLEAD(imt,jmt),licomqice(imt,jmt))
#ifdef SPMD
      allocate(qice_global(imt_global,jmt_global))
#endif
      allocate(su(imt,jmt),sv(imt,jmt),psa(imt,jmt),tsa(imt,jmt))
      allocate(sss(imt,jmt),swv(imt,jmt),uva(imt,jmt),qar(imt,jmt),cld(imt,jmt))
      allocate(ddd(imt,jmt),qqq(imt,jmt),sst(imt,jmt),nswv(imt,jmt))
!lhl
       allocate(dqdt(imt,jmt),chloro(imt,jmt),lwv(imt,jmt),seaice(imt,jmt))
       allocate(rain(imt,jmt),snow(imt,jmt),fresh(imt,jmt),runoff(imt,jmt))
       allocate(lthf(imt,jmt),sshf(imt,jmt)) !only for output
!linpf091126

!lhl1204
      allocate(USTAR(imt,jmt),BUOYTUR(imt,jmt),BUOYSOL(imt,jmt))
!lhl
!#endif
!
!
#if (defined FRC_DAILY)
      allocate(su_in_io(imt_global,jmt_global),sv_in_io(imt_global,jmt_global))
      allocate(su_io(imt_global,jmt_global),sv_io(imt_global,jmt_global))
#endif
#if (defined BOUNDARY)
      allocate(restore(imt,jmt,km,ntra))
      allocate(restore_io(imt,jmt_global,km,ntra))
      ! lihuimin 2012.6.18
      ! modi end
#endif
!
#ifdef SPMD
       allocate( tsf_global(imt_global,jmt_global),ssf_global(imt_global,jmt_global),&
                 su_global(imt_global,jmt_global),sv_global(imt_global,jmt_global),&
                 swv_global(imt_global,jmt_global), mius_global(imt_global,jmt_global),&
                 fresh1_global(imt_global,jmt_global),u3_global(imt_global,jmt_global),&
                 sss1_global(imt_global,jmt_global),fresh11_global(imt_global,jmt_global))
#endif
       allocate(tsf(imt,jmt),ssf(imt,jmt),mius(imt,jmt),fresh1(imt,jmt),&
                u3(imt,jmt),sss1(imt,jmt),fresh11(imt,jmt))
#if(defined SOLARCHLORO)
       allocate(ztr(km))
       allocate(Tr(0:km-1,0:nsub))
#endif
      allocate(PXB(imt,jmt),PYB(imt,jmt),PAX(imt,jmt),PAY(imt,jmt),WHX(imt,jmt),WHY(imt,jmt),WGP(imt,jmt))
      allocate(wka(imt,jmt,km))
      allocate(work(imt,jmt),work1(imt,jmt),work2(imt,jmt))
      allocate( wki(imt))
      allocate(wkj(jmt))
      allocate(wkk(kmp1))

#ifdef SPMD
      allocate( wkj_global(jmt_global))
#endif

#if (defined ISO)
      allocate(fzisop(km))
      allocate(kisrpl(km))
      allocate(zt(km),dzw(0:km),dzwr(0:km),dzr(km),dytr(jmt),cstrdytr(jmt),dyur(jmt),tmask(imt,km,jmt))
      allocate(F3(jmt))
#endif
      allocate(psi(jmt_global,km+1,3))
      allocate(bsf(imt_global,jmt_global))
      
      end subroutine alloc_init
