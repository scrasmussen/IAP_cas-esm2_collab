  deallocate(ip2mem)
  deallocate(TERRAIN, RAINCON, RAINNON, T2, SWDOWN, &
             LATITCRS, LONGICRS, LAND_USE, UST, U10, V10,RMOL,&
             PBL_HGT,NPBL,clflo,clfmi,clfhi,HGT1,CLDOPD)
  deallocate(ip2memGas)
  deallocate(EmtaGas,EmtpGas,EmttGas,EmtbGas,DryVelGas)
  DEALLOCATE(WETDEP,WETDEP2)

!  deallocate(DryDGas,WetDGas)

  deallocate(ip3mem)
  deallocate(u,v,t,h,w,dx,dy,dz,rnw,clw,rh1,Plev,TAUCLDI,TAUCLDC)
!--------------lijie modify--------------------
 deallocate(globalno2,globalo3,globalco)
!----------------finish------------------------
  ! Zifa 2004/09/02
  if(imasskeep==1)then
  deallocate(kpmass_m1, RatioMass ,kpmass_m2)
  endif

  deallocate(ip4mem)
  deallocate(gas)
!  deallocate(gastmp) !acm2
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
  deallocate(gasOLD)
  deallocate(GasTermBal,ipGasTermBal)
  deallocate(IGGPOS)
!++++++++++++++++++++++++for process by lijie+++++++++++++++++++++++
! *** DUST and SEA SALT
  DEALLOCATE(GRAVEL,DUSTK,DUSTHGTF,TOTALDUST)
  DEALLOCATE(FICE,FSNOW,SOILT,SOILRH,FVEG,EMITFACT)
  DEALLOCATE(DUSTEMISS,DUSTDRY,DUSTWET,DUSTGRAV)
  DEALLOCATE(DUSTDRYSO4,DUSTDRYNO3,DUSTDRYFEII,DUSTDRYFEIII)
  DEALLOCATE(DUSTWETSO4,DUSTWETNO3,DUSTWETFEII,DUSTWETFEIII)
  DEALLOCATE(DUSTGRAVSO4,DUSTGRAVNO3,DUSTGRAVFEII,DUSTGRAVFEIII)
  DEALLOCATE(SEAEMISS)
  DEALLOCATE(RK_HETSO2_DUST,RK_HETHNO3_DUST)
! ***  AQUEOUS CHEMISTRY ******
  deallocate(CPH)
  DEALLOCATE (OPE) ! OPE
!=========================cloud and convection=========================
!  deallocate(ppp,ttn,ffn,conc,ip2mem2dconv,&
!                             ip3mem3dconv)
                              
!========================finished======================================
  deallocate(ktop,tropp) 
  deallocate(ASO4,ANO3,ACL,ANA,ANH4)
  deallocate(ip5mem)
  deallocate (ip5memc,ip5memcs) ! DUST and SEA SALT
  deallocate(aer,aer_src)
  deallocate(jo1d,jno2)
  deallocate(DUSTEXT,EXT,VISIB,SSA,AOD,PBLAOD,DUSTAOD)
  deallocate(EXTASO4,EXTANO3,EXTANH4,EXTBC,EXTOC)
  deallocate(DUSO2,DUNO2,DUO3)
  deallocate(UVB,UVBS,UVA,VIS)
  deallocate(syc,sxc,exc,eyc)
  deallocate(atestR,atestS)
  deallocate(atestR0,atestS0) ! juanxiong he
  !!!!!!!!!!!!!!!!!!!!!!
  ! for Source Mark
  if(ifsmt>0)then
  deallocate(igMark,iaMarkAer,iaMarkSiz)
  deallocate(MapSource,tmpMarkCon)
  deallocate(SourceMark)
  endif
  !!!!!!!!!!!!!!!!!!!!!!


  if(allocated(STAMARK)) deallocate(STAMARK)











