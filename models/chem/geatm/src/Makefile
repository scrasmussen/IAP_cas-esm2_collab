#
.SUFFIXES:.F90 .f90 .f .o

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f
SFC             =       ifort 
DM_FC           =       mpiifort -f90=$(SFC) -mcmodel=medium 
FC              =       $(DM_FC)
FORMAT_FREE     =       -FR
BYTESWAPIO      =      -convert big_endian 
# FCDEBUG       =  -O3 -no-prec-div -traceback -fp-model source -assume byterecl -ftz -safe-cray-ptr -fno-alias -fno-fnalias -m64  $(BYTESWAPIO)
# FCDEBUG       =    -w -ftz -fp-model precise -IPF_fp_relaxed -ip -fno-alias -fno-fnalias -m64 \
                     -g -traceback -check all -check noarg_temp_created -fstack-protector -assume protect_parens \
                     $(BYTESWAPIO)
#FCBASEOPTS      = -O3 -traceback -w -ftz -fp-model fast=2 -IPF_fp_relaxed -ip -fno-alias -fno-fnalias -safe-cray-ptr $(BYTESWAPIO)
FCBASEOPTS      = -O3 -traceback -w -ftz -fp-model precise -IPF_fp_relaxed -ip -fno-alias -fno-fnalias -g -safe-cray-ptr $(BYTESWAPIO)
AR              =      ar
ARFLAGS         =      ru

naq_moddir = ./module_main
naq_dyndir = ./naqpms_dyn
smpbox_dir = ./smpgas_chem
GAS_INC = ./include
GAS_SRC = ./gaschem
AEROSOL_SRC = ./aerosol
LIQ_SRC = ./aqeous
MODEL   = ./model
MARK   = ./mark
APM_SRC = ./apm_phys_box/src
APM_INC = ./apm_phys_box/src
APM_DYN = ./apm_dyn
APM_MEM = ./apm_var
AQCHEM_DIR = ./aqchem
EMSFRC_DIR = ./emit.sen

INCS = /public/software/mpi/intelmpi/5.0.2.044/include64 \
	-I$(GAS_INC) -I$(naq_moddir)\
	-I$(APM_INC) -I$(APM_DYN)

NETCDF = /public/software/mathlib/netcdf/4.3.2/intel

LIBTARGET = libgeatmlib.a

VPATH = model:model_main:include: \
	aerosol:mark:apm_var:apm_dyn:\
	apm_phys_box:aqchem:aqeous: \
	naqpms_dyn:gaschem:smpgas_chem: \
	emit.sen:main  

OBJS    = 
	   
MODULES =  \
	$(naq_moddir)/naqpms_var.o \
	$(naq_moddir)/smp_sulfchem_var.o \
	$(naq_moddir)/met_fields.o \
	$(naq_moddir)/work_vars.o \
	$(naq_moddir)/naqpms_grid.o \
	$(naq_moddir)/cputime.o \
	$(naq_moddir)/flags4wdep.o \
	$(naq_dyndir)/adv1d_comv.o \
	$(smpbox_dir)/sulfur_smp_gaschem.o \
	$(MODEL)/dep_parm.o \
	$(MODEL)/wrf_micromet.o \
	$(MODEL)/kvcalc_ysu.o \
	$(MODEL)/kvcalc_myj.o \
	$(MODEL)/kvcalc_cmaq.o \
	$(MODEL)/kvcalc_ob70.o \
	$(APM_DYN)/sub_max_ij.o \
	$(APM_DYN)/print_time.o \
	$(APM_DYN)/get_funit.o \
	$(APM_SRC)/apm_opti_mod.o \
	$(APM_SRC)/apm_nucl_mod.o \
	$(APM_SRC)/apm_init_mod.o \
	$(APM_SRC)/apm_coag_mod.o \
	$(APM_SRC)/apm_grow_mod.o \
	$(APM_SRC)/apm_phys_mod.o \
	$(APM_MEM)/apm_var.o \
	$(AQCHEM_DIR)/aqchem_var.o \
	$(AQCHEM_DIR)/sub_get_cldflag.o \
	$(AQCHEM_DIR)/aqchem_driver.o \
	$(AQCHEM_DIR)/sub_cal_inorganic.o \
	$(APM_DYN)/sub_cal_yspgf.o \
	$(APM_DYN)/v3_apm_inout.o \
	$(APM_DYN)/test_apm_io4nucl.o \
	$(APM_DYN)/apm_inout-aq.o \
	$(APM_DYN)/initialize_apm.o \
	$(APM_DYN)/apm_rst.o \
	$(APM_DYN)/apm_5bdy.o \
	$(APM_DYN)/apm_snstvy_emit.o \
	$(APM_DYN)/apm_zero_emit.o \
	$(APM_DYN)/apm_bcoc_sulf_emit.o \
	$(APM_DYN)/apm_salt_emit.o \
	$(APM_DYN)/apm_dust_emit.o \
	$(APM_DYN)/apm_put_emit.o \
	$(APM_DYN)/new-apm_h_adv.o \
	$(APM_DYN)/apm_3d_adv.o \
	$(APM_DYN)/apm_v_adv.o \
	$(APM_DYN)/new-apm_h_dif.o \
	$(APM_DYN)/apm_v_dif.o \
	$(APM_DYN)/apm_cld_convect.o \
	$(APM_DYN)/sub_wetdep_apm.o \
	$(APM_DYN)/sub_drydep_apm.o \
	$(APM_DYN)/sub_gradep_apm.o \
	$(APM_DYN)/apm_wet_dep.o \
	$(APM_DYN)/v2_apm_wet_dep.o \
	$(APM_DYN)/v3_apm_wet_dep.o \
	$(APM_DYN)/apm_gra_vel.o \
	$(APM_DYN)/apm_gra_dep.o \
	$(APM_DYN)/v2_apm_gra_dep.o \
	$(APM_DYN)/v3_apm_gra_dep.o \
	$(APM_DYN)/apm_dry_dep.o \
	$(APM_DYN)/v1.0-apm_dry_dep.o \
	$(APM_DYN)/sub_dry_sink.o \
	$(APM_DYN)/v2_sub_drydep_apm.o \
	$(APM_DYN)/v2_apm_dry_dep.o \
	$(APM_DYN)/v3_apm_dry_dep.o \
	$(APM_DYN)/apm_org_trsf.o \
	$(APM_DYN)/apm_acid.o \
	$(APM_DYN)/apm_h2so4_cbmz.o \
	$(APM_DYN)/apm_pso4_so2.o \
	$(APM_DYN)/apm_h2so4_gas.o \
	$(APM_DYN)/apm_so4_00dis.o \
	$(APM_DYN)/apm_isor_so4.o \
	$(APM_DYN)/apm_isor_tracer.o \
	$(APM_DYN)/sub_cal_lv.o \
	$(APM_DYN)/apm_clean_sp.o \
	$(APM_DYN)/apm_boxphy_driver.o \
	$(APM_DYN)/test_print.o \
	$(APM_DYN)/sub_time.o \
	$(APM_DYN)/print_gas.o \
	$(APM_DYN)/keep_dust_zero.o \
	$(APM_DYN)/apm_tracer_check.o \
	$(APM_DYN)/sub_check_value.o \
	$(APM_DYN)/sub_check_rgf.o \
	$(APM_DYN)/sub_check_matrix.o \
	$(APM_DYN)/force_rgf.o \
	$(APM_DYN)/apm_cal_opt.o \
	$(naq_moddir)/get_funitnaqpms.o \
	$(naq_moddir)/naqpms_ic.o \
	$(naq_moddir)/wr_naq_init.o \
	$(naq_moddir)/rd_global.o \
	$(naq_moddir)/rd_met.o \
	$(naq_dyndir)/hadvppm.o \
	$(naq_dyndir)/tridiag_zadv.o \
	$(naq_dyndir)/vrtslv.o \
	$(naq_dyndir)/calc_zrates.o \
	$(naq_dyndir)/naqpms_trop_pbltop.o \
	$(naq_dyndir)/naqpms_dms_emit.o \
	$(naq_dyndir)/naqpms_rd_emit.o \
	$(naq_dyndir)/naqpms_ppb2ugom3.o \
	$(naq_dyndir)/naqpms_put_emit_v2.o \
	$(naq_dyndir)/naqpms_putdust.o \
	$(naq_dyndir)/naqpms_putsalt.o \
	$(naq_dyndir)/naqpms_h_adv.o \
	$(naq_dyndir)/naqpms_v_adv.o \
	$(naq_dyndir)/advec1d_walcek.o \
	$(naq_dyndir)/naqpms_3d_adv.o \
	$(naq_dyndir)/module_gtiming.o \
	$(naq_dyndir)/naqpms_cld_convect.o \
	$(naq_dyndir)/naqpms_h_dif.o \
	$(naq_dyndir)/naqpms_gra_dep.o \
	$(naq_dyndir)/naqpms_dry_dep.o \
	$(naq_dyndir)/naqpms_drv_gaschem.o \
	$(naq_dyndir)/naqpms_bcoc_agt.o \
	$(naq_dyndir)/naqpms_wet_dep.o \
	$(naq_dyndir)/naqpms_drv_isorropia.adj.o \
	$(naq_dyndir)/naqpms_cal_diag.o \
	$(naq_dyndir)/naqpms_drv_soap.o \
	$(naq_dyndir)/naqpms_v_dif.o \
	$(naq_dyndir)/naqpms_output.o \
	$(naq_dyndir)/naqpms_aveload.o \
	$(naq_dyndir)/naqpms_output_v2.o \
	$(naq_dyndir)/naqpms_out_wet.o \
	$(naq_dyndir)/naqpms_out_dry.o \
	$(naq_dyndir)/naqpms_out_dust.o \
	$(naq_dyndir)/naqpms_out_salt.o \
	$(naq_dyndir)/naqpms_out_smark.o \
	$(naq_dyndir)/naqpms_out_term.o \
	$(naq_dyndir)/check_naqpms_tracer.o \
	$(MODEL)/cal_hgfac.o \
	$(MODEL)/caldate.o \
	$(MODEL)/juldate.o \
	$(MODEL)/clddiag.o \
	$(MODEL)/getznth.o \
	$(MODEL)/micromet_el.o \
	$(MODEL)/vd_gas_zhang.o \
	$(MODEL)/vd_aer_zhang.o \
	$(MODEL)/getland.o \
	$(MODEL)/getland_usgs2zhang.o \
	./main/geatm_vartype.o \
	$(MODEL)/parall.o   \
	$(MODEL)/inout.o \
	$(MODEL)/inout3.o \
	$(MODEL)/model.o \
	$(MODEL)/sub_rd_gridinfo.o \
	$(MODEL)/sub_rd_monthly_oxidants.o \
	$(MODEL)/adv_hori.o \
	$(MODEL)/adv_hori_ds.o \
	$(MODEL)/adv_vert2.o \
	$(MODEL)/dif_hori.o \
	$(MODEL)/DryVel.o \
	$(MODEL)/setboundary.o \
	$(MODEL)/setwindbound.o \
	$(MODEL)/setboundtop1.o \
	$(MODEL)/termbal.o   \
	$(MODEL)/termballi.o   \
	$(MODEL)/getboundary.o\
	$(MODEL)/clddiagf.o\
	$(MODEL)/convmix.o\
	$(MODEL)/strtmix.o\
	$(MODEL)/diffus.o\
	$(MODEL)/diffus_ds.o\
	$(MODEL)/vdiffimp.o\
	$(MODEL)/vdiffimp_ds.o\
	$(MODEL)/trdiag.o\
	$(MODEL)/trdiag_ds.o\
	$(MODEL)/tropause.o\
	$(MODEL)/calfconv.o\
	$(MODEL)/calmol.o\
	$(MODEL)/eddyz.o\
	$(MODEL)/convect43c.o\
	$(MODEL)/drydep_gas.o\
	$(MODEL)/drydep_aer.o\
	$(MODEL)/drydep_aer_old.o\
	$(MODEL)/drydep_dust.o\
	$(MODEL)/GETDU.o\
	$(MODEL)/MODIS.o\
	$(MODEL)/ALLD.o\
	$(LIQ_SRC)/GETCLOUDDEPTH.o\
	$(LIQ_SRC)/AQUEOUS.o\
	$(LIQ_SRC)/raqchem.o\
	$(LIQ_SRC)/hlconst.o\
	$(LIQ_SRC)/hlindex.o\
	$(LIQ_SRC)/WETDEP_GAS.o\
	$(LIQ_SRC)/scavrat.o\
	$(LIQ_SRC)/WETDEP_AER.o\
	$(LIQ_SRC)/WETDEP_AER_old.o\
	$(MARK)/diffus_mark.o\
	$(MARK)/vdiffimp_mark.o\
	$(MARK)/trdiag_mark.o\
	$(MARK)/clddiagf_mark.o\
	$(MARK)/convmix_mark.o\
	$(MARK)/strtmix_mark.o\
	$(MARK)/inout2.o \
	$(MARK)/getSMarkChange.o \
	$(MARK)/adv_hori_mark2.o \
	$(MARK)/adv_vert2_mark.o \
	$(MARK)/dif_hori_mark.o \
	$(MARK)/chemprod.o\
	$(GAS_SRC)/Jday.o \
	$(GAS_SRC)/CALCLF.o\
	$(GAS_SRC)/CALCLDOPD.o\
	$(GAS_SRC)/cbmz.o\
	$(GAS_SRC)/chemope.o\
	$(GAS_SRC)/Hetero.o\
	$(GAS_SRC)/Feevolution.o\
	$(AEROSOL_SRC)/ISOCOM.o\
	$(AEROSOL_SRC)/ISOFWD.o\
	$(AEROSOL_SRC)/ISOREV.o\
	$(AEROSOL_SRC)/ISRPINTR.o\
	$(AEROSOL_SRC)/SOAP.o\
	$(AEROSOL_SRC)/EXT.o\
	$(AEROSOL_SRC)/GETVIS.o \
	./main/geatm_comp_mct.o


$(LIBTARGET) :  $(MODULES)
	$(AR) $(ARFLAGS) $(LIBTARGET) *.o $(OBJS)

.f90.o :
	$(FC) $(FCBASEOPTS) -I$(NETCDF)/include -I$(OBJROOT)/mct/mct -I$(OBJROOT)/csm_share -I$(INCS)\
        -I$(OBJROOT)/lib/include -I./ -L$(NETCDF)/lib -L$(LIBROOT) -lmct -lnetcdf -lnetcdff -lpio -lcsm_share \
        -L$(OBJROOT)/mct/mpeu -lmpeu  -c -module ./ $*.f90      
        
.f.o :
	$(FC)  $(FCBASEOPTS) -I$(NETCDF)/include -I$(OBJROOT)/mct/mct -I$(OBJROOT)/csm_share -I$(INCS)\
        -I$(OBJROOT)/lib/include -I./ -L$(NETCDF)/lib -L$(LIBROOT) -lmct -lnetcdf -lnetcdff -lpio -lcsm_share \
        -L$(OBJROOT)/mct/mpeu -lmpeu -c -module ./ $*.f

.F90.o :
	$(FC)  $(FCBASEOPTS) -I$(NETCDF)/include -I$(OBJROOT)/mct/mct -I$(OBJROOT)/csm_share -I$(INCS)\
        -I$(OBJROOT)/lib/include -I./ -L$(NETCDF)/lib -L$(LIBROOT) -lmct -lnetcdf -lnetcdff -lpio -lcsm_share \
        -L$(OBJROOT)/mct/mpeu -lmpeu -c -module ./ $*.F90

clean:
	rm -f *.o *.mod

.PHONY: geatm clean
		
ALLD.o :
AQUEOUS.o :
CALCLDOPD.o :
CALCLF.o :
CALCLF2.o :
DryVel.o :
EXT.o :
Feevolution.o :
GETCLOUDDEPTH.o :
GETDU.o :
GETVIS.o :
Hetero.o :
ISOCOM.o :
ISOFWD.o :
ISOREV.o :
ISRPINTR.o :
Jday.o :
MODIS.o :
SOAP.o :
WETDEP_AER.o :
WETDEP_AER_old.o :
WETDEP_GAS.o :
adv1d_2016.o :
adv1d_comv.o :
adv_hori.o :
adv_hori_ds.o :
adv_hori_mark.o :
adv_hori_mark2.o :
adv_vert.o :
adv_vert2.o :
adv_vert2_ds.o :
adv_vert2_mark.o :
adv_vert_mark.o :
advec1d_walcek.o :
aerv2_ch_bdy.o :
aerv2_receive_ip.o :
aerv2_receive_jp.o :
aerv2_send_bdy.o :
allosea.o :
apm_3d_adv.o : met_fields.o adv1d_comv.o
apm_5bdy.o :
apm_acid.o :
apm_albd_mod.o :
apm_area.o :
apm_bcoc_sulf_emit.o :
apm_bins2bulk.o :
apm_boxphy_driver.ageing_v1.o :
apm_boxphy_driver.o :
apm_bulk2bins.o :
apm_cal_opt.o :
apm_ch_bdy.o :
apm_cld_convect.o :
apm_clean_sp.o :
apm_coated_species.o :
apm_dry_dep.o :
apm_dust_emit.o :
apm_eff_rad.o :
apm_eff_rad_sub.o :
apm_emit.o :
apm_gra_dep.o :
apm_gra_vel.o :
apm_h2so4_cbmz.o :
apm_h2so4_gas.o :
apm_h_adv.o :
apm_h_dif.o :
apm_init_parm.o :
apm_inout-aq.o :
apm_inout.o :
apm_isor_so4.o :
apm_isor_tracer.o :
apm_opti_mod.o :
apm_org_trsf.o :
apm_polar_trspt_receive.o :
apm_polar_trspt_send.o :
apm_pso4_so2.o :
apm_put_emit-erro.o :
apm_put_emit.o :
apm_receive_ip.o :
apm_receive_jp.o :
apm_rrtmg_sw.o :
apm_rst.o :
apm_salt_emit.o :
apm_send_bdy.o :
apm_snstvy_emit.o :
apm_so4_00dis.o :
apm_tracer_check.o : cputime.o
apm_v_adv.o :
apm_v_dif.o :
apm_var.o :
apm_wet_dep.o :
apm_zero_emit.o :
aqchem_driver.o :
aqchem_var.o :
back-aqchem_driver.o :
back_wetdep.o :
bcoc_bin_info.o :
box_apm_physics.o :
box_aqchem.o :
cal_hgfac.o : met_fields.o 
calc_air_density.o : met_fields.o 
calc_zrates.o : met_fields.o 
calfconv.o :
call.o :
calpv.o :
caltheta.o :
check_naqpms_tracer.o :
chemope.o :
chemprod.o :
clddiag.o : 
cnn.o :
cputime.o :
dbl-apm_var.o :
deallocate_inc.o :
dep_parm.o :
dif_hori.o :
dif_hori_ds.o :
dif_hori_mark.o :
diffus.o : vdiffimp.o
diffus_ds.o : vdiffimp_ds.o
vdiffimp.o : trdiag.o
vdiffimp_ds.o : trdiag_ds.o
dms.o :
drv.o : 
drydep_aer.o : dep_parm.o 
drydep_aer_old.o : dep_parm.o 
drydep_dust.o : dep_parm.o 
drydep_gas.o : dep_parm.o 
eddyz.o :
em4-apm_bcoc_sulf_emit.o :
em5-apm_bcoc_sulf_emit.o :
fine-aqchem_driver.o :
flags4wdep.o :
force_rgf.o :
getSMarkChange.o :
get_funit.o :
get_funitnaqpms.o :
getboundary.o :
getland.o :
getland_usgs2zhang.o :
h2o2.o :
inert.naqpms_drv_isorropia.o : met_fields.o 
initialize_apm.o : 
inout.o : 
inout2.o :
inout3.o :
isorropia_var.o : 
keep_dust_zero.o :
kvcalc_cmaq.o :
kvcalc_myj.o :
kvcalc_ob70.o : 
kvcalc_ysu.o :   
lijie-AQUEOUS.o :   
met_fields.o : 
model.o : 
module_gtiming.o :
naqpms_aveload.o : met_fields.o 
naqpms_bcoc_agt.o :   
naqpms_cal_diag.o : met_fields.o 
naqpms_cld_convect.o : met_fields.o 
naqpms_dms_emit.o : met_fields.o 
naqpms_drv_gaschem.o : ork_vars.o met_fields.o
naqpms_drv_isorropia.adj.o : met_fields.o 
naqpms_drv_isorropia.o : met_fields.o 
naqpms_drv_soap.o : met_fields.o 
naqpms_dry_dep.o : met_fields.o 
naqpms_gra_dep.o :   
naqpms_grid.o : 
naqpms_h_adv.o : met_fields.o 
naqpms_ic.o : work_vars.o met_fields.o 
naqpms_out_dry.o :   
naqpms_out_dust.o : met_fields.o 
naqpms_out_salt.o : met_fields.o 
naqpms_out_smark.o : met_fields.o 
naqpms_out_term.o :   
naqpms_out_wet.o : met_fields.o
naqpms_output.o : met_fields.o 
naqpms_output_v2.o : met_fields.o 
naqpms_ppb2ugom3.o : met_fields.o 
naqpms_put_emit.o :   
naqpms_put_emit_v2.o : met_fields.o 
naqpms_putdust.o : work_vars.o met_fields.o 
naqpms_putsalt.o : work_vars.o met_fields.o 
naqpms_rd_emit.o :   
naqpms_trop_pbltop.o : met_fields.o 
naqpms_v_adv.o : met_fields.o 
naqpms_var.o :   
naqpms_var_idx.o : 
naqpms_wet_dep.o : 
naqpms_v_dif.o : diffus.o diffus_ds.o
naqpms_h_dif.o :dif_hori.o dif_hori_mark.o
naqpms_3d_adv.o : advec1d_walcek.o hadvppm.o vrtslv.o 
hadvppm.o :
vrtslv.o : tridiag_zadv.o get_funitnaqpms.o
new-apm_h_adv.o : 
new-apm_h_dif.o : 
note.o : 
parall.o : 
photolysis_cld_adjust.o : 
print_gas.o : 
print_time.o : 
rd_global.o : 
rd_met.o : met_fields.o 
reaction.o :  
setboundary.o :  
setboundtop1.o :  
setwindbound.o :   
smp_sulfchem_var.o : 
sola_apm_boxphy_driver.o : 
sub_apm_ddepmass.o : 
sub_cal_cosz.o : 
sub_cal_inorganic.o : 
sub_cal_lv.o : 
sub_cal_yspgf.o : met_fields.o 
sub_check_isrp.o : 
sub_check_matrix.o : sub_check_matrix.f90  
sub_check_rgf.o : sub_check_rgf.f90  
sub_check_value.o : sub_check_value.f90  
sub_dry_sink.o : sub_dry_sink.f90  
sub_drydep_apm.o :
sub_get_apmsulf.o : 
sub_get_cldflag.o :
sub_gradep_apm.o : 
sub_julianday.o : 
sub_max_ij.o :
sub_naqpms_ddepmass.o :
sub_naqpms_emitsen.o :
sub_rd_gridinfo.o :
sub_rd_monthly_oxidants.o :
sub_time.o :
sub_wetdep_apm.o : 
sulfur_smp_gaschem.o :
termbal.o :
termballi.o :
test_apm_io4nucl.o :
test_print.o :
tropause.o :
v1.0-apm_dry_dep.o :
v1.naqpms_rd_emit.o :
v2_apm_dry_dep.o :
v2_apm_gra_dep.o :
v2_apm_inout.o :
v2_apm_wet_dep.o :
v2_sub_drydep_apm.o : dep_parm.o
v3_apm_dry_dep.o :
v3_apm_gra_dep.o :
v3_apm_inout.o :
v3_apm_wet_dep.o :
var_description.o :
wet-apm_dry_dep.o :
work_vars.o :
wr_naq_init.o : wr_naq_init.o work_vars.o met_fields.o
wrf_micromet.o :
geatm_vartype.o :
geatm_comp_mct.o : naqpms_var.o \
	smp_sulfchem_var.o \
	met_fields.o \
	work_vars.o \
	naqpms_grid.o \
	cputime.o \
	flags4wdep.o \
	adv1d_comv.o \
	sulfur_smp_gaschem.o \
	dep_parm.o \
	wrf_micromet.o \
	kvcalc_ysu.o \
	kvcalc_myj.o \
	kvcalc_cmaq.o \
	kvcalc_ob70.o \
	sub_max_ij.o \
	print_time.o \
	get_funit.o \
	apm_opti_mod.o \
	apm_nucl_mod.o \
	apm_init_mod.o \
	apm_coag_mod.o \
	apm_grow_mod.o \
	apm_phys_mod.o \
	apm_var.o \
	aqchem_var.o \
	sub_get_cldflag.o \
	aqchem_driver.o \
	sub_cal_inorganic.o \
	sub_cal_yspgf.o \
	v3_apm_inout.o \
	test_apm_io4nucl.o \
	apm_inout-aq.o \
	initialize_apm.o \
	apm_rst.o \
	apm_5bdy.o \
	apm_snstvy_emit.o \
	apm_zero_emit.o \
	apm_bcoc_sulf_emit.o \
	apm_salt_emit.o \
	apm_dust_emit.o \
	apm_put_emit.o \
	new-apm_h_adv.o \
	apm_3d_adv.o \
	apm_v_adv.o \
	new-apm_h_dif.o \
	apm_v_dif.o \
	apm_cld_convect.o \
	sub_wetdep_apm.o \
	sub_drydep_apm.o \
	sub_gradep_apm.o \
	apm_wet_dep.o \
	v2_apm_wet_dep.o \
	v3_apm_wet_dep.o \
	apm_gra_vel.o \
	apm_gra_dep.o \
	v2_apm_gra_dep.o \
	v3_apm_gra_dep.o \
	apm_dry_dep.o \
	v1.0-apm_dry_dep.o \
	sub_dry_sink.o \
	v2_sub_drydep_apm.o \
	v2_apm_dry_dep.o \
	v3_apm_dry_dep.o \
	apm_org_trsf.o \
	apm_acid.o \
	apm_h2so4_cbmz.o \
	apm_pso4_so2.o \
	apm_h2so4_gas.o \
	apm_so4_00dis.o \
	apm_isor_so4.o \
	apm_isor_tracer.o \
	sub_cal_lv.o \
	apm_clean_sp.o \
	apm_boxphy_driver.o \
	test_print.o \
	sub_time.o \
	print_gas.o \
	keep_dust_zero.o \
	apm_tracer_check.o \
	sub_check_value.o \
	sub_check_rgf.o \
	sub_check_matrix.o \
	force_rgf.o \
	apm_cal_opt.o \
	get_funitnaqpms.o \
	naqpms_ic.o \
	wr_naq_init.o \
	rd_global.o \
	rd_met.o \
	hadvppm.o \
	tridiag_zadv.o \
	vrtslv.o \
	calc_zrates.o \
	module_gtiming.o \
	naqpms_trop_pbltop.o \
	naqpms_dms_emit.o \
	naqpms_rd_emit.o \
	naqpms_ppb2ugom3.o \
	naqpms_put_emit_v2.o \
	naqpms_putdust.o \
	naqpms_putsalt.o \
	naqpms_h_adv.o \
	naqpms_v_adv.o \
	advec1d_walcek.o \
	naqpms_3d_adv.o \
	naqpms_cld_convect.o \
	naqpms_h_dif.o \
	naqpms_gra_dep.o \
	naqpms_dry_dep.o \
	naqpms_drv_gaschem.o \
	naqpms_bcoc_agt.o \
	naqpms_wet_dep.o \
	naqpms_drv_isorropia.adj.o \
	naqpms_cal_diag.o \
	naqpms_drv_soap.o \
	naqpms_v_dif.o \
	naqpms_output.o \
	naqpms_aveload.o \
	naqpms_output_v2.o \
	naqpms_out_wet.o \
	naqpms_out_dry.o \
	naqpms_out_dust.o \
	naqpms_out_salt.o \
	naqpms_out_smark.o \
	naqpms_out_term.o \
	check_naqpms_tracer.o \
	cal_hgfac.o \
	caldate.o \
	juldate.o \
	clddiag.o \
	getznth.o \
	micromet_el.o \
	vd_gas_zhang.o \
	vd_aer_zhang.o \
	getland.o \
	getland_usgs2zhang.o \
	geatm_vartype.o \
	parall.o   \
	inout.o \
	model.o \
	sub_rd_gridinfo.o \
	sub_rd_monthly_oxidants.o \
	adv_hori.o \
	adv_hori_ds.o \
	adv_vert2.o \
	dif_hori.o \
	DryVel.o \
	setboundary.o \
	setwindbound.o \
	setboundtop1.o \
	termbal.o   \
	termballi.o   \
	getboundary.o\
	clddiagf.o\
	convmix.o\
	strtmix.o\
	diffus.o\
	diffus_ds.o\
	vdiffimp.o\
	vdiffimp_ds.o\
	trdiag.o\
	trdiag_ds.o\
	tropause.o\
	calfconv.o\
	calmol.o\
	eddyz.o\
	convect43c.o\
	drydep_gas.o\
	drydep_aer.o\
	drydep_aer_old.o\
	drydep_dust.o\
	GETDU.o\
	MODIS.o\
	ALLD.o\
	GETCLOUDDEPTH.o\
	AQUEOUS.o\
	raqchem.o\
	hlconst.o\
	hlindex.o\
	WETDEP_GAS.o\
	scavrat.o\
	WETDEP_AER.o\
	WETDEP_AER_old.o\
	diffus_mark.o\
	vdiffimp_mark.o\
	trdiag_mark.o\
	clddiagf_mark.o\
	convmix_mark.o\
	strtmix_mark.o\
	inout2.o \
	inout3.o \
	getSMarkChange.o \
	adv_hori_mark2.o \
	adv_vert2_mark.o \
	dif_hori_mark.o \
	chemprod.o\
	Jday.o \
	CALCLF.o\
	CALCLDOPD.o\
	cbmz.o\
	chemope.o\
	Hetero.o\
	Feevolution.o\
	ISOCOM.o\
	ISOFWD.o\
	ISOREV.o\
	ISRPINTR.o\
	SOAP.o\
	EXT.o\
	GETVIS.o
