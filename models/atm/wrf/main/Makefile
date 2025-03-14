#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f

MODULES =  wrf_comp_mct.o

OBJS    = 

LIBPATHS =

include ../configure.wrf


$(SOLVER)_wrf : wrf.o ../main/wrf_comp_mct.o
	$(RANLIB) $(LIBWRFLIB)
	$(LD) -o wrf.exe $(LDFLAGS) wrf.o ../main/wrf_comp_mct.o $(LIBWRFLIB) $(LIB)

$(SOLVER)_wrf_SST_ESMF : wrf_ESMFMod.o wrf_SST_ESMF.o ../main/wrf_comp_mct.o
	$(RANLIB) $(LIBWRFLIB)
	$(LD) -o wrf_SST_ESMF.exe $(LDFLAGS) wrf_SST_ESMF.o wrf_ESMFMod.o ../main/wrf_comp_mct.o $(LIBWRFLIB) $(LIB)

$(SOLVER)_ideal : module_initialize ideal.o
	$(RANLIB) $(LIBWRFLIB)
	$(LD) -o ideal.exe $(LDFLAGS) ideal.o ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o $(LIBWRFLIB) $(LIB)

$(SOLVER)_real : module_initialize real_$(SOLVER).o ndown_$(SOLVER).o nup_$(SOLVER).o tc_$(SOLVER).o
	$(RANLIB) $(LIBWRFLIB)
	$(LD) -o ndown.exe $(LDFLAGS) ndown_$(SOLVER).o  ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o $(LIBWRFLIB) $(LIB)
	$(LD) -o nup.exe $(LDFLAGS) nup_$(SOLVER).o  ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o $(LIBWRFLIB) $(LIB)
	$(LD) -o tc.exe $(LDFLAGS) tc_$(SOLVER).o  ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o $(LIBWRFLIB) $(LIB)

convert_em : convert_em.o
	$(RANLIB) $(LIBWRFLIB)
	$(LD) -o convert_em.exe $(LDFLAGS) convert_em.o $(LIBWRFLIB) $(LIB)

convert_nmm : convert_nmm.o
	$(RANLIB) $(LIBWRFLIB)
	$(FC) -o convert_nmm.exe $(LDFLAGS) convert_nmm.o $(LIBWRFLIB) $(LIB)

real_nmm : real_nmm.o
	( cd ../dyn_nmm ;  $(MAKE) module_initialize_real.o )
	$(RANLIB) $(LIBWRFLIB)
	$(FC) -o real_nmm.exe $(LDFLAGS) real_nmm.o ../dyn_nmm/module_initialize_real.o ../share/module_optional_input.o ../share/input_wrf.o ../share/module_io_domain.o $(LIBWRFLIB) $(LIB)

module_initialize :
	( cd ../dyn_$(SOLVER) ;  $(MAKE) module_initialize_$(IDEAL_CASE).o )

## prevent real being compiled for OMP -- only for regtesting
#$(SOLVER)_real : module_initialize real_$(SOLVER).o
#	$(RANLIB) $(LIBWRFLIB)
#	if [ -z "$(OMP)" ] ; then $(FC) -o real.exe $(LDFLAGS) real_$(SOLVER).o ../dyn_$(SOLVER)/module_initialize_$(IDEAL_CASE).o $(LIBWRFLIB) $(LIB) ; fi
#
## prevent module_initialize being compiled for OMP --remove after IBM debugging
#module_initialize :
#	if [ -z "$(OMP)" ] ; then ( cd ../dyn_$(SOLVER) ;  $(MAKE) module_initialize_$(IDEAL_CASE).o ) ; fi
# end of regtest changes

clean:
	@ echo 'use the clean script'

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

convert_nmm.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	$(ESMF_MOD_DEPENDENCE)

convert_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	$(ESMF_MOD_DEPENDENCE)

ideal.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../share/module_io_domain.o \
	../dyn_$(SOLVER)/$(CASE_MODULE) \
	$(ESMF_MOD_DEPENDENCE)

ndown_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../frame/module_wrf_error.o \
	../frame/module_integrate.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	../share/module_get_file_names.o \
	../share/module_soil_pre.o \
	../dyn_em/module_initialize_$(IDEAL_CASE).o \
	../dyn_em/module_big_step_utilities_em.o \
	$(ESMF_MOD_DEPENDENCE)

nup_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_streams.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../frame/module_wrf_error.o \
	../frame/module_integrate.o \
	../share/module_bc.o \
	../share/module_io_domain.o \
	../share/module_get_file_names.o \
	../share/module_soil_pre.o \
	../dyn_em/module_initialize_real.o \
	../dyn_em/module_big_step_utilities_em.o \
	$(ESMF_MOD_DEPENDENCE)

# this already built above :../dyn_em/module_initialize.real.o 
tc_em.o: \
	../frame/module_machine.o \
	../frame/module_domain.o \
	../frame/module_driver_constants.o \
	../frame/module_configure.o \
	../frame/module_timing.o \
	../frame/module_dm.o \
	../dyn_em/module_big_step_utilities_em.o \
	../share/module_io_domain.o \
	../share/module_date_time.o \
	../share/module_optional_input.o \
	../share/module_bc_time_utilities.o \
	$(ESMF_MOD_DEPENDENCE)

wrf.o:  ../main/wrf_comp_mct.o

wrf_ESMFMod.o:  ../main/wrf_comp_mct.o

wrf_SST_ESMF.o:  wrf_ESMFMod.o

wrf_comp_mct.o: ../frame/module_machine.o \
		../frame/module_domain.o \
		../frame/module_integrate.o \
		../frame/module_driver_constants.o \
		../frame/module_configure.o \
		../frame/module_timing.o \
		../frame/module_wrf_error.o \
		../frame/module_atm_communicator.o \
		../frame/module_state_description.o\
		../frame/module_dm.o \
		../frame/module_comm_dm.o \
		../share/module_io_domain.o \
		../share/module_io_wrf.o \
		../share/module_date_time.o \
		../share/module_check_a_mundo.o \
		../frame/module_io.o \
		../share/module_llxy.o \
		../prep/module_real.o \
		../prep/module_geogrid.o \
		../prep/llxy_module.o \
		../prep/gridinfo_module.o \
		../prep/list_module.o \
		../prep/parallel_module.o \
		../prep/source_data_module.o \
		../prep/module_camdomain.o \
		../prep/module_geodomain.o \
		../frame/module_nesting.o \
		$(OBJROOT)/csm_share/seq_flds_indices.o \
		$(OBJROOT)/csm_share/seq_flds_mod.o \
		$(OBJROOT)/csm_share/seq_infodata_mod.o \
		$(OBJROOT)/csm_share/seq_timemgr_mod.o \
		$(OBJROOT)/csm_share/seq_cdata_mod.o \
		$(OBJROOT)/csm_share/seq_comm_mct.o \
		$(OBJROOT)/csm_share/shr_file_mod.o \
		$(OBJROOT)/csm_share/shr_sys_mod.o \
		$(OBJROOT)/csm_share/mct_mod.o \
		$(OBJROOT)/csm_share/esmf_mod.o \
		$(OBJROOT)/utils/timing/perf_mod.o \
		$(ESMF_MOD_DEPENDENCE) 

# DO NOT DELETE
