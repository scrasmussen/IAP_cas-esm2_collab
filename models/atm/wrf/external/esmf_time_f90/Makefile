# To build this by itself, use the make target esmf_time_f90_only 
# from the top-level WRF Makefile.  
# > cd ../..
# configure
# make esmf_time_f90_only

.SUFFIXES: .F90 .o .f

# get rid of single quotes after comments
# WARNING:  This will break if a quoted string is followed by a comment that has
#           a single quote.
SED_FTN = sed -e "/\!.*'/s/'//g"

RM = /bin/rm -f
AR = ar
#RANLIB	= ranlib
RANLIB	= echo

OBJS =  \
	module_symbols_util.o \
	module_utility.o 
	
default: libesmf_time.a

tests: Test1_ESMF.exe Test1_WRFU.exe

libesmf_time.a : $(OBJS)
	$(RM) libesmf_time.a
	if [ "$(AR)" != "lib.exe" ] ; then \
          $(AR) ru libesmf_time.a $(OBJS) ; \
        else \
          $(AR) /out:libesmf_time.a $(OBJS) ; \
        fi
	$(RANLIB) libesmf_time.a

Test1_ESMF.f : Test1.F90
	$(RM) Test1_ESMF.b Test1_ESMF.f
	cp Test1.F90 Test1_ESMF.b
	$(CPP) -C -P -I. Test1_ESMF.b > Test1_ESMF.f

Test1_ESMF.exe : libesmf_time.a Test1_ESMF.o
	$(FC) -o Test1_ESMF.exe Test1_ESMF.o libesmf_time.a

Test1_WRFU.f : Test1.F90
	$(RM) Test1_WRFU.b Test1_WRFU.f
	sed -e "s/ESMF_Mod/module_utility/g" -e "s/ESMF_/WRFU_/g" Test1.F90 > Test1_WRFU.b
	$(CPP) -C -P -I. Test1_WRFU.b > Test1_WRFU.f

Test1_WRFU.exe : libesmf_time.a Test1_WRFU.o
	$(FC) -o Test1_WRFU.exe Test1_WRFU.o libesmf_time.a

.F90.o :
	$(RM) $@
	$(SED_FTN) $*.F90 > $*.b
	$(CPP)  -I$(LIBROOT)/include  -C -P -I. $*.b > $*.f
	$(RM) $*.b
	$(FC)  -I$(LIBROOT)/include -o $@ -c $*.f

.F90.f :
	$(RM) $@
	$(SED_FTN) $*.F90 > $*.b
	$(CPP) -C -P -I. $*.b > $*.f
	$(RM) $*.b

.f.o :
	$(RM) $@
	$(RM) $*.b
	$(FC) -c $*.f

clean : testclean

testclean: 
	$(RM) *.b *.f *.o *.obj *.i libesmf_time.a *.mod Test1*.exe

superclean: testclean
	$(RM) Test1*.out make_tests.out

# DEPENDENCIES : only dependencies after this line 

#$$$ update dependencies!  


module_utility.o :  module_symbols_util.o
Test1.o :  module_utility.o

