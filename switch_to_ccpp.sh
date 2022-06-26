#!/bin/bash

set -e

export BASEDIR=/glade/work/heinzell/iap/IAP-AGCM_CCPP

sed -i 's/#CPPDEFS += -DCCPP/CPPDEFS += -DCCPP/g' Macros.cheyenne
grep -e 'CPPDEFS += -DCCPP' Macros.cheyenne

cd $BASEDIR/IAP_cas-esm2/models/atm/iap
./ccpp/framework/scripts/ccpp_prebuild.py --config ./ccpp/config/ccpp_prebuild_config.py --debug --verbose 2>&1 | tee ccpp_prebuild.log

rsync -av bld/ccpp/physics/ src/dynamics/physics/
rm -fr bld/ccpp/physics

mv -v $BASEDIR/IAP_cas-esm2/models/atm/cam/src/dynamics/iap/physics/ccpp_static_api.F90 $BASEDIR/IAP_cas-esm2/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/IAP_cas-esm2/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_cap.F90 $BASEDIR/IAP_cas-esm2/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/IAP_cas-esm2/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_test1_cap.F90 $BASEDIR/IAP_cas-esm2/models/atm/cam/src/dynamics/iap/

cp -av ccpp/framework/src/ccpp_types.F90 src/dynamics

cd $BASEDIR/IAP_cas-esm2/models/atm/iap/ccpp/physics/physics/
cp -av  cldwat_ccpp.F90 dummy_loop.F90 iap_* zm_conv*F90 ../../../src/physics/

cd $BASEDIR/IAP_cas-esm2/scripts/FAMIPC5_FD14

./FAMIPC5_FD14.cheyenne.clean_build
./FAMIPC5_FD14.cheyenne.build
