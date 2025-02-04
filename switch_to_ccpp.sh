#!/bin/bash

set -e -x

export BASEDIR=$(pwd)

cd $BASEDIR/scripts/ccsm_utils/Machines/
sed -i 's/#CPPDEFS += -DCCPP/CPPDEFS += -DCCPP/g' Macros.derecho
sed -i 's/CPPDEFS += -DCCPP -DCCPP_SASAS/CPPDEFS += -DCCPP/g' Macros.derecho
grep -e 'CPPDEFS += -DCCPP' Macros.derecho

cd $BASEDIR/models/atm/iap
./ccpp/framework/scripts/ccpp_prebuild.py --config ./ccpp/config/ccpp_prebuild_config.py --debug --verbose 2>&1 | tee ccpp_prebuild.log

rsync -av bld/ccpp/physics/ src/dynamics/physics/
rm -rf bld/ccpp/physics

mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_static_api.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_sasasr_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_test1_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_test2_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/

cp -av ccpp/framework/src/ccpp_types.F90 src/dynamics

cd $BASEDIR/models/atm/iap/ccpp/physics/physics/
# cp -av  cldwat_ccpp.F90 dummy_loop.F90 iap_* zm_conv*F90 ../../../src/physics/
cp -av  cldwat_ccpp.F90 iap_* zm_conv*F90 ../../../src/physics/

echo "Now go to scripts/ directory to setup and build case, see README"
# cd $BASEDIR/scripts/FAMIPC5_FD14
# ./FAMIPC5_FD14.derecho.clean_build
# ./FAMIPC5_FD14.derecho.build
