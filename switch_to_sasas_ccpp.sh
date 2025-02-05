#!/bin/bash

set -e -x

export BASEDIR=$(pwd)

# for file scripts/ccsm_utils/Machines/Macros.derecho
cd $BASEDIR/scripts/ccsm_utils/Machines/
sed -i 's/CPPDEFS += -DCCPP$/CPPDEFS += -DCCPP -DCCPP_SASAS/g' Macros.derecho
sed -i 's/#CPPDEFS += -DCCPP$/CPPDEFS += -DCCPP -DCCPP_SASAS/g' Macros.derecho
# grep -e -c 'CPPDEFS += -DCCPP -DCCPP_SASAS' Macros.derecho

cd $BASEDIR/models/atm/iap
./ccpp/framework/scripts/ccpp_prebuild.py --config ./ccpp/config/ccpp_prebuild_config.py --debug --verbose 2>&1 | tee ccpp_prebuild.log

rsync -av bld/ccpp/physics/ src/dynamics/physics/
rm -rf bld/ccpp/physics
cp -av ccpp/framework/src/ccpp_types.F90 src/dynamics

mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_static_api.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_sasasr_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_test1_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/
mv -v $BASEDIR/models/atm/cam/src/dynamics/iap/physics/ccpp_IAP_test_test2_cap.F90 $BASEDIR/models/atm/cam/src/dynamics/iap/


cd $BASEDIR/models/atm/iap/src/physics
ln -sf ../../ccpp/physics/physics/sascnvnr.F .
ln -sf ../../ccpp/physics/physics/sascnvnr.meta .
ln -sf ../../ccpp/physics/physics/funcphys.f90 .
ln -sf ../../ccpp/physics/physics/physcons.F90 .
ln -sf ../../ccpp/physics/physics/physcons.F90 .
ln -sf ../../ccpp/physics/physics/machine.F .

cd $BASEDIR/models/atm/iap/ccpp/physics/physics/
cp -av  cldwat_ccpp.F90 iap_* zm_conv*F90 ../../../src/physics/
# cd $BASEDIR
# # rm -vf models/atm/cam/src/dynamics/iap/ccpp_static_api.F90
# rm -vf models/atm/cam/src/dynamics/iap/ccpp_types.F90
# # rm -vf models/atm/cam/src/dynamics/iap/ccpp_IAP_test_cap.F90
# rm -vf models/atm/cam/src/dynamics/iap/ccpp_data.F90
# # rm -vf models/atm/cam/src/dynamics/iap/ccpp_IAP_test_test1_cap.F90

echo "Now go to scripts/ directory to setup and build case, see README"
# cd $BASEDIR/scripts/FAMIPC5_FD14
# ./FAMIPC5_FD14.derecho.clean_build
# ./FAMIPC5_FD14.derecho.build
