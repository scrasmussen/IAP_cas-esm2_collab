#!/bin/bash

set -e

export BASEDIR=$(pwd)


sed -i 's/CPPDEFS += -DCCPP/#CPPDEFS += -DCCPP/g' Macros.derecho
grep -e '#CPPDEFS += -DCCPP' Macros.derecho

rm -vf models/atm/cam/src/dynamics/iap/ccpp_static_api.F90
rm -vf models/atm/cam/src/dynamics/iap/ccpp_types.F90
rm -vf models/atm/cam/src/dynamics/iap/ccpp_IAP_test_cap.F90
rm -vf models/atm/cam/src/dynamics/iap/ccpp_data.F90
rm -vf models/atm/cam/src/dynamics/iap/ccpp_IAP_test_test1_cap.F90

./FAMIPC5_FD14.derecho.clean_build
./FAMIPC5_FD14.derecho.build
