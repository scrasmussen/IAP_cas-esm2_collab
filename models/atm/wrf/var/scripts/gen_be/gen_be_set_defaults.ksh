#!/bin/ksh
#-----------------------------------------------------------------------
# Script: gen_be_set_defaults.ksh
#
# Purpose: This scripts sets the environment variables used within the 
# entire scripts system to values corresponding to a "standard case".
# The standard case currently used is the con200 application.
# The namelist parameters specified here is that sub-set of the entire
# range of parameters for all namelists that we have found necessary
# to change through experience of the applications tested so far. As
# new applications are tests, additional environment valiables may
# be added. 

#-----------------------------------------------------------------------
# [1] Set defaults for required environment variables:
#-----------------------------------------------------------------------

# Decide which stages to run (run if true):

export RUN_GEN_BE_STAGE0=${RUN_GEN_BE_STAGE0:-false} # Run stage 0 (create perturbation files).
export RUN_GEN_BE_STAGE1=${RUN_GEN_BE_STAGE1:-false} # Run stage 1 (Remove mean, split variables).
export RUN_GEN_BE_STAGE2=${RUN_GEN_BE_STAGE2:-false} # Run stage 2 (Regression coefficients).
export RUN_GEN_BE_STAGE2A=${RUN_GEN_BE_STAGE2A:-false} # Run stage 2a (Regression coefficients).
export RUN_GEN_BE_STAGE3=${RUN_GEN_BE_STAGE3:-false} # Run stage 3 (Vertical Covariances).
export RUN_GEN_BE_STAGE4=${RUN_GEN_BE_STAGE4:-false} # Run stage 4 (Horizontal Covariances).
export RUN_GEN_BE_DIAGS=${RUN_GEN_BE_DIAGS:-false}   # Run gen_be diagnostics.
export RUN_GEN_BE_DIAGS_READ=${RUN_GEN_BE_DIAGS_READ:-false}   # Run gen_be diagnostics_read.
export RUN_GEN_BE_MULTICOV=${RUN_GEN_BE_MULTICOV:-false} # Set to calculate chi/T/ps regression diags. 

export DOMAIN=${DOMAIN:-01}                          # domain id.
export START_DATE=${START_DATE:-2003010200}          # Time of first perturbation.
export END_DATE=${END_DATE:-2003012812}              # Time of last perturbation.
export FCST_RANGE=${FCST_RANGE:-24}                  # Forecast range of forecast (hours).
export FCST_RANGE1=${FCST_RANGE1:-24}                # Forecast range of forecast 1 (hours).
export FCST_RANGE2=${FCST_RANGE2:-12}                # Forecast range of forecast 2 (hours).
export INTERVAL=${INTERVAL:-12}                      # Period between files (hours).
export BE_METHOD=${BE_METHOD:-NMC}                   # NMC (NMC-method), ENS (Ensemble-Method).
export NE=${NE:-1}                                   # Number of ensemble members (for ENS).
export BIN_TYPE=${BIN_TYPE:-5}                       # 0=None, 1=1:ni, 2=latitude, ....
export LAT_MIN=${LAT_MIN:--90.0}                     # Used if BIN_TYPE = 2.
export LAT_MAX=${LAT_MAX:-90.0}                      # Used if BIN_TYPE = 2.
export BINWIDTH_LAT=${BINWIDTH_LAT:-10.0}            # Used if BIN_TYPE = 2.
export BINWIDTH_HGT=${BINWIDTH_HGT:-1000.0}          # Used if BIN_TYPE = 2.
export HGT_MIN=${HGT_MIN:-0.0}                       # Used if BIN_TYPE = 2.
export HGT_MAX=${HGT_MAX:-20000.0}                   # Used if BIN_TYPE = 2.
export REMOVE_MEAN=${REMOVE_MEAN:-.true.}            # Remove time/ensemble/area mean.
export GAUSSIAN_LATS=${GAUSSIAN_LATS:-.false.}       # Set if Gaussian latitudes used (global only). 
export TESTING_EOFS=${TESTING_EOFS:-.true.}          # True if performing EOF tests.
export NUM_PASSES=${NUM_PASSES:-0}                   # Number of passes of recursive filter.
export RF_SCALE=${RF_SCALE:-1.0}                     # Recursive filter scale.
export USE_GLOBAL_EOFS=${USE_GLOBAL_EOFS:-.true.}    # Use domain-averaged EOFS for stage3.
export DATA_ON_LEVELS=${DATA_ON_LEVELS:-.false.}     # False if fields projected onto modes.
export GLOBAL=${GLOBAL:-false}                       # Global or regional models
export NUM_LEVELS=${NUM_LEVELS:-27}                  # Hard-wired for now....
export N_SMTH_SL=${N_SMTH_SL:-2}                     # Amount of lengthscale smoothing (0=none).
export STRIDE=${STRIDE:-1}                           # Calculate correlation evert STRIDE point (stage4 regional).
export NBINS=${NBINS:-1}                             # Number of latitude bins for length scale computation
export IBIN=${IBIN:-1}                               # Index of latitude bin to compute length scale for
export TESTING_SPECTRAL=${TESTING_SPECTRAL:-.false.} # True if performing spectral tests.
export LOCAL=${LOCAL:-true}                          # True if local machine.
export NUM_JOBS=${NUM_JOBS:-1}                       # Number of jobs to run (stage4 regional)).
export MACHINES=${MACHINES:-" node1 node1 node2 node2 node3 node3 node4 node4 "\
                            " node5 node5 node6 node6 node7 node7 node8 node8"}
export REGION=${REGION:-con200}
export EXPT=${EXPT:-noobs}
export ID=${ID:-gen_be}
export ID1=${ID1:-${BE_METHOD}.bin_type${BIN_TYPE}}
export VARIABLE1=${VARIABLE1:-chi_u}              # For cov3d
export VARIABLE2=${VARIABLE2:-chi}                # For cov3d
export CLEAN=false
export FILE_TYPE=${FILE_TYPE:-wrfout}

# Directories:
export REL_DIR=${REL_DIR:-$HOME/trunk}            # Directory containing codes.
export WRFVAR_DIR=${WRFVAR_DIR:-$REL_DIR/wrfvar}  # WRF-Var code directory.
export BUILD_DIR=${BUILD_DIR:-$WRFVAR_DIR/var/da} # WRF-Var code build directory.
export DAT_DIR=${DAT_DIR:-${HOME}/data}           # Top-level data directory.
export REG_DIR=${REG_DIR:-$DAT_DIR/$REGION}       # Region-specific data dir.
export EXP_DIR=${EXP_DIR:-$REG_DIR/$EXPT}         # Experiment-specific data dir.
export FC_DIR=${FC_DIR:-$EXP_DIR/fc}              # Forecast directory
export RUN_DIR=${RUN_DIR:-$EXP_DIR/gen_be$BIN_TYPE} # Run dir.
export WORK_DIR=${WORK_DIR:-$RUN_DIR/working}     # Working directory
export STAGE0_DIR=${STAGE0_DIR:-$WORK_DIR/stage0} # Output for stage0.

if $GLOBAL; then
   export UH_METHOD=power
else
   export UH_METHOD=scale
fi

export CONTROL_VARIABLES=${CONTROL_VARIABLES:-" psi chi_u t_u rh ps_u "}
export DELETE_DIRS=${DELETE_DIRS:-" "}

