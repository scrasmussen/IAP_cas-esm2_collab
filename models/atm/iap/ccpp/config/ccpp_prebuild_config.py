#!/usr/bin/env python

# CCPP prebuild config for IAP AGCM


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "IAP"

# Add all files with metadata tables on the host model side and in CCPP,
# relative to basedir = top-level directory of host model. This includes
# kind and type definitions used in CCPP physics. Also add any internal
# dependencies of these files to the list.
VARIABLE_DEFINITION_FILES = [
    # actual variable definition files
    'ccpp/framework/src/ccpp_types.F90',
    '../../csm_share/shr/shr_kind_mod.F90',
    'src/physics/physics_types.F90',
    'src/physics/ppgrid.F90',
    '../cam/src/dynamics/iap/pmgrid.F90',
    '../cam/src/control/hycoef.F90',
    '../cam/src/utils/spmd_utils.F90',
    '../cam/src/control/cam_logfile.F90',
    '../cam/src/control/camsrfexch_types.F90',
    '../cam/src/physics/cam/constituents.F90',
    '../cam/src/control/physconst.F90',
    '../cam/src/dynamics/iap/ccpp_data.F90',
    'src/physics/buffer.F90',
    ]

TYPEDEFS_NEW_METADATA = {
    'ccpp_types' : {
        'ccpp_types' : '',
        'ccpp_t' : 'cdata',
        },
    'ccpp_data' : {
        'ccpp_data' : ''
        },
    'physics_types' : {
        'physics_state' : 'phys_state(cdata%blk_no)',
        'physics_int_ephem' : 'phys_int_ephem(cdata%blk_no)',
        'physics_int_pers' : 'phys_int_pers(cdata%blk_no)',
        'physics_global' : 'phys_global',
        'physics_types' : '',
        },
     'camsrfexch_types' : {
        'cam_in_t' : 'cam_in(cdata%blk_no)',
        'camsrfexch_types' : '',
        },
    'buffer' : {
        'pblht(:,ccpp_block_number)'  : 'pblht(:,cdata%blk_no)',
        'tpert(:,ccpp_block_number)'  : 'tpert(:,cdata%blk_no)',
        'tpert2(:,ccpp_block_number)' : 'tpert2(:,cdata%blk_no)',
        'qpert2(:,ccpp_block_number)' : 'qpert2(:,cdata%blk_no)',
        'buffer' : '',
        },
    }

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'ccpp/physics/physics/zm_convr.F90'             ,
    'ccpp/physics/physics/zm_convr_post.F90'        ,
    'ccpp/physics/physics/zm_conv_evap.F90'         ,
    'ccpp/physics/physics/zm_conv_evap_post.F90'    ,
    'ccpp/physics/physics/zm_conv_momtran.F90'      ,
    'ccpp/physics/physics/zm_conv_momtran_post.F90' ,
    'ccpp/physics/physics/zm_conv_convtran.F90'     ,
    'ccpp/physics/physics/zm_conv_all_post.F90'     ,
    ]

# Default build dir, relative to current working directory,
# if not specified as command-line argument
DEFAULT_BUILD_DIR = 'bld'

# Auto-generated makefile/cmakefile snippets that contain all type definitions
TYPEDEFS_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.mk'
TYPEDEFS_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.cmake'
TYPEDEFS_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_TYPEDEFS.sh'

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_SCHEMES.cmake'
SCHEMES_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_SCHEMES.sh'

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE   = '{build_dir}/ccpp/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE  = '{build_dir}/ccpp/physics/CCPP_CAPS.cmake'
CAPS_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_CAPS.sh'

# Directory where to put all auto-generated physics caps
CAPS_DIR = '{build_dir}/ccpp/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'ccpp/suites'

# Optional arguments - only required for schemes that use
# optional arguments. ccpp_prebuild.py will throw an exception
# if it encounters a scheme subroutine with optional arguments
# if no entry is made here. Possible values are: 'all', 'none',
# or a list of standard_names: [ 'var1', 'var3' ].
OPTIONAL_ARGUMENTS = {
    #'subroutine_name_1' : 'all',
    #'subroutine_name_2' : 'none',
    #'subroutine_name_2' : [ 'var1', 'var3'],
    }

# Directory where to write static API to
STATIC_API_DIR = '{build_dir}/ccpp/physics'
STATIC_API_CMAKEFILE = '{build_dir}/ccpp/physics/CCPP_STATIC_API.cmake'
STATIC_API_SOURCEFILE = '{build_dir}/ccpp/physics/CCPP_STATIC_API.sh'

# Directory for writing HTML pages generated from metadata files
METADATA_HTML_OUTPUT_DIR = '{build_dir}/ccpp/physics/physics/docs'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = '{build_dir}/ccpp/physics/CCPP_VARIABLES_IAP.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = '{build_dir}/ccpp/framework/doc/DevelopersGuide/CCPP_VARIABLES_IAP.tex'
