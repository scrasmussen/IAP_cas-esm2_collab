#!/bin/csh -f

#===============================================================================
#  ocn.iage.setup.csh : perform setup tasks for ideal age module
#
#  recognized commands, possibly with arguments, are
#    set_nt         nt_filename
#    namelist       licom_in_filename
#    set_tavg_nml   
#    tavg_contents  tavg_contents_filename my_stream
#    prestage       res_dpt_dir res_indpt_dir
#    document       licom_document_files
#    ccsm_prestage  ccsm_prestage_file
#
#  SVN:$Id: ocn.iage.setup.csh 22408 2010-04-15 20:53:19Z njn01 $
#
#===============================================================================

if ($#argv < 1) then
   echo ocn.iage.setup.csh : command argument missing
   exit 1
endif

#===============================================================================
#  set module name, which is required for tavg_nml
#  module name must match the name of this setup script
#===============================================================================
set module = iage

set command = $1

if ($command == set_nt) then

   echo ocn.iage.setup.csh : setting nt                                          >> $LICOM_BLDNML
   echo ------------------------------------------------------------------------ >> $LICOM_BLDNML

   if ($#argv < 2) then
      echo nt_filename argument missing
      exit 3
   endif

   set nt_filename = $2

   if !(-f $nt_filename) then
      echo nt_filename = $nt_filename does not exist
      exit 3
   endif

   @ nt_in = `cat $nt_filename`
   @ nt_in += 1
   echo $nt_in >! $nt_filename

else if ($command == namelist) then

   echo ocn.iage.setup.csh : setting namelist options                            >> $LICOM_BLDNML
   echo ------------------------------------------------------------------------ >> $LICOM_BLDNML

   if ($#argv < 2) then
      echo licom_in_filename argument missing
      exit 4
   endif

   set licom_in_filename = $2

   if !(-f $licom_in_filename) then
      echo licom_in_filename = $pop_in_filename does not exist
      exit 4
   endif

   set init_iage_option = ccsm_$runtype

   cat >> $licom_in_filename << EOF

&iage_nml
   init_iage_option = '$init_iage_option'
   init_iage_init_file = 'same_as_TS'
/
EOF

else if ($command == set_tavg_nml) then

  #-------------------------------------------------------------------------------------
  # if there is no module-related tavg output, set n_tavg_streams_tracer = 0
  #-------------------------------------------------------------------------------------
    set n_tavg_streams_tracer = 0
cat >&! $LICOM_DOCDIR/$module.tavg << EOF
n_tavg_streams_tracer =  $n_tavg_streams_tracer
EOF

else if ($command == tavg_contents) then

   echo ocn.iage.setup.csh : setting tavg_contents variables                     >> $LICOM_BLDNML
   echo ------------------------------------------------------------------------ >> $LICOM_BLDNML

   if ($#argv < 3) then
      echo tavg_contents_filename argument missing
      exit 5
   endif

   set tavg_contents_filename = $2

   if !(-f $tavg_contents_filename) then
      echo tavg_contents_filename = $tavg_contents_filename does not exist
      exit 5
   endif

   @ my_stream = $3
   if ($my_stream < 1) then
      echo invalid my_stream number  ($my_stream)
      exit 5
   endif

   #------------------------------------------------------------------------------------
   # For now, set streams manually. You must only set as many streams as are declared
   #  in the tavg_nml section. For example, if there are three streams:
   #  @ s1 = $my_stream
   #  @ s2 = $s1 + 1
   #  @ s3 = $s2 + 1
   #------------------------------------------------------------------------------------

   @ s1 = 1   # use base-model stream 1

cat >> $tavg_contents_filename << EOF
$s1  IAGE
EOF
   #-------------------------------------------------------------------------------------
   # Add optional tracer budget terms
   #-------------------------------------------------------------------------------------
if ($OCN_TAVG_TRACER_BUDGET == TRUE) then
cat >> $tavg_contents_filename << EOF
$s1  IAGE_RESET_TEND
$s1  DIA_IMPVF_IAGE
$s1  HDIFE_IAGE
$s1  HDIFN_IAGE
$s1  HDIFB_IAGE
$s1  UE_IAGE
$s1  VN_IAGE
$s1  WT_IAGE
EOF
endif

#  disable the following until they are computed correctly
#  IAGE_SQR 
#  UE_IAGE
#  VN_IAGE
#  WT_IAGE
#  ADV_IAGE
#  J_IAGE
#  Jint_IAGE
#  STF_IAGE
#  RESID_IAGE
#  FvPER_IAGE
#  FvICE_IAGE

else if ($command == prestage) then

   echo ocn.iage.setup.csh : prestaging data files                               >> $LICOM_BLDNML
   echo ------------------------------------------------------------------------ >> $LICOM_BLDNML

else if ($command == document) then

   echo ocn.iage.setup.csh : documenting inputdata files                         >> $LICOM_BLDNML
   echo ------------------------------------------------------------------------ >> $LICOM_BLDNML

else if ($command == ccsm_prestage) then

  #echo ocn.iage.setup.csh : writing ccsm prestaging information                 >> $LICOM_BLDNML

else

   echo ocn.iage.setup.csh : unrecognized command argument $command              >> $LICOM_BLDNML
   echo ------------------------------------------------------------------------ >> $LICOM_BLDNML

   exit 2

endif
