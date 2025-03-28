#!/usr/bin/env perl
#===============================================================================
# SVN $Id: taskmaker.pl 24672 2010-09-10 17:30:01Z fischer $
# SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/scripts/branch_tags/cesm1_0_rel_tags/cesm1_0_3_n02_scripts4_110531b/ccsm_utils/Tools/taskmaker.pl $
#===============================================================================
# This is a program to derive task and thread geometry info
# based on env var values set in a ccsm4 env_pes file
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# parse arg vector: select option
#-------------------------------------------------------------------------------

$sumpeflag    = -1; # total number of pes: (mpi tasks)x(threads)
$sumtasks     = -1; # total number of mpi tasks
$maxthrds     = -1; # max threads over all mpi tasks
$taskgeomflag = -1; # task geometry string for IBM
$thrdgeomflag = -1; # thread geometry string for IBM
$aprunflag    = -1; # aprun options for Cray XT
$psbrsflag    = -1; # pbs resources option for NAS pleiades
$document     = -1; # document the layout
$removedeadtasks = 0; # remove dead tasks or reset to 1

if ($#ARGV  < 0 ){ # no arguments
   $taskgeomflag = 1;  # default option
}
elsif ($#ARGV eq 0 ){ # one argument
#  $opt=shift(@ARGV);
   $opt=$ARGV[0];
   if     ($opt eq "-sumonly"){
      $sumpeflag = 1;
   }
   elsif ($opt eq "-sumpes"){
      $sumpeflag = 1;
   }
   elsif ($opt eq "-sumtasks"){
      $sumtaskflag = 1;
   }
   elsif ($opt eq "-maxthrds"){
      $maxthrdflag = 1;
   }
   elsif ($opt eq "-taskgeom"){
      $taskgeomflag = 1;
   }
   elsif ($opt eq "-thrdgeom"){
      $thrdgeomflag = 1;
   }
   elsif ($opt eq "-aprun"){
      $aprunflag = 1;
   }
   elsif ($opt eq "-pbsrs"){
      $pbsrsflag = 1;
   }
   elsif ($opt eq "-document"){
      $document = 1;
   }
   else {
      print "(taskmaker.pl) Usage: taskmaker.pl [-taskgeom|-thrdgeom|-sumpes|-sumtasks|-maxthrds|-aprun|-pbsrs] \n";
      exit;
   }
}
else { # more than one argument
   print "(taskmaker.pl) Usage: taskmaker.pl [-taskgeom|-thrdgeom|-sumpes|-sumtasks|-maxthrds|-aprun|-pbsrs] \n";
   exit;
}

#-------------------------------------------------------------------------------
# get task/thread layout data via env vars
#-------------------------------------------------------------------------------

$COMP_CPL     = $ENV{'COMP_CPL'};
$NTASKS_CPL   = $ENV{'NTASKS_CPL'};
$NTHRDS_CPL   = $ENV{'NTHRDS_CPL'};
$ROOTPE_CPL   = $ENV{'ROOTPE_CPL'};
$PSTRID_CPL   = $ENV{'PSTRID_CPL'};
$COMP_ATM     = $ENV{'COMP_ATM'};
$NTASKS_ATM   = $ENV{'NTASKS_ATM'};
$NTHRDS_ATM   = $ENV{'NTHRDS_ATM'};
$ROOTPE_ATM   = $ENV{'ROOTPE_ATM'};
$PSTRID_ATM   = $ENV{'PSTRID_ATM'};
$COMP_LND     = $ENV{'COMP_LND'};
$NTASKS_LND   = $ENV{'NTASKS_LND'};
$NTHRDS_LND   = $ENV{'NTHRDS_LND'};
$ROOTPE_LND   = $ENV{'ROOTPE_LND'};
$PSTRID_LND   = $ENV{'PSTRID_LND'};
$COMP_ICE     = $ENV{'COMP_ICE'};
$NTASKS_ICE   = $ENV{'NTASKS_ICE'};
$NTHRDS_ICE   = $ENV{'NTHRDS_ICE'};
$ROOTPE_ICE   = $ENV{'ROOTPE_ICE'};
$PSTRID_ICE   = $ENV{'PSTRID_ICE'};
$COMP_OCN     = $ENV{'COMP_OCN'};
$NTASKS_OCN   = $ENV{'NTASKS_OCN'};
$NTHRDS_OCN   = $ENV{'NTHRDS_OCN'};
$ROOTPE_OCN   = $ENV{'ROOTPE_OCN'};
$PSTRID_OCN   = $ENV{'PSTRID_OCN'};
$COMP_GLC     = $ENV{'COMP_GLC'};
$NTASKS_GLC   = $ENV{'NTASKS_GLC'};
$NTHRDS_GLC   = $ENV{'NTHRDS_GLC'};
$ROOTPE_GLC   = $ENV{'ROOTPE_GLC'};
$PSTRID_GLC   = $ENV{'PSTRID_GLC'};
$MAXTPN       = $ENV{'MAX_TASKS_PER_NODE'};
$MACH         = $ENV{'MACH'};

if ($MACH eq 'pleiades_wes') {
    $NAS_NODES = 'wes';
} elsif ($MACH eq 'pleiades') {
    $NAS_NODES = 'har';
}

#$NAS_NODES = substr($MACH,9);


if ($MAXTPN < 1) {$MAXTPN = 1 ;}

@mcomps = (  $COMP_CPL,   $COMP_ATM,   $COMP_LND,   $COMP_ICE,   $COMP_OCN,   $COMP_GLC);
@ntasks = ($NTASKS_CPL, $NTASKS_ATM, $NTASKS_LND, $NTASKS_ICE, $NTASKS_OCN, $NTASKS_GLC);
@nthrds = ($NTHRDS_CPL, $NTHRDS_ATM, $NTHRDS_LND, $NTHRDS_ICE, $NTHRDS_OCN, $NTHRDS_GLC);
@rootpe = ($ROOTPE_CPL, $ROOTPE_ATM, $ROOTPE_LND, $ROOTPE_ICE, $ROOTPE_OCN, $ROOTPE_GLC);
@pstrid = ($PSTRID_CPL, $PSTRID_ATM, $PSTRID_LND, $PSTRID_ICE, $PSTRID_OCN, $PSTRID_GLC);

#print "ntasks = @ntasks \n";
#print "nthrds = @nthrds \n";
#print "rootpe = @rootpe \n";
#print "pstrid = @pstrid \n";

#-------------------------------------------------------------------------------
# compute total number of mpi tasks
#-------------------------------------------------------------------------------

$tottasks = 0;
for ($c1=0; $c1 <= $#ntasks; $c1++){
    $n = $ntasks[$c1];
    $t = $nthrds[$c1];
    $r = $rootpe[$c1];
    $p = $pstrid[$c1];

    $tt = $r + ($n - 1) * $p + 1;
    if ($tt > $tottasks) {$tottasks = $tt ;}
}

#-------------------------------------------------------------------------------
# compute max threads for each mpi task
#-------------------------------------------------------------------------------

# initialize maxt, max threads for each task
for ($c1=0; $c1 < $tottasks; $c1++){
    $maxt[$c1] = 0;
}

# compute maxt array (max threads for each task)
for ($c1=0; $c1 <= $#ntasks; $c1++){
    $n = $ntasks[$c1];
    $t = $nthrds[$c1];
    $r = $rootpe[$c1];
    $p = $pstrid[$c1];

    $c2 = 0;
    while ($c2 < $n) {
       $s = $r + $c2 * $p;
       if ($t > $maxt[$s]) {$maxt[$s] = $t;}
       $c2 = $c2 + 1;
    }
}

# remove tasks with zero threads if requested
if ($removedeadtasks > 0) {
  $alltasks = $tottasks;
  for ($c1=0; $c1 < $alltasks; $c1++){
    if ($c1 < $tottasks && $maxt[$c1] < 1) {
      for ($c2=$c1; $c2 < $tottasks-1; $c2++){
        $maxt[$c2] = $maxt[$c2+1];
      }
      $maxt[$tottasks] = 0;
      $tottasks = $tottasks - 1;
    }
  }
}

# compute min/max threads over all mpi tasks and sum threads
# also reset maxt values from zero to one after checking for min values
# but before checking for max and summing
$minthrds = $maxt[0];
$maxthrds = $maxt[0];
$sumt[0] = 0;
for ($c1=1; $c1 < $tottasks; $c1++){ 
   if ($maxt[$c1] < $minthrds) {$minthrds = $maxt[$c1] ;}
   if ($maxt[$c1] < 1) {$maxt[$c1] = 1;}
   if ($maxt[$c1] > $maxthrds) {$maxthrds = $maxt[$c1] ;}
   $sumt[$c1] = $sumt[($c1-1)] + $maxt[($c1-1)];
}

#-------------------------------------------------------------------------------
# compute task & thread geometry strings for use on IBM machines
#-------------------------------------------------------------------------------

$fullsum = 0;     # sum of all tasks on all nodes
$sum = $maxt[0];  # sum of all tasks on one node
$taskgeom = "(0";
$thrdgeom = " $maxt[0]";
$taskcnt = 1;
$thrdcnt = $maxt[0];
$aprun = "";
$pbsrs = "";

for ($c1=1; $c1 < $tottasks; $c1++){     # assign each task to a node
    $sum = $sum + $maxt[$c1];
    if ($sum > $MAXTPN) {
        $fullsum = $fullsum + $MAXTPN;
        $sum = $maxt[$c1];
        $taskgeom = $taskgeom.")($c1";   # this is 1st task on a new node
    }
    else {
        $taskgeom = $taskgeom.",$c1";    # append this task to current node
    }
    $thrdgeom = $thrdgeom.":$maxt[$c1]"; # number of threads assigned to this task
#    $t = $maxt[$c1];
    if ($maxt[$c1] != $thrdcnt) {
      $taskpernode = $MAXTPN / $thrdcnt;
      $aprun = $aprun." -n $taskcnt -N $taskpernode -d $thrdcnt ./ccsm.exe :";
      $nodecnt = $taskcnt / $taskpernode ;
      $pbsrs = $pbsrs."${nodecnt}:ncpus=${MAXTPN}:mpiprocs=${taskpernode}:ompthreads=${thrdcnt}:model=${NAS_NODES}+";
      $thrdcnt = $maxt[$c1];
      $taskcnt = 1;
    }
    else {
      $taskcnt = $taskcnt + 1;
    }
}
$fullsum = $fullsum + $sum;
$taskgeom = $taskgeom.")";
$taskpernode = $MAXTPN / $thrdcnt;
$aprun = $aprun." -n $taskcnt -N $taskpernode -d $thrdcnt ./ccsm.exe";

$nodecnt = $taskcnt / $taskpernode ;
$pbsrs = $pbsrs."${nodecnt}:ncpus=${MAXTPN}:mpiprocs=${taskpernode}:ompthreads=${thrdcnt}:model=${NAS_NODES}";

#print "taskgeom = $taskgeom \n";

#-------------------------------------------------------------------------------
# output what was asked for
#-------------------------------------------------------------------------------

#print "test output1 = $fullsum $tottasks $maxthrds \n";
#print "test output2 = $taskgeom \n";
#print "test output3 = $thrdgeom \n";
#print "test output4 = $aprun \n";

if    ($sumpeflag    > 0) {
    print "$fullsum";
}
elsif ($sumtaskflag  > 0) {
    print " $tottasks";
}
elsif ($maxthrdflag  > 0) {
    print "$maxthrds";
}
elsif ($taskgeomflag > 0) {
    print "$taskgeom";
}
elsif ($thrdgeomflag > 0) {
    print "$thrdgeom";
}
elsif ($aprunflag > 0) {
    print "$aprun";
}
elsif ($pbsrsflag > 0) {
    print "$pbsrs";
}
elsif ($document > 0) {
    print "# ---------------------------------------- \n";
    print "# PE LAYOUT: \n";
    print "#   total number of tasks  = $tottasks \n";
    print "#   maximum threads per task = $maxthrds \n";
  for ($c1=0; $c1 <= $#ntasks; $c1++) {
    $n = $ntasks[$c1];
    $t = $nthrds[$c1];
    $r = $rootpe[$c1];
    $p = $pstrid[$c1];
    $tt = $r + ($n - 1) * $p;
    print "#   $mcomps[$c1] ntasks=$n  nthreads=$t rootpe=$r \n";
  }
    print "#   \n";
    print "#   total number of hw pes = $fullsum \n";
  for ($c1=0; $c1 <= $#ntasks; $c1++) {
    $n = $ntasks[$c1];
    $t = $nthrds[$c1];
    $r = $rootpe[$c1];
    $p = $pstrid[$c1];
    $tt = $r + ($n - 1) * $p;
    $tm = $sumt[$tt] + $t - 1;
    print "#     $mcomps[$c1] hw pe range ~ from $sumt[$r] to $tm \n";
  }
  if ($minthrds < 1) {
    print "#   \n";
    print "#   WARNING there appear to be some IDLE hw pes \n";
    print "#   Please consider reviewing your env_mach_pes file \n";
    }
    print "# ---------------------------------------- \n";
}
else {
    print "(taskmaker.pl) internal output selection error";
}

exit;

#===============================================================================
