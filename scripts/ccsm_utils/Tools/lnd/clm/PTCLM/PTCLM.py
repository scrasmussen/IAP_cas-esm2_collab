#!/usr/bin/env python
#########################################################################################
#
# PTCLM.py
#
# Python script to create cases to run point simulations of CLM4
# using Tower Datasets for Ameriflux tower sites, using the CESM1
# framework.
#
# Python script originally created by: 
#
#  Daniel M. Riccciuto, Dali Wang, Peter E. Thornton, Wilfred M. Poist
#  of Environmental Sciences Division, Oak Ridge National Lab.
#
#  Quinn Thomas
#  of Cornell University
#
#  Modified by Erik Kluzek (NCAR) to be incorporated as a standard part of CLM4.
#
#  For help on PTCLM.py type:
#
#   PTCLM.py --help
#
#  Also see the README file
#
#  Requirements:
#
#  python, UNIX shell, NCL (NCAR Command Language)
#  To create tools: GNU make, Fortran compiler, C compiler
#
# NOTE:  mkgriddata mksurfdata,surf mkdatadomain must be compiled!
#           You should only have to compile them once.
#           you must also have ncl installed.
#
#########################################################################################
description = 'Python script to create cases to run single point simulations with tower site data.'
import os, csv, time, re, sys
from   xml.sax.handler import ContentHandler
from   xml.sax         import make_parser 

######  THE ERROR FUNCTION
##############################################################

def error( desc ):
     "error function"
     print "ERROR("+sys.argv[0]+"):: "+desc
     os.abort()

######  SET SOME VARIABLES ##############################################################

#configure case options          
myres="pt1_pt1"           #single-point mode (don't change)
#run time defaults
defmyrun_units="default"  #default time units to run
defmyrun_n=-999           #default number of time to run
defSitesGroup = "EXAMPLE" #default site group name
defmyscratch  = "defscr"  #default scratch root directory

ccsm_input="default"

######  GET VERSION INFORMATION #########################################################

if sys.version_info < (2, 4):
    error( "The version of Python being used is too old for PTCLM" )


svnurl="$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/PTCLM/branch_tags/cesm1_0_rel_tags/cesm1_0_3_n01_PTCLM1_110504/PTCLM.py $"
if   ( svnurl.split('/')[4] == "trunk"       ):
   svnvers="scripts_trunk"
elif ( svnurl.split('/')[4] == "trunk_tags"  ):
   svnvers=svnurl.split('/')[5]
elif ( svnurl.split('/')[4] == "branches"    ):
   svnvers="scripts_branch_"+svnurl.split('/')[5]
elif ( svnurl.split('/')[4] == "branch_tags" ):
   svnvers="scripts_brnchtag_"+svnurl.split('/')[6]
else:
   print( "Error getting version from: "+svnurl)
   os.abort()
version="PTCLM"+str(0.5)+"_"+svnvers

### PARSE THE COMMAND LINE INPUT ########################################################

from optparse import OptionParser, OptionGroup

#parse arguments
cmdline = ""
for arg in sys.argv:
    cmdline = cmdline+arg+" "
parser = OptionParser( usage="%prog [options] -d inputdatadir -m machine -s sitename", description=description, version=version )
required = OptionGroup( parser, "Required Options" )
required.add_option("-d", "--csmdata", dest="ccsm_input", default=" ", \
                  help="Location of CCSM input data")
required.add_option("-m", "--machine", dest="mymachine", default="none" \
                  ,help="Machine, valid CESM script machine (-m list to list valid machines)")
required.add_option("-s", "--site", dest="mysite", default="none", \
                  help="Site-code to run, FLUXNET code or CLM1PT name (-s list to list valid names)")
parser.add_option_group(required)
options  = OptionGroup( parser, "Configure and Run Options" )
options.add_option("-c", "--compset", dest="mycompset", default="ICN", \
                  help="Compset for CCSM simulation (Must be a valid 'I' compset [other than IG compsets], use -c list to list valid compsets)")
options.add_option("--coldstart", dest="coldstart", action="store_true", default=False, \
                  help="Do a coldstart with arbitrary initial conditions")
options.add_option("--ad_spinup", action="store_true", \
                  dest="ad_spinup", default=False, \
                  help="Run accelerated decomposition spinup (CN compset only)")
options.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
options.add_option("--cesm_root", dest="base_cesm", \
                  default=" ", help = \
                  "Root CESM directory (top level directory with models and scripts subdirs)")
options.add_option("--debug", dest="debug", action="store_true", default=False, \
                  help="Flag to turn on debug mode so won't run, but display what would happen")
options.add_option("--exit_spinup", action="store_true", \
                  dest="exit_spinup", default=False, \
                  help="Exit from accelerated decomposition spinup (CN compset only)")
options.add_option("--final_spinup", action="store_true", \
                  dest="final_spinup", default=False, \
                  help="Use for one more spinup for at least 50 years in normal mode")
options.add_option("--finidat", dest="finidat", default=" ", \
                  help="Name of finidat initial conditions file to start CLM from")
options.add_option("--list", dest="list", default=False, action="store_true", \
                  help="List all valid: sites, compsets, and machines")
options.add_option("--namelist", dest="namelist", default=" " \
                  ,help="List of namelist items to add to CLM namelist "+ \
                        "(example: --namelist=\"hist_fincl1='TG',hist_nhtfrq=-1\"" )
options.add_option("--QIAN_tower_yrs",action="store_true",\
                  dest="QIAN_tower_yrs",default=False,\
                  help="Use the QIAN forcing data year that correspond to the tower years")
options.add_option("--rmold", dest="rmold", action="store_true", default=False, \
                  help="Remove the old case directory before starting")
options.add_option("--run_n", dest="myrun_n", default=defmyrun_n, \
                  help="Number of time units to run simulation" )
options.add_option("--run_units", dest="myrun_units", default=defmyrun_units, \
                  help="Time units to run simulation (steps,days,years, etc.)")
options.add_option("--scratchroot", dest="scratchroot", default=defmyscratch, \
                  help="Directory name of scratch space to build/run model in (can only be set if using a generic machine)")
options.add_option("--quiet", action="store_true", \
                  dest="quiet", default=False, \
                  help="Print minimul information on what the script is doing")
options.add_option("--sitegroupname", dest="sitegroup", default=defSitesGroup, \
                  help="Name of the group of sites to search for you selected site in "+ \
                  "(look for prefix group names in the PTCLM_sitedata directory)")
options.add_option("--stdurbpt", dest="stdurbpt", default=False, action="store_true", \
                  help="If you want to setup for standard urban namelist settings")
options.add_option("--useQIAN", dest="useQIAN", help= \
                  "use QIAN input forcing data instead of tower site meterology data", \
                  default=False, action="store_true")
options.add_option("--verbose", action="store_true", \
                  dest="verbose", default=False, \
                  help="Print out extra information on what the script is doing")
parser.add_option_group(options)

suprtclm1ptSettings="For supported CLM1PT single-point datasets, you MUST run with the "+ \
              "following settings:" + \
              " --nopointdata" + \
              " --ndepgrid" + \
              " And you must NOT set any of these:" + \
              " --soilgrid" + \
              " --pftgrid"+ \
              " --aerdepgrid"+ \
              " --owritesrfaer"
indatgengroup = OptionGroup( parser, "Input data generation options", \
                  "These are options having to do with generation of input datasets.  " + \
                  "Note: When running for supported CLM1PT single-point datasets you can NOT generate new datasets.  "+ \
                  suprtclm1ptSettings )
parser.add_option_group(indatgengroup)
indatgengroup.add_option("--aerdepgrid", action="store_true", \
                  dest="aerdepgrid", help="Do NOT regrid aerosol deposition data (use standard global grid data)", \
                  default=False)
indatgengroup.add_option("--ndepgrid", action="store_true", \
                  dest="ndepgrid", help="Do NOT regrid Nitrogen deposition data (use standard global grid data)", \
                  default=False)
indatgengroup.add_option("--nopointdata", action="store_true", \
                  dest="nopointdata", help="Do NOT make point data (use data already created)", \
                  default=False)
indatgengroup.add_option("--owritesrfaer", action="store_true", \
                  dest="owritesrfaer", help=\
                  "Overwrite the existing surface/aerosol datasets if they exist (normally do NOT recreate them)", \
                  default=False)
indatgengroup.add_option("--pftgrid", dest="pftgrid", help = \
                  "Use pft information from global gridded file (rather than site data)", \
                  action="store_true", default=False)
indatgengroup.add_option("--soilgrid", dest="soilgrid", help = \
                  "Use soil information from global gridded file (rather than site data)",\
                   action="store_true", default=False)
versiongroup  = OptionGroup( parser, "Main Script Version Id: $Id: PTCLM.py 28806 2011-06-04 05:08:01Z erik $ Scripts URL: "+svnurl )
parser.add_option_group(versiongroup)

(options, args) = parser.parse_args()
if len(args) != 0:
    parser.error("incorrect number of arguments")

### END PARSE THE COMMAND LINE INPUT ####################################################

### SOME FUNCTIONS    ###################################################################

def system( cmd ):
     "system function with error checking and debug prining"
     if plev>0: print "Run command: "+cmd
     if ( not options.debug or cmd.startswith( "./create_newcase" ) ):
        rcode = os.system( cmd )
     else: rcode = 0
     if ( rcode != 0 ):
        error( "Error running command"+cmd )

def queryFilename( queryopts, filetype ):
    "query the XML database to get a filename"
    query = abs_base_cesm+"/models/lnd/clm/bld/queryDefaultNamelist.pl -silent " \
             +"-justvalue "
    if ( ccsm_input != "default" ): 
       query = query + " -csmdata "+ccsm_input
    cmd = query+queryopts+" -var "+filetype
    file = os.popen( cmd )
    filename = file.read() 
    if ( file.close() != None ):
       print "Query = "+cmd
       error( "Error getting file from XML database" )
    # Remove the trailing new line from the filename
    if ( (filename == None) or (filename == "") ): 
       print "Query = "+cmd
       error( "Trouble finding file from XML database: "+filetype )
    return( filename.replace( "\n", "" ) )

def Get_envconf_Value( var ):
     'Function to get the value of a variable from the env_conf.xml file'
     pwd = os.getcwd()
     os.chdir( mycase)
     env_val = re.search('value="([^"]*)"', os.popen("grep "+var+" env_conf.xml").read() ).group(1)
     if ( env_val == None ): error( "Did NOT find envconf value for: "+var )
     os.chdir(pwd)
     return env_val

def xmlchange_envconf_value( var, value ):
     'Function to set the value of a variable in the env_conf.xml file'
     cmd = "./xmlchange -file env_conf.xml -id "+var+" -val "+value
     system( cmd )

def xmlchange_envrun_value( var, value ):
     'Function to set the value of a variable in the env_run.xml file'
     cmd = "./xmlchange -file env_run.xml -id "+var+" -val "+value
     system( cmd )

if sys.version_info < (2, 5):
   def rpartition( string, sep ):
       'Reverse order of dividing string by seperator'
       before = string[0:string.rfind(sep)];
       after  = string[before.count(""):];
       return ( before, sep, after )
#
# Some classes
#

#
# List the machines in the config_machines.xml file
#
class MachineList( ContentHandler ):

   def startDocument(self):
     self.list = [];
   
   def startElement(self, name, attrs):
     if name == 'machine':     
       self.list.append( str( attrs.get('MACH',"") ) )

   def endDocument(self):
     print "\nValid Machines: "+str(self.list)+"\n\n";
#
# List the compsets  in the config_compsets.xml file
#
class ICompSetsList( ContentHandler ):

   def startDocument(self):
     self.list = [];
   
   def startElement(self, name, attrs):
     if name == 'compset':
       name = str( attrs.get('NAME',"") )
       if ( name.startswith( "I" ) and not name.endswith( "_GLC" ) ): 
          nlen = len(self.list)
          if (   nlen == 0                 ): self.list.append( name )
          elif ( name != self.list[nlen-1] ): self.list.append( name )

   def endDocument(self):
     print "\nValid Compsets: "+str(self.list)+"\n\n";


###### SET OPTIONS BASED ON INPUT FROM PARSER  ##########################################

mymachine  = options.mymachine
mysite     = options.mysite
mycompset  = options.mycompset
SitesGroup = options.sitegroup
infohelp   = "\n\n Use --help option for help on usage.\n";
if(options.list):
    mycompset = "list"
    mymachine = "list"
    mysite    = "list"
if ( mysite == "none" ): parser.error("sitename is a required argument, set it to a valid value"+infohelp )
if ( (options.ad_spinup and options.exit_spinup) or (options.exit_spinup and options.final_spinup) or \
     (options.ad_spinup and options.final_spinup) ):
    parser.error( "More than one spinup option is selected"+infohelp )
if ( options.stdurbpt and (options.ad_spinup or options.exit_spinup or \
     options.final_spinup) ):
    parser.error( "Standard urban point namelist setup option is incompatible with ANY spinup option"+infohelp )
if ( options.verbose and options.quiet ):
    parser.error( "options quiet and verbose are mutually exclusive"+infohelp )
if(options.ad_spinup and options.finidat != " "):
    parser.error( "ad_spinup and setting the finidat file are mutually exclusive"+infohelp )
if(options.coldstart and options.finidat != " "):
    parser.error( "coldstart and setting the finidat file are mutually exclusive"+infohelp )

if (   options.verbose ): plev = 2
elif ( options.quiet   ): plev = 0
else:                     plev = 1

sitedata=SitesGroup+"_sitedata.txt"
soildata=SitesGroup+"_soildata.txt"
pftdata=SitesGroup+"_pftdata.txt"


if plev>0: print "---------------- PTCLM version "+str(version)+"-----------------------------\n"
if plev>0: print "   "+cmdline+"\n"
if plev>0: print "   OPTIONS:\n" 
if plev>0: print "Site name:\t\t\t\t\t\t"+mysite+"\n"

#Set case name based on site and other information

if ( options.mycaseid == "" ):
   mycasename=mysite
else:
   mycasename=options.mycaseid+"_"+mysite

mycasename=mycasename+"_"+mycompset
  
if options.useQIAN:
   mycasename+="_QIAN"
  
if options.ad_spinup:
   mycasename+="_ad_spinup"
elif options.exit_spinup:
   mycasename+="_exit_spinup"
elif options.final_spinup:
   mycasename+="_final_spinup"

if plev>0: print "CESM Component set:\t\t\t\t\t"+mycompset
if plev>0: print "CESM machine:\t\t\t\t\t\t"+options.mymachine

base_cesm = options.base_cesm
if base_cesm == " ":
    #assume base directory is one level up from where script
    #  is executed, if not specified
    stdout    = os.popen("cd ../../../../../..; pwd")
    base_cesm = os.path.abspath( str.strip(stdout.read()) )
    ptclm_dir = base_cesm+"/scripts/ccsm_utils/Tools/lnd/clm/PTCLM"
else:
    stdout    = os.popen("pwd")
    ptclm_dir = str.strip(stdout.read())

abs_base_cesm = os.path.abspath( base_cesm )
if plev>0: print "Root CLM directory:\t\t\t\t\t"+abs_base_cesm

# Get the case directory name
if mycasename.startswith("/"): mycase = mycasename
else:                          mycase = abs_base_cesm+"/scripts/"+mycasename
if sys.version_info < (2, 5):
   mycasedir  = rpartition(mycase,"/")[0]
   mycasename = rpartition(mycase,"/")[2]
else:
   mycasedir  = mycase.rpartition("/")[0]
   mycasename = mycase.rpartition("/")[2]
if plev>0: print "Case name:\t\t\t\t\t\t"+mycasename
if plev>0: print "Case directory:\t\t\t\t\t"+mycasedir

ad_spinup  = options.ad_spinup
if plev>0: print "Accelerated Decomposition mode:\t\t\t\t"+str(ad_spinup)

exit_spinup = options.exit_spinup
if plev>0: print "Exit spinup mode:\t\t\t\t\t"+str(exit_spinup)

final_spinup=options.final_spinup
if plev>0: print "Final spin up:\t\t\t\t\t\t"+str(final_spinup)

useQIAN=options.useQIAN
if plev>0: print "Using QIAN climate inputs\t\t\t\t"+str(useQIAN)

finidat    = options.finidat
if (finidat == " "):
    if plev>0: print "Finidat file:\t\t\t\t\t\t<none>"
else:
    if plev>0: print "Finidat file:\t\t\t\t\t\t"+finidat
if plev>0:     print "Sites group name:\t\t\t\t\t"+SitesGroup

if plev>0: print "Use preexisting point data:\t\t\t\t"+str(options.nopointdata)
if (options.nopointdata):
    makeptfiles=False
else:
    makeptfiles=True
    if plev>0: print "** Surface data file will be built using site-level data " + \
          "when available unless otherwise specified ** \n"
    if plev>0: print "\tExtract PFT data from gridded files:\t\t"+str(options.pftgrid)
    if plev>0: print "\tExtract soil data from gridded files:\t\t"+str(options.soilgrid)

if(mymachine == "list"):
    machXML = make_parser()
    machXML.setContentHandler(MachineList()) 
    machXML.parse( abs_base_cesm+"/scripts/ccsm_utils/Machines/config_machines.xml" )

if(mycompset == "list"):
    compXML = make_parser()
    compXML.setContentHandler(ICompSetsList()) 
    compXML.parse( abs_base_cesm+"/scripts/ccsm_utils/Case.template/config_compsets.xml" )


###### END SET OPTIONS BASED ON INPUT FROM PARSER  ######################################

########## GET SITE LAT, LON, AND TOWER MET YEARS #######################################

siteDir = "PTCLM_sitedata"
os.chdir(ptclm_dir+"/"+siteDir)
#get lat/lon, start/end years from sitedata file
suprtclm1pt = True
if plev>0: print "\nOpen Site data file: "+siteDir+"/"+sitedata+"\n"
AFdatareader = csv.reader(open(sitedata, "rb"))
if ( mysite == "list" ):  plev = 2
for row in AFdatareader:
    if plev>1: print " site = %9s name: %-55s Region: %s" % ( row[0], row[1], row[2] )
    if row[0] == mysite:
        suprtclm1pt = False
        lon=float(row[3])
        if (lon < 0):
            lon=360.0+float(row[3]) 
        lat=float(row[4])
        startyear=int(row[6])
        endyear=int(row[7])
        alignyear = int(row[8])

if ( mysite == "list" ): 
  print "\nSupported CLM1PT name dataset names are:\n";
  resList = queryFilename( " -res list", "none" )
  resols = resList.split(" ");
  for i, res in enumerate(resols):
     if ( res.startswith( "1x1_" ) ):
        print " site = %9s " % ( res );

# Exit early for list options
if ( mysite == "list" or mycompset == "list" or mymachine == "list" ): 
  exit()

# inputdata directory -- set after list options
ccsm_input=options.ccsm_input
if ccsm_input == " ":
   parser.error( "inputdatadir is a required argument, set it to the directory where you have your inputdata"+infohelp )
if plev>0: print "CCSM input data directory:\t\t\t\t"+ccsm_input
#define data and utility directroies
clm_tools  = abs_base_cesm+'/models/lnd/clm/tools'
clm_input  = ccsm_input+'/lnd/clm2'
datm_input = ccsm_input+'/atm/datm7'

mask = "navy"
if ( suprtclm1pt ):
   if plev>0: print "Did NOT find input sitename:"+mysite+" in sitedata:"+sitedata
   if plev>0: print "Assuming that this is a supported CLM1PT single-point dataset"
   if ( not options.nopointdata or options.soilgrid or options.pftgrid or \
        not options.ndepgrid or options.aerdepgrid or options.owritesrfaer ):
      error( suprtclm1ptSettings )
   clmusrdatname    = ""
   clmres           = mysite
   pft_phys_out     = ""
   clmusrdat        = ""
else:
   clmusrdatname    = "1x1pt_"+mysite
   clmusrdat        = " -usrname "+clmusrdatname
   clmres           = clmusrdatname
   pft_phys_file    = queryFilename( " ", "fpftcon" )
   pft_phys_out     = re.search( "(.+)\.nc$", pft_phys_file ).group(1)+"."+mysite+".nc"

if plev>0: print "----------------------------------------------------------------\n"

######   CREATE NEW CASE and GET USE_CASE ###############################################

os.chdir(ptclm_dir)      
if plev>0: print "Creating new case\n"

os.chdir(abs_base_cesm+"/scripts")
if ( mymachine.find( "generic" ) == 0 ):
    if ( options.scratchroot == defmyscratch ): options.scratchroot = abs_base_cesm+"/run"
    opt = " -scratchroot "+options.scratchroot+" -max_tasks_per_node 1 " \
         +" -din_loc_root_csmdata "+ccsm_input
else:
    opt = " "
    if ( mymachine == "none" ): parser.error( "machine is a required argument, set it to a valid value"+infohelp )
    if ( options.scratchroot != defmyscratch ): parser.error( "scratchroot can only be set for a generic machine"+infohelp )

if ( options.rmold ): system( "/bin/rm -rf "+mycase )

cmd = "./create_newcase -case "+mycase+" -mach "+mymachine+" -compset "+mycompset \
      +" -res "+myres+opt
system( cmd )

clmnmlusecase    = Get_envconf_Value( "CLM_NML_USE_CASE" )

# Get any options already set in CLM_CONFIG_OPTS and check for consistency ##############
clmconfigopts    = Get_envconf_Value( "CLM_CONFIG_OPTS" )
datmpresaero     = Get_envconf_Value( "DATM_PRESAERO" )

bgctypeCN = re.search('-bgc (cn[a-z]*)', clmconfigopts )

if ( bgctypeCN == None ):
   if (ad_spinup):  parser.error( "ad_spinup only works with CN compsets"+infohelp   )
   if (exit_spinup): parser.error( "exit_spinup only works with CN compsets"+infohelp )

filen = mycase+"/README.PTCLM"
if plev>0: print "Write "+filen+" with command line"
output = open( filen,'w')
output.write(cmdline+"\n")
output.close

############# GET SIM_YEAR, RCP and SIM_YEAR_RANGE based on USE-CASE ####################
############# CLM configure ensures naming conventions are followed  ####################
############# And setup Query options based on them #####################################

if (   clmnmlusecase.endswith("_transient") ):
     transient = re.search('^([0-9]+-[0-9]+)_*(.*)_(transient$)',   clmnmlusecase )
     if ( transient ):
        sim_year_range = transient.group(1)
        sim_year       = re.search( '^([0-9]+)-',    transient.group(1) ).group(1)
        rcpcase        = re.search( '^rcp([0-9.]+)', transient.group(2) )
        if ( rcpcase == None ): rcp = -999.9
        else:                   rcp = rcpcase.group(1)
     elif ( clmnmlusecase.startswith("20thC_") ):
        sim_year_range = "1850-2000"
        sim_year       = "1850"
        rcp            = "-999.9"
     else:
        error( "Can not parse use-case name, does not follow conventions:"+clmnmlusecase )

     if ( sim_year_range == "1850-2000" ): actual_year_range = "1849-2006"
     else:                                 actual_year_range = sim_year_range
elif ( clmnmlusecase.endswith("_control") ):
          sim_year       = re.search( '^([0-9]+)_', clmnmlusecase ).group(1)
          if ( sim_year == None ): error( "Trouble finding sim_year from:"+clmnmlusecase )
          sim_year       = str(sim_year)
          sim_year_range = "constant"
          rcp            = str(-999.9)
elif ( clmnmlusecase.endswith("_pd") or clmnmlusecase == "UNSET" ):
          sim_year       = "2000"   
          sim_year_range = "constant"
          rcp            = str(-999.9)
else:
          error( "Can not parse use-case name:, does not follow conventions"+clmnmlusecase )

qoptionsbase   = " -options mask="+mask+",rcp="+rcp+",datm_presaero="+datmpresaero
   
if ( ad_spinup    ): 
   qoptionsbase += ",bgc=cn,ad_spinup=on"
if ( exit_spinup  ): 
   qoptionsbase += ",bgc=cn,exit_spinup=on"
if ( final_spinup ): 
   qoptionsbase += ",final_spinup=on"

qoptions       = qoptionsbase+",sim_year="+sim_year+",sim_year_range="+sim_year_range;
queryOpts      = " -onlyfiles -res "+clmres+clmusrdat+qoptions
queryOptsNousr = qoptions
queryOptsNavy  = " -res 0.33x0.33 "+qoptions

if ( suprtclm1pt ):
    supqryOpts = queryOptsNousr+" -namelist default_settings"
    startyear  = queryFilename( supqryOpts, "datm_cycle_beg_year" )
    endyear    = queryFilename( supqryOpts, "datm_cycle_end_year" )
    alignyear  = startyear

####### ANY OTHER LAST SETTINGS BEFORE CREATING DATASETS ################################

myrun_n     = options.myrun_n
myrun_units = options.myrun_units

#default simulation length for different types of runs
qryOpts  = queryOptsNousr + " -namelist seq_timemgr_inparm"
if (  myrun_units == defmyrun_units ):
   myrun_units = queryFilename( qryOpts, "stop_option" )
if (  myrun_n == defmyrun_n ):
   myrun_n     = queryFilename( qryOpts, "stop_n"      )

if plev>0: print "Number of simulation "+myrun_units+" to run:\t\t\t\t"+str(myrun_n)
############# BEGIN CREATE POINT DATASETS ###############################################


if makeptfiles:
    if plev>0: print("Making input files for the point (this may take a while if creating transient datasets)")

    gridfile   = queryFilename( queryOpts, "fatmgrid"   )
    fracfile   = queryFilename( queryOpts, "fatmlndfrc" )
    surffile   = queryFilename( queryOpts, "fsurdat"    )
    domainfile = queryFilename( queryOpts+" -namelist shr_strdata_nml", "domainfile" )

    mksrf_fnavyoro  = queryFilename( queryOptsNavy+" -namelist clmexp", "mksrf_fnavyoro" )

    os.chdir(ptclm_dir)
    #make and move grid and frac files ##################################################
    if plev>0: print "Creating grid data"
    #write mkgriddata namelist with site lat/lon and correct input
    input  = open("usr_files/mkgriddata.TEMPLATE")
    output = open(clm_tools+"/mkgriddata/namelist",'w')
    line=0
    for s in input:
        line +=1
        if line == 1:
            output.write(s)
        if line == 2:
            output.write(" mksrf_fnavyoro = '"+ mksrf_fnavyoro+"'\n")
        if line == 3 or line == 4:
            output.write(s)
        if line == 5:
            output.write(s.replace("SITEE",str(lon+0.05)))
        if line == 6:
            output.write(s.replace("SITEW",str(lon-0.05)))
        if line == 7:
            output.write(s.replace("SITES",str(lat-0.05)))
        if line == 8:
            output.write(s.replace("SITEN",str(lat+0.05)))
        if line > 8:
            output.write(s)
    output.close()
    input.close()
    if plev>1: os.system( "cat "+clm_tools+"/mkgriddata/namelist" )
    system(clm_tools+"/mkgriddata/mkgriddata < "+clm_tools+"/mkgriddata/namelist > mkgriddata.log")
    system("/bin/mv -f ./fracdata_0001x0001.nc "+fracfile )
    system("/bin/mv -f ./griddata_0001x0001.nc "+gridfile )

    if (sim_year_range != "constant"):
       pftdynfile = queryFilename( queryOpts, "fpftdyn" )
    else:
       pftdynfile = None
    #make surface data and dynpft #######################################################
    if ( (not options.owritesrfaer) and os.path.exists( surffile) and \
         ((pftdynfile == None) or os.path.exists( pftdynfile ) ) ):
        print "\n\nWARNING: Use existing surface file rather than re-creating it:\t"+surffile
    else:
        if plev>0: print "\n\nRe-create surface dataset:\t"+surffile
        if ( os.path.exists( surffile ) ): print "Over write file: "+surffile
        if ( sim_year_range == "constant" ):
           mksrfyears = sim_year
        else:
           mksrfyears = sim_year_range

        # --- use site-level data for mksurfdata when available ----
        #PFT information for the site
        if (options.pftgrid == False):
            if plev>0: print "Replacing PFT information in surface data file"
            os.chdir(ptclm_dir+"/PTCLM_sitedata")
            AFdatareader = csv.reader(open(pftdata, "rb"))
            os.chdir(ptclm_dir)
            pft_frac=[0,0,0,0,0]
            pft_code=[0,0,0,0,0]
            found=0
            for row in AFdatareader:
                if plev>1: print " site = %9s" % row[0]
                if row[0] == mysite:
                    found=1
                    output=open("./tempsitePFT.txt","w")      
                    output.write(' '.join(row[1:11]))
                    output.close()
                    for thispft in range(0,5):
                        pft_frac[thispft]=float(row[1+2*thispft])
                        pft_code[thispft]=int(row[2+2*thispft])
            if ( found == 0 ):
               error( "Did NOT find input sitename:"+mysite+" in pftdata:"+pftdata+ \
                      " run with pftgrid instead")
            # Find index of first zero
            for i in range(0,len(pft_frac)):
               if ( pft_frac[i] == 0.0 ):
                  nzero = i
                  break
            pftopts=" -pft_frc \""+str(pft_frac[0:nzero])+'"' \
                       " -pft_idx \""+str(pft_code[0:nzero])+'"'
        else: 
            pftopts=""
   
        #Read in the soil conditions for the site #######################################
        if (options.soilgrid == False):

            #soil information
            os.chdir(ptclm_dir+"/PTCLM_sitedata")
            if plev>0: print "Replacing soil information in surface data file"
            AFdatareader = csv.reader(open(soildata, "rb"))
            os.chdir(ptclm_dir)
            found=0
            for row in AFdatareader:
                if plev>1: print " site = %9s" % row[0]
                if row[0] == mysite:
                    found=1
                    output=open("./tempsitesoil.txt","w")
                    output.write(' '.join(row[1:7]))
                    output.close()
                    # The first three items are NOT used
                    soil_depth = float(row[1])  # This is ignored
                    n_layers   = int(row[2])    # This is ignored
                    layer_depth = float(row[3]) # This is ignored
                    sandpct     = float(row[4])
                    claypct     = float(row[5])
            if ( found == 0 ):
               error( "Did NOT find input sitename:"+mysite+" in soildata:"+soildata+ \
                      " run with soilgrid instead")
            if plev>0: print " sandpct="+str(sandpct)+" claypct="+str(claypct)
            soilopts=" -soil_cly "+str(claypct)+" -soil_snd "+str(sandpct)
        else: soilopts=""
        #----- create dynamic pft input file --------------- ############################
        if (options.pftgrid == False) and (sim_year_range != "constant"):

            if plev>0: print "Creating site-specific dynamics PFTs and harvesting"

            pftdyn_site_filename = ptclm_dir + "/PTCLM_sitedata/" + \
                                   mysite + "_dynpftdata.txt"

            # only set dynpft file if the file exists
            if ( os.path.exists( pftdyn_site_filename ) ):
               if plev>0: print "Transition PFT file exists, so using it for changes in PFT"
               # Convert the file from transition years format to mksurfdata pftdyn format
               cnv = ptclm_dir + \
                     "/PTCLM_sitedata/cnvrt_trnsyrs2_pftdyntxtfile.pl "+sim_year_range
               pftdynoutfile = mycase+"/pftdyn_"+mycasename+".txt"
               system( cnv+pftdyn_site_filename+" > "+pftdynoutfile )
               dynpftopts = " -dynpft "+pftdynoutfile
            else:
               if plev>0: print "Transition PFT file did NOT exist, so proceeding with constant PFT"
               dynpftopts = ""
               
        else: 
            dynpftopts = ""

        # Now run mksurfdata  ###########################################################
        mksurfopts = "-nomv -res "+clmres+clmusrdat+ \
                     " -dinlc "+ccsm_input+" -y "+mksrfyears+" -rcp "+rcp+\
                     soilopts+pftopts+dynpftopts
        system(clm_tools+"/mksurfdata/mksurfdata.pl "+mksurfopts+" > mksurfdata.log")

        #move surface data and pftdyn file to correct location
        system("/bin/mv -f ./surfdata_"+clmres+"_*_c*.nc  "+surffile )
        system("/bin/mv -f ./surfdata_"+clmres+"_*_c*.log "+clm_input+"/surfdata" )
        if (sim_year_range != "constant"):
            system("/bin/mv -f surfdata.pftdyn_"+clmres+"_*.nc "+pftdynfile )

    #make data domain file ##############################################################
    if plev>0: print "Creating data domain"
    #write mkdatadomain namelist with correct inputs
    input  = open("usr_files/mkdatadomain.TEMPLATE")
    output = open(clm_tools+"/mkdatadomain/namelist",'w')
    line=0
    for s in input:
        line +=1
        if line < 3:
            output.write(s)
        if line == 3:
            output.write(" f_fracdata = '"+fracfile+"'\n")
        if line == 4:
            output.write(" f_griddata = '"+gridfile+"'\n")
        if line == 5:
            output.write(" f_domain = '"+domainfile+"'\n")
        if line > 5:
            output.write(s)
    output.close()
    input.close()

    if plev>1: os.system( "cat "+clm_tools+"/mkdatadomain/namelist" )
    system(clm_tools+"/mkdatadomain/mkdatadomain < "+clm_tools+"/mkdatadomain/namelist > mkdatadomain.log")

    #
    #make and move ndep and aerosol deposition files ####################################
    #
    #Always make them for a simulation year range of 1850-2000, and 
    #then later in namelists you will narrow down which year (or years) to focus on
    os.environ["RES"]     = clmres
    os.environ["GRDFIL"]  = gridfile
    os.environ["CSMDATA"] = ccsm_input
    os.environ["SIM_YR"]  = "1850"
    os.environ["RCP"]     = rcp
    if ( sim_year_range == "constant" ):
        os.environ["SIM_YR_RNG"]="1850-2000"
    else:
        os.environ["SIM_YR_RNG"]=sim_year_range


    queryOptsDep = " -onlyfiles -res "+clmres+clmusrdat+qoptionsbase
    queryOptsDep += ",sim_year="+os.environ["SIM_YR"]+",sim_year_range="+os.environ["SIM_YR_RNG"]
    if ( not options.aerdepgrid ):

       os.chdir(clm_tools+"/ncl_scripts")
       aerfile  = queryFilename( queryOptsDep+" -namelist datm_internal", "datm_file_aero" )
       if plev>0: print "Creating aerosol deposition data"
       if plev>0: print "RES="+os.environ["RES"]+" GRDFIL="+os.environ["GRDFIL"]+ \
             " CSMDATA="+os.environ["CSMDATA"]+" SIM_YR="+os.environ["SIM_YR"]+ \
             " SIM_YR_RNG="+os.environ["SIM_YR_RNG"]+" RCP="+os.environ["RCP"]
       if ( os.path.exists( aerfile ) and (not options.owritesrfaer) ):
          print "WARNING: do NOT overwrite existing file: "+aerfile
       else:
          lfile = "aerosoldep_*_mean_"+clmres+"_c*.nc"
          system("/bin/rm -f "+lfile )
          system("ncl "+clm_tools+"/ncl_scripts/aerdepregrid.ncl > aerdepregrid.log")
          system("/bin/mv -f "+lfile+"  "+aerfile )

    if ( not options.ndepgrid and (bgctypeCN != None) ):

       os.chdir(clm_tools+"/ncl_scripts")
       queryOptsNDep = queryOptsDep+",bgc="+bgctypeCN.group(1)+" -namelist ndepdyn_nml"
       if plev>0: print "Creating Nitrogen deposition data"
       ndepfile  = queryFilename( queryOptsNDep, "stream_fldfilename_ndep" )
       if plev>0: print "RES="+os.environ["RES"]+" GRDFIL="+os.environ["GRDFIL"]+ \
             " CSMDATA="+os.environ["CSMDATA"]+" SIM_YR="+os.environ["SIM_YR"]+ \
             " SIM_YR_RNG="+os.environ["SIM_YR_RNG"]+" RCP="+os.environ["RCP"]
       lfile = "fndep_clm_*simyr*_"+clmres+"_c*.nc"
       system("/bin/rm -f "+lfile )
       system("ncl "+clm_tools+"/ncl_scripts/ndepregrid.ncl   > ndepregrid.log")
       system("/bin/mv -f "+lfile+"  "+ndepfile )
  
    # Default PFT-physiology file used to make site-level file ##########################
    pft_phys_file  = queryFilename( queryOptsNousr, "fpftcon" )

    #create site-specific pft-physiology file (only if does NOT exist)
    if ( not os.path.exists( pft_phys_out ) ):
       system("/bin/cp "+pft_phys_file+" "+pft_phys_out )
    if plev>0: print "pft_phys_file = "+pft_phys_out+"\n"

else:
    print "WARNING: nopointdata option was selected.  Model will crash if the site level data have not been created\n"    
   
####### END CREATE POINT DATASETS #######################################################


#####  ENV XML CHANGES ##################################################################
os.chdir(mycase)

xmlchange_envconf_value( "CLM_PT1_NAME",    clmres )
if ( clmusrdatname != "" ):
   xmlchange_envconf_value( "CLM_USRDAT_NAME", clmusrdatname )

if(useQIAN):
    xmlchange_envconf_value(        "DATM_MODE",             "CLM_QIAN" )
    if(options.QIAN_tower_yrs):
       xmlchange_envconf_value(     "DATM_CLMNCEP_YR_START", str(startyear) )
       if(endyear < 2005):
           xmlchange_envconf_value( "DATM_CLMNCEP_YR_END",   str(endyear) )
       else:
           xmlchange_envconf_value( "DATM_CLMNCEP_YR_END",   "2004" )
else:
    xmlchange_envconf_value( "DATM_MODE",             "CLM1PT" )
    xmlchange_envconf_value( "DATM_CLMNCEP_YR_START", str(startyear) )
    xmlchange_envconf_value( "DATM_CLMNCEP_YR_END",   str(endyear) )

xmlchange_envconf_value( "CLM_BLDNML_OPTS", "'-mask "+mask+"'" )

# If MPISERIAL support is available use it
if ( Get_envconf_Value( "MPISERIAL_SUPPORT" ) == "TRUE" ):
   xmlchange_envconf_value( "USE_MPISERIAL", "TRUE"   )

###### SET Spinup and ENV_RUN.XML VALUES ################################################
   
hist_nhtfrq = 0
if (ad_spinup):
   hist_mfilt  = 100
   hist_nhtfrq = -8760
   xmlchange_envconf_value( "CLM_CONFIG_OPTS",     "'"+clmconfigopts+" -ad_spinup on'" )
   xmlchange_envconf_value( "CLM_FORCE_COLDSTART", "'on'"  )
   xmlchange_envrun_value( "STOP_DATE",   "06010101"       )
elif (exit_spinup):
   hist_mfilt  = 12
   xmlchange_envconf_value( "CLM_CONFIG_OPTS",     "'"+clmconfigopts+" -exit_spinup on'" )
   xmlchange_envconf_value( "RUN_TYPE",            "branch" )
   xmlchange_envconf_value( "RUN_REFCASE",         \
                            mycase.replace("_exit_spinup","_ad_spinup" ) )
   xmlchange_envconf_value( "RUN_REFDATE",         "0601-01-01" )
   xmlchange_envconf_value( "GET_REFCASE",         "FALSE" )
elif (final_spinup):
   hist_mfilt  = 12*int(myrun_n)
elif(options.stdurbpt):
   hist_mfilt  = str(myrun_n)+", "+str(myrun_n)+", "+str(myrun_n)
   hist_nhtfrq = "-1,-1,-1"
   if ( clmnmlusecase != "UNSET" and clmnmlusecase != "2000_control" ):
      error( "Option stdurbpt is incompatible with this compset" )
   xmlchange_envconf_value( "CLM_NML_USE_CASE", "stdurbpt_pd" )
   xmlchange_envconf_value( "ATM_NCPL",  str(24) )
else:
   hist_mfilt  = 1200

if ( suprtclm1pt ):
   clmconfigopts = Get_envconf_Value( "CLM_CONFIG_OPTS" )
   xmlchange_envconf_value( "CLM_CONFIG_OPTS", "'"+clmconfigopts+" -sitespf_pt "+ \
                             clmres+"'")
   run_startdate = queryFilename( queryOptsNousr+" -namelist default_settings", "run_startdate" )
   xmlchange_envconf_value( "RUN_STARTDATE",   run_startdate  )
   starttod   = queryFilename( queryOptsNousr+" -namelist seq_timemgr_inparm", "start_tod" )
   xmlchange_envrun_value( "START_TOD", starttod  )
   xmlchange_envconf_value( "DATM_PRESAERO", "pt1_pt1" )

xmlchange_envrun_value( "STOP_N",      str(myrun_n) )
xmlchange_envrun_value( "STOP_OPTION", myrun_units )
xmlchange_envrun_value( "REST_OPTION", str(myrun_units) )
rest_n = max( 1, int(myrun_n) // 5 )
xmlchange_envrun_value( "REST_N",      str(rest_n)     )

xmlchange_envrun_value( "DIN_LOC_ROOT_CSMDATA", ccsm_input   )

if ( options.coldstart ):
   xmlchange_envconf_value( "CLM_FORCE_COLDSTART", "on" )
   if ( Get_envconf_Value( "RUN_TYPE" ) == "hybrid" ): 
      xmlchange_envconf_value( "RUN_TYPE", "startup" )

####  SET NAMELIST OPTIONS ##############################################################

output = open("user_nl_clm",'w')
output.write("&clm_inparm\n")
output.write(   " hist_nhtfrq = "+str(hist_nhtfrq)+"\n" )
output.write(   " hist_mfilt  = "+str(hist_mfilt)+"\n" )
if( pft_phys_out != "" ):
   output.write(" fpftcon = '"+pft_phys_out+"'\n" )
if(options.namelist != " "):
   output.write(options.namelist+"\n")
if(finidat != " "):
   output.write(" finidat = '"+finidat+     "'\n")
output.write("/\n")
output.close()
if plev>1: os.system( "cat user_nl_clm" )

###### END SET Spinup and ENV_RUN.XML VALUES ############################################

if plev>0: print "Scripts created successfully\n"
if plev>0: print "cd "+mycase+" and then..."
if plev>0: print "Configure, build and run your case as normal\n"

###   END PTCLM SCRIPT ####################################################################

