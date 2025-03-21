c----CAMx v5.40 111010
c
c     CAMx.PRM defines model parameters
c
c     Copyright 1996 - 2010
c     ENVIRON International Corporation
c
c     Modifications:
c        7/5/02      Added parameters for IRON-PIG
c        9/3/02      removed MXVECCL
c        1/10/03     Added MXVECDP
c        9/12/03     Revised parameters for IRON-PiG
c        7/20/05     Moved sampling grid parameters from TRACER.COM
c       10/20/05     Standardized the temperature/pressure ranges for
c                    chemistry rates
c        7/1/08      Reconfigured for CAMx/MPI dynamic memory allocation
c
c-----------------------------------------------------------------------
c      MXCELLS -- maximum number of cells (used for local arrays)
c      MXLAYER -- maximum number of layers (used for local arrays)
c      MXSPEC  -- maximum number of species (used for local arrays)
c      MXREACT -- maximum number of reactions (used for local arrays)
c      MXGRID  -- maximum number of grids (used for local arrays)
c      MXPTSRC -- maximum number of point sources (used for local arrays)
c      MXCHDRN -- maximum number of children grids for any parent grid
c      NOVERLP -- number of overlap rows/columns for each slice
c-----------------------------------------------------------------------
c
      integer   MXCELLS
      integer   MXLAYER
      integer   MXSPEC
      integer   MXREACT
      integer   MXGRID
      integer   MXPTSRC
      integer   MXCHDRN
      integer   NOVERLP
c
      parameter ( MXCELLS =      200 )
      parameter ( MXLAYER =       30 )
      parameter ( MXSPEC  =      138 )
      parameter ( MXREACT =      288 )
      parameter ( MXGRID  =       10 )
      parameter ( MXPTSRC =   100000 )
      parameter ( MXCHDRN =        3 )
      parameter ( NOVERLP =        5 )
c
c***********************************************************************
c    >>>>>>>>> MXTRSP Parameter for number of tracer species.
c              The value of this parameter greatly affects the
c              amount of memory required by the model.
c
c    OSAT>>>>> For OSAT, use the following formula to determine the
c              appropriate value:
c
c       MXTRSP = 2*ntime*nreg*(ndays+1) + 4*(1+nbound+ngroup*nreg)
c
c             where:
c
c        ntime  = number of timing releases per day (usually 0)
c        nreg   = number of source regions
c        ndays  = number of days in complete simulation
c        nbound = 1 if "stratify boundary" is off, 5 otherwise
c        ngroup = number of source groupings
c
c     For a simulation without source apportionment, set MXTRSP = 1.
c     
c   DDM>>>>>>> In the DDM, MXTRSP represents the number of species
c              times the number of parameters (MXFDDM) for which
c              sensitivities are calculated.  Use the following formula
c              for calculating MXFDDM and MXTRSP:
c
c       MXFDDM = nicddm + (nbcddm * nbound) + (nemddm * ngroup * nreg)
c                                           + nrateddm + nhddm
c
c       MXTRSP = nspec * MXFDDM
c
c             where:
c
c       nspec  = total number of modeled species
c       nicddm = number of IC DDM species groups
c       nbcddm = number of BC DDM species groups
c       nbound = 1 if "stratify boundary" is off, 5 otherwise
c       nemddm = number of emission DDM species groups
c       nreg   = number of source regions
c       ngroup = number of source groupings
c       nrateddm = number of rate constant groups
c       nhddm  = number of HDDM sensitivity groups
c
c   PA>>>>>>>> For Process Analysis, MXTRSP must be greater than or
c              equal to the parameter MXCPA which is set in procan.com
c              The current value for MXCPA in CAMx is 99
c
c***********************************************************************
c
c   MXTRSP   I   maximum number of total tracer species in simulation
c
      integer   MXTRSP
c
      parameter( MXTRSP =    1 )
c
c-----------------------------------------------------------------------
c      Parameters for input met files:
c
c      MXVAR     -- max number of variables on file
c      NUM3DMET  -- number of 3D met variables
c      NUM2DMET  -- number of 2D met variables
c      NUMKVMET  -- number of 3d Kv variables
c      NUMCLDMET -- number of 3D cloud variables
c
c-----------------------------------------------------------------------
c
      integer MXVAR
      integer NUM3DMET
      integer NUM2DMET
      integer NUMKVMET
      integer NUMCLDMET
c
      parameter ( MXVAR     = 100 )
      parameter ( NUM3DMET  = 6 )
      parameter ( NUM2DMET  = 1 )
      parameter ( NUMKVMET  = 1 )
      parameter ( NUMCLDMET = 5 )
c
c-----------------------------------------------------------------------
c      Parameters for Chemistry:
c
c      MXPHT1  -- number of primary photolysis reactions
c      MXPHT2  -- number of secondary photolysis reactions
c      NTEMPR  -- number of temperature bins for chemical reactions
c      TEMPLO  -- minimum temperature for chemical reactions
c      TEMPHI  -- maximum temperature for chemical reactions
c      NPRESR  -- number of pressure bins for chemical reactions
c      PRESLO  -- minimum pressure for chemical reactions
c      PRESHI  -- maximum pressure for chemical reactions
c      MXBIN   -- number of aerosol size bins (CF)
c      MXSECT  -- maximum number of particle sections (CMU)
c
c-----------------------------------------------------------------------
c    
      integer   MXPHT1
      integer   MXPHT2
      integer   NTEMPR
      real      TEMPLO
      real      TEMPHI
      integer   NPRESR
      real      PRESLO
      real      PRESHI
      integer   MXBIN
      integer   MXSECT
c
      parameter ( MXPHT1  =    27 )
      parameter ( MXPHT2  =    31 )
      parameter ( NTEMPR  =   141 )
      parameter ( TEMPLO  = 180.0 )
      parameter ( TEMPHI  = 330.0 )
      parameter ( NPRESR  =    11 )
      parameter ( PRESLO  =  60.0 )
      parameter ( PRESHI  = 1100. )
      parameter ( MXBIN   =     1 )
      parameter ( MXSECT  =    10 )
c
c-----------------------------------------------------------------------
c     Parameters for various lookup tables:
c
c     NOZN   -- number of ozone column values
c     NALB   -- number of albedo values
c     NHAZE  -- number of haze values
c     NRUF   -- number of roughness values
c     NLUW89 -- number of land use categories for Wesely 89 scheme
c     NLUZ03 -- number of land use categories for Zhang 03 schem
c     NHGHT  -- number of height values
c     NZEN   -- number of zenith values
c
c-----------------------------------------------------------------------
c
      integer   NOZN
      integer   NALB
      integer   NHAZE
      integer   NRUF
      integer   NLUW89
      integer   NLUZ03
      integer   NHGHT
      integer   NZEN
c
      parameter ( NOZN  =  5 )
      parameter ( NALB  =  5 )
      parameter ( NHAZE =  3 )
      parameter ( NRUF  =  9 )
      parameter ( NLUW89   = 11 )
      parameter ( NLUZ03   = 26 )
      parameter ( NHGHT = 11 )
      parameter ( NZEN  = 10 )
c
c-----------------------------------------------------------------------
c     Parameters for PiG
c
c     GRESPIG  -- flag for determining if GREASD-PiG is used
c     IRONPIG  -- flag for determining if IRON-PiG is used
c     MXPIG    -- maximum number of PiG puffs allowed in the model
c     MXRECTR  -- maximum number of IRON puff reactors
c     FLEAK    -- fraction of horizontal scale to start IRON puff leakage
c     LEAKON   -- toggle to enable/disable IRON puff leakage
c     LVISPIG  -- toggle for IRON PiG visualization in average file
c     OVERLAP  -- toggle to enable/disable puff overlap chemistry
c     DXYMAX   -- maximum allowed IRON puff size (m) before dumping
c                 If zero, the grid scale is used
c     PIGMXDT  -- maximum allowed duration of release (s)
c     AGEMAX   -- maximum puff age (s), after which puff will be slaughtered
c     MXSAMPLE -- maximum number of PiG sampling grids
c     MXCOLSMP -- maximum number of PiG sampling grid columns
c     MXROWSMP -- maximum number of PiG sampling grid rows
c----------------------------------------------------------------------
c    
      integer   GRESPIG
      integer   IRONPIG
      integer   MXPIG
      integer   MXRECTR
      integer   MXSAMPLE
      integer   MXCOLSMP
      integer   MXROWSMP
      logical   LEAKON
      logical   LVISPIG
      logical   OVERLAP
      real      FLEAK
      real      DXYMAX
      real      AGEMAX
      real      PIGMXDT
c     
      parameter ( GRESPIG  = 1 )
      parameter ( IRONPIG  = 2 )
      parameter ( MXPIG    = 50000 )
      parameter ( MXRECTR  = 1 )
      parameter ( FLEAK    = .25 )
      parameter ( LEAKON   = .FALSE. )
      parameter ( LVISPIG  = .FALSE. )
      parameter ( OVERLAP  = .FALSE. )
      parameter ( DXYMAX   = 0. )
      parameter ( AGEMAX   = 18.*3600. )
      parameter ( PIGMXDT  = 300. )
      parameter ( MXSAMPLE = 1 )
      parameter ( MXCOLSMP = 1 )
      parameter ( MXROWSMP = 1 )
c
c-----------------------------------------------------------------------
c     Parameters for Edge cell definition:
c
c   IDNORTH  -- index of the North edge in boundary condition array
c   IDSOUTH  -- index of the South edge in boundary condition array
c   IDWEST   -- index of the West edge in boundary condition array
c   IDEAST   -- index of the East edge in boundary condition array
c
      integer IDNORTH
      integer IDSOUTH
      integer IDWEST
      integer IDEAST
c
      parameter( IDNORTH = 1 )
      parameter( IDSOUTH = 2 )
      parameter( IDWEST  = 3 )
      parameter( IDEAST  = 4 )
