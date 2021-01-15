#include <define.h>

module paramodel

!----------------------------------------------------------------------
! Define the dimension of model array
!----------------------------------------------------------------------

      integer nl_soil       ! number of soil layers
      integer maxsnl        ! max number of snow layers
      integer nfcon_col     ! number of time constant variables
      integer nfcon_pft     ! number of time constant variables
      integer nftune        ! number of clm tunable constants
      integer nfvar_col     ! number of time varying variables
      integer nfvar_pft     ! number of time varying variables
      integer nforc         ! number of forcing variables
      integer nflai         ! number of leaf time varying variables
      integer nflux         ! number of LSM flux variables
      integer maxpatch      ! maximum number of patches in model grid
      integer nsoilcateg    ! number of soil texture categories
      integer nsoilcolor    ! number of soil color categories
      integer nlandcateg    ! number of land cover categories
      integer grasscateg    ! land cover category of grass type
      integer oceancateg    ! land cover category of ocean type
      integer watercateg    ! land cover category of inland water type
      integer baresoilcateg
      integer lakecateg
      integer wetlandcateg
      integer glaciercateg
      integer urbancateg

#ifdef DGVM
      integer numpft        ! number of PFTs (excluding soil)
      integer numpft_nat    ! number of natural PFTs (excluding soil & crop)
#endif

      parameter(nl_soil    = 10)
      parameter(maxsnl     = -5)
      parameter(nfcon_col  = 87)

#ifndef DGVM
      parameter(nfcon_pft  = 35+8)
      parameter(nfvar_col  = 81)
      parameter(nfvar_pft  = 45)
#else
      parameter(numpft     = 16)
      parameter(numpft_nat = 14)
#ifndef DyN
#ifdef BNUDGVM
      parameter(nfcon_pft  = 34+40+8)    ! BNUDGVM jidy@09/Sep/2014
#else
      parameter(nfcon_pft  = 34+40+17+8) ! IAPDGVM jidy@09/Sep/2014
#endif
      parameter(nfvar_col  = 81+3+1)     !IAPDGVM wliq6mon
      parameter(nfvar_pft  = 45+68+2)    ! fire
#else 
#ifdef BNUDGVM
      parameter(nfcon_pft  = 34+40+2+8)    ! BNUDGVM jidy@09/Sep/2014
#else
      parameter(nfcon_pft  = 34+40+2+17+8) ! IAPDGVM jidy@09/Sep/2014
#endif
      parameter(nfvar_col  = 81+3+6+1)
      parameter(nfvar_pft  = 45+68+17+2)
#endif
#endif

      parameter(nftune     = 14)
      parameter(nforc      = 18)
      parameter(nflai      =  4)

      parameter(grasscateg    = 13)
      parameter(baresoilcateg = 17)
      parameter(lakecateg     = 18)
      parameter(wetlandcateg  = 19)
      parameter(glaciercateg  = 20)
      parameter(urbancateg    = 21)
      parameter(oceancateg    = 22)
      parameter(nlandcateg    = 22)

      parameter(maxpatch   = nlandcateg)
      parameter(nsoilcateg = 17)
      parameter(nsoilcolor = 20)

end module paramodel
