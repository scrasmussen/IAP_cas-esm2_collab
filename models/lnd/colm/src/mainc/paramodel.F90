#include <define.h>

module paramodel

!----------------------------------------------------------------------
! Define the dimension of model array
!----------------------------------------------------------------------

      integer nl_soil       ! number of soil layers
      integer nl_lake       ! number of lake layers
      integer nl_csoil      ! number of soil layers for solving carbon
      integer maxsnl        ! max number of snow layers
      integer nftune        ! number of clm tunable constants
      integer nforc         ! number of forcing variables
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

      integer npftpara      ! number of parameters for each PFT
      integer numlandc
#ifdef IAPDGVM
      parameter(npftpara   = 49)
#else
      parameter(npftpara   = 32)
#endif    
      parameter(numlandc   = 21)

      integer numpft        ! number of PFTs (excluding soil)
      integer numpft_nat    ! number of natural PFTs (excluding soil & crop)

#ifdef SOIL10
      parameter(nl_soil    = 10)   !10L
#endif
#ifdef SOIL15
      parameter(nl_soil    = 15)   !15L
#endif
      parameter(nl_lake    = 10)
      parameter(maxsnl     = -5)

      parameter(nl_csoil   = 10)

      parameter(numpft     = 16)
      parameter(numpft_nat = 14)

      parameter(nftune     = 14)
      parameter(nforc      = 18)

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
