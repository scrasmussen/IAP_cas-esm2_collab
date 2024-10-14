!  CVS: $Id: setbcx.F90,v 1.1.1.1 2004/04/29 06:22:39 lhl Exp $
#include <def-undef.h>
 
#if (defined ISO)
!     ===================================
      SUBROUTINE setbcx (a, imttmp, jmtorkmtmp)
!     ===================================
use precision_mod
use param_mod
      IMPLICIT NONE
      INTEGER :: imttmp,jmtorkmtmp,ktmp
      REAL(r8) :: a(imttmp,jmtorkmtmp)
      if (nx_proc==1) then
      DO ktmp = 1,jmtorkmtmp
         a (1,ktmp) = a (imttmp -1,ktmp)
         a (imttmp,ktmp) = a (2,ktmp)
      END DO
      else
      DO ktmp=1,jmtorkmtmp
      call exchange_boundary(a(1,ktmp),1)
      END DO
      end if
 
      RETURN
      END SUBROUTINE setbcx
 
#else
      SUBROUTINE setbcx ()
      RETURN
      END SUBROUTINE setbcx
#endif 
 
