
subroutine apm_area(flag)

use apm_varlist
implicit none
include 'apm_parm.inc'

character :: flag*4

if(    flag.eq.'sulf') then
 return
elseif(flag.eq.'salt') then
 return
elseif(flag.eq.'dust') then
 return
elseif(flag.eq.'bc') then
 return
elseif(flag.eq.'oc') then
 return
endif

end
