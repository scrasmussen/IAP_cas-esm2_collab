
subroutine apm_pso4_so2( myid, lapm, dt, i03, ixy, sulfate_aqp )

 use apm_varlist
 implicit none
 include 'apm_parm.inc'
 integer :: myid
 real    :: dt
 logical :: lapm
 integer :: ixy,i02,i03,iapm
 real    :: sulfate_aqp

 !IF(lapm) THEN ! apm flag

  pso4_so2(i03+ixy) = sulfate_aqp

 !ENDIF ! apm flag


end subroutine apm_pso4_so2





