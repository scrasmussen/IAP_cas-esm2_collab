
subroutine dry_dep_apm_sink(myid,deltz,dtddep,vdep,con)

 implicit none
 integer :: myid
 real    :: deltz,dtddep,vdep,con

 con=con-vdep*con*dtddep/deltz
 con=amax1(con,0.0)

end subroutine dry_dep_apm_sink
