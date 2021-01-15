
 do is=1,NSO4
    itsp=itsp+1
    iapm = ip_sulf(k,is,ne)
!    call getvalue_apm( myid,apm_sulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                      &,i,j,atestS(k,itsp) )
    call getvalue_apm( myid,apm_sulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
                      &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )
 
 enddo


 do is=1,NSEA
    itsp=itsp+1
    iapm = ip_salt(k,is,ne)
    call getvalue_apm( myid,apm_salt(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                      &,i,j,atestS(k,itsp) )
                      &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )

 enddo


 do is=1,NDSTB
    itsp=itsp+1
    iapm = ip_dust(k,is,ne)
    call getvalue_apm( myid,apm_dust(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                      &,i,j,atestS(k,itsp) )
                      &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )
 enddo


 do is=1,NBCOCT
    itsp=itsp+1
    iapm = ip_bcoc(k,is,ne)
    call getvalue_apm( myid,apm_bcoc(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                      &,i,j,atestS(k,itsp) )
                      &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )
 enddo


 itsp=itsp+1
 iapm=ip3mem(k,ne)
 call getvalue_apm( myid,msltsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                  &,i,j,atestS(k,itsp) )
                  &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )

 itsp=itsp+1
 iapm=ip3mem(k,ne)
 call getvalue_apm( myid,mdstsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                  &,i,j,atestS(k,itsp) )
                  &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )


 itsp=itsp+1
 iapm=ip3mem(k,ne)
 call getvalue_apm( myid,mbcsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                  &,i,j,atestS(k,itsp) )
                  &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )

 itsp=itsp+1
 iapm=ip3mem(k,ne)
 call getvalue_apm( myid,mocsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                  &,i,j,atestS(k,itsp) )
                  &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )

 itsp=itsp+1
 iapm=ip3mem(k,ne)
 call getvalue_apm( myid,h2so4_gas(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                  &,i,j,atestS(k,itsp) )
                  &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )



 do is=1,nbincb
    itsp=itsp+1
    iapm = ip_cbbin(k,is,ne)
    call getvalue_apm( myid,apm_binbc(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                      &,i,j,atestS(k,itsp) )
                  &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )

 enddo


 do is=1,nbincb
    itsp=itsp+1
    iapm = ip_cbbin(k,is,ne)
    call getvalue_apm( myid,apm_binoc(iapm),sx(ne),ex(ne),sy(ne),ey(ne) &
!                      &,i,j,atestS(k,itsp) )
                  &,IPOLARMRK(4,IPS),IPOLARMRK(5,IPS),ATESTS(K,ITSP) )
 enddo







