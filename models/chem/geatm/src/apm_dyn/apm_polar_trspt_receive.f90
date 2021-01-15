
do is=1,NSO4
    ips1 = ips1 + 1
    iapm = ip_sulf(k,is,ne)
    call putvalue_apm(myid,apm_sulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
enddo

do is=1,NSEA
    ips1 = ips1 + 1
    iapm = ip_salt(k,is,ne)
    call putvalue_apm(myid,apm_salt(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
enddo

do is=1,NDSTB
    ips1 = ips1 + 1
    iapm = ip_dust(k,is,ne)
    call putvalue_apm(myid,apm_dust(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
enddo

do is=1,NBCOCT
    ips1 = ips1 + 1
    iapm = ip_bcoc(k,is,ne)
    call putvalue_apm(myid,apm_bcoc(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
enddo



ips1 = ips1 + 1
iapm = ip3mem(k,ne)
call putvalue_apm(myid,msltsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))


ips1 = ips1 + 1
iapm = ip3mem(k,ne)
call putvalue_apm(myid,mdstsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))

ips1 = ips1 + 1
iapm = ip3mem(k,ne)
call putvalue_apm(myid,mbcsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))

ips1 = ips1 + 1
iapm = ip3mem(k,ne)
call putvalue_apm(myid,mocsulf(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))

ips1 = ips1 + 1
iapm = ip3mem(k,ne)
call putvalue_apm(myid,h2so4_gas(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))



do is=1,nbincb
    ips1 = ips1 + 1
    iapm = ip_cbbin(k,is,ne)
    call putvalue_apm(myid,apm_binbc(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
enddo


do is=1,nbincb
    ips1 = ips1 + 1
    iapm = ip_cbbin(k,is,ne)
    call putvalue_apm(myid,apm_binoc(iapm),sx(ne),ex(ne),sy(ne),ey(ne), &
            &    IPOLARMRK(3,IPS),IPP,atestR(k,ips1))
enddo

