

! these codes are in do k loop in main.f90

if(lfor_sulf) then
 do is=1,NSO4
  iapm=ip_sulf(k,is,ne)
  call exchng2( myid, apm_sulf(iapm) &
               ,sx(ne), ex(ne), sy(ne), ey(ne) &
               ,comm2d(ne), stride(ne) &
               ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
 enddo
endif

if(lfor_salt) then
 do is=1,NSEA
  iapm=ip_salt(k,is,ne)
  call exchng2( myid, apm_salt(iapm) &
               ,sx(ne), ex(ne), sy(ne), ey(ne) &
               ,comm2d(ne), stride(ne) &
               ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
 enddo
endif

if(lfor_dust) then
 do is=1,NDSTB
  iapm=ip_dust(k,is,ne)
  call exchng2( myid, apm_dust(iapm) &
               ,sx(ne), ex(ne), sy(ne), ey(ne) &
               ,comm2d(ne), stride(ne) &
               ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
 enddo
endif

if(lfor_bcoc) then
 do is=1,NBCOCT
  iapm=ip_bcoc(k,is,ne)
  call exchng2( myid, apm_bcoc(iapm) &
               ,sx(ne), ex(ne), sy(ne), ey(ne) &
               ,comm2d(ne), stride(ne) &
               ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
 enddo
endif


!================
!> coated sulfate
if(lcoated_dyn) then

 iapm=ip3mem(k,ne)

 call exchng2( myid, msltsulf(iapm) &
              ,sx(ne), ex(ne), sy(ne), ey(ne) &
              ,comm2d(ne), stride(ne) &
              ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
 
 call exchng2( myid, mdstsulf(iapm) &
              ,sx(ne), ex(ne), sy(ne), ey(ne) &
              ,comm2d(ne), stride(ne) &
              ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )

 call exchng2( myid, mbcsulf(iapm)  &
              ,sx(ne), ex(ne), sy(ne), ey(ne) &
              ,comm2d(ne), stride(ne) &
              ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )

 call exchng2( myid, mocsulf(iapm)  &
              ,sx(ne), ex(ne), sy(ne), ey(ne) &
              ,comm2d(ne), stride(ne) &
              ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
endif
!!!!!


if(lfor_h2so4) then
 iapm=ip3mem(k,ne)
 call exchng2( myid, h2so4_gas(iapm) &
              ,sx(ne), ex(ne), sy(ne), ey(ne) &
              ,comm2d(ne), stride(ne) &
              ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
endif


 do is=1,nbincb
  iapm=ip_cbbin(k,is,ne)
  call exchng2( myid, apm_binbc(iapm) &
               ,sx(ne), ex(ne), sy(ne), ey(ne) &
               ,comm2d(ne), stride(ne) &
               ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
 enddo

 do is=1,nbincb
  iapm=ip_cbbin(k,is,ne)
  call exchng2( myid, apm_binoc(iapm) & 
               ,sx(ne), ex(ne), sy(ne), ey(ne) &
               ,comm2d(ne), stride(ne) &
               ,nbrleft(ne), nbrright(ne), nbrtop(ne), nbrbottom(ne) )
 enddo





