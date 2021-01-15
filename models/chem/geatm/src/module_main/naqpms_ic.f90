
subroutine initialize_naqpms( myid,iyear,imonth,iday,ihour,iminute &
                  ,nest,nzz,nx,ny,nz,sx,ex,sy,ey,ne &
                  ,igas,isize,iaer &
                  ,ndustcom &
                  ,ifsmt,ifsm,idmSet,ismMax,igMark )

use naqpms_varlist
use met_fields
use naqpms_gridinfo

use  work_vars, only: tropp

implicit none

integer :: myid
integer :: ne,nest
integer :: nx(5),ny(5),nz(5),nzz
integer :: sy(5),ey(5),sx(5),ex(5)

integer :: iyear,imonth,iday,ihour,iminute


integer :: ifsm(5)

integer :: idmSet,ismMax

integer :: igMark(idmSet)

integer :: igas,isize,iaer

integer :: ndustcom

integer :: ifsmt


integer :: k,i0,ig,ia,is,idug,iduc

integer :: idm,ism

character :: cdnum*1
character :: date*10

integer :: irec
integer :: funit

logical :: lexist

character :: fname*100

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour


!-------------------- to read initial data----------------------------------
loop_nest : do ne=1,nest

 write(cdnum(1:1),'(i1)') ne

 fname='init/testd'//cdnum//'.'//date(1:10)//'.grd'

 inquire(file=fname,exist=lexist)

 if(lexist) then ! initial data file exist 

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',access='direct' &
            ,recl=nx(ne)*ny(ne),status='old')

  irec=1

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,dz(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
         
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,w(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,t(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
 
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,rh1(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,Plev(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                   nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
   i0=ip2mem(ne) 
   call read2d(myid,tropp(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                      nx(ne),ny(ne),irec,funit)
  enddo                   
  
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,h(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo

  do ig=1,igas
     do k=1,nzz
     i0=ip4mem(k,ig,ne)
     call read2d(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
     enddo
  enddo

if(1==2) then
  do ig=1,naersp
  do is=1,naerbin
  do k=1,nzz
     i0=ip4mem_aer(k,is,ig,ne)
     call read2d(myid,aerom(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
  enddo
  enddo
endif

!cycle ! shun : only read gas tracers in NAQPMS, goto next domain

!stop

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,jo1d(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo
  
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,jno2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo
  
  
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,EXT(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo              

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,EXTASO4(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,EXTANO3(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,EXTANH4(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,EXTBC(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,EXTOC(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,VISIB(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo


  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,UVB(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo
  
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,UVBS(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo
  
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,UVA(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,VIS(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo  
 
  do k=1,nzz 
  i0=ip3mem(k,ne)
  call read2d(myid,SSA(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo  

  i0=ip2mem(ne)
  call read2d(myid,AOD(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
       nx(ne),ny(ne),irec,funit)

  i0=ip2mem(ne)
  call read2d(myid,CLDOPD(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
       nx(ne),ny(ne),irec,funit)
       
  i0=ip2mem(ne)
  call read2d(myid,DUSO2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
       nx(ne),ny(ne),irec,funit)

  i0=ip2mem(ne)
  call read2d(myid,DUO3(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
       nx(ne),ny(ne),irec,funit)
        
  i0=ip2mem(ne)
  call read2d(myid,DUNO2(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
       nx(ne),ny(ne),irec,funit)

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,ANA(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo
  
  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,ASO4(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,ANH4(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo


  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,ANO3(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo


  do k=1,nzz
  i0=ip3mem(k,ne)
  call read2d(myid,ACL(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo
             
  do K = 1, NZZ
  i0=ip3mem(k,ne)
  call read2d(myid,CPH(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo

  do K = 1, NZZ
  i0=ip3mem(k,ne)
  call read2d(myid,OPE(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
        nx(ne),ny(ne),irec,funit)
  enddo
  
       
  do ia=1,iaer
  do is=1,isize
  do k=1,nzz
  i0=ip5mem(k,is,ia,ne)
     call read2d(myid,aer(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
  enddo
  enddo

  do ig=1,naersp
  do is=1,naerbin
  do k=1,nzz
     i0=ip4mem_aer(k,is,ig,ne)
     call read2d(myid,aerom(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
  enddo
  enddo

if(myid.eq.0)  print*,'naq_ic irec=',irec

!stop

  close(funit)

  cycle

! ---- read dust init data
  fname='init/dustd'//cdnum//'.'//date(1:10)//'.grd'

  inquire(file=fname,exist=lexist)

  if(lexist) then ! initial data file exist 

  call get_funitnaqpms(funit)

  open(funit,file=trim(fname),form='unformatted',access='direct' &
            ,recl=nx(ne)*ny(ne),status='old')

  irec=1
  
  do k=1,nzz
   i0=ip3mem(k,ne)
   call read2d(myid,u(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
             nx(ne),ny(ne),irec,funit)
  enddo

  do k=1,nzz
   i0=ip3mem(k,ne)
   call read2d(myid,v(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
             nx(ne),ny(ne),irec,funit)
  enddo

  do ia=2,2
  do is=1,isize
  do k=1,nzz
     i0=ip5mem(k,is,ia,ne)
     call read2d(myid,aer(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),irec,funit)
  enddo
  enddo
  enddo

  do iduc = 1, ndustcom
  do is = 1, isize
  do k =1 ,nzz
     i0=ip5memc(k,is,iduc,ne)
   call read2d(myid,dustcomp(i0),sx(ne),ex(ne),sy(ne),ey(ne),&
           nx(ne),ny(ne),irec,funit)
  enddo
  enddo
  enddo

  close(funit)

  endif

 else  ! initial data file does NOT exist

 ! to ozone as 60ppb vertically
  do ig=11,11
     do k=1,nzz
     i0=ip4mem(k,ig,ne)
      if(k==1.or.k==2.or.k==3) &
        call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                   nx(ne),ny(ne),20.)
      if(k==4.or.k==5.or.k==6) &
        call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                   nx(ne),ny(ne),10.)
      if(k==7.or.k==8.or.k==9) &
        call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),20.)
      if(k>=10)  &
        call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),30.)

     enddo
  enddo

if(1==2) then  !winne@2017.07.12,CO emit = 0
  do ig=17,17  !to CO as 200ppbv
     do k=1,nzz
     i0=ip4mem(k,ig,ne)
     call puto3(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),300.)
     enddo
  enddo
endif

  DO IA = 1, IAER  ! TO DUST AND SEA SALT 5 UG/M3
  DO IS = 1, ISIZE
  DO K =1 ,NZZ
     I0 = IP5MEM(K,IS,IA,NE)
   CALL PUTO3 (MYID, AER(I0), SX(NE), EX(NE), SY(NE), EY(NE),&
               NX(NE),NY(NE),0.1)
  ENDDO ! K
  ENDDO ! IS
  ENDDO ! IA

 endif   ! if no exist of the init file, no input

enddo loop_nest ! do nest loop

!!!!!!!!!!!!!!!!!!!!!
! For Source Mark
if(ifsmt>0)then    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do ne=1,nest

 if(ifsm(ne)==1)then !mmmmm

  write(cdnum(1:1),'(i1)') ne
  
  fname='initSM/resmd'//cdnum//'.'//date(1:10)//'.grd' 
 
  inquire(file=fname,exist=lexist)

  irec=1

  if(lexist) then !xxxxxxxxxxxxxxx

     call get_funitnaqpms(funit)

     open(funit,file=trim(fname),form='unformatted',access='direct' &
               ,recl=nx(ne)*ny(ne),status='old')     

     do idm=1,idmSet

      do ism=1,ismMax  !ism
      do k=1,nzz
      i0=ipSMmem(k,ism,idm,ne)
     call read2d(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irec,funit)
      enddo
      enddo            !ism

      do k=1,nzz
      ig=igMark(idm)
      i0=ip4mem(k,ig,ne)
      call read2d(myid,gas(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
               nx(ne),ny(ne),irec,funit)
      enddo

     enddo

  else ! irecSM(ne) .gt. 0

     do idm=1,idmSet
     do ism=1,ismMax
       do k=1,nzz-1
        i0=ipSMmem(k,ism,idm,ne)
           if(ism==3)then   ! initial 100% (1 boundary,2 strato, 3 init)
     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),1.)
           else
     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),0.)
           endif
       enddo !k
! set top conditions 100%
       do k=nzz,nzz
           i0=ipSMmem(k,ism,idm,ne)
          if(ism==2) then
     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
                    nx(ne),ny(ne),1.)
          else
     call puto3(myid,SourceMark(i0),sx(ne),ex(ne),sy(ne),ey(ne), &
              nx(ne),ny(ne),0.)
          endif
       enddo !k
     enddo !ism
     enddo !idm

  endif                   ! lexist

  endif                   ! ifsm(ne).eq.1

  enddo

endif ! if sm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine initialize_naqpms
