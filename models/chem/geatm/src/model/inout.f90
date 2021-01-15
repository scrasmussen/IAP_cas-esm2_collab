  subroutine write2d( myid, a, sx, ex, sy, ey, k, id)
  integer myid, sx, ex, sy, ey, k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(10+id)myid,k,sx,ex,sy,ey
  write(10+id)((a(i,j), i = sx,ex),j=sy,ey)
  return
  end

  subroutine write2d_v2( myid, a, sx, ex, sy, ey, k, id,funit)
  integer myid, sx, ex, sy, ey, k,funit
  real a(sx-1:ex+1,sy-1:ey+1)

  write(funit)myid,k,sx,ex,sy,ey
  write(funit)((a(i,j), i = sx,ex),j=sy,ey)
  return
  end



 subroutine writeterm( myid, a, sx, ex, sy, ey, k, ip, ig, id)
  integer myid, sx, ex, sy, ey, k ,ig, ip
  real a(sx-1:ex+1,sy-1:ey+1)

  write(70+id)myid,k,ip,ig,sx,ex,sy,ey
  write(70+id)((a(i,j), i = sx,ex),j=sy,ey)
  return
  end

 subroutine writeterm_v2( myid, a, sx, ex, sy, ey, k, ip, ig,id,funit)
  integer myid, sx, ex, sy, ey, k ,ig, ip,funit
  real a(sx-1:ex+1,sy-1:ey+1)

  write(funit)myid,k,ip,ig,sx,ex,sy,ey
  write(funit)((a(i,j), i = sx,ex),j=sy,ey)
  return
  end





! CCC AQUEOUS WET

  subroutine writewet( myid, a, sx, ex, sy, ey, k,id)
  integer myid, sx, ex, sy, ey ,k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(290+id) myid,k,sx,ex,sy,ey
  write(290+id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end

  subroutine writewet_v2( myid, a, sx, ex, sy, ey, k,id,funit)
  integer myid, sx, ex, sy, ey ,k,funit
  real a(sx-1:ex+1,sy-1:ey+1)

  write(funit) myid,k,sx,ex,sy,ey
  write(funit)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end


  subroutine writedry_v2( myid, a, sx, ex, sy, ey,isp,id,funit)
  integer myid, sx, ex, sy, ey ,k,funit,isp
  real a(sx-1:ex+1,sy-1:ey+1)

  write(funit) myid,isp,sx,ex,sy,ey
  write(funit)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end


!CCC END


! CCC DUST

  subroutine writedust( myid, a, sx, ex, sy, ey, k,id)
  integer myid, sx, ex, sy, ey ,k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(400+id) myid,k,sx,ex,sy,ey
  write(400+id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end

  subroutine writedust_v2( myid, a, sx, ex, sy, ey, k,id,funit)
  integer myid, sx, ex, sy, ey ,k,funit
  real a(sx-1:ex+1,sy-1:ey+1)

  write(funit) myid,k,sx,ex,sy,ey
  write(funit)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end


!CCC END


! CCC SEA SALT

  subroutine writesea( myid, a, sx, ex, sy, ey, k,id)
  integer myid, sx, ex, sy, ey ,k
  real a(sx-1:ex+1,sy-1:ey+1)

  write(500+id) myid,k,sx,ex,sy,ey
  write(500+id)((a(i,j), i = sx,ex),j=sy,ey)

  return
  end

!CCC END


  subroutine read2d( myid, a, sx, ex, sy, ey, &
                     nx,ny,irec,id)
  integer myid, sx, ex, sy, ey, id
  real a(sx-1:ex+1,sy-1:ey+1), v0(nx,ny)
!  print*,irec,id
  read(id,rec=irec)((v0(i,j),i=1,nx),j=1,ny)
  irec=irec+1

!print*,'kk-irec=',irec

  do i= sx,ex
  do j= sy,ey
!  do i= sx-1,ex+1
!  do j= sy-1,ey+1
  !a(i,j)=((i-nx/2)**2+(j-ny/2)**2)*k
   a(i,j)=v0(i,j)
!   print*,a(i,j),i,j
  enddo
  enddo
! print*,a(47,33)
  return
  end


  subroutine puto3( myid, a, sx, ex, sy, ey, &
                     nx,ny,value)
  integer myid, sx, ex, sy, ey
  real a(sx-1:ex+1,sy-1:ey+1), value

  do i= sx,ex
  do j= sy,ey
   a(i,j)=value
  enddo
  enddo

  return
  end


  subroutine openfil1( myid,nx,ny,irec80,irec60,id, &
                        iyear,imonth,iday, ihour)
  integer myid,id, irec80,irec60,iyear,imonth,ihour
  character*1 cdnum
  character*20 fname
  character*10 date
  logical :: lexist

  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  if ( myid .lt. 0 ) then
      close( 11 )
      return
  endif
  write(cdnum(1:1),'(i1)')id
  write(fname(1:3),'(i3.3)') myid

 inquire(file=fname,exist=lexist)
 if(lexist) then ! initial data file exist 
  open( 80+id, file='init/testd'//cdnum//'.'// &
        date(1:10)//'.grd', form='unformatted',  &
        access='direct',recl=nx*ny, status='old')
  if(myid.eq.0)  print *, 80+id,'init/testd'//cdnum//'.'//date//'.grd', ' opened'
  irec80=1

  open (60+id,file='init/dustd'//cdnum//'.'//&
         date(1:10)//'.grd',form='unformatted',&
         access = 'direct',recl=nx*ny,status='old')
  if(myid.eq.0)  print*,60+id, 'init/dustd'//cdnum//'.'//date//'.grd',' opened'
  irec60 = 1  
  else

  irec80 = -1 
  irec60 = -1
  end if
  return   
  end

 subroutine openfileEmit( myid,nx,ny,irecg,id, &
                        iyear,imonth,iday, ihour)
  integer myid,irecg,id
  character*1 cdnum
  character*20 fname
  character*10 date
  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  if ( myid .lt. 0 ) then
      close( 11 )
      return
  endif
  write(cdnum(1:1),'(i1)')id
  write (fname(1:3),'(i3.3)') myid
  close(50+id)

   open(50+id,file='emit/data.emit/emitgrid_'//date(1:4)//'-'//date(5:6)//'-'//&
        date(7:8)//'_'//date(9:10)//'.d'//cdnum,form='unformatted', &
        access='direct',recl=nx*ny)

  if(myid==0)print *,  &
      'opening emit/emitgrid_'//date(1:4)//'-'//date(5:6)//'-'//&
                  date(7:8)//'.d'//cdnum
  irecg=1
  return
  end


  subroutine openfil2( myid,nx,ny,irec,irec_as,irecglobal,irechgt,id,&! irectop added by lijie 05-06-3
                       iyear,imonth,iday,ihour,ikosaline)          ! ireckv added by lijie 050712
                                                           ! irecpbl and  irecxmol added by lijie 071203 ACM2
  integer myid,irec,irec_as,id,irecgtop,ireckv,irecglobal,irecpbl
  character*1 cdnum
  character*10 date
  
  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  write(cdnum(1:1),'(i1)')id
  irec=1

!  juanxiong he
!  if(ikosaline==2)then
!  close(40+id)
!  open(40+id, file='/name8/kosa/estimate/data/'// &
!       'emitkosa.'//date(1:10), &
!        form='unformatted',access='direct',recl=nx*ny)
!  endif
!  irec_as = 1

 !---------------------global condition----------------
  close(300+id)
  open(300+id,file='global/g2m/g2md0'//cdnum//'_'//date(1:4)//'-'//date(5:6)//&
       '-'//date(7:8)// '_'//date(9:10)//'.DAT',form='unformatted',&
       access='direct',recl=nx*ny )
  irecglobal=1
!-----------------------------------------------------
                                
                      
 close(1000)

! juanxiong he
! open(1000,file='datagrid/wrfd0'//cdnum//'.dat',form='unformatted',&
!                     access='direct',recl=nx*ny)
! irechgt=2               


  return
  end
  
subroutine openfil3(myid,nx,ny,irecg,id,&
                      iyear,imonth,iday, ihour)
  integer myid,irec,id, irec80
  character*1 cdnum
  character*20 fname
  character*10 date
  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  if ( myid .lt. 0 ) then
      close( 11 )
      return
  endif
  write(cdnum(1:1),'(i1)')id
  write (fname(1:3),'(i3.3)') myid

  if(1==2) then
  close(10+id)
  call system("mkdir -p out/tmp")
  call system("mkdir -p out/tmp/"//date(1:8))
  open( 10+id, file='out/tmp/'//date(1:8)// &
        '/food'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
  !lijie change tmp to tmp1 050622 for run 2 apops.new at the same time
  endif

  close(70+id)
  call system("mkdir -p term/tmp")
  call system("mkdir -p term/tmp/"//date(1:8))
  open( 70+id, file='term/tmp/'//date(1:8)// &
         '/term'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
  return
  end

subroutine openfil4(myid,nx,ny,irecg,id,&
                      iyear,imonth,iday, ihour)
  integer irec,id, irec80
  character*1 cdnum
  character*20 fname
  character*10 date
  write(date(1:4),'(i4)')iyear
  write(date(5:6),'(i2.2)')imonth
  write(date(7:8),'(i2.2)')iday
  write(date(9:10),'(i2.2)')ihour

  if ( myid .lt. 0 ) then
      close( 11 )
      return
  endif
  write(cdnum(1:1),'(i1)')id
  write (fname(1:3),'(i3.3)') myid
  close(10+id)
  call system("mkdir -p out/tmp")
  call system("mkdir -p out/tmp/"//date(1:8))
  open( 10+id, file='out/tmp/'//date(1:8)// &
        '/food'//cdnum//'.'//fname(1:3)//'.'//date,&
        form='unformatted',access='direct',recl=nx*ny)
  irecg=1
  return
  end


! CCCC DUST CONCENTRATION AND EMISSIONS, DRY, WET and GRAVITY DEPOSITION

  subroutine openfildust(myid,nx,ny,irecdust,id,iyear,imonth,iday,ihour)
   integer myid,irecdust,id
   character*1 cdnum
   character*20 fname
   character*10 date
    write(date(1:4),'(i4)')iyear
    write(date(5:6),'(i2.2)')imonth
    write(date(7:8),'(i2.2)')iday
    write(date(9:10),'(i2.2)')ihour

    write(cdnum(1:1),'(i1)')id
    write (fname(1:3),'(i3.3)') myid
    close(400+id)
    call system("mkdir -p dust/tmp")
    call system("mkdir -p dust/tmp/"//date(1:8))
    open( 400+id, file='dust/tmp/'//date(1:8)// &
        '/dustd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
    irecdust = 1
    return
  end

! CCC END CCC


! CCCC SEA CONCENTRATION AND EMISSIONS, DRY, WET and GRAVITY DEPOSITION

  subroutine openfilsea(myid,nx,ny,irecsea,id,iyear,imonth,iday,ihour)
   integer myid,irecsea,id
   character*1 cdnum
   character*20 fname
   character*10 date
    write(date(1:4),'(i4)')iyear
    write(date(5:6),'(i2.2)')imonth
    write(date(7:8),'(i2.2)')iday
    write(date(9:10),'(i2.2)')ihour

    write(cdnum(1:1),'(i1)')id
    write (fname(1:3),'(i3.3)') myid
    close(500+id)
    call system("mkdir -p seasalt/tmp")
    call system("mkdir -p seasalt/tmp/"//date(1:8))
    open( 500+id, file='seasalt/tmp/'//date(1:8)// &
        '/seasaltd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
    irecsea = 1
    return
  end

! CCC END CCC



!CCCC   AQUEOUS DEPOSITION  CCCCCCCCCC
  subroutine openfilwet(myid,nx,ny,irecwet,id,iyear,imonth,iday,ihour)
   integer myid,irecwet,id
   character*1 cdnum
   character*20 fname
   character*10 date
    write(date(1:4),'(i4)')iyear
    write(date(5:6),'(i2.2)')imonth
    write(date(7:8),'(i2.2)')iday
    write(date(9:10),'(i2.2)')ihour

    write(cdnum(1:1),'(i1)')id
    write (fname(1:3),'(i3.3)') myid
    close(290+id)
    call system("mkdir -p wetdep/tmp")
    call system("mkdir -p wetdep/tmp/"//date(1:8))
    open( 290+id, file='wetdep/tmp/'//date(1:8)// &
        '/wetdepd'//cdnum//'.'//fname(1:3)//'.'//date,&
          form='unformatted',status='UNKNOWN' )
    irecwet = 1
    return
  end
! CCC END CCC

    
  subroutine needSNDwrite( a, sx, ex, sy, ey, k,b,ntotalnum,nzz)
  integer sx, ex, sy, ey, k,  ntotalnum
  real a(sx-1:ex+1,sy-1:ey+1),b(ntotalnum,nzz)

  if(ntotalnum .ne. (ex-sx+1)*(ey-sy+1))then
    print *,'number setting error'
   stop
  endif

  ii=1
  do i=sx,ex
  do j=sy,ey
    b(ii,k)=a(i,j)
    ii=ii+1
    enddo
    enddo
  return
  end


  subroutine needGETwrite(a,nx,ny,sx,ex,sy,ey,k,b,ntotalnum,nzz)
  integer sx, ex, sy, ey, k, ntotalnum
  real a(nx,ny,nzz),b(ntotalnum,nzz)

  if(ntotalnum .ne. (ex-sx+1)*(ey-sy+1))then
    print *,'number setting error'
   stop
  endif

  ii=1
  do i=sx,ex
  do j=sy,ey
    a(i,j,k)=b(ii,k)
    ii=ii+1
    enddo
    enddo
  return
  end


  subroutine putvalue( myid, a, sx, ex, sy, ey, &
                     ix,iy,value)
  integer myid, sx, ex, sy, ey,ix,iy
  real a(sx-1:ex+1,sy-1:ey+1), value
   a(ix,iy)=value
  return
  end
  

  subroutine getvalue( myid, a, sx, ex, sy, ey, &
                     ix,iy,value)
  integer myid, sx, ex, sy, ey, ix, iy
  real a(sx-1:ex+1,sy-1:ey+1), value
  
  value=a(ix,iy)

  return
  end

  subroutine write2dconv( myid, a, sx, ex, sy, ey, k, id)
  integer myid, sx, ex, sy, ey, k
  integer a(sx-1:ex+1,sy-1:ey+1)
     
  write(10+id)myid,k,sx,ex,sy,ey
  write(10+id)((float(a(i,j)), i = sx,ex),j=sy,ey)
  return
  end

  subroutine read2dconv( myid, a, sx, ex, sy, ey, &
                       nx,ny,irec,id)
  integer myid, sx, ex, sy, ey
  integer a(sx-1:ex+1,sy-1:ey+1)
  real v0(nx,ny)
   read(id,rec=irec)((v0(i,j),i=1,nx),j=1,ny)
   irec=irec+1
   do j= sy,ey
   do i= sx,ex
     a(i,j)=int(v0(i,j))
   enddo
   enddo
     return
  end

