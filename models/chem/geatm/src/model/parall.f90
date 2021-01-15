  subroutine exchng2( myid, v, sx, ex, sy, ey,  &
                      comm2d, stride,  &
                      nbrleft, nbrright, nbrtop, nbrbottom  )
  include "mpif.h"
  integer myid, sx, ex, sy, ey, stride
  real v(sx-1:ex+1,sy-1:ey+1)
  integer nbrleft, nbrright, nbrtop, nbrbottom, comm2d
  integer status(MPI_STATUS_SIZE), ierr, nx
!
  nx = ex - sx + 1
!  These are just like the 1-d versions, except for less data
  call MPI_SENDRECV( v(sx,ey),  nx, MPI_REAL,  &
                    nbrtop, 0,  &
                    v(sx,sy-1), nx, MPI_REAL,  &
                    nbrbottom, 0, comm2d, status, ierr )
  call MPI_SENDRECV( v(sx,sy),  nx, MPI_REAL,  &
                    nbrbottom, 1,  &
                    v(sx,ey+1), nx, MPI_REAL,  &
                    nbrtop, 1, comm2d, status, ierr )
! This uses the "strided" datatype
!       v(ex,sy-1) = -100 - myid
  call MPI_SENDRECV( v(ex,sy-1),  1, stride, nbrright,  2,  &
                     v(sx-1,sy-1), 1, stride, nbrleft,  2,  &
                     comm2d, status, ierr )
!       v(sx,sy-1) = -200 - myid
  call MPI_SENDRECV( v(sx,sy-1),  1, stride, nbrleft,   3,  &
                     v(ex+1,sy-1), 1, stride, nbrright, 3,  &
                     comm2d, status, ierr )
  return
  end

 subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
  integer n, numprocs, myid, s, e
  integer nlocal
  integer deficit
!
  nlocal  = n / numprocs
  s       = myid * nlocal + 1
  deficit = mod(n,numprocs)
  s       = s + min(myid,deficit)
  if (myid .lt. deficit) then
      nlocal = nlocal + 1
  endif
  e = s + nlocal - 1
  if (e .gt. n .or. myid .eq. numprocs-1) e = n
  return
  end


  SUBROUTINE EXINFO(MYID,NUMPROCS,AA1,AA2,AA3,AA4,&
                   SMARK,STADOM)
  INCLUDE "mpif.h"
  INTEGER MYID,NUMPROCS,STADOM
  INTEGER AA1, AA2, AA4, AA3
  INTEGER STATUS(MPI_STATUS_SIZE),IERR

  INTEGER :: SMARK(STADOM)

  !DO ISTD=1,(4*NUMPROCS+4)  !! need to modify
  DO ISTD=1,(4*NUMPROCS)  !! modified by chenhs
    SMARK(ISTD)=0
  ENDDO

  SMARK(4*MYID+1)=AA1
  SMARK(4*MYID+2)=AA2
  SMARK(4*MYID+3)=AA3
  SMARK(4*MYID+4)=AA4

  IDSTD=4*MYID+1

  IF(MYID>0)THEN
  DO IDETC=0,(MYID-1)
  CALL MPI_SEND(SMARK(IDSTD),4,MPI_INTEGER,IDETC,0,MPI_COMM_WORLD,IERR)
  ENDDO
  ENDIF

  IF(MYID<(NUMPROCS-1))THEN
  DO IDETC=(MYID+1),(NUMPROCS-1)
  CALL MPI_SEND(SMARK(IDSTD),4,MPI_INTEGER,IDETC,0,MPI_COMM_WORLD,IERR)
  ENDDO
  ENDIF

  IF(NUMPROCS .GT. 1)CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )

  IF(MYID>0)THEN
  DO IDETC=0,(MYID-1)
  IDSTD=IDETC*4+1
  CALL MPI_RECV(SMARK(IDSTD),4,MPI_INTEGER,IDETC,0,MPI_COMM_WORLD,STATUS,IERR)
  ENDDO
  ENDIF
  IF(MYID<(NUMPROCS-1))THEN
  DO IDETC=(MYID+1),(NUMPROCS-1)
  IDSTD=IDETC*4+1
  CALL MPI_RECV(SMARK(IDSTD),4,MPI_INTEGER,IDETC,0,MPI_COMM_WORLD,STATUS,IERR)
  ENDDO
  ENDIF

  IF(NUMPROCS .GT. 1)CALL MPI_BARRIER( MPI_COMM_WORLD, IERR )

  IF(AA1<0)THEN
    ISID=SMARK(MYID*4+2)
    DO
      IF(SMARK(ISID*4+2)<0)THEN
        AA1=ISID
        GOTO 1200
      ELSE
        ISID=SMARK(ISID*4+2)
      ENDIF
    ENDDO
  ENDIF
  1200 CONTINUE

  IF(AA2<0)THEN
    ISID=SMARK(MYID*4+1)
    DO
      IF(SMARK(ISID*4+1)<0)THEN
        AA2=ISID
        GOTO 1201
      ELSE
        ISID=SMARK(ISID*4+1)
      ENDIF
    ENDDO
  ENDIF
  1201 CONTINUE

  RETURN
 END SUBROUTINE EXINFO

!---------------------------------------------------------------------------------------
 subroutine gather_whole_field(patch_array, ps1, pe1, ps2, pe2, ps3, pe3, &
                                domain_array, ds1, de1, ds2, de2, ds3, de3)
      use geatm_vartype, only : procs, local_com2d, myid2d, dims
      implicit none
      include "mpif.h"
      ! arguments
      integer, intent(in) :: ps1, pe1, ps2, pe2, ps3, pe3, &
                             ds1, de1, ds2, de2, ds3, de3
      real(8), dimension(ps1:pe1,ps2:pe2,ps3:pe3), intent(in) :: patch_array
      real(8), dimension(ds1:de1,ds2:de2,ds3:de3), intent(inout) :: domain_array
  
      ! local variables
      integer :: i, ii, j, jj, k, kk, m
      integer, dimension(2) :: idims, jdims
      integer :: mpi_ierr
      integer, dimension(mpi_status_size) :: mpi_stat
      real(8), dimension(:),allocatable ::  av
      
      if (myid2d .eq. 0) then
  
          do i=0,dims(2,1)-1
            do j=0,dims(1,1)-1
               if (procs(i,j) .ne. 0) then

                  call mpi_recv(jdims, 2, mpi_integer, procs(i,j),mpi_any_tag,local_com2d, mpi_stat, mpi_ierr)
                  call mpi_recv(idims, 2, mpi_integer, procs(i,j),mpi_any_tag,local_com2d, mpi_stat, mpi_ierr)

                  allocate(av((idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de2-ds2+1))) 
                  call mpi_recv(av,(idims(2)-idims(1)+1)*(jdims(2)-jdims(1)+1)*(de2-ds2+1), &
                                MPI_DOUBLE_PRECISION, procs(i,j), mpi_any_tag,local_com2d, mpi_stat, mpi_ierr)

                  m=0
                  do jj=jdims(1), jdims(2)
                    do kk=ds2, de2
                      do ii=idims(1), idims(2)
                        m=m+1
                        domain_array(ii,kk,jj)=av(m)
                       enddo
                    enddo
                  enddo
                  deallocate(av)

               else
                  domain_array(ps1:pe1,ps2:pe2,ps3:pe3) =patch_array(ps1:pe1,ps2:pe2,ps3:pe3)
               end if
            end do
           end do
  
      else
  
         jdims(1) = ps3
         jdims(2) = pe3
         call mpi_send(jdims, 2, mpi_integer, 0, myid2d, local_com2d, mpi_ierr)
         idims(1) = ps1
         idims(2) = pe1
         call mpi_send(idims, 2, mpi_integer, 0, myid2d, local_com2d, mpi_ierr)

         allocate(av((pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1)))
         m=0
         do jj=ps3, pe3
          do kk=ps2, pe2
           do ii=ps1, pe1 
             m = m+1
             av(m)= patch_array(ii,kk,jj)
           enddo
          enddo
         enddo
         call mpi_send(av, (pe1-ps1+1)*(pe2-ps2+1)*(pe3-ps3+1), &
                       MPI_DOUBLE_PRECISION, 0, myid2d, local_com2d, mpi_ierr)
         deallocate(av)
      end if
 
   end subroutine gather_whole_field
