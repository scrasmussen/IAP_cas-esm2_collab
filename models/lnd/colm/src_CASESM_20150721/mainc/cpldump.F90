
#include <define.h>

#ifdef COUP_CSM
#ifdef CPL_DUMP

module cpldump

   use precision
   use paramodel,   only: nforc
   use colm_varMod, only: numgrid, numgrid_glob
   use colm_cplMod, only: forcg
   use spmd
   use spmd_decomp

   implicit none

   private

   integer, parameter :: nrecv  = nforc
   integer, parameter :: nsend  = 19

   integer, parameter :: iurecv = 31
   integer, parameter :: iusend = 32

   character(len=255) :: fnrecv = "recv.dat"
   character(len=255) :: fnsend = "send.dat"

   real(r8), pointer  :: vrecv(:,:)
   real(r8), pointer  :: vsend(:,:)

   integer, save      :: irecv = 0
   integer, save      :: isend = 0

   logical, save      :: readmode = .false.

   interface cpldump_init
      module procedure cpldump_init
   end interface cpldump_init

   interface cpldump_send
      module procedure cpldump_send
   end interface cpldump_send

   interface cpldump_recv
      module procedure cpldump_recv
   end interface cpldump_recv

   interface cpldump_read
      module procedure cpldump_read
   end interface cpldump_read

   interface cpldump_exit
      module procedure cpldump_exit
   end interface cpldump_exit

   public cpldump_init, cpldump_send, cpldump_recv, cpldump_read, cpldump_exit

contains

   subroutine cpldump_init(lread)

      logical,intent(in) :: lread

      character(len=12) :: fstat
      integer length

      readmode = lread

      if (readmode) then
         fstat = "old"
      else
         fstat = "unknown"
      end if

      allocate (vrecv(nforc,numgrid_glob))
      allocate (vsend(nsend,numgrid_glob))

      if (p_master) then
         inquire(iolength=length) vrecv
         open(iurecv,file=trim(fnrecv),form='unformatted',access='direct',recl=length,status=trim(fstat))

         inquire(iolength=length) vsend
         open(iusend,file=trim(fnsend),form='unformatted',access='direct',recl=length,status='unknown')
      end if

   end subroutine cpldump_init

   subroutine cpldump_read(forcg_loc)

      real(r8), intent(out) :: forcg_loc(nforc,numgrid)

      integer i

      if (.not. readmode) return

      irecv = irecv+1

      if (p_master) read(iurecv,rec=irecv) vrecv

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do

      p_counts = p_counts*nforc
      p_displs = p_displs*nforc

      call mpi_scatterv(vrecv,p_counts,p_displs,mpi_real8,&
                        forcg_loc,p_counts(p_iam),mpi_real8,0,p_comm,p_err)
#else
      forcg_loc = vrecv
#endif

   end subroutine cpldump_read

   subroutine cpldump_recv

      integer i

      if (readmode) return

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do

      p_counts = p_counts*nrecv
      p_displs = p_displs*nrecv

      call mpi_gatherv(forcg,size(forcg),mpi_real8,&
                       vrecv,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      vrecv = forcg
#endif

      irecv = irecv+1

      if (p_master) write(iurecv,rec=irecv) vrecv

   end subroutine cpldump_recv

   subroutine cpldump_send(trad,scv,avsdr,avsdf,anidr,anidf,tref,qref,u10m,v10m,sublim,taux,tauy,lhflx,shflx,lwup,qflx,swabs,nee)

      real(r8), intent(in) :: trad(numgrid)
      real(r8), intent(in) :: scv(numgrid)
      real(r8), intent(in) :: avsdr(numgrid)
      real(r8), intent(in) :: avsdf(numgrid)
      real(r8), intent(in) :: anidr(numgrid)
      real(r8), intent(in) :: anidf(numgrid)
      real(r8), intent(in) :: tref(numgrid)
      real(r8), intent(in) :: qref(numgrid)
      real(r8), intent(in) :: u10m(numgrid)
      real(r8), intent(in) :: v10m(numgrid)
      real(r8), intent(in) :: sublim(numgrid)
      real(r8), intent(in) :: taux(numgrid)
      real(r8), intent(in) :: tauy(numgrid)
      real(r8), intent(in) :: lhflx(numgrid)
      real(r8), intent(in) :: shflx(numgrid)
      real(r8), intent(in) :: lwup(numgrid)
      real(r8), intent(in) :: qflx(numgrid)
      real(r8), intent(in) :: swabs(numgrid)
      real(r8), intent(in) :: nee(numgrid)

      real(r8), pointer ::  vsend_loc(:,:)

      integer i

      allocate (vsend_loc(nsend,numgrid))

      vsend_loc( 1,:) = trad(:)
      vsend_loc( 2,:) = scv(:)
      vsend_loc( 3,:) = avsdr(:)
      vsend_loc( 4,:) = avsdf(:)
      vsend_loc( 5,:) = anidr(:)
      vsend_loc( 6,:) = anidf(:)
      vsend_loc( 7,:) = tref(:)
      vsend_loc( 8,:) = qref(:)
      vsend_loc( 9,:) = u10m(:)
      vsend_loc(10,:) = v10m(:)
      vsend_loc(11,:) = sublim(:)
      vsend_loc(12,:) = taux(:)
      vsend_loc(13,:) = tauy(:)
      vsend_loc(14,:) = lhflx(:)
      vsend_loc(15,:) = shflx(:)
      vsend_loc(16,:) = lwup(:)
      vsend_loc(17,:) = qflx(:)
      vsend_loc(18,:) = swabs(:)
      vsend_loc(19,:) = nee(:)

#ifdef SPMD
      p_counts(:) = numgrid_proc(:)
      p_displs(0) = 0
      do i = 1, p_nprocs-1
         p_displs(i) = sum(numgrid_proc(0:i-1))
      end do

      p_counts = p_counts*nsend
      p_displs = p_displs*nsend

      call mpi_gatherv(vsend_loc,size(vsend_loc),mpi_real8,&
                       vsend,p_counts,p_displs,mpi_real8,0,p_comm,p_err)
#else
      vsend = vsend_loc
#endif

      isend = isend+1

      if (p_master) write(iusend,rec=isend) vsend

      deallocate (vsend_loc)

   end subroutine cpldump_send

   subroutine cpldump_exit

      if (p_master) then
         close(iurecv)
         close(iusend)
      end if

      deallocate (vrecv)
      deallocate (vsend)

   end subroutine cpldump_exit

end module cpldump

#endif
#endif
