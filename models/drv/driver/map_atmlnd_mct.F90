module map_atmlnd_mct

!---------------------------------------------------------------------
!
! Purpose:
!
! Collect coupling routines for sequential coupling of LND-ATM.
!       
!
! Author: R. Jacob, M. Vertenstein
! Revision History:
! 30Mar06 - P. Worley - added optional arguments to MCT_Rearrange call
! 13Apr06 - M. Vertenstein - cleaned up interfaces 
!
!---------------------------------------------------------------------

  use shr_kind_mod      ,only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_sys_mod
  use shr_const_mod
  use shr_mct_mod, only: shr_mct_sMatPInitnc
  use mct_mod

  use seq_comm_mct, only : logunit, loglevel
  use seq_cdata_mod
  use seq_flds_indices
  use seq_infodata_mod
  use seq_map_mod
  use m_die

  implicit none
  save
  private  ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: map_lnd2atm_init_mct
  public :: map_atm2lnd_init_mct
  public :: map_lnd2atm_mct
  public :: map_atm2lnd_mct

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(mct_rearr), private :: Re_lnd2atm
  type(mct_rearr), private :: Re_atm2lnd
  type(mct_sMatp), private :: sMatp_Fa2l
  type(mct_sMatp), private :: sMatp_Sa2l
  type(mct_sMatp), private :: sMatp_Fl2a
  type(mct_sMatp), private :: sMatp_Sl2a

#ifdef CPP_VECTOR
    logical :: usevector = .true.
#else
    logical :: usevector = .false.
#endif

#ifdef SYSUNICOS
    logical :: usealltoall = .true.
#else
    logical :: usealltoall = .false.
#endif
  logical, private :: samegrid_mapa2l

!=======================================================================
   contains
!=======================================================================

  subroutine map_atm2lnd_init_mct( cdata_a, cdata_l)

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_a
    type(seq_cdata),intent(in) :: cdata_l
    ! 
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_a           ! atm gsMap
    type(mct_gsMap), pointer :: gsMap_l           ! lnd gsMap
    type(mct_gGrid), pointer :: dom_a             ! atm domain
    type(mct_gGrid), pointer :: dom_l             ! lnd domain
    integer                  :: mpicom            ! communicator spanning atm and lnd
    integer                  :: ka, km            ! indices
    integer                  :: lsize             ! size of attribute vector
    type(mct_aVect)          :: areasrc           ! atm areas from mapping file
    type(mct_aVect)          :: areadst           ! lnd areas set to atm areas
    character(*),parameter :: subName = '(map_atm2lnd_init_mct) '
    !--------------------------------------------------

    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a, dom=dom_a)
    call seq_cdata_setptrs(cdata_l, gsMap=gsMap_l, dom=dom_l)
    call seq_cdata_setptrs(cdata_a, mpicom=mpicom, infodata=infodata)

    call seq_infodata_GetData( infodata, samegrid_al=samegrid_mapa2l)

    lsize = mct_gsMap_lsize(gsMap_a, mpicom)
    call mct_aVect_init( areasrc, rList="aream", lsize=lsize )
       
    lsize = mct_gsMap_lsize(gsMap_l, mpicom)
    call mct_aVect_init( areadst, rList="aream", lsize=lsize )

    if (samegrid_mapa2l) then

       call mct_rearr_init(gsMap_a, gsMap_l, mpicom, Re_atm2lnd)

       ! --- copy atm area into land aream

       ka = mct_aVect_indexRA(areasrc   , "aream")
       km = mct_aVect_indexRa(dom_a%data, "aream" )
       areasrc%rAttr(ka,:) = dom_a%data%rAttr(km,:)
       call mct_rearr_rearrange(areasrc, areadst, Re_atm2lnd, VECTOR=usevector, &
          ALLTOALL=usealltoall)

    else

       call shr_mct_sMatPInitnc(sMatp_Fa2l,gsMap_a,gsMap_l,"seq_maps.rc", &
          "atm2lndFmapname:","atm2lndFmaptype:",mpicom, &
          areasrc=areasrc, areadst=areadst)

       call shr_mct_sMatPInitnc(sMatp_Sa2l,gsMap_a,gsMap_l,"seq_maps.rc", &
          "atm2lndSmapname:","atm2lndSmaptype:",mpicom)

    endif

    ka = mct_aVect_indexRA(areadst   ,"aream")
    km = mct_aVect_indexRA(dom_l%data,"aream")
    dom_l%data%rAttr(km,:) = areadst%rAttr(ka,:)

    call mct_aVect_clean(areasrc)      
    call mct_aVect_clean(areadst)      

  end subroutine map_atm2lnd_init_mct

!=======================================================================

  subroutine map_lnd2atm_init_mct( cdata_l, cdata_a )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata),intent(in) :: cdata_l
    type(seq_cdata),intent(in) :: cdata_a
    ! 
    ! Local variables
    !
    type(seq_infodata_type), pointer :: infodata
    type(mct_gsMap), pointer :: gsMap_a           ! atm gsMap
    type(mct_gsMap), pointer :: gsMap_l           ! lnd gsMap
    type(mct_gGrid), pointer :: dom_l             ! lnd domain
    integer                  :: kf,iam,ierr,lsize
    integer                  :: mpicom            ! communicator spanning atm and lnd
    character(*),parameter :: subName = '(map_lnd2atm_init_mct) '
    !--------------------------------------------------

    call seq_cdata_setptrs(cdata_l, gsMap=gsMap_l, dom=dom_l)
    call seq_cdata_setptrs(cdata_a, gsMap=gsMap_a)
    call seq_cdata_setptrs(cdata_a, mpicom=mpicom,infodata=infodata)
    call mpi_comm_rank(mpicom,iam,ierr)
    
    call seq_infodata_GetData( infodata, samegrid_al=samegrid_mapa2l)

    ! Initialize lnd -> atm mapping or rearranger

    if (samegrid_mapa2l) then

       call mct_rearr_init(gsMap_l, gsMap_a, mpicom, Re_lnd2atm)

    else

       call shr_mct_sMatPInitnc(sMatp_Fl2a, gsMap_l, gsMap_a, "seq_maps.rc", &
            "lnd2atmFmapname:", "lnd2atmFmaptype:", mpicom)

       call shr_mct_sMatPInitnc(sMatp_Sl2a, gsMap_l, gsMap_a, "seq_maps.rc", &
            "lnd2atmSmapname:", "lnd2atmSmaptype:", mpicom)

    endif

  end subroutine map_lnd2atm_init_mct

!=======================================================================

  subroutine map_atm2lnd_mct( cdata_a, av_a, cdata_l, av_l, fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata), intent(in)  :: cdata_a
    type(mct_aVect), intent(in)  :: av_a
    type(seq_cdata), intent(in)  :: cdata_l
    type(mct_aVect), intent(out) :: av_l
    character(len=*),intent(in), optional :: statelist
    character(len=*),intent(in), optional :: fluxlist
    !
    ! Local
    ! 
    integer                :: lsize
    type(mct_aVect)        :: av_l_f     ! temporary flux attribute vector
    type(mct_aVect)        :: av_l_s     ! temporary state attribute vector
    character(*),parameter :: subName = '(map_atm2lnd_mct) '
    !--------------------------------------------------

    if (samegrid_mapa2l) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_a, av_l, Re_atm2lnd, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
              call mct_rearr_rearrange_fldlist(av_a, av_l, Re_atm2lnd, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_a, av_l, Re_atm2lnd, VECTOR=usevector, ALLTOALL=usealltoall)
       end if

    else

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
            call seq_map_avNorm(av_a, av_l, sMatp_Fa2l, rList=fluxlist, donorm=.false.)
          end if
          if (present(statelist)) then
            call seq_map_avNorm(av_a, av_l, sMatp_Sa2l, rList=statelist, donorm=.false.)
          end if
       else
          !--- default is flux mapping
          call seq_map_avNorm(av_a, av_l, sMatp_Fa2l, donorm=.false.)
       endif
    endif
       
  end subroutine map_atm2lnd_mct

!=======================================================================

  subroutine map_lnd2atm_mct( cdata_l, av_l, cdata_a, av_a, &
                              fractions_l, fractions_a, &
                              fluxlist, statelist )

    !--------------------------------------------------
    ! 
    ! Arguments
    !
    type(seq_cdata) ,intent(in)  :: cdata_l
    type(mct_aVect) ,intent(in)  :: av_l
    type(seq_cdata) ,intent(in)  :: cdata_a
    type(mct_aVect) ,intent(out) :: av_a
    type(mct_aVect) ,intent(in), optional :: fractions_l
    type(mct_aVect) ,intent(in), optional :: fractions_a
    character(len=*),intent(in), optional :: fluxlist
    character(len=*),intent(in), optional :: statelist
    !
    ! Local
    !
    integer  :: i,j,kl,lsize,numats,ier
    type(mct_aVect)          :: av_a_f     ! temporary flux attribute vector
    type(mct_aVect)          :: av_a_s     ! temporary state attribute vector
    character(*),parameter :: subName = '(map_lnd2atm_mct) '
    !--------------------------------------------------
!----------------------------- zhh debug 2013.07.23 ------------------------
!    print*, '====== at map_lnd2atm_mct ======: samegrid_mapa2l =', samegrid_mapa2l
!----------------------------- zhh debug 2013.07.23 ------------------------

    if (samegrid_mapa2l) then

       if (present(fluxlist) .or. present(statelist)) then
          if (present(fluxlist)) then
             call mct_rearr_rearrange_fldlist(av_l, av_a, Re_lnd2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=fluxlist)
          endif
          if (present(statelist)) then
             call mct_rearr_rearrange_fldlist(av_l, av_a, Re_lnd2atm, VECTOR=usevector, &
                  ALLTOALL=usealltoall, fldlist=statelist)
          endif
       else
          call mct_rearr_rearrange(av_l, av_a, Re_lnd2atm, VECTOR=usevector, ALLTOALL=usealltoall)
       end if

    else

       ! Normalize input data with fraction of land on land grid

!----------------------------- zhh debug 2013.07.23 ------------------------
!    print*, ' fractions_l =', present(fractions_l)
!----------------------------- zhh debug 2013.07.23 ------------------------

       if (present(fractions_l)) then

          if (present(fluxlist) .or. present(statelist)) then
             if (present(fluxlist)) then
                call seq_map_avNorm(av_l, av_a, sMatp_Fl2a, fractions_l, 'lfrin', fractions_a, 'lfrin', rList=fluxlist)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_l, av_a, sMatp_Sl2a, fractions_l, 'lfrin', fractions_a, 'lfrin', rList=statelist)
             end if
          else
             call seq_map_avNorm(av_l, av_a, sMatp_Fl2a, fractions_l, 'lfrin', fractions_a, 'lfrin')
          endif

       else

!----------------------------- zhh debug 2013.07.23 ------------------------
!          print*, ' fluxlist =', present(fluxlist)
!          print*, ' statelist =', present(statelist)
!----------------------------- zhh debug 2013.07.23 ------------------------
          if (present(fluxlist) .or. present(statelist)) then
             if (present (fluxlist)) then
                call seq_map_avNorm(av_l, av_a, sMatp_Fl2a, rList=fluxlist, donorm=.false.)
             end if
             if (present(statelist)) then
                call seq_map_avNorm(av_l, av_a, sMatp_Sl2a, rList=statelist, donorm=.false.)
             end if
          else
             ! --- default is flux mapping
             call seq_map_avNorm(av_l, av_a, sMatp_Fl2a, donorm=.false.)
          endif

       endif

    endif

  end subroutine map_lnd2atm_mct

end module map_atmlnd_mct
