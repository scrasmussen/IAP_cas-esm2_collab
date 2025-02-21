module hycoef
!> \section arg_table_hycoef Argument Table
!! \htmlinclude hycoef.html
!!
    use shr_kind_mod, only: r8 => shr_kind_r8
    use spmd_utils,   only: masterproc
    use pmgrid,       only: plev, plevp
    use cam_logfile,  only: iulog
    use pio,          only: file_desc_t, var_desc_t, &
                            pio_inq_dimid, pio_inq_dimlen, pio_inq_varid, &
                            pio_double, pio_def_dim, pio_def_var, &
                            pio_put_var, pio_get_var
!wangty modify
#ifdef wrf
    use Dyn_const,    only: PMTOP, SIG, SIGL ! juanxiong he, 2010.08
#endif
    implicit none
    private

!-----------------------------------------------------------------------
!
! Purpose: Hybrid level definitions: p = a*p0 + b*ps
!          interfaces   p(k) = hyai(k)*ps0 + hybi(k)*ps
!          midpoints    p(k) = hyam(k)*ps0 + hybm(k)*ps
!
!-----------------------------------------------------------------------

    real(r8), public :: hyai(plevp)       ! ps0 component of hybrid coordinate - interfaces
    real(r8), public :: hyam(plev)        ! ps0 component of hybrid coordinate - midpoints
    real(r8), public :: hybi(plevp)       ! ps component of hybrid coordinate - interfaces
    real(r8), public :: hybm(plev)        ! ps component of hybrid coordinate - midpoints
    real(r8), public :: etamid(plev)      ! hybrid coordinate - midpoints

    real(r8), public :: hybd(plev)        ! difference  in b (hybi) across layers
    real(r8), public :: hypi(plevp)       ! reference pressures at interfaces
    real(r8), public :: hypm(plev)        ! reference pressures at midpoints
    real(r8), public :: hypd(plev)        ! reference pressure layer thickness

    real(r8), public, parameter :: ps0 = 1.0e5_r8    ! Base state surface pressure (pascals)
    real(r8), public, parameter :: psr = 1.0e5_r8    ! Reference surface pressure (pascals)

    integer, public :: nprlev       ! number of pure pressure levels at top

    public hycoef_init

    type(var_desc_t) :: hyam_desc, hyai_desc, hybm_desc, hybi_desc
    public init_restart_hycoef, write_restart_hycoef, read_restart_hycoef

contains

  subroutine hycoef_init

    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Defines the locations of model interfaces from input data in the
    ! hybrid coordinate scheme.  Actual pressure values of model level
    ! interfaces are determined elsewhere from the fields set here.
    !
    ! Method:
    ! the following fields are set:
    ! hyai     fraction of reference pressure used for interface pressures
    ! hyam     fraction of reference pressure used for midpoint pressures
    ! hybi     fraction of surface pressure used for interface pressures
    ! hybm     fraction of surface pressure used for midpoint pressures
    ! hybd     difference of hybi's
    ! hypi     reference state interface pressures
    ! hypm     reference state midpoint pressures
    ! hypd     reference state layer thicknesses
    ! hypdln   reference state layer thicknesses (log p)
    ! hyalph   distance from interface to level (used in integrals)
    ! prsfac   log pressure extrapolation factor (used to compute psl)
    !
    ! Author: B. Boville
    ! revised by juanxiong he, 2010.08
    !            ZhangHe, 2011-12-09
        !-----------------------------------------------------------------------
    !
    ! $Id$
    ! $Author$
    !
    !-----------------------------------------------------------------------

    use shr_kind_mod, only: r8 => shr_kind_r8
    use spmd_utils,   only: masterproc
    use pmgrid, only : plev, plevp
!    use Dyn_const,    only: PMTOP, SIG, SIGL   ! juanxiong he, 2010.08
!wangty modify
#ifdef wrf 
    use dycore,        only: dycore_is
#endif
    !---------------------------Local workspace-----------------------------
    integer  :: k        ! Level index
    real(r8) :: amean, bmean, atest, btest, eps
    !-----------------------------------------------------------------------

    eps    = 1.e-05_r8
    nprlev = 0

    ! Set layer locations
    !
    do k=1,plev
       !
       ! Interfaces. Set nprlev to the interface above, the first time a
       ! nonzero surface pressure contribution is found. "nprlev"
       ! identifies the lowest pure pressure interface.
       !
       if (nprlev==0 .and. hybi(k).ne.0.0_r8) nprlev = k - 1
    end do
    !
    ! Set nprlev if no nonzero b's have been found. All interfaces are
    ! pure pressure. A pure pressure model requires other changes as well.
    !
    if (nprlev==0) nprlev = plev + 2
    !
    ! Set delta sigma part of layer thickness and reference state midpoint
    ! pressures
    !

    do k=1,plev
!wangty modify 20180905
#ifdef wrf 
     if(dycore_is('IAP')) then
       hyai(k) = 0.0      !Zhang He, 2011-12-09
       hybi(k) = sig(k)   !juanxiong he
       hyam(k) = 0.0      !Zhang He, 2011-12-09
       hybm(k) = sigl(k)  !juanxiong he
     end if
#endif
       hybd(k) = hybi(k+1) - hybi(k)
       hypm(k) = hyam(k)*ps0 + hybm(k)*psr
       etamid(k) = hyam(k) + hybm(k)
    end do
    !
    ! Reference state interface pressures
    !
    do k=1,plevp
!wangty modify
#ifdef wrf 
     if(dycore_is('IAP')) then
       hyai(k) = 0.0      !Zhang He, 2011-12-09
       hybi(k) = sig(k)  !juanxiong he
     end if
#endif
       hypi(k) = hyai(k)*ps0 + hybi(k)*psr
    end do
    !
    ! Reference state layer thicknesses
    !
    do k=1,plev
       hypd(k) = hypi(k+1) - hypi(k)
    end do
    !
    ! Test that A's and B's at full levels are arithmetic means of A's and
    ! B's at interfaces
    !
    do k = 1,plev
       amean = ( hyai(k+1) + hyai(k) )*0.5_r8
       bmean = ( hybi(k+1) + hybi(k) )*0.5_r8
       if(amean == 0._r8 .and. hyam(k) == 0._r8) then
          atest = 0._r8
       else
          atest = abs( amean - hyam(k) )/ ( 0.5_r8*( abs(amean + hyam(k)) ) )
       endif
       if(bmean == 0._r8 .and. hybm(k) == 0._r8) then
          btest = 0._r8
       else
          btest = abs( bmean - hybm(k) )/ ( 0.5_r8*( abs(bmean + hybm(k)) ) )
       endif
       if (atest > eps) then
          if (masterproc) then
             write(iulog,9850)
             write(iulog,*)'k,atest,eps=',k,atest,eps
          end if
          !        call endrun
       endif

       if (btest > eps) then
          if (masterproc) then
             write(iulog,9850)
             write(iulog,*)'k,btest,eps=',k,btest,eps
          end if
          !        call endrun
       endif
    end do

    if (masterproc) then
       write(iulog,'(a)')'1 Layer Locations (*1000) '
       do k=1,plev
          write(iulog,9800)k,hyai(k),hybi(k),hyai(k)+hybi(k)
          write(iulog,9810) hyam(k), hybm(k), hyam(k)+hybm(k)
       end do

       write(iulog,9800)plevp,hyai(plevp),hybi(plevp),hyai(plevp)+hybi(plevp)
       write(iulog,9820)
       do k=1,plev
          write(iulog,9830) k, hypi(k)
          write(iulog,9840) hypm(k), hypd(k)
       end do
       write(iulog,9830) plevp, hypi(plevp)
    end if

9800 format( 1x, i3, 3p, 3(f10.4,10x) )
9810 format( 1x, 3x, 3p, 3(10x,f10.4) )
9820 format(1x,'reference pressures (Pa)')
9830 format(1x,i3,f15.4)
9840 format(1x,3x,15x,2f15.4)
9850 format('HYCOEF: A and/or B vertical level coefficients at full',/, &
         ' levels are not the arithmetic mean of half-level values')

  end subroutine hycoef_init

  subroutine init_restart_hycoef(File, vdimids)

    type(file_desc_t), intent(inout) :: File
    integer, intent(out) :: vdimids(:)

!
!   PIO traps errors internally, no need to check ierr
!
    integer :: ierr

    ierr = PIO_Def_Dim(File,'lev',plev,vdimids(1))
    ierr = PIO_Def_Dim(File,'ilev',plevp,vdimids(2))


    ierr = pio_def_var(File, 'hyai', pio_double, vdimids(2:2), hyai_desc)
    ierr = pio_def_var(File, 'hyam', pio_double, vdimids(1:1), hyam_desc)
    ierr = pio_def_var(File, 'hybi', pio_double, vdimids(2:2), hybi_desc)
    ierr = pio_def_var(File, 'hybm', pio_double, vdimids(1:1), hybm_desc)

    return

  end subroutine init_restart_hycoef

  subroutine write_restart_hycoef(file)

    type(file_desc_t), intent(inout) :: File


!
!   PIO traps errors internally, no need to check ierr
!
    integer :: ierr


    ierr = pio_put_var(File, hyai_desc, hyai)
    ierr = pio_put_var(File, hyam_desc, hyam)
    ierr = pio_put_var(File, hybi_desc, hybi)
    ierr = pio_put_var(File, hybm_desc, hybm)




  end subroutine write_restart_hycoef
!
! This code is used both for initial and restart reading.
!
  subroutine read_restart_hycoef(File)

    use abortutils,   only: endrun

    type(file_desc_t), intent(inout) :: File

    integer :: flev, filev, lev_dimid, ierr
!
!   PIO traps errors internally, no need to check ierr
!
    ierr = PIO_Inq_DimID(File, 'lev', lev_dimid)
    ierr = PIO_Inq_dimlen(File, lev_dimid, flev)
    if(plev/=flev) then
       write(iulog,*) 'Restart file lev does not match model. lev (file, model):',flev, plev
       call endrun()
    end if

    ierr = PIO_Inq_DimID(File, 'ilev', lev_dimid)
    ierr = PIO_Inq_dimlen(File, lev_dimid, filev)


    if(plevp /= filev) then
       write(iulog,*) 'Restart file ilev does not match model plevp (file, model):',filev, plevp
       call endrun()
    end if


    ierr = pio_inq_varid(File, 'hyai', hyai_desc)

    ierr = pio_inq_varid(File, 'hyam', hyam_desc)

    ierr = pio_inq_varid(File, 'hybi', hybi_desc)

    ierr = pio_inq_varid(File, 'hybm', hybm_desc)

    ierr = pio_get_var(File, hyai_desc, hyai)
    ierr = pio_get_var(File, hybi_desc, hybi)
    ierr = pio_get_var(File, hyam_desc, hyam)
    ierr = pio_get_var(File, hybm_desc, hybm)

#if ( defined OFFLINE_DYN )
!
! make sure top interface is non zero for fv dycore
!
   if (hyai(1) .eq. 0._r8) then
      if (hybm(1) .ne. 0.0_r8) then
         hyai(1) = hybm(1)*1.e-2_r8
      else if (hyam(1) .ne. 0.0_r8) then
         hyai(1) = hyam(1)*1.e-2_r8
      else
         call endrun('Not able to set hyai(1) to non-zero.')
      endif
   endif
#endif


  end subroutine read_restart_hycoef


end module hycoef
