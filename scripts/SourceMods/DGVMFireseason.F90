
  subroutine Fireseason (dtime,afirefrac1,nfireg1,forc_i,forc_rh,forc_wind,&
                         fm_col,fm_col1,agb_col,fa_col,parea)
                        
!-----------------------------------------------------------------------
! Calculate length of fire season in a year
! called every tstep.
! CALLED FROM: CLMMAIN
!-----------------------------------------------------------------------
    use precision
    implicit none

! INTENT IN VARIABLES:
    real(r8), INTENT(in) :: dtime        ! model step size
!    integer,  INTENT(in) :: ivt          ! vegetation type for this pft
!    real(r8), INTENT(in) :: litterag    ! above ground litter (gC/m2 veget'd area)
!    real(r8), INTENT(in) :: lm_ind
!    real(r8), INTENT(in) :: sm_ind
!    real(r8), INTENT(in) :: hm_ind
!    real(r8), INTENT(in) :: nind
!    real(r8), INTENT(in) :: wf
!    real(r8), INTENT(in) :: fpcgrid
 !   real(r8), INTENT(in) :: fpc_lf
    real(r8), INTENT(in) :: forc_i
    real(r8), INTENT(in) :: forc_rh
    real(r8), INTENT(in) :: forc_wind
    real(r8), INTENT(in) :: agb_col
    real(r8), INTENT(in) :: fm_col
    real(r8), INTENT(in) :: fm_col1
    real(r8), INTENT(in) :: fa_col
    real(r8), INTENT(in) :: parea
!    real(r8), INTENT(in) :: btran
    real(r8), INTENT(inout) :: afirefrac1
    real(r8), INTENT(inout) :: nfireg1
!    real(r8), INTENT(out) :: baf

! OTHER LOCAL VARIABLES:
!    integer  :: j   ! indices
    real(r8), parameter :: PI = 3.14159265358979323846
    
    real(r8) :: nfire         ! the number of fires in km2 per time step
    real(r8) :: fb            !available of fuel 
    real(r8) :: spread_m      ! moisture limitation for fire spread
    real(r8) :: ign           !presence of a source of  ignition per time step
    real(r8) :: Lb_lf         !length-to-breadth ratio   
    real(r8) :: fd_lf         !duration time
    real(r8) :: baf_lf        !burnt area for a PFT
    real(r8),parameter :: lfuel=110     ! lower threshold of fuel mass (gC/m2) for ignition
    real(r8),parameter  :: ufuel=1050     ! upper threshold of fuel mass(gC/m2) for ignition 
    real(r8),parameter  :: g0=0.05 
!**********************************************************************************************
    fb=0.
    ign=0.

 ! ********fire occurrence**********
       ign=forc_i/3600.0*dtime
       fb=max(0.0,min(1.0,(agb_col-lfuel)/(ufuel-lfuel)))
      nfire=ign*fb*fm_col*(1.0-max(0.0,min(1.0,(min(max(forc_rh,0.0),1.0)-0.3)/(0.7-0.3))))
 
      ! ********fire spread **********
      spread_m=fm_col1*(1.0-max(0.0,min(1.0,(min(max(forc_rh,0.0),1.0)-0.3)/(0.7-0.3))))
      fd_lf=24*3600
      Lb_lf=1.0+10.0*(1.0-EXP(-0.06*forc_wind))
      baf_lf=(g0*spread_m*fa_col*fd_lf/1000)**2*nfire*PI*Lb_lf
      afirefrac1=afirefrac1+baf_lf
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (afirefrac1.gt.1E2 )then
         write(6,*) 'imposible afirefrac1 (Fireseason) = ', afirefrac1
         write(6,*) 'imposible fm_col1 (Fireseason) = ', fm_col1
         write(6,*) 'imposible fm_col (Fireseason) = ',  fm_col
         write(6,*) 'imposible agb_col (Fireseason) = ', agb_col
         write(6,*) 'imposible fa_col (Fireseason) = ',  fa_col
         write(6,*) 'imposible forc_rh (Fireseason) = ', forc_rh
         write(6,*) 'imposible forc_i (Fireseason) = ',  forc_i
         write(6,*) 'imposible forc_wind (Fireseason)=', forc_wind
      end if    
  !burned area fraction from column to pft
 !     if(ivt(p) == 17)then
 !        baf(p)=0.0
  !    else
  !      baf(p)=baf_lf
  !     end if
      nfireg1=nfireg1+nfire*parea
 end subroutine Fireseason

