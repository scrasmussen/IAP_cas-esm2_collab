#include <define.h>

module colm_cplMod

#ifdef COUP_CSM

   use precision
   use colm_varctl

   implicit none

   real(r8), pointer :: lnd_taux(:)
   real(r8), pointer :: lnd_tauy(:)
   real(r8), pointer :: lnd_shflx(:)
   real(r8), pointer :: lnd_lhflx(:)
   real(r8), pointer :: lnd_qflx(:)
   real(r8), pointer :: lnd_swabs(:)
   real(r8), pointer :: lnd_lwup(:)
   real(r8), pointer :: lnd_scv(:)
   real(r8), pointer :: lnd_avsdr(:)
   real(r8), pointer :: lnd_avsdf(:)
   real(r8), pointer :: lnd_anidr(:)
   real(r8), pointer :: lnd_anidf(:)
   real(r8), pointer :: lnd_trad(:)
   real(r8), pointer :: lnd_tref(:)
   real(r8), pointer :: lnd_qref(:)
   real(r8), pointer :: lnd_u10m(:)
   real(r8), pointer :: lnd_v10m(:)
   real(r8), pointer :: lnd_sublim(:)
   real(r8), pointer :: lnd_nee(:)

   type forc_grid_type
      real(r8), pointer :: pco2m(:)
      real(r8), pointer :: po2m(:)
      real(r8), pointer :: us(:)
      real(r8), pointer :: vs(:)
      real(r8), pointer :: tm(:)
      real(r8), pointer :: qm(:)
      real(r8), pointer :: prc(:)
      real(r8), pointer :: prl(:)
      real(r8), pointer :: rain(:)
      real(r8), pointer :: snow(:)
      real(r8), pointer :: pbot(:)
      real(r8), pointer :: psrf(:)
      real(r8), pointer :: sols(:)
      real(r8), pointer :: soll(:)
      real(r8), pointer :: solsd(:)
      real(r8), pointer :: solld(:)
      real(r8), pointer :: frl(:)
      real(r8), pointer :: hu(:)
      real(r8), pointer :: ht(:)
      real(r8), pointer :: hq(:)
   end type forc_grid_type

   type(forc_grid_type) :: forcg ! forcing variables received from Coupler

   real(r8) :: eccen     !Earth's eccentricity factor (unitless) (typically 0 to 0.1)
   real(r8) :: obliq     !Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
   real(r8) :: mvelp     !Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)

  ! Orbital information after call to routine shr_orbit_params
  
   real(r8) :: obliqr    !Earth's obliquity in radians
   real(r8) :: lambm0    !Mean longitude (radians) of perihelion at the vernal equinox
   real(r8) :: mvelpp    !Earth's moving vernal equinox longitude

   interface colm_cpl_init
      module procedure colm_cpl_init
   end interface

   interface colm_cpl_l2a
      module procedure colm_cpl_l2a
   end interface

   interface colm_cpl_a2l
      module procedure colm_cpl_a2l
   end interface

   interface colm_cpl_exit
      module procedure colm_cpl_exit
   end interface

CONTAINS

   subroutine colm_cpl_init

      use colm_varMod, only: numgrid

      implicit none

! routine:

      allocate (lnd_taux   (numgrid))
      allocate (lnd_tauy   (numgrid))
      allocate (lnd_shflx  (numgrid))
      allocate (lnd_lhflx  (numgrid))
      allocate (lnd_qflx   (numgrid))
      allocate (lnd_swabs  (numgrid))
      allocate (lnd_lwup   (numgrid))
      allocate (lnd_scv    (numgrid))
      allocate (lnd_avsdr  (numgrid))
      allocate (lnd_avsdf  (numgrid))
      allocate (lnd_anidr  (numgrid))
      allocate (lnd_anidf  (numgrid))
      allocate (lnd_trad   (numgrid))
      allocate (lnd_tref   (numgrid))
      allocate (lnd_qref   (numgrid))
      allocate (lnd_u10m   (numgrid))
      allocate (lnd_v10m   (numgrid))
      allocate (lnd_sublim (numgrid))
      allocate (lnd_nee    (numgrid))

      allocate (forcg%pco2m(numgrid))
      allocate (forcg%po2m (numgrid))
      allocate (forcg%us   (numgrid))
      allocate (forcg%vs   (numgrid))
      allocate (forcg%tm   (numgrid))
      allocate (forcg%qm   (numgrid))
      allocate (forcg%prc  (numgrid))
      allocate (forcg%prl  (numgrid))
      allocate (forcg%rain (numgrid))
      allocate (forcg%snow (numgrid))
      allocate (forcg%pbot (numgrid))
      allocate (forcg%psrf (numgrid))
      allocate (forcg%sols (numgrid))
      allocate (forcg%soll (numgrid))
      allocate (forcg%solsd(numgrid))
      allocate (forcg%solld(numgrid))
      allocate (forcg%frl  (numgrid))
      allocate (forcg%hu   (numgrid))
      allocate (forcg%ht   (numgrid))
      allocate (forcg%hq   (numgrid))

   end subroutine colm_cpl_init

   subroutine colm_cpl_l2a

      use colm_varMod
      use timemgr, only: get_curr_year
      use landuse, only: landuse_cflux
      use spmd_decomp, only: gxmap, gymap

      implicit none

! local variables:

      real(r8) lucflux(lon_points,lat_points)
      integer curr_yr, g, i, j

! routine:

      curr_yr = get_curr_year()

      call landuse_cflux(curr_yr,lucflux)

    ! from g(C)/m2/s to kg(CO2)/m2/s
      lucflux(:,:) = lucflux(:,:)/12.0_r8*44.0_r8*1.0e-3

      do g = 1,numgrid
         lnd_taux(g)  = fldv%taux(g)
         lnd_tauy(g)  = fldv%tauy(g)
         lnd_shflx(g) = fldv%fsena(g)
         lnd_lhflx(g) = fldv%lfevpa(g)
         lnd_qflx(g)  = fldv%fevpa(g)
         lnd_swabs(g) = fldv%sabvsun(g) + fldv%sabvsha(g) + fldv%sabg(g)
         lnd_lwup(g)  = fldv%olrg(g)
         lnd_avsdr(g) = fldv%avsdr(g)
         lnd_avsdf(g) = fldv%avsdf(g)
         lnd_anidr(g) = fldv%anidr(g)
         lnd_anidf(g) = fldv%anidf(g)
         lnd_trad(g)  = fldv%trad(g)
         lnd_tref(g)  = fldv%tref(g)
         lnd_qref(g)  = fldv%qref(g)
         lnd_scv(g)   = fldv%scv(g)/1000._r8    ! mm(H2O) => m(H2O)
         lnd_u10m(g)  = fldv%u10m(g)
         lnd_v10m(g)  = fldv%v10m(g)
         lnd_sublim(g)= fldv%qsubl(g)

       ! nep=assim-respc-fmicr
       ! assim:photosynthesis rate[mol m-2 s-1], respc:plant respiration[mol m-2 s-1], fmicr:soil respiration[umol m-2 s-1]
         lnd_nee(g)   = fldv%assim(g)-fldv%respc(g)-fldv%fmicr(g)

         if(lnep_adjust) then
       !when dmcf=28.5
       !*28.5*0.5: convert molCO2/m2/s to gC/m2/s; 
       !/12: convert gC/m2/s to molC/m2/s;    zjw20180720
           lnd_nee(g) = lnd_nee(g)*28.5*0.5/12.0  
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           lnd_nee(g) = lnd_nee(g)-nep_residual(g)
         end if

         if(curr_yr.ge.lnd_cflux_year) then
          ! convert from molC/m2/s to kgCO2/m2/s
            lnd_nee(g) = lnd_nee(g)*44.0_r8*1.0e-3_r8

          ! add land use carbon emission
            i = gxmap(g)
            j = gymap(g)
            lnd_nee(g) = lnd_nee(g)-lucflux(i,j)
         else
            lnd_nee(g) = 0._r8
         end if
      end do

      if(lnep_adjust) lnep_adjust = .false.

#ifdef MYBUG
      write(6,*) 'lnd_taux' , minval(lnd_taux)    , maxval(lnd_taux)
      write(6,*) 'lnd_tauy' , minval(lnd_tauy)    , maxval(lnd_tauy)
      write(6,*) 'lnd_shflx', minval(lnd_shflx)   , maxval(lnd_shflx)
      write(6,*) 'lnd_lhflx', minval(lnd_lhflx)   , maxval(lnd_lhflx)
      write(6,*) 'lnd_qflx' , minval(lnd_qflx)    , maxval(lnd_qflx)
      write(6,*) 'lnd_swabs', minval(lnd_swabs)   , maxval(lnd_swabs)
      write(6,*) 'lnd_lwup' , minval(lnd_lwup)    , maxval(lnd_lwup)
      write(6,*) 'lnd_scv'  , minval(lnd_scv)     , maxval(lnd_scv)
      write(6,*) 'lnd_avsdr', minval(lnd_avsdr(:)), maxval(lnd_avsdr(:))
      write(6,*) 'lnd_avsdf', minval(lnd_avsdf(:)), maxval(lnd_avsdf(:))
      write(6,*) 'lnd_anidr', minval(lnd_anidr(:)), maxval(lnd_anidr(:))
      write(6,*) 'lnd_anidf', minval(lnd_anidf(:)), maxval(lnd_anidf(:))
      write(6,*) 'lnd_trad' , minval(lnd_trad)    , maxval(lnd_trad)
      write(6,*) 'lnd_tref' , minval(lnd_tref)    , maxval(lnd_tref)
      write(6,*) 'lnd_qref' , minval(lnd_qref)    , maxval(lnd_qref)
#endif

   end subroutine colm_cpl_l2a

   subroutine colm_cpl_a2l

      use colm_varMod, only: numcolumn, numgrid, forc
      use spmd_decomp, only: cgmap

      implicit none

      integer c, g

! routine:

      lsnowfrac = .true.

      do c = 1,numcolumn
         g = cgmap(c)
         forc%pco2m(c) = forcg%pco2m(g)
         forc%po2m (c) = forcg%po2m (g)
         forc%us   (c) = forcg%us   (g)
         forc%vs   (c) = forcg%vs   (g)
         forc%tm   (c) = forcg%tm   (g)
         forc%qm   (c) = forcg%qm   (g)
         forc%prc  (c) = forcg%prc  (g)
         forc%prl  (c) = forcg%prl  (g)
         forc%rain (c) = forcg%rain (g)
         forc%snow (c) = forcg%snow (g)
         forc%pbot (c) = forcg%pbot (g)
         forc%psrf (c) = forcg%psrf (g)
         forc%sols (c) = forcg%sols (g)
         forc%soll (c) = forcg%soll (g)
         forc%solsd(c) = forcg%solsd(g)
         forc%solld(c) = forcg%solld(g)
         forc%frl  (c) = forcg%frl  (g)
         forc%hu   (c) = forcg%hu   (g)
         forc%ht   (c) = forcg%ht   (g)
         forc%hq   (c) = forcg%hq   (g)
      end do

   end subroutine colm_cpl_a2l

   subroutine colm_cpl_exit

! routine:

      deallocate (forcg%pco2m)
      deallocate (forcg%po2m )
      deallocate (forcg%us   )
      deallocate (forcg%vs   )
      deallocate (forcg%tm   )
      deallocate (forcg%qm   )
      deallocate (forcg%prc  )
      deallocate (forcg%prl  )
      deallocate (forcg%rain )
      deallocate (forcg%snow )
      deallocate (forcg%pbot )
      deallocate (forcg%psrf )
      deallocate (forcg%sols )
      deallocate (forcg%soll )
      deallocate (forcg%solsd)
      deallocate (forcg%solld)
      deallocate (forcg%frl  )
      deallocate (forcg%hu   )
      deallocate (forcg%ht   )
      deallocate (forcg%hq   )

      deallocate (lnd_taux )
      deallocate (lnd_tauy )
      deallocate (lnd_shflx)
      deallocate (lnd_lhflx)
      deallocate (lnd_qflx )
      deallocate (lnd_swabs)
      deallocate (lnd_lwup )
      deallocate (lnd_scv  )
      deallocate (lnd_avsdr)
      deallocate (lnd_avsdf)
      deallocate (lnd_anidr)
      deallocate (lnd_anidf)
      deallocate (lnd_trad )
      deallocate (lnd_tref )
      deallocate (lnd_qref )
      deallocate (lnd_u10m )
      deallocate (lnd_v10m )
      deallocate (lnd_sublim )
      deallocate (lnd_nee  )

   end subroutine colm_cpl_exit

#endif

end module colm_cplMod
