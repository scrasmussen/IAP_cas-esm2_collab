subroutine da_transfer_wrftoxb(xbx, grid, config_flags)

   !---------------------------------------------------------------------------
   ! Purpose: Transfers fields from WRF to first guess structure.
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !---------------------------------------------------------------------------

   implicit none
   
   type (xbx_type), intent(inout)     :: xbx        ! Header & non-gridded vars.

   type(domain), intent(inout)        :: grid
   type(grid_config_rec_type), intent(in) :: config_flags

   integer :: i, j, k, ij

   real    :: theta, tmpvar

   real, dimension(ims:ime,jms:jme) :: rgh_fac

   character(len=19) :: current_date

   real :: loc_psac_mean

   real, dimension(jds:jde) :: loc_latc_mean

   integer :: size2d

   real, dimension(kms:kme) :: DDT

   real   :: qvf1, cvpm, cpovcv, ppb, ttb, albn, aln, height, temp
   real, allocatable :: arrayglobal(:,:)
#ifdef A2C
   real   :: uu, vv
#endif

   !---------------------------------------------------------------------------
   ! Set xb array range indices for processor subdomain.
   !---------------------------------------------------------------------------

   if (trace_use) call da_trace_entry("da_transfer_wrftoxb")

   grid%xb % map  = grid%map_proj
   grid%xb % ds   = grid%dx

   grid%xb % mix = grid%xp % ide - grid%xp % ids + 1
   grid%xb % mjy = grid%xp % jde - grid%xp % jds + 1
   grid%xb % mkz = grid%xp % kde - grid%xp % kds + 1

   ! WHY?
   !rizvi   xbx%big_header%bhi(5,5) = grid%start_year
   !rizvi   xbx%big_header%bhi(6,5) = grid%start_month
   !rizvi   xbx%big_header%bhi(7,5) = grid%start_day
   !rizvi   xbx%big_header%bhi(8,5) = grid%start_hour

   !---------------------------------------------------------------------------
   ! WRF-specific fitelds:
   !---------------------------------------------------------------------------

   ptop = grid%p_top

   grid%xb%sigmaf(kte+1) = grid%znw(kte+1)

   grid%xb%znw(kte+1) = grid%znw(kte+1)
   grid%xb%znu(kte+1) = 0.0
 
   do k=kts,kte
      grid%xb%sigmah(k) = grid%znu(k)
      grid%xb%sigmaf(k) = grid%znw(k)

      grid%xb%znu(k) = grid%znu(k)
      grid%xb%znw(k) = grid%znw(k)
      grid%xb%dn(k)  = grid%dn(k)
      grid%xb%dnw(k) = grid%dnw(k)
   end do

   grid%xb % ptop = ptop
      
   !---------------------------------------------------------------------------
   ! Convert WRF fitelds to xb:
   !---------------------------------------------------------------------------

   if (print_detail_xb) then
      write(unit=stdout, fmt='(3a, i8)') &
         'file:', __FILE__, ', line:', __LINE__

      write(unit=stdout, fmt=*) 'its,ite=', its,ite
      write(unit=stdout, fmt=*) 'jts,jte=', jts,jte
      write(unit=stdout, fmt=*) 'kts,kte=', kts,kte

      write(unit=stdout, fmt='(/5a/)') &
         'lvl         dnw                dn             rdnw       rdn'

      do k=kts,kte+1
         write(unit=stdout, fmt='(i3,8f16.8)') k, &
            grid%dnw(k), grid%dn(k), grid%rdnw(k), grid%rdn(k)
      end do

      write(unit=stdout, fmt='(/5a/)') &
         'lvl         znu                 znw           rdnw       rdn'

      do k=kts,kte+1
         write(unit=stdout, fmt='(i3,8f16.8)') k, &
            grid%xb%sigmah(k), grid%xb%sigmaf(k), grid%rdnw(k), &
            grid%rdn(k)
      end do

      write(unit=stdout, fmt='(/5a/)') &
         'lvl         phb                 ph_2'

      do k=kts,kte
         write(unit=stdout, fmt='(i3,8e20.12)') k, &
               grid%phb(its,jts,k), grid%ph_2(its,jts,k)
      end do

      write(unit=stdout, fmt=*) 'simple variables:'

      if (jte == jde) then
         write(unit=stdout, fmt=*) ' '

         do k=kts+5,kte,10
            do i=its,ite,10
               write(unit=stdout, fmt=*) &
                    '  grid%v_2(', i, ',', jde+1, ',', k, ')=', &
                       grid%v_2(i, jde+1,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      if (ite == ide) then
         write(unit=stdout, fmt=*) ' '

         do k=kts+5,kte,10
            do j=jts,jte,10
               write(unit=stdout, fmt=*) &
                  '  grid%u_2(', ide+1, ',', j, ',', k, ')=', &
                  grid%u_2(ide+1,j,k)
            end do
            write(unit=stdout, fmt=*) ' '
         end do
      end if

      write(unit=stdout, fmt=*) 'simple variables:'

      write(unit=stdout,fmt=*) &
         '  grid%u_1(its,jts,kts)=', grid%u_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%v_1(its,jts,kts)=', grid%v_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%w_1(its,jts,kts)=', grid%w_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%t_1(its,jts,kts)=', grid%t_1(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%ph_1(its,jts,kts)=',grid%ph_1(its,jts,kts)


      write(unit=stdout,fmt=*) &
         '  grid%u_2(its,jte,kts)=', grid%u_2(its,jte,kts)
      write(unit=stdout,fmt=*) &
         '  grid%v_2(ite,jts,kts)=', grid%v_2(ite,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%w_2(its,jts,kts)=', grid%w_2(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%t_2(its,jts,kts)=', grid%t_2(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%ph_2(its,jts,kts)=',grid%ph_2(its,jts,kts)
      write(unit=stdout,fmt=*) &
         '  grid%phb(its,jts,kts)=', grid%phb(its,jts,kts)

      write(unit=stdout, fmt=*) &
         '  grid%sm31,grid%em31,grid%sm32,grid%em32, grid%sm33,grid%em33=', &
         grid%sm31,grid%em31,grid%sm32,grid%em32,grid%sm33,grid%em33

      write(unit=stdout, fmt=*) '  grid%p_top=', grid%p_top
      write(unit=stdout, fmt=*) '  grid%znu(kts)=', grid%znu(kts)
      write(unit=stdout, fmt=*) '  grid%mub(its,jts)=', grid%mub(its,jts)
      write(unit=stdout, fmt=*) '  grid%mu_2(its,jts)=', &
         grid%mu_2(its,jts)

      write(unit=stdout, fmt=*) '  hbot(its,jts)=', grid%hbot(its,jts)
      write(unit=stdout, fmt=*) '  htop(its,jts)=', grid%htop(its,jts)

      write(unit=stdout, fmt=*) '  grid%p_top=', grid%p_top
      write(unit=stdout, fmt=*) '  num_moist=', num_moist
      write(unit=stdout, fmt=*) '  P_QV=', P_QV

      write(unit=stdout, fmt=*) '  moist(its,jts,kts,2)=', &
         grid%moist(its,jts,kts,2)
      write(unit=stdout, fmt=*) ' '
   end if

   !---------------------------------------------------------------
   ! Need this to exchange values in the halo region.
   ! grid%xa%u and grid%xa%v are used as temporary arrays and so
   ! it is easy to use the existing exchange scheme.
   !
   ! Note, this is needed as u_2 and v_2 has no guarantee
   ! the most east column, and the most north row are
   ! properly initailized for each tile.
   !---------------------------------------------------------------
#ifdef A2C

   grid%xa%u(its:ite+1,jts:jte,kts:kte) = grid%u_2(its:ite+1,jts:jte,kts:kte)
   grid%xa%v(its:ite,jts:jte+1,kts:kte) = grid%v_2(its:ite,jts:jte+1,kts:kte)
   grid%xb%u(its:ite+1,jts:jte,kts:kte) = grid%u_2(its:ite+1,jts:jte,kts:kte)
   grid%xb%v(its:ite,jts:jte+1,kts:kte) = grid%v_2(its:ite,jts:jte+1,kts:kte)

!rizvi's update
  if ((fg_format==fg_format_wrf_arw_regional  .or. &
       fg_format==fg_format_wrf_arw_global  ) .and. ide == ipe ) then
     ipe = ipe + 1
     ide = ide + 1
  end if

  if ((fg_format==fg_format_wrf_arw_regional  .or. &
       fg_format==fg_format_wrf_arw_global  ) .and. jde == jpe ) then
     jpe = jpe + 1
     jde = jde + 1
  end if
!rizvi's update
#ifdef DM_PARALLEL
#include "HALO_XB_UV.inc"
#endif

!rizvi's update
  if ((fg_format==fg_format_wrf_arw_regional  .or. &
       fg_format==fg_format_wrf_arw_global  ) .and. ide == ipe ) then
     ipe = ipe - 1
     ide = ide - 1
  end if

  if ((fg_format==fg_format_wrf_arw_regional  .or. &
       fg_format==fg_format_wrf_arw_global  ) .and. jde == jpe ) then
     jpe = jpe - 1
     jde = jde - 1
  end if

#else

   do j=jts,jte
      do k=kts,kte
         do i=its,ite+1
            grid%xa%u(i,j,k) = grid%u_2(i,j,k)
         end do
      end do
   end do

   do j=jts,jte+1
      do k=kts,kte
         do i=its,ite
            grid%xa%v(i,j,k) = grid%v_2(i,j,k)
         end do
      end do
   end do

   ! Fill the halo region for u and v.

#ifdef DM_PARALLEL
#include "HALO_PSICHI_UV_ADJ.inc"
#endif

#endif

   if (print_detail_xb) then
      write(unit=stdout, fmt=*) &
         ' ids,ide,jds,jde,kds,kde=', ids,ide,jds,jde,kds,kde
      write(unit=stdout, fmt=*) &
         ' its,ite,jts,jte,kts,kte=', its,ite,jts,jte,kts,kte
      write(unit=stdout, fmt=*) &
          ' ims,ime,jms,jme,kms,kme=', ims,ime,jms,jme,kms,kme
         
      write(unit=stdout, fmt=*) &
         ' lbound(grid%xb%u)=',   lbound(grid%xb%u)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%xb%v)=',   lbound(grid%xb%v)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%u_2)=', lbound(grid%u_2)
      write(unit=stdout, fmt=*) &
         ' lbound(grid%v_2)=', lbound(grid%v_2)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%xb%u)=',   ubound(grid%xb%u)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%xb%v)=',   ubound(grid%xb%v)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%u_2)=', ubound(grid%u_2)
      write(unit=stdout, fmt=*) &
         ' ubound(grid%v_2)=', ubound(grid%v_2)
   end if

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j, k, cvpm, cpovcv, ppb, temp, ttb ) &
   !$OMP PRIVATE ( albn, qvf1, aln, theta )
   do ij = 1 , grid%num_tiles

   do j=grid%j_start(ij), grid%j_end(ij)
      k = kte+1

      do i=its,ite
         grid%p(i,j,k) = 0.0
         grid%xb%map_factor(i,j) = grid%msft(i,j)
         grid%xb%cori(i,j) = grid%f(i,j)
         grid%xb%tgrn(i,j) = grid%sst(i,j)
         if (grid%xb%tgrn(i,j) < 100.0) &
            grid%xb%tgrn(i,j) = grid%tmn(i,j)
         grid%xb%lat(i,j) = grid%xlat(i,j)
         grid%xb%lon(i,j) = grid%xlong(i,j)
         grid%xb%terr(i,j) = grid%ht(i,j)
         grid%xb%snow(i,j) = grid%snowc(i,j)
         grid%xb%lanu(i,j) = grid%lu_index(i,j)
         grid%xb%landmask(i,j) = grid%landmask(i,j)
         grid%xb%xland(i,j) = grid%xland(i,j)
         ! Z. Liu below are variables used by RTTOV
         grid%xb%tsk(i,j) = grid%tsk(i,j)
         grid%xb%smois(i,j) = grid%smois(i,j,1)
         grid%xb%tslb(i,j) = grid%tslb(i,j,1)
         grid%xb%xice(i,j) = grid%xice(i,j)
         grid%xb%ivgtyp(i,j) = grid%ivgtyp(i,j)
         grid%xb%isltyp(i,j) = grid%isltyp(i,j)
         grid%xb%vegfra(i,j) = grid%vegfra(i,j)
         grid%xb%snowh(i,j) = grid%snowh(i,j)*1000.0 ! meter to mm    
      end do

      ! WHY?
      ! Adapted the code from "real.init.code" by Y.-R. Guo 05/13/2004:

      ! do i=its,ite
      !    k = kte
      !    qvf1 = 0.5*(grid%moist(i,j,k,P_QV)+grid%moist(i,j,k,P_QV))
      !    qvf2 = 1.0/(1.0+qvf1)
      !    qvf1 = qvf1*qvf2
      !    grid%xb%p(i,j,k) = -0.5*(grid%mu_2(i,j)+qvf1* &
      !       grid%mub(i,j))/grid%rdnw(k)/qvf2

      !    do k = kte-1,1,-1
      !       qvf1 = 0.5*(grid%moist(i,j,k,P_QV)+grid%moist(i,j,k+1,P_QV))
      !       qvf2 = 1.0/(1.0+qvf1)
      !       qvf1 = qvf1*qvf2
      !       grid%p(i,j,k) = grid%p(i,j,k+1) - &
      !          (grid%mu_2(i,j)+qvf1*grid%mub(i,j))/qvf2/rdn(k+1)
      !    end do
      ! end do

      ! Adapted the code from WRF module_big_step_utilitites_em.F ----
      !         subroutine calc_p_rho_phi      Y.-R. Guo (10/20/2004)

      ! NOTE: as of V3.1, P (pressure perturbation) and PB (base state pressure)
      ! are included in the wrfinput file. However, P and PB are still
      ! re-calculated here.

      cvpm =  - (1.0 - gas_constant/cp)
      cpovcv = cp / (cp - gas_constant)

      do k=kts,kte
         do i=its,ite
            ! The base specific volume (from real.init.code)
            ppb  = grid%znu(k) * grid%mub(i,j) + ptop
            grid%pb(i,j,k) = ppb
            temp = MAX ( iso_temp, base_temp + base_lapse*log(ppb/base_pres) )
            ttb  = temp * (base_pres/ppb)**kappa
            ! ttb  = (base_temp + base_lapse*log(ppb/base_pres)) * &
            !   (base_pres/ppb)**kappa
            albn = (gas_constant/base_pres) * ttb * (ppb/base_pres)**cvpm

            qvf1 = 1.0 + grid%moist(i,j,k,P_QV) / rd_over_rv
            aln  = -1.0 / (grid%mub(i,j)+grid%mu_2(i,j)) * &
               (albn*grid%mu_2(i,j) + grid%rdnw(k) * &
               (grid%ph_2(i,j,k+1) - grid%ph_2(i,j,k)))
            ! total pressure:
            grid%xb%p(i,j,k) = base_pres * &
               ((gas_constant*(t0+grid%t_2(i,j,k))*qvf1) / &
               (base_pres*(aln+albn)))**cpovcv
            ! total density
            grid%xb%rho(i,j,k)= 1.0 / (albn+aln)
            ! pressure purtubation:
            grid%p(i,j,k) = grid%xb%p(i,j,k) - ppb
         end do
      end do

      do k=kts,kte+1
         do i=its,ite
            grid%xb%hf(i,j,k) = (grid%phb(i,j,k)+grid%ph_2(i,j,k))/gravity
            grid%xa%w (i,j,k) = grid%w_2(i,j,k)
            grid%xb%w (i,j,k) = grid%w_2(i,j,k)
         end do
      end do

      do k=kts,kte
         do i=its,ite
            grid%xb%u(i,j,k) = 0.5*(grid%xa%u(i,j,k)+grid%xa%u(i+1,j,k))
            grid%xb%v(i,j,k) = 0.5*(grid%xa%v(i,j,k)+grid%xa%v(i,j+1,k))
            grid%xb%wh(i,j,k)= 0.5*(grid%xb%w(i,j,k)+grid%xb%w(i,j,k+1))
            grid%xb%h(i,j,k) = 0.5*(grid%xb%hf(i,j,k)+grid%xb%hf(i,j,k+1))

            grid%xb%q(i,j,k) = grid%moist(i,j,k,P_QV)

            theta = t0 + grid%t_2(i,j,k)
            grid%xb%t(i,j,k) = theta*(grid%xb%p(i,j,k)/base_pres)**kappa

            ! Convert to specific humidity from mixing ratio of water vapor:
            grid%xb%q(i,j,k)=grid%xb%q(i,j,k)/(1.0+grid%xb%q(i,j,k))
            if (num_pseudo == 0 .and. grid%xb%q(i,j,k) < 1.0e-9) &
               grid%xb%q(i,j,k) = 1.0e-9
   
            ! Background qrn needed for radar radial velocity assmiilation:

            if (size(grid%moist,dim=4) >= 4) then
               grid%xb%qcw(i,j,k) = grid%moist(i,j,k,p_qc)
               grid%xb%qrn(i,j,k) = grid%moist(i,j,k,p_qr)
               grid%xb%qt (i,j,k) = grid%xb%q(i,j,k) + grid%xb%qcw(i,j,k) + &
                  grid%xb%qrn(i,j,k)
            end if

            if (size(grid%moist,dim=4) >= 6) then
               grid%xb%qci(i,j,k) = grid%moist(i,j,k,p_qi)
               grid%xb%qsn(i,j,k) = grid%moist(i,j,k,p_qs)
            end if

            if (size(grid%moist,dim=4) >= 7) then
               grid%xb%qgr(i,j,k) = grid%moist(i,j,k,p_qg)
            end if

            if ( config_flags%mp_physics == 3 ) then   ! WSM3-class scheme
               if ( grid%xb%t(i,j,k) <= t_kelvin ) then
                  grid%xb%qci(i,j,k) = grid%xb%qcw(i,j,k)
                  grid%xb%qcw(i,j,k) = 0.0
                  grid%xb%qsn(i,j,k) = grid%xb%qrn(i,j,k)
                  grid%xb%qrn(i,j,k) = 0.0
               end if
            end if

         end do
      end do

      do i=its,ite
         grid%xb%psac(i,j) = grid%mub(i,j)+grid%mu_2(i,j)
         grid%xb%psfc(i,j) = grid%mub(i,j)+grid%p(i,j,kts)+grid%p_top

         if (grid%xb%tgrn(i,j) < 100.0) &    
            grid%xb%tgrn(i,j) = grid%xb%t(i,j,kts)+ &
            0.0065*(grid%xb%h(i,j,kts)-grid%xb%hf(i,j,kts))
      end do
   end do

   end do
   !$OMP END PARALLEL DO

   grid%xb%ztop = grid%xb%hf(its,jts,kte+1)

   if (print_detail_xb) then
      write(unit=stdout, fmt=*) ' '
      if (print_detail_xb) then
         write(unit=stdout, fmt='(/5a/)') &
            'lvl         h                 p                t'

         do k=kts,kte
            write(unit=stdout, fmt='(i3,8e20.12)') k, &
               grid%xb%h(its,jts,k), grid%xb%p(its,jts,k), grid%xb%t(its,jts,k)
         end do
      end if

      write(unit=stdout,fmt=*) ' '
      write(unit=stdout,fmt=*) 'grid%xb%u(its,jte,kte)=', grid%xb%u(its,jte,kte)
      write(unit=stdout,fmt=*) 'grid%xb%v(ite,jts,kte)=', grid%xb%v(ite,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%w(its,jts,kte)=', grid%xb%w(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%t(its,jts,kte)=', grid%xb%t(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%p(its,jts,kte)=', grid%xb%p(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%q(its,jts,kte)=', grid%xb%q(its,jts,kte)
      write(unit=stdout,fmt=*) 'grid%xb%h(its,jts,kte)=', grid%xb%h(its,jts,kte)
      write(unit=stdout,fmt=*) &
         'grid%xb%hf(its,jts,kte)=', grid%xb%hf(its,jts,kte)
      write(unit=stdout,fmt=*) &
         'grid%xb%map_factor(its,jts)=', grid%xb%map_factor(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%cori(its,jts)=', grid%xb%cori(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%tgrn(its,jts)=', grid%xb%tgrn(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lat(its,jts)=', grid%xb%lat(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lon(its,jts)=', grid%xb%lon(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%terr(its,jts)=', grid%xb%terr(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%snow(its,jts)=', grid%xb%snow(its,jts)
      write(unit=stdout,fmt=*) 'grid%xb%lanu(its,jts)=', grid%xb%lanu(its,jts)
      write(unit=stdout,fmt=*) &
         'grid%xb%landmask(its,jts)=', grid%xb%landmask(its,jts)
      write(unit=stdout,fmt=*) '(ite,jte)=', ite,jte                   
      write(unit=stdout,fmt=*) 'grid%xb%lat(ite,jte)=', grid%xb%lat(ite,jte)
      write(unit=stdout,fmt=*) 'grid%xb%lon(ite,jte)=', grid%xb%lon(ite,jte)
      write(unit=stdout,fmt=*) ' '
   end if

   !---------------------------------------------------------------------------
   ! [3.0] Calculate vertical inner product for use in vertical transform:
   !---------------------------------------------------------------------------
      
   if (vertical_ip == vertical_ip_sqrt_delta_p) then
      ! Vertical inner product is sqrt(Delta p):
      do k=kts,kte
         grid%xb % vertical_inner_product(its:ite,jts:jte,k) = &
            sqrt(grid%xb % psac(its:ite,jts:jte) * grid%xb%sigmah(k))
      end do 
   else if (vertical_ip == vertical_ip_delta_p) then

      ! Vertical inner product is Delta p:
      do k=1,grid%xb%mkz
         grid % xb % vertical_inner_product(its:ite,jts:jte,k) = &
         grid % xb % psac(its:ite,jts:jte) * grid%xb%sigmah(k)
      end do
   end if

   !---------------------------------------------------------------------------
   ! Roughness
   !---------------------------------------------------------------------------

   current_date = 'yyyy-mm-dd_hh:mm:ss'

   write(current_date(1:19), fmt='(i4.4, 5(a1, i2.2))') &
      grid%start_year, '-', &
      grid%start_month, '-', &
      grid%start_day, '_', &
      grid%start_hour, ':', &
      grid%start_minute, ':', &
      grid%start_second

   xbx % mminlu = 'USGS'

   call da_roughness_from_lanu(19, xbx % mminlu, current_date, &
      grid%xb % lanu, grid%xb % rough)

   !---------------------------------------------------------------------------
   ! Calculate 1/grid box areas:
   !---------------------------------------------------------------------------

   if (print_detail_xb) then
      write(unit=stdout, fmt='(/a, e24.16)') &
         'grid%xb % ds=', grid%xb % ds

      write(unit=stdout, fmt='(a, e24.16/)') &
           'grid%xb % map_factor(its,jts)=', grid%xb % map_factor(its,jts)
   end if

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, i, j, tmpvar, height, message )
   do ij = 1 , grid%num_tiles

   do j=grid%j_start(ij),grid%j_end(ij)
      do i=its,ite
         if (grid%xb%ztop < grid%xb%hf(i,j,kte+1)) &
             grid%xb%ztop = grid%xb%hf(i,j,kte+1)

         tmpvar = grid%xb%ds / grid%xb%map_factor(i,j)

         grid%xb % grid_box_area(i,j) = tmpvar*tmpvar

         ! Calculate surface variable(wind, moisture, temperature)
         ! sfc variables: 10-m wind, and 2-m T, Q, at cross points

         height = grid%xb%h(i,j,kts) - grid%xb%terr(i,j)
         if (height <= 0.0) then
            message(1) = "Negative height found"
            write (unit=message(2),FMT='(2I6,A,F10.2,A,F10.2)') &
               i,j,' ht = ',grid%xb%h(i,j,kts) ,' terr =  ',grid%xb%terr(i,j)
            call da_error(__FILE__,__LINE__, message(1:2))
         end if

#ifdef A2C
         uu = 0.5*(grid%xb%u(i,j,kts)+grid%xb%u(i+1,j,kts) )
         vv = 0.5*(grid%xb%v(i,j,kts)+grid%xb%v(i,j+1,kts) )
#endif
         call da_sfc_wtq(grid%xb%psfc(i,j), grid%xb%tgrn(i,j), &
#ifdef A2C
            grid%xb%p(i,j,kts), grid%xb%t(i,j,kts), grid%xb%q(i,j,kts),uu,vv, &
#else
            grid%xb%p(i,j,kts), grid%xb%t(i,j,kts), grid%xb%q(i,j,kts), &
            grid%xb%u(i,j,kts), grid%xb%v(i,j,kts), &
#endif
            grid%xb%p(i,j,kts+1), grid%xb%t(i,j,kts+1), grid%xb%q(i,j,kts+1), &
            height,  grid%xb%rough(i,j),grid%xb%xland(i,j), &
            grid%xb%u10(i,j), grid%xb%v10(i,j), grid%xb%t2(i,j), &
            grid%xb%q2(i,j), grid%xb%regime(i,j))
      end do
   end do

   end do
   !$OMP END PARALLEL DO

   !---------------------------------------------------------------------------
   ! Calculate saturation vapour pressure and relative humidity:
   !---------------------------------------------------------------------------

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, k, j, i )
   do ij = 1 , grid%num_tiles
      do k=kts,kte
         do j=grid%j_start(ij),grid%j_end(ij)
            do i=its,ite
               call da_tpq_to_rh(grid%xb % t(i,j,k), grid%xb % p(i,j,k), &
                  grid%xb % q(i,j,k), grid%xb %es(i,j,k), grid%xb %qs(i,j,k), &
                  grid%xb %rh(i,j,k))
            end do
         end do
      end do
   end do
   !$OMP END PARALLEL DO 

   !---------------------------------------------------------------------------
   ! Calculate dew point temperature:
   !---------------------------------------------------------------------------

   call da_trh_to_td (grid)

   if (print_detail_xb) then
      i=its; j=jts; k=kts

      write(unit=stdout, fmt=*) 'i,j,k=', i,j,k
      write(unit=stdout, fmt=*) 'grid%xb % td(i,j,k)=', grid%xb % td(i,j,k)
      write(unit=stdout, fmt=*) 'grid%xb % es(i,j,k)=', grid%xb % es(i,j,k)
      write(unit=stdout, fmt=*) 'grid%xb % rh(i,j,k)=', grid%xb % rh(i,j,k)
      write(unit=stdout, fmt=*) 'grid%xb % qs(i,j,k)=', grid%xb % qs(i,j,k)
      write(unit=stdout, fmt=*) ' '
   end if

   !---------------------------------------------------------------------------
   ! Sea level pressure and total precipitable water
   !---------------------------------------------------------------------------

   call da_wrf_tpq_2_slp (grid)

   ! WHY?
   ! do j = jts,jte
   !    do i = its,ite
   !       call da_tpq_to_slp(grid%xb%t(i,j,:), grid%xb%q(i,j,:), &
   !          grid%xb%p(i,j,:), grid%xb%terr(i,j), &
   !          grid%xb%psfc(i,j), grid%xb%slp(i,j), grid%xp)
   !    end do
   ! end do

   call da_integrat_dz(grid)

   !---------------------------------------------------------------------------
   ! Surface wind speed
   !---------------------------------------------------------------------------

   tmpvar = log(10.0/0.0001)

   !$OMP PARALLEL DO &
#ifndef A2C
   !$OMP PRIVATE (ij, i, j, height)
#else
   !$OMP PRIVATE (ij, i, j, height, uu, vv)
#endif
   do ij = 1, grid%num_tiles

   do j=grid%j_start(ij), grid%j_end(ij)
      do i=its,ite
         height = grid%xb%h(i,j,kts) - grid%xb%terr(i,j)
         rgh_fac(i,j) = 1.0/log(height/0.0001)
#ifndef A2C
         grid%xb%speed(i,j) = sqrt(grid%xb%u(i,j,kts)*grid%xb%u(i,j,kts) &
                         + grid%xb%v(i,j,kts)*grid%xb%v(i,j,kts) + 1.0e-6) &
                    *tmpvar*rgh_fac(i,j)
#else
         uu = 0.5*(grid%xb%u(i,j,kts)+grid%xb%u(i+1,j,kts) )
         vv = 0.5*(grid%xb%v(i,j,kts)+grid%xb%v(i,j+1,kts) )
         grid%xb%speed(i,j) = sqrt(uu*uu + vv*vv + 1.0e-6)*tmpvar*rgh_fac(i,j)
#endif
      end do
   end do

   end do
   !$OMP END PARALLEL DO

   !---------------------------------------------------------------------------
   ! Brightness temperature SH Chen
   !---------------------------------------------------------------------------

   if (use_ssmitbobs)   &
      call da_transform_xtotb(grid)

   !---------------------------------------------------------------------------
   ! GPS Refractivity linked by Y.-R. Guo 05/28/2004
   !---------------------------------------------------------------------------

   call da_transform_xtogpsref(grid)

   !---------------------------------------------------------------------------
   ! Ground-based GPS ZTD must follow the GPS Refractivity calculation.
   !---------------------------------------------------------------------------

   ! WHY? For certain computation method, not current one.
   if (use_gpsztdobs) then
      call da_transform_xtoztd(grid)
      if (print_detail_xb) then
        i=its; j=jts
        write(unit=stdout, fmt=*) 'grid%xb % tpw(i,j)=', grid%xb % tpw(i,j)
        write(unit=stdout, fmt=*) 'grid%xb % ztd(i,j)=', grid%xb % ztd(i,j)
        write(unit=stdout, fmt=*) ' '
      end if
   end if

   !---------------------------------------------------------------------------
   ! Calculate means for later use in setting up background errors.
   !---------------------------------------------------------------------------

   ! WHY?
   ! if (.not. associated(xbx % latc_mean)) then
   allocate (xbx % latc_mean(jds:jde))
   if (trace_use) call da_trace("da_transfer_wrftoxb",&
      message="allocated xbx%latc_mean")
   ! end if

   size2d = (ide-ids+1)*(jde-jds+1)

   tmpvar = 1.0/real(size2d)

   ! Bitwitse-exact reduction preserves operation order of serial code for
   ! testing, at the cost of much-increased run-time.  Turn it off when not
   ! testing.  Thits will always be .false. for a serial or 1-process MPI run.
   if (test_dm_exact) then
      allocate(arrayglobal(ids:ide, jds:jde))
      call da_patch_to_global(grid,grid%xb%psac, arrayglobal)
      loc_psac_mean = tmpvar*sum(arrayglobal(ids:ide,jds:jde))
      deallocate(arrayglobal)
   else
      loc_psac_mean = tmpvar*sum(grid%xb % psac(its:ite,jts:jte))
   end if

   tmpvar = 1.0/real(ide-ids+1)

   loc_latc_mean = 0.0

   if (test_dm_exact) then
      allocate(arrayglobal(ids:ide, jds:jde))
      call da_patch_to_global(grid,grid%xb%lat, arrayglobal)
      do j=jds,jde
         loc_latc_mean(j) = tmpvar*sum(arrayglobal(ids:ide, j))
      end do
      deallocate(arrayglobal)
   else
      do j=jts,jte
         loc_latc_mean(j) = tmpvar*sum(grid%xb % lat(its:ite, j))
      end do
   end if

   if (test_dm_exact) then
      ! Broadcast result from monitor to other tasks.
      call wrf_dm_bcast_real(loc_psac_mean, 1)
      xbx % psac_mean = loc_psac_mean
      ! Broadcast result from monitor to other tasks.
      call wrf_dm_bcast_real(loc_latc_mean, (jde-jds+1))
      xbx % latc_mean = loc_latc_mean
   else
      xbx % psac_mean = wrf_dm_sum_real(loc_psac_mean)
      call wrf_dm_sum_reals(loc_latc_mean, xbx % latc_mean)
   end if

   if (print_detail_xb) then
      ! write(unit=stdout, fmt=*) 'loc_psac_mean  =', loc_psac_mean
      write(unit=stdout, fmt=*) 'xbx % psac_mean=', xbx % psac_mean

      ! write(unit=stdout, fmt=*) 'loc_latc_mean  =', loc_latc_mean(jts)
      write(unit=stdout, fmt=*) 'xbx % latc_mean=', xbx % latc_mean(jts)
   end if


   ! Fill the halo region for xb        

#ifdef DM_PARALLEL
#include "HALO_XB.inc"
#endif

   ! Calculate time step from one dimensional cloud model parameterization

   if (dt_cloud_model) then
      do j = jts, jte
         do i = its, ite
            call da_cloud_model (grid%xb%t(I,J,:),  grid%xb%p(I,J,:), &
               grid%xb%q(I,J,:), grid%xb%qcw(I,J,:), grid%xb%qrn(I,J,:), &
               grid%xb%h(I,J,:), grid%xb%hf(I,J,:), ddt, kts, kte)

            do k = kts, kte
               grid%xb%delt(i,j,k) = DDT(k)
            end do
         end do
      end do
   end if

   if (trace_use) call da_trace_exit("da_transfer_wrftoxb")

end subroutine da_transfer_wrftoxb

