#include <define.h>

module DustEmis

#ifdef DUST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   2014-05-20 Wu Chenglai
!   Dust Emission scheme from Shao's Scheme (Shao, 2004).
!   Related programs are copied from module_qf03.F in WRF_chem_dust.
!   The schemes needs the full-dispersed PSD and mini-dispersed PSD.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use precision
  use phycon_module, only : grav 
  implicit none

  contains

  subroutine DustEmission (npatch,num_filterp,filterp,wt_patch,wt_column,&
                    vegcov,lai,fsno,smois,ustar_in,rhoair,dustsource,isltyp,soil_top_cat,&
                    dustemis_bin_1,dustemis_bin_2,dustemis_bin_3,dustemis_bin_4,dustemis_total)
 
!---------------------Argument------------------------------------------

  integer, INTENT(in) ::&
        npatch             ,&! number of pfts in a column
        num_filterp       ,&! number of filtered patches in a column
        filterp(npatch)      ! array stored patch index in a column

  real(r8), INTENT(in) ::&
        wt_patch(npatch)   ,&! weight of pft relative to grid
        wt_column           ! weight of column relative to grid

  real(r8), INTENT(in) ::&
        vegcov            ,&! fraction of vegetation cover  old version: use fveg(npatch)                        
        lai(npatch)        ,&! leaf-area index
        fsno              ,&! fractional snow cover
        smois             ,&! soil moisture [m3/m3]
        ustar_in (npatch)  ,&! u* in similarity theory [m/s]
        rhoair              ! density air [kg/m3]

  integer,INTENT(in)  ::&
               dustsource ,&! index for potential dust source
                 isltyp     ! dominant soil type

  real(r8),INTENT(in) :: soil_top_cat(1:12) ! fraction for 12-categories soil

  real(r8),INTENT(out) ::&
        dustemis_bin_1(npatch) ,& ! dust emission flux for four size bins [kg/m2/s]
        dustemis_bin_2(npatch) ,&
        dustemis_bin_3(npatch) ,&
        dustemis_bin_4(npatch) ,&
        dustemis_total(npatch)    ! total dust emission flux [kg/m2/s] 
  
  ! 
  ! Wu Chenglai [2012-02-19] new varibles (the same to tlai)
  !   real(r8), pointer :: tvegcov(:)         
  ! 
  !  local variables
  !   real(r8) :: soil_tex(12)       

  !   integer i                 ! loop counter
  integer p,fp              ! pft index


  !  Special variables for dust emission scheme
  !   Variables:
  real(r8)          :: rho,  ustar, w, cf_in
  real(r8)          :: ust_min, rough_cor_in, smois_cor_in, smois_correc, g_in
  !Wu chenglai {2011/07/30} modify the value of calpha
  !       real(8) :: calpha
  integer :: i, j, j0, idst, kk, ij, imod, n   !not n, which has been assigned.
  real(r8) :: total, qtotal, ftotal, cf
  real(r8) :: ftotalb, ftotalp
  !
  !       real(8), parameter :: calpha = 5.d0, cbeta = 1.37d0
  !

  ! 2012-02-19 Wu Chenglai
 
  real(r8)                    :: tot_soilc       
  real(r8),dimension(16)      :: soilc       
  integer                     :: domsoilc,isl, cc
  REAL(r8)                    :: cont_clay
  character(4)                :: s_type
  !
  ! 2011-12-02 Wu Chenglai
  integer, parameter :: imax=113 ! No. of particle size intervals for psd
  integer, parameter :: jmax=3   ! No. of log-normal distributions for constructing psd
  real(r8), dimension(3, jmax) :: soilpsd  ! Coefs for sand minimally dispersed
 
  real(r8), dimension(0:imax) ::     d0
  real(r8), dimension(imax)   :: dd, psdm, dpsdm, ppsdm
  real(r8), dimension(imax)   ::     psdf, dpsdf, ppsdf
  real(r8), dimension(imax)   ::     psds, dpsds, ppsds
  !
  real(r8), dimension(imax)   :: beta_d, beta_s  ! beta1 and beta2 used by Shao et al. (1996) dust emission
  real(r8), dimension(imax)   :: ustart, q, ffq, fff
  real(r8), dimension(imax,imax) :: f
  !
  real(r8), parameter :: c_lambda=0.35
  real(r8) :: lambda, fc, phl, ccc, ghl
  ! 2011-12-02 Wu Chenglai
  real(r8), parameter :: dcut=10.e0_r8          ! dust cutoff particle size
  real(r8) :: cys, u0, al0, sx, ppr, rhos, smass, omega, rys
  real(r8) :: ddm, a1, a2, a3, a3b, a3p
  real(r8) :: zeta, sigma_m                 ! u*sqrt(rhos/p), bombardment coefficient
  !
  !real(8) :: r_c, h_c, qwhite, f_mb, f_hlys, pmass, ustart0, vhlys
  !external r_c, h_c, qwhite, f_mb, f_hlys, pmass, ustart0, vhlys
  real(r8) :: ustart0_out, qwhite_out, f_mb_out, f_hlys_out, pmass_out, vhlys_out

  integer, parameter :: nbins=4
  real(r8) :: sigma
  real(r8) :: fclay
  real(r8) :: dbin(nbins), fbin(nbins)
  integer :: ibin(nbins)
  ! 2011-12-02 Wu Chenglai
  data dbin /1.0_r8, 2.5_r8,5._r8,10._r8/    !size cut diameter (um)
  real(r8) :: rhop
  real     :: new_ftotal, new_test
  !kang [2009/01/07] 
  real(r8) :: cell_fbin(nbins), cell_ftotal


  integer, parameter :: mmax=4
  real(r8), dimension(3, mmax) :: csandm         ! Coefs for sand minimally dispersed
  data csandm /0._r8,     0._r8,     0._r8,     &
      &             0.0329_r8, 4.3733_r8, 0.8590_r8, &
      &             0.9671_r8, 5.7689_r8, 0.2526_r8, &
      &             0._r8,     0._r8,     0._r8      /

  real(r8), dimension(3, mmax) :: cloamm         ! Coefs for loam minimally dispersed
  data cloamm /0.1114_r8, 4.3565_r8, 0.4257_r8, &
      &             0.4554_r8, 5.1674_r8, 0.3824_r8, &
      &             0.4331_r8, 5.4092_r8, 1.0000_r8, &
      &             0._r8,     0._r8,     0._r8      /
 
  real(r8), dimension(3, mmax) :: csloamm        ! Coefs for sandy clay loam minimally dispersed, very dusty
  data csloamm/0.0722_r8, 2.2675_r8, 1.0000_r8, &
      &             0.6266_r8, 4.9654_r8, 0.3496_r8, &
      &             0.3012_r8, 5.5819_r8, 0.5893_r8, &
      &             0._r8,     0._r8,     0._r8      /
  real(r8), dimension(3, mmax) :: cclaym         ! Coefs for clay minimally dispersed
  data cclaym /0.3902_r8, 3.5542_r8, 1.0000_r8, &
      &             0.2813_r8, 4.2239_r8, 0.2507_r8, & 
      &             0.3286_r8, 5.1638_r8, 0.4632_r8, &
      &             0._r8,     0._r8,     0._r8      /
  real(r8), dimension(3, mmax) :: csandf         ! Coefs for sand fully dispersed
  data csandf /0._r8,     0._r8,     0._r8,     &
      &             0.0338_r8, 0.6931_r8, 1.0000_r8, &
      &             0.9662_r8, 5.6300_r8, 0.2542_r8, &
      &             0._r8,     0._r8,     0._r8      /
  real(r8), dimension(3, mmax) :: cloamf         ! Coefs for loam fully dispersed
  data cloamf /0.5844_r8, 4.6079_r8, 0.6141_r8, &
      &             0.3304_r8, 5.2050_r8, 0.2897_r8, &
      &             0.0522_r8, 7.0553_r8, 1.0000_r8, &
      &             0.0330_r8, 0.6931_r8, 1.0000_r8  /
  real(r8), dimension(3, mmax) :: csloamf         ! Coefs for sandy clay loam fully dispersed
  data csloamf/0.2344_r8, 1.8079_r8, 0.6141_r8, &
      &             0.3634_r8, 4.2050_r8, 0.2897_r8, &
      &             0.4022_r8, 5.6553_r8, 1.0000_r8, &
      &             0._r8,     0._r8,     0._r8      /
  real(r8), dimension(3, mmax) :: cclayf         ! Coefs for clay fully dispersed
  data cclayf /0.0872_r8, 0.6931_r8, 1.0000_r8, &
      &             0.4464_r8, 3.9323_r8, 0.9181_r8, &
      &             0.4665_r8, 5.4486_r8, 0.3916_r8, &
      &             0._r8,     0._r8,     0._r8      /
!+++++++++++++++++++++++++
!------------------------------------------------------------------------

  ! Loop through pfts

  ! initialize variables which get passed to the atmosphere
  dustemis_bin_1(:)=0.0
  dustemis_bin_2(:)=0.0
  dustemis_bin_3(:)=0.0
  dustemis_bin_4(:)=0.0
  dustemis_total(:)=0.0

  !Scheme select:
  imod=5

  if(dustsource.eq.1)then
    do fp = 1, num_filterp
       p = filterp(fp)
    

       !Parameters:
       rhop = 2560.e0_r8      ! particle density [kg/m3]
       rhos = 1000.e0_r8      ! bulk density of soil [kg/m3] ???
       phl = 30000._r8        ! plastic pressure [N/m2]
       cys = 0.00001_r8       ! cys : parameter

       
       ! only perform the following calculations if lnd_frc_mbl is non-zero 
       
       !if (dustsource(c).eq.0)cycle        ! Wu Chenglai 2012-01-12
       !if (lnd_frc_mbl(p) > 0.0_r8) then
       !+++++++++++++++++++++++++++++++++++++++++
       ! 2011-11-30 Wu Chenglai
       ! Begin of new lines:        
       ! Preparation:
       rho   = rhoair
       ustar = ustar_in(p)   
       sigma = rhop/rho    ! particle-air density ratio
       g_in = grav       ! change "g" to "grav"  
       !cf=tlai(p)
       !lambda=tlai(p)/4._r8   ! Test here !!!
       !cf  = tvegcov(p)/100._r8  !vegetation cover
       !cf=0.                      !currently just set cf=0. !for test version
       cf =vegcov/100._r8
       lambda=lai(p)/4.

       soilc(1:12)=soil_top_cat(1:12)  
       w  = smois        !vosoil moisture 
 
       ! slevis: adding liqfrac here, because related to effects from soil water

       ! liqfrac = max( 0.0_r8, min( 1.0_r8, h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_r8) ) )
 
       !-------------------
       ! frontal area index
       !-------------------
       !Wu[2010/08/05]

       ! lambda=max(0._r8,lambda)
       ! call r_c(lambda, rough_cor_in)

       !Soil moisture:

       ! isl = cc
       ! call h_c(w, isl, smois_correc)
                  
       !**************************************************
       ! 2012-02-19 Wu Chenglai
       !The codes below are copied from modulde_qf03.F

       !*************************************************
       !kang [2009/01/07] initialization
       cell_ftotal = 0._r8
       do n = 1, nbins
          cell_fbin(n) = 0._r8
       enddo

       !DO cc = 1, 16   ! soil category
       DO cc = 1, 12   ! new soil category

          if(soilc(cc).eq.0._r8) then
             go to 103
          endif

          if(cc.eq.1.or.cc.eq.2) then
             s_type = 'sand'
          else if (cc.eq.8.or.cc.eq.4.or.cc.eq.5.or.cc.eq.6.or.cc.eq.7.or.cc.eq.9) then
             s_type = 'loam'
          else if (cc.eq.3) then
             !Wu Chenglai, 2013-03-01
             !Wu Chenglai, 2013-03-19
             !To make the emission in East Asia(i.e, China and Mongolia) larger
             !s_type = 'sloa'
             s_type = 'loam'
          else if (cc.eq.10.or.cc.eq.11.or.cc.eq.12) then
             s_type = 'clay'
          else 
             go to 103
          endif

          if(s_type .eq. 'sand') then 
             call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax, csandm, jmax)      ! psd for minimally dispersed
             call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax, csandf, jmax)      ! psd for fully dispersed
          elseif(s_type .eq. 'sloa') then                                        ! sandy clay loam
             call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax, csloamm, jmax)      
             call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax, csloamf, jmax)      
          elseif (s_type .eq. 'loam') then 
             call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax, cloamm, jmax)      
             call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax, cloamf, jmax)      
          elseif (s_type .eq. 'clay') then 
             call psd_create(d0, dd, psdm, dpsdm, ppsdm, imax, cclaym, jmax)      
             call psd_create(d0, dd, psdf, dpsdf, ppsdf, imax, cclayf, jmax)      
          else
             !print *, "not yet programmed."
#ifdef MYBUG
             write(6,*)"call psd_create"
#endif
             write(6,*) "not yet programmed."
             stop
          endif

          !-------------------
          ! frontal area index
          !-------------------
          !Wu[2010/08/05]
          !if(cf>=(0.-1.e-9).and.cf<1.)then
          !cf=max(0.,cf)
          !lambda = - c_lambda*dlog( 1.d0 - cf )
          !r = r_c(lambda)

          lambda=max(0._r8,lambda)     
          call r_c(lambda, rough_cor_in)
          !else
          !lambda=999.
          !rough_cor_in=999.
          !endif
          !h = h_c(w)
          ! kang[2009/01/07] add soil category for moisture correction of Fecan
          !call h_c(w, smois_cor_in)
          ! kang[2009/02/25] Now, we're using new soil texture class. So we don't have to 
          ! match the class.
          ! Matching WRF_soil class with Shao's class (cc:WRF_soil class, isl:Shao's class)
          ! cc
          ! 1:sand, 2:loamy sand, 3:sandy loam, 4:silt loam, 5:silt, 6:loam, 7:sandy clay loam,
          ! 8:silty clay loam, 9:clay loam, 10:sandy clay, 11:silty clay, 12:clay
          ! isl
          ! 1:sand, 2:loamy sand, 3:sandy loam, 4:loam, 5:silt loam, 6:sandy clay loam, 
          ! 7:clay loam, 8:silty clay loam, 9:sandy clay, 10:silty clay, 11:clay, 12:heavy clay
          if(cc.eq.1.or.cc.eq.2.or.cc.eq.3.or.cc.eq.8) then
             isl = cc
          else if(cc.eq.4) then
             isl = 5
          else if(cc.eq.5) then
             isl = 5
          elseif(cc.eq.6) then
             isl = 4
          else if(cc.eq.7) then
             isl = 6
          else if(cc.eq.9) then
             isl = 7
          else if(cc.eq.10) then
             isl = 9
          else if(cc.eq.11) then
             isl = 10
          else if(cc.eq.12) then
             isl = 11
          else
             isl = 12  !!!test
          endif

          isl = cc
          call h_c(w, isl, smois_correc)

          !print*, 'lambda,r,h', lambda, r, h

          !----------------------------------------------
          ! for each particle size group, estimate ustart
          !----------------------------------------------
          !Threshold friction velocity:
          ust_min = 999.0_r8
          do i = 1, imax
             call ustart0(dd(i), sigma, g_in, rho, ustart0_out)
             ! write(iulog,*) i, sigma, g_in, rho, ustart0_out
             ustart(i) = ustart0_out
             ustart(i) = rough_cor_in*smois_correc*ustart(i)
             !write(iulog,*) 'smois cor', smois_correc
             !write(iulog,*) 'ust_t0', ustart0_out, 'correction', ustart(i)
             ust_min = dmin1(ust_min, ustart(i))

             !Streamwise sand flux:

             call qwhite(ustart(i), ustar, rho, g_in, qwhite_out)
             q(i) = qwhite_out
          enddo
          !write(iulog,*) 'ust_t_min ', ust_min
          !write(iulog,*) 'rough_cor: ', rough_cor_in
          !write(iulog,*) 'smois_cor: ', smois_correc
          !if (cc.eq.domsoilc) then
          ! smois_cor_in = smois_correc
          !endif

          !Vertical dust flux:
          if( ustar .le. ust_min ) then    ! no erosion goto 102
             q = 0.e0_r8
             ffq = 0.e0_r8
             qtotal = 0.e0_r8
             fff = 0.e0_r8
             ftotal = 0.e0_r8
             fbin = 0.e0_r8
             !goto 102
          else	
             ! print*, 'u*', ustar, 'u*t', ust_min
             ghl = dexp( -(ustar - ust_min)**3.e0_r8 )
             dpsds = ghl*dpsdm + (1._r8-ghl)*dpsdf
             psds  = ghl*psdm + (1._r8-ghl)*psdf
             ppsds = ghl*ppsdm + (1._r8-ghl)*ppsdf

             ffq = q*dpsds
             qtotal = sum(ffq)

             !--------------
             ! dust emission
             !--------------
             ! Before calculating dust emission, some parameters should be set.

             ! for Lu & Shao (1999) scheme, mass fraction of dust is needed.
             j = 0
1            j = j+1
             if( dd(j) .le. dcut ) then
                idst = j
                goto 1
             endif
             fc = 0._r8
             do j = 1, idst
                fc = fc + dpsds(j)              ! mass fraction of dust
             enddo

             ! for MB95 scheme, clay contents are needed.
             fclay=0._r8
             do i = 1, imax
                if(d0(i).le.2.1_r8) fclay=fclay+dpsds(i)
             enddo 

             !write(iulog,*) 'fclay',  fclay
             if(imod.eq.2.and.fclay.eq.0._r8) stop 'wrong clay content'

             ! size bin
             do n=1,nbins
                ibin(n)=0
                do i=imax,1,-1
                   if(d0(i).ge.dbin(n)) ibin(n)=i
                enddo
                if(ibin(n).eq.0) stop 'wrong dust classes'
                !write(iulog,*) n, ibin(n),dbin(n),d0(ibin(n))
             enddo

             !write(iulog,*) 'ustar : ', ustar
             !------------------------------------------------------------------------
             ! Shao et al. (1996) dust emission model, for dust smaller than 20 micron
             !------------------------------------------------------------------------
             IF(imod .eq. 1 ) THEN 
                !print*,'imod = 1, come here?'
                beta_d=0._r8
                beta_s=0._r8
                do j = 1, idst
                   beta_d(j) = 1.e-1_r8*dexp(-140.7_r8*dd(j)*1.e-3_r8+0.37_r8)     ! dd in [um] 
                enddo

                do i = idst+1, imax
                   beta_s(i)= 0.125e-4_r8*dlog(dd(i)*1.e-3_r8) + 0.328e-4_r8      ! dd in [um]
                   beta_s(i)= dmax1(beta_s(i),0.e0_r8)
                enddo
          
                f=0._r8
                ccc=2._r8*rhop*2.5_r8*g_in/(3._r8*rho)*0.01_r8
                do i = idst+1, imax
                do j = 1, idst
                   f(i,j) = ccc*beta_d(j)*beta_s(i)*ffq(i)/ustart(j)**2._r8
                   !if(f(i,j)>10000._r8)then
                   !write(iulog,*)"test_01",ccc,beta_d(j),beta_s(i),ffq(i),ustart(j)
                   !endif              
                enddo
                enddo
           
                !write(iulog,*)"ffq:", ffq
                !write(iulog,*)"beta_d: ",beta_d
                !write(iulog,*)"beta_s: ",beta_s
                !write(iulog,*)"ustart: ",ustart
        
                !write(iulog,*)"ccc:",ccc         
                !write(iulog,*)"rho",rho          
 
                ftotal = 0.0_r8
                !new_ftotal=0.
                !new_test=0.
                do j = 1, idst
                   fff(j) = 0._r8
                   do i = idst+1, imax
                      fff(j) = fff(j) + f(i,j)*dpsds(j)
                      !new_test=max(new_test,f(i,j))
                   enddo
                   ftotal = ftotal + fff(j)
                   !new_ftotal=new_ftotal+fff(j)
                enddo
                            
                !write(iulog,*) 'ftotal', ftotal
                !write(iulog,*) 'new_ftotal', new_ftotal
                !write(iulog,*) 'new_test', new_test

                do n=1,nbins
                   j0=1
                   if(n.gt.1) j0=ibin(n-1)+1
                   fbin(n)=0._r8
                   do j=j0,ibin(n)
                      fbin(n)=fbin(n)+fff(j)
                      !write(iulog,*)'fbin', n,fbin(n)
                   enddo
                enddo

                !print*, 'fbin', n, fbin(n)
                !enddo
                !print*, 'ftotal', ftotal

             !--------------------------------
             ! Shao (2004) dust emission model
             !--------------------------------
             ELSEIF (imod .eq. 5) THEN 
                !print*,'imod = 5, come here??'
                !Wu [2010/01/29] add for experiment
                !Note that different soil type should have different parameters(Shao,2004)
                !e.g.,Cy and p should be larger and smaller for loose sandy soils respectly.
                if(s_type .eq. 'sand') then
                   !phl=1.0e+3_r8
                   !cys=5.0e-5_r8
                   !!Wu[2013-03-15] for test
                   if(cc .eq. 1)then
                      phl=5.0e+2_r8
                      !Wu Chenglai 2013-03-21
                      !v15: cys=5.0e-5_r8 probably too small
                      !v16: 
                      cys=6.0e-5_r8
                   else
                      phl=1.0e+4_r8
                      cys=1.0e-6_r8
                   endif
                endif
                if(s_type .eq. 'sloa') then
                   !Wu Chenglai 2013-03-21
                   !v15: phl=1.0e+4_r8 probably too large
                   !v16:
                   phl=2.0e+4_r8
                   cys=3.0e-5_r8
                endif
                if(s_type .eq. 'loam') then
                   if(cc.eq.3)then
                      phl=3.0e+4_r8
                      cys=1.0e-6_r8
                   else
                      phl=3.0e+4_r8
                      cys=1.0e-6_r8
                   endif
                endif
                if(s_type .eq. 'clay') then
                   phl=5.0e+4_r8
                   !2013-03-14, Wu Chenglai
                   !It is possible that soil type of clay produces two much dust emission;
                   ! so we scale the cy according to the comparison of observed and offline simulation
                   ! for test, we scale it by 0.01.
                   !cys=1.0e-5_r8
                   cys=1.0e-7_r8
                endif

                zeta = ustar*dsqrt( rhos/phl )
                sigma_m = 12.e0_r8*zeta*zeta*(1.e0_r8+14.e0_r8*zeta)


                do i = idst+1, imax
                do j = 1, idst
                   rys = psdm(j)/psdf(j)
                   a1 = cys*dpsdf(j)*( (1._r8-ghl) + ghl*rys )
                   a2 = (1._r8+sigma_m)
                   a3 = ffq(i)*g_in/ustar**2._r8
                   f(i,j) = a1*a2*a3
                enddo
                enddo
          
                ftotal = 0.0_r8
                do j = 1, idst
                   fff(j) = 0._r8
                   do i = idst+1, imax
                      fff(j) = fff(j) + f(i,j)
                   enddo
                   ftotal = ftotal + fff(j)
                enddo
 
                do n=1,nbins
                   j0=1
                   if(n.gt.1) j0=ibin(n-1)+1
                   fbin(n)=0._r8
                   do j=j0,ibin(n)
                      fbin(n)=fbin(n)+fff(j) 
                   enddo
                enddo
                !do n=1, nbins
                !print*, 'fbin', n, fbin(n)
                !enddo
                !print*, 'ftotal', ftotal
         
             ELSE
                write(6,*)"Dust emission not calculate" 
                stop
                ! print*,"Dust emission not calculate" 
             ENDIF       ! imod
          ENDIF   !ustar vs ust_min

!102   continue
  

          do n = 1, nbins
             cell_fbin(n) = cell_fbin(n) + soilc(cc)*fbin(n)

             !cell_fbin(n) = cell_fbin(n) + (soilc(cc)/tot_soilc)*fbin(n)
             !print*, 'fbin', cc, cell_fbin(n)
          enddo

          cell_ftotal = cell_ftotal + ftotal*soilc(cc)

          !cell_ftotal = cell_ftotal + ftotal*soilc(cc)/tot_soilc
          !print*, 'ftotal', cc, cell_ftotal

103   continue  

       END do      !cc
 
       !if(cf>=(0.-1.e-9).and.cf<1.)then
       cf=max(0._r8,cf)
       cf=min(1._r8,cf)
       !do n=1,nbins
       !cell_fbin(n)=cell_fbin(n)*(1.-cf)
       !enddo
       !cell_ftotal=cell_ftotal*(1.-cf)
       !else
       ! do n=1,nbins
       !    cell_fbin(n)=0.
       ! enddo
       ! cell_ftotal=0.
       ! endif

       ! if(cell_ftotal>8.e-6)then
       ! do n=1,nbins
       !    cell_fbin(n)=cell_fbin(n)*8.e-6/cell_ftotal
       ! enddo 
       ! cell_ftotal=8.e-6
       ! endif

       ! Accumulate the dust emission to 4 bins:
       ! do n=1,nbins
       !    flx_mss_vrt_dst(p,n) = cell_fbin(n)* & 
       !                           lnd_frc_mbl(p) * liqfrac*(1._r8-cf)
       !    flx_mss_vrt_dst_tot(p) = flx_mss_vrt_dst_tot(p) + flx_mss_vrt_dst(p,n)
       ! enddo
       ! Wu Chenglai [2012-02-28]
       ! flx_mss_vrt_dst_1(p) =  flx_mss_vrt_dst(p,1)
       ! flx_mss_vrt_dst_2(p) =  flx_mss_vrt_dst(p,2)
       ! flx_mss_vrt_dst_3(p) =  flx_mss_vrt_dst(p,3)
       ! flx_mss_vrt_dst_4(p) =  flx_mss_vrt_dst(p,4)
       ! end of new varibles Wu Chenglai
       ! end if   ! lnd_frc_mbl > 0.0
       dustemis_bin_1(p)  = cell_fbin(1)*(1._r8-cf)*(1._r8-fsno)
       dustemis_bin_2(p)  = cell_fbin(2)*(1._r8-cf)*(1._r8-fsno)
       dustemis_bin_3(p)  = cell_fbin(3)*(1._r8-cf)*(1._r8-fsno)
       dustemis_bin_4(p)  = cell_fbin(4)*(1._r8-cf)*(1._r8-fsno)
       dustemis_total(p)  = cell_ftotal*(1._r8-cf)*(1._r8-fsno)

    end do  !fp
 
   ENDIF       ! dustsource.eq.0
 
  end subroutine DustEmission


!+++++++++++++++++++++++++++++++++++++++++++
! 2011-11-30, Wu Chenglai
! new subroutines
!*****************************************************************************
      subroutine psd_create(d, dm, psd, dpsd, ppsd, imax, cmtrix, jmax)
!
!----------------------------------------------------------------------------
! Yaping Shao, 13 June 2000
!
! - Generate particle size distribution density function 
!   (both minimally-dispersed and fully-dispersed as the 
!   sum of four log-normal distributions.
!
! d(0,imax):	output, particle size at 0, 1, 2, ..., imax points	[um]
! dm(imax):	output, particle size at middle of 0-1, 1-2, etc        [um]
! psd(imax):    output, particle size distribution density at dm	[um^-1]
! dpsd(imax):   output, Delta P for sections 0-1, 1-2, etc.  		[ ]
! ppsd(imax):   output, P for sections 0-1, 1-2, etc.  			[ ]
! imax:         input, length dm, psd, dpsdm, ppsd, etc.
! cmtrix: 	jmaxx coefficient matrix 
!               e.g.
!		w1   = cmtrix(1, 1): weight for first log-normal distribution
!               dln1 = cmtrix(2, 1): mean log-particle size of first log-normal distribution 	
!               sig1 = cmtrix(3, 1): sigma of log-particle size for first log-normal distribution
!		etc. 
!		careful with the dimension of dln and sig
!----------------------------------------------------------------------------
!
      integer :: i, j, imax, jmax	
      real(r8), dimension(3, jmax) :: cmtrix
      real(r8) :: d(0:imax), dm(imax)                 
      real(r8) :: psd(imax), dpsd(imax), ppsd(imax)   ! for p(d), Delta P(d) and P(d)
      real(r8) :: p, pp, w, dln, sig
      real(r8) :: cn, sum_percent
      real(r8), parameter :: eps=1.e-7_r8
      real(r8), parameter :: dref=1000.e0_r8
      real(r8) :: fu, fd, phi
!
      cn = 1.e0_r8/dsqrt(2.e0_r8*3.14159e0_r8)
!
! initialise psd, dpsd, ppsd
!
      psd = 0.e0_r8
      dpsd = 0.e0_r8
      ppsd = 0.e0_r8
!
! Estimate d using phi scale. phi varies between from 9 to -1
! with increment 0.1. Reference particle size d0 = 1000 um
!
! 2011-12-02 Wu Chenglai
!      fu = 9.e0_r8
!      fu = 12.e0_r8
!      fu = 11.e0_r8
      fu = 10.3e0_r8
      fd = -1.e0_r8
      do i = 0, imax
	phi = fu - i*(fu-fd)/imax
	d(i) = dref/2.e0_r8**phi
      enddo 	
!
! 2011-12-01 Wu Chenglai  (only for %sand+%silt+%clay >%50)
   sum_percent=sum(cmtrix(1,1:jmax))
   if(sum_percent.gt.0.5_r8)then
 
      do i = 1, imax
        dm(i) = dexp( (dlog(d(i))+dlog(d(i-1)) )/2.e0_r8 )

        pp = 0.e0_r8
	do j = 1, jmax
	  w = cmtrix(1, j)
          dln = cmtrix(2, j)
          sig = cmtrix(3, j)
	  if ( (w.gt.eps) .and. (sig.ne.0._r8) ) then
	    p = w*cn/sig*dexp( -(dlog(dm(i))-dln)**2._r8/(2._r8*sig**2._r8) )
          else
            p=0.e0_r8	
          endif
	  pp = pp + p
	enddo
!
	dpsd(i) = pp*( dlog(d(i)) - dlog(d(i-1)) )        ! Delta P over i
	if (i.eq.1) then 
          ppsd(i) = 0.e0_r8 + dpsd(i)                        ! P(d), with P(0) = 0
	else
          ppsd(i) = ppsd(i-1) + dpsd(i)      
	endif
        psd(i) = pp/dm(i)                                 ! p(d), particle size distribution density

      enddo
!
! Renormalisation, in case ppsd(imax) is not 1
!
      dpsd = dpsd/ppsd(imax)
      psd  = psd/ppsd(imax)
      ppsd = ppsd/ppsd(imax)
   else
     dpsd = 0._r8
     psd  = 0._r8
     ppsd = 0._r8
   endif  !sum_percent

!
     end subroutine psd_create

!*****************************************************************************
!      real(8) function ustart0(dum, sigma, g_in, rho)
      subroutine ustart0(dum, sigma, g_in, rho, ustart0_out)
!
! Y. Shao, 13 June 2000
!
! Calculate ustar0(d) using Shao and Lu (2000) for uncovered
! dry surface
!
! dum:    particle diameter			[um]
! ustar0: threshold friction velocity   	[m/s]
!
!      real(8) :: dm, dum
      real(r8), intent(in) :: dum, sigma, g_in, rho
      real(r8), intent(out) :: ustart0_out
      real(r8) :: dm
      real(r8), parameter :: gamma = 1.65e-4_r8      ! a constant
      real(r8), parameter :: f = 0.0123_r8

      dm = dum*1e-6_r8

!      ustart0 = f*(sigma*g_in*dm + gamma/(rho*dm) )
!      ustart0 = dsqrt( ustart0 )
      ustart0_out = f*(sigma*g_in*dm + gamma/(rho*dm) )
      ustart0_out = dsqrt( ustart0_out )

!      end function
      end subroutine ustart0
!*****************************************************************************
!      real(8) function qwhite(ust, ustar, rho, g_in)
      subroutine qwhite(ust, ustar, rho, g_in, qwhite_out)
!
! Yaping Shao 17-07-99!
!
! White (1979) Sand Flux Equation
! Q = c*rho*u_*^3 over g (1 - u_*t over u_*)(1 + u_*t^2/u_*^2)
! qwhite: Streamwise Sand Flux;       [kg m-1 s-1]
! c     : 2.6
! ust   : threhold friction velocity  [m/s]
! ustar : friction velocity           [m/s]
!C
!      real(8) :: c, ust, ustar
      real(r8) :: c
      real(r8), intent(in) :: ust, ustar, rho, g_in
      real(r8), intent(out) :: qwhite_out
      real(r8) :: a, b

!!Wu[2010/02/04]change c to 2.6
      !c = 1.0e0_r8
      !c  = 10.e0_r8  !test_01
      !c = 50.e0_r8   !test_02
       c = 85.e0_r8   
!!Wu[2013-03-15] change c to 10. for global model
!      c = 1.0e1_r8

!      c=2.6d0
      a = rho/g_in
!      IF (ustar.lt.ust) THEN 
!        qwhite = 0.d0
!      ELSE
!        b = ust/ustar 
!        qwhite = c*a*ustar**3.*(1.-b)*(1.+b*b)
!      ENDIF
      IF (ustar.lt.ust) THEN 
        qwhite_out = 0.e0_r8
      ELSE
        b = ust/ustar 
        qwhite_out = c*a*ustar**3._r8*(1._r8-b)*(1._r8+b*b)
      ENDIF

!      END function
      END subroutine qwhite
!
!*****************************************************************************
!      real(8) FUNCTION r_c (x)
      subroutine r_c (x,r)
!
!   Y. Shao 17-07-92
!   CORRECTION FUNCTION FOR UST(D) BASED ON Raupach et al. (1992)
!   x = frontal area index
!
!   R_C = (1 - sig m x)^{1/2} (1 + m beta x)^{1/2}
!         Note I deife R_C = u_{*tR}/u_{*tS}
!         While Raupach et al. defined 
!                      R_C = u_{*tS}/u_{*tR} and their R function is
!   R_C = (1 - sig m x)^{-1/2} (1 + m beta x)^{-1/2}
!                      
!   sig   : basal to frontal area; about 1
!   m     : parameter less than 1; about 0.5
!   beta  : a ratio of drag coef.; about 90.
!
!      real(8) :: x, xc
      real(r8) :: xc
      real(r8), intent(in) :: x
      real(r8), intent(out) :: r
      real(r8), parameter :: sig=1._r8, m=0.5_r8, beta=90._r8
!
      xc = 1._r8/(sig*m)
      IF (x.ge.xc) THEN
!        R_C = 999.           ! Full covered surface
        r   = 999._r8           ! Full covered surface
      ELSE
!        R_C = dsqrt(1.-sig*m*x)*dsqrt(1.+m*beta*x)
        r   = dsqrt(1._r8-sig*m*x)*dsqrt(1._r8+m*beta*x)
      ENDIF
!
!      END FUNCTION
      END subroutine r_c

! kang[2009/01/07] This is Fecan's correction from Dr. Jung
!----------------------------------------------------------------------
! A routine for correction of ust for soil moisture content
!
! w : volumetric soil moisture
! isl: soil texture type, ranging from 1 to 12
!
! Author: Yaping Shao, 5/05/2001
! Reference: Fecan et al. (1999), Ann. Geophysicae,17,149-157
!
! Data based on Shao and Jung, 2000, unpublished manuscript
! Data invented for silty loam, sandy clay loam, silty clay loamn sandy
! claym silty clay
!                        isl=5, 6, 8, 9, 10, 12
!----------------------------------------------------------------------
      subroutine h_c (w, isl, h)

!      real(8) :: w, wr(12), a(12), b(12)
      real(r8) :: wr(12), a(12), b(12)
      real(r8), intent(in)  :: w
      real(r8), intent(out) :: h
      integer, intent(in)  :: isl

      data wr/0.005_r8, 0.01_r8, 0.037_r8, 0.049_r8, 0.059_r8, 0.075_r8, 0.095_r8, 0.110_r8, 0.125_r8, 0.140_r8, 0.156_r8, 0.171_r8/
      data a /21.19_r8, 30.0_r8, 44.87_r8, 17.79_r8, 21.79_r8, 25.79_r8, 29.86_r8, 27.50_r8, 25.20_r8, 22.90_r8, 20.47_r8, 20.47/
      data b /0.68_r8,  0.90_r8, 0.85_r8,  0.61_r8,  0.67_r8,  0.74_r8,  0.8_r8,   0.75_r8, 0.70_r8,  0.65_r8,  0.59_r8,  0.54/

      if ( w.lt.0._r8 ) then
!        print *, 'h_c, w = ', w, ' <  0'
         write(6,*) 'h_c, w = ', w, ' <  0'
        stop
      endif

      if ( w.le.wr(isl) ) then
         h = 1.0_r8
      else
         h = sqrt( 1._r8 + a(isl)*( w-wr(isl) )**b(isl) )
      endif

      END subroutine h_c

#endif

end module DustEmis
