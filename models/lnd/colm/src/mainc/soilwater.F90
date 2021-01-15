
 subroutine soilwater ( nl_soil,dtime,wimp,smpmin, &
                        porsl,phi0,bsw,hksati, &
                        z,dz,zi,tss,vol_liq,vol_ice,eff_porosity, &
                        qinfl,etr,rootr,dwat,hk,dhkdw,fice)     !niu

!=======================================================================
! Original author : Yongjiu Dai, September 15, 1999
!
! Soil moisture is predicted from a 10-layer model (as with soil temperatures),
! in which the vertical soil moisture transport is governed by infiltration,
! runoff, gradient diffusion, gravity, and root extraction by canopy transpiration;
! the net water applied to the surface layer is by snowmelt, precipitation,
! the throughfall of canopy dew, and minus surface runoff and evaporation.
! 
! The vertical water flow in an unsaturated porous media is described by Darcy's law,
! and the hydraulic conductivity and the soil negative potential vary with
! soil water content and soil texture based on the work of Clapp and Hornberger (1978)
! and Cosby et al. (1984). The equation is integrated over the layer thickness,
! in which the time rate of change in water mass must equal the net flow across
! the bounding interface, plus the rate of internal source or sink.
! The terms of water flow across the layer interfaces are linearly expanded
! by using first-order Taylor expansion. The equations result in a tridiagonal
! system equation.
!
!             Note: length units here are all millimeter 
!                   (in temperature subroutine uses same  soil layer 
!                   structure required but lengths are m)
!
!                   Richards equation:
!
!                   d wat     d     d wat d psi
!                   ----- =  -- [ k(----- ----- - 1) ] + S
!                     dt     dz       dz  d wat
!
!                   where: wat = volume of water per volume of soil (mm**3/mm**3)
!                          psi = soil matrix potential (mm)
!                          dt  = time step (s)
!                          z   = depth (mm)
!                          dz  = thickness (mm)
!                          qin = inflow at top (mm h2o /s) 
!                          qout= outflow at bottom (mm h2o /s)
!                          s   = source/sink flux (mm h2o /s) 
!                          k   = hydraulic conductivity (mm h2o /s)
!
!                                 d qin                 d qin
!           qin[n+1] = qin[n] +  -------- d wat(j-1) + --------- d wat(j)
!                                 d wat(j-1)            d wat(j)
!                          ==================|================= 
!                                            < qin 
!
!                          d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j) 
!
!                                            > qout
!                          ==================|================= 
!                                    d qout               d qout
!             qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!                                    d wat(j)             d wat(j+1)
!
!
!                   Solution: linearize k and psi about d wat and use tridiagonal 
!                             system of equations to solve for d wat, 
!                             where for layer j
!
!                   r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
!=======================================================================

  use precision
  use debug
  use phycon_module, only : tfrz
  implicit none
!     
!------------------------- dummy argument ------------------------------
  integer, INTENT(in) :: nl_soil  ! soil layers

  real(r8), INTENT(in) :: &
        dtime,            &! time step [s]
        wimp,             &!water impremeable if porosity less than wimp
        smpmin,           &!restriction for min of soil poten. (mm)

        phi0(1:nl_soil),  &! saturated soil suction [mm]
        bsw(1:nl_soil),   &! Clapp-Hornberger "B"
        porsl(1:nl_soil), &! saturated volumetric soil water content(porosity)
        hksati(1:nl_soil),&! hydraulic conductivity at saturation [mm h2o/s]
        z(1:nl_soil),     &! layer depth [mm]
        dz(1:nl_soil),    &! layer thickness [mm]
        zi(0:nl_soil),    &! interface level below a "z" level [mm]
 
        qinfl,            &! net water input from top [mm h2o/s]
        etr,              &! actual transpiration [mm h2o/s]
        rootr(1:nl_soil), &! root resistance of a layer, all layers add to 1.0
 eff_porosity(1:nl_soil), &! effective porosity
      vol_liq(1:nl_soil), &! soil water per unit volume [mm/mm]
      vol_ice(1:nl_soil), &! soil water per unit volume [mm/mm]
         fice(1:nl_soil), &! impermeable area due to ice (-)            !niu
      tss(1:nl_soil)       ! soil temperature

  real(r8), INTENT(out) :: &
        dwat (1:nl_soil), &! change of soil water [m3/m3]
        hk   (1:nl_soil), &! hydraulic conductivity [mm h2o/s]
        dhkdw(1:nl_soil)   ! d(hk)/d(vol_liq)
                   
!-------------------------- local variables -----------------------------
  real(r8) :: &
       amx(1:nl_soil),    &! "a" left off diagonal of tridiagonal matrix
       bmx(1:nl_soil),    &! "b" diagonal column for tridiagonal matrix
       cmx(1:nl_soil),    &! "c" right off diagonal tridiagonal matrix
       den,               &! used in calculating qin, qout
       dqidw0,            &! d(qin)/d(vol_liq(i-1))
       dqidw1,            &! d(qin)/d(vol_liq(i))
       dqodw1,            &! d(qout)/d(vol_liq(i))
       dqodw2,            &! d(qout)/d(vol_liq(i+1))
       dsmpdw(1:nl_soil), &! d(smp)/d(vol_liq)
       num,               &! used in calculating qin, qout
       qin,               &! flux of water into soil layer [mm h2o/s]
       qout,              &! flux of water out of soil layer [mm h2o/s]
       rmx(1:nl_soil),    &! "r" forcing term of tridiagonal matrix
       s_node,            &! soil wetness
       s1,                &! "s" at interface of layer
       s2,                &! k*s**(2b+2)
       smp(1:nl_soil)      ! soil matrix potential [mm]

  real(r8) :: &
       liqvm,             &! two-layer mean volumetric liq (m3/m3) 
       icevm,             &! two-layer mean volumetric ice (m3/m3)
       porem,             &! two-layer mean porosity (m3/m3)
       ficem               ! two-layer mean impermeable fraction (m3/m3)

  integer j                ! do loop indices

!=======================================================================
! Set zero to hydraulic conductivity if effective porosity 5% in any of 
! two neighbour layers or liquid content (theta) less than 0.001

!niu      do j = 1, nl_soil
!niu         if((eff_porosity(j) < wimp) &
!niu             .OR. (eff_porosity(min(nl_soil,j+1)) < wimp) &
!niu             .OR. (vol_liq(j) <= 1.e-3))then
!niu           hk(j) = 0.
!niu           dhkdw(j) = 0.
!niu         else
!niu           s1 = 0.5*(vol_liq(j)/porsl(j)&
!niu                   + vol_liq(min(nl_soil,j+1))/porsl(min(nl_soil,j+1)))
!niu           s2 = hksati(j)*s1**(2.*bsw(j)+2.)
!niu
!niu           hk(j) = s1*s2  
!niu
!niu           dhkdw(j) = (2.*bsw(j)+3.)*s2*0.5/porsl(j)
!niu           if(j == nl_soil) dhkdw(j) = dhkdw(j) * 2.
!niu         endif
!niu           
!niu      end do
! Evaluate hydraulic conductivity, soil matrix potential,
!     d(smp)/d(vol_liq), and d(hk)/d(vol_liq).

!niu      do j = 1, nl_soil
!niu         if(tss(j)>tfrz) then
!niu            if(porsl(j)<1.e-6)then     ! bed rock
!niu               s_node = 0.001
!niu               smp(j) = -phi0(j)
!niu               dsmpdw(j) = 0.
!niu            else
!niu               s_node = min(1.,vol_liq(j)/porsl(j))
!niu               s_node = max(s_node, 0.001)
!niu               smp(j) = -phi0(j)*s_node**(-bsw(j))
!niu               smp(j) = max(smpmin, smp(j))        ! Limit soil suction
!niu               dsmpdw(j) = -bsw(j)*smp(j)/(s_node*porsl(j))
!niu            endif
!niu         else
! When ice is present, the matric potential is only related to temperature
! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
! Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm 
!niu            smp(j) = 1.e3 * 0.3336e6/9.80616*(tss(j)-tfrz)/tss(j)
!niu            smp(j) = max(smpmin, smp(j))        ! Limit soil suction
!niu            dsmpdw(j) = 0.
!niu         endif
!niu
!niu      end do

!niuin (see Niu&Yang, 2006 for explanations)

      do j = 1, nl_soil
           liqvm = 0.5_r8*(vol_liq(j) + vol_liq(min(nl_soil,j+1)) )
           icevm = 0.5_r8*(vol_ice(j) + vol_ice(min(nl_soil,j+1)) )
           porem = 0.5_r8*(porsl  (j) + porsl  (min(nl_soil,j+1)) )
           ficem = 0.5_r8*(fice   (j) + fice   (min(nl_soil,j+1)) )

           if(porem > 0.) then   !jidy
              s1 = min(1._r8,(liqvm+icevm)/porem)
           else
              s1 = 0.            !jidy
           end if
           s2 = hksati(j)*s1**(2._r8*bsw(j)+2._r8)

           hk(j) = (1._r8-ficem) * s1*s2  

           if(porsl(j) > 0.) then  !jidy
              dhkdw(j) = (1._r8-ficem) * (2._r8*bsw(j)+3._r8)*s2*0.5/porsl(j)
           else
              dhkdw(j) = 0.
           end if

         ! if(j == nl_soil) dhkdw(j) = dhkdw(j) * 2._r8        !comment by jidy
           if(j > 1) then
              if (porsl(j-1) > 0. .and. porsl(j) < 1.e-6) then !add by jidy
                 dhkdw(j) = dhkdw(j) * 2._r8                   !add by jidy
              end if
           end if

           if(porsl(j) > 0.) then
              s_node = min(1._r8,(vol_liq(j)+vol_ice(j))/porsl(j))
           else
              s_node = 0.
           end if
           s_node = max(s_node, 0.001_r8)
           smp(j) = -phi0(j)*s_node**(-bsw(j))
           smp(j) = max(smpmin, smp(j))        ! Limit soil suction

           if(porsl(j) > 0.) then                              !add by jidy
              dsmpdw(j) = -bsw(j)*smp(j)/(s_node*porsl(j))
           else
              dsmpdw(j) = 0.                                   !add by jidy
           end if
      end do

!niuout
!-----------------------------------------------------------------------
! Set up r, a, b, and c vectors for tridiagonal solution
! Node j=1

      j      = 1
      qin    = qinfl

      den    = (z(j+1)-z(j))
      num    = (smp(j+1)-smp(j)) - den
      qout   = -hk(j)*num/den
      dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
      dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den

      rmx(j) =  qin - qout - etr*rootr(j)
      amx(j) =  0.
      bmx(j) =  dz(j)/dtime + dqodw1
      cmx(j) =  dqodw2
 
!-----------------------------------------------------------------------
! Nodes j=2 to j=nl_soil-1

      do j = 2, nl_soil - 1
         den    = (z(j) - z(j-1))
         num    = (smp(j)-smp(j-1)) - den
         qin    = -hk(j-1)*num/den
         dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
         dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den

         den    = (z(j+1)-z(j))
         num    = (smp(j+1)-smp(j)) - den
         qout   = -hk(j)*num/den
         dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den
         dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den

         rmx(j) =  qin - qout - etr*rootr(j)
         amx(j) = -dqidw0
         bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1
         cmx(j) =  dqodw2
      end do
 
!-----------------------------------------------------------------------
! node j=nl_soil
      j      = nl_soil
      den    = (z(j) - z(j-1))
      num    = (smp(j)-smp(j-1)) - den
      qin    = -hk(j-1)*num/den
      dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den
      dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den

!niu      qout   =  hk(j)
!niu      dqodw1 =  dhkdw(j)
      qout   =  0.      !zero-flux bottom boundary condition
      dqodw1 =  0.      !zero-flux bottom boundary condition

      rmx(j) =  qin - qout - etr*rootr(j)
      amx(j) = -dqidw0
      bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1
      cmx(j) =  0.


!-----------------------------------------------------------------------
! Solve for dwat

      call tridia (nl_soil ,amx ,bmx ,cmx ,rmx ,dwat)

 end subroutine soilwater
