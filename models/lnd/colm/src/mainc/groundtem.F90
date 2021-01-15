#include <define.h>
 subroutine groundtem (itypwat,lb,nl_soil,dtime,n_pft,num_filterp,filterp,wt_patch, &
                       wt_column,capr,cnfac,csol,porsl,dkmg,dkdry,dksatu, &
                       sigf,dz,z,zi,tss,wice,wliq,scv,snowdp, &
                       frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg,&
#if (defined FHNP) && (defined FTF)
                       frostdp, thawdp, frostdp0,D_temperature,N_time,&
                       frost_day, thaw_day, &!liruichaoadd          
#endif
                       imelt,sm,xmf,fact,phi0,bsw)
! Snow and soil temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
!   the formulation used in SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction
!   in 10 soil layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two neighbor layers
!   (j, j+1) are derived from an assumption that the flux across the interface
!   is equal to that from the node j to the interface and the flux from the
!   interface to the node j+1. The equation is solved using the Crank-Nicholson
!   method and resulted in a tridiagonal system equation.
!
! Phase change (see meltf.F90)
! 
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!=======================================================================

  use precision
  use phycon_module, only : stefnc,tfrz
!use debug
  implicit none
  integer, INTENT(in) :: lb             !lower bound of array
  integer, INTENT(in) :: nl_soil        !upper bound of array
  integer, INTENT(in) :: n_pft          !number of pfts in a single column
  integer, INTENT(in) :: itypwat        !land water type (0=soil,1=urban or built-up,2=wetland,
                                        !3=land ice, 4=deep lake, 5=shallow lake)
  integer, INTENT(in) :: num_filterp    !
  integer, INTENT(in) :: filterp(n_pft) !
  real(r8), INTENT(in) :: wt_patch(n_pft) ! relative weight of pfts to grid area
  real(r8), INTENT(in) :: wt_column     ! relative weight of column to grid area

  real(r8), INTENT(in) :: dtime    !model time step [second]
  real(r8), INTENT(in) :: capr     !tuning factor to turn first layer T into surface T
  real(r8), INTENT(in) :: cnfac    !Crank Nicholson factor between 0 and 1

  real(r8), INTENT(in) :: csol(1:nl_soil)  !heat capacity of soil solids [J/(m3 K)]
  real(r8), INTENT(in) :: porsl(1:nl_soil) !soil porosity [-]
  real(r8), INTENT(in) :: dkmg(1:nl_soil)  !thermal conductivity of soil minerals [W/m-K]
  real(r8), INTENT(in) :: dkdry(1:nl_soil) !thermal conductivity of dry soil [W/m-K]
  real(r8), INTENT(in) :: dksatu(1:nl_soil)!thermal conductivity of saturated soil [W/m-K]
  real(r8), INTENT(in) :: phi0(1:nl_soil)  !soil water suction, negative potential [mm] !niu
  real(r8), INTENT(in) :: bsw(1:nl_soil)   !clapp and hornbereger "b" parameter [-]     !niu

  real(r8), INTENT(in) :: sigf(n_pft)      !fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: dz(lb:nl_soil)   !layer thickiness [m]
  real(r8), INTENT(in) :: z (lb:nl_soil)   !node depth [m]
  real(r8), INTENT(in) :: zi(lb-1:nl_soil) !interface depth [m]

  real(r8), INTENT(in) :: sabg(n_pft)     !solar radiation absorbed by ground [W/m2]
  real(r8), INTENT(in) :: frl             !atmospheric infrared (longwave) radiation [W/m2]
  real(r8), INTENT(in) :: dlrad(n_pft)    !downward longwave radiation blow the canopy [W/m2]
  real(r8), INTENT(in) :: fseng(n_pft)    !sensible heat flux from ground [W/m2]
  real(r8), INTENT(in) :: fevpg(n_pft)    !evaporation heat flux from ground [mm/s]
  real(r8), INTENT(in) :: cgrnd(n_pft)    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
  real(r8), INTENT(in) :: htvp            !latent heat of vapor of water (or sublimation) [j/kg]
  real(r8), INTENT(in) :: emg             !ground emissivity (0.97 for snow,

  real(r8), INTENT(inout) :: tss (lb:nl_soil) !soil temperature [K]
  real(r8), INTENT(inout) :: wice(lb:nl_soil) !ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq(lb:nl_soil) !liqui water [kg/m2]
  real(r8), INTENT(inout) :: scv              !snow cover, water equivalent [mm, kg/m2]
  real(r8), INTENT(inout) :: snowdp           !snow depth [m]

  real(r8), INTENT(out) :: sm               !rate of snowmelt [kg/(m2 s)]
  real(r8), INTENT(out) :: xmf              !total latent heat of phase change of ground water
  real(r8), INTENT(out) :: fact(lb:nl_soil) !used in computing tridiagonal matrix
  integer, INTENT(out) :: imelt(lb:nl_soil)    !flag for melting or freezing [-]
!------------------------ local variables ------------------------------
  real(r8) cv(lb:nl_soil)     ! heat capacity [J/(m2 K)]
  real(r8) tk(lb:nl_soil)     ! thermal conductivity [W/(m K)]

  real(r8) at(lb:nl_soil)     !"a" vector for tridiagonal matrix
  real(r8) bt(lb:nl_soil)     !"b" vector for tridiagonal matrix
  real(r8) ct(lb:nl_soil)     !"c" vector for tridiagonal matrix
  real(r8) rt(lb:nl_soil)     !"r" vector for tridiagonal solution

  real(r8) fn  (lb:nl_soil)   ! heat diffusion through the layer interface [W/m2]
  real(r8) fn1 (lb:nl_soil)   ! heat diffusion through the layer interface [W/m2]
  real(r8) dzm                ! used in computing tridiagonal matrix
  real(r8) dzp                ! used in computing tridiagonal matrix

  real(r8) tssbef(lb:nl_soil) ! soil/snow temperature before update
  real(r8) hs_pft(n_pft)      ! net energy flux into the surface on pfts(w/m2)
  real(r8) hs                 ! net energy flux into the surface on column(w/m2)
  real(r8) dhsdT_pft(n_pft)   ! d(hs)/dT on pft
  real(r8) dhsdT              ! d(hs)/dT on column
  real(r8) brr(lb:nl_soil)    ! temporay set
  integer i,j,p,fp

#if (defined FHNP) && (defined FTF)
  !liruichao add
  real(r8), INTENT(inout) :: frostdp         !frost depth
  real(r8), INTENT(inout) :: thawdp          !thaw deppth
  real(r8), INTENT(inout) :: frostdp0        !initial frost depth
  real(r8), INTENT(inout) :: D_temperature   !frost or thaw index
  real(r8), INTENT(inout) :: N_time          !step counter
  real(r8), INTENT(inout) :: frost_day       !frost days
  real(r8), INTENT(inout) :: thaw_day        !thaw days
  real(r8)  frostdp_last      !the frostdp of last step time
  real(r8)  thawdp_last       !the thawdp of last step time
  real(r8)  frostdp0_last     !initial frost duration days
  real(r8)  R_downward(1:10)  !dz/tk
  real(r8)  N_downward(1:10)  !The frost or thaw index of each layer
  real(r8)  sumR_downward(1:10)!The sum of each layer dz/tk
  real(r8)  sumN_downward(0:10)!The sum of frost or thaw index of each layer
  real(r8)  D_upper            !quadratic formula
  real(r8)  dzf                !frost depth in soil layer
  real(r8)  dzt                !thaw depth in soil layer
  integer, parameter :: Day=4  !determine to enter frost or thaw period
  real(r8)  sm_vfraction(1:10)   ! volumetric fraction of soil moisture
  real(r8) fact1(lb:nl_soil+1)   !used in computing tridiagonal matrix
  real(r8) at1(lb:nl_soil+1)     !"a" vector for tridiagonal matrix
  real(r8) bt1(lb:nl_soil+1)     !"b" vector for tridiagonal matrix
  real(r8) ct1(lb:nl_soil+1)     !"c" vector for tridiagonal matrix
  real(r8) rt1(lb:nl_soil+1)     !"r" vector for tridiagonal solution
  real(r8) fn2(lb:nl_soil+1)     ! heat diffusion through the layerinterface[W/m2]
  real(r8) fn2_1(lb:nl_soil+1)
  real(r8) tss1(lb:nl_soil+1)    ! soil temperature [K]
  real(r8) deltaz1               ! the thickness between frost depth and frontsoil layer
  real(r8) deltaz2               ! the thickness between frost depth and lastsoil layer
  integer j1
  real(r8) cv1(lb:nl_soil+1)     ! heat capacity [J/(m2 K)]
  real(r8) tk1(lb:nl_soil+1)     ! thermal conductivity [W/(m K)]
!==liruichao add for derivation of the 3 layers tk because of FTF==
  real(r8) thk1(lb:nl_soil+1)     ! thermal conductivity [W/(m K)]
  real(r8) thkhc(lb:nl_soil)     ! thermal conductivity [W/(m K)]
  real(r8) z1 (lb:nl_soil+1)   !node depth [m] !real(r8) fn2(lb:nl_soil+1)
  real(r8) zi1 (lb-1:nl_soil+1)   !node depth [m] !real(r8) fn2(lb:nl_soil+1)
!==liruichao add for derivation of the 3 layers tk because of FTF==
  real(r8) tssbef1(lb:nl_soil+1)   ! heat diffusion through the layer interface[W/m2]

  real(r8)  h1
  real(r8)  h2
  real(r8)  h3
  real(r8)  h4


#endif
!=======================================================================

! heat capacity 
      call hCapacity (itypwat,lb,nl_soil,csol,porsl,wice,wliq,scv,dz,cv)

! thermal conductivity
#if (defined FHNP) && (defined FTF)
!====liruichao added output for thkhc=====
      call hConductivity (itypwat,lb,nl_soil,&
                          dkmg,dkdry,dksatu,porsl,dz,z,zi,tss,wice,wliq,tk,thkhc)
#else
!====liruichao added output for thkhc=====
      call hConductivity (itypwat,lb,nl_soil,&
                          dkmg,dkdry,dksatu,porsl,dz,z,zi,tss,wice,wliq,tk)
#endif


#if (defined FHNP) && (defined FTF)
!liruichao  add
!begin

         if(N_time<=71)then
                N_time=N_time+1
                D_temperature=D_temperature+(tss(lb)-tfrz)*dtime
        end if
        if(N_time==72)then
          N_time=0
          frostdp_last=frostdp
          thawdp_last=thawdp        !the frostdp and thawdp of last step time

          do i=1,10
                            sm_vfraction(i)=wliq(i)/(1000*dz(i))+wice(i)/(917*dz(i))
          end do

          do i=1,10
            R_downward(i)=dz(i)/tk(i)
          end do
          N_downward(1:10)=0
          N_downward(1)=335000*917*sm_vfraction(1)*dz(1)*dz(1)/tk(1)/2!Li=335000*917
          do i=2,10
            do j=1,i-1
              N_downward(i)=335000*917*sm_vfraction(i)*dz(i)*R_downward(j)+N_downward(i)
            end do
            N_downward(i)=N_downward(i)+335000*917*sm_vfraction(i)*dz(i)*dz(i)/tk(i)/2
          end do

          sumR_downward(1:10)=0
            do i=2,10
              do j=1,i-1
                sumR_downward(i)=sumR_downward(i)+R_downward(j)
              end do
            end do
            sumN_downward(0:10)=0
              do i=1,10
                do j=1,i
                  sumN_downward(i)=sumN_downward(i)+N_downward(j)
                end do
              end do

            if(D_temperature<=0)then
                if(frost_day<Day)then
                        frost_day=frost_day+1
                        D_upper=0
                end if

                if(frost_day==Day)then
                        thaw_day=0

                        if(thawdp==0)then
                                if(frostdp==0)then
                                        D_upper=D_temperature
                                end if 
                                if(frostdp>=zi(10))then
                                        D_upper=0
                                end if

                                if(frostdp>0.and.frostdp<zi(10))then
                                        do j=1,10
                      if(frostdp>=zi(j-1).and.frostdp<zi(j))then
                        D_upper=((frostdp-zi(j-1)+tk(j)*sumR_downward(j))**2-tk(j)**2*sumR_downward(j)**2)*(335000*917*sm_vfraction(j))/(2*tk(j))+sumN_downward(j-1)-D_temperature
                        D_upper=-D_upper
                        exit
                      end if
                    end do
                                 

                                end if
                        end if
                        if(thawdp>0.and.thawdp<=zi(10))then
                                do j=1,10
                  
                     if(thawdp>=zi(j-1).and.thawdp<zi(j))then
                    D_upper=((thawdp-zi(j-1)+tk(j)*sumR_downward(j))**2-tk(j)**2*sumR_downward(j)**2)*(335000*917*sm_vfraction(j))/(2*tk(j))+sumN_downward(j-1)+D_temperature
                     ! D_upper=-D_upper
                      exit
                    end if
                  end do
                              

                if(D_upper<=0)then
                    D_upper=0
                    thawdp=0     ! D_upper=
                  end if
                end if
            end if
        end if
              if(D_temperature>0)then
                if(thaw_day<Day)then
                        thaw_day=thaw_day+1
                        D_upper=0
                end if

                if(thaw_day==Day)then
                        frost_day=0

                        if(thawdp==0)then
                                if(frostdp==0)then
                                        D_upper=0
                                end if
                      if(frostdp>0.and.frostdp<=zi(10))then
                                        do j=1,10
                      if(thawdp>=zi(j-1).and.thawdp<zi(j))then
                        D_upper=((thawdp-zi(j-1)+tk(j)*sumR_downward(j))**2-tk(j)**2*sumR_downward(j)**2)*(335000*917*sm_vfraction(j))/(2*tk(j))+sumN_downward(j-1)+D_temperature
                        exit
                      end if
                    end do

                                end if
                        end if

                        if(thawdp>0.and.thawdp<zi(10))then
  
                                  do j=1,10
                      if(thawdp>=zi(j-1).and.thawdp<zi(j))then
                        D_upper=((thawdp-zi(j-1)+tk(j)*sumR_downward(j))**2-tk(j)**2*sumR_downward(j)**2)*(335000*917*sm_vfraction(j))/(2*tk(j))+sumN_downward(j-1)+D_temperature
                        exit
                      end if
                    end do

                  end if
                                     end if
                
          end if
           if(D_upper<0)then
              D_upper=-D_upper

                do j=1,10
                  if(D_upper>sumN_downward(j-1).and.D_upper<=sumN_downward(j))then
                    dzf=-tk(j)*sumR_downward(j)+sqrt(tk(j)**2*sumR_downward(j)**2+2*tk(j)*(D_upper-sumN_downward(j-1))/(335000*917*sm_vfraction(j)))
                    frostdp=dzf+zi(j-1)
                    end if
    
              end do
  
                    if(D_upper>sumN_downward(10))then
                      frostdp=zi(10)
                    end if
                 

                if(frostdp_last/=thawdp_last.and.((frostdp_last-thawdp_last)*(frostdp-thawdp)<=0))then
                  frostdp=0
                  thawdp=0
                !end if
              end if
              D_upper=-D_upper
            end if

           if(D_upper>0)then
              do j=1,10
                if(D_upper>sumN_downward(j-1).and.D_upper<=sumN_downward(j))then
                  dzt=-tk(j)*sumR_downward(j)+sqrt(tk(j)**2*sumR_downward(j)**2+2*tk(j)*(D_upper-sumN_downward(j-1))/(335000*917*sm_vfraction(j)))
                  thawdp=dzt+zi(j-1)

          !   exit
                   end if
           !    exit
             end do

           

                  if(D_upper>sumN_downward(10))then
                     thawdp=zi(10)
                  end if
                !end if
              !end do

              if(frostdp_last/=thawdp_last.and.((frostdp_last-thawdp_last)*(frostdp-thawdp)<=0))then
                frostdp=0
                thawdp=0
              end if
            end if
             D_temperature=0
        end if


!liruichao end
#endif  


! net ground heat flux into the surface and its temperature derivative
! Added a pfts loop here to get the average of hs and dhsdT over all 
! PFTs on the column  add by zhq. 07/20/2009
      hs = 0.
      dhsdT = 0.

      do fp = 1, num_filterp
          p = filterp(fp)

        hs_pft(p) = sabg(p) + dlrad(p) &
           + (1.-sigf(p))*emg*frl - emg*stefnc*tss(lb)**4 &
           - (fseng(p)+fevpg(p)*htvp) 

        dhsdT_pft(p) = - cgrnd(p) - 4.*emg * stefnc * tss(lb)**3
        hs = hs + hs_pft(p) * wt_patch(p) / wt_column
        dhsdT = dhsdT + dhsdT_pft(p) *wt_patch(p) / wt_column
      end do

      tssbef(lb:) = tss(lb:)
      j       = lb
      fact(j) = dtime / cv(j) &
              * dz(j) / (0.5*(z(j)-zi(j-1)+capr*(z(j+1)-zi(j-1))))

      do j = lb + 1, nl_soil
        fact(j) = dtime/cv(j)
      enddo

      do j = lb, nl_soil - 1
        fn(j) = tk(j)*(tss(j+1)-tss(j))/(z(j+1)-z(j))
      enddo
      fn(nl_soil) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
      j     = lb
      dzp   = z(j+1)-z(j)
      at(j) = 0.
      bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
      ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = tss(j) + fact(j)*( hs - dhsdT*tss(j) + cnfac*fn(j) )


      do j = lb + 1, nl_soil - 1
         dzm   = (z(j)-z(j-1))
         dzp   = (z(j+1)-z(j))
         at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
         bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
         ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
         rt(j) = tss(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
      end do

      j     =  nl_soil
      dzm   = (z(j)-z(j-1))
      at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
      ct(j) = 0.
      rt(j) = tss(j) - cnfac*fact(j)*fn(j-1)

! solve for tss
      i = size(at)
      call tridia (i ,at ,bt ,ct ,rt ,tss) 


!=======================================================================
! melting or freezing 
!=======================================================================
!#if 0
!      do j = lb, nl_soil - 1
!         fn1(j) = tk(j)*(tss(j+1)-tss(j))/(z(j+1)-z(j))
!      enddo
!      fn1(nl_soil) = 0.
!
!      j = lb
!      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)
!
!      do j = lb + 1, nl_soil
!         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
!      enddo
!#endif

#if (defined FHNP) && (defined FTF)
!liruichao add Calculating soil temperature

        deltaz1=0.0
        deltaz2=0.0

        j1=0.0
        at1=0.0_r8
        bt1=0.0_r8
        ct1=0.0_r8
        rt1=0.0_r8

        z1=0.0_r8
        zi1=0.0_r8
        tssbef1  =0.0_r8
        tk1 =0.0_r8
        thk1 =0.0_r8
        cv1 =0.0_r8
        fact1 =0.0_r8
        fn2 =0.0_r8

        h1=0.0_r8
        h2=0.0_r8
        h3=0.0_r8
        h4=0.0_r8
        tss1=0.0_r8
 
    if(frostdp>=z(2).and.frostdp<z(9))then
         do j=2,9
           if(frostdp>=z(j).and.frostdp<z(j+1))then
                  j1=j
           endif
         enddo

     if((tssbef(j1)-tfrz)*(tssbef(j1+1)-tfrz)<0)then
               deltaz1=frostdp-z(j1)
               deltaz2=z(j1+1)-frostdp

           if(deltaz1>=0.01.and.deltaz2>=0.01)then
                z1(lb:j1)=z(lb:j1)
                z1(j1+1)=frostdp
                z1(j1+2:nl_soil+1)=z(j1+1:nl_soil)
                !h1,h2,h3,h4 4 sections that occupie the new layer
                !h1,h4 are the layers of neighbor layer
                !h2,h3 are the layers of the frostdp Jinbo Xie
                h1=0.5*(frostdp+z(j1))-zi(j1-1)
                h2=zi(j1)-0.5*(frostdp+z(j1))
                h3=0.5*(z(j1+1)+frostdp)-zi(j1)
                h4=zi(j1+1)-0.5*(z(j1+1)+frostdp)
                !set prepare temperature
                tssbef1(lb:j1-1)=tssbef(lb:j1-1)
                tssbef1(j1)=(tssbef(j1)*(h1+h2)-tfrz*h2)/h1
                tssbef1(j1+1)=tfrz
                tssbef1(j1+2)=(tssbef(j1+1)*(h3+h4)-tfrz*h3)/h4
                tssbef1(j1+3:nl_soil+1)=tssbef(j1+2:nl_soil)
!====liruichao add for 3 layers tk because of FTF====
!====this uses the contraint that the flux pass into the interface====
!====equals the flux that pass out the interface====
!====the tk is thus derived from upper and lower node tk as a function====
!====of depth dz====
                !thermal conductivity
                zi1(lb-1:j1-1)=zi(lb-1:j1-1)
                zi1(j1)=zi(j1)-h2
                zi1(j1+1)=zi(j1)+h3
                zi1(j1+2:nl_soil+1)=zi(j1+1:nl_soil)

                thk1(lb:j1-1)=thkhc(lb:j1-1)
                thk1(j1)=h1*thkhc(j1)/(h2+h1)
                thk1(j1+1)=h2*thkhc(j1)/(h2+h1)+h3*thkhc(j1+1)/(h3+h4)
                thk1(j1+2)=h4*thkhc(j1+1)/(h3+h4)
                thk1(j1+3:nl_soil+1)=thkhc(j1+2:nl_soil)
!====Jinbo Xie====
!====the linear1 separation of thk1 is a simple estimate since the====
!====presence of dkdry term in the estimation of thk not fully linear====
!====may be further modified later when tk and hk are further====
!====in the effect of ftf====
                tk1(lb:j1-2)=tk(lb:j1-2)
                do j=j1-1,j1+2
                tk1(j) = thk1(j)*thk1(j+1)*(z1(j+1)-z1(j)) &
                             /(thk1(j)*(z1(j+1)-zi1(j))+thk1(j+1)*(zi1(j)-z1(j)))
                enddo
                tk1(j1+3:nl_soil+1)=tk(j1+2:nl_soil)
!====liruichao add for 3 layers tk because of FTF====
                !heat capacity
                cv1(lb:j1-1)=cv(lb:j1-1)
                cv1(j1)=((z(j1)+frostdp)/2-zi(j1-1))&
                      *cv(j1)/dz(j1)
                cv1(j1+1)=cv(j1)/dz(j1)* (zi(j1)&
                        -(z(j1)+frostdp)/2)&
                        +cv(j1+1)/dz(j1+1)&
                       *((z(j1+1)+frostdp)/2-zi(j1))
                cv1(j1+2)=(zi(j1+1)&
                               -(z(j1+1)+frostdp)/2)&
                               *cv(j1+1)/dz(j1+1)
                cv1(j1+3:nl_soil+1)=cv(j1+2:nl_soil)
!====Jinbo Xie====
!====the linear separation of cv implicitly====
!====assumes a linear separation of wliq and ====
!====wice in the new 3 layers separated by the FTF====
!====================================================
                j       = lb
                  fact1(j) = dtime / cv1(j) &
                             * dz(j) / (0.5*(z(j)-zi(j-1)+capr*(z(j+1)-zi(j-1))))
             do j = lb+1, nl_soil+1
                fact1(j) = dtime/cv1(j)
             enddo

             do j = lb, nl_soil
                fn2(j) = tk1(j)*(tssbef1(j+1)-tssbef1(j))/(z1(j+1)-z1(j))
             enddo
                fn2(nl_soil+1) = 0.

             j     = lb
             dzp   = z(j+1)-z(j)
             at1(j) = 0.
             bt1(j) = 1+(1.-cnfac)*fact1(j)*tk1(j)/dzp-fact1(j)*dhsdT
             ct1(j) =  -(1.-cnfac)*fact1(j)*tk1(j)/dzp
             rt1(j) = tssbef1(j) + fact1(j)*( hs - dhsdT*tssbef1(j)&
                                 + cnfac*fn2(j) )
         do  j =  lb+1, nl_soil
         dzm   = (z1(j)-z1(j-1))
         dzp   = (z1(j+1)-z1(j))
         at1(j) =   - (1.-cnfac)*fact1(j)* tk1(j-1)/dzm
         bt1(j) = 1.+ (1.-cnfac)*fact1(j)*(tk1(j)/dzp + tk1(j-1)/dzm)
         ct1(j) =   - (1.-cnfac)*fact1(j)* tk1(j)/dzp
         rt1(j) = tssbef1(j) + cnfac*fact1(j)*( fn2(j) - fn2(j-1) )
         end do

              j   =  nl_soil+1
           dzm    = (z1(j)-z1(j-1))
           at1(j) =   - (1.-cnfac)*fact1(j)*tk1(j-1)/dzm
           bt1(j) = 1.+ (1.-cnfac)*fact1(j)*tk1(j-1)/dzm
           ct1(j) = 0.
           rt1(j) = tssbef1(j) - cnfac*fact1(j)*fn2(j-1)

      i = size(at1)
      call tridia (i ,at1 ,bt1 ,ct1 ,rt1 ,tss1)

      do j=lb, j1-1
        tss(j) = tss1(j)
      enddo
      
      tss(j1)=tss1(j1)*h1/dz(j1) +tss1(j1+1)*h2/dz(j1)
      tss(j1+1)=tss1(j1+1)*h3/dz(j1+1) +tss1(j1+2)*h4/dz(j1+1)
            
      do j = j1+3, nl_soil+1
        tss(j-1) = tss1(j)
      enddo


                   endif
                    
        endif
          
endif

!end
!-----------------------------------------------------------------------------------------------------
#endif
!=======================================================================
! melting or freezing 
!=======================================================================
!#if 0
#if (defined FHNP) && (defined FTF)
!no need for brr
!Jinbo Xie
#else
      do j = lb, nl_soil - 1
         fn1(j) = tk(j)*(tss(j+1)-tss(j))/(z(j+1)-z(j))
      enddo
      fn1(nl_soil) = 0.

      j = lb
      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

      do j = lb +1, nl_soil
         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
      enddo
#endif
      call meltf ( lb, nl_soil, dtime, &
                   fact(lb:), brr(lb:), hs, dhsdT, &
                   tssbef(lb:), tss(lb:), wliq(lb:), wice(lb:), imelt(lb:), &
                   scv, snowdp, sm, xmf, porsl(1:), phi0(1:), bsw(1:), dz(1:) ) !niu

!-----------------------------------------------------------------------

 end subroutine groundtem
