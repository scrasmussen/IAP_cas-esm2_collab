#include <define.h>
#if (defined FHNP) && (defined FTF)
 subroutine hConductivity (itypwat,lb,nl_soil,& 
                           dkmg,dkdry,dksatu,porsl,dz,z,zi,tss,wice,wliq,tk,thkhc)
#else
subroutine hConductivity (itypwat,lb,nl_soil,&
                          dkmg,dkdry,dksatu,porsl,dz,z,zi,tss,wice,wliq,tk)
#endif

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, September 15, 1999
!
! calculation of thermal conductivities of snow / soil layers
! The thermal conductivity of soil is computed from
! the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
! the formulation used in SNTHERM (Jordan 1991).
!
! The thermal conductivities at the interfaces between two neighbor layers
! (j, j+1) are derived from an assumption that the flux across the interface
! is equal to that from the node j to the interface and the flux from the
! interface to the node j+1.
!-----------------------------------------------------------------------

  use precision
  use phycon_module, only : denh2o,denice,tfrz,tkwat,tkice,tkair
  implicit none

  integer, INTENT(in) :: lb            ! lower bound of array
  integer, INTENT(in) :: nl_soil       ! upper bound of array
  integer, INTENT(in) :: itypwat       ! land water type (0=soil, 1=urban, 2=wetland,
                                       ! 3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) ::   dkmg(1:nl_soil)  ! dkm**dmvol, where dkm is the mineral cond.
  real(r8), INTENT(in) ::  dkdry(1:nl_soil)  ! thermal conductivity for dry soil [W/m-K]
  real(r8), INTENT(in) :: dksatu(1:nl_soil)  ! Thermal conductivity of saturated soil [W/m-K]
  real(r8), INTENT(in) ::  porsl(1:nl_soil)  ! fractional volume between soil grains=1.-dmvol
  real(r8), INTENT(in) ::   dz(lb:nl_soil)   ! layer thickiness [m]
  real(r8), INTENT(in) ::    z(lb:nl_soil)   ! node depth [m]
  real(r8), INTENT(in) ::   zi(lb-1:nl_soil) ! interface depth [m]
  real(r8), INTENT(in) ::  tss(lb:nl_soil)   ! Nodal temperature [K]
  real(r8), INTENT(in) :: wice(lb:nl_soil)   ! ice lens [kg/m2]
  real(r8), INTENT(in) :: wliq(lb:nl_soil)   ! liqui water [kg/m2]
  real(r8), INTENT(out) :: tk(lb:nl_soil)    ! thermal conductivity [W/(m K)]
!====liruichao add====
#if (defined FHNP) && (defined FTF)
  real(r8), INTENT(out) :: thkhc(lb:nl_soil)
#endif
!====liruichao add====
! local
  real(r8) rhosnow      ! partitial density of water (ice + liquid)
  real(r8) dksat        ! thermal conductivity for saturated soil (j/(k s m))
  real(r8) dke          ! kersten number
  real(r8) fl           ! fraction of liquid or unfrozen water to total water
  real(r8) satw         ! relative total water content of soil.
  real(r8) thk(lb:nl_soil)  ! thermal conductivity of layer

  integer i

!-----------------------------------------------------------------------
! Thermal conductivity of soil from Farouki (1981),
      do i = 1, nl_soil

         if(itypwat<=1)then         !soil ground
            thk(i) = dkdry(i)       !rock or dry soil

            if(porsl(i)>1.e-05)then
               satw = (wliq(i)/denh2o+wice(i)/denice)/(dz(i)*porsl(i))
               satw = min(1., satw)
               if(satw>.1e-6)then          
                  fl = wliq(i)/(wice(i)+wliq(i))
                  if(tss(i) >= tfrz) then       ! Unfrozen soil
                     dke = log10(satw) + 1.0
                     dke = max(dke, 0.)
                     dksat = dksatu(i)
                  else                          ! Frozen soil
                     dke = satw
                     dksat = dkmg(i)*0.249**(fl*porsl(i))*2.29**porsl(i)
                  end if
                  thk(i) = dke*dksat + (1.-dke)*dkdry(i)
               endif
            endif
         else                       ! wetland or glacier
            thk(i) = tkwat
            if(tss(i)<tfrz) thk(i) = tkice
         endif

      enddo

! Thermal conductivity of snow, which from Jordan (1991) pp. 18
      if(lb < 1)then
        do i = lb, 0
          rhosnow = (wice(i)+wliq(i))/dz(i)
          thk(i) = tkair+(7.75e-5*rhosnow+1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
        enddo
      endif

! Thermal conductivity at the layer interface
      do i = lb, nl_soil-1

! the following consideration is try to avoid the snow conductivity 
! to be dominant in the thermal conductivity of the interface. 
! Because when the distance of bottom snow node to the interfacee 
! is larger than that of interface to top soil node,
! the snow thermal conductivity will be dominant, and the result is that 
! lees heat tranfer between snow and soil 
         if((i==0) .AND. (z(i+1)-zi(i)<zi(i)-z(i)))then
            tk(i) = 2.*thk(i)*thk(i+1)/(thk(i)+thk(i+1))
            tk(i) = max(0.5*thk(i+1),tk(i))
         else
            tk(i) = thk(i)*thk(i+1)*(z(i+1)-z(i)) &
                  /(thk(i)*(z(i+1)-zi(i))+thk(i+1)*(zi(i)-z(i)))
         endif
      enddo
!=====liruichao add====
!for derivation of the tk of new interface for the 3 layers adjusted later
#if (defined FHNP) && (defined FTF)
      thkhc(lb:nl_soil) = 0.0_r8
      thkhc(lb:nl_soil)=thk(lb:nl_soil)
#endif
!=====liruichao add====
 end subroutine hConductivity

 subroutine hConductivity_lakesoil (itypwat,lb,nl_soil,&
          dkdry,dksatu,porsl,dz_soisno,z_soisno,zi_soisno,t_soisno,wice_soisno,wliq_soisno,tk,tktopsoil)
!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, September 15, 1999
!
! calculation of thermal conductivities of snow / soil layers
! The thermal conductivity of soil is computed from
! the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
! the formulation used in SNTHERM (Jordan 1991).
!
! The thermal conductivities at the interfaces between two neighbor layers
! (j, j+1) are derived from an assumption that the flux across the interface
! is equal to that from the node j to the interface and the flux from the
! interface to the node j+1.
!
! ________________
! REVISION HISTORY:
! 07/19/2014, Yongjiu Dai: treat the wetland as soil column instead of
!                          water body.
! 08/16/2014, Nan Wei: recalculate the heat conductivity of soil layers 
!                      underneath the lake
!-----------------------------------------------------------------------

  use precision
  use phycon_module, only : denh2o,denice,tfrz,tkwat,tkice,tkair
  implicit none

  integer, INTENT(in) :: lb       ! lower bound of array
  integer, INTENT(in) :: nl_soil  ! upper bound of array
  integer, INTENT(in) :: itypwat  ! land water type (0=soil, 1=urban,2=wetland,
                                  ! 3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) ::  dkdry(1:nl_soil)  ! thermal conductivity for dry soil [W/m-K]
  real(r8), INTENT(in) :: dksatu(1:nl_soil)  ! Thermal conductivity of saturated soil [W/m-K]
  real(r8), INTENT(in) ::  porsl(1:nl_soil)  ! fractional volume between soil grains=1.-dmvol
  real(r8), INTENT(in) ::   dz_soisno(lb:nl_soil)   ! layer thickiness [m]
  real(r8), INTENT(in) ::    z_soisno(lb:nl_soil)   ! node depth [m]
  real(r8), INTENT(in) ::   zi_soisno(lb-1:nl_soil) ! interface depth [m]
  real(r8), INTENT(in) ::  t_soisno(lb:nl_soil)   ! Nodal temperature [K]
  real(r8), INTENT(in) :: wice_soisno(lb:nl_soil)   ! ice lens [kg/m2]
  real(r8), INTENT(in) :: wliq_soisno(lb:nl_soil)   ! liqui water [kg/m2]

  real(r8), INTENT(out) :: tk(lb:nl_soil)    ! thermal conductivity [W/(m K)]
  real(r8), optional, INTENT(out) :: tktopsoil

! local
  real(r8) rhosnow  ! partitial density of water (ice + liquid)
  real(r8) dksat    ! thermal conductivity for saturated soil (j/(k s m))
  real(r8) dke      ! kersten number
  real(r8) fl       ! fraction of liquid or unfrozen water to total water
  real(r8) satw     ! relative total water content of soil.
  real(r8) thk(lb:nl_soil)  ! thermal conductivity of layer
  real(r8)  xicevol

  integer i

!-----------------------------------------------------------------------
! Thermal conductivity of soil from Farouki (1981),
      do i = 1, nl_soil
         thk(i) = dkdry(i)       !rock or dry soil
         if(porsl(i)>1.e-05 .and. (wice_soisno(i)+wliq_soisno(i)) > 0.0)then
            satw = 1.
            fl = wliq_soisno(i)/(wice_soisno(i)+wliq_soisno(i))
            if(t_soisno(i) >= tfrz) then  ! Unfrozen soil
               dke = log10(satw) + 1.0
               dke = max(dke, 0.)
               dksat = dksatu(i)
            else                          ! Frozen soil
               dke = satw
               dksat = dksatu(i)*(2.29/0.57)**((1.-fl)*porsl(i))
            end if                                                
            thk(i) = dke*dksat + (1.-dke)*dkdry(i)
            satw = (wliq_soisno(i)/denh2o+wice_soisno(i)/denice)/(dz_soisno(i)*porsl(i))
            if(satw > 1.0)then
               xicevol = (satw-1.0)*porsl(i)
               thk(i) = (thk(i) + xicevol*tkice)/(1.0 + xicevol)/(1.0 + xicevol)
            end if
            if(present(tktopsoil))tktopsoil = thk(1)
         end if
      enddo

! Thermal conductivity of snow
      if(lb < 1)then
         do i = lb, 0
            rhosnow = (wice_soisno(i)+wliq_soisno(i))/dz_soisno(i)

         ! presently option [1] is the default option
         ! [1] Jordan (1991) pp. 18
            thk(i) = tkair+(7.75e-5*rhosnow+1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
         end do
      end if

      do i = lb, nl_soil-1
         if (i /= 0) then
             tk(i) = thk(i)*thk(i+1)*(z_soisno(i+1)-z_soisno(i)) &
             /(thk(i)*(z_soisno(i+1)-zi_soisno(i))+thk(i+1)*(zi_soisno(i)-z_soisno(i)))
         else if (i == 0 .and. i>=lb) then
             tk(i) = thk(i)
         end if
      enddo
      tk(nl_soil) = 0.

 end subroutine hConductivity_lakesoil
