#include <define.h>

 subroutine meltf (lb, nl_soil, dtime, &
                   fact, brr, hs, dhsdT, &
                   tssbef, tss, wliq, wice, imelt, &
                   scv, snowdp, sm, xmf, porsl, phi0, bsw, dz)   !niu

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, September 15, 1999
!
! calculation of the phase change within snow and soil layers:
! 
! (1) check the conditions which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e., melting),
!     or layer temperature is less than the freezing point
!     and the liquid water mass is not equal to zero (i.e., freezing);
! (2) assess the rate of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point;
! (3) re-adjust the ice and liquid mass, and the layer temperature
!
!-----------------------------------------------------------------------

  use precision
  use debug
  use phycon_module, only : tfrz, hfus, grav  !niu
  implicit none

!-----------------------------------------------------------------------

   integer, INTENT(in) :: nl_soil             ! upper bound of array (i.e., soil layers)
   integer, INTENT(in) :: lb                  ! lower bound of array (i.e., snl +1)
  real(r8), INTENT(in) :: dtime               ! time step [second]
  real(r8), INTENT(in) :: tssbef(lb:nl_soil)  ! temperature at previous time step [K]
  real(r8), INTENT(in) :: brr   (lb:nl_soil)
  real(r8), INTENT(in) :: fact  (lb:nl_soil)  ! temporary variables
  real(r8), INTENT(in) :: hs                  ! net ground heat flux into the surface
  real(r8), INTENT(in) :: dhsdT               ! temperature derivative of "hs"

  real(r8), INTENT(in) :: porsl(1:nl_soil)    ! soil porosity [-]                           !niu
  real(r8), INTENT(in) :: phi0 (1:nl_soil)    ! soil water suction, negative potential [mm] !niu
  real(r8), INTENT(in) :: bsw  (1:nl_soil)    ! clapp and hornbereger "b" parameter [-]     !niu
  real(r8), INTENT(in) :: dz   (1:nl_soil)    ! layer thickiness [m]                        !niu

  real(r8), INTENT(inout) :: tss (lb:nl_soil) ! temperature at current time step [K]
  real(r8), INTENT(inout) :: wice(lb:nl_soil) ! ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq(lb:nl_soil) ! liquid water [kg/m2]
  real(r8), INTENT(inout) :: scv              ! snow mass [kg/m2]
  real(r8), INTENT(inout) :: snowdp           ! snow depth [m]

  real(r8), INTENT(out) :: sm                 ! rate of snowmelt [mm/s, kg/(m2 s)]
  real(r8), INTENT(out) :: xmf                ! total latent heat of phase change
   integer, INTENT(out) :: imelt(lb:nl_soil)  ! flag for melting or freezing [-]

! Local 
  real(r8) :: hm(lb:nl_soil)                  ! energy residual [W/m2]
  real(r8) :: xm(lb:nl_soil)                  ! metling or freezing within a time step [kg/m2]
  real(r8) :: heatr                           ! energy residual or loss after melting or freezing
  real(r8) :: temp1                           ! temporary variables [kg/m2]
  real(r8) :: temp2                           ! temporary variables [kg/m2]

  real(r8) :: smp
  real(r8) :: supercool(1:nl_soil)

  real(r8), dimension(lb:nl_soil) :: wmass0, wice0, wliq0
#if (defined FHNP) && (defined FTF)
  real(r8) :: tinc2(lb:nl_soil)
  !Jinbo Xie new tinc calculation following clm4.5
#endif
  real(r8) :: propor,tinc, we, scvold  
  integer j

!-----------------------------------------------------------------------

  sm = 0.
  xmf = 0.
  do j = lb, nl_soil
     imelt(j) = 0
     hm(j) = 0.
     xm(j) = 0.
     wice0(j) = wice(j)
     wliq0(j) = wliq(j)
     wmass0(j) = wice(j) + wliq(j)
#if (defined FHNP) && (defined FTF)
     tinc2=0.!Jinbo Xie new tinc calculation following clm4.5
#endif
  enddo


  scvold=scv
  we=0.
  if(lb<=0) we = sum(wice(lb:0)+wliq(lb:0))

  do j = lb, 0       !snow layers      !niu

     ! Melting identification
     ! if ice exists above melt point, melt some to liquid.

     if(wice(j) > 0. .AND. tss(j) > tfrz)then
#if (defined FHNP) && (defined FTF)
        tinc2(j)=tfrz-tss(j)!Jinbo Xie new tinc calculation following clm4.5
#endif
        imelt(j) = 1
        tss(j) = tfrz
     endif

     ! Freezing identification
     ! if liquid exists below melt point, freeze some to ice.

     if(wliq(j) > 0. .AND. tss(j) < tfrz) then
#if (defined FHNP) && (defined FTF)
        tinc2(j)=tfrz-tss(j)!Jinbo Xie new tinc calculation following clm4.5
#endif
        imelt(j) = 2
        tss(j) = tfrz
     endif
  enddo

!niu in

  do j = 1, nl_soil       ! soil layers

     ! Melting identification
     ! if ice exists above melt point, melt some to liquid.

     if(wice(j) > 0. .AND. tss(j) > tfrz)then
#if (defined FHNP) && (defined FTF)
        tinc2(j)=tfrz-tss(j)!Jinbo Xie new tinc calculation following clm4.5
#endif
        imelt(j) = 1
        tss(j) = tfrz
     endif

     ! supercooling water
     supercool(j) = 0.0
     if(tss(j) < tfrz) then
        smp = hfus * (tfrz-tss(j))/(grav*tss(j)) * 1000.     ! mm
        if (porsl(j) > 0.) then                              ! jidy
           supercool(j) = porsl(j)*(smp/phi0(j))**(-1.0/bsw(j))
        else
           supercool(j) = 0.                                 ! jidy
        end if
        supercool(j) = supercool(j)*dz(j)*1000.              ! mm
     end if

     ! Freezing identification
     ! if liquid exists below melt point, freeze some to ice.

     if(wliq(j) > supercool(j) .AND. tss(j) < tfrz) then
#if (defined FHNP) && (defined FTF)
        tinc2(j)=tfrz-tss(j)!Jinbo Xie new tinc calculation following clm4.5
#endif
        imelt(j) = 2
        tss(j) = tfrz
     endif
  enddo

!niu out


! If snow exists, but its thickness less than the critical value (0.01 m)

  if(lb == 1 .AND. scv > 0.)then
     if(tss(1) > tfrz)then
#if (defined FHNP) && (defined FTF)
        tinc2(1)=tfrz-tss(1)!Jinbo Xie new tinc calculation following clm4.5
#endif
        imelt(1) = 1
        tss(1) = tfrz
     endif
  endif




! Calculate the energy surplus and loss for melting and freezing

  do j = lb, nl_soil
     if(imelt(j) > 0)then
        tinc = tss(j)-tssbef(j)
        if(j > lb)then
#if (defined FHNP) && (defined FTF)
           hm(j) =    - tinc2(j)/fact(j)
           !Jinbo Xie new tinc calculation following clm4.5
#else
           hm(j) = brr(j) - tinc/fact(j) 
#endif
        else
#if (defined FHNP) && (defined FTF)
           hm(j) =      dhsdT*tinc2(j)  - tinc2(j)/fact(j)
           !Jinbo Xie new tinc calculation following clm4.5
#else
           hm(j) = hs + dhsdT*tinc + brr(j) - tinc/fact(j) 
#endif
        endif
     endif
  enddo



  do j = lb, nl_soil
     if(imelt(j) == 1 .AND. hm(j) < 0.) then
       hm(j) = 0.
       imelt(j) = 0
     endif

! this error was checked carefully, it results from the the computed error
! of "Tridiagonal-Matrix" in subroutine "thermal".

     if(imelt(j) == 2 .AND. hm(j) > 0.) then
       hm(j) = 0.
       imelt(j) = 0
     endif
  enddo

! The rate of melting and freezing

  do j = lb, nl_soil

     if(imelt(j) > 0 .AND. abs(hm(j)) > .0) then

        xm(j) = hm(j)*dtime/hfus                        ! kg/m2

        ! if snow exists, but its thickness less than the critical value (1 cm)
        ! Note: more work is need on how to tune the snow depth at this case
        if(j == 1)then
        if ((lb == 1) .AND. (scv > 0.) .AND. (xm(j) > 0.))then
           temp1 = scv                                 ! kg/m2
           scv = max(0.,temp1-xm(j))
           propor = scv/temp1
           snowdp = propor * snowdp
           heatr = hm(j) - hfus*(temp1-scv)/dtime       ! W/m2
           if(heatr > 0.) then
              xm(j) = heatr*dtime/hfus                  ! kg/m2
              hm(j) = heatr                            ! W/m2
           else
              xm(j) = 0.
              hm(j) = 0.
           endif
           sm = max(0.,(temp1-scv))/dtime              ! kg/(m2 s)
           xmf = hfus*sm
        endif
        endif

        heatr = 0.
        if(xm(j) > 0.) then
           wice(j) = max(0., wice0(j)-xm(j))
           heatr = hm(j) - hfus*(wice0(j)-wice(j))/dtime
!niu    else if(xm(j) < 0.) then
        else if(xm(j) <= 0.) then
!niu       wice(j) = min(wmass0(j), wice0(j)-xm(j))
           if(j <= 0) then  ! snow                                        !niu
               wice(j) = min(wmass0(j), wice0(j)-xm(j))                   !niu
           else             ! soil                                        !niu
               if(wmass0(j) < supercool(j)) then                          !niu
                   wice(j) = 0.                                           !niu
               else                                                       !niu
                   wice(j) = min(wmass0(j)-supercool(j), wice0(j)-xm(j))  !niu
               endif                                                      !niu
           endif
           heatr = hm(j) - hfus*(wice0(j)-wice(j))/dtime  
        endif

        wliq(j) = max(0.,wmass0(j)-wice(j))

        if(abs(heatr) > 0.)then
           if(j > lb)then
              tss(j) = tss(j) + fact(j)*heatr
           else
              tss(j) = tss(j) + fact(j)*heatr/(1.-fact(j)*dhsdT)
           endif
!niu       if(wliq(j)*wice(j) > 0.) tss(j) = tfrz
           if(j <= 0) then ! snow                                         !niu
               if(wliq(j)*wice(j) > 0.) tss(j) = tfrz                     !niu
           end if                                                         !niu
        endif

        xmf = xmf + hfus * (wice0(j)-wice(j))/dtime

        if(imelt(j) == 1 .AND. j < 1) &
        sm = sm + max(0.,(wice0(j)-wice(j)))/dtime  

     endif

  enddo

  !scvold=scv
  if(lb<=0) then
  we = sum(wice(lb:0)+wliq(lb:0))-we
  if(abs(we)>1.e-6) stop 'meltf error'
  endif

 end subroutine meltf
