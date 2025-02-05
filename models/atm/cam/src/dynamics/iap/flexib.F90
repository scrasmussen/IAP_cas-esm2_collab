module flexib
!------------------------------------------------------------------------------------------------
! Purpose:    Set constants & variables which is flexible or alterable 
! Author :    ZhangHe
! Completed : 2005.9.1
! Update    : 2006.12.9, Zhanghe, set three kinds of vertical levels;
!            2006.12.30, ZhangHe, use 'parame2.inc' instead of statement in this module
!            2007.04.16, ZhangHe, about 'DTHDFS'
!            2007.04.20, ZhangHe, added 'NPHF', 'NENG' & 'NFST'
!            2007.05.23, ZhangHe, update 'PHALFC26' according to initial data of CAM3.1
!            2007.12.22, ZhangHe, reset DTIMEQ
!            2008.04.07, ZhangHe, added statement of Kstd
!            2008.04.10, ZhangHe, added statement of Istar
!            2010.08.10, ZhangHe, added Ndq and modify DTIMEQ
!            2011.05.03, ZhangHe, let Kstd = 0 if adiabatic run
!            2011.11.14, ZhangHe, added PHALFC30
!            2011.12.17, ZhangHe, added PHALFC(K) = PHALFC30(K)
!            2012.10.25, Jiang Jinrong, removed ksa to Dyn_const.F90
!            2013.03.29, ZhangHe, let NENG = 0
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,     only: NX, NY, NL, NZ

   implicit none
   save
   public

!begin czy20180506==========================================================
!czy20180327   integer , parameter :: Ndt    = 3                                              !wjp 2007.06
!   integer , parameter :: Ndt    = 1                                              !wjp 2007.06 czy20180327
!czy20180327   integer , parameter :: Ndq    = 3                                              !zhh 2010.08.10
!   integer , parameter :: Ndq    = 1                                              !zhh 2010.08.10 czy20180327
!   real(r8), parameter :: DTDY   = 200.0D0    ! timestep of dynamical integration !wjp 2007.06
!   real(r8), parameter :: DTDY   = 100.0D0    ! timestep of dynamical integration !czy 2017.08.19
!   real(r8), parameter :: DTDY   = 50.0D0    ! timestep of dynamical integration !czy 2018.03.26
!   real(r8), parameter :: DTlin  = DTDY       ! timestep of linear term
!   real(r8), parameter :: DTadv  = Ndt*DTlin  ! timestep of advection term
!   real(r8), parameter :: DTIMEQ = Ndq*DTlin   ! timestep of Q
   integer  :: Ndt            !wjp 2007.06 czy20180327
   integer  :: Ndq
   real(r8) :: DTDY
   real(r8) :: DTlin
   real(r8) :: DTadv
   real(r8) :: DTIMEQ
!end czy20180506==========================================================
   real(r8) :: DTHDFS         ! timestep of diffusion
   real(r8) :: PHALFC(NZ)     ! initialized pressure at sigma level 
   real(r8) :: PHALFC19(20)   ! initialized pressure at sigma level 
   real(r8) :: PHALFC26(27)   ! initialized pressure at sigma level 
   real(r8) :: PHALFC30(31)   ! initialized pressure at sigma level 
   real(r8) :: PHALFC31(32)   ! initialized pressure at sigma level 
   !czy20170715
   real(r8) :: PHALFC51(52)   ! initialized pressure at sigma level 
   !czy20180627
   real(r8) :: PHALFC35(36)   ! initialized pressure at sigma level 
!jjr   real(r8) :: ksa(NX,NY)     ! coefficient of diffusion of P'sa
   real(r8) :: EPS
   real(r8) :: a11, a12, a13  ! flexible coefficients of L(VT), for splitting-up algorithms
   real(r8) :: a21, a22, a23  ! flexible coefficients of L(UT)
   real(r8) :: a31, a32, a33  ! flexible coefficients of L(TT)
   real(r8) :: b1             ! flexible coefficients of [PX1,PY1,OO1]
   real(r8) :: b21, b22       ! flexible coefficients of [PY2,OO2Y] & [PX2,OO2X]
   real(r8) :: c10            ! flexible coefficients of [FSV, FSU]
   real(r8) :: BETA(3)        ! parameter of SHAPIRO SMOOTHER
   real(r8) :: DFS0           ! Karman conts which controled the coefficient of horizon diffusion
   integer  :: kappa          ! denotative constant for computing Dpsa
   integer  :: kappa2         ! for boundary condition coupled to ocean model
   integer  :: IDSOSS         ! switch of filter
   integer  :: IBCFFT         ! switch of FFT/shapiro smooth
   integer  :: NBCFFT         ! switch of FFT/shapiro smooth for physics
   integer  :: NPHF           ! switch of filter for physics
   integer  :: NENG           ! switch to compute energy budget 
   integer  :: NFST           ! switch to do smoothing 
   integer  :: Kstd           ! switch of standard stratification approximation
   integer  :: Istar          ! index of the methods to defeine flexible substitute
!zmh 2018.08.25 new flags
   logical  :: zmh_cfl_drag  = .false.          ! these two can be set in namelist
   real(r8) :: zmh_cfl_u0r   = 1.0_r8  !u interval m/s

   logical  :: zmh_nlitti1 = .false. !false. !true. !true.  combine forth and back integration in tend_lin
   logical  :: zmh_nlitti2 = .false. !false. !true.  combine forth and back integration in tend_adv
   logical  :: zmh_Coriolis= .false. !false.   ! Modified Coriolis algorithm
   !logical  :: zmh_tend1   = .false.  ! embedded substepping in tend_lin
   !logical  :: zmh_tend2   = .false.  ! embedded substepping in tend_adv
!wbd 2018.11.16 new flags
   logical  :: wbd_nlitti1 = .false. !true.  added the leaping format in tend_lin
   logical  :: wbd_nlitti2 = .false. !true.  added the leaping format in tend_adv


!  Set the data of the constants
! ------------------------------- according to ECWMF ---------------------------------
   data PHALFC19 / 0.000000, 20.00000, 40.00000, 60.83326, 86.24240       &   
                 , 119.7071, 163.6012, 219.2266, 286.8024, 365.4415       &
                 , 453.1477, 546.8523, 642.5411, 735.4829, 820.6007       &
                 , 893.0138, 948.7813, 985.8773, 1005.429, 1013.250 / 
! -------------------------------- according to CAM3 ----------------------------------
   data PHALFC26 / 2.194067,  4.895209,  9.882418,  18.05201,  29.83724   &    
                 , 44.62334,  61.60587,  78.51243,  92.36580,  108.66359  & 
                 , 127.83708, 150.39371, 176.93043, 208.14944, 244.87709  &                                                       
                 , 288.08522, 338.91731, 398.71865, 469.07180, 551.83871  &  														 
                 , 649.20969, 744.38289, 831.02123, 903.30029, 955.99746  &												
     	         , 985.11220, 1000.0 /
! -------------------------------- according to CAM5 ----------------------------------
   data PHALFC30 / 2.255240,  5.031692,  10.15795,  18.55532,  30.66912   &    
                 , 45.86748,  63.32348,  80.70142,  94.94104,  111.69321  & 
                 , 131.40127, 154.58681, 181.86335, 213.95282, 251.70442  &                                                       
                 , 296.11722, 348.36659, 409.83522, 482.14993, 567.22442  &  														 
                 , 652.33297, 730.44589, 796.36307, 845.35367, 873.71587  &												
     	         , 900.32463, 924.96446, 947.43233, 967.53862, 985.11219  &
                 , 1000.0   /
! ------------------------------- according to ECWMF ---------------------------------
   data PHALFC31 / 0.000000, 20.00000, 40.00000, 60.00000, 80.00000       &   
                 , 100.1574, 121.1638, 143.6299, 167.9520, 194.3569       &
                 , 222.9420, 253.7096, 286.5963, 321.4966, 358.2815       &
                 , 396.8118, 436.9459, 478.5425, 521.4575, 565.5346       &
                 , 610.6005, 656.4283, 702.7315, 749.1255, 795.0949       &
                 , 839.9533, 882.7980, 922.4593, 957.4447, 985.8772       &
                 , 1005.429, 1013.250 / 
! ------------------------------- according to ZhangMH & ZhangH ---------------------------------
!czy20170715 unit: hPa
   data PHALFC51 / 0.5000000,   1.0000000,   2.4079030,  6.1620226,   14.0680704        &
                , 25.7996941,  38.5396194,  50.9282112, 62.9261475,   74.2591858        &
                , 85.0522766,  95.7094650, 106.5234909                                  & 
                ,117.8194809, 129.9549408, 143.3197632, 158.3362427, 175.3836823        &
                ,194.7931061, 216.6529388, 241.0215454, 267.9503479, 297.4174194        & 
                ,329.3274231, 364.0086975, 401.4205933, 441.3315735, 483.3193054        &
                ,526.7706299, 570.8815308, 615.0315552, 658.2039185, 699.5576782        &
                ,738.0511475, 772.7767334, 802.8096924, 828.1500244, 849.7794189        &
                ,868.5499878, 885.1893921, 900.3826904, 914.4606323, 927.1274414        &
                ,938.5086060, 948.7559814, 958.0477295, 966.5883789, 974.5658569        &
                ,981.6986694, 988.2683716, 994.3688354, 1000.000000 /

! ------------------------------- according to ZhangMH & ZhangH ---------------------------------
!czy20180627 unit: hPa
   data PHALFC35 /2.25523950000000e+000, 5.03169190000000e+000, 10.1579470000000e+000, &
                18.5553170000000e+000, 30.6691230000000e+000, 45.8674770000000e+000, &
                63.3234830000000e+000, 80.7014180000000e+000, 94.9410420000000e+000, &
                111.693210000000e+000, 131.401270000000e+000, 154.586810000000e+000, &
                181.863350000000e+000, 213.952820000000e+000, 250.204930000000e+000, &
                291.983140000000e+000, 338.018650000000e+000, 392.278030000000e+000, &
                454.250020000000e+000, 521.832790000000e+000, 598.035750000000e+000, &
                673.895580000000e+000, 744.750690000000e+000, 787.179780000000e+000, &
                816.636130000000e+000, 841.753390000000e+000, 863.740960000000e+000, &
                883.406670000000e+000, 902.194170000000e+000, 920.004150000000e+000, &
                936.776980000000e+000, 952.486880000000e+000, 967.718710000000e+000, &
                981.766660000000e+000, 994.271540000000e+000, 1.00000000000000e+003 /
!
   data EPS    / 0.010D0 /
   data a11, a12, a13 / 1.0, 1.0, 1.0 /
   data a21, a22, a23 / 1.0, 1.0, 1.0 /
   data a31, a32, a33 / 1.0, 1.0, 1.0 /
   data b1            / 1.0 /
   data c10           / 1.0 /
   data b21, b22      / 1.0, 1.0 /
   data BETA   / 0.02D0, 0.05D0, 0.1D0 /
   data DFS0   / 0.1D0 /
   data kappa  / 1 /
   data kappa2 / 1 /
   data IDSOSS / 0 /
   data IBCFFT / 1 /
   data NBCFFT / 1 /
   data NPHF   / 1 /
   data NENG   / 0 /
   data NFST   / 1 /
   data Kstd   / 1 /   ! Kstd = 0; use standard stratification approximation
   data Istar  / 1 /      

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   subroutine setfle
!-----------------------------------------------------------------------------------------------
!    Purpose:  Set some consts & variables in module flexib
!    Author :    ZhangHe
!    Completed : 2005.9.20 
!    Update    : 2007.8.25, Zhanghe, change sub. name from 'dytest' to 'setfle'
!                2007.11.11, Zhanghe, add set PHALFC
!-----------------------------------------------------------------------------------------------
      use pmgrid,     only: beglatdyn, endlatdyn
      use cam_control_mod, only: adiabatic        !zhh, 2011-11-18
      implicit none
!----------------------------------Local workspace----------------------------------------------
      integer :: i,j,k
!-----------------------------------------------------------------------------------------------
      do K = 1, NZ                       
         if ( NL == 19 ) then
            PHALFC(K) = PHALFC19(K)
         else if ( NL == 26 ) then
            PHALFC(K) = PHALFC26(K)
         else if ( NL == 30 ) then
            PHALFC(K) = PHALFC30(K)
         else if ( NL == 31 ) then
            PHALFC(K) = PHALFC31(K)
            !czy20170715
         else if ( NL == 51 ) then
            PHALFC(K) = PHALFC51(K)
            !czy20180625
         else if ( NL == 35 ) then
            PHALFC(K) = PHALFC35(K)
         else
            print*, 'the value of NLAY is wrong'
         end if
      end do
!
      if (adiabatic) then
         Kstd = 0 
         kappa2 = 0
         kappa = 0
      end if

      return
   end subroutine

end module
