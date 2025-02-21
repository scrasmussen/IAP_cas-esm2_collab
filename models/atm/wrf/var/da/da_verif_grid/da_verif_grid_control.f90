module da_verif_grid_control 
   !----------------------------------------------------------------------------
   ! History:
   !
   ! Author:   Syed RH Rizvi     NCAR/MMM       10/08/2007
   !
   ! Abstract: Main module for defining and initializing various constants, 
   !   namelist variables etc. 
   !
   !----------------------------------------------------------------------------

  implicit none

  integer, parameter                     :: max_3d_variables = 6
  integer, parameter                     :: max_2d_variables = 7
  integer, parameter                     :: num_vert_levels = 20 
  integer, parameter                     :: max_num_scores = 3
  integer, parameter                     :: max_num_exps = 10
  real, parameter                        :: missing=1.0E35
!-----------------------------------------------------------------
  integer, parameter                     :: namelist_unit = 7
  integer, parameter                     :: input_file_unit = 9
  integer, parameter                     :: first_unit = 11
  integer, parameter                     :: second_unit = 12
  integer, parameter                     :: time_series_unit = 20
  integer, parameter                     :: time_average_unit = 21
!-----------------------------------------------------------------
  character (len=512)                    :: profile_time_series_3d
  character (len=512)                    :: time_series_2d
  character (len=50)                     :: filename
!-----------------------------------------------------------------
  integer                                :: stime(6), etime(6)
  character (len=4)                      :: year
  character (len=2)                      :: month, day, hour
  character (len=19)                     :: hstart, hend, hdate
  character (len=10)                     :: date,pdate
!-----------------------------------------------------------------
  real, dimension( num_vert_levels )     :: vert_levels
  integer                                :: nx, ny, nz, number_of_levels 
!-----------------------------------------------------------------
  integer                                :: io_status
!-----------------------------------------------------------------
  logical                                :: debug1, debug2
!
  logical                                :: verify_its_own_analysis
  integer                                :: num_verifying_experiments
  integer                                :: verify_forecast_hour
  integer                                :: domain
  character (len=512)                    :: control_exp_dir
  character (len=512), dimension (max_num_exps) :: verif_dirs
  character (len=512), dimension (max_num_exps) :: out_dirs
  integer                                :: start_year, end_year
  integer                                :: start_month, end_month
  integer                                :: start_day, end_day
  integer                                :: start_hour, end_hour
  integer                                :: start_minutes=0, end_minutes=0
  integer                                :: start_seconds=0, end_seconds=0
  integer                                :: interval_hour
  integer                                :: num3dvar
  integer                                :: num2dvar
  character(len=20)                      :: var3d(max_3d_variables)
  character(len=20)                      :: var2d(max_2d_variables)
  integer                                :: num_scores
  character(len=20)                      :: score_names(max_num_scores)
  character (len=1)                      :: vertical_type
  character (len=20)                     :: verification_file_string

!-----------------------------------------------------------------

  namelist /control_main/ verify_its_own_analysis, num_verifying_experiments, &
                verification_file_string,control_exp_dir,verif_dirs, out_dirs, &
                verify_forecast_hour, domain, vertical_type
  namelist /control_times/ start_year, start_month, start_day, start_hour, &
                end_year, end_month, end_day, end_hour, interval_hour
  namelist /control_vars/ num3dvar, var3d, num2dvar, var2d


  data vert_levels  /1000.0, 925.0, 850.0, 700.0, 500.0, 400.0, 300.0, 250.0, &
                      200.0, 150.0, 100.0,  70.0,  50.0,  30.0,  20.0,  10.0, &
                        5.0,   3.0,   2.0,   1.0/
end module da_verif_grid_control 
