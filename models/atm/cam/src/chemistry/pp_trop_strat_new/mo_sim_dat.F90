      module mo_sim_dat
      private
      public :: set_sim_dat
      contains
      subroutine set_sim_dat
      use chem_mods, only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass
      use chem_mods, only : diag_map
      use chem_mods, only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods, only : pht_alias_lst, pht_alias_mult
      use chem_mods, only : het_lst, extfrc_lst, inv_lst, slvd_lst
      use abortutils, only : endrun
      use mo_tracname, only : solsym
      use chem_mods, only : frc_from_dataset
      use shr_kind_mod,only : r8 => shr_kind_r8
      use cam_logfile, only : iulog
      implicit none
!--------------------------------------------------------------
! ... local variables
!--------------------------------------------------------------
      integer :: ios
      clscnt(:) = (/ 5, 0, 0, 50, 0 /)
      cls_rxt_cnt(:,1) = (/ 2, 11, 0, 5 /)
      cls_rxt_cnt(:,4) = (/ 4, 34, 61, 50 /)
      solsym(: 55) = (/ 'O3      ','O       ','O1D     ','N2O     ','N       ', &
                        'NO      ','NO2     ','NO3     ','HNO3    ','HO2NO2  ', &
                        'N2O5    ','CH4     ','CH3O2   ','CH3OOH  ','CH3OH   ', &
                        'CH2O    ','CO      ','H       ','OH      ','HO2     ', &
                        'H2O2    ','CL      ','CL2     ','CLO     ','OCLO    ', &
                        'CL2O2   ','HCL     ','HOCL    ','CLONO2  ','ISOP    ', &
                        'CFC11   ','CFC12   ','CO2     ','SO2     ','DMS     ', &
                        'H2SO4   ','SOAG    ','SOA     ','so4_a1  ','pom_a1  ', &
                        'soa_a1  ','bc_a1   ','dst_a1  ','ncl_a1  ','num_a1  ', &
                        'so4_a2  ','soa_a2  ','ncl_a2  ','num_a2  ','dst_a3  ', &
                        'ncl_a3  ','so4_a3  ','num_a3  ','ncl_a4  ','num_a4  ' /)
      adv_mass(: 55) = (/ 47.99820_r8, 15.99940_r8, 15.99940_r8, 44.01288_r8, 14.00674_r8, &
                          30.00614_r8, 46.00554_r8, 62.00494_r8, 63.01234_r8, 79.01174_r8, &
                          108.0105_r8, 16.04060_r8, 47.03200_r8, 48.03940_r8, 32.04000_r8, &
                          30.02520_r8, 28.01040_r8, 1.007400_r8, 17.00680_r8, 33.00620_r8, &
                          34.01360_r8, 35.45270_r8, 70.90540_r8, 51.45210_r8, 67.45150_r8, &
                          102.9042_r8, 36.46010_r8, 52.45950_r8, 97.45764_r8, 68.11420_r8, &
                          137.3675_r8, 120.9132_r8, 44.00980_r8, 64.06480_r8, 62.13240_r8, &
                          98.07840_r8, 12.01100_r8, 144.1320_r8, 115.1073_r8, 12.01100_r8, &
                          12.01100_r8, 12.01100_r8, 135.0640_r8, 58.44247_r8, 1.007400_r8, &
                          115.1073_r8, 12.01100_r8, 58.44247_r8, 1.007400_r8, 135.0640_r8, &
                          58.44247_r8, 115.1073_r8, 1.007400_r8, 58.44247_r8, 1.007400_r8 /)
      fix_mass(: 4) = (/ 0.00000000_r8, 28.0134792_r8, 31.9988003_r8, 18.0142002_r8 /)
      clsmap(: 5,1) = (/ 12, 4, 31, 32, 33 /)
      clsmap(: 50,4) = (/ 1, 2, 3, 17, 5, 6, 7, 19, 8, 9, &
                            10, 11, 13, 14, 16, 18, 20, 21, 22, 23, &
                            24, 25, 26, 27, 28, 29, 30, 15, 34, 35, &
                            38, 36, 37, 39, 40, 41, 42, 43, 44, 45, &
                            46, 47, 48, 49, 50, 51, 52, 53, 54, 55 /)
      permute(: 50,4) = (/ 49, 47, 34, 38, 1, 50, 48, 45, 43, 2, &
                             35, 28, 39, 36, 42, 40, 44, 33, 41, 3, &
                             46, 4, 26, 5, 29, 32, 37, 31, 27, 30, &
                              6, 7, 8, 9, 10, 11, 12, 13, 14, 15, &
                             16, 17, 18, 19, 20, 21, 22, 23, 24, 25 /)
      diag_map(: 50) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
                            11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                            21, 22, 23, 24, 25, 26, 29, 32, 35, 40, &
                            44, 48, 53, 57, 67, 72, 77, 84, 90, 97, &
                           103, 111, 124, 141, 165, 181, 195, 209, 225, 238 /)
      het_lst(: 55) = (/ 'O3      ','O       ','O1D     ','N2O     ','N       ', &
                         'NO      ','NO2     ','NO3     ','HNO3    ','HO2NO2  ', &
                         'N2O5    ','CH4     ','CH3O2   ','CH3OOH  ','CH3OH   ', &
                         'CH2O    ','CO      ','H       ','OH      ','HO2     ', &
                         'H2O2    ','CL      ','CL2     ','CLO     ','OCLO    ', &
                         'CL2O2   ','HCL     ','HOCL    ','CLONO2  ','ISOP    ', &
                         'CFC11   ','CFC12   ','CO2     ','SO2     ','DMS     ', &
                         'H2SO4   ','SOAG    ','SOA     ','so4_a1  ','pom_a1  ', &
                         'soa_a1  ','bc_a1   ','dst_a1  ','ncl_a1  ','num_a1  ', &
                         'so4_a2  ','soa_a2  ','ncl_a2  ','num_a2  ','dst_a3  ', &
                         'ncl_a3  ','so4_a3  ','num_a3  ','ncl_a4  ','num_a4  ' /)
      extfrc_lst(: 10) = (/ 'NO      ','NO2     ','CO      ','SO2     ','so4_a1  ', &
                            'so4_a2  ','pom_a1  ','bc_a1   ','num_a1  ','num_a2  ' /)
      frc_from_dataset(: 10) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true., .true. /)
      inv_lst(: 4) = (/ 'M       ', 'N2      ', 'O2      ', 'H2O     ' /)
      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(:rxt_tag_cnt) = (/ 'jo2_b           ', 'jo3_a           ', 'jo3_b           ', 'jn2o            ', &
                                     'jno2            ', 'jno3_a          ', 'jno3_b          ', 'jho2no2_a       ', &
                                     'jho2no2_b       ', 'jch3ooh         ', 'jch2o_a         ', 'jch2o_b         ', &
                                     'jh2o2           ', 'jhocl           ', 'jclono2_a       ', 'jclono2_b       ', &
                                     'jcfcl3          ', 'jcf2cl2         ', 'usr_O_O2        ', 'O_O3            ', &
                                     'usr_O_O         ', 'O1D_N2          ', 'O1D_O2b         ', 'ox_l1           ', &
                                     'O1D_N2Oa        ', 'O1D_N2Ob        ', 'O1D_O3          ', 'O1D_CFC11       ', &
                                     'O1D_CFC12       ', 'O1D_CH4a        ', 'O1D_CH4b        ', 'O1D_CH4c        ', &
                                     'H_O2            ', 'H_O3            ', 'H_HO2a          ', 'H_HO2b          ', &
                                     'H_HO2c          ', 'OH_O            ', 'ox_l2           ', 'OH_HO2          ', &
                                     'OH_OH           ', 'OH_OH_M         ', 'OH_H2O2         ', 'HO2_O           ', &
                                     'ox_l3           ', 'usr_HO2_HO2     ', 'H2O2_O          ', 'NO_O_M          ', &
                                     'ox_p1           ', 'NO_O3           ', 'NO2_O           ', 'NO2_O_M         ', &
                                     'NO2_O3          ', 'tag_NO2_NO3     ', 'usr_N2O5_M      ', 'tag_NO2_OH      ', &
                                     'NO3_NO          ', 'NO3_O           ', 'NO3_OH          ', 'NO3_HO2         ', &
                                     'tag_NO2_HO2     ', 'HO2NO2_OH       ', 'usr_HO2NO2_M    ', 'CL_O3           ', &
                                     'CL_HO2a         ', 'CL_HO2b         ', 'CLO_O           ', 'CLO_OHa         ', &
                                     'CLO_OHb         ', 'CLO_HO2         ', 'CLO_NO          ', 'CLO_NO2_M       ', &
                                     'CLO_CLOa        ', 'CLO_CLOb        ', 'CLO_CLOc        ', 'tag_CLO_CLO_M   ', &
                                     'usr_CL2O2_M     ', 'CH4_OH          ', 'usr_CO_OH_b     ', 'CO_OH_M         ', &
                                     'CH2O_NO3        ', 'CH2O_OH         ', 'CH2O_O          ', 'ox_p2           ', &
                                     'CH3O2_HO2       ', 'CH3O2_CH3O2a    ', 'CH3O2_CH3O2b    ', 'CH3OH_OH        ', &
                                     'CH3OOH_OH       ', 'usr_N2O5_aer    ', 'usr_NO3_aer     ', 'usr_NO2_aer     ', &
                                     'usr_SO2_OH      ', 'DMS_OHb         ', 'usr_DMS_OH      ', 'DMS_NO3         ', &
                                     'usr_HO2_aer     ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, &
                                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, &
                                       31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                                       41, 42, 43, 44, 45, 46, 47, 48, 49, 50, &
                                       51, 52, 53, 54, 55, 56, 57, 58, 59, 60, &
                                       61, 62, 63, 64, 65, 66, 67, 68, 69, 70, &
                                       71, 72, 73, 74, 75, 76, 77, 78, 79, 80, &
                                       81, 82, 83, 84, 85, 86, 87, 88, 89, 92, &
                                       93, 94, 95, 96, 97, 98, 99 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ 'userdefined     ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_lst(:,2) = (/ '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ', '                ', '                ', &
                              '                ', '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8, 1._r8, 1._r8, &
                          1._r8, 1._r8, 1._r8 /)
      end subroutine set_sim_dat
      end module mo_sim_dat
