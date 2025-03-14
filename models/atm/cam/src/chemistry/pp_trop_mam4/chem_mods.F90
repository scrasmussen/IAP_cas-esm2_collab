      module chem_mods
!--------------------------------------------------------------
! ... Basic chemistry parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod, only : r8 => shr_kind_r8
      implicit none
      save
      integer, parameter :: hetcnt = 22, & ! number of heterogeneous processes
                            phtcnt = 1, & ! number of photolysis reactions
                            rxntot = 7, & ! number of total reactions
                            gascnt = 6, & ! number of gas phase reactions
                            nabscol = 2, & ! number of absorbing column densities
                            gas_pcnst = 22, & ! number of "gas phase" species
                            nfs = 8, & ! number of "fixed" species
                            relcnt = 0, & ! number of relationship species
                            grpcnt = 0, & ! number of group members
                            nzcnt = 24, & ! number of non-zero matrix entries
                            extcnt = 7, & ! number of species with external forcing
                            clscnt1 = 0, & ! number of species in explicit class
                            clscnt2 = 0, & ! number of species in hov class
                            clscnt3 = 0, & ! number of species in ebi class
                            clscnt4 = 22, & ! number of species in implicit class
                            clscnt5 = 0, & ! number of species in rodas class
                            indexm = 1, & ! index of total atm density in invariant array
                            indexh2o = 4, & ! index of water vapor density
                            clsze = 1, & ! loop length for implicit chemistry
                            rxt_tag_cnt = 4, &
                            nslvd = 0
      integer :: clscnt(5) = 0
      integer :: cls_rxt_cnt(4,5) = 0
      integer :: clsmap(gas_pcnst,5) = 0
      integer :: permute(gas_pcnst,5) = 0
      integer :: diag_map(clscnt4) = 0
      real(r8) :: adv_mass(gas_pcnst) = 0._r8
      real(r8) :: fix_mass(max(1,nfs))
      integer, allocatable :: rxt_tag_map(:)
      real(r8), allocatable :: pht_alias_mult(:,:)
      character(len=16), allocatable :: rxt_tag_lst(:)
      character(len=16), allocatable :: pht_alias_lst(:,:)
      character(len=8) :: het_lst(max(1,hetcnt))
      character(len=8) :: inv_lst(max(1,nfs))
      character(len=8) :: extfrc_lst(max(1,extcnt))
      logical :: frc_from_dataset(max(1,extcnt))
      character(len=8) :: slvd_lst(max(1,nslvd))
      end module chem_mods
