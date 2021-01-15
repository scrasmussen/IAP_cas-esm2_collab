      module mo_nln_matrix
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: nlnmat
      contains
      subroutine nlnmat01( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(225) = -(rxt(20)*y(2) + rxt(27)*y(3) + rxt(34)*y(18) + rxt(39)*y(19) &
                      + rxt(45)*y(20) + rxt(50)*y(6) + rxt(53)*y(7) + rxt(64)*y(22) &
                      + rxt(91)*y(30))
         mat(197) = -rxt(20)*y(1)
         mat(65) = -rxt(27)*y(1)
         mat(101) = -rxt(34)*y(1)
         mat(169) = -rxt(39)*y(1)
         mat(146) = -rxt(45)*y(1)
         mat(237) = -rxt(50)*y(1)
         mat(210) = -rxt(53)*y(1)
         mat(107) = -rxt(64)*y(1)
         mat(83) = -rxt(91)*y(1)
         mat(195) = -(rxt(20)*y(1) + 4._r8*rxt(21)*y(2) + rxt(38)*y(19) + rxt(44) &
                      *y(20) + rxt(47)*y(21) + rxt(48)*y(6) + (rxt(51) + rxt(52) &
                      ) * y(7) + rxt(58)*y(8) + rxt(67)*y(24) + rxt(83)*y(16))
         mat(223) = -rxt(20)*y(2)
         mat(167) = -rxt(38)*y(2)
         mat(144) = -rxt(44)*y(2)
         mat(56) = -rxt(47)*y(2)
         mat(235) = -rxt(48)*y(2)
         mat(208) = -(rxt(51) + rxt(52)) * y(2)
         mat(127) = -rxt(58)*y(2)
         mat(182) = -rxt(67)*y(2)
         mat(115) = -rxt(83)*y(2)
         mat(167) = mat(167) + 2.000_r8*rxt(41)*y(19)
         mat(100) = rxt(37)*y(20)
         mat(144) = mat(144) + rxt(37)*y(18)
         mat(57) = -(rxt(27)*y(1))
         mat(212) = -rxt(27)*y(3)
         mat(84) = -((rxt(79) + rxt(80)) * y(19))
         mat(158) = -(rxt(79) + rxt(80)) * y(17)
         mat(214) = .050_r8*rxt(91)*y(30)
         mat(187) = rxt(83)*y(16)
         mat(158) = mat(158) + rxt(82)*y(16)
         mat(121) = rxt(81)*y(16)
         mat(109) = rxt(83)*y(2) + rxt(82)*y(19) + rxt(81)*y(8)
         mat(78) = .050_r8*rxt(91)*y(1)
         mat(238) = -(rxt(48)*y(2) + rxt(49)*y(20) + rxt(50)*y(1) + rxt(57)*y(8) &
                      + rxt(71)*y(24) + rxt(84)*y(13))
         mat(198) = -rxt(48)*y(6)
         mat(147) = -rxt(49)*y(6)
         mat(226) = -rxt(50)*y(6)
         mat(130) = -rxt(57)*y(6)
         mat(185) = -rxt(71)*y(6)
         mat(96) = -rxt(84)*y(6)
         mat(198) = mat(198) + rxt(51)*y(7)
         mat(211) = rxt(51)*y(2)
         mat(209) = -((rxt(51) + rxt(52)) * y(2) + rxt(53)*y(1) + rxt(54)*y(8) + rxt(56) &
                      *y(19) + rxt(61)*y(20) + rxt(72)*y(24))
         mat(196) = -(rxt(51) + rxt(52)) * y(7)
         mat(224) = -rxt(53)*y(7)
         mat(128) = -rxt(54)*y(7)
         mat(168) = -rxt(56)*y(7)
         mat(145) = -rxt(61)*y(7)
         mat(183) = -rxt(72)*y(7)
         mat(224) = mat(224) + rxt(50)*y(6)
         mat(196) = mat(196) + rxt(48)*y(6) + rxt(58)*y(8)
         mat(236) = rxt(50)*y(1) + rxt(48)*y(2) + 2.000_r8*rxt(57)*y(8) + rxt(84) &
                      *y(13) + rxt(49)*y(20) + rxt(71)*y(24)
         mat(168) = mat(168) + rxt(59)*y(8) + rxt(62)*y(10)
         mat(128) = mat(128) + rxt(58)*y(2) + 2.000_r8*rxt(57)*y(6) + rxt(59)*y(19) &
                      + rxt(60)*y(20)
         mat(71) = rxt(62)*y(19)
         mat(95) = rxt(84)*y(6)
         mat(145) = mat(145) + rxt(49)*y(6) + rxt(60)*y(8)
         mat(183) = mat(183) + rxt(71)*y(6)
         mat(165) = -(rxt(38)*y(2) + rxt(39)*y(1) + rxt(40)*y(20) + (4._r8*rxt(41) &
                      + 4._r8*rxt(42)) * y(19) + rxt(43)*y(21) + rxt(56)*y(7) + rxt(59) &
                      *y(8) + rxt(62)*y(10) + (rxt(68) + rxt(69)) * y(24) + (rxt(79) &
                      + rxt(80)) * y(17) + rxt(82)*y(16) + rxt(88)*y(15) + rxt(89) &
                      *y(14) + rxt(90)*y(30) + rxt(95)*y(34) + (rxt(96) + rxt(97) &
                      ) * y(35))
         mat(193) = -rxt(38)*y(19)
         mat(221) = -rxt(39)*y(19)
         mat(142) = -rxt(40)*y(19)
         mat(55) = -rxt(43)*y(19)
         mat(206) = -rxt(56)*y(19)
         mat(126) = -rxt(59)*y(19)
         mat(70) = -rxt(62)*y(19)
         mat(180) = -(rxt(68) + rxt(69)) * y(19)
         mat(87) = -(rxt(79) + rxt(80)) * y(19)
         mat(114) = -rxt(82)*y(19)
         mat(47) = -rxt(88)*y(19)
         mat(76) = -rxt(89)*y(19)
         mat(82) = -rxt(90)*y(19)
         mat(30) = -rxt(95)*y(19)
         mat(43) = -(rxt(96) + rxt(97)) * y(19)
         mat(221) = mat(221) + rxt(34)*y(18) + rxt(45)*y(20)
         mat(193) = mat(193) + rxt(83)*y(16) + rxt(44)*y(20) + rxt(47)*y(21)
         mat(233) = rxt(49)*y(20)
         mat(165) = mat(165) + .300_r8*rxt(89)*y(14)
         mat(126) = mat(126) + rxt(60)*y(20)
         mat(76) = mat(76) + .300_r8*rxt(89)*y(19)
         mat(114) = mat(114) + rxt(83)*y(2)
         mat(99) = rxt(34)*y(1) + 2.000_r8*rxt(35)*y(20)
         mat(142) = mat(142) + rxt(45)*y(1) + rxt(44)*y(2) + rxt(49)*y(6) + rxt(60) &
                      *y(8) + 2.000_r8*rxt(35)*y(18) + rxt(66)*y(22)
         mat(55) = mat(55) + rxt(47)*y(2)
         mat(105) = rxt(66)*y(20)
         mat(124) = -(rxt(54)*y(7) + rxt(57)*y(6) + rxt(58)*y(2) + rxt(59)*y(19) &
                      + rxt(60)*y(20) + rxt(81)*y(16) + rxt(98)*y(35))
         mat(204) = -rxt(54)*y(8)
         mat(231) = -rxt(57)*y(8)
         mat(191) = -rxt(58)*y(8)
         mat(163) = -rxt(59)*y(8)
         mat(140) = -rxt(60)*y(8)
         mat(112) = -rxt(81)*y(8)
         mat(41) = -rxt(98)*y(8)
         mat(219) = rxt(53)*y(7)
         mat(191) = mat(191) + rxt(52)*y(7)
         mat(204) = mat(204) + rxt(53)*y(1) + rxt(52)*y(2)
         mat(199) = rxt(56)*y(19)
         mat(148) = rxt(56)*y(7)
         mat(117) = rxt(81)*y(16) + rxt(98)*y(35)
         mat(108) = rxt(81)*y(8)
         mat(38) = rxt(98)*y(8)
         mat(67) = -(rxt(62)*y(19))
         mat(155) = -rxt(62)*y(10)
         mat(202) = rxt(61)*y(20)
         mat(134) = rxt(61)*y(7)
         mat(200) = rxt(54)*y(8)
         mat(119) = rxt(54)*y(7)
         mat(90) = -(rxt(84)*y(6) + rxt(85)*y(20) + (4._r8*rxt(86) + 4._r8*rxt(87) &
                      ) * y(13))
         mat(227) = -rxt(84)*y(13)
         mat(136) = -rxt(85)*y(13)
         mat(215) = 1.860_r8*rxt(91)*y(30)
         mat(159) = .700_r8*rxt(89)*y(14)
         mat(73) = .700_r8*rxt(89)*y(19)
         mat(79) = 1.860_r8*rxt(91)*y(1)
         mat(72) = -(rxt(89)*y(19))
         mat(156) = -rxt(89)*y(14)
         mat(89) = rxt(85)*y(20)
         mat(135) = rxt(85)*y(13)
         mat(111) = -(rxt(81)*y(8) + rxt(82)*y(19) + rxt(83)*y(2))
         mat(123) = -rxt(81)*y(16)
         mat(162) = -rxt(82)*y(16)
         mat(190) = -rxt(83)*y(16)
         mat(218) = .870_r8*rxt(91)*y(30)
         mat(230) = rxt(84)*y(13)
         mat(162) = mat(162) + .300_r8*rxt(89)*y(14) + rxt(88)*y(15)
         mat(92) = rxt(84)*y(6) + (4.000_r8*rxt(86)+2.000_r8*rxt(87))*y(13)
         mat(75) = .300_r8*rxt(89)*y(19)
         mat(80) = .870_r8*rxt(91)*y(1)
         mat(45) = rxt(88)*y(19)
         mat(97) = -(rxt(34)*y(1) + (rxt(35) + rxt(36) + rxt(37)) * y(20))
         mat(216) = -rxt(34)*y(18)
         mat(137) = -(rxt(35) + rxt(36) + rxt(37)) * y(18)
         mat(188) = rxt(38)*y(19)
         mat(85) = rxt(79)*y(19)
         mat(160) = rxt(38)*y(2) + rxt(79)*y(17) + rxt(82)*y(16)
         mat(110) = rxt(82)*y(19)
         mat(141) = -((rxt(35) + rxt(36) + rxt(37)) * y(18) + rxt(40)*y(19) + rxt(44) &
                      *y(2) + rxt(45)*y(1) + 4._r8*rxt(46)*y(20) + rxt(49)*y(6) + rxt(60) &
                      *y(8) + rxt(61)*y(7) + (rxt(65) + rxt(66)) * y(22) + rxt(70) &
                      *y(24) + rxt(85)*y(13))
         mat(98) = -(rxt(35) + rxt(36) + rxt(37)) * y(20)
         mat(164) = -rxt(40)*y(20)
         mat(192) = -rxt(44)*y(20)
         mat(220) = -rxt(45)*y(20)
         mat(232) = -rxt(49)*y(20)
         mat(125) = -rxt(60)*y(20)
         mat(205) = -rxt(61)*y(20)
         mat(104) = -(rxt(65) + rxt(66)) * y(20)
         mat(179) = -rxt(70)*y(20)
         mat(93) = -rxt(85)*y(20)
         mat(220) = mat(220) + rxt(39)*y(19) + .060_r8*rxt(91)*y(30)
         mat(192) = mat(192) + rxt(83)*y(16) + rxt(47)*y(21)
         mat(86) = rxt(80)*y(19)
         mat(232) = mat(232) + rxt(84)*y(13)
         mat(164) = mat(164) + rxt(39)*y(1) + rxt(80)*y(17) + rxt(59)*y(8) + rxt(43) &
                      *y(21) + rxt(68)*y(24) + rxt(88)*y(15) + .500_r8*rxt(97)*y(35)
         mat(125) = mat(125) + rxt(59)*y(19) + rxt(81)*y(16)
         mat(93) = mat(93) + rxt(84)*y(6) + 4.000_r8*rxt(86)*y(13)
         mat(113) = rxt(83)*y(2) + rxt(81)*y(8)
         mat(54) = rxt(47)*y(2) + rxt(43)*y(19)
         mat(179) = mat(179) + rxt(68)*y(19)
         mat(81) = .060_r8*rxt(91)*y(1)
         mat(46) = rxt(88)*y(19)
         mat(42) = .500_r8*rxt(97)*y(19)
         mat(53) = -(rxt(43)*y(19) + rxt(47)*y(2))
         mat(154) = -rxt(43)*y(21)
         mat(186) = -rxt(47)*y(21)
         mat(154) = mat(154) + 2.000_r8*rxt(42)*y(19)
         mat(133) = 2.000_r8*rxt(46)*y(20)
         mat(103) = -(rxt(64)*y(1) + (rxt(65) + rxt(66)) * y(20))
         mat(217) = -rxt(64)*y(22)
         mat(138) = -(rxt(65) + rxt(66)) * y(22)
         mat(189) = rxt(67)*y(24)
         mat(229) = rxt(71)*y(24)
         mat(161) = rxt(68)*y(24)
         mat(177) = rxt(67)*y(2) + rxt(71)*y(6) + rxt(68)*y(19) + (4.000_r8*rxt(73) &
                       +2.000_r8*rxt(75))*y(24)
         mat(171) = 2.000_r8*rxt(74)*y(24)
         mat(181) = -(rxt(67)*y(2) + (rxt(68) + rxt(69)) * y(19) + rxt(70)*y(20) &
                      + rxt(71)*y(6) + rxt(72)*y(7) + (4._r8*rxt(73) + 4._r8*rxt(74) &
                      + 4._r8*rxt(75) + 4._r8*rxt(76)) * y(24))
         mat(194) = -rxt(67)*y(24)
         mat(166) = -(rxt(68) + rxt(69)) * y(24)
         mat(143) = -rxt(70)*y(24)
         mat(234) = -rxt(71)*y(24)
         mat(207) = -rxt(72)*y(24)
         mat(222) = rxt(64)*y(22)
         mat(143) = mat(143) + rxt(66)*y(22)
         mat(106) = rxt(64)*y(1) + rxt(66)*y(20)
      end subroutine nlnmat01
      subroutine nlnmat02( mat, y, rxt )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat(172) = 2.000_r8*rxt(75)*y(24)
         mat(174) = 2.000_r8*rxt(76)*y(24)
         mat(149) = rxt(69)*y(24)
         mat(131) = rxt(65)*y(22)
         mat(102) = rxt(65)*y(20)
         mat(173) = rxt(69)*y(19)
         mat(132) = rxt(70)*y(24)
         mat(175) = rxt(70)*y(20)
         mat(201) = rxt(72)*y(24)
         mat(176) = rxt(72)*y(7)
         mat(77) = -(rxt(90)*y(19) + rxt(91)*y(1))
         mat(157) = -rxt(90)*y(30)
         mat(213) = -rxt(91)*y(30)
         mat(44) = -(rxt(88)*y(19))
         mat(153) = -rxt(88)*y(15)
         mat(88) = 2.000_r8*rxt(87)*y(13)
         mat(29) = -(rxt(95)*y(19))
         mat(151) = -rxt(95)*y(34)
         mat(151) = mat(151) + (rxt(96)+.500_r8*rxt(97))*y(35)
         mat(118) = rxt(98)*y(35)
         mat(39) = (rxt(96)+.500_r8*rxt(97))*y(19) + rxt(98)*y(8)
         mat(40) = -((rxt(96) + rxt(97)) * y(19) + rxt(98)*y(8))
         mat(152) = -(rxt(96) + rxt(97)) * y(35)
         mat(120) = -rxt(98)*y(35)
         mat(150) = rxt(95)*y(34)
         mat(28) = rxt(95)*y(19)
      end subroutine nlnmat02
      subroutine nlnmat_finit( mat, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(inout) :: mat(nzcnt)
!----------------------------------------------
! ... local variables
!----------------------------------------------
!----------------------------------------------
! ... complete matrix entries implicit species
!----------------------------------------------
         mat( 1) = lmat( 1)
         mat( 2) = lmat( 2)
         mat( 3) = lmat( 3)
         mat( 4) = lmat( 4)
         mat( 5) = lmat( 5)
         mat( 6) = lmat( 6)
         mat( 7) = lmat( 7)
         mat( 8) = lmat( 8)
         mat( 9) = lmat( 9)
         mat( 10) = lmat( 10)
         mat( 11) = lmat( 11)
         mat( 12) = lmat( 12)
         mat( 13) = lmat( 13)
         mat( 14) = lmat( 14)
         mat( 15) = lmat( 15)
         mat( 16) = lmat( 16)
         mat( 17) = lmat( 17)
         mat( 18) = lmat( 18)
         mat( 19) = lmat( 19)
         mat( 20) = lmat( 20)
         mat( 21) = lmat( 21)
         mat( 22) = lmat( 22)
         mat( 23) = lmat( 23)
         mat( 24) = lmat( 24)
         mat( 25) = lmat( 25)
         mat( 26) = lmat( 26)
         mat( 27) = lmat( 27)
         mat( 29) = mat( 29) + lmat( 29)
         mat( 31) = lmat( 31)
         mat( 32) = lmat( 32)
         mat( 33) = lmat( 33)
         mat( 34) = lmat( 34)
         mat( 35) = lmat( 35)
         mat( 36) = lmat( 36)
         mat( 37) = lmat( 37)
         mat( 40) = mat( 40) + lmat( 40)
         mat( 44) = mat( 44) + lmat( 44)
         mat( 48) = lmat( 48)
         mat( 49) = lmat( 49)
         mat( 50) = lmat( 50)
         mat( 51) = lmat( 51)
         mat( 52) = lmat( 52)
         mat( 53) = mat( 53) + lmat( 53)
         mat( 55) = mat( 55) + lmat( 55)
         mat( 57) = mat( 57) + lmat( 57)
         mat( 58) = lmat( 58)
         mat( 59) = lmat( 59)
         mat( 60) = lmat( 60)
         mat( 61) = lmat( 61)
         mat( 62) = lmat( 62)
         mat( 63) = lmat( 63)
         mat( 64) = lmat( 64)
         mat( 66) = lmat( 66)
         mat( 67) = mat( 67) + lmat( 67)
         mat( 68) = lmat( 68)
         mat( 69) = lmat( 69)
         mat( 70) = mat( 70) + lmat( 70)
         mat( 71) = mat( 71) + lmat( 71)
         mat( 72) = mat( 72) + lmat( 72)
         mat( 74) = lmat( 74)
         mat( 75) = mat( 75) + lmat( 75)
         mat( 76) = mat( 76) + lmat( 76)
         mat( 77) = mat( 77) + lmat( 77)
         mat( 84) = mat( 84) + lmat( 84)
         mat( 90) = mat( 90) + lmat( 90)
         mat( 97) = mat( 97) + lmat( 97)
         mat( 98) = mat( 98) + lmat( 98)
         mat( 103) = mat( 103) + lmat( 103)
         mat( 109) = mat( 109) + lmat( 109)
         mat( 110) = mat( 110) + lmat( 110)
         mat( 111) = mat( 111) + lmat( 111)
         mat( 117) = mat( 117) + lmat( 117)
         mat( 124) = mat( 124) + lmat( 124)
         mat( 127) = mat( 127) + lmat( 127)
         mat( 128) = mat( 128) + lmat( 128)
         mat( 130) = mat( 130) + lmat( 130)
         mat( 133) = mat( 133) + lmat( 133)
         mat( 141) = mat( 141) + lmat( 141)
         mat( 159) = mat( 159) + lmat( 159)
         mat( 165) = mat( 165) + lmat( 165)
         mat( 181) = mat( 181) + lmat( 181)
         mat( 195) = mat( 195) + lmat( 195)
         mat( 197) = mat( 197) + lmat( 197)
         mat( 199) = mat( 199) + lmat( 199)
         mat( 206) = mat( 206) + lmat( 206)
         mat( 208) = mat( 208) + lmat( 208)
         mat( 209) = mat( 209) + lmat( 209)
         mat( 211) = mat( 211) + lmat( 211)
         mat( 212) = mat( 212) + lmat( 212)
         mat( 223) = mat( 223) + lmat( 223)
         mat( 225) = mat( 225) + lmat( 225)
         mat( 238) = mat( 238) + lmat( 238)
         mat( 91) = 0._r8
         mat( 94) = 0._r8
         mat( 116) = 0._r8
         mat( 122) = 0._r8
         mat( 129) = 0._r8
         mat( 139) = 0._r8
         mat( 170) = 0._r8
         mat( 178) = 0._r8
         mat( 184) = 0._r8
         mat( 203) = 0._r8
         mat( 228) = 0._r8
         mat( 1) = mat( 1) - dti
         mat( 2) = mat( 2) - dti
         mat( 3) = mat( 3) - dti
         mat( 4) = mat( 4) - dti
         mat( 5) = mat( 5) - dti
         mat( 6) = mat( 6) - dti
         mat( 7) = mat( 7) - dti
         mat( 8) = mat( 8) - dti
         mat( 9) = mat( 9) - dti
         mat( 10) = mat( 10) - dti
         mat( 11) = mat( 11) - dti
         mat( 12) = mat( 12) - dti
         mat( 13) = mat( 13) - dti
         mat( 14) = mat( 14) - dti
         mat( 15) = mat( 15) - dti
         mat( 16) = mat( 16) - dti
         mat( 17) = mat( 17) - dti
         mat( 18) = mat( 18) - dti
         mat( 19) = mat( 19) - dti
         mat( 20) = mat( 20) - dti
         mat( 21) = mat( 21) - dti
         mat( 22) = mat( 22) - dti
         mat( 23) = mat( 23) - dti
         mat( 24) = mat( 24) - dti
         mat( 25) = mat( 25) - dti
         mat( 26) = mat( 26) - dti
         mat( 29) = mat( 29) - dti
         mat( 32) = mat( 32) - dti
         mat( 35) = mat( 35) - dti
         mat( 40) = mat( 40) - dti
         mat( 44) = mat( 44) - dti
         mat( 48) = mat( 48) - dti
         mat( 53) = mat( 53) - dti
         mat( 57) = mat( 57) - dti
         mat( 67) = mat( 67) - dti
         mat( 72) = mat( 72) - dti
         mat( 77) = mat( 77) - dti
         mat( 84) = mat( 84) - dti
         mat( 90) = mat( 90) - dti
         mat( 97) = mat( 97) - dti
         mat( 103) = mat( 103) - dti
         mat( 111) = mat( 111) - dti
         mat( 124) = mat( 124) - dti
         mat( 141) = mat( 141) - dti
         mat( 165) = mat( 165) - dti
         mat( 181) = mat( 181) - dti
         mat( 195) = mat( 195) - dti
         mat( 209) = mat( 209) - dti
         mat( 225) = mat( 225) - dti
         mat( 238) = mat( 238) - dti
      end subroutine nlnmat_finit
      subroutine nlnmat( mat, y, rxt, lmat, dti )
      use chem_mods, only : gas_pcnst, rxntot, nzcnt
      implicit none
!----------------------------------------------
! ... dummy arguments
!----------------------------------------------
      real(r8), intent(in) :: dti
      real(r8), intent(in) :: lmat(nzcnt)
      real(r8), intent(in) :: y(gas_pcnst)
      real(r8), intent(in) :: rxt(rxntot)
      real(r8), intent(inout) :: mat(nzcnt)
      call nlnmat01( mat, y, rxt )
      call nlnmat02( mat, y, rxt )
      call nlnmat_finit( mat, lmat, dti )
      end subroutine nlnmat
      end module mo_nln_matrix
