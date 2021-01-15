      module mo_prod_loss
      use shr_kind_mod, only : r8 => shr_kind_r8
      private
      public :: exp_prod_loss
      public :: imp_prod_loss
      contains
      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )
      use ppgrid, only : pver
      implicit none
!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:,:,:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:,:,:)
      real(r8), intent(in) :: rxt(:,:,:)
      real(r8), intent(in) :: het_rates(:,:,:)
!--------------------------------------------------------------------
! ... loss and production for Explicit method
!--------------------------------------------------------------------
         loss(:,:,1) = ((rxt(:,:,30) +rxt(:,:,31) +rxt(:,:,32))* y(:,:,3) +rxt(:,:,78) &
                 * y(:,:,19) + het_rates(:,:,12))* y(:,:,12)
         prod(:,:,1) = 0._r8
         loss(:,:,2) = ((rxt(:,:,25) +rxt(:,:,26))* y(:,:,3) + rxt(:,:,4) &
                  + het_rates(:,:,4))* y(:,:,4)
         prod(:,:,2) = 0._r8
         loss(:,:,3) = (rxt(:,:,28)* y(:,:,3) + rxt(:,:,17) + het_rates(:,:,31)) &
                 * y(:,:,31)
         prod(:,:,3) = 0._r8
         loss(:,:,4) = (rxt(:,:,29)* y(:,:,3) + rxt(:,:,18) + het_rates(:,:,32)) &
                 * y(:,:,32)
         prod(:,:,4) = 0._r8
         loss(:,:,5) = ( + het_rates(:,:,33))* y(:,:,33)
         prod(:,:,5) = 0._r8
      end subroutine exp_prod_loss
      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )
      use ppgrid, only : pver
      implicit none
!--------------------------------------------------------------------
! ... dummy args
!--------------------------------------------------------------------
      real(r8), dimension(:), intent(out) :: &
            prod, &
            loss
      real(r8), intent(in) :: y(:)
      real(r8), intent(in) :: rxt(:)
      real(r8), intent(in) :: het_rates(:)
!--------------------------------------------------------------------
! ... loss and production for Implicit method
!--------------------------------------------------------------------
         loss(49) = (rxt(20)* y(2) +rxt(27)* y(3) +rxt(50)* y(6) +rxt(53)* y(7) &
                  +rxt(34)* y(18) +rxt(39)* y(19) +rxt(45)* y(20) +rxt(64)* y(22) &
                  +rxt(91)* y(30) + rxt(2) + rxt(3) + het_rates(1))* y(1)
         prod(49) =rxt(19)*y(2)
         loss(47) = (rxt(20)* y(1) + 2._r8*rxt(21)* y(2) +rxt(48)* y(6) + (rxt(51) + &
                 rxt(52))* y(7) +rxt(58)* y(8) +rxt(83)* y(16) +rxt(38)* y(19) &
                  +rxt(44)* y(20) +rxt(47)* y(21) +rxt(67)* y(24) + rxt(19) &
                  + het_rates(2))* y(2)
         prod(47) = (rxt(22) +rxt(23))*y(3) +rxt(3)*y(1) +rxt(5)*y(7) +rxt(6)*y(8) &
                  +rxt(37)*y(20)*y(18) +rxt(41)*y(19)*y(19)
         loss(34) = (rxt(27)* y(1) + (rxt(25) +rxt(26))* y(4) + (rxt(30) +rxt(31) + &
                 rxt(32))* y(12) +rxt(28)* y(31) +rxt(29)* y(32) + rxt(22) + rxt(23) &
                  + rxt(24) + het_rates(3))* y(3)
         prod(34) =rxt(2)*y(1)
         loss(38) = ((rxt(79) +rxt(80))* y(19) + het_rates(17))* y(17)
         prod(38) = (rxt(11) +rxt(12) +rxt(81)*y(8) +rxt(82)*y(19) +rxt(83)*y(2)) &
                 *y(16) +.050_r8*rxt(91)*y(30)*y(1)
         loss(1) = ( + het_rates(5))* y(5)
         prod(1) = 0._r8
         loss(50) = (rxt(50)* y(1) +rxt(48)* y(2) +rxt(57)* y(8) +rxt(84)* y(13) &
                  +rxt(49)* y(20) +rxt(71)* y(24) + het_rates(6))* y(6)
         prod(50) = (rxt(5) +.500_r8*rxt(94) +rxt(51)*y(2))*y(7) &
                  +2.000_r8*rxt(25)*y(4)*y(3) +rxt(7)*y(8)
         loss(48) = (rxt(53)* y(1) + (rxt(51) +rxt(52))* y(2) +rxt(54)* y(8) +rxt(56) &
                 * y(19) +rxt(61)* y(20) +rxt(72)* y(24) + rxt(5) + rxt(94) &
                  + het_rates(7))* y(7)
         prod(48) = (rxt(48)*y(2) +rxt(49)*y(20) +rxt(50)*y(1) + &
                 2.000_r8*rxt(57)*y(8) +rxt(71)*y(24) +rxt(84)*y(13))*y(6) + (rxt(6) + &
                 rxt(58)*y(2) +rxt(59)*y(19) +rxt(60)*y(20))*y(8) + (rxt(9) +rxt(63) + &
                 rxt(62)*y(19))*y(10) +rxt(55)*y(11) +rxt(16)*y(29)
         loss(45) = (rxt(39)* y(1) +rxt(38)* y(2) +rxt(56)* y(7) +rxt(59)* y(8) &
                  +rxt(62)* y(10) +rxt(78)* y(12) +rxt(89)* y(14) +rxt(88)* y(15) &
                  +rxt(82)* y(16) + (rxt(79) +rxt(80))* y(17) + 2._r8*(rxt(41) + &
                 rxt(42))* y(19) +rxt(40)* y(20) +rxt(43)* y(21) + (rxt(68) +rxt(69)) &
                 * y(24) +rxt(90)* y(30) +rxt(95)* y(34) + (rxt(96) +rxt(97))* y(35) &
                  + het_rates(19))* y(19)
         prod(45) = (2.000_r8*rxt(35)*y(18) +rxt(44)*y(2) +rxt(45)*y(1) + &
                 rxt(49)*y(6) +rxt(60)*y(8) +rxt(66)*y(22))*y(20) + (rxt(47)*y(21) + &
                 rxt(83)*y(16))*y(2) + (2.000_r8*rxt(24) +rxt(30)*y(12))*y(3) &
                  + (rxt(10) +.300_r8*rxt(89)*y(19))*y(14) +rxt(34)*y(18)*y(1) &
                  +.500_r8*rxt(94)*y(7) +rxt(8)*y(10) +2.000_r8*rxt(13)*y(21) +rxt(14) &
                 *y(28)
         loss(43) = (rxt(58)* y(2) +rxt(57)* y(6) +rxt(54)* y(7) +rxt(81)* y(16) &
                  +rxt(59)* y(19) +rxt(60)* y(20) +rxt(98)* y(35) + rxt(6) + rxt(7) &
                  + rxt(93) + het_rates(8))* y(8)
         prod(43) = (rxt(52)*y(2) +rxt(53)*y(1))*y(7) +rxt(8)*y(10) +rxt(55)*y(11) &
                  +rxt(15)*y(29)
         loss(2) = ( + het_rates(9))* y(9)
         prod(2) = (rxt(93) +rxt(81)*y(16) +rxt(98)*y(35))*y(8) + (.500_r8*rxt(94) + &
                 rxt(56)*y(19))*y(7) +2.000_r8*rxt(92)*y(11)
         loss(35) = (rxt(62)* y(19) + rxt(8) + rxt(9) + rxt(63) + het_rates(10)) &
                 * y(10)
         prod(35) =rxt(61)*y(20)*y(7)
         loss(28) = ( + rxt(55) + rxt(92) + het_rates(11))* y(11)
         prod(28) =rxt(54)*y(8)*y(7)
         loss(39) = (rxt(84)* y(6) + 2._r8*(rxt(86) +rxt(87))* y(13) +rxt(85)* y(20) &
                  + het_rates(13))* y(13)
         prod(39) = (rxt(30)*y(3) +rxt(78)*y(19))*y(12) +1.860_r8*rxt(91)*y(30)*y(1) &
                  +.700_r8*rxt(89)*y(19)*y(14)
         loss(36) = (rxt(89)* y(19) + rxt(10) + het_rates(14))* y(14)
         prod(36) =rxt(85)*y(20)*y(13)
         loss(42) = (rxt(83)* y(2) +rxt(81)* y(8) +rxt(82)* y(19) + rxt(11) + rxt(12) &
                  + het_rates(16))* y(16)
         prod(42) = (rxt(84)*y(6) +2.000_r8*rxt(86)*y(13) +rxt(87)*y(13))*y(13) &
                  + (rxt(31)*y(12) +rxt(32)*y(12))*y(3) + (rxt(10) + &
                 .300_r8*rxt(89)*y(19))*y(14) +.870_r8*rxt(91)*y(30)*y(1) &
                  +rxt(88)*y(19)*y(15)
         loss(40) = (rxt(34)* y(1) + (rxt(35) +rxt(36) +rxt(37))* y(20) + rxt(33) &
                  + het_rates(18))* y(18)
         prod(40) = (rxt(38)*y(2) +rxt(79)*y(17) +rxt(82)*y(16))*y(19) +rxt(31)*y(12) &
                 *y(3) +rxt(10)*y(14) +2.000_r8*rxt(11)*y(16)
         loss(44) = (rxt(45)* y(1) +rxt(44)* y(2) +rxt(49)* y(6) +rxt(61)* y(7) &
                  +rxt(60)* y(8) +rxt(85)* y(13) + (rxt(35) +rxt(36) +rxt(37))* y(18) &
                  +rxt(40)* y(19) + 2._r8*rxt(46)* y(20) + (rxt(65) +rxt(66))* y(22) &
                  +rxt(70)* y(24) + rxt(99) + het_rates(20))* y(20)
         prod(44) = (rxt(39)*y(1) +rxt(43)*y(21) +rxt(59)*y(8) +rxt(68)*y(24) + &
                 rxt(80)*y(17) +rxt(88)*y(15) +.500_r8*rxt(97)*y(35))*y(19) &
                  + (rxt(47)*y(21) +rxt(83)*y(16))*y(2) + (rxt(9) +rxt(63))*y(10) &
                  + (rxt(84)*y(6) +2.000_r8*rxt(86)*y(13))*y(13) &
                  +.060_r8*rxt(91)*y(30)*y(1) +rxt(31)*y(12)*y(3) +rxt(81)*y(16)*y(8) &
                  +rxt(33)*y(18)
         loss(33) = (rxt(47)* y(2) +rxt(43)* y(19) + rxt(13) + het_rates(21))* y(21)
         prod(33) = (.500_r8*rxt(99) +rxt(46)*y(20))*y(20) +rxt(42)*y(19)*y(19)
         loss(41) = (rxt(64)* y(1) + (rxt(65) +rxt(66))* y(20) + het_rates(22))* y(22)
         prod(41) = (rxt(67)*y(2) +rxt(68)*y(19) +rxt(71)*y(6) + &
                 2.000_r8*rxt(73)*y(24) +rxt(75)*y(24))*y(24) &
                  + (3.000_r8*rxt(28)*y(31) +2.000_r8*rxt(29)*y(32))*y(3) +rxt(14) &
                 *y(28) +rxt(15)*y(29)
         loss(3) = ( + het_rates(23))* y(23)
         prod(3) =rxt(74)*y(24)*y(24)
         loss(46) = (rxt(67)* y(2) +rxt(71)* y(6) +rxt(72)* y(7) + (rxt(68) +rxt(69)) &
                 * y(19) +rxt(70)* y(20) + 2._r8*(rxt(73) +rxt(74) +rxt(75) +rxt(76)) &
                 * y(24) + het_rates(24))* y(24)
         prod(46) = (rxt(64)*y(1) +rxt(66)*y(20))*y(22) +2.000_r8*rxt(77)*y(26) &
                  +rxt(16)*y(29)
         loss(4) = ( + het_rates(25))* y(25)
         prod(4) =rxt(75)*y(24)*y(24)
         loss(26) = ( + rxt(77) + het_rates(26))* y(26)
         prod(26) =rxt(76)*y(24)*y(24)
         loss(5) = ( + het_rates(27))* y(27)
         prod(5) =rxt(69)*y(24)*y(19) +rxt(65)*y(22)*y(20)
         loss(29) = ( + rxt(14) + het_rates(28))* y(28)
         prod(29) =rxt(70)*y(24)*y(20)
         loss(32) = ( + rxt(15) + rxt(16) + het_rates(29))* y(29)
         prod(32) =rxt(72)*y(24)*y(7)
         loss(37) = (rxt(91)* y(1) +rxt(90)* y(19) + het_rates(30))* y(30)
         prod(37) = 0._r8
         loss(31) = (rxt(88)* y(19) + het_rates(15))* y(15)
         prod(31) =rxt(87)*y(13)*y(13)
         loss(27) = (rxt(95)* y(19) + het_rates(34))* y(34)
         prod(27) = (rxt(96)*y(19) +.500_r8*rxt(97)*y(19) +rxt(98)*y(8))*y(35)
         loss(30) = (rxt(98)* y(8) + (rxt(96) +rxt(97))* y(19) + het_rates(35))* y(35)
         prod(30) = 0._r8
         loss(6) = ( + het_rates(38))* y(38)
         prod(6) = 0._r8
         loss(7) = ( + het_rates(36))* y(36)
         prod(7) =rxt(95)*y(34)*y(19)
         loss(8) = ( + het_rates(37))* y(37)
         prod(8) = 0._r8
         loss(9) = ( + het_rates(39))* y(39)
         prod(9) = 0._r8
         loss(10) = ( + het_rates(40))* y(40)
         prod(10) = 0._r8
         loss(11) = ( + het_rates(41))* y(41)
         prod(11) = 0._r8
         loss(12) = ( + het_rates(42))* y(42)
         prod(12) = 0._r8
         loss(13) = ( + het_rates(43))* y(43)
         prod(13) = 0._r8
         loss(14) = ( + het_rates(44))* y(44)
         prod(14) = 0._r8
         loss(15) = ( + het_rates(45))* y(45)
         prod(15) = 0._r8
         loss(16) = ( + het_rates(46))* y(46)
         prod(16) = 0._r8
         loss(17) = ( + het_rates(47))* y(47)
         prod(17) = 0._r8
         loss(18) = ( + het_rates(48))* y(48)
         prod(18) = 0._r8
         loss(19) = ( + het_rates(49))* y(49)
         prod(19) = 0._r8
         loss(20) = ( + het_rates(50))* y(50)
         prod(20) = 0._r8
         loss(21) = ( + het_rates(51))* y(51)
         prod(21) = 0._r8
         loss(22) = ( + het_rates(52))* y(52)
         prod(22) = 0._r8
         loss(23) = ( + het_rates(53))* y(53)
         prod(23) = 0._r8
         loss(24) = ( + het_rates(54))* y(54)
         prod(24) = 0._r8
         loss(25) = ( + het_rates(55))* y(55)
         prod(25) = 0._r8
      end subroutine imp_prod_loss
      end module mo_prod_loss
