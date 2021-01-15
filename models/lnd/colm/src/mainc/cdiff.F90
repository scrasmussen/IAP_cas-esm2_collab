
subroutine cdiff(dtime, z, dz, tss, cpool_fast, cpool_slow)

  use precision
  use paramodel, only : nl_soil, nl_csoil
  implicit none

  real(r8), parameter     :: D = 0.001/(86400*365)  !diffusion constant 

  real(r8), intent(in)    :: dtime                  !model time step [second]
  real(r8), intent(in)    ::  z(nl_soil)            !node depth [m]
  real(r8), intent(in)    :: dz(nl_soil)            !layer thickiness [m]
  real(r8), intent(in)    :: tss(nl_soil)           !soil temperature [K]
  real(r8), intent(inout) :: cpool_fast(nl_soil)    !fast soil carbon pool (kg m-2)
  real(r8), intent(inout) :: cpool_slow(nl_soil)    !slow soil carbon pool (kg m-2)

 !local variables

  real(r8) cf(1:nl_csoil)   !density of fast carbon pools
  real(r8) cs(1:nl_csoil)   !density of slow carbon pools

  real(r8) at(1:nl_csoil)   !"a" vector for tridiagonal matrix
  real(r8) bt(1:nl_csoil)   !"b" vector for tridiagonal matrix
  real(r8) ct(1:nl_csoil)   !"c" vector for tridiagonal matrix
  real(r8) rt(1:nl_csoil)   !"r" vector for tridiagonal solution

  real(r8) dzm              !used in computing tridiagonal matrix
  real(r8) dzp              !used in computing tridiagonal matrix
  real(r8) dzx              !used in computing tridiagonal matrix
  real(r8) cnfac            !Crank Nicholson factor between 0 and 1

  real(r8) cfmass0, csmass0, cfmass, csmass

  integer j

!=======================================================================

  cfmass0 = sum(cpool_fast(1:nl_csoil))
  csmass0 = sum(cpool_slow(1:nl_csoil))

  cf = cpool_fast(1:nl_csoil)/dz(1:nl_csoil)
  cs = cpool_slow(1:nl_csoil)/dz(1:nl_csoil)

! set up vector r and vectors a, b, c that define tridiagonal matrix
! a: subdiagonal elements
! b: diagonal elements
! c: superdiagonal elements

  j     = 1
  dzp   = (z(j+1)-z(j))*(z(j+1)+z(j))
  at(j) = 0.
  bt(j) = 1.+ 2*D*dtime/dzp
  ct(j) =   - 2*D*dtime/dzp
  rt(j) = cf(j)

  do j = 2, nl_csoil - 1
     dzm   = (z(j)-z(j-1))
     dzp   = (z(j+1)-z(j))
     dzx   = (z(j+1)-z(j-1))
     at(j) =   - 2*D*dtime/(dzm*dzx)
     bt(j) = 1.+ 2*D*dtime/(dzm*dzp)
     ct(j) =   - 2*D*dtime/(dzp*dzx)
     rt(j) = cf(j)
  end do

!***If nl_csoil==nl_soil***!
! j     =  nl_csoil
! dzm   = (z(j)-z(j-1))*(z(j)-z(j-1))
! at(j) =   - D*dtime/dzm
! bt(j) = 1.+ D*dtime/dzm
! ct(j) = 0.
! rt(j) = cf(j)

  j     = nl_csoil
  dzm   = (z(j)-z(j-1))*(z(j+1)-z(j-1))
  at(j) =   - 2*D*dtime/dzm
  bt(j) = 1.+ 2*D*dtime/dzm
  ct(j) = 0.
  rt(j) = cf(j)

! solve for soil fast carbon pool diffusion
  rt(1:nl_csoil) = cf(1:nl_csoil)
  call tridia (nl_csoil ,at ,bt ,ct ,rt ,cf(1:nl_csoil))

! solve for soil slow carbon pool diffusion
  rt(1:nl_csoil) = cs(1:nl_csoil)
  call tridia (nl_csoil ,at ,bt ,ct ,rt ,cs(1:nl_csoil))

  cpool_fast(1:nl_csoil) = cf(1:nl_csoil)*dz(1:nl_csoil)
  cpool_slow(1:nl_csoil) = cs(1:nl_csoil)*dz(1:nl_csoil)

  cfmass = sum(cpool_fast(1:nl_csoil))
  csmass = sum(cpool_slow(1:nl_csoil))

  if(abs(cfmass-cfmass0)>1.e-6 .or. abs(csmass-csmass0)>1.e-6) then
     write(6,*) 'Mass conservation violation:', cfmass0, cfmass, csmass0, csmass
     call abort
  end if

end subroutine cdiff
