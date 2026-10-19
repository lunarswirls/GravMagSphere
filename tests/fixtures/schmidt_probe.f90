subroutine schmidt_probe(lm, theta, phi, p, dp, p_over_sin, basis) bind(c)
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use gravmag_sphere_harmonics, only: schmidt_legendre, schmidt_basis
  implicit none
  integer(c_int), value :: lm
  real(c_double), value :: theta, phi
  real(c_double), intent(out) :: p(0:lm,0:lm), dp(0:lm,0:lm), p_over_sin(0:lm,0:lm)
  real(c_double), intent(out) :: basis(3,lm*(lm+2))
  integer(c_int) :: ig(0:lm,0:lm), ih(0:lm,0:lm), l, m, index

  ig = 0
  ih = 0
  index = 0
  do l=1, lm
    do m=0, l
      index = index + 1
      ig(l,m) = index
      if (m > 0) then
        index = index + 1
        ih(l,m) = index
      end if
    end do
  end do
  call schmidt_legendre(lm, theta, p, dp, p_over_sin)
  call schmidt_basis(theta, phi, lm, ig, ih, basis(1,:), basis(2,:), basis(3,:))
end subroutine schmidt_probe
