module gravmag_sphere_harmonics
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none
  private
  public :: schmidt_legendre, schmidt_basis

contains

  ! real schmidt semi-normalization without the condon-shortley phase
  ! differentiate the recurrence in colatitude and recur p/sin(theta) for m > 0
  ! the latter remains finite at the poles without a clamped denominator
  pure subroutine schmidt_legendre(lm, theta, p, dp, p_over_sin)
    integer(int32), intent(in) :: lm
    real(real64), intent(in) :: theta
    real(real64), intent(out) :: p(0:lm,0:lm), dp(0:lm,0:lm), p_over_sin(0:lm,0:lm)
    integer(int32) :: l, m
    real(real64) :: st, ct, factor, alpha, beta, denominator

    st = sin(theta)
    ct = cos(theta)
    p = 0.0_real64
    dp = 0.0_real64
    p_over_sin = 0.0_real64
    p(0,0) = 1.0_real64
    if (lm >= 1) then
      p(1,1) = st
      dp(1,1) = ct
      p_over_sin(1,1) = 1.0_real64
    end if
    do m=2, lm
      factor = sqrt(real(2*m-1,real64)/real(2*m,real64))
      p(m,m) = factor*st*p(m-1,m-1)
      dp(m,m) = factor*(ct*p(m-1,m-1) + st*dp(m-1,m-1))
      p_over_sin(m,m) = factor*p(m-1,m-1)
    end do
    do m=0, lm-1
      factor = sqrt(real(2*m+1,real64))
      p(m+1,m) = factor*ct*p(m,m)
      dp(m+1,m) = factor*(ct*dp(m,m) - st*p(m,m))
      if (m > 0) p_over_sin(m+1,m) = factor*ct*p_over_sin(m,m)
    end do
    do m=0, lm
      do l=m+2, lm
        denominator = real(l-m,real64)*real(l+m,real64)
        alpha = real(2*l-1,real64)/sqrt(denominator)
        beta = sqrt(real(l-m-1,real64)*real(l+m-1,real64)/denominator)
        p(l,m) = alpha*ct*p(l-1,m) - beta*p(l-2,m)
        dp(l,m) = alpha*(ct*dp(l-1,m) - st*p(l-1,m)) - beta*dp(l-2,m)
        if (m > 0) p_over_sin(l,m) = alpha*ct*p_over_sin(l-1,m) - beta*p_over_sin(l-2,m)
      end do
    end do
  end subroutine schmidt_legendre

  ! shared vector basis at the coefficient reference radius
  pure subroutine schmidt_basis(theta, phi, lm, ig, ih, basis_br, basis_bt, basis_bp)
    real(real64), intent(in) :: theta, phi
    integer(int32), intent(in) :: lm, ig(0:lm,0:lm), ih(0:lm,0:lm)
    real(real64), intent(out) :: basis_br(:), basis_bt(:), basis_bp(:)
    real(real64) :: p(0:lm,0:lm), dp(0:lm,0:lm), p_over_sin(0:lm,0:lm), cm, sm
    integer(int32) :: l, m, idx

    call schmidt_legendre(lm, theta, p, dp, p_over_sin)
    basis_br = 0.0_real64
    basis_bt = 0.0_real64
    basis_bp = 0.0_real64
    do l=1, lm
      do m=0, l
        cm = cos(real(m,real64)*phi)
        sm = sin(real(m,real64)*phi)
        idx = ig(l,m)
        basis_br(idx) = real(l+1,real64)*p(l,m)*cm
        basis_bt(idx) = -dp(l,m)*cm
        if (m > 0) basis_bp(idx) = real(m,real64)*p_over_sin(l,m)*sm
        if (m > 0) then
          idx = ih(l,m)
          basis_br(idx) = real(l+1,real64)*p(l,m)*sm
          basis_bt(idx) = -dp(l,m)*sm
          basis_bp(idx) = -real(m,real64)*p_over_sin(l,m)*cm
        end if
      end do
    end do
  end subroutine schmidt_basis

end module gravmag_sphere_harmonics
