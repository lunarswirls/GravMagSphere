module gravmag_sphere_harmonics
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none
  private
  public :: schmidt_legendre, schmidt_basis, schmidt_seeds, schmidt_start, schmidt_advance
  public :: schmidt_values, schmidt_point, coefficient_index
  integer, parameter :: scale_step = 512
  real(real64), parameter :: scale_low = 2.0_real64**(-256), scale_high = 2.0_real64**256

contains

  ! fukushima extended-range arithmetic with radix 2**512 and double mantissas
  ! sectoral quotients are retained before multiplication by the small sine
  pure subroutine schmidt_seeds(lm, theta, seeds, powers)
    integer(int32), intent(in) :: lm
    real(real64), intent(in) :: theta(:)
    real(real64), allocatable, intent(out) :: seeds(:,:)
    integer(int32), allocatable, intent(out) :: powers(:,:)
    real(real64) :: st(size(theta)), factor
    integer(int32) :: m, i, st_power(size(theta))
    allocate(seeds(size(theta),0:lm), powers(size(theta),0:lm))
    st = sin(theta)
    st_power = 0
    do i=1, size(theta)
      do while (abs(st(i)) < scale_low .and. st(i) /= 0.0_real64)
        st(i) = scale(st(i), scale_step)
        st_power(i) = st_power(i)-1
      end do
    end do
    seeds(:,0) = 1.0_real64
    powers = 0
    if (lm >= 1) seeds(:,1) = 1.0_real64
    do m=2, lm
      factor = sqrt(real(2*m-1,real64)/real(2*m,real64))
      seeds(:,m) = factor*st*seeds(:,m-1)
      powers(:,m) = powers(:,m-1)+st_power
      do i=1, size(theta)
        do while (abs(seeds(i,m)) < scale_low .and. seeds(i,m) /= 0.0_real64)
          seeds(i,m) = scale(seeds(i,m), scale_step)
          powers(i,m) = powers(i,m)-1
        end do
        ! the scaled sine can exceed one at extremely small polar angles
        if (abs(seeds(i,m)) > scale_high) then
          seeds(i,m) = scale(seeds(i,m), -scale_step)
          powers(i,m) = powers(i,m)+1
        end if
      end do
    end do
  end subroutine schmidt_seeds

  pure subroutine schmidt_start(m, st, ct, seed, seed_power, previous, current, powers)
    integer(int32), intent(in) :: m, seed_power(:)
    real(real64), intent(in) :: st(:), ct(:), seed(:)
    real(real64), intent(out) :: previous(:,:), current(:,:)
    integer(int32), intent(out) :: powers(:)
    previous = 0.0_real64
    powers = seed_power
    if (m == 0) then
      current(:,1) = 1.0_real64
      current(:,2:3) = 0.0_real64
    else
      current(:,1) = st*seed
      current(:,2) = real(m,real64)*ct*seed
      current(:,3) = seed
    end if
  end subroutine schmidt_start

  ! both recurrence states and all three quantities share a per-latitude exponent
  pure subroutine schmidt_advance(l, m, st, ct, previous, current, powers)
    integer(int32), intent(in) :: l, m
    real(real64), intent(in) :: st(:), ct(:)
    real(real64), intent(inout) :: previous(:,:), current(:,:)
    integer(int32), intent(inout) :: powers(:)
    real(real64) :: next(size(st),3), denominator, alpha, beta, largest
    integer(int32) :: i
    denominator = real(l-m,real64)*real(l+m,real64)
    alpha = real(2*l-1,real64)/sqrt(denominator)
    beta = sqrt(real(l-m-1,real64)*real(l+m-1,real64)/denominator)
    next(:,1) = alpha*ct*current(:,1)-beta*previous(:,1)
    next(:,2) = alpha*(ct*current(:,2)-st*current(:,1))-beta*previous(:,2)
    next(:,3) = alpha*ct*current(:,3)-beta*previous(:,3)
    previous = current
    current = next
    do i=1, size(st)
      largest = max(maxval(abs(previous(i,:))), maxval(abs(current(i,:))))
      if (largest > scale_high) then
        previous(i,:) = scale(previous(i,:), -scale_step)
        current(i,:) = scale(current(i,:), -scale_step)
        powers(i) = powers(i)+1
      end if
    end do
  end subroutine schmidt_advance

  pure subroutine schmidt_values(current, powers, values)
    real(real64), intent(in) :: current(:,:)
    integer(int32), intent(in) :: powers(:)
    real(real64), intent(out) :: values(:,:)
    integer(int32) :: i
    values = current
    do i=1, size(powers)
      if (powers(i) < -2) then
        values(i,:) = 0.0_real64
      else if (powers(i) /= 0) then
        values(i,:) = scale(current(i,:), scale_step*powers(i))
      end if
    end do
  end subroutine schmidt_values

  pure integer(int32) function coefficient_index(l, m) result(index)
    integer(int32), intent(in) :: l, m
    index = l*l + max(0,2*m-1)
  end function coefficient_index

  ! schmidt semi-normalization without the condon-shortley phase
  pure subroutine schmidt_legendre(lm, theta, p, dp, p_over_sin)
    integer(int32), intent(in) :: lm
    real(real64), intent(in) :: theta
    real(real64), intent(out) :: p(0:lm,0:lm), dp(0:lm,0:lm), p_over_sin(0:lm,0:lm)
    integer(int32) :: l, m, powers(1)
    integer(int32), allocatable :: seed_powers(:,:)
    real(real64), allocatable :: seeds(:,:)
    real(real64) :: st(1), ct(1), previous(1,3), current(1,3), values(1,3)
    st = sin(theta)
    ct = cos(theta)
    call schmidt_seeds(lm, [theta], seeds, seed_powers)
    p = 0.0_real64
    dp = 0.0_real64
    p_over_sin = 0.0_real64
    do m=0, lm
      call schmidt_start(m, st, ct, seeds(:,m), seed_powers(:,m), previous, current, powers)
      do l=m, lm
        if (l > m) call schmidt_advance(l, m, st, ct, previous, current, powers)
        call schmidt_values(current, powers, values)
        p(l,m) = values(1,1)
        dp(l,m) = values(1,2)
        p_over_sin(l,m) = values(1,3)
      end do
    end do
  end subroutine schmidt_legendre

  ! point synthesis streams orders without a square basis workspace
  pure subroutine schmidt_point(lm, theta, phi, coefficients, field)
    integer(int32), intent(in) :: lm
    real(real64), intent(in) :: theta, phi, coefficients(:)
    real(real64), intent(out) :: field(3)
    real(real64), allocatable :: seeds(:,:)
    integer(int32), allocatable :: seed_powers(:,:)
    real(real64) :: st(1), ct(1), previous(1,3), current(1,3), values(1,3), cm, sm, g, h, amplitude
    integer(int32) :: l, m, index, powers(1)
    st = sin(theta)
    ct = cos(theta)
    call schmidt_seeds(lm, [theta], seeds, seed_powers)
    field = 0.0_real64
    do m=0, lm
      call schmidt_start(m, st, ct, seeds(:,m), seed_powers(:,m), previous, current, powers)
      cm = cos(real(m,real64)*phi)
      sm = sin(real(m,real64)*phi)
      do l=m, lm
        if (l > m) call schmidt_advance(l, m, st, ct, previous, current, powers)
        if (l == 0) cycle
        call schmidt_values(current, powers, values)
        index = coefficient_index(l,m)
        g = coefficients(index)
        h = 0.0_real64
        if (m > 0) h = coefficients(index+1)
        amplitude = g*cm+h*sm
        field(1) = field(1)+real(l+1,real64)*values(1,1)*amplitude
        field(2) = field(2)-values(1,2)*amplitude
        field(3) = field(3)+real(m,real64)*values(1,3)*(g*sm-h*cm)
      end do
    end do
  end subroutine schmidt_point

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
