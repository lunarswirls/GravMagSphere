module gravmag_sphere_transform
  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  use gravmag_sphere_harmonics, only: schmidt_seeds, schmidt_start, schmidt_advance, schmidt_values, coefficient_index
  implicit none
  private
  public :: transform_grid, gauss_nodes, fft_roots, native_fft, sample_fourier, analyze_fourier
  public :: synthesize_fourier, inverse_ring, filter_coefficients, spectral_norm, transform_memory
  real(real64), parameter :: pi = acos(-1.0_real64)
  complex(real64), parameter :: imaginary = (0.0_real64,1.0_real64)

  abstract interface
    subroutine ring_sampler(theta, field)
      import :: real64
      real(real64), intent(in) :: theta
      real(real64), intent(out) :: field(:,:)
    end subroutine ring_sampler
  end interface
contains

  subroutine transform_grid(lm, ntheta, nphi)
    integer(int32), intent(in) :: lm
    integer(int32), intent(inout) :: ntheta, nphi
    integer(int32) :: requested
    if (lm < 1 .or. lm > 46339) error stop 'ERROR: degree exceeds the supported integer indexing range'
    ntheta = max(ntheta,lm+1)
    requested = max(nphi,2*lm+1)
    nphi = 1
    do while (nphi < requested)
      if (nphi > huge(nphi)/2) error stop 'ERROR: longitude count exceeds the integer range'
      nphi = 2*nphi
    end do
  end subroutine transform_grid

  ! symmetric newton iteration for legendre roots and exact quadrature weights
  subroutine gauss_nodes(n, theta, weights)
    integer(int32), intent(in) :: n
    real(real64), intent(out) :: theta(n), weights(n)
    integer(int32) :: i, k, iteration
    real(real64) :: z, correction, p0, p1, p2, derivative
    !$omp parallel do default(shared) private(i,k,iteration,z,correction,p0,p1,p2,derivative)
    do i=1, (n+1)/2
      z = cos(pi*(real(i,real64)-0.25_real64)/(real(n,real64)+0.5_real64))
      do iteration=1, 30
        p0 = 1.0_real64
        p1 = z
        do k=2, n
          p2 = (real(2*k-1,real64)*z*p1-real(k-1,real64)*p0)/real(k,real64)
          p0 = p1
          p1 = p2
        end do
        derivative = real(n,real64)*(p0-z*p1)/(1.0_real64-z*z)
        correction = p1/derivative
        z = z-correction
        if (abs(correction) <= 2.0_real64*epsilon(z)*max(abs(z),0.1_real64)) exit
      end do
      if (iteration > 30) error stop 'ERROR: gauss-legendre root iteration did not converge'
      ! reevaluate the derivative at the final root
      p0 = 1.0_real64
      p1 = z
      do k=2, n
        p2 = (real(2*k-1,real64)*z*p1-real(k-1,real64)*p0)/real(k,real64)
        p0 = p1
        p1 = p2
      end do
      derivative = real(n,real64)*(p0-z*p1)/(1.0_real64-z*z)
      theta(i) = acos(z)
      theta(n+1-i) = pi-theta(i)
      weights(i) = 2.0_real64/((1.0_real64-z*z)*derivative*derivative)
      weights(n+1-i) = weights(i)
    end do
    !$omp end parallel do
  end subroutine gauss_nodes

  subroutine fft_roots(n, roots)
    integer(int32), intent(in) :: n
    complex(real64), intent(out) :: roots(0:n/2-1)
    integer(int32) :: i
    do i=0, n/2-1
      roots(i) = cmplx(cos(2*pi*real(i,real64)/real(n,real64)), &
                      -sin(2*pi*real(i,real64)/real(n,real64)),kind=real64)
    end do
  end subroutine fft_roots

  ! native radix-two fft with normalized forward and unnormalized inverse transforms
  subroutine native_fft(values, roots, inverse)
    complex(real64), intent(inout) :: values(0:)
    complex(real64), intent(in) :: roots(0:)
    logical, intent(in) :: inverse
    integer(int32) :: n, i, j, bit, width, half, first, k, stride
    complex(real64) :: temporary, twiddle
    n = size(values)
    if (n < 2 .or. iand(n,n-1) /= 0 .or. size(roots) /= n/2) error stop 'ERROR: invalid native fft dimensions'
    j = 0
    do i=0, n-1
      if (i < j) then
        temporary = values(i)
        values(i) = values(j)
        values(j) = temporary
      end if
      bit = n/2
      do while (bit > 0)
        if (iand(j,bit) == 0) exit
        j = ieor(j,bit)
        bit = bit/2
      end do
      j = ieor(j,bit)
    end do
    width = 2
    do
      half = width/2
      stride = n/width
      do first=0, n-1, width
        do k=0, half-1
          twiddle = roots(k*stride)
          if (inverse) twiddle = conjg(twiddle)
          temporary = twiddle*values(first+k+half)
          values(first+k+half) = values(first+k)-temporary
          values(first+k) = values(first+k)+temporary
        end do
      end do
      if (width == n) exit
      width = 2*width
    end do
    if (.not. inverse) values = values/real(n,real64)
  end subroutine native_fft

  subroutine sample_fourier(lm, nphi, theta, weights, joint, sampler, fourier, energy)
    integer(int32), intent(in) :: lm, nphi
    real(real64), intent(in) :: theta(:), weights(:), joint
    procedure(ring_sampler) :: sampler
    complex(real64), intent(out) :: fourier(:,0:,:)
    real(real64), intent(out) :: energy
    complex(real64), allocatable :: roots(:), work(:)
    real(real64), allocatable :: field(:,:), energies(:)
    integer(int32) :: i, component
    allocate(roots(0:nphi/2-1),energies(size(theta)))
    call fft_roots(nphi,roots)
    !$omp parallel default(shared) private(i,component,work,field)
    allocate(work(0:nphi-1),field(nphi,3))
    !$omp do schedule(dynamic)
    do i=1, size(theta)
      call sampler(theta(i),field)
      energies(i) = 2*pi*weights(i)*sum(field(:,1)**2+joint*(field(:,2)**2+field(:,3)**2))/real(nphi,real64)
      do component=1, 3
        work = cmplx(field(:,component),0.0_real64,kind=real64)
        call native_fft(work,roots,.false.)
        fourier(i,:,component) = work(0:lm)
      end do
    end do
    !$omp end do
    deallocate(work,field)
    !$omp end parallel
    energy = sum(energies)
  end subroutine sample_fourier

  pure real(real64) function spectral_norm(l,joint) result(norm)
    integer(int32), intent(in) :: l
    real(real64), intent(in) :: joint
    norm = 4*pi*(real(l+1,real64)**2+joint*real(l,real64)*real(l+1,real64))/real(2*l+1,real64)
  end function spectral_norm

  ! orthogonal vector projection with one shared exterior-potential coefficient set
  subroutine analyze_fourier(lm, theta, weights, fourier, joint, coefficients, degree_energy)
    integer(int32), intent(in) :: lm
    real(real64), intent(in) :: theta(:), weights(:), joint
    complex(real64), intent(in) :: fourier(:,0:,:)
    real(real64), intent(out) :: coefficients(lm*(lm+2)), degree_energy(lm)
    real(real64), allocatable :: seeds(:,:), st(:), ct(:), previous(:,:), current(:,:), values(:,:)
    integer(int32), allocatable :: seed_powers(:,:), powers(:)
    complex(real64) :: projection
    integer(int32) :: l, m, index, ntheta
    ntheta = size(theta)
    allocate(st(ntheta),ct(ntheta))
    st = sin(theta)
    ct = cos(theta)
    call schmidt_seeds(lm,theta,seeds,seed_powers)
    coefficients = 0.0_real64
    !$omp parallel default(shared) private(l,m,index,projection,previous,current,values,powers)
    allocate(previous(ntheta,3),current(ntheta,3),values(ntheta,3),powers(ntheta))
    !$omp do schedule(dynamic)
    do m=0, lm
      call schmidt_start(m,st,ct,seeds(:,m),seed_powers(:,m),previous,current,powers)
      do l=m, lm
        if (l > m) call schmidt_advance(l,m,st,ct,previous,current,powers)
        if (l == 0) cycle
        call schmidt_values(current,powers,values)
        projection = 2*pi*sum(weights*(real(l+1,real64)*values(:,1)*fourier(:,m,1) &
                     -joint*values(:,2)*fourier(:,m,2)+imaginary*joint*real(m,real64)*values(:,3)*fourier(:,m,3))) &
                     /spectral_norm(l,joint)
        index = coefficient_index(l,m)
        coefficients(index) = real(projection,real64)
        if (m > 0) coefficients(index+1) = -aimag(projection)
      end do
    end do
    !$omp end do
    deallocate(previous,current,values,powers)
    !$omp end parallel
    do l=1, lm
      degree_energy(l) = spectral_norm(l,joint)*sum(coefficients(l*l:(l+1)*(l+1)-1)**2)
    end do
  end subroutine analyze_fourier

  ! synthesize longitude fourier coefficients without allocating a spatial grid
  subroutine synthesize_fourier(lm, theta, coefficients, fourier)
    integer(int32), intent(in) :: lm
    real(real64), intent(in) :: theta(:), coefficients(:)
    complex(real64), intent(out) :: fourier(:,0:,:)
    real(real64), allocatable :: seeds(:,:), st(:), ct(:), previous(:,:), current(:,:), values(:,:)
    integer(int32), allocatable :: seed_powers(:,:), powers(:)
    integer(int32) :: l, m, index, ntheta
    complex(real64) :: amplitude
    ntheta = size(theta)
    allocate(st(ntheta),ct(ntheta))
    st = sin(theta)
    ct = cos(theta)
    call schmidt_seeds(lm,theta,seeds,seed_powers)
    fourier = (0.0_real64,0.0_real64)
    !$omp parallel default(shared) private(l,m,index,amplitude,previous,current,values,powers)
    allocate(previous(ntheta,3),current(ntheta,3),values(ntheta,3),powers(ntheta))
    !$omp do schedule(dynamic)
    do m=0, lm
      call schmidt_start(m,st,ct,seeds(:,m),seed_powers(:,m),previous,current,powers)
      do l=m, lm
        if (l > m) call schmidt_advance(l,m,st,ct,previous,current,powers)
        if (l == 0) cycle
        call schmidt_values(current,powers,values)
        index = coefficient_index(l,m)
        amplitude = cmplx(coefficients(index),0.0_real64,kind=real64)
        if (m > 0) amplitude = 0.5_real64*cmplx(coefficients(index),-coefficients(index+1),kind=real64)
        fourier(:,m,1) = fourier(:,m,1)+real(l+1,real64)*values(:,1)*amplitude
        fourier(:,m,2) = fourier(:,m,2)-values(:,2)*amplitude
        fourier(:,m,3) = fourier(:,m,3)-imaginary*real(m,real64)*values(:,3)*amplitude
      end do
    end do
    !$omp end do
    deallocate(previous,current,values,powers)
    !$omp end parallel
  end subroutine synthesize_fourier

  subroutine inverse_ring(fourier, roots, field)
    complex(real64), intent(in) :: fourier(0:,:), roots(0:)
    real(real64), intent(out) :: field(:,:)
    complex(real64) :: work(0:size(field,1)-1)
    integer(int32) :: component, m, lm, nphi
    nphi = size(field,1)
    lm = size(fourier,1)-1
    do component=1, 3
      work = (0.0_real64,0.0_real64)
      work(0:lm) = fourier(:,component)
      do m=1, lm
        work(nphi-m) = conjg(work(m))
      end do
      call native_fft(work,roots,.true.)
      field(:,component) = real(work,real64)
    end do
  end subroutine inverse_ring

  subroutine filter_coefficients(lm, raw, degree_energy, data_energy, joint, damping, power, &
                                 coefficients, residual, penalty, degrees_of_freedom)
    integer(int32), intent(in) :: lm
    real(real64), intent(in) :: raw(:), degree_energy(:), data_energy, joint, damping, power
    real(real64), intent(out) :: coefficients(lm*(lm+2)), residual, penalty, degrees_of_freedom
    real(real64) :: floor, norm, multiplier, weight, removed
    integer(int32) :: l
    floor = 0.0_real64
    do l=1, lm
      floor = floor+real(2*l+1,real64)*spectral_norm(l,joint)
    end do
    floor = epsilon(1.0_real64)*floor/real(lm*(lm+2),real64)
    penalty = 0.0_real64
    degrees_of_freedom = 0.0_real64
    ! separate the unresolved residual from damping to avoid subtracting nearly equal totals repeatedly
    residual = max(0.0_real64,data_energy-sum(degree_energy(1:lm)))
    do l=1, lm
      norm = spectral_norm(l,joint)
      weight = (real(l,real64)/real(lm,real64))**power
      multiplier = 1.0_real64/(1.0_real64+damping*weight+floor/norm)
      coefficients(l*l:(l+1)*(l+1)-1) = multiplier*raw(l*l:(l+1)*(l+1)-1)
      removed = 1.0_real64-multiplier
      residual = residual+degree_energy(l)*removed*removed
      penalty = penalty+degree_energy(l)*weight*multiplier*multiplier
      degrees_of_freedom = degrees_of_freedom+real(2*l+1,real64)*multiplier
    end do
  end subroutine filter_coefficients

  pure integer(int64) function transform_memory(lm,ntheta,nphi,correction) result(bytes)
    integer(int32), intent(in) :: lm, ntheta, nphi
    logical, intent(in) :: correction
    integer(int64) :: degree, latitudes, longitudes
    degree = int(lm,int64)
    latitudes = int(ntheta,int64)
    longitudes = int(nphi,int64)
    ! peak global workspaces plus conservative ring and recurrence storage for 32 threads
    bytes = 60_int64*latitudes*(degree+1)+16_int64*degree*(degree+2) &
          +32_int64*(80_int64*latitudes+80_int64*longitudes)
    if (correction) bytes = bytes+40_int64*latitudes*longitudes+48_int64*latitudes*(degree+1)
  end function transform_memory
end module gravmag_sphere_transform
