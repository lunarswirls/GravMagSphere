subroutine native_roundtrip_probe(lm, ntheta, nphi, joint, coefficients, recovered, statistics) bind(c)
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use gravmag_sphere_transform
  implicit none
  integer(c_int), value :: lm,ntheta,nphi
  real(c_double), value :: joint
  real(c_double), intent(in) :: coefficients(lm*(lm+2))
  real(c_double), intent(out) :: recovered(lm*(lm+2)),statistics(3)
  real(c_double), allocatable :: theta(:),weights(:),degree_energy(:),field(:,:)
  complex(c_double), allocatable :: fourier(:,:,:),roots(:),work(:)
  integer(c_int) :: i,component
  allocate(theta(ntheta),weights(ntheta),degree_energy(lm),fourier(ntheta,0:lm,3),roots(0:nphi/2-1))
  call gauss_nodes(ntheta,theta,weights)
  call fft_roots(nphi,roots)
  call synthesize_fourier(lm,theta,coefficients,fourier)
  !$omp parallel default(shared) private(i,component,field,work)
  allocate(field(nphi,3),work(0:nphi-1))
  !$omp do schedule(dynamic)
  do i=1, ntheta
    call inverse_ring(fourier(i,:,:),roots,field)
    do component=1, 3
      work = cmplx(field(:,component),0.0_c_double,kind=c_double)
      call native_fft(work,roots,.false.)
      fourier(i,:,component) = work(0:lm)
    end do
  end do
  !$omp end do
  deallocate(field,work)
  !$omp end parallel
  call analyze_fourier(lm,theta,weights,fourier,joint,recovered,degree_energy)
  statistics(1) = sum(weights)
  statistics(2) = maxval(abs(recovered-coefficients))/maxval(abs(coefficients))
  statistics(3) = sqrt(sum((recovered-coefficients)**2)/sum(coefficients**2))
end subroutine native_roundtrip_probe

subroutine native_selected_probe(lm,m,theta,values) bind(c)
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use gravmag_sphere_harmonics
  implicit none
  integer(c_int), value :: lm,m
  real(c_double), value :: theta
  real(c_double), intent(out) :: values(3)
  real(c_double), allocatable :: seeds(:,:)
  integer(c_int), allocatable :: seed_powers(:,:)
  real(c_double) :: st(1),ct(1),previous(1,3),current(1,3),result(1,3)
  integer(c_int) :: l,powers(1)
  st = sin(theta)
  ct = cos(theta)
  call schmidt_seeds(lm,[theta],seeds,seed_powers)
  call schmidt_start(m,st,ct,seeds(:,m),seed_powers(:,m),previous,current,powers)
  do l=m+1, lm
    call schmidt_advance(l,m,st,ct,previous,current,powers)
  end do
  call schmidt_values(current,powers,result)
  values = result(1,:)
end subroutine native_selected_probe

subroutine native_fft_probe(n,values,inverse) bind(c)
  use, intrinsic :: iso_c_binding, only: c_int,c_double_complex,c_double
  use gravmag_sphere_transform
  implicit none
  integer(c_int), value :: n,inverse
  complex(c_double_complex), intent(inout) :: values(0:n-1)
  complex(c_double) :: roots(0:n/2-1)
  call fft_roots(n,roots)
  call native_fft(values,roots,inverse /= 0)
end subroutine native_fft_probe

subroutine native_nodes_probe(n,theta,weights) bind(c)
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use gravmag_sphere_transform
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(out) :: theta(n),weights(n)
  call gauss_nodes(n,theta,weights)
end subroutine native_nodes_probe

subroutine native_point_probe(lm,theta,phi,coefficients,field) bind(c)
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use gravmag_sphere_harmonics, only: schmidt_point
  implicit none
  integer(c_int), value :: lm
  real(c_double), value :: theta,phi
  real(c_double), intent(in) :: coefficients(lm*(lm+2))
  real(c_double), intent(out) :: field(3)
  call schmidt_point(lm,theta,phi,coefficients,field)
end subroutine native_point_probe
