module gravmag_orbital
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none
contains
  subroutine orbital_field(nobs, nbody, nq, radius_km, xyz_km, bodies, nodes, weights, field_nt, status) &
      bind(c, name="orbital_field")
    ! integrate constant cartesian magnetization over spherical block volumes
    ! coordinates and source depths use km, magnetization a/m, field output nt
    integer(c_int), value, intent(in) :: nobs, nbody, nq
    real(c_double), value, intent(in) :: radius_km
    real(c_double), intent(in) :: xyz_km(3,nobs), bodies(9,nbody), nodes(nq), weights(nq)
    real(c_double), intent(out) :: field_nt(3,nobs)
    integer(c_int), intent(out) :: status
    integer :: ibody, ilat, ilon, ir, iobs
    real(c_double) :: pi, radians, lat_half, lon_half, radius_half, radius_mid
    real(c_double) :: lat, lon, radial, volume, source(3), moment(3), delta(3), distance2, inverse3, dot

    pi = acos(-1.0_c_double)
    radians = pi / 180.0_c_double
    status = 0
    field_nt = 0.0_c_double
    if (nq < 1 .or. nobs < 1 .or. nbody < 1 .or. radius_km <= 0.0_c_double) then
      status = 1
      return
    end if
    do ibody = 1, nbody
      ! body order: latitude, longitude, latitude width, longitude width,
      ! top depth, thickness, mx, my, mz
      if (any(bodies(3:4,ibody) <= 0.0_c_double) .or. bodies(5,ibody) < 0.0_c_double .or. &
          bodies(6,ibody) <= 0.0_c_double .or. sum(bodies(5:6,ibody)) >= radius_km) then
        status = 1
        return
      end if
      lat_half = bodies(3,ibody) * radians * 0.5_c_double
      lon_half = bodies(4,ibody) * radians * 0.5_c_double
      radius_half = bodies(6,ibody) * 500.0_c_double
      radius_mid = (radius_km - bodies(5,ibody) - bodies(6,ibody) * 0.5_c_double) * 1000.0_c_double
      do ilat = 1, nq
        lat = bodies(1,ibody) * radians + lat_half * nodes(ilat)
        do ilon = 1, nq
          lon = bodies(2,ibody) * radians + lon_half * nodes(ilon)
          do ir = 1, nq
            radial = radius_mid + radius_half * nodes(ir)
            source = radial * [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]
            volume = radial**2 * cos(lat) * lat_half * lon_half * radius_half &
                     * weights(ilat) * weights(ilon) * weights(ir)
            moment = bodies(7:9,ibody) * volume
            do iobs = 1, nobs
              delta = xyz_km(:,iobs) * 1000.0_c_double - source
              distance2 = sum(delta**2)
              if (distance2 <= 1.0e-12_c_double) then
                status = 2
                return
              end if
              inverse3 = 1.0_c_double / (distance2 * sqrt(distance2))
              dot = sum(moment * delta)
              ! mu0/(4*pi) times tesla-to-nt is 100
              field_nt(:,iobs) = field_nt(:,iobs) + 100.0_c_double * inverse3 &
                                * (3.0_c_double * dot * delta / distance2 - moment)
            end do
          end do
        end do
      end do
    end do
  end subroutine orbital_field

  subroutine dipole_field(nobs, nsource, xyz_km, source_xyz_km, moments_am2, field_nt, status) &
      bind(c, name="dipole_field")
    ! evaluate fitted point-dipole moments at new map or orbital positions
    ! use the same dipole kernel and si moment units as the equivalent-grid fitter
    integer(c_int), value, intent(in) :: nobs, nsource
    real(c_double), intent(in) :: xyz_km(3,nobs), source_xyz_km(3,nsource), moments_am2(3,nsource)
    real(c_double), intent(out) :: field_nt(3,nobs)
    integer(c_int), intent(out) :: status
    integer :: iobs, isource
    real(c_double) :: delta(3), distance2, inverse3, dot

    status = 0
    field_nt = 0.0_c_double
    if (nobs < 1 .or. nsource < 1) then
      status = 1
      return
    end if
    do isource = 1, nsource
      do iobs = 1, nobs
        delta = (xyz_km(:,iobs) - source_xyz_km(:,isource)) * 1000.0_c_double
        distance2 = sum(delta**2)
        if (distance2 <= 1.0e-12_c_double) then
          status = 2
          return
        end if
        inverse3 = 1.0_c_double / (distance2 * sqrt(distance2))
        dot = sum(moments_am2(:,isource) * delta)
        field_nt(:,iobs) = field_nt(:,iobs) + 100.0_c_double * inverse3 &
                          * (3.0_c_double * dot * delta / distance2 - moments_am2(:,isource))
      end do
    end do
  end subroutine dipole_field
  subroutine monopole_kernel(nobs, nsource, xyz_km, source_xyz_km, kernel, status) &
      bind(c, name="monopole_kernel")
    ! Equivalent magnetic charges q in nT km^2: B=q*(x-s)/|x-s|^3.
    ! Kernel layout matches a C-contiguous Python (nobs,3,nsource) array.
    integer(c_int), value, intent(in) :: nobs, nsource
    real(c_double), intent(in) :: xyz_km(3,nobs), source_xyz_km(3,nsource)
    real(c_double), intent(out) :: kernel(nsource,3,nobs)
    integer(c_int), intent(out) :: status
    integer :: i, j
    real(c_double) :: delta(3), distance2
    status = 0
    kernel = 0.0_c_double
    if (nobs < 1 .or. nsource < 1) then
      status = 1
      return
    end if
    do i = 1, nobs
      do j = 1, nsource
        delta = xyz_km(:,i) - source_xyz_km(:,j)
        distance2 = sum(delta**2)
        if (distance2 <= 1.0e-18_c_double) then
          status = 2
          return
        end if
        kernel(j,:,i) = delta / (distance2 * sqrt(distance2))
      end do
    end do
  end subroutine monopole_kernel
end module gravmag_orbital
