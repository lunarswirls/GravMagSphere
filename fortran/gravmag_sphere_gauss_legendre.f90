program gravmag_sphere_gauss_legendre
  use, intrinsic :: iso_fortran_env, only: real64, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use gravmag_paths, only: default_output_path
  implicit none

  ! restore the nested longitude, latitude, radius gauss-legendre volume method
  ! legacy reference: gravmag_sphere_brtp.f at git revision 71d13a1
  ! von frese et al (1981), journal of geophysics 49, 234-242; ravat (1989)
  ! use current cards, si properties, planet-fixed xyz vectors and real64
  ! polygon longitude slabs end at vertices so boundary slopes stay continuous

  real(real64), parameter :: pi = acos(-1.0_real64), radians = pi / 180.0_real64
  real(real64), parameter :: gravity = 6.67430e-6_real64
  character(len=4096) :: infile, outfile, line
  character(len=32) :: diag_env
  integer :: input_unit, output_unit, ios, body, narg, index, subdivisions
  integer :: nr, ntheta, nphi, nblim, ifield, nrem, nlat, nlon, npts
  integer :: radial_order, latitude_order, longitude_order
  integer :: ilat, ilon, ir, slab, span, ip, jp, kp, i, j, count, next
  integer :: iprint, nfile, nopt
  real(real64) :: radius_km, lat0, lon0, dlat, dlon, altitude, amplitude, inc, dec
  real(real64) :: unused(5), first, conint, top, bottom, latlo, lathi, lonlo, lonhi
  real(real64) :: lon, lat, radial, lonhalf, lonmid, lathalf, latmid, rhalf, rmid
  real(real64) :: volume, total_volume, weight, moment(3), position(3), delta(3), distance2
  real(real64) :: contribution(3), corrected(3), updated(3), inverse3, factor, swap
  real(real64) :: start_time, finish_time
  real(real64), allocatable :: poly_lat(:), poly_lon(:), breaks(:), crossings(:)
  real(real64), allocatable :: rn(:), rw(:), tn(:), tw(:), pn(:), pw(:)
  real(real64), allocatable :: observation(:,:,:), field(:,:,:), compensation(:,:,:)
  integer(int64) :: source_nodes
  logical :: diagnostics

  narg = command_argument_count()
  if (narg < 2 .or. narg > 7) then
    write(*,'(a)') 'usage: gravmag_sphere_gauss_legendre radius_km input.in [output.txt] '// &
      '[radial_order] [latitude_order] [longitude_order] [subdivisions]'
    stop 2
  end if
  call get_command_argument(1, line)
  read(line,*,iostat=ios) radius_km
  if (ios /= 0) error stop 'error: invalid radius'
  if (.not. ieee_is_finite(radius_km) .or. radius_km <= 0) error stop 'error: invalid radius'
  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)
  radial_order = 0
  latitude_order = 0
  longitude_order = 0
  subdivisions = 1
  do index = 4, narg
    call get_command_argument(index, line)
    read(line,*,iostat=ios) i
    if (ios /= 0) error stop 'error: invalid quadrature option'
    if (i < 0 .or. i > 256) error stop 'error: quadrature options must be between 0 and 256'
    select case (index)
    case (4)
      radial_order = i
    case (5)
      latitude_order = i
    case (6)
      longitude_order = i
    case (7)
      subdivisions = i
    end select
  end do
  if (subdivisions < 1) error stop 'error: subdivisions must be positive'
  open(newunit=input_unit, file=trim(infile), status='old', iostat=ios)
  if (ios /= 0) error stop 'error: cannot open input'
  call default_output_path(trim(infile), outfile, '_quadrature.txt')
  open(newunit=output_unit, file=trim(outfile), status='replace', iostat=ios)
  if (ios /= 0) error stop 'error: cannot open output'
  write(output_unit,'(a)') '# solver=gauss_legendre'
  diag_env = ''
  call get_environment_variable('GRAVMAG_DIAGNOSTICS', diag_env)
  diagnostics = len_trim(diag_env) > 0 .and. trim(diag_env) /= '0'
  body = 0

  do
    call read_line(.true.)
    if (ios < 0) exit
    body = body + 1
    call cpu_time(start_time)
    call read_line()
    read(line,*,iostat=ios) lat0, lon0, dlat, dlon, altitude, nlat, nlon
    if (ios /= 0) error stop 'error: invalid card 2'
    if (.not. all(ieee_is_finite([lat0, lon0, dlat, dlon, altitude]))) error stop 'error: nonfinite grid'
    if (nlat < 1 .or. nlon < 1 .or. altitude < 0) error stop 'error: invalid observation grid'
    if (max(abs(lat0), abs(lat0 + (nlat-1)*dlat)) > 90) error stop 'error: invalid observation latitude'
    call read_line()
    read(line,*,iostat=ios) nr, ntheta, nphi, nblim
    if (ios /= 0) error stop 'error: invalid card 3'
    if (radial_order > 0) nr = radial_order
    if (latitude_order > 0) ntheta = latitude_order
    if (longitude_order > 0) nphi = longitude_order
    if (min(nr, ntheta, nphi) < 1 .or. max(nr, ntheta, nphi) > 256) &
      error stop 'error: quadrature orders must be between 1 and 256'
    call read_line()
    read(line,*,iostat=ios) unused
    if (ios /= 0) error stop 'error: invalid card 4'
    call read_line()
    read(line,*,iostat=ios) ifield, nrem, amplitude, inc, dec
    if (ios /= 0) error stop 'error: invalid card 5'
    if (.not. all(ieee_is_finite([amplitude, inc, dec]))) error stop 'error: nonfinite material'
    if (ifield /= 1 .and. ifield /= 2) error stop 'error: unsupported field mode'
    if (ifield == 2 .and. nrem /= 1) error stop 'error: magnetic mode requires nrem=1'
    moment = amplitude * [cos(inc*radians)*cos(dec*radians), cos(inc*radians)*sin(dec*radians), sin(inc*radians)]
    call read_line()
    read(line,*,iostat=ios) iprint, nfile, nopt, first, conint
    if (ios /= 0) error stop 'error: invalid card 6'
    call read_line()
    if (nblim == 1) then
      read(line,*,iostat=ios) lathi, latlo, lonhi, lonlo, top, bottom
      if (ios /= 0) error stop 'error: invalid block limits'
      if (lathi <= latlo .or. lonhi <= lonlo) error stop 'error: block limits must increase'
      npts = 4
      allocate(poly_lat(npts), poly_lon(npts))
      poly_lat = [latlo, latlo, lathi, lathi]
      poly_lon = [lonlo, lonhi, lonhi, lonlo]
    else if (nblim == 0) then
      read(line,*,iostat=ios) npts, top, bottom
      if (ios /= 0) error stop 'error: invalid polygon header'
      if (npts < 3) error stop 'error: polygon requires at least three vertices'
      allocate(poly_lat(npts), poly_lon(npts))
      do i = 1, npts
        call read_line()
        read(line,*,iostat=ios) poly_lat(i), poly_lon(i)
        if (ios /= 0) error stop 'error: invalid polygon vertex'
        if (abs(poly_lat(i)) > 90 .and. abs(poly_lon(i)) <= 90) then
          swap = poly_lat(i)
          poly_lat(i) = poly_lon(i)
          poly_lon(i) = swap
        end if
      end do
    else
      error stop 'error: nblim must be 0 or 1'
    end if
    if (.not. all(ieee_is_finite([top, bottom, poly_lat, poly_lon]))) error stop 'error: nonfinite geometry'
    if (top < 0 .or. bottom <= top .or. bottom >= radius_km) error stop 'error: invalid depths'
    if (altitude + top <= 0) error stop 'error: observation sphere must be above the source top'
    if (any(abs(poly_lat) >= 90)) error stop 'error: source must avoid the poles'
    if (nblim == 0) then
      do i = 2, npts
        poly_lon(i) = poly_lon(i-1) + modulo(poly_lon(i)-poly_lon(i-1)+180.0_real64, 360.0_real64)-180.0_real64
      end do
    end if
    if (maxval(poly_lon)-minval(poly_lon) >= 180) error stop 'error: split sources spanning 180 degrees longitude'
    ! match straight edges in the unwrapped longitude/latitude card coordinates
    allocate(breaks(npts), crossings(npts))
    breaks = poly_lon
    call sort_values(breaks)
    allocate(rn(nr), rw(nr), tn(ntheta), tw(ntheta), pn(nphi), pw(nphi))
    call gauss_legendre(rn, rw)
    call gauss_legendre(tn, tw)
    call gauss_legendre(pn, pw)
    allocate(observation(3,nlon,nlat), field(3,nlon,nlat), compensation(3,nlon,nlat))
    do i = 1, nlat
      lat = (lat0 + (i-1)*dlat)*radians
      do j = 1, nlon
        lon = (lon0 + (j-1)*dlon)*radians
        observation(:,j,i) = (radius_km+altitude)*1000 * [cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)]
      end do
    end do
    field = 0
    compensation = 0
    source_nodes = 0
    total_volume = 0
    ! integrate each longitude slab and each pair of polygon crossings
    do slab = 1, npts-1
      if (breaks(slab+1) <= breaks(slab)) cycle
      lonhalf = (breaks(slab+1)-breaks(slab))/(2*subdivisions)
      do ip = 1, subdivisions
        lonmid = breaks(slab) + (2*ip-1)*lonhalf
        do ilon = 1, nphi
          lon = lonmid + lonhalf*pn(ilon)
          count = 0
          do i = 1, npts
            next = mod(i,npts)+1
            if ((poly_lon(i) > lon) .neqv. (poly_lon(next) > lon)) then
              count = count+1
              crossings(count) = poly_lat(i) + (poly_lat(next)-poly_lat(i)) * &
                (lon-poly_lon(i))/(poly_lon(next)-poly_lon(i))
            end if
          end do
          if (mod(count,2) /= 0) error stop 'error: unpaired polygon crossings'
          call sort_values(crossings(:count))
          do span = 1, count, 2
            lathalf = (crossings(span+1)-crossings(span))/(2*subdivisions)
            if (lathalf <= 0) cycle
            do jp = 1, subdivisions
              latmid = crossings(span) + (2*jp-1)*lathalf
              do ilat = 1, ntheta
                lat = (latmid+lathalf*tn(ilat))*radians
                position = [cos(lat)*cos(lon*radians), cos(lat)*sin(lon*radians), sin(lat)]
                weight = lonhalf*lathalf*radians**2*pw(ilon)*tw(ilat)*cos(lat)
                rhalf = (bottom-top)*500/subdivisions
                do kp = 1, subdivisions
                  rmid = (radius_km-bottom)*1000 + (2*kp-1)*rhalf
                  do ir = 1, nr
                    radial = rmid+rhalf*rn(ir)
                    volume = radial**2*rhalf*rw(ir)*weight
                    total_volume = total_volume+volume
                    source_nodes = source_nodes+1
                    do i = 1, nlat
                      do j = 1, nlon
                        delta = observation(:,j,i)-radial*position
                        distance2 = dot_product(delta,delta)
                        if (distance2 <= 1e-12_real64) error stop 'error: singular quadrature node'
                        inverse3 = 1/(distance2*sqrt(distance2))
                        if (ifield == 2) then
                          factor = 3*dot_product(moment,delta)/distance2
                          contribution = 100*volume*inverse3*(factor*delta-moment)
                        else
                          contribution = -gravity*amplitude*volume*inverse3*delta
                        end if
                        ! compensated accumulation protects cancellation between volume elements
                        corrected = contribution-compensation(:,j,i)
                        updated = field(:,j,i)+corrected
                        compensation(:,j,i) = (updated-field(:,j,i))-corrected
                        field(:,j,i) = updated
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    if (source_nodes == 0) error stop 'error: polygon has zero volume'
    if (.not. all(ieee_is_finite(field))) error stop 'error: nonfinite field'
    if (ifield == 2) then
      write(output_unit,'(a)') '# body_id lon_deg lat_deg bx_nT by_nT bz_nT btot_nT'
    else
      write(output_unit,'(a)') '# body_id lon_deg lat_deg gx_mGal gy_mGal gz_mGal gtot_mGal'
    end if
    do i = 1, nlat
      do j = 1, nlon
        lon = modulo(lon0+(j-1)*dlon+180.0_real64,360.0_real64)-180.0_real64
        write(output_unit,'(i0,6(1x,es24.16))') body, lon, lat0+(i-1)*dlat, &
          field(:,j,i), sqrt(sum(field(:,j,i)**2))
      end do
    end do
    call cpu_time(finish_time)
    if (diagnostics) then
      write(*,'(a,i0)') 'DIAG|gauss_legendre|meta|body_id=', body
      write(*,'(a,i0)') 'DIAG|gauss_legendre|meta|source_nodes=', source_nodes
      write(*,'(a,i0)') 'DIAG|gauss_legendre|meta|radial_order=', nr
      write(*,'(a,i0)') 'DIAG|gauss_legendre|meta|latitude_order=', ntheta
      write(*,'(a,i0)') 'DIAG|gauss_legendre|meta|longitude_order=', nphi
      write(*,'(a,i0)') 'DIAG|gauss_legendre|meta|subdivisions=', subdivisions
      write(*,'(a,es24.16)') 'DIAG|gauss_legendre|meta|volume_m3=', total_volume
      write(*,'(a,f14.6)') 'DIAG|gauss_legendre|time|body_total_s=', finish_time-start_time
    end if
    deallocate(poly_lat, poly_lon, breaks, crossings, rn, rw, tn, tw, pn, pw, observation, field, compensation)
  end do
  if (body == 0) error stop 'error: no bodies in input'
  close(input_unit)
  close(output_unit)

contains

  subroutine read_line(allow_eof)
    logical, optional, intent(in) :: allow_eof
    do
      read(input_unit,'(a)',iostat=ios) line
      if (ios /= 0) then
        if (present(allow_eof)) then
          if (allow_eof .and. ios < 0) return
        end if
        error stop 'error: incomplete input cards'
      end if
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) /= '#' .and. line(1:1) /= '!') return
    end do
  end subroutine read_line

  subroutine sort_values(values)
    real(real64), intent(inout) :: values(:)
    real(real64) :: value
    integer :: a, b
    do a = 2, size(values)
      value = values(a)
      b = a-1
      do while (b >= 1)
        if (values(b) <= value) exit
        values(b+1) = values(b)
        b = b-1
      end do
      values(b+1) = value
    end do
  end subroutine sort_values

  subroutine gauss_legendre(nodes, weights)
    real(real64), intent(out) :: nodes(:), weights(:)
    real(real64) :: root, step, p0, p1, p2, derivative
    integer :: n, a, b, iteration
    n = size(nodes)
    do a = 1, (n+1)/2
      root = cos(pi*(a-0.25_real64)/(n+0.5_real64))
      do iteration = 1, 100
        p0 = 1
        p1 = root
        do b = 2, n
          p2 = ((2*b-1)*root*p1-(b-1)*p0)/b
          p0 = p1
          p1 = p2
        end do
        derivative = n*(root*p1-p0)/(root**2-1)
        step = p1/derivative
        if (abs(step) < 4*epsilon(root)) exit
        root = root-step
      end do
      if (iteration > 100) error stop 'error: quadrature roots did not converge'
      nodes(a) = -root
      nodes(n+1-a) = root
      weights(a) = 2/((1-root**2)*derivative**2)
      weights(n+1-a) = weights(a)
    end do
  end subroutine gauss_legendre
end program gravmag_sphere_gauss_legendre
