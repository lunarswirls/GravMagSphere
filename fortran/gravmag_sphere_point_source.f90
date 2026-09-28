program gravmag_sphere_point_source
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use gravmag_sphere_subs
  use gravmag_sphere_physics
  use gravmag_paths, only: default_output_path
  implicit none

  !***********************************************************************
  ! gravmag_sphere_point_source
  !
  ! this program computes either magnetic or gravity field components on a
  ! spherical observation grid due to a 3-d body below the surface and
  ! output field is saved in the cartesian xyz basis
  !
  ! adapted from the legacy "sphere" family comments/instructions:
  ! - repeated body blocks are supported to build total multi-body response
  ! - card 2 (observation grid) should be identical across bodies when
  !     summing bodies in one output file
  ! - numerical integration controls determine cost and accuracy
  !
  ! numerical method:
  ! - magnetic fields sum point charges on discretized top, bottom and side surfaces
  ! - gravity fields sum point masses from midpoint volume cells
  ! - nr controls radial sampling and refine_factor scales horizontal sampling
  !
  ! mode selection (card 5):
  ! - ifield=2: magnetic mode
  !   fields: bx/by/bz [nT], btot [nT]
  !   card5: ifield nrem m_amp[a/m] inc[deg] dec[deg], with nrem=1
  ! - ifield=1: gravity mode
  !   fields: gx/gy/gz [mgal], gtot [mgal]
  !   card5: ifield nrem density_contrast[kg/m^3] dummy dummy
  !
  ! geometry (card 7):
  ! - fixed limits: lat_max lat_min lon_max lon_min depth_top depth_bot
  ! - polygon: npts depth_top depth_bot + npts vertices
  !   vertices may be provided as (lat lon), or (lon lat) when one value
  !   is clearly outside latitude limits.
  !
  ! REFERENCES: (in publications, please give appropriate credit to following
  !             individuals for their
  !             efforts in developing these programs)
  !
  ! von Frese, R.R.B., W.J. Hinze, L.W. Braile, and A.J. Luca, 1981,
  !     Spherical-sphere Gravity and Magnetic Anomaly Modeling by Gauss-
  !     Legendre Quadrature Integration, Journal of Geophysics.,
  !     vol.49,pp.234-242.
  !
  ! Ravat, D., 1989, Magsat Investigations over the Greater African Region,
  !     Ph.D. Dissertation, Purdue University, West Lafayette, IN, 234p.
  !***********************************************************************

  ! ---- constants ----
  real(real32), parameter :: pi      = 3.14159265358979323846_real32
  real(real32), parameter :: deg2rad = pi/180.0_real32
  real(real32), parameter :: T_to_nT = 1.0e9_real32
  real(real32), parameter :: mps2_to_mgal = 1.0e5_real32

  ! ---- CLI ----
  real(real32) :: rsphere_km
  character(len=256) :: infile, outfile, tmp
  integer(int32) :: ios, narg

  ! ---- per-body cards ----
  character(len=256) :: title
  integer(int32) :: knk

  ! Card 2
  real(real32) :: lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km
  integer(int32) :: nlat, nlon

  ! Card 3:
  ! legacy: nr, ntheta, nphi gaussian nodes + nblim geometry flag.
  ! current in this code:
  ! - nr is used as radial midpoint count
  ! - ntheta, nphi are read for compatibility and provenance
  ! - nblim controls fixed limits (1) vs polygon (0)
  integer(int32) :: nr, ntheta, nphi, nblim

  ! Card 4 (kept for backwards compatibility, but unused now)
  real(real32) :: htheta, hphi, phi1, phi2, rho_dummy

  ! Card 5
  integer(int32) :: ifield, nrem
  real(real32) :: M_amp_Apm, inc_deg, dec_deg
  real(real32) :: rho_kgm3

  ! Card 6
  integer(int32) :: iprint, nfile, nopt
  real(real32) :: first, conint

  ! Card 7 (source geometry)
  integer(int32) :: npts
  real(real32) :: depth_top_km, depth_bot_km
  real(real32) :: lat_max, lat_min, lon_max, lon_min
  real(real32), allocatable :: poly_lat(:), poly_lon(:)

  ! ---- derived geometry ----
  real(real32) :: rsphere_m, ro_m
  real(real32) :: r_top_m, r_bot_m, dr_m

  ! magnetization vector (planet-fixed XYZ)
  real(real32) :: cosI, sinI, cosD, sinD
  real(real32) :: Mx, My, Mz

  ! observation grid xyz (meters)
  real(real32), allocatable :: xo(:,:), yo(:,:), zo(:,:)

  ! outputs (Tesla)
  real(real32), allocatable :: Bx(:,:), By(:,:), Bz(:,:)

  ! ---- rasterization controls ----
  ! more nodes improve accuracy but cost more
  ! if convergence is poor, increase resolution or split bodies
  integer(int32) :: refine_factor
  integer(int32) :: src_nlat_cli, src_nlon_cli, src_nr_cli
  integer(int32) :: nlat_s, nlon_s, nlat_seed, nlon_seed, nr_s
  real(real32) :: dlat_s_deg, dlon_s_deg
  real(real32) :: lat_lo, lat_hi, lon_lo, lon_hi, lon_ref

  ! source list (dipole moments)
  integer(int32) :: ns, cap
  real(real32), allocatable :: xs_list(:), ys_list(:), zs_list(:)
  real(real32), allocatable :: qmag_list(:)
  real(real32), allocatable :: mass_list(:)

  ! loop indices
  integer(int32) :: i, j, is, it, ir
  integer(int32) :: ie, ie2, nedge, nseg, ks
  real(real32) :: lat_deg, lon_deg, lat_rad, lon_rad
  real(real32) :: xs, ys, zs, dV, r_mid, w_r
  real(real32) :: edge_deg, denom_deg, winding
  real(real64) :: signed_area
  real(real32) :: t0, t1, lat_a, lon_a, lat_b, lon_b, r0, r1
  real(real32) :: erx, ery, erz, sigma
  real(real32) :: dBx, dBy, dBz
  real(real32) :: btot_nt
  real(real32) :: gtot_mgal
  real(real32) :: p1, p2

  ! per-thread accumulators
  real(real32) :: bx_acc, by_acc, bz_acc
  logical :: diag_enabled
  character(len=32) :: diag_env
  real :: t_body0, t_body1, t_phase0, t_phase1
  real :: t_obsgrid_s, t_source_s, t_accum_s, t_output_s
  integer(int64) :: obs_cells64, source_cap64
  integer(int64) :: mem_obs_bytes, mem_source_bytes, mem_total_bytes

  ! ----------------------------
  ! CLI
  ! ----------------------------
  narg = command_argument_count()
  if (narg < 2) then
    write(*,*) "Usage: gravmag_sphere_point_source <R_sphere_km> <input.in> [output.txt] "// &
               "[refine_factor] [source_nlat] [source_nlon] [source_nr]"
    stop 2
  end if
  call get_command_argument(1, tmp)
  read(tmp,*,iostat=ios) rsphere_km
  if (ios /= 0) stop "ERROR: could not parse R_sphere_km"

  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)

  open(4, file=trim(infile), status='old', iostat=ios)
  if (ios /= 0) stop "ERROR: cannot open input file"

  call default_output_path(trim(infile), outfile, '.txt')
  open(2, file=trim(outfile), status='replace', iostat=ios)
  if (ios /= 0) stop "ERROR: cannot open output file"

  rsphere_m = rsphere_km * 1000.0_real32
  knk = 0

  refine_factor = 2   ! default
  src_nlat_cli = 0
  src_nlon_cli = 0
  src_nr_cli = 0
  if (narg >= 4) then
    call get_command_argument(4, tmp)
    read(tmp,*,iostat=ios) refine_factor
    if (ios /= 0 .or. refine_factor < 1) stop "ERROR: invalid refine_factor CLI argument"
  end if
  if (narg >= 5) then
    call get_command_argument(5, tmp)
    read(tmp,*,iostat=ios) src_nlat_cli
    if (ios /= 0 .or. src_nlat_cli < 0) stop "ERROR: invalid source_nlat CLI argument"
  end if
  if (narg >= 6) then
    call get_command_argument(6, tmp)
    read(tmp,*,iostat=ios) src_nlon_cli
    if (ios /= 0 .or. src_nlon_cli < 0) stop "ERROR: invalid source_nlon CLI argument"
  end if
  if (narg >= 7) then
    call get_command_argument(7, tmp)
    read(tmp,*,iostat=ios) src_nr_cli
    if (ios /= 0 .or. src_nr_cli < 0) stop "ERROR: invalid source_nr CLI argument"
  end if

  diag_env = ''
  call get_environment_variable('GRAVMAG_DIAGNOSTICS', diag_env)
  diag_enabled = (len_trim(diag_env) > 0 .and. trim(diag_env) /= '0')

  do
    ! ---- read title (skip blanks/comments) ----
    title = ''
    do
      read(4,'(A)', iostat=ios) title
      if (ios /= 0) exit
      if (len_trim(title) == 0) cycle
      if (title(1:1) == '!' .or. title(1:1) == '#') cycle
      exit
    end do
    if (ios /= 0) exit

    knk = knk + 1
    write(*,'(A,I0,2A)') 'Body ', knk, ': ', trim(title)
    t_obsgrid_s = 0.0
    t_source_s = 0.0
    t_accum_s = 0.0
    t_output_s = 0.0
    call cpu_time(t_body0)

    ! ---- Cards 2-6 ----
    read(4,*,iostat=ios) lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km, nlat, nlon
    if (ios /= 0) stop "ERROR reading Card 2"

    read(4,*,iostat=ios) nr, ntheta, nphi, nblim
    if (ios /= 0) stop "ERROR reading Card 3"

    read(4,*,iostat=ios) htheta, hphi, phi1, phi2, rho_dummy
    if (ios /= 0) stop "ERROR reading Card 4"

    read(4,*,iostat=ios) ifield, nrem, M_amp_Apm, inc_deg, dec_deg
    if (ios /= 0) stop "ERROR reading Card 5"
    if (ifield /= 1 .and. ifield /= 2) stop "ERROR: Card 5 ifield must be 1(gravity) or 2(magnetic)"
    if (ifield == 2 .and. nrem /= 1) stop "ERROR: magnetics require nrem=1"

    read(4,*,iostat=ios) iprint, nfile, nopt, first, conint
    if (ios /= 0) stop "ERROR reading Card 6"

    if (ifield == 2) then
      ! Magnetic mode:
      !   Card 5 values: ifield nrem M_amp[A/m] inc[deg] dec[deg]
      cosI = cos(inc_deg*deg2rad)
      sinI = sin(inc_deg*deg2rad)
      cosD = cos(dec_deg*deg2rad)
      sinD = sin(dec_deg*deg2rad)

      Mx = M_amp_Apm * (cosI*cosD)
      My = M_amp_Apm * (cosI*sinD)
      Mz = M_amp_Apm * (sinI)
      rho_kgm3 = 0.0_real32
    else
      ! Gravity mode:
      !   Card 5 values: ifield nrem density_contrast[kg/m^3] dummy dummy
      Mx = 0.0_real32
      My = 0.0_real32
      Mz = 0.0_real32
      rho_kgm3 = M_amp_Apm
    end if

    ! ---- geometry Card 7 ----
    if (allocated(poly_lat)) deallocate(poly_lat, poly_lon)

    if (nblim == 1) then
      ! lat_max lat_min lon_max lon_min depth_top_km depth_bot_km
      read(4,*,iostat=ios) lat_max, lat_min, lon_max, lon_min, depth_top_km, depth_bot_km
      if (ios /= 0) stop "ERROR reading fixed limits (Card 7) - expected 6 values"
      npts = 5
      allocate(poly_lat(npts), poly_lon(npts))
      poly_lat = [lat_min, lat_max, lat_max, lat_min, lat_min]
      poly_lon = [lon_min, lon_min, lon_max, lon_max, lon_min]
    else
      read(4,*,iostat=ios) npts, depth_top_km, depth_bot_km
      if (ios /= 0) stop "ERROR reading polygon header: npts depth_top depth_bot"
      if (npts < 3) stop "ERROR: polygon needs >=3 vertices"
      allocate(poly_lat(npts), poly_lon(npts))
      do i=1,npts
        read(4,*,iostat=ios) p1, p2
        if (ios /= 0) stop "ERROR reading polygon vertex"
        ! Accept both (lat lon) and (lon lat):
        ! - if one value is outside valid latitude range, treat that as longitude
        ! - if both are within [-90,90], default to legacy (lat lon)
        if (abs(p1) <= 90.0_real32 .and. abs(p2) > 90.0_real32) then
          poly_lat(i) = p1
          poly_lon(i) = p2
        else if (abs(p1) > 90.0_real32 .and. abs(p2) <= 90.0_real32) then
          poly_lon(i) = p1
          poly_lat(i) = p2
        else
          poly_lat(i) = p1
          poly_lon(i) = p2
        end if
      end do
    end if

    ! unwrap polygon longitudes near mean
    lon_ref = sum(poly_lon) / real(size(poly_lon), real32)
    call wrap180_scalar(lon_ref)
    call unwrap_lons_inplace(poly_lon, lon_ref)

    lat_lo = minval(poly_lat); lat_hi = maxval(poly_lat)
    lon_lo = minval(poly_lon); lon_hi = maxval(poly_lon)

    if (depth_top_km < 0.0_real32 .or. depth_bot_km < 0.0_real32) stop "ERROR: depth must be >= 0"
    if (depth_bot_km <= depth_top_km) stop "ERROR: need depth_bot > depth_top"

    r_top_m = rsphere_m - depth_top_km*1000.0_real32
    r_bot_m = rsphere_m - depth_bot_km*1000.0_real32
    if (r_top_m <= 0.0_real32 .or. r_bot_m <= 0.0_real32) stop "ERROR: depth too large (radius <= 0)"

    if (r_bot_m > r_top_m) then
      ! swap
      dr_m   = r_bot_m
      r_bot_m = r_top_m
      r_top_m = dr_m
    end if
    dr_m = r_top_m - r_bot_m

    ro_m = (rsphere_km + elvo_km) * 1000.0_real32

    ! ---- allocate observation arrays ----
    call cpu_time(t_phase0)
    if (.not. allocated(xo)) then
      allocate(xo(nlat,nlon), yo(nlat,nlon), zo(nlat,nlon))
      allocate(Bx(nlat,nlon), By(nlat,nlon), Bz(nlat,nlon))
    else
      if (size(xo,1) /= nlat .or. size(xo,2) /= nlon) then
        deallocate(xo,yo,zo,Bx,By,Bz)
        allocate(xo(nlat,nlon), yo(nlat,nlon), zo(nlat,nlon))
        allocate(Bx(nlat,nlon), By(nlat,nlon), Bz(nlat,nlon))
      end if
    end if

    ! Build obs xyz from grid definition (lon in [-180,180))
    do i=1,nlat
      do j=1,nlon
        lat_deg = lat0_deg + real(i-1,real32)*dlat_deg
        lon_deg = lon0_deg + real(j-1,real32)*dlon_deg
        call wrap180_scalar(lon_deg)

        lat_rad = lat_deg*deg2rad
        lon_rad = lon_deg*deg2rad

        xo(i,j) = ro_m * cos(lat_rad) * cos(lon_rad)
        yo(i,j) = ro_m * cos(lat_rad) * sin(lon_rad)
        zo(i,j) = ro_m * sin(lat_rad)
      end do
    end do
    call cpu_time(t_phase1)
    t_obsgrid_s = t_obsgrid_s + (t_phase1 - t_phase0)

    ! ------------------------------------------
    ! rasterize polygon on fine grid
    ! ------------------------------------------
    ! Source mesh is now decoupled from Card 2 output grid:
    ! - default source controls come from Card 3 (ntheta, nphi, nr)
    ! - optional CLI overrides can replace them globally
    call cpu_time(t_phase0)
    nlat_seed = ntheta
    nlon_seed = nphi
    nr_s = nr
    if (src_nlat_cli > 0) nlat_seed = src_nlat_cli
    if (src_nlon_cli > 0) nlon_seed = src_nlon_cli
    if (src_nr_cli > 0) nr_s = src_nr_cli
    if (nlat_seed <= 0) nlat_seed = max(1_int32, nlat)
    if (nlon_seed <= 0) nlon_seed = max(1_int32, nlon)
    if (nr_s <= 0) nr_s = 1

    nlat_s = max(1_int32, refine_factor * nlat_seed)
    nlon_s = max(1_int32, refine_factor * nlon_seed)

    dlat_s_deg = (lat_hi - lat_lo) / real(nlat_s, real32)
    dlon_s_deg = (lon_hi - lon_lo) / real(nlon_s, real32)

    ! source list capacity estimate (expanded as needed)
    cap = max(4096_int32, 4*nlat_s*nlon_s + 32*nr_s*max(1_int32,refine_factor))
    ns  = 0
    if (allocated(xs_list)) deallocate(xs_list, ys_list, zs_list, qmag_list, mass_list)
    allocate(xs_list(cap), ys_list(cap), zs_list(cap))
    allocate(qmag_list(cap), mass_list(cap))

    write(*,'(A,I0,A,I0,A,I0)') '  Source mesh cells (nlat,nlon,nr): ', nlat_s, ', ', nlon_s, ', ', nr_s

    if (ifield == 2) then
      ! ------------------------------------------------------------------
      ! Magnetic mode: build equivalent magnetic surface charges
      ! q = (M · n_out) dS, and dB = mu0/(4pi) * q * R / |R|^3
      !
      ! Surfaces:
      ! 1) top spherical cap (r = r_top)
      ! 2) bottom spherical cap (r = r_bot, outward normal inward)
      ! 3) polygon side wall (triangulated edge x radial strips)
      ! ------------------------------------------------------------------
      ! Top/bottom charge patches over polygon interior
      do it=1, nlat_s
        lat_deg = lat_lo + (real(it,real32)-0.5_real32)*dlat_s_deg
        if (lat_deg <= -90.0_real32 .or. lat_deg >= 90.0_real32) cycle

        do is=1, nlon_s
          lon_deg = lon_lo + (real(is,real32)-0.5_real32)*dlon_s_deg
          if (.not. point_in_poly_lonlat(lon_deg, lat_deg, poly_lon, poly_lat)) cycle

          lat_rad = lat_deg*deg2rad
          lon_rad = lon_deg*deg2rad

          erx = cos(lat_rad) * cos(lon_rad)
          ery = cos(lat_rad) * sin(lon_rad)
          erz = sin(lat_rad)

          ! Top surface (outward = +er)
          xs = r_top_m * erx
          ys = r_top_m * ery
          zs = r_top_m * erz
          dV = (r_top_m*r_top_m) * cos(lat_rad) * (dlat_s_deg*deg2rad) * (dlon_s_deg*deg2rad)
          sigma = Mx*erx + My*ery + Mz*erz
          call append_source(cap, ns, xs_list, ys_list, zs_list, qmag_list, mass_list, xs, ys, zs, sigma*dV, 0.0_real32)

          ! Bottom surface (outward = -er)
          xs = r_bot_m * erx
          ys = r_bot_m * ery
          zs = r_bot_m * erz
          dV = (r_bot_m*r_bot_m) * cos(lat_rad) * (dlat_s_deg*deg2rad) * (dlon_s_deg*deg2rad)
          sigma = -(Mx*erx + My*ery + Mz*erz)
          call append_source(cap, ns, xs_list, ys_list, zs_list, qmag_list, mass_list, xs, ys, zs, sigma*dV, 0.0_real32)
        end do
      end do

      if (ns == 0) stop "ERROR: no polygon interior cells found; increase refine_factor or check polygon"

      ! Side wall: triangulate edge x radial strips
      nedge = size(poly_lat)
      if (nedge > 1) then
        if (abs(poly_lat(1)-poly_lat(nedge)) < 1.0e-6_real32 .and. &
            abs(poly_lon(1)-poly_lon(nedge)) < 1.0e-6_real32) nedge = nedge - 1
      end if

      ! winding in unwrapped lon-lat coordinates fixes outward normals even at concave edges
      ! translate the shoelace sum to the first vertex to avoid cancellation for small footprints
      signed_area = 0.0_real64
      do ie=1, nedge
        ie2 = mod(ie, nedge) + 1
        signed_area = signed_area + &
          real(poly_lon(ie)-poly_lon(1),real64)*real(poly_lat(ie2)-poly_lat(1),real64) - &
          real(poly_lon(ie2)-poly_lon(1),real64)*real(poly_lat(ie)-poly_lat(1),real64)
      end do
      if (signed_area == 0.0_real64) stop "ERROR: polygon has zero signed area"
      winding = real(sign(1.0_real64, signed_area), real32)

      denom_deg = max(1.0e-4_real32, min(abs(dlat_s_deg), abs(dlon_s_deg)))

      do ie=1, nedge
        ie2 = ie + 1
        if (ie == nedge) ie2 = 1

        edge_deg = sqrt( (poly_lat(ie2)-poly_lat(ie))**2 + &
                         ((poly_lon(ie2)-poly_lon(ie))*cos(0.5_real32*(poly_lat(ie2)+poly_lat(ie))*deg2rad))**2 )
        nseg = max(1_int32, int(ceiling(edge_deg/denom_deg)))

        do ks=1, nseg
          t0 = real(ks-1,real32)/real(nseg,real32)
          t1 = real(ks,real32)/real(nseg,real32)
          lat_a = poly_lat(ie) + (poly_lat(ie2)-poly_lat(ie))*t0
          lon_a = poly_lon(ie) + (poly_lon(ie2)-poly_lon(ie))*t0
          lat_b = poly_lat(ie) + (poly_lat(ie2)-poly_lat(ie))*t1
          lon_b = poly_lon(ie) + (poly_lon(ie2)-poly_lon(ie))*t1

          do ir=1, nr_s
            r0 = r_bot_m + (real(ir-1,real32)/real(nr_s,real32))*dr_m
            r1 = r_bot_m + (real(ir,real32)/real(nr_s,real32))*dr_m

            call add_side_quad_sources(r0, r1, lat_a, lon_a, lat_b, lon_b, winding, &
                                       mx, my, mz, cap, ns, xs_list, ys_list, zs_list, qmag_list, mass_list)
          end do
        end do
      end do

      write(*,'(A,I0)') '  Source elements (mag surface charges): ', ns

    else
      ! Gravity mode: volume midpoint masses
      do it=1, nlat_s
        lat_deg = lat_lo + (real(it,real32)-0.5_real32)*dlat_s_deg
        if (lat_deg <= -90.0_real32 .or. lat_deg >= 90.0_real32) cycle

        do is=1, nlon_s
          lon_deg = lon_lo + (real(is,real32)-0.5_real32)*dlon_s_deg
          if (.not. point_in_poly_lonlat(lon_deg, lat_deg, poly_lon, poly_lat)) cycle

          lat_rad = lat_deg*deg2rad
          lon_rad = lon_deg*deg2rad

          do ir=1, nr_s
            w_r  = 1.0_real32/real(nr_s,real32)
            r_mid = r_bot_m + ( (real(ir,real32)-0.5_real32)/real(nr_s,real32) ) * dr_m

            dV = (r_mid*r_mid) * cos(lat_rad) * (dlat_s_deg*deg2rad) * (dlon_s_deg*deg2rad) * (dr_m*w_r)
            if (dV <= 0.0_real32) cycle

            xs = r_mid * cos(lat_rad) * cos(lon_rad)
            ys = r_mid * cos(lat_rad) * sin(lon_rad)
            zs = r_mid * sin(lat_rad)

            call append_source(cap, ns, xs_list, ys_list, zs_list, qmag_list, mass_list, &
                               xs, ys, zs, 0.0_real32, rho_kgm3*dV)
          end do
        end do
      end do

      write(*,'(A,I0)') '  Source elements (grav masses): ', ns
    end if
    call cpu_time(t_phase1)
    t_source_s = t_source_s + (t_phase1 - t_phase0)

    if (diag_enabled) then
      obs_cells64 = int(nlat, int64) * int(nlon, int64)
      source_cap64 = int(cap, int64)
      mem_obs_bytes = 6_int64 * obs_cells64 * 4_int64
      mem_source_bytes = source_cap64 * 5_int64 * 4_int64
      mem_total_bytes = mem_obs_bytes + mem_source_bytes
      write(*,'(A,I0)') 'DIAG|point_source|meta|source_capacity=', cap
      write(*,'(A,I0)') 'DIAG|point_source|meta|source_elements=', ns
      write(*,'(A,F12.6)') 'DIAG|point_source|memory|obs_arrays_mib=', bytes_to_mib(mem_obs_bytes)
      write(*,'(A,F12.6)') 'DIAG|point_source|memory|source_arrays_mib=', bytes_to_mib(mem_source_bytes)
      write(*,'(A,F12.6)') 'DIAG|point_source|memory|total_est_mib=', bytes_to_mib(mem_total_bytes)
    end if

    ! ------------------------------------------
    ! Field accumulation: OpenMP over obs grid
    ! ------------------------------------------
    call cpu_time(t_phase0)
    Bx(:,:) = 0.0_real32
    By(:,:) = 0.0_real32
    Bz(:,:) = 0.0_real32

    if (ifield == 2) then
      !$omp parallel do collapse(2) default(none) &
      !$omp shared(nlat,nlon,xo,yo,zo,Bx,By,Bz,ns,xs_list,ys_list,zs_list,qmag_list) &
      !$omp private(i,j,is,dBx,dBy,dBz,bx_acc,by_acc,bz_acc) schedule(static)
      do i=1,nlat
        do j=1,nlon
          bx_acc = 0.0_real32
          by_acc = 0.0_real32
          bz_acc = 0.0_real32

          do is=1, ns
            call add_mcharge_field_cart(xo(i,j), yo(i,j), zo(i,j), &
                                        xs_list(is), ys_list(is), zs_list(is), &
                                        qmag_list(is), dBx, dBy, dBz)
            bx_acc = bx_acc + dBx
            by_acc = by_acc + dBy
            bz_acc = bz_acc + dBz
          end do

          Bx(i,j) = bx_acc
          By(i,j) = by_acc
          Bz(i,j) = bz_acc
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do collapse(2) default(none) &
      !$omp shared(nlat,nlon,xo,yo,zo,Bx,By,Bz,ns,xs_list,ys_list,zs_list,mass_list) &
      !$omp private(i,j,is,dBx,dBy,dBz,bx_acc,by_acc,bz_acc) schedule(static)
      do i=1,nlat
        do j=1,nlon
          bx_acc = 0.0_real32
          by_acc = 0.0_real32
          bz_acc = 0.0_real32

          do is=1, ns
            call add_pointmass_field_cart(xo(i,j), yo(i,j), zo(i,j), &
                                          xs_list(is), ys_list(is), zs_list(is), &
                                          mass_list(is), dBx, dBy, dBz)
            bx_acc = bx_acc + dBx
            by_acc = by_acc + dBy
            bz_acc = bz_acc + dBz
          end do

          Bx(i,j) = bx_acc
          By(i,j) = by_acc
          Bz(i,j) = bz_acc
        end do
      end do
      !$omp end parallel do
    end if
    call cpu_time(t_phase1)
    t_accum_s = t_accum_s + (t_phase1 - t_phase0)

    ! ---- output ----
    call cpu_time(t_phase0)
    if (nfile /= 0) then
      if (ifield == 2) then
        write(2,'(A)') '# body_id lon_deg lat_deg Bx_nT By_nT Bz_nT Btot_nT'
      else
        write(2,'(A)') '# body_id lon_deg lat_deg gx_mGal gy_mGal gz_mGal gtot_mGal'
      end if

      do i=1,nlat
        do j=1,nlon
          lat_deg = lat0_deg + real(i-1,real32)*dlat_deg
          lon_deg = lon0_deg + real(j-1,real32)*dlon_deg
          call wrap180_scalar(lon_deg)

          if (ifield == 2) then
            btot_nt = sqrt( (Bx(i,j)*T_to_nT)**2 + (By(i,j)*T_to_nT)**2 + (Bz(i,j)*T_to_nT)**2 )
            write(2,'(I6,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
                 knk, lon_deg, lat_deg, Bx(i,j)*T_to_nT, By(i,j)*T_to_nT, Bz(i,j)*T_to_nT, btot_nt
          else
            gtot_mgal = sqrt( (Bx(i,j)*mps2_to_mgal)**2 + (By(i,j)*mps2_to_mgal)**2 + (Bz(i,j)*mps2_to_mgal)**2 )
            write(2,'(I6,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
                 knk, lon_deg, lat_deg, Bx(i,j)*mps2_to_mgal, By(i,j)*mps2_to_mgal, Bz(i,j)*mps2_to_mgal, gtot_mgal
          end if
        end do
      end do
    end if
    call cpu_time(t_phase1)
    t_output_s = t_output_s + (t_phase1 - t_phase0)

    if (diag_enabled) then
      call cpu_time(t_body1)
      write(*,'(A,F12.6)') 'DIAG|point_source|time|obs_grid_s=', real(t_obsgrid_s, real64)
      write(*,'(A,F12.6)') 'DIAG|point_source|time|source_mesh_s=', real(t_source_s, real64)
      write(*,'(A,F12.6)') 'DIAG|point_source|time|field_accum_s=', real(t_accum_s, real64)
      write(*,'(A,F12.6)') 'DIAG|point_source|time|output_write_s=', real(t_output_s, real64)
      write(*,'(A,F12.6)') 'DIAG|point_source|time|body_total_s=', real(t_body1 - t_body0, real64)
    end if

  end do

  close(4); close(2)

contains

  pure real(real64) function bytes_to_mib(nbytes)
    integer(int64), intent(in) :: nbytes
    bytes_to_mib = real(nbytes, real64) / (1024.0_real64 * 1024.0_real64)
  end function bytes_to_mib

  !---------------------------------------------------------------------
  ! append_source
  !
  ! centralized append helper for either:
  ! magnetic charge elements (qmag) or gravity mass elements (mass)
  !---------------------------------------------------------------------
  subroutine append_source(cap, ns, xs, ys, zs, qmag, mass, x, y, z, qv, mv)
    integer(int32), intent(inout) :: cap, ns
    real(real32), allocatable, intent(inout) :: xs(:), ys(:), zs(:)
    real(real32), allocatable, intent(inout) :: qmag(:), mass(:)
    real(real32), intent(in) :: x, y, z, qv, mv

    if (ns >= cap) call grow_sources(cap, xs, ys, zs, qmag, mass)

    ns = ns + 1
    xs(ns) = x
    ys(ns) = y
    zs(ns) = z
    qmag(ns) = qv
    mass(ns) = mv
  end subroutine append_source

  !---------------------------------------------------------------------
  ! sph_to_xyz
  !
  ! convert spherical (radius, latitude, longitude) to cartesian xyz
  ! all distances are meters, all angles are degrees on entry
  !---------------------------------------------------------------------
  subroutine sph_to_xyz(r_m, lat_deg, lon_deg, x, y, z)
    real(real32), intent(in) :: r_m, lat_deg, lon_deg
    real(real32), intent(out) :: x, y, z
    real(real32) :: lat_rad_l, lon_rad_l

    lat_rad_l = lat_deg * deg2rad
    lon_rad_l = lon_deg * deg2rad

    x = r_m * cos(lat_rad_l) * cos(lon_rad_l)
    y = r_m * cos(lat_rad_l) * sin(lon_rad_l)
    z = r_m * sin(lat_rad_l)
  end subroutine sph_to_xyz

  !---------------------------------------------------------------------
  ! add_side_quad_sources
  !
  ! triangulates one side-wall quadrilateral strip and delegates each
  ! triangle to add_triangle_charge_source. this is where side-wall
  ! equivalent magnetic charge is built in magnetic mode
  !---------------------------------------------------------------------
  subroutine add_side_quad_sources(r0, r1, lat_a, lon_a, lat_b, lon_b, winding, &
                                   mx, my, mz, cap, ns, xs, ys, zs, qmag, mass)
    real(real32), intent(in) :: r0, r1
    real(real32), intent(in) :: lat_a, lon_a, lat_b, lon_b
    real(real32), intent(in) :: winding, mx, my, mz
    integer(int32), intent(inout) :: cap, ns
    real(real32), allocatable, intent(inout) :: xs(:), ys(:), zs(:)
    real(real32), allocatable, intent(inout) :: qmag(:), mass(:)

    real(real32) :: ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz

    call sph_to_xyz(r0, lat_a, lon_a, ax, ay, az)
    call sph_to_xyz(r0, lat_b, lon_b, bx, by, bz)
    call sph_to_xyz(r1, lat_b, lon_b, cx, cy, cz)
    call sph_to_xyz(r1, lat_a, lon_a, dx, dy, dz)

    ! split the quad into two triangles: (a,b,c) and (a,c,d)
    call add_triangle_charge_source(ax, ay, az, bx, by, bz, cx, cy, cz, &
                                    winding, mx, my, mz, cap, ns, xs, ys, zs, qmag, mass)
    call add_triangle_charge_source(ax, ay, az, cx, cy, cz, dx, dy, dz, &
                                    winding, mx, my, mz, cap, ns, xs, ys, zs, qmag, mass)
  end subroutine add_side_quad_sources

  !---------------------------------------------------------------------
  ! add_triangle_charge_source
  !
  ! computes one equivalent magnetic charge element from a triangle:
  !   q = (M dot n_out) * area
  ! the element is stored at the triangle centroid
  !---------------------------------------------------------------------
  subroutine add_triangle_charge_source(ax, ay, az, bx, by, bz, cx, cy, cz, &
                                        winding, mx, my, mz, cap, ns, xs, ys, zs, qmag, mass)
    real(real32), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz
    real(real32), intent(in) :: winding, mx, my, mz
    integer(int32), intent(inout) :: cap, ns
    real(real32), allocatable, intent(inout) :: xs(:), ys(:), zs(:)
    real(real32), allocatable, intent(inout) :: qmag(:), mass(:)

    real(real32) :: ux, uy, uz, vx, vy, vz, nx, ny, nz
    real(real32) :: nmag, area, gx, gy, gz, qel

    ux = bx - ax
    uy = by - ay
    uz = bz - az
    vx = cx - ax
    vy = cy - ay
    vz = cz - az

    nx = uy*vz - uz*vy
    ny = uz*vx - ux*vz
    nz = ux*vy - uy*vx

    nmag = sqrt(nx*nx + ny*ny + nz*nz)
    if (nmag <= 1.0e-20_real32) return

    area = 0.5_real32 * nmag
    nx = nx / nmag
    ny = ny / nmag
    nz = nz / nmag

    gx = (ax + bx + cx) / 3.0_real32
    gy = (ay + by + cy) / 3.0_real32
    gz = (az + bz + cz) / 3.0_real32

    ! edge x radial cross products point outward for counterclockwise boundaries
    ! reverse their sign for clockwise boundaries when assigning magnetic charge
    qel = winding * (mx*nx + my*ny + mz*nz) * area
    call append_source(cap, ns, xs, ys, zs, qmag, mass, gx, gy, gz, qel, 0.0_real32)
  end subroutine add_triangle_charge_source

  !---------------------------------------------------------------------
  ! grow_sources
  !
  ! dynamic array growth for source element buffers. growth factor is
  ! aggressive to avoid frequent reallocations in large meshes but $$
  !---------------------------------------------------------------------
  subroutine grow_sources(cap, xs, ys, zs, qmag, mass)
    integer(int32), intent(inout) :: cap
    real(real32), allocatable, intent(inout) :: xs(:), ys(:), zs(:)
    real(real32), allocatable, intent(inout) :: qmag(:)
    real(real32), allocatable, intent(inout) :: mass(:)

    integer(int32) :: newcap
    real(real32), allocatable :: xs2(:), ys2(:), zs2(:), qmag2(:), mass2(:)

    ! grow by 2x or by a fixed large block, whichever is larger
    newcap = max(2*cap, cap + 100000_int32)
    allocate(xs2(newcap), ys2(newcap), zs2(newcap), qmag2(newcap), mass2(newcap))

    xs2(1:cap) = xs; ys2(1:cap) = ys; zs2(1:cap) = zs
    qmag2(1:cap) = qmag
    mass2(1:cap) = mass

    call move_alloc(xs2, xs)
    call move_alloc(ys2, ys)
    call move_alloc(zs2, zs)
    call move_alloc(qmag2, qmag)
    call move_alloc(mass2, mass)

    cap = newcap
  end subroutine grow_sources

end program gravmag_sphere_point_source
