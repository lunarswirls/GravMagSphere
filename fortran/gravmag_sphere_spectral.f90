program gravmag_sphere_spectral
  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  use gravmag_paths, only: default_output_path
  use gravmag_sphere_harmonics, only: schmidt_basis
  implicit none

  !***********************************************************************
  ! gravmag_sphere_spectral
  !
  ! spherical-harmonic field solver:
  ! - build elemental dipoles (magnetic) or point masses (gravity)
  ! - fit shared spherical-harmonic coefficients to vector samples on a fitting sphere
  ! - reconstruct field components on the observation grid
  !
  ! resolution note:
  ! - this file implements a regularized least-squares spherical-harmonic fit
  ! - for compact/sharp sources, low lmax can smooth and spread anomalies
  ! - gravmag_sphere_direct_sum and gravmag_sphere_gauss_legendre avoid harmonic truncation
  !***********************************************************************

  integer, parameter :: wp = real64
  real(wp), parameter :: pi = 3.1415926535897932384626433832795_wp
  real(wp), parameter :: deg2rad = pi/180.0_wp
  real(wp), parameter :: mu0_over_4pi = 1.0e-7_wp
  real(wp), parameter :: grav_G = 6.67430e-11_wp
  real(wp), parameter :: T_to_nT = 1.0e9_wp
  real(wp), parameter :: mps2_to_mgal = 1.0e5_wp

  character(len=256) :: infile, outfile, tmp
  integer(int32) :: ios, argc

  real(wp) :: rsphere_km, rsphere_m

  ! solver controls (optional CLI overrides)
  integer(int32) :: lmax, refine_factor, ntheta_fit, nphi_fit
  integer(int32) :: lmax_use
  integer(int32) :: src_nlat_cli, src_nlon_cli, src_nr_cli
  real(wp) :: reg_lambda, reg_power
  real(wp) :: reg_lambda_use, reg_power_use
  integer(int32) :: auto_mode_i, edge_corr_i, hybrid_mode_i, complex_vertex_threshold
  logical :: auto_mode, edge_corr_enabled
  real(wp) :: joint_strength
  real(wp) :: edge_sigma_deg, edge_band_deg, edge_lambda
  real(wp) :: hybrid_band_deg, hybrid_transition_deg
  integer(int32) :: edge_centers_per_edge

  ! per-body input cards
  character(len=256) :: title
  integer(int32) :: knk
  integer(int32) :: nbody_total

  real(wp) :: lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km
  integer(int32) :: nlat, nlon
  real(wp) :: lat0_ref_deg, lon0_ref_deg, dlat_ref_deg, dlon_ref_deg, elvo_ref_km
  integer(int32) :: nlat_ref, nlon_ref
  logical :: have_ref_grid

  integer(int32) :: nr, ntheta, nphi, nblim

  real(wp) :: htheta, hphi, phi1, phi2, rho_dummy

  integer(int32) :: ifield, nrem
  integer(int32) :: ifield_ref
  real(wp) :: M_amp_Apm, inc_deg, dec_deg
  real(wp) :: rho_kgm3

  integer(int32) :: iprint, nfile, nopt
  integer(int32) :: nfile_ref
  real(wp) :: first, conint

  integer(int32) :: npts
  real(wp) :: depth_top_km, depth_bot_km
  real(wp) :: lat_max, lat_min, lon_max, lon_min
  real(wp), allocatable :: poly_lat(:), poly_lon(:)

  ! geometry & arrays
  real(wp) :: ro_m, r_top_m, r_bot_m, dr_m
  real(wp) :: cosI, sinI, cosD, sinD
  real(wp) :: Mx, My, Mz

  real(wp), allocatable :: xo(:,:), yo(:,:), zo(:,:)
  real(wp), allocatable :: Bx(:,:), By(:,:), Bz(:,:)

  real(wp) :: lat_lo, lat_hi, lon_lo, lon_hi, lon_ref
  integer(int32) :: nlat_s, nlon_s, nlat_seed, nlon_seed, nr_s
  real(wp) :: dlat_s_deg, dlon_s_deg

  integer(int32) :: ns, cap
  real(wp), allocatable :: xs_list(:), ys_list(:), zs_list(:)
  real(wp), allocatable :: mx_list(:), my_list(:), mz_list(:)
  real(wp), allocatable :: mass_list(:)
  integer(int32), allocatable :: body_poly_start(:), body_poly_n(:)
  real(wp), allocatable :: body_poly_lon(:), body_poly_lat(:), body_lon_ref(:)
  integer(int32), allocatable :: body_complex_flag(:), body_nblim_arr(:), body_vertex_n(:)
  integer(int32) :: body_meta_cap, body_poly_cap, body_poly_size

  integer(int32) :: i, j, is, it, ir, l, m
  real(wp) :: lat_deg, lon_deg, lat_rad, lon_rad
  real(wp) :: xs, ys, zs, dV, r_mid

  ! Gauss coefficient fit + sampled fitting data
  integer(int32), allocatable :: idx_g(:,:), idx_h(:,:)
  integer(int32) :: ncoef
  real(wp), allocatable :: coef(:)
  integer(int32) :: nsamp, samp_idx, ncorr
  real(wp), allocatable :: samp_theta(:), samp_phi(:), samp_w(:)
  real(wp), allocatable :: samp_lon(:), samp_lat(:)
  real(wp), allocatable :: samp_br(:), samp_bt(:), samp_bp(:)
  real(wp), allocatable :: fit_res_br(:), fit_res_bt(:), fit_res_bp(:)
  real(wp), allocatable :: corr_lon(:), corr_lat(:), corr_coef_br(:), corr_coef_bt(:), corr_coef_bp(:)
  real(wp), allocatable :: corr_center_lon(:), corr_center_lat(:)
  real(wp), allocatable :: corr_coef_br_all(:), corr_coef_bt_all(:), corr_coef_bp_all(:)
  integer(int32), allocatable :: corr_start(:), corr_n(:)
  integer(int32) :: corr_cap, corr_size, ncorr_body
  real(wp) :: fit_rss, fit_regnorm

  real(wp) :: theta, phi, st, ct, cp, sp
  real(wp) :: xb, yb, zb
  real(wp) :: dBx, dBy, dBz
  real(wp) :: br_fit, bt_fit, bp_fit, w
  real(wp) :: br_monopole
  integer(int32) :: itheta, iphi

  real(wp) :: br, bt, bp, bxv, byv, bzv, btot_nt, gtot_mgal
  real(wp) :: reg_w, lnorm
  real(wp) :: corr_br, corr_bt, corr_bp
  real(wp) :: corr_br_body, corr_bt_body, corr_bp_body
  integer(int32) :: out_body_id, ibody, npts_body, p0, c0
  logical :: body_complex, input_is_complex, hybrid_active
  integer(int32) :: hybrid_cells, hybrid_mode_used
  real(wp) :: d_edge, d_edge_min, wd, ws, tau
  real(wp) :: bx_dir, by_dir, bz_dir
  logical :: diag_enabled
  character(len=32) :: diag_env
  real :: t_body0, t_body1, t_phase0, t_phase1
  real :: t_obsgrid_s, t_source_s, t_fit_s, t_solve_s, t_eval_s, t_hybrid_s, t_output_s
  integer(int64) :: obs_cells64, lside64, ncoef64
  integer(int64) :: mem_obs_bytes, mem_source_bytes, mem_solver_bytes, mem_total_bytes

  ! ----------------------------
  ! CLI
  ! ----------------------------
  argc = command_argument_count()
  if (argc < 2) then
    write(*,*) 'Usage: gravmag_sphere_spectral <R_sphere_km> <input.in> [output.txt] '// &
               '[lmax] [refine_factor] [ntheta_fit] [nphi_fit] [reg_lambda] [reg_power] '// &
               '[source_nlat] [source_nlon] [source_nr] [auto_mode] [joint_strength] [edge_correction] '// &
               '[hybrid_mode] [hybrid_band_deg] [complex_vertex_threshold] [hybrid_transition_deg]'
    stop 2
  end if

  call get_command_argument(1, tmp)
  read(tmp,*,iostat=ios) rsphere_km
  if (ios /= 0) stop 'ERROR: could not parse R_sphere_km'

  call get_command_argument(2, infile)
  call get_command_argument(3, outfile)

  lmax = 24
  refine_factor = 2
  ntheta_fit = 72
  nphi_fit = 144
  reg_lambda = 0.2_wp
  reg_power = 4.0_wp
  src_nlat_cli = 0_int32
  src_nlon_cli = 0_int32
  src_nr_cli = 0_int32
  auto_mode_i = 1_int32
  edge_corr_i = 1_int32
  hybrid_mode_i = 1_int32
  complex_vertex_threshold = 12_int32
  auto_mode = .true.
  edge_corr_enabled = .true.
  joint_strength = 1.0_wp
  edge_sigma_deg = 0.6_wp
  edge_band_deg = 1.5_wp
  edge_lambda = 5.0e-2_wp
  edge_centers_per_edge = 2_int32
  hybrid_band_deg = 1.5_wp
  hybrid_transition_deg = 0.75_wp

  if (argc >= 4) then
    call get_command_argument(4, tmp)
    read(tmp,*,iostat=ios) lmax
    if (ios /= 0 .or. lmax < 1) stop 'ERROR: invalid lmax'
  end if
  if (argc >= 5) then
    call get_command_argument(5, tmp)
    read(tmp,*,iostat=ios) refine_factor
    if (ios /= 0 .or. refine_factor < 1) stop 'ERROR: invalid refine_factor'
  end if
  if (argc >= 6) then
    call get_command_argument(6, tmp)
    read(tmp,*,iostat=ios) ntheta_fit
    if (ios /= 0 .or. ntheta_fit < 4) stop 'ERROR: invalid ntheta_fit'
  end if
  if (argc >= 7) then
    call get_command_argument(7, tmp)
    read(tmp,*,iostat=ios) nphi_fit
    if (ios /= 0 .or. nphi_fit < 8) stop 'ERROR: invalid nphi_fit'
  end if
  if (argc >= 8) then
    call get_command_argument(8, tmp)
    read(tmp,*,iostat=ios) reg_lambda
    if (ios /= 0 .or. reg_lambda < 0.0_wp) stop 'ERROR: invalid reg_lambda'
  end if
  if (argc >= 9) then
    call get_command_argument(9, tmp)
    read(tmp,*,iostat=ios) reg_power
    if (ios /= 0 .or. reg_power < 0.0_wp) stop 'ERROR: invalid reg_power'
  end if
  if (argc >= 10) then
    call get_command_argument(10, tmp)
    read(tmp,*,iostat=ios) src_nlat_cli
    if (ios /= 0 .or. src_nlat_cli < 0) stop 'ERROR: invalid source_nlat'
  end if
  if (argc >= 11) then
    call get_command_argument(11, tmp)
    read(tmp,*,iostat=ios) src_nlon_cli
    if (ios /= 0 .or. src_nlon_cli < 0) stop 'ERROR: invalid source_nlon'
  end if
  if (argc >= 12) then
    call get_command_argument(12, tmp)
    read(tmp,*,iostat=ios) src_nr_cli
    if (ios /= 0 .or. src_nr_cli < 0) stop 'ERROR: invalid source_nr'
  end if
  if (argc >= 13) then
    call get_command_argument(13, tmp)
    read(tmp,*,iostat=ios) auto_mode_i
    if (ios /= 0 .or. (auto_mode_i /= 0 .and. auto_mode_i /= 1)) stop 'ERROR: invalid auto_mode (use 0/1)'
    auto_mode = (auto_mode_i == 1)
  end if
  if (argc >= 14) then
    call get_command_argument(14, tmp)
    read(tmp,*,iostat=ios) joint_strength
    if (ios /= 0 .or. joint_strength < 0.0_wp) stop 'ERROR: invalid joint_strength'
  end if
  if (argc >= 15) then
    call get_command_argument(15, tmp)
    read(tmp,*,iostat=ios) edge_corr_i
    if (ios /= 0 .or. (edge_corr_i /= 0 .and. edge_corr_i /= 1)) stop 'ERROR: invalid edge_correction (use 0/1)'
    edge_corr_enabled = (edge_corr_i == 1)
  end if
  if (argc >= 16) then
    call get_command_argument(16, tmp)
    read(tmp,*,iostat=ios) hybrid_mode_i
    if (ios /= 0 .or. (hybrid_mode_i /= 0 .and. hybrid_mode_i /= 1 .and. hybrid_mode_i /= 2)) then
      stop 'ERROR: invalid hybrid_mode (use 0=off,1=auto,2=force)'
    end if
  end if
  if (argc >= 17) then
    call get_command_argument(17, tmp)
    read(tmp,*,iostat=ios) hybrid_band_deg
    if (ios /= 0 .or. hybrid_band_deg < 0.0_wp) stop 'ERROR: invalid hybrid_band_deg'
  end if
  if (argc >= 18) then
    call get_command_argument(18, tmp)
    read(tmp,*,iostat=ios) complex_vertex_threshold
    if (ios /= 0 .or. complex_vertex_threshold < 3_int32) stop 'ERROR: invalid complex_vertex_threshold'
  end if
  if (argc >= 19) then
    call get_command_argument(19, tmp)
    read(tmp,*,iostat=ios) hybrid_transition_deg
    if (ios /= 0 .or. hybrid_transition_deg < 0.0_wp) stop 'ERROR: invalid hybrid_transition_deg'
  end if

  diag_env = ''
  call get_environment_variable('GRAVMAG_DIAGNOSTICS', diag_env)
  diag_enabled = (len_trim(diag_env) > 0 .and. trim(diag_env) /= '0')

  open(4, file=trim(infile), status='old', iostat=ios)
  if (ios /= 0) stop 'ERROR: cannot open input file'

  call default_output_path(trim(infile), outfile, '_spectral.txt')
  open(2, file=trim(outfile), status='replace', iostat=ios)
  if (ios /= 0) stop 'ERROR: cannot open output file'

  rsphere_m = rsphere_km * 1000.0_wp
  knk = 0
  have_ref_grid = .false.
  nbody_total = 0
  t_obsgrid_s = 0.0
  t_source_s = 0.0
  t_fit_s = 0.0
  t_solve_s = 0.0
  t_eval_s = 0.0
  t_hybrid_s = 0.0
  t_output_s = 0.0
  body_meta_cap = 0
  body_poly_cap = 0
  body_poly_size = 0
  corr_cap = 0
  corr_size = 0
  ncorr = 0
  cap = 1024
  ns = 0
  allocate(xs_list(cap), ys_list(cap), zs_list(cap))
  allocate(mx_list(cap), my_list(cap), mz_list(cap))
  allocate(mass_list(cap))
  call cpu_time(t_body0)

  do
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

    read(4,*,iostat=ios) lat0_deg, lon0_deg, dlat_deg, dlon_deg, elvo_km, nlat, nlon
    if (ios /= 0) stop 'ERROR reading Card 2'

    read(4,*,iostat=ios) nr, ntheta, nphi, nblim
    if (ios /= 0) stop 'ERROR reading Card 3'

    read(4,*,iostat=ios) htheta, hphi, phi1, phi2, rho_dummy
    if (ios /= 0) stop 'ERROR reading Card 4'

    read(4,*,iostat=ios) ifield, nrem, M_amp_Apm, inc_deg, dec_deg
    if (ios /= 0) stop 'ERROR reading Card 5'
    if (ifield /= 1 .and. ifield /= 2) stop 'ERROR: Card 5 ifield must be 1(gravity) or 2(magnetic)'
    if (ifield == 2 .and. nrem /= 1) stop 'ERROR: magnetics require nrem=1'

    read(4,*,iostat=ios) iprint, nfile, nopt, first, conint
    if (ios /= 0) stop 'ERROR reading Card 6'

    if (.not. have_ref_grid) then
      lat0_ref_deg = lat0_deg
      lon0_ref_deg = lon0_deg
      dlat_ref_deg = dlat_deg
      dlon_ref_deg = dlon_deg
      elvo_ref_km = elvo_km
      nlat_ref = nlat
      nlon_ref = nlon
      ifield_ref = ifield
      nfile_ref = nfile
      ro_m = (rsphere_km + elvo_km) * 1000.0_wp
      have_ref_grid = .true.
    else
      if (nlat /= nlat_ref .or. nlon /= nlon_ref) stop 'ERROR: collective mode requires identical Card 2 nlat/nlon for all bodies'
      if (abs(lat0_deg - lat0_ref_deg) > 1.0e-9_wp) stop 'ERROR: collective mode requires identical Card 2 lat0 for all bodies'
      if (abs(lon0_deg - lon0_ref_deg) > 1.0e-9_wp) stop 'ERROR: collective mode requires identical Card 2 lon0 for all bodies'
      if (abs(dlat_deg - dlat_ref_deg) > 1.0e-12_wp) stop 'ERROR: collective mode requires identical Card 2 dlat for all bodies'
      if (abs(dlon_deg - dlon_ref_deg) > 1.0e-12_wp) stop 'ERROR: collective mode requires identical Card 2 dlon for all bodies'
      if (abs(elvo_km - elvo_ref_km) > 1.0e-9_wp) stop 'ERROR: collective mode requires identical Card 2 elevation for all bodies'
      if (ifield /= ifield_ref) stop 'ERROR: collective mode requires identical field type for all bodies'
      if (nfile /= nfile_ref) nfile_ref = max(nfile_ref, nfile)
    end if

    if (ifield == 2) then
      cosI = cos(inc_deg*deg2rad)
      sinI = sin(inc_deg*deg2rad)
      cosD = cos(dec_deg*deg2rad)
      sinD = sin(dec_deg*deg2rad)
      Mx = M_amp_Apm * (cosI*cosD)
      My = M_amp_Apm * (cosI*sinD)
      Mz = M_amp_Apm * (sinI)
      rho_kgm3 = 0.0_wp
    else
      Mx = 0.0_wp
      My = 0.0_wp
      Mz = 0.0_wp
      rho_kgm3 = M_amp_Apm
    end if

    if (allocated(poly_lat)) deallocate(poly_lat, poly_lon)
    if (nblim == 1) then
      read(4,*,iostat=ios) lat_max, lat_min, lon_max, lon_min, depth_top_km, depth_bot_km
      if (ios /= 0) stop 'ERROR reading fixed limits (Card 7)'
      npts = 5
      allocate(poly_lat(npts), poly_lon(npts))
      poly_lat = [lat_min, lat_max, lat_max, lat_min, lat_min]
      poly_lon = [lon_min, lon_min, lon_max, lon_max, lon_min]
    else
      read(4,*,iostat=ios) npts, depth_top_km, depth_bot_km
      if (ios /= 0) stop 'ERROR reading polygon header (Card 7)'
      if (npts < 3) stop 'ERROR: polygon needs >=3 vertices'
      allocate(poly_lat(npts), poly_lon(npts))
      do i=1,npts
        read(4,*,iostat=ios) poly_lat(i), poly_lon(i)
        if (ios /= 0) stop 'ERROR reading polygon vertex'
      end do
    end if

    lon_ref = sum(poly_lon) / real(size(poly_lon), wp)
    call wrap180_scalar(lon_ref)
    call unwrap_lons_inplace(poly_lon, lon_ref)

    call ensure_body_meta_capacity(knk, body_meta_cap, body_poly_start, body_poly_n, body_lon_ref, &
                                   body_complex_flag, body_nblim_arr, body_vertex_n)
    body_poly_start(knk) = body_poly_size + 1
    body_poly_n(knk) = npts
    body_lon_ref(knk) = lon_ref
    body_nblim_arr(knk) = nblim
    body_vertex_n(knk) = polygon_vertex_count(poly_lon, poly_lat)
    body_complex = classify_source_complexity(nblim, poly_lon, poly_lat, complex_vertex_threshold)
    body_complex_flag(knk) = merge(1_int32, 0_int32, body_complex)
    if (body_complex) then
      write(*,'(A,I0,A,I0,A)') '  Complexity check: body ', knk, ' classified as COMPLEX (vertices=', &
                                body_vertex_n(knk), ')'
    else
      write(*,'(A,I0,A,I0,A)') '  Complexity check: body ', knk, ' classified as SIMPLE (vertices=', &
                                body_vertex_n(knk), ')'
    end if
    do i=1, npts
      body_poly_size = body_poly_size + 1
      if (body_poly_size > body_poly_cap) then
        call grow_poly_store(body_poly_cap, body_poly_lon, body_poly_lat)
      end if
      body_poly_lon(body_poly_size) = poly_lon(i)
      body_poly_lat(body_poly_size) = poly_lat(i)
    end do

    lat_lo = minval(poly_lat)
    lat_hi = maxval(poly_lat)
    lon_lo = minval(poly_lon)
    lon_hi = maxval(poly_lon)

    if (depth_top_km < 0.0_wp .or. depth_bot_km < 0.0_wp) stop 'ERROR: depth must be >= 0'
    if (depth_bot_km <= depth_top_km) stop 'ERROR: need depth_bot > depth_top'

    r_top_m = rsphere_m - depth_top_km*1000.0_wp
    r_bot_m = rsphere_m - depth_bot_km*1000.0_wp
    if (r_top_m <= 0.0_wp .or. r_bot_m <= 0.0_wp) stop 'ERROR: depth too large'
    if (r_bot_m > r_top_m) then
      dr_m = r_bot_m
      r_bot_m = r_top_m
      r_top_m = dr_m
    end if
    dr_m = r_top_m - r_bot_m

    ! Source mesh is decoupled from Card 2 output grid:
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
    dlat_s_deg = (lat_hi - lat_lo) / real(nlat_s, wp)
    dlon_s_deg = (lon_hi - lon_lo) / real(nlon_s, wp)

    do it=1, nlat_s
      lat_deg = lat_lo + (real(it,wp)-0.5_wp)*dlat_s_deg
      if (lat_deg <= -90.0_wp .or. lat_deg >= 90.0_wp) cycle

      do is=1, nlon_s
        lon_deg = lon_lo + (real(is,wp)-0.5_wp)*dlon_s_deg
        if (.not. point_in_poly_lonlat(lon_deg, lat_deg, poly_lon, poly_lat)) cycle

        lat_rad = lat_deg*deg2rad
        lon_rad = lon_deg*deg2rad
        do ir=1, nr_s
          r_mid = r_bot_m + ((real(ir,wp)-0.5_wp)/real(nr_s,wp))*dr_m
          dV = (r_mid*r_mid) * cos(lat_rad) * (dlat_s_deg*deg2rad) * (dlon_s_deg*deg2rad) * (dr_m/real(nr_s,wp))
          if (dV <= 0.0_wp) cycle

          xs = r_mid * cos(lat_rad) * cos(lon_rad)
          ys = r_mid * cos(lat_rad) * sin(lon_rad)
          zs = r_mid * sin(lat_rad)

          ns = ns + 1
          if (ns > cap) then
            call grow_sources(cap, xs_list, ys_list, zs_list, mx_list, my_list, mz_list, mass_list)
          end if

          xs_list(ns) = xs
          ys_list(ns) = ys
          zs_list(ns) = zs
          if (ifield == 2) then
            mx_list(ns) = Mx * dV
            my_list(ns) = My * dV
            mz_list(ns) = Mz * dV
            mass_list(ns) = 0.0_wp
          else
            mx_list(ns) = 0.0_wp
            my_list(ns) = 0.0_wp
            mz_list(ns) = 0.0_wp
            mass_list(ns) = rho_kgm3 * dV
          end if
        end do
      end do
    end do

    write(*,'(A,I0,A,I0,A,I0)') '  Source mesh cells (nlat,nlon,nr): ', nlat_s, ', ', nlon_s, ', ', nr_s
    write(*,'(A,I0)') '  Cumulative source elements: ', ns
    call cpu_time(t_phase1)
    t_source_s = t_source_s + (t_phase1 - t_phase0)
  end do

  if (.not. have_ref_grid) stop 'ERROR: no bodies found in input file'
  nbody_total = knk
  lat0_deg = lat0_ref_deg
  lon0_deg = lon0_ref_deg
  dlat_deg = dlat_ref_deg
  dlon_deg = dlon_ref_deg
  elvo_km = elvo_ref_km
  nlat = nlat_ref
  nlon = nlon_ref
  ifield = ifield_ref
  nfile = nfile_ref
  input_is_complex = (nbody_total > 1)
  do ibody=1, nbody_total
    if (body_complex_flag(ibody) /= 0) then
      input_is_complex = .true.
      exit
    end if
  end do
  select case (hybrid_mode_i)
  case (0_int32)
    hybrid_active = .false.
    hybrid_mode_used = 0_int32
  case (1_int32)
    hybrid_active = input_is_complex
    hybrid_mode_used = 1_int32
  case default
    hybrid_active = .true.
    hybrid_mode_used = 2_int32
  end select
  hybrid_cells = 0_int32
  if (hybrid_active) then
    write(*,'(A,L1,A,F6.2,A,F6.2,A,I0)') 'Hybrid mode active. input_complex=', input_is_complex, &
                                          ' band_deg=', hybrid_band_deg, ' transition_deg=', hybrid_transition_deg, &
                                          ' vertex_threshold=', complex_vertex_threshold
  else
    write(*,'(A,L1,A,I0)') 'Hybrid mode inactive. input_complex=', input_is_complex, &
                           ' mode=', hybrid_mode_used
  end if
  if (nbody_total > 1) then
    write(*,'(A,I0,A)') 'Collective solve mode active: solving ', nbody_total, ' bodies together.'
  else
    write(*,'(A)') 'Single-body input detected; using the collective solver path.'
  end if

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
  do i=1,nlat
    do j=1,nlon
      lat_deg = lat0_deg + real(i-1,wp)*dlat_deg
      lon_deg = lon0_deg + real(j-1,wp)*dlon_deg
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
  ! Build fitting samples on sphere (Br/Btheta/Bphi)
  ! ------------------------------------------
  call cpu_time(t_phase0)
  nsamp = ntheta_fit * nphi_fit
  if (allocated(samp_theta)) deallocate(samp_theta, samp_phi, samp_w, samp_lon, samp_lat, samp_br, samp_bt, samp_bp)
  allocate(samp_theta(nsamp), samp_phi(nsamp), samp_w(nsamp))
  allocate(samp_lon(nsamp), samp_lat(nsamp))
  allocate(samp_br(nsamp), samp_bt(nsamp), samp_bp(nsamp))

  ! fix the gravity degree-zero coefficient from the collective source mass
  ! fitting only the remainder keeps regularization from changing the total mass
  br_monopole = 0.0_wp
  if (ifield == 1) br_monopole = -grav_g * sum(mass_list(1:ns)) / ro_m**2

  samp_idx = 0
  do itheta=1, ntheta_fit
    theta = (real(itheta,wp)-0.5_wp) * pi / real(ntheta_fit,wp)
    st = sin(theta)
    ct = cos(theta)

    do iphi=1, nphi_fit
      phi = -pi + (real(iphi-1,wp) * 2.0_wp*pi/real(nphi_fit,wp))
      cp = cos(phi)
      sp = sin(phi)
      xb = ro_m * st * cp
      yb = ro_m * st * sp
      zb = ro_m * ct

      dBx = 0.0_wp
      dBy = 0.0_wp
      dBz = 0.0_wp
      do is=1, ns
        if (ifield == 2) then
          call add_dipole_field_cart64(xb, yb, zb, xs_list(is), ys_list(is), zs_list(is), &
                                       mx_list(is), my_list(is), mz_list(is), br, bt, bp)
        else
          call add_pointmass_field_cart64(xb, yb, zb, xs_list(is), ys_list(is), zs_list(is), &
                                          mass_list(is), br, bt, bp)
        end if
        dBx = dBx + br
        dBy = dBy + bt
        dBz = dBz + bp
      end do

      br_fit = dBx*(st*cp) + dBy*(st*sp) + dBz*ct
      bt_fit = dBx*(ct*cp) + dBy*(ct*sp) - dBz*st
      bp_fit = -dBx*sp + dBy*cp

      samp_idx = samp_idx + 1
      samp_theta(samp_idx) = theta
      samp_phi(samp_idx) = phi
      samp_w(samp_idx) = max(st, 1.0e-10_wp)
      samp_lon(samp_idx) = phi / deg2rad
      call wrap180_scalar(samp_lon(samp_idx))
      samp_lat(samp_idx) = 90.0_wp - theta/deg2rad
      samp_br(samp_idx) = br_fit - br_monopole
      samp_bt(samp_idx) = bt_fit
      samp_bp(samp_idx) = bp_fit
    end do
  end do
  call cpu_time(t_phase1)
  t_fit_s = t_fit_s + (t_phase1 - t_phase0)

  ! ------------------------------------------
  ! Auto-select (lmax, reg_lambda, reg_power) from a pilot sweep
  ! ------------------------------------------
  lmax_use = lmax
  reg_lambda_use = reg_lambda
  reg_power_use = reg_power
  if (auto_mode) then
    call cpu_time(t_phase0)
    call auto_select_params(lmax, reg_lambda, reg_power, joint_strength, &
                            samp_theta, samp_phi, samp_w, samp_br, samp_bt, samp_bp, &
                            lmax_use, reg_lambda_use, reg_power_use)
    call cpu_time(t_phase1)
    t_fit_s = t_fit_s + (t_phase1 - t_phase0)
  end if

  ! ------------------------------------------
  ! Joint vector fit (Br/Btheta/Bphi) with shared coefficients
  ! ------------------------------------------
  call cpu_time(t_phase0)
  if (allocated(idx_g)) deallocate(idx_g, idx_h)
  if (allocated(coef)) deallocate(coef)
  call solve_coeffs_from_samples(lmax_use, reg_lambda_use, reg_power_use, joint_strength, &
                                 samp_theta, samp_phi, samp_w, samp_br, samp_bt, samp_bp, &
                                 idx_g, idx_h, coef, ncoef, fit_rss, fit_regnorm, ios)
  if (ios /= 0) stop 'ERROR: failed to solve normal equations for Gauss coefficients'
  call cpu_time(t_phase1)
  t_solve_s = t_solve_s + (t_phase1 - t_phase0)

  if (diag_enabled) then
    obs_cells64 = int(nlat, int64) * int(nlon, int64)
    ncoef64 = int(ncoef, int64)
    lside64 = int(lmax_use + 1, int64)
    mem_obs_bytes = 6_int64 * obs_cells64 * 8_int64
    mem_source_bytes = int(cap, int64) * 7_int64 * 8_int64
    mem_solver_bytes = (ncoef64*ncoef64 + 3_int64*ncoef64 + 3_int64*lside64*lside64) * 8_int64 + &
                       (2_int64*lside64*lside64) * 4_int64
    mem_total_bytes = mem_obs_bytes + mem_source_bytes + mem_solver_bytes
    write(*,'(A)') 'DIAG|spectral|meta|normalization=schmidt_semi_normalized'
    write(*,'(A)') 'DIAG|spectral|meta|condon_shortley_phase=excluded'
    write(*,'(A,I0)') 'DIAG|spectral|meta|source_capacity=', cap
    write(*,'(A,I0)') 'DIAG|spectral|meta|source_elements=', ns
    write(*,'(A,I0)') 'DIAG|spectral|meta|collective_bodies=', nbody_total
    write(*,'(A,I0)') 'DIAG|spectral|meta|lmax_used=', lmax_use
    write(*,'(A,ES12.4)') 'DIAG|spectral|meta|reg_lambda_used=', reg_lambda_use
    write(*,'(A,F8.3)') 'DIAG|spectral|meta|reg_power_used=', reg_power_use
    write(*,'(A,F8.3)') 'DIAG|spectral|meta|joint_strength=', joint_strength
    write(*,'(A,F12.6)') 'DIAG|spectral|memory|obs_arrays_mib=', bytes_to_mib(mem_obs_bytes)
    write(*,'(A,F12.6)') 'DIAG|spectral|memory|source_arrays_mib=', bytes_to_mib(mem_source_bytes)
    write(*,'(A,F12.6)') 'DIAG|spectral|memory|solver_arrays_mib=', bytes_to_mib(mem_solver_bytes)
    write(*,'(A,F12.6)') 'DIAG|spectral|memory|total_est_mib=', bytes_to_mib(mem_total_bytes)
  end if

  ! ------------------------------------------
  ! Local edge RBF correction (sum per-body edge models)
  ! ------------------------------------------
  ncorr = 0
  if (edge_corr_enabled) then
    call cpu_time(t_phase0)
    if (allocated(fit_res_br)) deallocate(fit_res_br, fit_res_bt, fit_res_bp)
    allocate(fit_res_br(nsamp), fit_res_bt(nsamp), fit_res_bp(nsamp))
    ! samples and fitted coefficients both exclude the fixed monopole
    ! their difference is also the residual of the complete field
    do samp_idx=1, nsamp
      call evaluate_from_coeffs_point(samp_theta(samp_idx), samp_phi(samp_idx), lmax_use, idx_g, idx_h, coef, br, bt, bp)
      fit_res_br(samp_idx) = samp_br(samp_idx) - br
      fit_res_bt(samp_idx) = samp_bt(samp_idx) - bt
      fit_res_bp(samp_idx) = samp_bp(samp_idx) - bp
    end do

    if (allocated(corr_start)) deallocate(corr_start, corr_n)
    allocate(corr_start(nbody_total), corr_n(nbody_total))
    corr_start(:) = 0
    corr_n(:) = 0

    corr_cap = max(32_int32, body_poly_size)
    corr_size = 0
    if (allocated(corr_center_lon)) then
      deallocate(corr_center_lon, corr_center_lat, corr_coef_br_all, &
                 corr_coef_bt_all, corr_coef_bp_all)
    end if
    allocate(corr_center_lon(corr_cap), corr_center_lat(corr_cap), &
             corr_coef_br_all(corr_cap), corr_coef_bt_all(corr_cap), &
             corr_coef_bp_all(corr_cap))

    do ibody=1, nbody_total
      p0 = body_poly_start(ibody)
      npts_body = body_poly_n(ibody)
      if (npts_body < 2) cycle
      if (allocated(corr_lon)) deallocate(corr_lon, corr_lat, corr_coef_br, corr_coef_bt, corr_coef_bp)
      call fit_edge_rbf_correction(samp_lon, samp_lat, fit_res_br, fit_res_bt, fit_res_bp, &
                                   body_lon_ref(ibody), body_poly_lon(p0:p0+npts_body-1), body_poly_lat(p0:p0+npts_body-1), &
                                   edge_sigma_deg, edge_band_deg, edge_lambda, edge_centers_per_edge, &
                                   corr_lon, corr_lat, corr_coef_br, corr_coef_bt, corr_coef_bp, ncorr_body, ios)
      if (ios /= 0 .or. ncorr_body <= 0) cycle
      if (corr_size + ncorr_body > corr_cap) then
        do while (corr_size + ncorr_body > corr_cap)
          call grow_correction_store(corr_cap, corr_center_lon, corr_center_lat, &
                                     corr_coef_br_all, corr_coef_bt_all, corr_coef_bp_all)
        end do
      end if
      c0 = corr_size + 1
      corr_start(ibody) = c0
      corr_n(ibody) = ncorr_body
      corr_center_lon(c0:c0+ncorr_body-1) = corr_lon(1:ncorr_body)
      corr_center_lat(c0:c0+ncorr_body-1) = corr_lat(1:ncorr_body)
      corr_coef_br_all(c0:c0+ncorr_body-1) = corr_coef_br(1:ncorr_body)
      corr_coef_bt_all(c0:c0+ncorr_body-1) = corr_coef_bt(1:ncorr_body)
      corr_coef_bp_all(c0:c0+ncorr_body-1) = corr_coef_bp(1:ncorr_body)
      corr_size = corr_size + ncorr_body
      ncorr = ncorr + ncorr_body
    end do
    call cpu_time(t_phase1)
    t_fit_s = t_fit_s + (t_phase1 - t_phase0)
  end if

  ! ------------------------------------------
  ! Evaluate field from coefficients on obs grid
  ! ------------------------------------------
  call cpu_time(t_phase0)
  do i=1, nlat
    do j=1, nlon
      lat_deg = lat0_deg + real(i-1,wp)*dlat_deg
      lon_deg = lon0_deg + real(j-1,wp)*dlon_deg
      call wrap180_scalar(lon_deg)

      theta = (90.0_wp - lat_deg) * deg2rad
      phi = lon_deg * deg2rad
      st = sin(theta)
      ct = cos(theta)
      cp = cos(phi)
      sp = sin(phi)

      call evaluate_from_coeffs_point(theta, phi, lmax_use, idx_g, idx_h, coef, br, bt, bp)
      br = br + br_monopole

      if (edge_corr_enabled .and. ncorr > 0) then
        corr_br = 0.0_wp
        corr_bt = 0.0_wp
        corr_bp = 0.0_wp
        do ibody=1, nbody_total
          if (corr_n(ibody) <= 0) cycle
          p0 = body_poly_start(ibody)
          npts_body = body_poly_n(ibody)
          c0 = corr_start(ibody)
          call evaluate_edge_rbf_correction(lon_deg, lat_deg, body_lon_ref(ibody), &
                                            body_poly_lon(p0:p0+npts_body-1), body_poly_lat(p0:p0+npts_body-1), &
                                            corr_center_lon(c0:c0+corr_n(ibody)-1), corr_center_lat(c0:c0+corr_n(ibody)-1), &
                                            corr_coef_br_all(c0:c0+corr_n(ibody)-1), corr_coef_bt_all(c0:c0+corr_n(ibody)-1), &
                                            corr_coef_bp_all(c0:c0+corr_n(ibody)-1), &
                                            edge_sigma_deg, edge_band_deg, corr_br_body, corr_bt_body, corr_bp_body)
          corr_br = corr_br + corr_br_body
          corr_bt = corr_bt + corr_bt_body
          corr_bp = corr_bp + corr_bp_body
        end do
        br = br + corr_br
        bt = bt + corr_bt
        bp = bp + corr_bp
      end if

      bxv = br*st*cp + bt*ct*cp - bp*sp
      byv = br*st*sp + bt*ct*sp + bp*cp
      bzv = br*ct - bt*st

      if (hybrid_active) then
        d_edge_min = huge(1.0_wp)
        do ibody=1, nbody_total
          p0 = body_poly_start(ibody)
          npts_body = body_poly_n(ibody)
          if (npts_body < 2) cycle
          d_edge = min_distance_to_poly_edges(unwrap_lon_to_ref(lon_deg, body_lon_ref(ibody)), lat_deg, &
                                              body_poly_lon(p0:p0+npts_body-1), body_poly_lat(p0:p0+npts_body-1))
          if (d_edge < d_edge_min) d_edge_min = d_edge
        end do
        if (d_edge_min <= hybrid_band_deg + hybrid_transition_deg) then
          call cpu_time(t_phase0)
          bx_dir = 0.0_wp
          by_dir = 0.0_wp
          bz_dir = 0.0_wp
          xb = xo(i,j)
          yb = yo(i,j)
          zb = zo(i,j)
          do is=1, ns
            if (ifield == 2) then
              call add_dipole_field_cart64(xb, yb, zb, xs_list(is), ys_list(is), zs_list(is), &
                                           mx_list(is), my_list(is), mz_list(is), dBx, dBy, dBz)
            else
              call add_pointmass_field_cart64(xb, yb, zb, xs_list(is), ys_list(is), zs_list(is), &
                                              mass_list(is), dBx, dBy, dBz)
            end if
            bx_dir = bx_dir + dBx
            by_dir = by_dir + dBy
            bz_dir = bz_dir + dBz
          end do
          call cpu_time(t_phase1)
          t_hybrid_s = t_hybrid_s + (t_phase1 - t_phase0)

          if (hybrid_transition_deg <= 0.0_wp) then
            if (d_edge_min <= hybrid_band_deg) then
              wd = 1.0_wp
            else
              wd = 0.0_wp
            end if
          else
            if (d_edge_min <= hybrid_band_deg) then
              wd = 1.0_wp
            else if (d_edge_min >= hybrid_band_deg + hybrid_transition_deg) then
              wd = 0.0_wp
            else
              tau = (d_edge_min - hybrid_band_deg) / hybrid_transition_deg
              ws = tau*tau*(3.0_wp - 2.0_wp*tau)
              wd = 1.0_wp - ws
            end if
          end if
          ws = 1.0_wp - wd
          bxv = wd*bx_dir + ws*bxv
          byv = wd*by_dir + ws*byv
          bzv = wd*bz_dir + ws*bzv
          hybrid_cells = hybrid_cells + 1
        end if
      end if
      Bx(i,j) = bxv
      By(i,j) = byv
      Bz(i,j) = bzv
    end do
  end do
  call cpu_time(t_phase1)
  t_eval_s = t_eval_s + (t_phase1 - t_phase0)

  call cpu_time(t_phase0)
  if (nfile /= 0) then
    out_body_id = 1
    if (hybrid_active) then
      write(2,'(A)') '# solver=spherical_harmonic_hybrid'
    else
      write(2,'(A)') '# solver=spherical_harmonic_spectral'
    end if
    if (nbody_total > 1) then
      write(2,'(A,I0)') '# collective_multi_body=1 n_bodies=', nbody_total
    end if
    write(2,'(A)') '# normalization=schmidt_semi_normalized condon_shortley_phase=excluded'
    write(2,'(A)') '# numerical_floor=epsilon_times_mean_data_diagonal'
    write(2,'(A,I0,A,ES13.5,A,F8.3)') '# lmax=', lmax_use, ' reg_lambda=', reg_lambda_use, ' reg_power=', reg_power_use
    write(2,'(A,I0,A,F8.3,A,I0,A,I0)') '# auto_mode=', merge(1,0,auto_mode), &
                                        ' joint_strength=', joint_strength, &
                                        ' edge_correction=', merge(1,0,edge_corr_enabled), &
                                        ' edge_centers=', ncorr
    write(2,'(A,I0,A,L1,A,I0,A,F7.3,A,F7.3,A,I0)') '# hybrid_mode=', hybrid_mode_used, &
                                        ' input_complex=', input_is_complex, &
                                        ' hybrid_active=', merge(1,0,hybrid_active), &
                                        ' hybrid_band_deg=', hybrid_band_deg, &
                                        ' hybrid_transition_deg=', hybrid_transition_deg, &
                                        ' hybrid_cells=', hybrid_cells
    if (ifield == 2) then
      write(2,'(A)') '# body_id lon_deg[-180,180) lat_deg Bx_nT By_nT Bz_nT Btot_nT'
    else
      write(2,'(A)') '# body_id lon_deg[-180,180) lat_deg gx_mGal gy_mGal gz_mGal gtot_mGal'
    end if
    do i=1,nlat
      do j=1,nlon
        lat_deg = lat0_deg + real(i-1,wp)*dlat_deg
        lon_deg = lon0_deg + real(j-1,wp)*dlon_deg
        call wrap180_scalar(lon_deg)
        if (ifield == 2) then
          btot_nt = sqrt((Bx(i,j)*T_to_nT)**2 + (By(i,j)*T_to_nT)**2 + (Bz(i,j)*T_to_nT)**2)
          write(2,'(I6,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
            out_body_id, lon_deg, lat_deg, Bx(i,j)*T_to_nT, By(i,j)*T_to_nT, Bz(i,j)*T_to_nT, btot_nt
        else
          gtot_mgal = sqrt((Bx(i,j)*mps2_to_mgal)**2 + (By(i,j)*mps2_to_mgal)**2 + (Bz(i,j)*mps2_to_mgal)**2)
          write(2,'(I6,1X,F12.6,1X,F12.6,1X,ES15.7,1X,ES15.7,1X,ES15.7,1X,ES15.7)') &
            out_body_id, lon_deg, lat_deg, Bx(i,j)*mps2_to_mgal, By(i,j)*mps2_to_mgal, Bz(i,j)*mps2_to_mgal, gtot_mgal
        end if
      end do
    end do
  end if
  call cpu_time(t_phase1)
  t_output_s = t_output_s + (t_phase1 - t_phase0)

  if (diag_enabled) then
    call cpu_time(t_body1)
    write(*,'(A,F12.6)') 'DIAG|spectral|time|obs_grid_s=', real(t_obsgrid_s, wp)
    write(*,'(A,F12.6)') 'DIAG|spectral|time|source_mesh_s=', real(t_source_s, wp)
    write(*,'(A,F12.6)') 'DIAG|spectral|time|fit_assembly_s=', real(t_fit_s, wp)
    write(*,'(A,F12.6)') 'DIAG|spectral|time|linear_solve_s=', real(t_solve_s, wp)
    write(*,'(A,F12.6)') 'DIAG|spectral|time|field_eval_s=', real(t_eval_s, wp)
    write(*,'(A,F12.6)') 'DIAG|spectral|time|hybrid_direct_eval_s=', real(t_hybrid_s, wp)
    write(*,'(A,F12.6)') 'DIAG|spectral|time|output_write_s=', real(t_output_s, wp)
    write(*,'(A,I0)') 'DIAG|spectral|meta|hybrid_cells=', hybrid_cells
    write(*,'(A,F12.6)') 'DIAG|spectral|time|body_total_s=', real(t_body1 - t_body0, wp)
  end if

  close(4)
  close(2)

contains

  pure real(wp) function bytes_to_mib(nbytes)
    integer(int64), intent(in) :: nbytes
    bytes_to_mib = real(nbytes, wp) / (1024.0_wp*1024.0_wp)
  end function bytes_to_mib

  ! keep longitudes on the canonical branch [-180,180)
  subroutine wrap180_scalar(lon_deg)
    real(wp), intent(inout) :: lon_deg
    lon_deg = modulo(lon_deg + 180.0_wp, 360.0_wp) - 180.0_wp
  end subroutine wrap180_scalar

  ! shift longitudes to a continuous branch near lon_ref_deg so polygon
  ! edges do not jump across the dateline during point-in-polygon tests
  subroutine unwrap_lons_inplace(lon_deg, lon_ref_deg)
    real(wp), intent(inout) :: lon_deg(:)
    real(wp), intent(in) :: lon_ref_deg
    integer(int32) :: i
    real(wp) :: x
    do i=1, size(lon_deg)
      x = lon_deg(i)
      call wrap180_scalar(x)
      do while (x - lon_ref_deg >= 180.0_wp)
        x = x - 360.0_wp
      end do
      do while (x - lon_ref_deg < -180.0_wp)
        x = x + 360.0_wp
      end do
      lon_deg(i) = x
    end do
  end subroutine unwrap_lons_inplace

  ! standard ray-casting polygon inclusion in lon-lat space
  ! polygon longitudes are assumed already unwrapped to a continuous branch
  logical function point_in_poly_lonlat(lon, lat, poly_lon, poly_lat)
    real(wp), intent(in) :: lon, lat
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)

    integer(int32) :: i, j, n
    real(wp) :: xi, yi, xj, yj, den, xint
    logical :: inside, intersect

    n = size(poly_lon)
    if (n < 3) then
      point_in_poly_lonlat = .false.
      return
    end if

    inside = .false.
    j = n
    do i=1, n
      xi = poly_lon(i); yi = poly_lat(i)
      xj = poly_lon(j); yj = poly_lat(j)

      intersect = ((yi > lat) .neqv. (yj > lat))
      if (intersect) then
        den = yj - yi
        if (abs(den) < 1.0e-30_wp) den = sign(1.0e-30_wp, den + 1.0e-30_wp)
        xint = (xj - xi) * (lat - yi) / den + xi
        if (lon < xint) inside = .not. inside
      end if

      j = i
    end do

    point_in_poly_lonlat = inside
  end function point_in_poly_lonlat

  ! dynamic growth helper for source-element arrays
  subroutine grow_sources(cap, xs, ys, zs, mx, my, mz, mass)
    integer(int32), intent(inout) :: cap
    real(wp), allocatable, intent(inout) :: xs(:), ys(:), zs(:)
    real(wp), allocatable, intent(inout) :: mx(:), my(:), mz(:)
    real(wp), allocatable, intent(inout) :: mass(:)

    integer(int32) :: newcap
    real(wp), allocatable :: xs2(:), ys2(:), zs2(:), mx2(:), my2(:), mz2(:), mass2(:)

    ! grow by 2x or by a fixed block to reduce reallocation churn
    newcap = max(2*cap, cap + 100000_int32)
    allocate(xs2(newcap), ys2(newcap), zs2(newcap), mx2(newcap), my2(newcap), mz2(newcap), mass2(newcap))

    xs2(1:cap) = xs; ys2(1:cap) = ys; zs2(1:cap) = zs
    mx2(1:cap) = mx; my2(1:cap) = my; mz2(1:cap) = mz
    mass2(1:cap) = mass

    call move_alloc(xs2, xs)
    call move_alloc(ys2, ys)
    call move_alloc(zs2, zs)
    call move_alloc(mx2, mx)
    call move_alloc(my2, my)
    call move_alloc(mz2, mz)
    call move_alloc(mass2, mass)

    cap = newcap
  end subroutine grow_sources

  subroutine ensure_body_meta_capacity(need, cap, poly_start, poly_n, lon_ref_arr, complex_arr, nblim_arr, vertex_arr)
    integer(int32), intent(in) :: need
    integer(int32), intent(inout) :: cap
    integer(int32), allocatable, intent(inout) :: poly_start(:), poly_n(:)
    real(wp), allocatable, intent(inout) :: lon_ref_arr(:)
    integer(int32), allocatable, intent(inout) :: complex_arr(:), nblim_arr(:), vertex_arr(:)

    integer(int32) :: newcap
    integer(int32), allocatable :: s2(:), n2(:), c2(:), b2(:), v2(:)
    real(wp), allocatable :: r2(:)

    if (need <= cap) return
    newcap = max(need, max(8_int32, 2_int32*max(1_int32, cap)))
    allocate(s2(newcap), n2(newcap), r2(newcap), c2(newcap), b2(newcap), v2(newcap))
    s2(:) = 0_int32
    n2(:) = 0_int32
    r2(:) = 0.0_wp
    c2(:) = 0_int32
    b2(:) = 0_int32
    v2(:) = 0_int32
    if (cap > 0) then
      s2(1:cap) = poly_start(1:cap)
      n2(1:cap) = poly_n(1:cap)
      r2(1:cap) = lon_ref_arr(1:cap)
      c2(1:cap) = complex_arr(1:cap)
      b2(1:cap) = nblim_arr(1:cap)
      v2(1:cap) = vertex_arr(1:cap)
    end if
    call move_alloc(s2, poly_start)
    call move_alloc(n2, poly_n)
    call move_alloc(r2, lon_ref_arr)
    call move_alloc(c2, complex_arr)
    call move_alloc(b2, nblim_arr)
    call move_alloc(v2, vertex_arr)
    cap = newcap
  end subroutine ensure_body_meta_capacity

  subroutine grow_poly_store(cap, poly_lon_arr, poly_lat_arr)
    integer(int32), intent(inout) :: cap
    real(wp), allocatable, intent(inout) :: poly_lon_arr(:), poly_lat_arr(:)

    integer(int32) :: newcap
    real(wp), allocatable :: lon2(:), lat2(:)

    if (cap <= 0) then
      newcap = 256_int32
      allocate(poly_lon_arr(newcap), poly_lat_arr(newcap))
      poly_lon_arr(:) = 0.0_wp
      poly_lat_arr(:) = 0.0_wp
      cap = newcap
      return
    end if

    newcap = max(2_int32*cap, cap + 1024_int32)
    allocate(lon2(newcap), lat2(newcap))
    lon2(:) = 0.0_wp
    lat2(:) = 0.0_wp
    lon2(1:cap) = poly_lon_arr(1:cap)
    lat2(1:cap) = poly_lat_arr(1:cap)
    call move_alloc(lon2, poly_lon_arr)
    call move_alloc(lat2, poly_lat_arr)
    cap = newcap
  end subroutine grow_poly_store

  subroutine grow_correction_store(cap, cen_lon, cen_lat, cbr, cbt, cbp)
    integer(int32), intent(inout) :: cap
    real(wp), allocatable, intent(inout) :: cen_lon(:), cen_lat(:), cbr(:), cbt(:), cbp(:)

    integer(int32) :: newcap
    real(wp), allocatable :: lon2(:), lat2(:), br2(:), bt2(:), bp2(:)

    if (cap <= 0) then
      newcap = 128_int32
      allocate(cen_lon(newcap), cen_lat(newcap), cbr(newcap), cbt(newcap), cbp(newcap))
      cen_lon(:) = 0.0_wp
      cen_lat(:) = 0.0_wp
      cbr(:) = 0.0_wp
      cbt(:) = 0.0_wp
      cbp(:) = 0.0_wp
      cap = newcap
      return
    end if

    newcap = max(2_int32*cap, cap + 512_int32)
    allocate(lon2(newcap), lat2(newcap), br2(newcap), bt2(newcap), bp2(newcap))
    lon2(:) = 0.0_wp
    lat2(:) = 0.0_wp
    br2(:) = 0.0_wp
    bt2(:) = 0.0_wp
    bp2(:) = 0.0_wp
    lon2(1:cap) = cen_lon(1:cap)
    lat2(1:cap) = cen_lat(1:cap)
    br2(1:cap) = cbr(1:cap)
    bt2(1:cap) = cbt(1:cap)
    bp2(1:cap) = cbp(1:cap)
    call move_alloc(lon2, cen_lon)
    call move_alloc(lat2, cen_lat)
    call move_alloc(br2, cbr)
    call move_alloc(bt2, cbt)
    call move_alloc(bp2, cbp)
    cap = newcap
  end subroutine grow_correction_store

  ! dipole field kernel in double precision
  subroutine add_dipole_field_cart64(xo, yo, zo, xs, ys, zs, mx, my, mz, dBx, dBy, dBz)
    real(wp), intent(in) :: xo, yo, zo
    real(wp), intent(in) :: xs, ys, zs
    real(wp), intent(in) :: mx, my, mz
    real(wp), intent(out) :: dBx, dBy, dBz

    real(wp) :: rx, ry, rz, r2, invr, invr3, invr5, mdotr, c

    rx = xo - xs
    ry = yo - ys
    rz = zo - zs

    r2 = rx*rx + ry*ry + rz*rz
    if (r2 <= 0.0_wp) then
      dBx = 0.0_wp; dBy = 0.0_wp; dBz = 0.0_wp
      return
    end if

    invr = 1.0_wp / sqrt(r2)
    invr3 = invr*invr*invr
    invr5 = invr3*invr*invr

    mdotr = mx*rx + my*ry + mz*rz
    c = 3.0_wp * mdotr * invr5

    dBx = mu0_over_4pi * (c*rx - mx*invr3)
    dBy = mu0_over_4pi * (c*ry - my*invr3)
    dBz = mu0_over_4pi * (c*rz - mz*invr3)
  end subroutine add_dipole_field_cart64

  ! point-mass gravity kernel in double precision
  subroutine add_pointmass_field_cart64(xo, yo, zo, xs, ys, zs, mass_kg, dgx, dgy, dgz)
    real(wp), intent(in) :: xo, yo, zo
    real(wp), intent(in) :: xs, ys, zs
    real(wp), intent(in) :: mass_kg
    real(wp), intent(out) :: dgx, dgy, dgz

    real(wp) :: rx, ry, rz, r2, invr, invr3, c

    rx = xs - xo
    ry = ys - yo
    rz = zs - zo

    r2 = rx*rx + ry*ry + rz*rz
    if (r2 <= 0.0_wp) then
      dgx = 0.0_wp; dgy = 0.0_wp; dgz = 0.0_wp
      return
    end if

    invr = 1.0_wp / sqrt(r2)
    invr3 = invr*invr*invr
    c = grav_G * mass_kg * invr3

    dgx = c * rx
    dgy = c * ry
    dgz = c * rz
  end subroutine add_pointmass_field_cart64

  ! map (l,m) coefficient pairs into packed 1-D vector indices
  ! g_lm and h_lm each get their own index stream except m=0
  ! where h_l0 is not used
  subroutine build_index_maps(lm, ig, ih, ncoef_out)
    integer(int32), intent(in) :: lm
    integer(int32), intent(out) :: ig(0:lm,0:lm), ih(0:lm,0:lm)
    integer(int32), intent(out) :: ncoef_out

    integer(int32) :: l, m, idx

    ig = 0
    ih = 0
    idx = 0
    do l=1, lm
      do m=0, l
        idx = idx + 1
        ig(l,m) = idx
        if (m > 0) then
          idx = idx + 1
          ih(l,m) = idx
        end if
      end do
    end do
    ncoef_out = idx
  end subroutine build_index_maps

  pure real(wp) function unwrap_lon_to_ref(lon_deg, lon_ref_deg)
    real(wp), intent(in) :: lon_deg, lon_ref_deg
    unwrap_lon_to_ref = modulo(lon_deg + 180.0_wp, 360.0_wp) - 180.0_wp
    do while (unwrap_lon_to_ref - lon_ref_deg >= 180.0_wp)
      unwrap_lon_to_ref = unwrap_lon_to_ref - 360.0_wp
    end do
    do while (unwrap_lon_to_ref - lon_ref_deg < -180.0_wp)
      unwrap_lon_to_ref = unwrap_lon_to_ref + 360.0_wp
    end do
  end function unwrap_lon_to_ref

  pure real(wp) function local_dist2_deg(lon1, lat1, lon2, lat2)
    real(wp), intent(in) :: lon1, lat1, lon2, lat2
    real(wp) :: dlon, dlat, c
    c = max(cos(0.5_wp*(lat1 + lat2)*deg2rad), 1.0e-6_wp)
    dlon = (lon1 - lon2) * c
    dlat = lat1 - lat2
    local_dist2_deg = dlon*dlon + dlat*dlat
  end function local_dist2_deg

  pure real(wp) function point_segment_distance_deg(lon, lat, lon1, lat1, lon2, lat2)
    real(wp), intent(in) :: lon, lat, lon1, lat1, lon2, lat2
    real(wp) :: c, px, py, ax, ay, bx, by, abx, aby, apx, apy, t, qx, qy
    c = max(cos(lat*deg2rad), 1.0e-6_wp)
    px = lon * c
    py = lat
    ax = lon1 * c
    ay = lat1
    bx = lon2 * c
    by = lat2
    abx = bx - ax
    aby = by - ay
    apx = px - ax
    apy = py - ay
    if (abx*abx + aby*aby <= 1.0e-18_wp) then
      qx = ax
      qy = ay
    else
      t = (apx*abx + apy*aby) / (abx*abx + aby*aby)
      t = max(0.0_wp, min(1.0_wp, t))
      qx = ax + t*abx
      qy = ay + t*aby
    end if
    point_segment_distance_deg = sqrt((px - qx)*(px - qx) + (py - qy)*(py - qy))
  end function point_segment_distance_deg

  pure real(wp) function min_distance_to_poly_edges(lon, lat, poly_lon, poly_lat)
    real(wp), intent(in) :: lon, lat
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)
    integer(int32) :: n, nseg, i, i1, i2
    logical :: closed
    real(wp) :: d
    n = size(poly_lon)
    if (n < 2) then
      min_distance_to_poly_edges = 1.0e30_wp
      return
    end if
    closed = (abs(poly_lon(1)-poly_lon(n)) + abs(poly_lat(1)-poly_lat(n)) < 1.0e-9_wp)
    if (closed) then
      nseg = n - 1
    else
      nseg = n
    end if
    min_distance_to_poly_edges = 1.0e30_wp
    do i=1, nseg
      i1 = i
      if (i < n) then
        i2 = i + 1
      else
        i2 = 1
      end if
      d = point_segment_distance_deg(lon, lat, poly_lon(i1), poly_lat(i1), poly_lon(i2), poly_lat(i2))
      if (d < min_distance_to_poly_edges) min_distance_to_poly_edges = d
    end do
  end function min_distance_to_poly_edges

  pure integer(int32) function polygon_vertex_count(poly_lon, poly_lat)
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)
    integer(int32) :: n
    n = size(poly_lon)
    if (n <= 0) then
      polygon_vertex_count = 0_int32
      return
    end if
    if (n >= 2 .and. (abs(poly_lon(1)-poly_lon(n)) + abs(poly_lat(1)-poly_lat(n)) < 1.0e-9_wp)) then
      polygon_vertex_count = max(0_int32, n - 1_int32)
    else
      polygon_vertex_count = n
    end if
  end function polygon_vertex_count

  pure logical function polygon_is_convex(poly_lon, poly_lat)
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)
    integer(int32) :: n, nv, i, i0, i1, i2
    real(wp) :: zcross, sign_ref
    logical :: have_sign

    nv = polygon_vertex_count(poly_lon, poly_lat)
    if (nv < 4) then
      polygon_is_convex = .true.
      return
    end if
    n = size(poly_lon)
    have_sign = .false.
    sign_ref = 0.0_wp
    do i=1, nv
      i0 = i
      i1 = modulo(i, nv) + 1
      i2 = modulo(i+1, nv) + 1
      zcross = (poly_lon(i1)-poly_lon(i0))*(poly_lat(i2)-poly_lat(i1)) - &
               (poly_lat(i1)-poly_lat(i0))*(poly_lon(i2)-poly_lon(i1))
      if (abs(zcross) <= 1.0e-12_wp) cycle
      if (.not. have_sign) then
        sign_ref = sign(1.0_wp, zcross)
        have_sign = .true.
      else
        if (sign_ref * zcross < 0.0_wp) then
          polygon_is_convex = .false.
          return
        end if
      end if
    end do
    polygon_is_convex = .true.
  end function polygon_is_convex

  pure real(wp) function polygon_compactness(poly_lon, poly_lat)
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)
    integer(int32) :: nv, i, i1, i2
    real(wp) :: lat_mean, c, x1, y1, x2, y2, area2, perim

    nv = polygon_vertex_count(poly_lon, poly_lat)
    if (nv < 3) then
      polygon_compactness = huge(1.0_wp)
      return
    end if
    lat_mean = 0.0_wp
    do i=1, nv
      lat_mean = lat_mean + poly_lat(i)
    end do
    lat_mean = lat_mean / real(nv, wp)
    c = max(cos(lat_mean*deg2rad), 1.0e-6_wp)

    area2 = 0.0_wp
    perim = 0.0_wp
    do i=1, nv
      i1 = i
      i2 = modulo(i, nv) + 1
      x1 = poly_lon(i1) * c
      y1 = poly_lat(i1)
      x2 = poly_lon(i2) * c
      y2 = poly_lat(i2)
      area2 = area2 + (x1*y2 - x2*y1)
      perim = perim + sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))
    end do
    if (abs(area2) <= 1.0e-12_wp) then
      polygon_compactness = huge(1.0_wp)
    else
      polygon_compactness = (perim*perim) / (2.0_wp*pi*abs(area2))
    end if
  end function polygon_compactness

  pure logical function classify_source_complexity(nblim_body, poly_lon, poly_lat, vertex_threshold)
    integer(int32), intent(in) :: nblim_body, vertex_threshold
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)
    integer(int32) :: nv
    logical :: is_convex
    real(wp) :: compactness

    if (nblim_body == 1_int32) then
      classify_source_complexity = .false.
      return
    end if

    nv = polygon_vertex_count(poly_lon, poly_lat)
    is_convex = polygon_is_convex(poly_lon, poly_lat)
    compactness = polygon_compactness(poly_lon, poly_lat)

    classify_source_complexity = .false.
    if (nv > vertex_threshold) classify_source_complexity = .true.
    if (.not. is_convex .and. nv >= max(5_int32, vertex_threshold - 2_int32)) classify_source_complexity = .true.
    if (compactness > 2.25_wp .and. nv >= max(6_int32, vertex_threshold - 4_int32)) classify_source_complexity = .true.
  end function classify_source_complexity

  subroutine evaluate_from_coeffs_point(theta, phi, lm, ig, ih, coef_in, br_out, bt_out, bp_out)
    real(wp), intent(in) :: theta, phi
    integer(int32), intent(in) :: lm
    integer(int32), intent(in) :: ig(0:lm,0:lm), ih(0:lm,0:lm)
    real(wp), intent(in) :: coef_in(:)
    real(wp), intent(out) :: br_out, bt_out, bp_out
    integer(int32) :: ncoef_local
    real(wp), allocatable :: row_br(:), row_bt(:), row_bp(:)

    ncoef_local = size(coef_in)
    allocate(row_br(ncoef_local), row_bt(ncoef_local), row_bp(ncoef_local))
    call schmidt_basis(theta, phi, lm, ig, ih, row_br, row_bt, row_bp)
    br_out = sum(row_br * coef_in)
    bt_out = sum(row_bt * coef_in)
    bp_out = sum(row_bp * coef_in)
    deallocate(row_br, row_bt, row_bp)
  end subroutine evaluate_from_coeffs_point

  subroutine solve_coeffs_from_samples(lm, reg_lam, reg_pow, joint_w, &
                                       theta_s, phi_s, w_s, br_s, bt_s, bp_s, &
                                       ig, ih, coef_out, ncoef_out, rss_out, regnorm_out, info)
    integer(int32), intent(in) :: lm
    real(wp), intent(in) :: reg_lam, reg_pow, joint_w
    real(wp), intent(in) :: theta_s(:), phi_s(:), w_s(:), br_s(:), bt_s(:), bp_s(:)
    integer(int32), allocatable, intent(out) :: ig(:,:), ih(:,:)
    real(wp), allocatable, intent(out) :: coef_out(:)
    integer(int32), intent(out) :: ncoef_out, info
    real(wp), intent(out) :: rss_out, regnorm_out

    integer(int32) :: nobs, i, j, s, l, m
    real(wp), allocatable :: ATA(:,:), ATb(:), row_br(:), row_bt(:), row_bp(:)
    real(wp) :: w, pred_br, pred_bt, pred_bp, rr, lnorm, reg_w, numerical_floor
    real(wp), allocatable :: normal_diagonal(:)

    nobs = size(theta_s)
    info = 0
    rss_out = 0.0_wp
    regnorm_out = 0.0_wp

    allocate(ig(0:lm,0:lm), ih(0:lm,0:lm))
    call build_index_maps(lm, ig, ih, ncoef_out)
    allocate(coef_out(ncoef_out))
    allocate(ATA(ncoef_out,ncoef_out), ATb(ncoef_out))
    allocate(row_br(ncoef_out), row_bt(ncoef_out), row_bp(ncoef_out))
    ATA(:,:) = 0.0_wp
    ATb(:) = 0.0_wp

    do s=1, nobs
      call schmidt_basis(theta_s(s), phi_s(s), lm, ig, ih, row_br, row_bt, row_bp)
      w = max(w_s(s), 1.0e-12_wp)

      if (joint_w > 0.0_wp) then
        do i=1, ncoef_out
          ATb(i) = ATb(i) + w * (row_br(i)*br_s(s) + joint_w*(row_bt(i)*bt_s(s) + row_bp(i)*bp_s(s)))
          do j=1, ncoef_out
            ATA(i,j) = ATA(i,j) + w * (row_br(i)*row_br(j) + joint_w*(row_bt(i)*row_bt(j) + row_bp(i)*row_bp(j)))
          end do
        end do
      else
        do i=1, ncoef_out
          ATb(i) = ATb(i) + w * row_br(i) * br_s(s)
          do j=1, ncoef_out
            ATA(i,j) = ATA(i,j) + w * row_br(i) * row_br(j)
          end do
        end do
      end if
    end do

    allocate(normal_diagonal(ncoef_out))
    do i=1, ncoef_out
      normal_diagonal(i) = ata(i,i)
    end do
    numerical_floor = epsilon(1.0_wp)*max(sum(normal_diagonal)/real(ncoef_out,wp), tiny(1.0_wp))
    if (reg_lam > 0.0_wp) then
      do l=1, lm
        lnorm = real(l,wp) / real(max(1_int32,lm),wp)
        reg_w = max(lnorm, 1.0e-6_wp)**reg_pow
        do m=0, l
          if (ig(l,m) > 0) ATA(ig(l,m),ig(l,m)) = ATA(ig(l,m),ig(l,m)) * (1.0_wp + reg_lam*reg_w)
          if (ih(l,m) > 0) ATA(ih(l,m),ih(l,m)) = ATA(ih(l,m),ih(l,m)) * (1.0_wp + reg_lam*reg_w)
        end do
      end do
    end if
    do i=1, ncoef_out
      ata(i,i) = ata(i,i) + numerical_floor
    end do

    call solve_linear_system(ATA, ATb, coef_out, ncoef_out, info)
    if (info /= 0) then
      deallocate(ata, atb, row_br, row_bt, row_bp, normal_diagonal)
      return
    end if

    rss_out = 0.0_wp
    do s=1, nobs
      call schmidt_basis(theta_s(s), phi_s(s), lm, ig, ih, row_br, row_bt, row_bp)
      pred_br = sum(row_br * coef_out)
      pred_bt = sum(row_bt * coef_out)
      pred_bp = sum(row_bp * coef_out)
      rr = pred_br - br_s(s)
      rss_out = rss_out + rr*rr
      if (joint_w > 0.0_wp) then
        rr = pred_bt - bt_s(s)
        rss_out = rss_out + joint_w * rr*rr
        rr = pred_bp - bp_s(s)
        rss_out = rss_out + joint_w * rr*rr
      end if
    end do

    ! use the same data-diagonal weighted norm as the fitted damping penalty
    regnorm_out = 0.0_wp
    do l=1, lm
      lnorm = real(l,wp) / real(max(1_int32,lm),wp)
      reg_w = max(lnorm, 1.0e-6_wp)**reg_pow
      do m=0, l
        if (ig(l,m) > 0) regnorm_out = regnorm_out + reg_w * normal_diagonal(ig(l,m)) * coef_out(ig(l,m))**2
        if (ih(l,m) > 0) regnorm_out = regnorm_out + reg_w * normal_diagonal(ih(l,m)) * coef_out(ih(l,m))**2
      end do
    end do

    deallocate(ata, atb, row_br, row_bt, row_bp, normal_diagonal)
  end subroutine solve_coeffs_from_samples

  subroutine auto_select_params(lmax_seed, reg_lam_seed, reg_pow_seed, joint_w, &
                                theta_s, phi_s, w_s, br_s, bt_s, bp_s, &
                                lmax_sel, reg_lam_sel, reg_pow_sel)
    integer(int32), intent(in) :: lmax_seed
    real(wp), intent(in) :: reg_lam_seed, reg_pow_seed, joint_w
    real(wp), intent(in) :: theta_s(:), phi_s(:), w_s(:), br_s(:), bt_s(:), bp_s(:)
    integer(int32), intent(out) :: lmax_sel
    real(wp), intent(out) :: reg_lam_sel, reg_pow_sel

    integer(int32) :: i, il, ilam, ipow, np, stride, nall, ncoef_t, info_t, lm_try
    integer(int32) :: lm_cand(3)
    real(wp) :: lam_cand(4), pow_cand(3)
    real(wp), allocatable :: th_p(:), ph_p(:), ww_p(:), br_p(:), bt_p(:), bp_p(:)
    integer(int32), allocatable :: ig_t(:,:), ih_t(:,:)
    real(wp), allocatable :: coef_t(:)
    real(wp) :: rss, regn, gcv, lcurve, best_gcv, best_lcurve
    real(wp) :: nobs_eff, dof_eff

    nall = size(theta_s)
    stride = max(1_int32, int(sqrt(real(max(1,nall),wp)/600.0_wp)))
    np = 0
    do i=1, nall, stride
      np = np + 1
    end do
    allocate(th_p(np), ph_p(np), ww_p(np), br_p(np), bt_p(np), bp_p(np))
    np = 0
    do i=1, nall, stride
      np = np + 1
      th_p(np) = theta_s(i)
      ph_p(np) = phi_s(i)
      ww_p(np) = w_s(i)
      br_p(np) = br_s(i)
      bt_p(np) = bt_s(i)
      bp_p(np) = bp_s(i)
    end do

    lm_cand(1) = max(6_int32, lmax_seed - 6_int32)
    lm_cand(2) = max(6_int32, lmax_seed)
    lm_cand(3) = max(6_int32, lmax_seed + 6_int32)

    if (reg_lam_seed <= 1.0e-12_wp) then
      lam_cand = [0.0_wp, 2.0e-2_wp, 1.0e-1_wp, 5.0e-1_wp]
    else
      lam_cand = [0.0_wp, max(0.0_wp, 0.25_wp*reg_lam_seed), max(0.0_wp, reg_lam_seed), max(0.0_wp, 1.5_wp*reg_lam_seed)]
    end if
    pow_cand = [max(1.0_wp, reg_pow_seed-2.0_wp), max(1.0_wp, reg_pow_seed), max(1.0_wp, reg_pow_seed+2.0_wp)]

    lmax_sel = lm_cand(2)
    reg_lam_sel = lam_cand(3)
    reg_pow_sel = pow_cand(2)
    best_gcv = huge(1.0_wp)
    best_lcurve = huge(1.0_wp)

    do il=1, 3
      lm_try = lm_cand(il)
      do ilam=1, 4
        do ipow=1, 3
          if (allocated(ig_t)) deallocate(ig_t, ih_t)
          if (allocated(coef_t)) deallocate(coef_t)
          call solve_coeffs_from_samples(lm_try, lam_cand(ilam), pow_cand(ipow), joint_w, &
                                         th_p, ph_p, ww_p, br_p, bt_p, bp_p, &
                                         ig_t, ih_t, coef_t, ncoef_t, rss, regn, info_t)
          if (info_t /= 0) cycle

          nobs_eff = real(size(th_p),wp)
          if (joint_w > 0.0_wp) nobs_eff = 3.0_wp * nobs_eff
          dof_eff = min(real(ncoef_t,wp), max(1.0_wp, nobs_eff-1.0_wp))
          gcv = rss / max(1.0_wp, (nobs_eff - dof_eff)*(nobs_eff - dof_eff))
          lcurve = sqrt(max(rss,1.0e-30_wp)) * sqrt(max(regn,1.0e-30_wp))

          if (gcv < best_gcv) then
            best_gcv = gcv
            best_lcurve = lcurve
            lmax_sel = lm_try
            reg_lam_sel = lam_cand(ilam)
            reg_pow_sel = pow_cand(ipow)
          else if (gcv <= 1.05_wp*best_gcv .and. lcurve < best_lcurve) then
            best_lcurve = lcurve
            lmax_sel = lm_try
            reg_lam_sel = lam_cand(ilam)
            reg_pow_sel = pow_cand(ipow)
          end if
        end do
      end do
    end do

    write(*,'(A,I0,A,ES12.4,A,F7.3)') '  Auto-selected spectral controls: lmax=', lmax_sel, &
                                       ' reg_lambda=', reg_lam_sel, ' reg_power=', reg_pow_sel
    deallocate(th_p, ph_p, ww_p, br_p, bt_p, bp_p)
  end subroutine auto_select_params

  subroutine fit_edge_rbf_correction(lon_s, lat_s, res_br, res_bt, res_bp, lon_ref_deg, poly_lon, poly_lat, &
                                     sigma_deg, band_deg, lambda_corr, centers_per_edge, &
                                     center_lon, center_lat, coef_br, coef_bt, coef_bp, ncorr, info)
    real(wp), intent(in) :: lon_s(:), lat_s(:), res_br(:), res_bt(:), res_bp(:)
    real(wp), intent(in) :: lon_ref_deg, poly_lon(:), poly_lat(:)
    real(wp), intent(in) :: sigma_deg, band_deg, lambda_corr
    integer(int32), intent(in) :: centers_per_edge
    real(wp), allocatable, intent(out) :: center_lon(:), center_lat(:), coef_br(:), coef_bt(:), coef_bp(:)
    integer(int32), intent(out) :: ncorr, info

    integer(int32) :: n, nseg, i, k, i1, i2, c, s, nper
    logical :: closed
    real(wp) :: t, lonu, d2, sigma2, band2, d_edge, w_edge, diag_scale
    real(wp), allocatable :: ATA(:,:), rhs_br(:), rhs_bt(:), rhs_bp(:), psi(:), Awrk(:,:), bwrk(:)

    info = 0
    ncorr = 0
    n = size(poly_lon)
    if (n < 2) then
      info = 1
      return
    end if

    closed = (abs(poly_lon(1)-poly_lon(n)) + abs(poly_lat(1)-poly_lat(n)) < 1.0e-9_wp)
    if (closed) then
      nseg = n - 1
    else
      nseg = n
    end if
    nper = max(1_int32, centers_per_edge)
    ncorr = max(1_int32, nseg * nper)

    allocate(center_lon(ncorr), center_lat(ncorr), coef_br(ncorr), coef_bt(ncorr), coef_bp(ncorr))
    c = 0
    do i=1, nseg
      i1 = i
      if (i < n) then
        i2 = i + 1
      else
        i2 = 1
      end if
      do k=1, nper
        t = (real(k,wp)-0.5_wp) / real(nper,wp)
        c = c + 1
        center_lon(c) = poly_lon(i1) + t*(poly_lon(i2)-poly_lon(i1))
        center_lat(c) = poly_lat(i1) + t*(poly_lat(i2)-poly_lat(i1))
      end do
    end do

    allocate(ATA(ncorr,ncorr), rhs_br(ncorr), rhs_bt(ncorr), rhs_bp(ncorr), psi(ncorr))
    ATA(:,:) = 0.0_wp
    rhs_br(:) = 0.0_wp
    rhs_bt(:) = 0.0_wp
    rhs_bp(:) = 0.0_wp
    sigma2 = max(sigma_deg, 1.0e-3_wp)**2
    band2 = max(band_deg, 1.0e-3_wp)**2

    do s=1, size(lon_s)
      lonu = unwrap_lon_to_ref(lon_s(s), lon_ref_deg)
      d_edge = min_distance_to_poly_edges(lonu, lat_s(s), poly_lon, poly_lat)
      w_edge = exp(-(d_edge*d_edge)/band2)
      if (w_edge < 1.0e-8_wp) cycle
      do c=1, ncorr
        d2 = local_dist2_deg(lonu, lat_s(s), center_lon(c), center_lat(c))
        psi(c) = exp(-0.5_wp*d2/sigma2)
      end do
      do i=1, ncorr
        rhs_br(i) = rhs_br(i) + w_edge*psi(i)*res_br(s)
        rhs_bt(i) = rhs_bt(i) + w_edge*psi(i)*res_bt(s)
        rhs_bp(i) = rhs_bp(i) + w_edge*psi(i)*res_bp(s)
        do k=1, ncorr
          ATA(i,k) = ATA(i,k) + w_edge*psi(i)*psi(k)
        end do
      end do
    end do

    diag_scale = 0.0_wp
    do i=1, ncorr
      diag_scale = diag_scale + ATA(i,i)
    end do
    diag_scale = max(diag_scale/real(max(1_int32,ncorr),wp), 1.0_wp)
    do i=1, ncorr
      ATA(i,i) = ATA(i,i) + max(0.0_wp, lambda_corr)*diag_scale + 1.0e-14_wp
    end do

    allocate(Awrk(ncorr,ncorr), bwrk(ncorr))
    Awrk(:,:) = ATA(:,:)
    bwrk(:) = rhs_br(:)
    call solve_linear_system(Awrk, bwrk, coef_br, ncorr, info)
    if (info /= 0) return
    Awrk(:,:) = ATA(:,:)
    bwrk(:) = rhs_bt(:)
    call solve_linear_system(Awrk, bwrk, coef_bt, ncorr, info)
    if (info /= 0) return
    Awrk(:,:) = ATA(:,:)
    bwrk(:) = rhs_bp(:)
    call solve_linear_system(Awrk, bwrk, coef_bp, ncorr, info)
    if (info /= 0) return

    deallocate(ATA, rhs_br, rhs_bt, rhs_bp, psi, Awrk, bwrk)
  end subroutine fit_edge_rbf_correction

  subroutine evaluate_edge_rbf_correction(lon_deg, lat_deg, lon_ref_deg, poly_lon, poly_lat, &
                                          center_lon, center_lat, coef_br, coef_bt, coef_bp, &
                                          sigma_deg, band_deg, corr_br, corr_bt, corr_bp)
    real(wp), intent(in) :: lon_deg, lat_deg, lon_ref_deg
    real(wp), intent(in) :: poly_lon(:), poly_lat(:)
    real(wp), intent(in) :: center_lon(:), center_lat(:), coef_br(:), coef_bt(:), coef_bp(:)
    real(wp), intent(in) :: sigma_deg, band_deg
    real(wp), intent(out) :: corr_br, corr_bt, corr_bp
    integer(int32) :: c
    real(wp) :: sigma2, band2, lonu, d2, d_edge, w_edge, psi

    corr_br = 0.0_wp
    corr_bt = 0.0_wp
    corr_bp = 0.0_wp
    if (size(center_lon) <= 0) return

    sigma2 = max(sigma_deg, 1.0e-3_wp)**2
    band2 = max(band_deg, 1.0e-3_wp)**2
    lonu = unwrap_lon_to_ref(lon_deg, lon_ref_deg)
    d_edge = min_distance_to_poly_edges(lonu, lat_deg, poly_lon, poly_lat)
    w_edge = exp(-(d_edge*d_edge)/band2)
    if (w_edge < 1.0e-8_wp) return

    do c=1, size(center_lon)
      d2 = local_dist2_deg(lonu, lat_deg, center_lon(c), center_lat(c))
      psi = exp(-0.5_wp*d2/sigma2)
      corr_br = corr_br + coef_br(c) * psi
      corr_bt = corr_bt + coef_bt(c) * psi
      corr_bp = corr_bp + coef_bp(c) * psi
    end do
    corr_br = w_edge * corr_br
    corr_bt = w_edge * corr_bt
    corr_bp = w_edge * corr_bp
  end subroutine evaluate_edge_rbf_correction

  ! dense Gaussian elimination with partial pivoting
  ! this keeps external dependencies minimal for a self-contained build
  subroutine solve_linear_system(A, b, x, n, info)
    integer(int32), intent(in) :: n
    real(wp), intent(inout) :: A(n,n)
    real(wp), intent(inout) :: b(n)
    real(wp), intent(out) :: x(n)
    integer(int32), intent(out) :: info

    integer(int32) :: i, j, k, piv
    real(wp) :: maxv, tmp, factor
    real(wp), allocatable :: rowtmp(:)

    allocate(rowtmp(n))
    info = 0

    do k=1, n-1
      piv = k
      maxv = abs(A(k,k))
      do i=k+1, n
        if (abs(A(i,k)) > maxv) then
          maxv = abs(A(i,k))
          piv = i
        end if
      end do

      ! near-singular pivot indicates ill-conditioned normal equations
      if (maxv < 1.0e-30_wp) then
        info = 1
        deallocate(rowtmp)
        return
      end if

      if (piv /= k) then
        rowtmp(:) = A(k,:)
        A(k,:) = A(piv,:)
        A(piv,:) = rowtmp(:)
        tmp = b(k)
        b(k) = b(piv)
        b(piv) = tmp
      end if

      do i=k+1, n
        factor = A(i,k) / A(k,k)
        A(i,k) = 0.0_wp
        do j=k+1, n
          A(i,j) = A(i,j) - factor*A(k,j)
        end do
        b(i) = b(i) - factor*b(k)
      end do
    end do

    if (abs(A(n,n)) < 1.0e-30_wp) then
      info = 1
      deallocate(rowtmp)
      return
    end if

    x(n) = b(n) / A(n,n)
    do i=n-1, 1, -1
      tmp = b(i)
      do j=i+1, n
        tmp = tmp - A(i,j)*x(j)
      end do
      if (abs(A(i,i)) < 1.0e-30_wp) then
        info = 1
        deallocate(rowtmp)
        return
      end if
      x(i) = tmp / A(i,i)
    end do

    deallocate(rowtmp)
  end subroutine solve_linear_system

end program gravmag_sphere_spectral
