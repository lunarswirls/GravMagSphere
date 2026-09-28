module gravmag_sphere_physics
  use, intrinsic :: iso_fortran_env, only: real32
  implicit none
  private
  public :: add_dipole_field_cart, add_pointmass_field_cart, add_mcharge_field_cart

  !***********************************************************************
  ! gravmag_sphere_physics
  !
  ! compact physics kernels used by gravmag_sphere_direct_sum
  !
  ! legacy context:
  ! - old sphere workflows used larger monolithic subroutine trees
  !     with gauss-legendre integration plumbing in many routines
  !
  ! current method:
  ! - keep only two direct cartesian kernels:
  !   1) magnetic dipole kernel from elemental dipole moments m*dV
  !   2) point-mass gravity kernel from elemental masses rho*dV
  ! - this keeps mode switching (mag vs grav) explicit
  !***********************************************************************

  real(real32), parameter :: pi  = 3.14159265358979323846_real32
  real(real32), parameter :: mu0 = 4.0_real32*pi*1.0e-7_real32     ! N/A^2 = T·m/A
  real(real32), parameter :: mu0_over_4pi = 1.0e-7_real32          ! mu0/(4pi)
  real(real32), parameter :: grav_G = 6.67430e-11_real32           ! m^3/(kg s^2)

contains

  !---------------------------------------------------------------------
  ! add_dipole_field_cart
  !
  ! classical magnetic dipole kernel:
  !   B = mu0/(4pi) * ( 3 (m.R) R / |R|^5 - m / |R|^3 )
  !---------------------------------------------------------------------
  subroutine add_dipole_field_cart(xo, yo, zo, xs, ys, zs, mx, my, mz, dBx, dBy, dBz)
    ! Dipole at source (xs,ys,zs) with moment (mx,my,mz) [A m^2]
    ! Field at observer (xo,yo,zo): returns dB [Tesla]
    real(real32), intent(in)  :: xo, yo, zo
    real(real32), intent(in)  :: xs, ys, zs
    real(real32), intent(in)  :: mx, my, mz
    real(real32), intent(out) :: dBx, dBy, dBz

    real(real32) :: Rxv, Ryv, Rzv
    real(real32) :: R2, invR, invR3, invR5
    real(real32) :: mdotR, c

    Rxv = xo - xs
    Ryv = yo - ys
    Rzv = zo - zs

    R2 = Rxv*Rxv + Ryv*Ryv + Rzv*Rzv
    ! singular self-term guard (observer exactly on source point)
    if (R2 <= 0.0_real32) then
      dBx = 0.0_real32; dBy = 0.0_real32; dBz = 0.0_real32
      return
    end if

    invR  = 1.0_real32 / sqrt(R2)
    invR3 = invR*invR*invR
    invR5 = invR3*invR*invR

    mdotR = mx*Rxv + my*Ryv + mz*Rzv
    c     = 3.0_real32 * mdotR * invR5

    dBx = mu0_over_4pi * (c*Rxv - mx*invR3)
    dBy = mu0_over_4pi * (c*Ryv - my*invR3)
    dBz = mu0_over_4pi * (c*Rzv - mz*invR3)
  end subroutine add_dipole_field_cart

  !---------------------------------------------------------------------
  ! add_pointmass_field_cart
  !
  ! newtonian point-mass acceleration:
  !   g = G m R / |R|^3
  ! with R from observer to source
  !---------------------------------------------------------------------
  subroutine add_pointmass_field_cart(xo, yo, zo, xs, ys, zs, mass_kg, dgx, dgy, dgz)
    ! Point-mass gravity contribution.
    ! Source mass at (xs,ys,zs) [kg], observer at (xo,yo,zo) [m]
    ! Returns acceleration components [m/s^2] in planet-fixed XYZ
    real(real32), intent(in)  :: xo, yo, zo
    real(real32), intent(in)  :: xs, ys, zs
    real(real32), intent(in)  :: mass_kg
    real(real32), intent(out) :: dgx, dgy, dgz

    real(real32) :: Rxv, Ryv, Rzv
    real(real32) :: R2, invR, invR3, c

    ! Vector from observer to source: gravity points toward source
    Rxv = xs - xo
    Ryv = ys - yo
    Rzv = zs - zo

    R2 = Rxv*Rxv + Ryv*Ryv + Rzv*Rzv
    if (R2 <= 0.0_real32) then
      dgx = 0.0_real32; dgy = 0.0_real32; dgz = 0.0_real32
      return
    end if

    invR  = 1.0_real32 / sqrt(R2)
    invR3 = invR*invR*invR
    c = grav_G * mass_kg * invR3

    dgx = c * Rxv
    dgy = c * Ryv
    dgz = c * Rzv
  end subroutine add_pointmass_field_cart

  !---------------------------------------------------------------------
  ! add_mcharge_field_cart
  !
  ! equivalent magnetic charge kernel used for the magnetic surface-
  ! charge formulation:
  !   B = mu0/(4pi) * q * R / |R|^3
  !---------------------------------------------------------------------
  subroutine add_mcharge_field_cart(xo, yo, zo, xs, ys, zs, qmag_A, dBx, dBy, dBz)
    ! Equivalent magnetic surface-charge contribution
    ! qmag_A is magnetic charge [A] = (M·n)*dS
    ! Returns dB [Tesla] in cartesian XYZ
    real(real32), intent(in)  :: xo, yo, zo
    real(real32), intent(in)  :: xs, ys, zs
    real(real32), intent(in)  :: qmag_A
    real(real32), intent(out) :: dBx, dBy, dBz

    real(real32) :: Rxv, Ryv, Rzv
    real(real32) :: R2, invR, invR3, c

    Rxv = xo - xs
    Ryv = yo - ys
    Rzv = zo - zs

    R2 = Rxv*Rxv + Ryv*Ryv + Rzv*Rzv
    if (R2 <= 0.0_real32) then
      dBx = 0.0_real32; dBy = 0.0_real32; dBz = 0.0_real32
      return
    end if

    invR  = 1.0_real32 / sqrt(R2)
    invR3 = invR*invR*invR
    c = mu0_over_4pi * qmag_A * invR3

    dBx = c * Rxv
    dBy = c * Ryv
    dBz = c * Rzv
  end subroutine add_mcharge_field_cart

end module gravmag_sphere_physics
