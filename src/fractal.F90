! Copyright (C) 2011-2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_fractal module.

!! This module includes all the conversions of different radii and
!! diameters to volume (and vice versa), such as
!! geometric radius/diameter, mobility equivalent radius/diameter, etc.
!! All equations used in this file are obtained from K.-H. Naumann,
!! COSIMA - a computer program simulating the dynamics
!! of fractal aerosols, Journal of Aerosol Science, Vol. 34, No. 10,
!! pp. 1371-1397, 2003. The equations are written in detail in the file
!! \c doc/fractal/fractal.tex.

!> The fractal_t structure and associated subroutines.
module pmc_fractal

  use pmc_spec_file
  use pmc_constants
  use pmc_netcdf

  !> Constant in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: A_SLIP = 1.142d0
  !> Constant in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: Q_SLIP = 0.588d0
  !> Constant in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: B_SLIP = 0.999d0
  !> Scaling factor in calculating accessible particle surface
  !> in Eq. 26 of Naumann [2003].
  real(kind=dp), parameter :: SCALE_FACTOR_S_ACC = 1d0
  !> Scaling exponent in calculating accessible particle surface
  !> in Eq. 26 of Naumann [2003].
  real(kind=dp), parameter :: SCALE_EXPONENT_S_ACC = 0.86d0

  !> Fractal data.
  !!
  !! The data in this structure is constant, as it represents physical
  !! quantities that cannot change over time.
  type fractal_t
     !> Volume fractal dimension (1).
     real(kind=dp) :: frac_dim
     !> Radius of primary particles (m).
     real(kind=dp) :: prime_radius
     !> Volume filling factor (1).
     real(kind=dp) :: vol_fill_factor
  end type fractal_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate fractal parameters.
  subroutine fractal_allocate(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(out) :: fractal

    fractal%frac_dim = 3d0
    fractal%prime_radius = 1d-8
    fractal%vol_fill_factor = 1d0

  end subroutine fractal_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine fractal_deallocate(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(inout) :: fractal

    fractal%frac_dim = 0d0
    fractal%prime_radius = 0d0
    fractal%vol_fill_factor = 0d0

  end subroutine fractal_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to geometric radius (m) for spherical particles.
  real(kind=dp) elemental function sphere_vol2rad(v)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v

    sphere_vol2rad = (3d0 * v / 4d0 / const%pi)**(1d0 / 3d0)

  end function sphere_vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to number of monomers in a fractal particle cluster.
  !> Based on Eq. 5 in Naumann [2003].
  real(kind=dp) elemental function vol2N(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2N = (sphere_vol2rad(v) / fractal%prime_radius)**3

  end function vol2N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to geometric radius (m).
  real(kind=dp) elemental function vol2rad(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2rad = fractal%prime_radius * (vol2N(v, fractal) &
         * fractal%vol_fill_factor)**(1d0 / fractal%frac_dim)

  end function vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to geometric diameter (m).
  real(kind=dp) elemental function vol2diam(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2diam = 2d0 * vol2rad(v, fractal)

  end function vol2diam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert radius (m) to diameter (m).
  real(kind=dp) elemental function rad2diam(r)

    !> Radius (m).
    real(kind=dp), intent(in) :: r

    rad2diam = 2d0 * r

  end function rad2diam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric radius (m) to volume (m^3) for spherical particles.
  real(kind=dp) elemental function sphere_rad2vol(r)

    !> Radius (m).
    real(kind=dp), intent(in) :: r

    sphere_rad2vol = 4d0 * const%pi * r**3 / 3d0

  end function sphere_rad2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric radius (m) to volume (m^3).
  real(kind=dp) elemental function rad2vol(r, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    rad2vol = sphere_rad2vol(fractal%prime_radius) &
         * (r / fractal%prime_radius)**fractal%frac_dim &
         / fractal%vol_fill_factor

  end function rad2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert diameter (m) to radius (m).
  real(kind=dp) elemental function diam2rad(d)

    !> Diameter (m).
    real(kind=dp), intent(in) :: d

    diam2rad = d / 2d0

  end function diam2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric diameter (m) to volume (m^3).
  real(kind=dp) elemental function diam2vol(d, fractal)

    !> Diameter (m).
    real(kind=dp), intent(in) :: d
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    diam2vol = rad2vol(diam2rad(d), fractal)

  end function diam2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to accessible particle surface (m^2).
  !> Based on Eq. 26 in Naumann [2003].
  real(kind=dp) function vol2S_acc(v, fractal)

    !> Volume (m^3)
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    ! Surface fractal dimension.
    real(kind=dp) :: ds

    if (fractal%frac_dim <= 2d0) then
       ds = 3d0
    elseif ((fractal%frac_dim > 2d0) .and. (fractal%frac_dim <= 3d0)) then
       ds = 6d0 / fractal%frac_dim
    end if

    vol2S_acc = 4d0 * const%pi * fractal%prime_radius**2 &
         * vol2N(v, fractal)**(ds / 3d0) * ((ds - 2d0) &
         * (SCALE_FACTOR_S_ACC / vol2N(v, fractal))**(1d0 &
         - SCALE_EXPONENT_S_ACC) - ds + 3d0)

  end function vol2S_acc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Kirkwood-Riseman ratio.
  !> Based on Eq. 21 in Naumann [2003].
  real(kind=dp) function h_KR(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    h_KR = -0.06483d0 * fractal%frac_dim**2 + 0.6353d0 * &
         fractal%frac_dim - 0.4898d0

  end function h_KR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to continuum regime mobility equivalent radius (m).
  !> Based on Eq. 21 in Naumann [2003].
  real(kind=dp) function vol2R_me_c(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2R_me_c = h_KR(fractal) * vol2rad(v, fractal)

  end function vol2R_me_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to particle effective radius (m).
  !> Based on Eq. 28 in Naumann [2003].
  real(kind=dp) function vol2R_eff(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2R_eff = vol2S_acc(v, fractal) / 4d0 / const%pi  &
         / vol2R_me_c(v, fractal)

  end function vol2R_eff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Slip correction function from continuum to free molecular regime.
  !> Based on Eq. 22 in Naumann [2003].
  real(kind=dp) function Slip_correct(r, tk, press, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    Slip_correct = 1d0 + A_SLIP * air_mean_free_path(tk, press) &
         / r + Q_SLIP * air_mean_free_path(tk, press) / r &
         * exp(-B_SLIP * r / air_mean_free_path(tk, press))

  end function Slip_correct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to mobility equivalent radius (m).
  !> Based on Eq. 5, 21 and 30 in Naumann [2003].
  real(kind=dp) function vol2Rme(v, tk, press, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: x, last_solution, solution
    real(kind=dp) :: Rmec, Reff, C_Reff, fp, f, df
    real(kind=dp) :: a1, a2, a3, a4, a5
    integer :: iter

    if (fractal%frac_dim == 3d0 .and. fractal%vol_fill_factor == 1d0) then
       vol2Rme = vol2rad(v, fractal)
       return
    end if

    Rmec = vol2R_me_c(v, fractal)
    Reff = vol2R_eff(v, fractal)
    C_Reff = Slip_correct(Reff, tk, press, fractal)
    fp = air_mean_free_path(tk, press)
    a1 = C_Reff
    a2 = -Rmec
    a3 = -Rmec * Q_SLIP * fp
    a4 = -B_SLIP / fp
    a5 = - Rmec * A_SLIP * fp

    x = vol2R_me_c(v, fractal)
    do iter = 1,7
       last_solution = x
       f = a1 * x**2 + a2 * x + a3 * exp(a4 * x) + a5
       df = 2d0 * a1 * x + a2 + a3 * a4 * exp(a4 * x)
       x = x - f / df
       solution = x
    end do
    call warn_assert_msg(397562310, &
         almost_equal(last_solution, solution), &
         "volume to Rme newton loop did not satisfy convergence tolerance")

    vol2Rme = x

  end function vol2Rme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function format of mobility equivalent radius.
  !> Based on Eq. 30 in Naumann [2003].
  function f_Rme(Rme, v, tk, press, fractal) result (y)

    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: Rme
    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: y

    real(kind=dp) :: R_me_c, R_eff, C_Reff
    R_me_c = vol2R_me_c(v, fractal)
    R_eff = vol2R_eff(v, fractal)
    C_Reff = Slip_correct(R_eff, tk, press, fractal)

    y = C_Reff * Rme**2 - R_me_c * Rme - R_me_c * Q_SLIP &
         * air_mean_free_path(tk, press) * exp(-B_SLIP * Rme &
         / air_mean_free_path(tk, press)) - R_me_c * A_SLIP &
         * air_mean_free_path(tk, press)

  end function f_Rme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> First derivative format of mobility equivalent radius function.
  function df_Rme(Rme, v, tk, press, fractal) result (y)

    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: Rme
    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    real(kind=dp) :: y

    real(kind=dp) :: R_me_c, R_eff, C_Reff
    R_me_c = vol2R_me_c(v, fractal)
    R_eff = vol2R_eff(v, fractal)
    C_Reff = Slip_correct(R_eff, tk, press, fractal)

    y = 2d0 * C_Reff * Rme - R_me_c + R_me_c * Q_SLIP * B_SLIP &
         * exp(-B_SLIP * Rme / air_mean_free_path(tk, press))

  end function df_Rme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius (m) to that in the continuum regime.
  !> Based on Eq. 30 in Naumann [2003].
  real(kind=dp) function Rme2R_me_c(Rme, tk, press, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: Rme
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: x, last_solution, solution
    real(kind=dp) :: C_Rme, fp, phi, ds, psi, c1, c2, f, df
    real(kind=dp) :: a1, a2, a3, a4, a5, a6, a7, a8
    integer :: iter

    C_Rme = Slip_correct(Rme, tk, press, fractal)
    fp = air_mean_free_path(tk, press)
    if (fractal%frac_dim <= 2d0) then
       ds = 3d0
    elseif ((fractal%frac_dim > 2d0) .and. (fractal%frac_dim <= 3d0)) then
       ds = 6d0 / fractal%frac_dim
    end if
    phi = fractal%prime_radius**2 / (fractal%vol_fill_factor &
         * h_KR(fractal)**fractal%frac_dim * &
         fractal%prime_radius**fractal%frac_dim)**(ds / 3d0)
    psi = 1d0 / (fractal%vol_fill_factor * h_KR(fractal)**fractal%frac_dim &
         * fractal%prime_radius**fractal%frac_dim)**(SCALE_EXPONENT_S_ACC &
         - 1d0)
    c1 = fractal%frac_dim * ds / 3d0 + (SCALE_EXPONENT_S_ACC - 1d0) &
         * fractal%frac_dim - 1d0
    c2 = fractal%frac_dim * ds / 3d0 - 1d0
    
    a1 = C_Rme / Rme * phi * psi * (ds - 2d0)
    a2 = C_Rme / Rme * phi * (3d0 - ds)
    a3 = -phi * psi * (ds - 2d0)
    a4 = phi * (ds - 3d0)
    a5 = -Q_SLIP * fp
    a6 = -B_SLIP / fp * phi * psi * (ds - 2d0)
    a7 = -B_SLIP / fp * phi * (3d0 - ds)
    a8 = -A_SLIP * fp

    x = Rme 
    do iter = 1,7
       last_solution = x
       f = a1 * x**(c1 + 1d0) + a2 * x**(c2 + 1d0) + a3 * x**c1 &
            + a4 * x**c2 + a5 * exp(a6 * x**c1 + a7 * x**c2) + a8
       df = a1 * (c1 + 1d0) * x**c1 + a2 * (c2 + 1d0) * x**c2 &
            + a3 * c1 * x**(c1 - 1d0) + a4 * c2 * x**(c2 - 1d0) &
            + a5 * (a6 * c1 * x**(c1 - 1d0) + a7 * c2 * x**(c2 - 1d0)) &
            * exp(a6 * x**c1 + a7 * x**c2)
       x = x - f / df
       solution = x
    end do
    call warn_assert_msg(397562311, &
            almost_equal(last_solution, solution), &
            "Rme to Rmec newton loop did not satisfy convergence tolerance")

    Rme2R_me_c = x

  end function Rme2R_me_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function format of mobility equivalent radius in continuum regime.
  !> Based on Eq. 30 in Naumann [2003].
  function f_Rmec(Rmec, Rme, tk, press, fractal) result (y)

    !> Mobility equivalent radius in continuum regime (m).
    real(kind=dp), intent(in) :: Rmec
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: Rme
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: y

    real(kind=dp) :: C_Rme, phi, ds, psi, c1, c2
    C_Rme = Slip_correct(Rme, tk, press, fractal)

    if (fractal%frac_dim <= 2d0) then
       ds = 3d0
    elseif ((fractal%frac_dim > 2d0) .and. (fractal%frac_dim <= 3d0)) then
       ds = 6d0 / fractal%frac_dim
    end if

    phi = fractal%prime_radius**2 / (fractal%vol_fill_factor &
         * h_KR(fractal)**fractal%frac_dim * &
         fractal%prime_radius**fractal%frac_dim)**(ds / 3d0)

    psi = 1d0 / (fractal%vol_fill_factor * h_KR(fractal)**fractal%frac_dim &
         * fractal%prime_radius**fractal%frac_dim)**(SCALE_EXPONENT_S_ACC &
         - 1d0)

    c1 = fractal%frac_dim * ds / 3d0 + (SCALE_EXPONENT_S_ACC - 1d0) &
         * fractal%frac_dim - 1d0

    c2 = fractal%frac_dim * ds / 3d0 - 1d0

    y = C_Rme / Rme * phi * psi * (ds - 2d0) * Rmec**(c1 + 1d0) &
         + C_Rme / Rme * phi * (3d0 - ds) * Rmec**(c2 + 1d0) &
         - phi * psi * (ds - 2d0) * Rmec**c1 + phi * (ds - 3d0) * Rmec**c2 &
         - Q_SLIP * air_mean_free_path(tk, press) * exp(-B_SLIP &
         / air_mean_free_path(tk, press) * (phi * psi * (ds - 2d0) &
         * Rmec**c1 + phi * (3d0 - ds) * Rmec**c2)) - A_SLIP &
         * air_mean_free_path(tk, press)

  end function f_Rmec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> First derivative format of mobility equivalent radius in continuum
  !> regime.
  function df_Rmec(Rmec, Rme, tk, press, fractal) result (y)

    !> Mobility equivalent radius in continuum regime (m).
    real(kind=dp), intent(in) :: Rmec
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: Rme
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: y

    real(kind=dp) :: C_Rme, phi, ds, psi, c1, c2
    C_Rme = Slip_correct(Rme, tk, press, fractal)

    if (fractal%frac_dim <= 2d0) then
       ds = 3d0
    elseif ((fractal%frac_dim > 2d0) .and. (fractal%frac_dim <= 3d0)) then
       ds = 6d0 / fractal%frac_dim
    end if

    phi = fractal%prime_radius**2 / (fractal%vol_fill_factor &
         * h_KR(fractal)**fractal%frac_dim * &
         fractal%prime_radius**fractal%frac_dim)**(ds / 3d0)

    psi = 1d0 / (fractal%vol_fill_factor * h_KR(fractal)**fractal%frac_dim &
         * fractal%prime_radius**fractal%frac_dim)**(SCALE_EXPONENT_S_ACC &
         - 1d0)

    c1 = fractal%frac_dim * ds / 3d0 + (SCALE_EXPONENT_S_ACC - 1d0) &
         * fractal%frac_dim - 1d0

    c2 = fractal%frac_dim * ds / 3d0 - 1d0

    y = C_Rme / Rme * phi * psi * (ds - 2d0) * (c1 + 1d0) * Rmec**c1 &
         + C_Rme / Rme * phi * (3d0 - ds) * (c2 + 1d0) * Rmec**c2 &
         - phi * psi * (ds - 2d0) * c1 * Rmec**(c1 - 1d0) &
         + phi * (ds - 3d0) * c2 * Rmec**(c2 - 1d0) &
         - Q_SLIP * air_mean_free_path(tk, press) * exp(-B_SLIP &
         / air_mean_free_path(tk, press) * (phi * psi * (ds - 2d0) * Rmec**c1 &
         + phi * (3d0 - ds) * Rmec**c2)) * (-B_SLIP &
         / air_mean_free_path(tk, press) * (phi * psi * (ds - 2d0) * c1 &
         * Rmec**(c1 - 1d0) + phi * (3d0 - ds) * c2* Rmec**(c2 - 1d0)))

  end function df_Rmec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius (m) to volume (m^3).
  !> Based on Eq. 5, 21 and 30 in Naumann [2003].
  real(kind=dp) function Rme2vol(r, tk, press, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: R_me_c, Rgeo
    if (fractal%frac_dim == 3d0 .and. fractal%vol_fill_factor == 1d0) then
       Rme2vol = rad2vol(r, fractal)
       return
    end if
    
    R_me_c = Rme2R_me_c(r, tk, press, fractal)
    Rgeo = R_me_c / h_KR(fractal)
    Rme2vol = rad2vol(Rgeo, fractal)

  end function Rme2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate air molecular mean free path (m).
  real(kind=dp) function air_mean_free_path(tk, press)

    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press

    real(kind=dp) :: boltz, avogad, mwair, rgas, rhoair, viscosd, &
         viscosk, gasspeed

    ! boltz   = boltzmann's constant (erg/K = g*cm^2/s/K)
    ! avogad  = avogadro's number (molecules/mol)
    ! mwair   = molecular weight of air (g/mol)
    ! rgas    = gas constant (atmos/(mol/liter)/K)
    ! rhoair  = air density (g/cm^3)
    ! viscosd = air dynamic viscosity (g/cm/s)
    ! viscosk = air kinematic viscosity (cm^2/s)
    ! gasspeed    = air molecule mean thermal velocity (cm/s)
    ! air_mean_free_path = air molecule mean free path (m)

    boltz = const%boltzmann * 1d7 ! J/K to erg/K
    avogad = const%avagadro
    mwair = const%air_molec_weight * 1d3 ! kg/mole to g/mole
    rgas = const%univ_gas_const * 1d-2 ! J/mole/K to atmos/(mol/liter)/K

    rhoair = 0.001d0 * ((press/1.01325d5) * mwair / (rgas * tk))

    viscosd = (1.8325d-04 * (296.16d0 + 120d0) / (tk + 120d0)) * (tk &
         / 296.16d0)**1.5d0
    viscosk = viscosd / rhoair
    gasspeed = sqrt(8d0 * boltz * tk * avogad / (const%pi * mwair))
    air_mean_free_path = 2d0 * viscosk / gasspeed * 1d-2

  end function air_mean_free_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the fractal dimension value to avoid an input value greater than 3.
  subroutine check_frac_dim(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(inout) :: fractal

    call assert_msg(801987241, fractal%frac_dim <= 3d0, &
         'fractal dimension greater than 3')

  end subroutine check_frac_dim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read fractal specification from a spec file.
  subroutine spec_file_read_fractal(file, fractal)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Fractal parameters.
    type(fractal_t), intent(inout) :: fractal

    logical :: do_fractal

    !> \page input_format_fractal Input File Format: Fractal Data
    !!
    !! The fractal parameters are divided into those specified at
    !! computed for the rest of the simulation.
    !!
    !! The fractal data file is specified by the parameters:
    !!   - \b frac_dim (real, dimensionless): the fractal dimension
    !!     (3 for spherical and less than 3 for agglomerate)
    !!   - \b prime_radius (real, unit m): radius of primary
    !!     particles
    !!   - \b vol_fill_factor (real, dimensionless): the volume filling
    !!     factor which accounts for the fact that even in a most closely
    !!     packed structure the spherical monomers can occupy only 74%
    !!     of the available volume (1 for compact structure)
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_fractal --- the corresponding output format

    call spec_file_read_logical(file, 'do_fractal', &
         do_fractal)
    if (do_fractal) then
       call spec_file_read_real(file, 'frac_dim', &
            fractal%frac_dim)
       call spec_file_read_real(file, 'prime_radius', &
            fractal%prime_radius)
       call spec_file_read_real(file, 'vol_fill_factor', &
            fractal%vol_fill_factor)
       call check_frac_dim(fractal)
    else
       fractal%frac_dim = 3d0
       fractal%prime_radius = 1d-8
       fractal%vol_fill_factor = 1d0
    end if

  end subroutine spec_file_read_fractal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine fractal_output_netcdf(fractal, ncid)

    !> \page output_format_fractal Output File Format: Fractal Data
    !!
    !! The fractal data NetCDF variables are:
    !!   - \b frac_dim (dimensionless): the fractal dimension
    !!   - \b prime_radius (unit m): radius of primary particles
    !!   - \b vol_fill_factor (dimensionless): volume filling factor
    !!
    !! See also:
    !!   - \ref input_format_fractal --- the corresponding input format

    !> Fractal parameters to write.
    type(fractal_t), intent(in) :: fractal
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_write_real(ncid, fractal%frac_dim, "fractal_dimension", &
         unit="1", standard_name="fractal_dimension")
    call pmc_nc_write_real(ncid, fractal%vol_fill_factor, &
         "volume_filling_factor", unit="1", &
         standard_name="volume_filling_factor")
    call pmc_nc_write_real(ncid, fractal%prime_radius, "prime_radius", &
         unit="m", standard_name="prime_radius")

  end subroutine fractal_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine fractal_input_netcdf(fractal, ncid)

    !> Fractal parameters to read.
    type(fractal_t), intent(inout) :: fractal
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_read_real(ncid, fractal%frac_dim, "fractal_dimension")
    call pmc_nc_read_real(ncid, fractal%vol_fill_factor, &
         "volume_filling_factor")
    call pmc_nc_read_real(ncid, fractal%prime_radius, "prime_radius")

  end subroutine fractal_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_fractal
