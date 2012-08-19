! Copyright (C) 2011-2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_fractal module.

!> The fractal_t structure and associated subroutines.
!!
!! This module includes all the conversions of different radii and
!! diameters to volume (and vice versa), such as
!! geometric radius/diameter, mobility equivalent radius/diameter, etc.
!! Here volume means the particle material volume (referred to "vol" in
!! code). Geometric radius (referred to "rad" in code) is defined as the
!! radius of the closest convex envelop.
!! This quantity provides the information about the collisional
!! cross section which is required to calculate coagulation rates.
!! Mobility equivalent radius (referred to "mobility_rad" in code)
!! will be the same as geometric radius for
!! spherical particles, but different for fractal particles. This is due
!! to hydrodynamic interactions between the primary particles, which
!! result in a decrease of the frictional forces acting upon the particles
!! accompanied by an increase of the translational diffusion coefficient.
!! Mobility equivalent radius is often used as output from experimental 
!! data such as SMPS measurement. Also use "mobility_rad_in_continuum"
!! to stand for mobility equivalent radius in continuum regime.
!!
!! All equations used in this file are obtained from
!! K.&nbsp;-H.&nbsp;Naumann (2003)
!! COSIMA - a computer program simulating the dynamics
!! of fractal aerosols, <i>Journal of Aerosol Science</i>, 34(10),
!! 1371-1397. The equations are written in detail in the file
!! \c doc/fractal/fractal.tex.
module pmc_fractal

  use pmc_spec_file
  use pmc_constants
  use pmc_netcdf

  !> Constant \f$A\f$ in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_A_SLIP = 1.142d0
  !> Constant \f$Q\f$ in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_Q_SLIP = 0.588d0
  !> Constant \f$b\f$ in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_B_SLIP = 0.999d0
  !> Scaling factor \f$z\f$ in calculating accessible particle surface
  !> in Eq. 26 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_SCALE_FACTOR_S_ACC = 1d0
  !> Scaling exponent \f$\gamma\f$ in calculating accessible particle surface
  !> in Eq. 26 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_SCALE_EXPONENT_S_ACC = 0.86d0

  !> Fractal data.
  !!
  !! The data in this structure is constant, as it represents physical
  !! quantities that cannot change over time.
  type fractal_t
     !> Volume fractal dimension \f$d_{\rm f}\f$ (1).
     real(kind=dp) :: frac_dim
     !> Radius of primary particles \f$R_{\rm 0}\f$ (m).
     real(kind=dp) :: prime_radius
     !> Volume filling factor \f$f\f$ (1).
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

  !> Convert volume \f$V\f$ (m^3) to geometric radius \f$R_{\rm geo}\f$ (m)
  !> for spherical particles.
  real(kind=dp) elemental function sphere_vol2rad(v)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v

    sphere_vol2rad = (3d0 * v / 4d0 / const%pi)**(1d0 / 3d0)

  end function sphere_vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume \f$V\f$ (m^3) to number of monomers \f$N\f$ in a
  !> fractal particle cluster.
  !!
  !! Based on Eq. 5 in Naumann [2003].
  real(kind=dp) elemental function vol_to_num_of_monomers(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol_to_num_of_monomers = (sphere_vol2rad(v) / fractal%prime_radius)**3

  end function vol_to_num_of_monomers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume \f$V\f$ (m^3) to geometric radius \f$R_{\rm geo}\f$ (m).
  real(kind=dp) elemental function vol2rad(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2rad = fractal%prime_radius * (vol_to_num_of_monomers(v, fractal) &
         * fractal%vol_fill_factor)**(1d0 / fractal%frac_dim)

  end function vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume \f$V\f$ (m^3) to geometric diameter \f$D_{\rm geo}\f$ (m).
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

  !> Convert geometric radius \f$R_{\rm geo}\f$ (m) to volume \f$V\f$ (m^3)
  !> for spherical particles.
  real(kind=dp) elemental function sphere_rad2vol(r)

    !> Radius (m).
    real(kind=dp), intent(in) :: r

    sphere_rad2vol = 4d0 * const%pi * r**3 / 3d0

  end function sphere_rad2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric radius \f$R_{\rm geo}\f$ (m) to volume \f$V\f$ (m^3).
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

  !> Convert geometric diameter \f$D_{\rm geo}\f$ (m) to volume \f$V\f$ (m^3).
  real(kind=dp) elemental function diam2vol(d, fractal)

    !> Diameter (m).
    real(kind=dp), intent(in) :: d
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    diam2vol = rad2vol(diam2rad(d), fractal)

  end function diam2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute surface fractal dimension \f$d_{\rm s}\f$.
  real(kind=dp) function fractal_cal_surface_frac_dim(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    if (fractal%frac_dim <= 2d0) then
       fractal_cal_surface_frac_dim = 3d0
    elseif ((fractal%frac_dim > 2d0) .and. (fractal%frac_dim <= 3d0)) then
       fractal_cal_surface_frac_dim = 6d0 / fractal%frac_dim
    else
       call die_msg(110248362, 'volume fractal dimension larger than 3')
    end if

  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume \f$V\f$ (m^3) to accessible particle surface
  !>  \f$S_{\rm acc}\f$ (m^2).
  !!
  !! Based on Eq. 26 in Naumann [2003].
  real(kind=dp) function vol_to_accessible_surface(v, fractal)

    !> Volume (m^3)
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    ! Surface fractal dimension.
    real(kind=dp) :: ds

    ds = fractal_cal_surface_frac_dim(fractal)

    vol_to_accessible_surface = 4d0 * const%pi * fractal%prime_radius**2 &
         * vol_to_num_of_monomers(v, fractal)**(ds / 3d0) &
         * ((ds - 2d0) * (FRACTAL_SCALE_FACTOR_S_ACC &
         / vol_to_num_of_monomers(v, fractal))**(1d0 &
         - FRACTAL_SCALE_EXPONENT_S_ACC) - ds + 3d0)

  end function vol_to_accessible_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Kirkwood-Riseman ratio \f$h_{\rm KR}\f$.
  !!
  !! Based on Eq. 21 in Naumann [2003].
  real(kind=dp) function fractal_kirkwood_riseman(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    fractal_kirkwood_riseman = -0.06483d0 * fractal%frac_dim**2 &
         + 0.6353d0 * fractal%frac_dim - 0.4898d0

  end function fractal_kirkwood_riseman

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume \f$V\f$ (m^3) to continuum regime mobility equivalent
  !> radius \f$R_{\rm me,c}\f$ (m).
  !!
  !! Based on Eq. 21 in Naumann [2003].
  real(kind=dp) function vol_to_mobility_rad_in_continuum(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol_to_mobility_rad_in_continuum = fractal_kirkwood_riseman(fractal) &
         * vol2rad(v, fractal)

  end function vol_to_mobility_rad_in_continuum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume \f$V\f$ (m^3) to particle effective radius
  !> \f$R_{\rm eff}\f$ (m).
  !!
  !! Based on Eq. 28 in Naumann [2003].
  real(kind=dp) function vol_to_effective_rad(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol_to_effective_rad = vol_to_accessible_surface(v, fractal) / 4d0 &
         / const%pi / vol_to_mobility_rad_in_continuum(v, fractal)

  end function vol_to_effective_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Slip correction function \f$C(R)\f$ from continuum to
  !> free molecular regime.
  !!
  !! Based on Eq. 22 in Naumann [2003].
  real(kind=dp) function fractal_slip_correct(r, tk, press, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: fp

    fp = air_mean_free_path(tk, press)
    fractal_slip_correct = 1d0 + FRACTAL_A_SLIP * fp / r &
         + FRACTAL_Q_SLIP * fp / r * exp(-FRACTAL_B_SLIP * r / fp)

  end function fractal_slip_correct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume \f$V\f$ (m^3) to mobility equivalent radius
  !> \f$R_{\rm me}\f$ (m).
  !!
  !! Based on Eq. 5, 21 and 30 in Naumann [2003].
  real(kind=dp) function vol_to_mobility_rad(v, tk, press, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    integer, parameter :: MAX_ITERATIONS = 10

    real(kind=dp) :: x
    integer :: iter

    if (fractal%frac_dim == 3d0 .and. fractal%vol_fill_factor == 1d0) then
       x = vol2rad(v, fractal)
    else
       x = vol_to_mobility_rad_in_continuum(v, fractal)
       do iter = 1, MAX_ITERATIONS
          x = x - fractal_f_mobility_rad(x, v, tk, press, fractal) &
              / fractal_df_mobility_rad(x, v, tk, press, fractal)
       end do
    end if
    vol_to_mobility_rad = x

  end function vol_to_mobility_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function format of mobility equivalent radius \f$R_{\rm me}\f$.
  !!
  !! Helper function. Do not call directly. To be solved in
  !! vol_to_mobility_rad.
  !!
  !! Based on Eq. 30 in Naumann [2003].
  real(kind=dp) function fractal_f_mobility_rad(mobility_rad, v, &
       tk, press, fractal)

    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: R_me_c, R_eff, C_Reff, fp

    R_me_c = vol_to_mobility_rad_in_continuum(v, fractal)
    R_eff = vol_to_effective_rad(v, fractal)
    C_Reff = fractal_slip_correct(R_eff, tk, press, fractal)
    fp = air_mean_free_path(tk, press)

    fractal_f_mobility_rad = C_Reff * mobility_rad**2 &
         - R_me_c * mobility_rad &
         - R_me_c * FRACTAL_Q_SLIP * fp &
         * exp(-FRACTAL_B_SLIP * mobility_rad / fp) &
         - R_me_c * FRACTAL_A_SLIP * fp

  end function fractal_f_mobility_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> First derivative format of mobility equivalent radius \f$R_{\rm me}\f$.
  !!
  !! Helper function. Do not call directly. To be solved in
  !! vol_to_mobility_rad.
  real(kind=dp) function fractal_df_mobility_rad(mobility_rad, v, &
       tk, press, fractal)

    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: R_me_c, R_eff, C_Reff

    R_me_c = vol_to_mobility_rad_in_continuum(v, fractal)
    R_eff = vol_to_effective_rad(v, fractal)
    C_Reff = fractal_slip_correct(R_eff, tk, press, fractal)

    fractal_df_mobility_rad = 2d0 * C_Reff * mobility_rad &
         - R_me_c &
         + R_me_c * FRACTAL_Q_SLIP * FRACTAL_B_SLIP &
         * exp(-FRACTAL_B_SLIP * mobility_rad / air_mean_free_path(tk, press))

  end function fractal_df_mobility_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius \f$R_{\rm me}\f$ (m) to that in the
  !> continuum regime \f$R_{\rm me,c}\f$ (m).
  !!
  !! Based on Eq. 30 in Naumann [2003].
  real(kind=dp) function mobility_rad_to_mobility_rad_in_continuum( &
       mobility_rad, tk, press, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) ::mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    integer, parameter :: MAX_ITERATIONS = 10

    real(kind=dp) :: x
    integer :: iter

    x = mobility_rad
    do iter = 1, MAX_ITERATIONS
      x = x - fractal_f_mobility_rad_in_continuum(x, mobility_rad, &
           tk, press, fractal) &
           / fractal_df_mobility_rad_in_continuum(x, mobility_rad, &
           tk, press, fractal)
    end do
    mobility_rad_to_mobility_rad_in_continuum = x

  end function mobility_rad_to_mobility_rad_in_continuum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function format of mobility equivalent radius
  !> in continuum regime \f$R_{\rm me,c}\f$.
  !!
  !! Helper function. Do not call directly. To be solved in
  !! mobility_rad_to_mobility_rad_in_continuum.
  !!
  !! Based on Eq. 30 in Naumann [2003].
  real(kind=dp) function fractal_f_mobility_rad_in_continuum( &
       mobility_rad_in_cont, mobility_rad, tk, press, fractal)

    !> Mobility equivalent radius in continuum regime (m).
    real(kind=dp), intent(in) :: mobility_rad_in_cont
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: C_Rme, fp, phi, ds, psi, c1, c2

    C_Rme = fractal_slip_correct(mobility_rad, tk, press, fractal)
    fp = air_mean_free_path(tk, press)
    ds = fractal_cal_surface_frac_dim(fractal)
    phi = fractal%prime_radius**2 / (fractal%vol_fill_factor &
         * fractal_kirkwood_riseman(fractal)**fractal%frac_dim &
         * fractal%prime_radius**fractal%frac_dim)**(ds / 3d0)
    psi = 1d0 / (fractal%vol_fill_factor &
         * fractal_kirkwood_riseman(fractal)**fractal%frac_dim &
         * fractal%prime_radius**fractal%frac_dim)&
         **(FRACTAL_SCALE_EXPONENT_S_ACC - 1d0)
    c1 = fractal%frac_dim * ds / 3d0 + (FRACTAL_SCALE_EXPONENT_S_ACC - 1d0) &
         * fractal%frac_dim - 1d0
    c2 = fractal%frac_dim * ds / 3d0 - 1d0

    fractal_f_mobility_rad_in_continuum = C_Rme / mobility_rad * phi &
         * psi * (ds - 2d0) * mobility_rad_in_cont**(c1 + 1d0) &
         + C_Rme / mobility_rad * phi * (3d0 - ds) &
         * mobility_rad_in_cont**(c2 + 1d0) &
         - phi * psi * (ds - 2d0) * mobility_rad_in_cont**c1 &
         + phi * (ds - 3d0) * mobility_rad_in_cont**c2 &
         - FRACTAL_Q_SLIP * fp * exp(-FRACTAL_B_SLIP / fp &
         * (phi * psi * (ds - 2d0) * mobility_rad_in_cont**c1 &
         + phi * (3d0 - ds) * mobility_rad_in_cont**c2)) &
         - FRACTAL_A_SLIP * fp

  end function fractal_f_mobility_rad_in_continuum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> First derivative format of mobility equivalent radius in continuum
  !> regime \f$R_{\rm me,c}\f$.
  !!
  !! Helper function. Do not call directly. To be solved in
  !! mobility_rad_to_mobility_rad_in_continuum.
  real(kind=dp) function fractal_df_mobility_rad_in_continuum( &
       mobility_rad_in_cont, mobility_rad, tk, press, fractal)

    !> Mobility equivalent radius in continuum regime (m).
    real(kind=dp), intent(in) :: mobility_rad_in_cont
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: C_Rme, fp, phi, ds, psi, c1, c2

    C_Rme = fractal_slip_correct(mobility_rad, tk, press, fractal)
    fp = air_mean_free_path(tk, press)
    ds = fractal_cal_surface_frac_dim(fractal)
    phi = fractal%prime_radius**2 / (fractal%vol_fill_factor &
         * fractal_kirkwood_riseman(fractal)**fractal%frac_dim * &
         fractal%prime_radius**fractal%frac_dim)**(ds / 3d0)
    psi = 1d0 / (fractal%vol_fill_factor &
         * fractal_kirkwood_riseman(fractal)**fractal%frac_dim &
         * fractal%prime_radius**fractal%frac_dim)&
         **(FRACTAL_SCALE_EXPONENT_S_ACC - 1d0)
    c1 = fractal%frac_dim * ds / 3d0 + (FRACTAL_SCALE_EXPONENT_S_ACC - 1d0) &
         * fractal%frac_dim - 1d0
    c2 = fractal%frac_dim * ds / 3d0 - 1d0

    fractal_df_mobility_rad_in_continuum = C_Rme / mobility_rad * phi * psi &
         * (ds - 2d0) * (c1 + 1d0) * mobility_rad_in_cont**c1 &
         + C_Rme / mobility_rad * phi * (3d0 - ds) * (c2 + 1d0) &
         * mobility_rad_in_cont**c2 &
         - phi * psi * (ds - 2d0) * c1 * mobility_rad_in_cont**(c1 - 1d0) &
         + phi * (ds - 3d0) * c2 * mobility_rad_in_cont**(c2 - 1d0) &
         - FRACTAL_Q_SLIP * fp * exp(-FRACTAL_B_SLIP / fp &
         * (phi * psi * (ds - 2d0) * mobility_rad_in_cont**c1 &
         + phi * (3d0 - ds) * mobility_rad_in_cont**c2)) &
         * (-FRACTAL_B_SLIP / fp * (phi * psi &
         * (ds - 2d0) * c1 * mobility_rad_in_cont**(c1 - 1d0) &
         + phi * (3d0 - ds) * c2* mobility_rad_in_cont**(c2 - 1d0)))

  end function fractal_df_mobility_rad_in_continuum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius \f$R_{\rm me}\f$ (m) to volume
  !> \f$V\f$ (m^3).
  !!
  !! Based on Eq. 5, 21 and 30 in Naumann [2003].
  real(kind=dp) function mobility_rad_to_vol(mobility_rad, tk, press, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    real(kind=dp) :: R_me_c, Rgeo

    if (fractal%frac_dim == 3d0 .and. fractal%vol_fill_factor == 1d0) then
       mobility_rad_to_vol = rad2vol(mobility_rad, fractal)
    else
       R_me_c = mobility_rad_to_mobility_rad_in_continuum(mobility_rad, tk, &
            press, fractal)
       Rgeo = R_me_c / fractal_kirkwood_riseman(fractal)
       mobility_rad_to_vol = rad2vol(Rgeo, fractal)
    end if

  end function mobility_rad_to_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate air molecular mean free path \f$l\f$ (m).
  real(kind=dp) function air_mean_free_path(tk, press)

    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press

    real(kind=dp) :: boltz, avogad, mwair, rgas, rhoair, viscosd, &
         viscosk, gasspeed

    boltz = const%boltzmann
    avogad = const%avagadro
    mwair = const%air_molec_weight
    rgas = const%univ_gas_const

    rhoair = (press * mwair) / (rgas * tk)

    viscosd = (1.8325d-05 * (296.16d0 + 120d0) / (tk + 120d0)) &
         * (tk / 296.16d0)**1.5d0
    viscosk = viscosd / rhoair
    gasspeed = sqrt(8d0 * boltz * tk * avogad / (const%pi * mwair))
    air_mean_free_path = 2d0 * viscosk / gasspeed

  end function air_mean_free_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       call assert_msg(801987241, fractal%frac_dim <= 3d0, &
            'fractal dimension greater than 3')
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
         "fractal_vol_fill_factor", unit="1", &
         standard_name="fractal_vol_fill_factor")
    call pmc_nc_write_real(ncid, fractal%prime_radius, &
         "fractal_prime_radius", unit="m", &
         standard_name="fractal_prime_radius")

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
         "fractal_vol_fill_factor")
    call pmc_nc_read_real(ncid, fractal%prime_radius, "fractal_prime_radius")

  end subroutine fractal_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_fractal
