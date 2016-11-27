! Copyright (C) 2011-2012, 2016 Jian Tian, Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_fractal module.

!> The fractal_t structure and associated subroutines.
!!
!! This module includes all the conversions of different radii and
!! diameters to volume (and vice versa), such as geometric radius/diameter,
!! mobility equivalent radius/diameter, etc. Here volume means the
!! particle mass-equivalent volume (referred to "vol" in code).
!! The mass-equivalent volume is the sum of the per-species volumes, which are
!! in turn the per-species masses divided by the per-species
!! densities in aero_data. Geometric radius
!! (referred to "rad" in code) is defined as the radius of the closest
!! convex envelop. This quantity provides the information about the
!! collisional cross section which is required to calculate coagulation rates.
!! Mobility equivalent radius (referred to "mobility_rad" in code)
!! will be the same as geometric radius for spherical particles, but
!! different for fractal particles. This is due to hydrodynamic interactions
!! between the primary particles, which result in a decrease of the frictional
!! forces acting upon the particles accompanied by an increase of the
!! translational diffusion coefficient. Mobility equivalent radius is often
!! used as output from experimental data such as SMPS measurements. We also
!! use "mobility_rad_in_continuum" for the mobility equivalent radius in
!! continuum regime.
!!
!! All equations used in this file are obtained from
!! K.&nbsp;-H.&nbsp;Naumann (2003) COSIMA - a computer program simulating
!! the dynamics of fractal aerosols, <i>Journal of Aerosol Science</i>,
!! 34(10), 1371-1397. The equations are written in detail in the file
!! \c doc/fractal/fractal.tex.
module pmc_fractal

  use pmc_spec_file
  use pmc_constants
  use pmc_netcdf
  use pmc_mpi

  !> Constant \f$A\f$ in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_A_SLIP = 1.142d0
  !> Constant \f$Q\f$ in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_Q_SLIP = 0.588d0
  !> Constant \f$b\f$ in slip correction equation in Eq. 22 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_B_SLIP = 0.999d0
  !> Scaling factor \f$z\f$ in calculating accessible particle surface area
  !> in Eq. 26 of Naumann [2003].
  real(kind=dp), parameter :: FRACTAL_SCALE_FACTOR_S_ACC = 1d0
  !> Scaling exponent \f$\gamma\f$ in calculating accessible particle surface
  !> area in Eq. 26 of Naumann [2003].
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

  !> Set fractal parameters for spherical particles.
  subroutine fractal_set_spherical(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(out) :: fractal

    fractal%frac_dim = 3d0
    fractal%prime_radius = 1d-8
    fractal%vol_fill_factor = 1d0

  end subroutine fractal_set_spherical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to number of 
  !> monomers \f$N\f$ in a fractal particle cluster.
  !!
  !! Based on Eq. 5 in Naumann [2003].
  real(kind=dp) elemental function fractal_vol_to_num_of_monomers(fractal, v)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    fractal_vol_to_num_of_monomers = (sphere_vol2rad(v) &
         / fractal%prime_radius)**3

  end function fractal_vol_to_num_of_monomers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to geometric radius
  !> \f$R_{\rm geo}\f$ (m).
  real(kind=dp) elemental function fractal_vol2rad(fractal, v)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    if (fractal_is_spherical(fractal)) then
       fractal_vol2rad = sphere_vol2rad(v)
    else
       fractal_vol2rad = fractal%prime_radius &
            * (fractal_vol_to_num_of_monomers(fractal, v) &
            * fractal%vol_fill_factor)**(1d0 / fractal%frac_dim)
    end if

  end function fractal_vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to geometric diameter
  !> \f$D_{\rm geo}\f$ (m).
  real(kind=dp) elemental function fractal_vol2diam(fractal, v)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    fractal_vol2diam = 2d0 * fractal_vol2rad(fractal, v)

  end function fractal_vol2diam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric radius \f$R_{\rm geo}\f$ (m) to mass-equivalent volume
  !> \f$V\f$ (m^3).
  real(kind=dp) elemental function fractal_rad2vol(fractal, r)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Geometric radius (m).
    real(kind=dp), intent(in) :: r

    if (fractal_is_spherical(fractal)) then
       fractal_rad2vol = sphere_rad2vol(r)
    else
       fractal_rad2vol = sphere_rad2vol(fractal%prime_radius) &
            * (r / fractal%prime_radius)**fractal%frac_dim &
            / fractal%vol_fill_factor
    end if

  end function fractal_rad2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric diameter \f$D_{\rm geo}\f$ (m) to
  !> mass-equivalent volume \f$V\f$ (m^3).
  real(kind=dp) elemental function fractal_diam2vol(fractal, d)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Geometric diameter (m).
    real(kind=dp), intent(in) :: d

    fractal_diam2vol = fractal_rad2vol(fractal, diam2rad(d))

  end function fractal_diam2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute surface fractal dimension \f$d_{\rm s}\f$.
  real(kind=dp) function fractal_surface_frac_dim(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    fractal_surface_frac_dim = 0d0 ! prevent uninitialized warnings
    if (fractal%frac_dim <= 2d0) then
       fractal_surface_frac_dim = 3d0
    elseif ((fractal%frac_dim > 2d0) .and. (fractal%frac_dim <= 3d0)) then
       fractal_surface_frac_dim = 6d0 / fractal%frac_dim
    else
       call die_msg(110248362, 'volume fractal dimension larger than 3')
    end if

  end function fractal_surface_frac_dim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to accessible
  !>  particle surface area \f$S_{\rm acc}\f$ (m^2).
  !!
  !! Based on Eq. 26 in Naumann [2003].
  real(kind=dp) function fractal_vol_to_accessible_surface(fractal, v)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Particle mass-equivalent volume (m^3)
    real(kind=dp), intent(in) :: v

    real(kind=dp) :: ds, N

    ds = fractal_surface_frac_dim(fractal)
    N = fractal_vol_to_num_of_monomers(fractal, v)

    fractal_vol_to_accessible_surface = 4d0 * const%pi &
         * fractal%prime_radius**2 * N**(ds / 3d0) * ((ds - 2d0) &
         * (FRACTAL_SCALE_FACTOR_S_ACC / N)&
         **(1d0 - FRACTAL_SCALE_EXPONENT_S_ACC) - ds + 3d0)

  end function fractal_vol_to_accessible_surface

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

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to continuum regime
  !> mobility equivalent radius \f$R_{\rm me,c}\f$ (m).
  !!
  !! Based on Eq. 21 in Naumann [2003].
  real(kind=dp) function fractal_vol_to_mobility_rad_in_continuum(fractal, v)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    fractal_vol_to_mobility_rad_in_continuum = &
         fractal_kirkwood_riseman(fractal) * fractal_vol2rad(fractal, v)

  end function fractal_vol_to_mobility_rad_in_continuum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to particle effective
  !> radius \f$R_{\rm eff}\f$ (m).
  !!
  !! Based on Eq. 28 in Naumann [2003].
  real(kind=dp) function fractal_vol_to_effective_rad(fractal, v)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    fractal_vol_to_effective_rad = &
         fractal_vol_to_accessible_surface(fractal, v) / 4d0 &
         / const%pi / fractal_vol_to_mobility_rad_in_continuum(fractal, v)

  end function fractal_vol_to_effective_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Slip correction function \f$C(R_eff)\f$ from continuum to
  !> free molecular regime.
  !!
  !! Based on Eq. 22 in Naumann [2003].
  real(kind=dp) function fractal_slip_correct(R_eff, temp, pressure)

    !> Effective radius (m).
    real(kind=dp), intent(in) :: R_eff
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp) :: fp

    fp = air_mean_free_path(temp, pressure)
    fractal_slip_correct = 1d0 + FRACTAL_A_SLIP * fp / R_eff &
         + FRACTAL_Q_SLIP * fp / R_eff * exp(-FRACTAL_B_SLIP * R_eff / fp)

  end function fractal_slip_correct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test whether a particle is spherical.
  elemental logical function fractal_is_spherical(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    if (fractal%frac_dim == 3d0 .and. fractal%vol_fill_factor == 1d0) then
       fractal_is_spherical = .true.
    else
       fractal_is_spherical = .false.
    end if

  end function fractal_is_spherical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to mobility equivalent
  !> radius \f$R_{\rm me}\f$ (m).
  !!
  !! Based on Eq. 5, 21 and 30 in Naumann [2003].
  real(kind=dp) function fractal_vol_to_mobility_rad(fractal, v, &
       temp, pressure)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp), parameter :: NEWTON_REL_TOL = 1d-14
    integer, parameter :: NEWTON_MAX_STEPS = 10

    real(kind=dp) :: x
    real(kind=dp) :: Rmec, Reff, C_Reff, fp, f, df
    real(kind=dp) :: a1, a2, a3, a4, a5
    integer ::newton_step

    if (fractal_is_spherical(fractal)) then
       fractal_vol_to_mobility_rad = fractal_vol2rad(fractal, v)
       return
    end if

    Rmec = fractal_vol_to_mobility_rad_in_continuum(fractal, v)
    Reff = fractal_vol_to_effective_rad(fractal, v)
    C_Reff = fractal_slip_correct(Reff, temp, pressure)
    fp = air_mean_free_path(temp, pressure)
    a1 = C_Reff
    a2 = -Rmec
    a3 = -Rmec * FRACTAL_Q_SLIP * fp
    a4 = -FRACTAL_B_SLIP / fp
    a5 = -Rmec * FRACTAL_A_SLIP * fp

    x = Rmec
    do newton_step = 1,NEWTON_MAX_STEPS
       f = a1 * x**2 + a2 * x + a3 * exp(a4 * x) + a5
       df = 2d0 * a1 * x + a2 + a3 * a4 * exp(a4 * x)
       x = x - f / df
       if (abs(f / df) / (abs(x + f / df) + abs(x)) &
            < NEWTON_REL_TOL) exit
    end do
    call assert_msg(728724732, newton_step < NEWTON_MAX_STEPS, &
         "fractal Newton loop failed to converge")
    fractal_vol_to_mobility_rad = x

  end function fractal_vol_to_mobility_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius \f$R_{\rm me}\f$ (m) to that in the
  !> continuum regime \f$R_{\rm me,c}\f$ (m).
  !!
  !! Based on Eq. 30 in Naumann [2003].
  real(kind=dp) function fractal_mobility_rad_to_mobility_rad_in_continuum( &
       fractal, mobility_rad, temp, pressure)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp), parameter :: NEWTON_REL_TOL = 1d-14
    integer, parameter :: NEWTON_MAX_STEPS = 10

    real(kind=dp) :: x
    real(kind=dp) :: C_Rme, fp, phi, ds, psi, c1, c2, f, df
    real(kind=dp) :: a1, a2, a3, a4, a5, a6, a7, a8
    integer :: newton_step

    C_Rme = fractal_slip_correct(mobility_rad, temp, pressure)
    fp = air_mean_free_path(temp, pressure)
    ds = fractal_surface_frac_dim(fractal)
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
    a1 = C_Rme / mobility_rad * phi * psi * (ds - 2d0)
    a2 = C_Rme / mobility_rad * phi * (3d0 - ds)
    a3 = -phi * psi * (ds - 2d0)
    a4 = phi * (ds - 3d0)
    a5 = -FRACTAL_Q_SLIP * fp
    a6 = -FRACTAL_B_SLIP / fp * phi * psi * (ds - 2d0)
    a7 = -FRACTAL_B_SLIP / fp * phi * (3d0 - ds)
    a8 = -FRACTAL_A_SLIP * fp

    x = mobility_rad
    do newton_step = 1,NEWTON_MAX_STEPS
       f = a1 * x**(c1 + 1d0) + a2 * x**(c2 + 1d0) + a3 * x**c1 &
            + a4 * x**c2 + a5 * exp(a6 * x**c1 + a7 * x**c2) + a8
       df = a1 * (c1 + 1d0) * x**c1 + a2 * (c2 + 1d0) * x**c2 &
            + a3 * c1 * x**(c1 - 1d0) + a4 * c2 * x**(c2 - 1d0) &
            + a5 * (a6 * c1 * x**(c1 - 1d0) + a7 * c2 * x**(c2 - 1d0)) &
            * exp(a6 * x**c1 + a7 * x**c2)
       x = x - f / df
       if (abs(f / df) / (abs(x + f / df) + abs(x)) &
            < NEWTON_REL_TOL) exit
    end do
    call assert_msg(875871349, newton_step < NEWTON_MAX_STEPS, &
         "fractal Newton loop failed to converge")
    fractal_mobility_rad_to_mobility_rad_in_continuum = x

  end function fractal_mobility_rad_to_mobility_rad_in_continuum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius \f$R_{\rm me}\f$ (m) to
  !> geometric radius \f$R_{\rm geo}\f$ (m^3).
  real(kind=dp) function fractal_mobility_rad_to_geometric_rad(fractal, &
       mobility_rad, temp, pressure)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp) :: Rmec

    if (fractal_is_spherical(fractal)) then
       fractal_mobility_rad_to_geometric_rad = mobility_rad
       return
    end if

    Rmec = fractal_mobility_rad_to_mobility_rad_in_continuum(fractal, &
         mobility_rad, temp, pressure)
    fractal_mobility_rad_to_geometric_rad &
         = Rmec / fractal_kirkwood_riseman(fractal)

  end function fractal_mobility_rad_to_geometric_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius \f$R_{\rm me}\f$ (m) to
  !> mass-equivalent volume \f$V\f$ (m^3).
  !!
  !! Based on Eq. 5, 21 and 30 in Naumann [2003].
  real(kind=dp) function fractal_mobility_rad_to_vol(fractal, mobility_rad, &
       temp, pressure)

    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    real(kind=dp) :: Rgeo

    Rgeo = fractal_mobility_rad_to_geometric_rad(fractal, mobility_rad, temp, &
         pressure)
    fractal_mobility_rad_to_vol = fractal_rad2vol(fractal, Rgeo)

  end function fractal_mobility_rad_to_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_fractal(val)

    !> Value to pack.
    type(fractal_t), intent(in) :: val

    pmc_mpi_pack_size_fractal = &
         pmc_mpi_pack_size_real(val%frac_dim) &
         + pmc_mpi_pack_size_real(val%prime_radius) &
         + pmc_mpi_pack_size_real(val%vol_fill_factor)

  end function pmc_mpi_pack_size_fractal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_fractal(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(fractal_t), intent(in) :: val

    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%frac_dim)
    call pmc_mpi_pack_real(buffer, position, val%prime_radius)
    call pmc_mpi_pack_real(buffer, position, val%vol_fill_factor)
    call assert(287864891, &
         position - prev_position <= pmc_mpi_pack_size_fractal(val))

  end subroutine pmc_mpi_pack_fractal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_fractal(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(fractal_t), intent(inout) :: val

    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%frac_dim)
    call pmc_mpi_unpack_real(buffer, position, val%prime_radius)
    call pmc_mpi_unpack_real(buffer, position, val%vol_fill_factor)
    call assert(294756245, &
         position - prev_position <= pmc_mpi_pack_size_fractal(val))

  end subroutine pmc_mpi_unpack_fractal

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
    !! The fractal parameters are all held constant for the simulation,
    !! and they are the same for all the particles.
    !!
    !! The fractal data file is specified by the parameters:
    !!   - \b frac_dim \f$d_{\rm f}\f$ (real, dimensionless): the fractal
    !!     dimension (3 for spherical and less than 3 for agglomerate)
    !!   - \b prime_radius \f$R_0\f$ (real, unit m): radius of primary
    !!     particles
    !!   - \b vol_fill_factor \f$f\f$ (real, dimensionless): the volume
    !!     filling factor which accounts for the fact that even
    !!     in a most closely packed structure the spherical monomers can
    !!     occupy only 74% of the available volume (1 for compact structure)
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
       call fractal_set_spherical(fractal)
    end if

  end subroutine spec_file_read_fractal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine fractal_output_netcdf(fractal, ncid)

    !> \page output_format_fractal Output File Format: Fractal Data
    !!
    !! The fractal data NetCDF variables are:
    !!   - \b fractal_dimension \f$d_{\rm f}\f$ (dimensionless):
    !!     particle volume fractal dimension
    !!   - \b fractal_prime_radius \f$R_0\f$ (unit m):
    !!     radius of primary particles
    !!   - \b fractal_vol_fill_factor \f$f\f$ (dimensionless):
    !!     volume filling factor
    !!
    !! See also:
    !!   - \ref input_format_fractal --- the corresponding input format

    !> Fractal parameters to write.
    type(fractal_t), intent(in) :: fractal
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_write_real(ncid, fractal%frac_dim, "fractal_dimension", &
         unit="1", description="particle volume fractal dimension")
    call pmc_nc_write_real(ncid, fractal%prime_radius, &
         "fractal_prime_radius", unit="m", &
         description="radius of primary particles")
    call pmc_nc_write_real(ncid, fractal%vol_fill_factor, &
         "fractal_vol_fill_factor", unit="1", &
         description="volume filling factor")

  end subroutine fractal_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine fractal_input_netcdf(fractal, ncid)

    !> Fractal parameters to read.
    type(fractal_t), intent(inout) :: fractal
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_read_real(ncid, fractal%frac_dim, "fractal_dimension")
    call pmc_nc_read_real(ncid, fractal%prime_radius, "fractal_prime_radius")
    call pmc_nc_read_real(ncid, fractal%vol_fill_factor, &
         "fractal_vol_fill_factor")

  end subroutine fractal_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_fractal
