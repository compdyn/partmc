! Copyright (C) 2011-2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_fractal module.

module pmc_fractal

  use pmc_spec_file
  use pmc_constants

  type fractal_t
     !> Whether to do fractal radii conversion.
     logical :: do_fractal 
     !> Constants in slip correction formula.
     real(kind=dp) :: A_slip
     real(kind=dp) :: Q_slip
     real(kind=dp) :: b_slip
     !> Scaling factor in calculating accessible particle surface.
     real(kind=dp) :: scale_factor_S_acc
     !> Scaling exponent in calculating accessible particle surface.
     real(kind=dp) :: scale_exponent_S_acc
     !> Volume fractal dimension.
     real(kind=dp) :: frac_dim
     !> Radius of primary particles (m).
     real(kind=dp) :: prime_radius
     !> Volume filling factor.
     real(kind=dp) :: vol_fill_factor
  end type fractal_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate fractal parameters.
  subroutine fractal_allocate(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(out) :: fractal

    fractal%do_fractal = .false.
    fractal%A_slip = 1.142d0
    fractal%Q_slip = 0.588d0
    fractal%b_slip = 0.999d0
    fractal%scale_factor_S_acc = 1d0
    fractal%scale_exponent_S_acc = 0.86d0
    fractal%frac_dim = 0d0
    fractal%prime_radius = 0d0
    fractal%vol_fill_factor = 0d0

  end subroutine fractal_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine fractal_deallocate(fractal)

    !> Fractal parameters.
    type(fractal_t), intent(inout) :: fractal

    fractal%do_fractal = .false.

  end subroutine fractal_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to radius (m).
  real(kind=dp) elemental function vol2rad(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters. 
    type(fractal_t), intent(in) :: fractal

    if (fractal%do_fractal) then
       vol2rad = vol2Rgeo(v, fractal)
    else
       vol2rad = (v / (4d0 / 3d0 * const%pi))**(1d0 / 3d0)
    end if

  end function vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to diameter (m).
  real(kind=dp) function vol2diam(v, fractal)

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

  !> Convert radius (m) to volume (m^3).
  real(kind=dp) elemental function rad2vol(r, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Fractal parameters. 
    type(fractal_t), intent(in) :: fractal
   
    if (fractal%do_fractal) then
       rad2vol = 4d0 * const%pi * fractal%prime_radius**3d0 * (r &
            / fractal%prime_radius)**fractal%frac_dim / 3d0 &
            / fractal%vol_fill_factor
    else
       rad2vol = 4d0 / 3d0 * const%pi * r**3d0
    end if

  end function rad2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Convert diameter (m) to radius (m).
  real(kind=dp) elemental function diam2rad(d)

    !> Diameter (m).
    real(kind=dp), intent(in) :: d

    diam2rad = d / 2d0

  end function diam2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert diameter (m) to volume (m^3).
  real(kind=dp) function diam2vol(d, fractal)

    !> Diameter (m).
    real(kind=dp), intent(in) :: d
    !> Fractal parameters. 
    type(fractal_t), intent(in) :: fractal

    diam2vol = rad2vol(diam2rad(d), fractal)

  end function diam2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the number of monomers in a fractal particle cluster.
  !> Based on Eq. 5 in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) elemental function vol2N(v, fractal)
   
    !> Volume (m^3)
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2N = 3d0 * v / 4d0 / const%pi / (fractal%prime_radius**3)

  end function vol2N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to geometric radius (m) for fractal particles
  real(kind=dp) elemental function vol2Rgeo(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters. 
    type(fractal_t), intent(in) :: fractal

    vol2Rgeo = fractal%prime_radius * (vol2N(v, fractal) &
         * fractal%vol_fill_factor)**(1d0 / fractal%frac_dim)

  end function vol2Rgeo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the accessible particle surface.
  !> Based on Eq. 26 in Naumann 2003 J. Aerosol. Sci.
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
         * vol2N(v, fractal)**(ds / 3d0)                 &
         * ((ds - 2d0) * (fractal%scale_factor_S_acc     &
         / vol2N(v, fractal))**(1d0                      &
         - fractal%scale_exponent_S_acc) - ds + 3d0)

  end function vol2S_acc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to continuum regime mobility equivalent radius
  real(kind=dp) function vol2R_me_c(v, fractal)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    ! Geometric radius (m).
    real(kind=dp) :: Rgeo

    Rgeo = vol2Rgeo(v, fractal)
    vol2R_me_c = (-0.06483d0 * fractal%frac_dim**2d0 + 0.6353d0 &
         * fractal%frac_dim - 0.4898d0) * Rgeo

  end function vol2R_me_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate effective radius.
  !> Based on Eq. 28 in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function vol2R_eff(v, fractal)

    !> Volume (m^3)
    real(kind=dp), intent(in) :: v
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal

    vol2R_eff = vol2S_acc(v, fractal) / 4d0 / const%pi  &
         / vol2R_me_c(v, fractal)

  end function vol2R_eff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Slip correction function from continuum to free molecular regime
  real(kind=dp) function Slip_correct(r, tk, press, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters. 
    type(fractal_t), intent(in) :: fractal
    
    Slip_correct = 1d0 + fractal%A_slip &
         * air_mean_free_path(tk, press) &
         / r + fractal%Q_slip * air_mean_free_path(tk, press) / r &
         * exp(-fractal%b_slip * r / air_mean_free_path(tk, press))

  end function Slip_correct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate particle mass equivalent radius.
  !> Based on Eq. 3 in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function vol2R_m(v)

    !> Volume (m^3).
    real(kind=dp), intent(in) :: v    

    vol2R_m = (3d0 * v / 4d0 / const%pi)**(1d0 / 3d0)

  end function vol2R_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert volume (m^3) to mobility equivalent radius (m)
  real(kind=dp) function vol2Rme(v, tk, press, fractal)

    !> Volume (m^3)
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters. 
    type(fractal_t), intent(in) :: fractal
    
    real(kind=dp) :: x
    integer, parameter :: MAX_ITERATIONS = 10
    integer :: iter
    iter = 1
    x = vol2R_me_c(v, fractal)

    !print *, 'analytical is ', df_Rme(x, v, tk, press, fractal)
    !print *, 'numerical is ', (f_Rme(x, v, tk, press, fractal) &
    !     - f_Rme(x-1d-10, v, tk, press, fractal)) / 1d-10
    do
      !print *, 'The solution for Rme is ', x
      x = x - f_Rme(x, v, tk, press, fractal) &
           / df_Rme(x, v, tk, press, fractal)
      if (iter > MAX_ITERATIONS) then
         exit
      end if
      iter= iter + 1
    end do

    vol2Rme = x

  end function vol2Rme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function f_Rme(x, v, tk, press, fractal) result (y)
    
    real(kind=dp), intent(in) :: x
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

    y = C_Reff * x**2 - R_me_c * x - R_me_c * fractal%Q_slip &
         * air_mean_free_path(tk, press) * exp(-fractal%b_slip * x &
         / air_mean_free_path(tk, press)) - R_me_c * fractal%A_slip &
         * air_mean_free_path(tk, press)
  end function f_Rme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function df_Rme(x, v, tk, press, fractal) result (y)
    
    real(kind=dp), intent(in) :: x
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

    y = 2d0 * C_Reff * x - R_me_c + R_me_c * fractal%Q_slip * fractal%b_slip &
         * exp(-fractal%b_slip * x / air_mean_free_path(tk, press))
  end function df_Rme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius to that in the continuum regime
  real(kind=dp) function Rme2R_me_c(r, tk, press, fractal)

    !> Radius (m)
    real(kind=dp), intent(in) :: r
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    real(kind=dp) :: x
    integer, parameter :: MAX_ITERATIONS = 10
    integer :: iter
    iter = 1
    x = r

    !print *, 'analytical is ', df_Rmec(x,r,tk,press,fractal)
    !print *, 'numerical is ', (f_Rmec(x,r,tk,press,fractal) &
    !     - f_Rmec(x-1d-10, r, tk, press, fractal))/1d-10
    do
    !  print *, 'The solution for Rme_c is ', x
      x = x - f_Rmec(x, r, tk, press, fractal) &
           / df_Rmec(x, r, tk, press, fractal)
      if(iter > MAX_ITERATIONS) then
         exit
      end if
      iter = iter + 1
    end do

    Rme2R_me_c = x

  end function Rme2R_me_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !> Convert fractal geometric radius (m) to volume (m^3).
  real(kind=dp) function Rgeo2vol(r, fractal)

    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Fractal parameters. 
    type(fractal_t), intent(in) :: fractal

    Rgeo2vol = 4d0 * const%pi * fractal%prime_radius**3d0 * (r &
         / fractal%prime_radius)**fractal%frac_dim / 3d0 &
         / fractal%vol_fill_factor

  end function Rgeo2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function f_Rmec(x, r, tk, press, fractal) result (y)
    
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: r
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    real(kind=dp) :: y

    real(kind=dp) :: C_Rme, phi, h_KR
    C_Rme = Slip_correct(r, tk, press, fractal)
    h_KR = -0.06483d0 * fractal%frac_dim**2 + 0.6353d0 * &
         fractal%frac_dim - 0.4898d0
    phi = fractal%prime_radius**(2d0 - fractal%frac_dim &
         * fractal%scale_exponent_S_acc) &
         / (fractal%vol_fill_factor**fractal%scale_exponent_S_acc &
         * h_KR**(fractal%frac_dim * fractal%scale_exponent_S_acc))
    
    y = C_Rme * x - fractal%A_slip * air_mean_free_path(tk, press) &
         * r / phi * x**(1d0 - fractal%frac_dim                    &
         * fractal%scale_exponent_S_acc)                           &
         - fractal%Q_slip * air_mean_free_path(tk, press) * r      &
         / phi * x**(1d0                                           &   
         - fractal%frac_dim * fractal%scale_exponent_S_acc)        &
         * exp(-fractal%b_slip * phi                               &
         / air_mean_free_path(tk, press)                           &
         * x**(fractal%frac_dim * fractal%scale_exponent_S_acc - 1d0)) - r
  
  end function f_Rmec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function df_Rmec(x, r, tk, press, fractal) result (y)
 
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: r
    !> Temperature (K).
    real(kind=dp), intent(in) :: tk
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    real(kind=dp) :: y

    real(kind=dp) :: C_Rme, phi, h_KR
    C_Rme = Slip_correct(r, tk, press, fractal)
    h_KR = -0.06483d0 * fractal%frac_dim**2 + 0.6353d0 * &
         fractal%frac_dim - 0.4898d0
    phi = fractal%prime_radius**(2d0 - fractal%frac_dim           &
         * fractal%scale_exponent_S_acc)                          & 
         / (fractal%vol_fill_factor**fractal%scale_exponent_S_acc &
         * h_KR**(fractal%frac_dim * fractal%scale_exponent_S_acc))

    y = C_Rme - air_mean_free_path(tk, press) * r / phi                 &
         * (1d0 - fractal%frac_dim * fractal%scale_exponent_S_acc)      &
         * x**(-fractal%frac_dim * fractal%scale_exponent_S_acc)        &
         * (fractal%A_slip + fractal%Q_slip * exp(-fractal%b_slip * phi &
         / air_mean_free_path(tk, press) * x**(fractal%frac_dim         &
         * fractal%scale_exponent_S_acc - 1d0))) - fractal%Q_slip       &
         * fractal%b_slip * r * (1d0 - fractal%frac_dim                 &
         * fractal%scale_exponent_S_acc) / x * exp(-fractal%b_slip      &
         * phi / air_mean_free_path(tk, press)                          &
         * x**(fractal%frac_dim * fractal%scale_exponent_S_acc - 1d0))

  end function df_Rmec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Convert mobility equivalent radius (m) to volume (m^3)
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
    R_me_c = Rme2R_me_c(r, tk, press, fractal)
    Rgeo = R_me_c / (-0.06483d0 * fractal%frac_dim**2 &
         + 0.6353d0 * fractal%frac_dim - 0.4898d0)
    Rme2vol = rad2vol(Rgeo, fractal)
  end function Rme2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read environment specification from a spec file.
  subroutine spec_file_read_fractal(file, fractal)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Fractal parameters.
    type(fractal_t), intent(inout) :: fractal

    !> \page input_format_fractal Input File Format: Fractal State
    !!
    !! The fractal parameters are divided into those specified at
    !! computed for the rest of the simulation.
    !!
    !! The fractal state is specified by the parameters:
    !! - \b frac_dim (real, dimensionless): the fractal dimension
    !!   (3 is spherical and less than 3 is agglomerate)
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_env_state --- the corresponding output
    !!     format
    !!   - \ref input_format_env_data --- the prescribed profiles of
    !!     other environment data

    call spec_file_read_logical(file, 'do_fractal', &
         fractal%do_fractal)
    if (fractal%do_fractal) then
       call spec_file_read_real(file, 'frac_dim', &
            fractal%frac_dim)
       call spec_file_read_real(file, 'prime_radius', &
            fractal%prime_radius)
       call spec_file_read_real(file, 'vol_fill_factor', &
            fractal%vol_fill_factor)
    else
       fractal%frac_dim = 3d0
       fractal%prime_radius = 1d-8 ! Can be set to any value
       fractal%vol_fill_factor = 1d0
    end if

  end subroutine spec_file_read_fractal

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
    air_mean_free_path = 2d0 * viscosk / gasspeed * 1d-02
  
  end function air_mean_free_path  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_fractal
