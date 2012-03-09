! Copyright (C) 2011-2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_chamber module.

module pmc_chamber

  use pmc_aero_data
  use pmc_aero_particle
  use pmc_aero_state
  use pmc_constants
  use pmc_spec_file
  use pmc_fractal

  type chamber_t
     !> Whether to do wall loss and sedimentation in the chamber.
     logical :: do_chamber
     !> Aerosol chamber volume (m^3).
     real(kind=dp) :: V_chamber
     !> Diffusional deposition area (m^2).
     real(kind=dp) :: A_diffuse
     !> Sedimentational deposition area (m^2).
     real(kind=dp) :: A_sedi
     !> Prefactor in dissusive boundary layer thickness (m).
     real(kind=dp) :: prefactor_BL
     !> Exponent in dissusive boundary layer thickness.
     real(kind=dp) :: exponent_BL
     !> Unit translational diffusion coefficient (m^2 s^{-1}).
     real(kind=dp) :: unit_diff_coef = 1d0
  end type chamber_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !> Based on Eq. 23 in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function chamber_diff_coef(aero_particle, fractal, temp, press)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Temperature (K).
    real(kind=dp) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press

    ! Particle effective radius.
    real(kind=dp) :: R_eff

    R_eff = vol2R_eff(aero_particle_volume(aero_particle), fractal)
    chamber_diff_coef = const%boltzmann * temp &
         * Slip_correct(R_eff, temp, press, fractal) / 6d0 / const%pi &
         / const%air_dyn_visc &
         / vol2R_me_c(aero_particle_volume(aero_particle), fractal)

  end function chamber_diff_coef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate diffusional boundary layer thickness.
  !> Based on Eq. 40 in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function chamber_diff_BL_thick(aero_particle, fractal, temp, &
       press)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Temperature (K).
    real(kind=dp) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press

    chamber_diff_BL_thick = chamber%prefactor_BL                      &
         * (chamber_diff_coef(aero_particle, fractal, temp, press)    &
         / chamber%unit_diff_coef)**(chamber%exponent_BL)

  end function chamber_diff_BL_thick

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the loss rate due to wall diffusion in chamber.
  !> Based on Eq. 37a in Naumann 2003 J. Aerosol. Sci. 
  real(kind=dp) function chamber_loss_wall(aero_particle, fractal, temp, press)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Temperature (K).
    real(kind=dp) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press

    chamber_loss_wall = chamber_diff_coef(aero_particle, fractal, temp, press) &
         * chamber%A_diffuse & 
         / chamber_diff_BL_thick(aero_particle, fractal, temp, press) &
         / chamber%V_chamber

  end function chamber_loss_wall

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the loss rate due to sedimentation in chamber.
  !> Based on Eq. 37b in Naumann 2003 J. Aerosol. Sci.
  real(kind=dp) function chamber_loss_sedi(aero_particle, aero_data, &
       fractal, temp, press)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data. 
    type(aero_data_t), intent(in) :: aero_data
    !> Fractal parameters.
    type(fractal_t), intent(in) :: fractal
    !> Temperature (K).
    real(kind=dp) :: temp 
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
   
    ! Particle density.
    real(kind=dp) :: rho

    rho = aero_particle_mass(aero_particle, aero_data)       &
         / aero_particle_volume(aero_particle)

    chamber_loss_sedi = 4d0 * const%pi * rho                        &
         * vol2R_m(aero_particle_volume(aero_particle)**3d0         &
         * const%grav_accel                                         &
         * chamber_diff_coef(aero_particle, fractal, temp, press)   &
         * chamber%A_sedi / 3d0 / const%boltzmann                   &
         / temp / chamber%V_chamber

  end function chamber_loss_sedi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read chamber specification from a spec file.
  subroutine spec_file_read_chamber(file, chamber)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Chamber data.
    type(chamber_t), intent(inout) :: chamber

    !> \page input_format_chamber Input File Format: Chamber
    !!
    !! The chamber parameters are divided into those specified at
    !! the start of the simulation and then either held constant or
    !! computed for the rest of the simulation.
    !! The variables below are for the first type.
    !!
    !! The chamber state is specified by the parameters:
    !! - \b rel_humidity (real, dimensionless): the relative humidity
    !!   (0 is completely unsaturated and 1 is fully saturated)
    !! - \b pressure (real, unit Pa): the atmospheric pressure
    !! - \b latitude (real, unit degrees_north): the latitude of the
    !!   simulation location
    !! - \b longitude (real, unit degrees_east): the longitude of the
    !!   simulation location
    !! - \b altitude (real, unit m): the altitude of the simulation
    !!   location
    !! - \b start_time (real, unit s): the time-of-day of the start of
    !!   the simulation (in seconds past midnight)
    !! - \b start_day (integer): the day-of-year of the start of the
    !!   simulation (starting from 1 on the first day of the year)
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_env_state --- the corresponding output
    !!     format
    !!   - \ref input_format_env_data --- the prescribed profiles of
    !!     other environment data

    call spec_file_read_logical(file, 'do_chamber', &
         chamber%do_chamber)
    if (chamber%do_chamber) then
       call spec_file_read_real(file, 'V_chamber', &
            chamber%V_chamber)
       call spec_file_read_real(file, 'A_diffuse', &
            chamber%A_diffuse)
       call spec_file_read_real(file, 'A_sedi', &
            chamber%A_sedi)
       call spec_file_read_real(file, 'prefactor_BL', &
            chamber%prefactor_BL)
       call spec_file_read_real(file, 'exponent_BL', &
            chamber%exponent_BL)
    end if
  end subroutine spec_file_read_chamber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_chamber
