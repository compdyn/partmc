! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_env_state module.

!> The env_state_t structure and associated subroutines.
module pmc_env_state

  use pmc_gas_state
  use pmc_aero_dist
  use pmc_constants
  use pmc_aero_data
  use pmc_aero_particle
  use pmc_util
  use pmc_gas_data
  use pmc_bin_grid
  use pmc_aero_state
  use pmc_aero_binned
  use pmc_inout
  use pmc_mpi
  use pmc_netcdf
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Current environment state.
  !!
  !! All quantities are instantaneous, describing the state at a
  !! particular instant of time. Constant data and other data not
  !! associated with the current environment state is store in
  !! env_data_t.
  !!
  !! The emissions and dilution are both described by pairs of a state
  !! and a rate. The product of these gives the actual emissions or
  !! dilution with units quantity per time. One way to think about
  !! this is to set the rate to 1/3600 and then regard the state as an
  !! amount per hour, etc.
  type env_state_t
     !> Temperature (K).
     real*8 :: temp
     !> Relative humidity (1).
     real*8 :: rel_humid
     !> Ambient pressure (Pa).
     real*8 :: pressure
     !> Longitude (degrees).
     real*8 :: longitude
     !> Latitude (degrees).
     real*8 :: latitude
     !> Altitude (m).
     real*8 :: altitude
     !> Start time (s since 00:00 UTC).
     real*8 :: start_time
     !> Start day of year (UTC).
     integer :: start_day
     !> Elapsed time since start_time (s).
     real*8 :: elapsed_time
     !> Box height (m).
     real*8 :: height
     !> Gas emissions.
     type(gas_state_t) :: gas_emissions
     !> Gas emisssion rate (s^{-1}).
     real*8 :: gas_emission_rate
     !> Background gas concentrations.
     type(gas_state_t) :: gas_background
     !> Gas-background dilution rate (s^{-1}).
     real*8 :: gas_dilution_rate
     !> Aerosol emissions.
     type(aero_dist_t) :: aero_emissions
     !> Aerosol emisssion rate (s^{-1}).
     real*8 :: aero_emission_rate
     !> Aerosol background.
     type(aero_dist_t) :: aero_background
     !> Aero-background dilute rate (s^{-1}).
     real*8 :: aero_dilution_rate
  end type env_state_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an empty environment.
  subroutine env_state_alloc(env_state)

    !> Environment.
    type(env_state_t), intent(out) :: env_state

    env_state%temp = 0d0
    env_state%rel_humid = 0d0
    env_state%pressure = 0d0
    env_state%longitude = 0d0
    env_state%latitude = 0d0
    env_state%altitude = 0d0
    env_state%start_time = 0d0
    env_state%start_day = 0
    env_state%elapsed_time = 0d0
    env_state%height = 0d0

    call gas_state_alloc(env_state%gas_emissions, 0)
    call gas_state_alloc(env_state%gas_background, 0)
    env_state%gas_emission_rate = 0d0
    env_state%gas_dilution_rate = 0d0
    call aero_dist_alloc(env_state%aero_emissions, 0, 0)
    call aero_dist_alloc(env_state%aero_background, 0, 0)
    env_state%aero_emission_rate = 0d0
    env_state%aero_dilution_rate = 0d0

  end subroutine env_state_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine env_state_free(env_state)

    !> Environment.
    type(env_state_t), intent(out) :: env_state

    call gas_state_free(env_state%gas_emissions)
    call gas_state_free(env_state%gas_background)
    call aero_dist_free(env_state%aero_emissions)
    call aero_dist_free(env_state%aero_background)

  end subroutine env_state_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> env_state += env_state_delta
  subroutine env_state_add(env_state, env_state_delta)

    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Increment.
    type(env_state_t), intent(in) :: env_state_delta

    env_state%temp = env_state%temp + env_state_delta%temp
    env_state%rel_humid = env_state%rel_humid + env_state_delta%rel_humid
    env_state%pressure = env_state%pressure + env_state_delta%pressure
    env_state%longitude = env_state%longitude + env_state_delta%longitude
    env_state%latitude = env_state%latitude + env_state_delta%latitude
    env_state%altitude = env_state%altitude + env_state_delta%altitude
    env_state%start_time = env_state%start_time + env_state_delta%start_time
    env_state%start_day = env_state%start_day + env_state_delta%start_day
    env_state%elapsed_time = env_state%elapsed_time &
         + env_state_delta%elapsed_time
    env_state%height = env_state%height + env_state_delta%height
    call gas_state_add(env_state%gas_emissions, env_state_delta%gas_emissions)
    env_state%gas_emission_rate = env_state%gas_emission_rate &
         + env_state_delta%gas_emission_rate
    call gas_state_add(env_state%gas_background, &
         env_state_delta%gas_background)
    env_state%gas_dilution_rate = env_state%gas_dilution_rate &
         + env_state_delta%gas_dilution_rate
    env_state%aero_emission_rate = env_state%aero_emission_rate &
         + env_state_delta%aero_emission_rate
    env_state%aero_dilution_rate = env_state%aero_dilution_rate &
         + env_state_delta%aero_dilution_rate
    
  end subroutine env_state_add
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> env_state *= alpha
  subroutine env_state_scale(env_state, alpha)

    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Scale factor.
    real*8, intent(in) :: alpha

    env_state%temp = env_state%temp * alpha
    env_state%rel_humid = env_state%rel_humid * alpha
    env_state%pressure = env_state%pressure * alpha
    env_state%longitude = env_state%longitude * alpha
    env_state%latitude = env_state%latitude * alpha
    env_state%altitude = env_state%altitude * alpha
    env_state%start_time = env_state%start_time * alpha
    env_state%start_day = nint(dble(env_state%start_day) * alpha)
    env_state%elapsed_time = env_state%elapsed_time * alpha
    env_state%height = env_state%height * alpha
    call gas_state_scale(env_state%gas_emissions, alpha)
    env_state%gas_emission_rate = env_state%gas_emission_rate * alpha
    call gas_state_scale(env_state%gas_background, alpha)
    env_state%gas_dilution_rate = env_state%gas_dilution_rate * alpha
    
  end subroutine env_state_scale
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> env_to = env_from
  subroutine env_state_copy(env_from, env_to)

    !> Original.
    type(env_state_t), intent(in) :: env_from
    !> Destination.
    type(env_state_t), intent(inout) :: env_to

    env_to%temp = env_from%temp
    env_to%rel_humid = env_from%rel_humid
    env_to%pressure = env_from%pressure
    env_to%longitude = env_from%longitude
    env_to%latitude = env_from%latitude
    env_to%altitude = env_from%altitude
    env_to%start_time = env_from%start_time
    env_to%start_day = env_from%start_day
    env_to%elapsed_time = env_from%elapsed_time
    env_to%height = env_from%height
    call gas_state_copy(env_from%gas_emissions, env_to%gas_emissions)
    env_to%gas_emission_rate = env_from%gas_emission_rate
    call gas_state_copy(env_from%gas_background, env_to%gas_background)
    env_to%gas_dilution_rate = env_from%gas_dilution_rate
    call aero_dist_copy(env_from%aero_emissions, env_to%aero_emissions)
    env_to%aero_emission_rate = env_from%aero_emission_rate
    call aero_dist_copy(env_from%aero_background, env_to%aero_background)
    env_to%aero_dilution_rate = env_from%aero_dilution_rate
    
  end subroutine env_state_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given water volume to the water vapor and updates all
  !> environment quantities.
  subroutine env_state_change_water_volume(env_state, aero_data, dv)
    
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Aero_data constants.
    type(aero_data_t), intent(in) :: aero_data
    !> Conc of water added (m^3/m^3).
    real*8, intent(in) :: dv
    
    real*8 pmv     ! ambient water vapor pressure (Pa)
    real*8 mv      ! ambient water vapor density (kg m^{-3})
                   ! pmv and mv are related by the factor molec_weight/(R*T)
    real*8 dmv     ! change of water density (kg m^{-3})
    
    dmv = dv * aero_data%density(aero_data%i_water)
    pmv = env_state_sat_vapor_pressure(env_state) * env_state%rel_humid
    mv = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const*env_state%temp) * pmv
    mv = mv - dmv    
    env_state%rel_humid = const%univ_gas_const * env_state%temp &
         / aero_data%molec_weight(aero_data%i_water) * mv &
         / env_state_sat_vapor_pressure(env_state)
    
  end subroutine env_state_change_water_volume
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the current saturation vapor pressure (Pa).
  real*8 function env_state_sat_vapor_pressure(env_state)
    
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    
    env_state_sat_vapor_pressure = const%water_eq_vap_press &
         * 10d0**(7.45d0 * (env_state%temp - const%water_freeze_temp) &
         / (env_state%temp - 38d0))
    
  end function env_state_sat_vapor_pressure
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the critical relative humidity from the kappa value (1).
  real*8 function aero_particle_kappa_rh(aero_particle, aero_data, &
       env_state)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real*8 :: kappa, diam, C, A
    
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    A = 4d0 * const%water_surf_eng * const%water_molec_weight &
         / (const%univ_gas_const * env_state%temp * const%water_density)
    C = sqrt(4d0 * A**3 / 27d0)
    diam = vol2diam(aero_particle_volume(aero_particle))
    aero_particle_kappa_rh = C / sqrt(kappa * diam**3) + 1d0

  end function aero_particle_kappa_rh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Air density (kg m^{-3}).
  real*8 function env_state_air_den(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_air_den = const%air_molec_weight &
         * env_state_air_molar_den(env_state)

  end function env_state_air_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Air molar density (mole m^{-3}).
  real*8 function env_state_air_molar_den(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_air_molar_den = env_state%pressure &
         / (const%univ_gas_const * env_state%temp)

  end function env_state_air_molar_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert (mole m^{-3}) to (ppb).
  subroutine gas_state_mole_dens_to_ppb(gas_state, env_state)

    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    
    gas_state%conc = gas_state%conc / env_state_air_molar_den(env_state) * 1d9
    
  end subroutine gas_state_mole_dens_to_ppb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment.
  subroutine env_state_update_gas_state(env_state, delta_t, &
       old_env_state, &
       gas_data, gas_state)

    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Time increment to update over.
    real*8, intent(in) :: delta_t
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state to update.
    type(gas_state_t), intent(inout) :: gas_state

    real*8 :: effective_dilution_rate
    type(gas_state_t) :: emission, dilution

    call gas_state_alloc(emission, gas_data%n_spec)
    call gas_state_alloc(dilution, gas_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env_state%gas_dilution_rate
    if (env_state%height > old_env_state%height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env_state%height - old_env_state%height) / delta_t / &
            old_env_state%height
    end if

    ! emission = delta_t * gas_emission_rate * gas_emissions
    ! but emissions are in (mole m^{-2} s^{-1})
    call gas_state_copy(env_state%gas_emissions, emission)
    call gas_state_scale(emission, 1d0 / env_state%height)
    call gas_state_mole_dens_to_ppb(emission, env_state)
    call gas_state_scale(emission, delta_t * env_state%gas_emission_rate)

    ! dilution = delta_t * gas_dilution_rate * (gas_background - gas_state)
    call gas_state_copy(env_state%gas_background, dilution)
    call gas_state_sub(dilution, gas_state)
    call gas_state_scale(dilution, delta_t * effective_dilution_rate)

    call gas_state_add(gas_state, emission)
    call gas_state_add(gas_state, dilution)

    !FIXME: ensure gas state is still positive

    call gas_state_free(emission)
    call gas_state_free(dilution)

  end subroutine env_state_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment for a
  !> particle aerosol distribution.
  subroutine env_state_update_aero_state(env_state, delta_t, old_env_state, &
       bin_grid, aero_data, aero_state, aero_binned)

    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Time increment to update over.
    real*8, intent(in) :: delta_t
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero state to update.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aero binned to update.
    type(aero_binned_t), intent(inout) :: aero_binned

    integer :: i
    real*8 :: sample_prop, effective_dilution_rate
    type(aero_state_t) :: aero_state_delta
    type(aero_binned_t) :: aero_binned_delta

    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state_delta)
    call aero_binned_alloc(aero_binned_delta, bin_grid%n_bin, aero_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env_state%aero_dilution_rate
    if (env_state%height > old_env_state%height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env_state%height - old_env_state%height) / delta_t / &
            old_env_state%height
    end if

    ! loss to background
    sample_prop = delta_t * effective_dilution_rate
    if (sample_prop > 1d0) then
       write(0,*) 'ERROR: effective dilution rate too high for this timestep'
       call exit(1)
    end if
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_state_sample(aero_state, aero_state_delta, sample_prop)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_binned_sub(aero_binned, aero_binned_delta)

    ! addition from background
    sample_prop = delta_t * effective_dilution_rate
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_state_add_aero_dist_sample(aero_state_delta, bin_grid, &
         aero_data, env_state%aero_background, sample_prop, &
         env_state%elapsed_time)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add_particles(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)
    
    ! emissions
    sample_prop = delta_t * env_state%aero_emission_rate / env_state%height
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_state_add_aero_dist_sample(aero_state_delta, bin_grid, &
         aero_data, env_state%aero_emissions, sample_prop, &
         env_state%elapsed_time)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add_particles(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)

    ! update computational volume
    aero_state%comp_vol = aero_state%comp_vol * env_state%temp &
         / old_env_state%temp

    call aero_state_free(aero_state_delta)
    call aero_binned_free(aero_binned_delta)

  end subroutine env_state_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment for a
  !> binned aerosol distribution.
  subroutine env_state_update_aero_binned(env_state, delta_t, & 
       old_env_state, &
       bin_grid, aero_data, aero_binned)

    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Time increment to update over.
    real*8, intent(in) :: delta_t
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero binned to update.
    type(aero_binned_t), intent(inout) :: aero_binned

    type(aero_binned_t) :: emission, dilution
    real*8 :: effective_dilution_rate

    call aero_binned_alloc(emission, bin_grid%n_bin, aero_data%n_spec)
    call aero_binned_alloc(dilution, bin_grid%n_bin, aero_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env_state%aero_dilution_rate
    if (env_state%height > old_env_state%height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env_state%height - old_env_state%height) / delta_t / &
            old_env_state%height
    end if

    ! emission = delta_t * aero_emission_rate * aero_emissions
    ! but emissions are #/m^2 so we need to divide by height
    call aero_binned_add_aero_dist(emission, bin_grid, aero_data, &
         env_state%aero_emissions)
    call aero_binned_scale(emission, &
         delta_t * env_state%aero_emission_rate / env_state%height)

    ! dilution = delta_t * aero_dilution_rate * (aero_background - aero_binned)
    call aero_binned_add_aero_dist(dilution, bin_grid, aero_data, &
         env_state%aero_background)
    call aero_binned_sub(dilution, aero_binned)
    call aero_binned_scale(dilution, delta_t * effective_dilution_rate)

    call aero_binned_add(aero_binned, emission)
    call aero_binned_add(aero_binned, dilution)

    call aero_binned_free(emission)
    call aero_binned_free(dilution)

  end subroutine env_state_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine inout_write_env_state(file, env_state)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Environment to write.
    type(env_state_t), intent(in) :: env_state
    
    call inout_write_comment(file, "begin env_state")
    call inout_write_real(file, "temp(K)", env_state%temp)
    call inout_write_real(file, "rel_humidity(1)", env_state%rel_humid)
    call inout_write_real(file, "pressure(Pa)", env_state%pressure)
    call inout_write_real(file, "longitude(deg)", env_state%longitude)
    call inout_write_real(file, "latitude(deg)", env_state%latitude)
    call inout_write_real(file, "altitude(m)", env_state%altitude)
    call inout_write_real(file, "start_time(s)", env_state%start_time)
    call inout_write_integer(file, "start_day(days)", env_state%start_day)
    call inout_write_real(file, "elapsed_time(s)", env_state%elapsed_time)
    call inout_write_real(file, "height(m)", env_state%height)
    call inout_write_gas_state(file, env_state%gas_emissions)
    call inout_write_real(file, "gas_emit_rate(1/s)", &
         env_state%gas_emission_rate)
    call inout_write_gas_state(file, env_state%gas_background)
    call inout_write_real(file, "gas_dilute_rate(1/s)", &
         env_state%gas_dilution_rate)
    call inout_write_aero_dist(file, env_state%aero_emissions)
    call inout_write_real(file, "aero_emit_rate(1/s)", &
         env_state%aero_emission_rate)
    call inout_write_aero_dist(file, env_state%aero_background)
    call inout_write_real(file, "aero_dilute_rat(1/s)", &
         env_state%aero_dilution_rate)
    call inout_write_comment(file, "end env_state")

  end subroutine inout_write_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine inout_read_env_state(file, env_state)
    
    !> File to read from.
    type(inout_file_t), intent(inout) :: file
    !> Environment to read.
    type(env_state_t), intent(out) :: env_state
    
    call inout_check_comment(file, "begin env_state")
    call inout_read_real(file, "temp(K)", env_state%temp)
    call inout_read_real(file, "rel_humidity(1)", env_state%rel_humid)
    call inout_read_real(file, "pressure(Pa)", env_state%pressure)
    call inout_read_real(file, "longitude(deg)", env_state%longitude)
    call inout_read_real(file, "latitude(deg)", env_state%latitude)
    call inout_read_real(file, "altitude(m)", env_state%altitude)
    call inout_read_real(file, "start_time(s)", env_state%start_time)
    call inout_read_integer(file, "start_day(days)", env_state%start_day)
    call inout_read_real(file, "elapsed_time(s)", env_state%elapsed_time)
    call inout_read_real(file, "height(m)", env_state%height)
    call inout_read_gas_state(file, env_state%gas_emissions)
    call inout_read_real(file, "gas_emit_rate(1/s)", &
         env_state%gas_emission_rate)
    call inout_read_gas_state(file, env_state%gas_background)
    call inout_read_real(file, "gas_dilute_rate(1/s)", &
         env_state%gas_dilution_rate)
    call inout_read_aero_dist(file, env_state%aero_emissions)
    call inout_read_real(file, "aero_emit_rate(1/s)", &
         env_state%aero_emission_rate)
    call inout_read_aero_dist(file, env_state%aero_background)
    call inout_read_real(file, "aero_dilute_rat(1/s)", &
         env_state%aero_dilution_rate)
    call inout_check_comment(file, "end env_state")

  end subroutine inout_read_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read environment specification from a inout file.
  subroutine spec_read_env_state(file, env_state)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Environment data.
    type(env_state_t), intent(out) :: env_state

    call env_state_alloc(env_state)
    call inout_read_real(file, 'rel_humidity', env_state%rel_humid)
    call inout_read_real(file, 'pressure', env_state%pressure)
    call inout_read_real(file, 'latitude', env_state%latitude)
    call inout_read_real(file, 'longitude', env_state%longitude)
    call inout_read_real(file, 'altitude', env_state%altitude)
    call inout_read_real(file, 'start_time', env_state%start_time)
    call inout_read_integer(file, 'start_day', env_state%start_day)

  end subroutine spec_read_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of an array of env_state.
  subroutine env_state_average(env_vec, env_avg)

    !> Array of env_state.
    type(env_state_t), intent(in) :: env_vec(:)
    !> Average of env_vec.
    type(env_state_t), intent(out) :: env_avg

    call average_real(env_vec%temp, env_avg%temp)
    call average_real(env_vec%rel_humid, env_avg%rel_humid)
    call average_real(env_vec%pressure, env_avg%pressure)
    call average_real(env_vec%longitude, env_avg%longitude)
    call average_real(env_vec%latitude, env_avg%latitude)
    call average_real(env_vec%altitude, env_avg%altitude)
    call average_real(env_vec%start_time, env_avg%start_time)
    call average_integer(env_vec%start_day, env_avg%start_day)
    call average_real(env_vec%elapsed_time, env_avg%elapsed_time)
    call average_real(env_vec%height, env_avg%height)
    call gas_state_average(env_vec%gas_emissions, env_avg%gas_emissions)
    call average_real(env_vec%gas_emission_rate, env_avg%gas_emission_rate)
    call gas_state_average(env_vec%gas_background, env_avg%gas_background)
    
  end subroutine env_state_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes.
  subroutine env_state_mix(val)

    !> Value to average.
    type(env_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(env_state_t) :: val_avg

    call env_state_alloc(val_avg)
    call pmc_mpi_allreduce_average_real(val%temp, val_avg%temp)
    call pmc_mpi_allreduce_average_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_allreduce_average_real(val%pressure, val_avg%pressure)
    val%temp = val_avg%temp
    val%rel_humid = val_avg%rel_humid
    val%pressure = val_avg%pressure
    call env_state_free(val_avg)
#endif

  end subroutine env_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_env_state(val)

    !> Value to pack.
    type(env_state_t), intent(in) :: val

    pmc_mpi_pack_size_env_state = &
         pmc_mpi_pack_size_real(val%temp) &
         + pmc_mpi_pack_size_real(val%rel_humid) &
         + pmc_mpi_pack_size_real(val%pressure) &
         + pmc_mpi_pack_size_real(val%longitude) &
         + pmc_mpi_pack_size_real(val%latitude) &
         + pmc_mpi_pack_size_real(val%altitude) &
         + pmc_mpi_pack_size_real(val%start_time) &
         + pmc_mpi_pack_size_integer(val%start_day) &
         + pmc_mpi_pack_size_real(val%elapsed_time) &
         + pmc_mpi_pack_size_real(val%height) &
         + pmc_mpi_pack_size_gas_state(val%gas_emissions) &
         + pmc_mpi_pack_size_real(val%gas_emission_rate) &
         + pmc_mpi_pack_size_gas_state(val%gas_background) &
         + pmc_mpi_pack_size_real(val%gas_dilution_rate) &
         + pmc_mpi_pack_size_aero_dist(val%aero_emissions) &
         + pmc_mpi_pack_size_real(val%aero_emission_rate) &
         + pmc_mpi_pack_size_aero_dist(val%aero_background) &
         + pmc_mpi_pack_size_real(val%aero_dilution_rate)

  end function pmc_mpi_pack_size_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_env_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(env_state_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real(buffer, position, val%temp)
    call pmc_mpi_pack_real(buffer, position, val%rel_humid)
    call pmc_mpi_pack_real(buffer, position, val%pressure)
    call pmc_mpi_pack_real(buffer, position, val%longitude)
    call pmc_mpi_pack_real(buffer, position, val%latitude)
    call pmc_mpi_pack_real(buffer, position, val%altitude)
    call pmc_mpi_pack_real(buffer, position, val%start_time)
    call pmc_mpi_pack_integer(buffer, position, val%start_day)
    call pmc_mpi_pack_real(buffer, position, val%elapsed_time)
    call pmc_mpi_pack_real(buffer, position, val%height)
    call pmc_mpi_pack_gas_state(buffer, position, val%gas_emissions)
    call pmc_mpi_pack_real(buffer, position, val%gas_emission_rate)
    call pmc_mpi_pack_gas_state(buffer, position, val%gas_background)
    call pmc_mpi_pack_real(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_pack_aero_dist(buffer, position, val%aero_emissions)
    call pmc_mpi_pack_real(buffer, position, val%aero_emission_rate)
    call pmc_mpi_pack_aero_dist(buffer, position, val%aero_background)
    call pmc_mpi_pack_real(buffer, position, val%aero_dilution_rate)
    call assert(464101191, &
         position - prev_position == pmc_mpi_pack_size_env_state(val))
#endif

  end subroutine pmc_mpi_pack_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_env_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(env_state_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real(buffer, position, val%temp)
    call pmc_mpi_unpack_real(buffer, position, val%rel_humid)
    call pmc_mpi_unpack_real(buffer, position, val%pressure)
    call pmc_mpi_unpack_real(buffer, position, val%longitude)
    call pmc_mpi_unpack_real(buffer, position, val%latitude)
    call pmc_mpi_unpack_real(buffer, position, val%altitude)
    call pmc_mpi_unpack_real(buffer, position, val%start_time)
    call pmc_mpi_unpack_integer(buffer, position, val%start_day)
    call pmc_mpi_unpack_real(buffer, position, val%elapsed_time)
    call pmc_mpi_unpack_real(buffer, position, val%height)
    call pmc_mpi_unpack_gas_state(buffer, position, val%gas_emissions)
    call pmc_mpi_unpack_real(buffer, position, val%gas_emission_rate)
    call pmc_mpi_unpack_gas_state(buffer, position, val%gas_background)
    call pmc_mpi_unpack_real(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_emissions)
    call pmc_mpi_unpack_real(buffer, position, val%aero_emission_rate)
    call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_background)
    call pmc_mpi_unpack_real(buffer, position, val%aero_dilution_rate)
    call assert(205696745, &
         position - prev_position == pmc_mpi_pack_size_env_state(val))
#endif

  end subroutine pmc_mpi_unpack_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_env_state(val, val_avg)

    !> Value to average.
    type(env_state_t), intent(in) :: val
    !> Result.
    type(env_state_t), intent(out) :: val_avg

    call env_state_alloc(val_avg)
    call env_state_copy(val, val_avg)
    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)

  end subroutine pmc_mpi_reduce_avg_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine env_state_output_netcdf(env_state, ncid)
    
    !> Environment state to write.
    type(env_state_t), intent(in) :: env_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_write_real(ncid, env_state%temp, "temperature", "K")
    call pmc_nc_write_real(ncid, env_state%rel_humid, "relative_humidity", "1")
    call pmc_nc_write_real(ncid, env_state%pressure, "pressure", "Pa")
    call pmc_nc_write_real(ncid, env_state%longitude, "longitude", "degrees")
    call pmc_nc_write_real(ncid, env_state%latitude, "latitude", "degrees")
    call pmc_nc_write_real(ncid, env_state%altitude, "altitude", "m")
    call pmc_nc_write_real(ncid, env_state%start_time, &
         "start_time_of_day", "s")
    call pmc_nc_write_integer(ncid, env_state%start_day, &
         "start_day_of_year", "1")
    call pmc_nc_write_real(ncid, env_state%elapsed_time, "elapsed_time", "s")
    call pmc_nc_write_real(ncid, env_state%height, "height", "m")

  end subroutine env_state_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env_state
