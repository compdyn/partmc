! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Environment parameters.
!
! The temperature profile is proscribed as a function of time by
! giving a number of times and the corresponding temperatures. Linear
! interpolation is used between the times, with constant interpolation
! outside of the range of times.

module pmc_env

  use pmc_gas_state
  use pmc_aero_dist
  
  type env_t
     real*8 :: temp                     ! temperature (K)
     real*8 :: rel_humid                ! relative humidity (1)
     real*8 :: pressure                 ! ambient pressure (Pa)
     real*8 :: longitude                ! longitude (degrees)
     real*8 :: latitude                 ! latitude (degrees)
     real*8 :: altitude                 ! altitude (m)
     real*8 :: start_time               ! start time (s since 00:00 UTC)
     integer :: start_day               ! start day of year (UTC)
     real*8 :: height                   ! box height (m)
     type(gas_state_t) :: gas_emissions ! gas emissions
     real*8 :: gas_emission_rate        ! gas emisssion rate (s^{-1})
     type(gas_state_t) :: gas_background ! background gas concentrations
     real*8 :: gas_dilution_rate        ! gas-background dilution rate (s^{-1})
     type(aero_dist_t) :: aero_emissions ! aerosol emissions
     real*8 :: aero_emission_rate       ! aerosol emisssion rate (s^{-1})
     type(aero_dist_t) :: aero_background ! aerosol background
     real*8 :: aero_dilution_rate       ! aero-background dilution rate (s^{-1})
  end type env_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_alloc(env)

    ! Allocate an empty environment.

    type(env_t), intent(out) :: env   ! environment

    env%temp = 0d0
    env%rel_humid = 0d0
    env%pressure = 0d0
    env%longitude = 0d0
    env%latitude = 0d0
    env%altitude = 0d0
    env%start_time = 0d0
    env%start_day = 0
    env%height = 0d0

    call gas_state_alloc(env%gas_emissions, 0)
    call gas_state_alloc(env%gas_background, 0)
    env%gas_emission_rate = 0d0
    env%gas_dilution_rate = 0d0
    call aero_dist_alloc(env%aero_emissions, 0, 0, 0)
    call aero_dist_alloc(env%aero_background, 0, 0, 0)
    env%aero_emission_rate = 0d0
    env%aero_dilution_rate = 0d0

  end subroutine env_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_free(env)

    ! Free all storage.

    type(env_t), intent(out) :: env   ! environment

    call gas_state_free(env%gas_emissions)
    call gas_state_free(env%gas_background)
    call aero_dist_free(env%aero_emissions)
    call aero_dist_free(env%aero_background)

  end subroutine env_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_add(env, env_delta)
    
    ! env += env_delta

    use pmc_gas_state
    use pmc_aero_dist
    
    type(env_t), intent(inout) :: env   ! environment
    type(env_t), intent(in) :: env_delta ! increment

    env%temp = env%temp + env_delta%temp
    env%rel_humid = env%rel_humid + env_delta%rel_humid
    env%pressure = env%pressure + env_delta%pressure
    env%longitude = env%longitude + env_delta%longitude
    env%latitude = env%latitude + env_delta%latitude
    env%altitude = env%altitude + env_delta%altitude
    env%start_time = env%start_time + env_delta%start_time
    env%start_day = env%start_day + env_delta%start_day
    env%height = env%height + env_delta%height
    call gas_state_add(env%gas_emissions, env_delta%gas_emissions)
    env%gas_emission_rate = env%gas_emission_rate + env_delta%gas_emission_rate
    call gas_state_add(env%gas_background, env_delta%gas_background)
    env%gas_dilution_rate = env%gas_dilution_rate + env_delta%gas_dilution_rate
    call aero_dist_add(env%aero_emissions, env_delta%aero_emissions)
    env%aero_emission_rate = env%aero_emission_rate &
         + env_delta%aero_emission_rate
    call aero_dist_add(env%aero_background, env_delta%aero_background)
    env%aero_dilution_rate = env%aero_dilution_rate &
         + env_delta%aero_dilution_rate
    
  end subroutine env_add
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_scale(env, alpha)
    
    ! env *= alpha

    use pmc_gas_state
    use pmc_aero_dist
    
    type(env_t), intent(inout) :: env   ! environment
    real*8, intent(in) :: alpha         ! scale factor

    env%temp = env%temp * alpha
    env%rel_humid = env%rel_humid * alpha
    env%pressure = env%pressure * alpha
    env%longitude = env%longitude * alpha
    env%latitude = env%latitude * alpha
    env%altitude = env%altitude * alpha
    env%start_time = env%start_time * alpha
    env%start_day = nint(dble(env%start_day) * alpha)
    env%height = env%height * alpha
    call gas_state_scale(env%gas_emissions, alpha)
    env%gas_emission_rate = env%gas_emission_rate * alpha
    call gas_state_scale(env%gas_background, alpha)
    env%gas_dilution_rate = env%gas_dilution_rate * alpha
    call aero_dist_scale(env%aero_emissions, alpha)
    env%aero_emission_rate = env%aero_emission_rate * alpha
    call aero_dist_scale(env%aero_background, alpha)
    env%aero_dilution_rate = env%aero_dilution_rate * alpha
    
  end subroutine env_scale
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_copy(env_from, env_to)
    
    ! env_to = env_from

    use pmc_gas_state
    use pmc_aero_dist
    
    type(env_t), intent(in) :: env_from ! original
    type(env_t), intent(inout) :: env_to ! destination

    env_to%temp = env_from%temp
    env_to%rel_humid = env_from%rel_humid
    env_to%pressure = env_from%pressure
    env_to%longitude = env_from%longitude
    env_to%latitude = env_from%latitude
    env_to%altitude = env_from%altitude
    env_to%start_time = env_from%start_time
    env_to%start_day = env_from%start_day
    env_to%height = env_from%height
    call gas_state_copy(env_from%gas_emissions, env_to%gas_emissions)
    env_to%gas_emission_rate = env_from%gas_emission_rate
    call gas_state_copy(env_from%gas_background, env_to%gas_background)
    env_to%gas_dilution_rate = env_from%gas_dilution_rate
    call aero_dist_copy(env_from%aero_emissions, env_to%aero_emissions)
    env_to%aero_emission_rate = env_from%aero_emission_rate
    call aero_dist_copy(env_from%aero_background, env_to%aero_background)
    env_to%aero_dilution_rate = env_from%aero_dilution_rate
    
  end subroutine env_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine env_change_water_volume(env, aero_data, dv)
    
    ! Adds the given water volume to the water vapor and updates all
    ! environment quantities.
    
    use pmc_constants
    use pmc_aero_data
    
    type(env_t), intent(inout) :: env   ! environment state to update
    type(aero_data_t), intent(in) :: aero_data ! aero_data constants
    real*8, intent(in) :: dv            ! conc of water added (m^3/m^3)
    
    real*8 pmv     ! ambient water vapor pressure (Pa)
    real*8 mv      ! ambient water vapor density (kg m^{-3})
                   ! pmv and mv are related by the factor molec_weight/(R*T)
    real*8 dmv     ! change of water density (kg m^{-3})
    
    dmv = dv * aero_data%density(aero_data%i_water)
    pmv = env_sat_vapor_pressure(env) * env%rel_humid
    mv = aero_data%molec_weight(aero_data%i_water)/(const%R*env%temp) * pmv
    mv = mv - dmv    
    env%rel_humid = const%R * env%temp &
         / aero_data%molec_weight(aero_data%i_water) * mv &
         / env_sat_vapor_pressure(env)
    
  end subroutine env_change_water_volume
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real*8 function env_sat_vapor_pressure(env) ! (Pa)

    ! Computes the current saturation vapor pressure.
    
    use pmc_constants
    
    type(env_t), intent(in) :: env      ! environment state
    
    env_sat_vapor_pressure = const%p00 * 10d0**(7.45d0 * (env%temp - const%T0) &
         / (env%temp - 38d0)) ! Pa
    
  end function env_sat_vapor_pressure
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_kappa_rh(aero_particle, aero_data, env) ! (1)

    ! Returns the critical relative humidity from the kappa value.

    use pmc_aero_particle
    use pmc_constants
    use pmc_aero_data
    use pmc_util

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(env_t), intent(in) :: env      ! environment state

    real*8 :: kappa, diam, C, A
    
    kappa = aero_particle_solute_kappa(aero_data, aero_particle)
    A = 4d0 * const%sig * const%water_molec_weight &
         / (const%R * env%temp * const%water_density)
    C = sqrt(4d0 * A**3 / 27d0)
    diam = vol2diam(aero_particle_volume(aero_particle))
    aero_particle_kappa_rh = C / sqrt(kappa * diam**3) + 1d0

  end function aero_particle_kappa_rh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function env_air_den(env) ! (kg m^{-3})

    use pmc_constants

    type(env_t), intent(in) :: env      ! environment state

    env_air_den = const%M_a * env_air_molar_den(env)

  end function env_air_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function env_air_molar_den(env) ! (mole m^{-3})

    use pmc_constants

    type(env_t), intent(in) :: env      ! environment state

    env_air_molar_den = env%pressure / (const%R * env%temp)

  end function env_air_molar_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_mole_dens_to_ppb(gas_state, env)
    
    ! Convert (mole m^{-3}) to (ppb).

    use pmc_gas_state
    
    type(gas_state_t), intent(inout) :: gas_state ! gas state
    type(env_t), intent(in) :: env      ! environment state
    
    gas_state%conc = gas_state%conc / env_air_molar_den(env) * 1d9
    
  end subroutine gas_state_mole_dens_to_ppb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_gas_state(env, delta_t, old_height, &
       gas_data, gas_state)

    ! Do emissions and background dilution from the environment.

    use pmc_gas_data
    use pmc_gas_state

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    real*8, intent(in) :: old_height    ! previous height (m)
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(gas_state_t), intent(inout) :: gas_state ! gas state to update

    real*8 :: effective_dilution_rate
    type(gas_state_t) :: emission, dilution

    call gas_state_alloc(emission, gas_data%n_spec)
    call gas_state_alloc(dilution, gas_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env%gas_dilution_rate
    if (env%height > old_height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env%height - old_height) / delta_t / old_height
    end if

    ! emission = delta_t * gas_emission_rate * gas_emissions
    ! but emissions are in (mole m^{-2} s^{-1})
    call gas_state_copy(env%gas_emissions, emission)
    call gas_state_scale(emission, 1d0 / env%height)
    call gas_state_mole_dens_to_ppb(emission, env)
    call gas_state_scale(emission, delta_t * env%gas_emission_rate)

    ! dilution = delta_t * gas_dilution_rate * (gas_background - gas_state)
    call gas_state_copy(env%gas_background, dilution)
    call gas_state_sub(dilution, gas_state)
    call gas_state_scale(dilution, delta_t * effective_dilution_rate)

    call gas_state_add(gas_state, emission)
    call gas_state_add(gas_state, dilution)

    call gas_state_free(emission)
    call gas_state_free(dilution)

  end subroutine env_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_aero_state(env, delta_t, old_height, bin_grid, &
       aero_data, aero_state, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! particle aerosol distribution.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_aero_binned

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    real*8, intent(in) :: old_height    ! previous height (m)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_state_t), intent(inout) :: aero_state ! aero state to update
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    integer :: i
    real*8 :: sample_vol, sample_prop, effective_dilution_rate
    type(aero_state_t) :: aero_state_delta
    type(aero_binned_t) :: aero_binned_delta

    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state_delta)
    call aero_binned_alloc(aero_binned_delta, bin_grid%n_bin, aero_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env%aero_dilution_rate
    if (env%height > old_height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env%height - old_height) / delta_t / old_height
    end if

    ! loss to background
    sample_prop = delta_t * effective_dilution_rate
    if (sample_prop > 1d0) then
       write(0,*) 'ERROR: effective dilution rate is too high for this timestep'
       call exit(1)
    end if
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_state_sample(aero_state, aero_state_delta, sample_prop)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_binned_sub(aero_binned, aero_binned_delta)

    ! addition from background
    sample_vol = delta_t * effective_dilution_rate * aero_state%comp_vol
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_dist_sample(bin_grid, aero_data, env%aero_background, &
         sample_vol, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add_particles(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)
    
    ! emissions
    sample_vol = delta_t * env%aero_emission_rate &
         * aero_state%comp_vol / env%height
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_dist_sample(bin_grid, aero_data, env%aero_emissions, &
         sample_vol, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add_particles(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)

    call aero_state_free(aero_state_delta)
    call aero_binned_free(aero_binned_delta)

  end subroutine env_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_aero_binned(env, delta_t, old_height, bin_grid, &
       aero_data, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! binned aerosol distribution.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_binned

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    real*8, intent(in) :: old_height    ! previous height (m)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    type(aero_binned_t) :: emission, dilution
    real*8 :: effective_dilution_rate

    call aero_binned_alloc(emission, bin_grid%n_bin, aero_data%n_spec)
    call aero_binned_alloc(dilution, bin_grid%n_bin, aero_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env%aero_dilution_rate
    if (env%height > old_height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env%height - old_height) / delta_t / old_height
    end if

    ! emission = delta_t * aero_emission_rate * aero_emissions
    ! but emissions are #/m^2 so we need to divide by height
    call aero_binned_add_aero_dist(emission, bin_grid, env%aero_emissions)
    call aero_binned_scale(emission, &
         delta_t * env%aero_emission_rate / env%height)

    ! dilution = delta_t * aero_dilution_rate * (aero_background - aero_binned)
    call aero_binned_add_aero_dist(dilution, bin_grid, env%aero_background)
    call aero_binned_sub(dilution, aero_binned)
    call aero_binned_scale(dilution, delta_t * effective_dilution_rate)

    call aero_binned_add(aero_binned, emission)
    call aero_binned_add(aero_binned, dilution)

    call aero_binned_free(emission)
    call aero_binned_free(dilution)

  end subroutine env_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_env(file, env)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(env_t), intent(in) :: env      ! environment to write
    
    call inout_write_comment(file, "begin env")
    call inout_write_real(file, "temp(K)", env%temp)
    call inout_write_real(file, "rel_humidity(1)", env%rel_humid)
    call inout_write_real(file, "pressure(Pa)", env%pressure)
    call inout_write_real(file, "longitude(deg)", env%longitude)
    call inout_write_real(file, "latitude(deg)", env%latitude)
    call inout_write_real(file, "altitude(m)", env%altitude)
    call inout_write_real(file, "start_time(s)", env%start_time)
    call inout_write_integer(file, "start_day(days)", env%start_day)
    call inout_write_real(file, "height(m)", env%height)
    call inout_write_gas_state(file, env%gas_emissions)
    call inout_write_real(file, "gas_emit_rate(1/s)", env%gas_emission_rate)
    call inout_write_gas_state(file, env%gas_background)
    call inout_write_real(file, "gas_dilute_rate(1/s)", env%gas_dilution_rate)
    call inout_write_aero_dist(file, env%aero_emissions)
    call inout_write_real(file, "aero_emit_rate(1/s)", env%aero_emission_rate)
    call inout_write_aero_dist(file, env%aero_background)
    call inout_write_real(file, "aero_dilute_rat(1/s)", env%aero_dilution_rate)
    call inout_write_comment(file, "end env")

  end subroutine inout_write_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_env(file, env)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(env_t), intent(out) :: env      ! environment to read
    
    call inout_check_comment(file, "begin env")
    call inout_read_real(file, "temp(K)", env%temp)
    call inout_read_real(file, "rel_humidity(1)", env%rel_humid)
    call inout_read_real(file, "pressure(Pa)", env%pressure)
    call inout_read_real(file, "longitude(deg)", env%longitude)
    call inout_read_real(file, "latitude(deg)", env%latitude)
    call inout_read_real(file, "altitude(m)", env%altitude)
    call inout_read_real(file, "start_time(s)", env%start_time)
    call inout_read_integer(file, "start_day(days)", env%start_day)
    call inout_read_real(file, "height(m)", env%height)
    call inout_read_gas_state(file, env%gas_emissions)
    call inout_read_real(file, "gas_emit_rate(1/s)", env%gas_emission_rate)
    call inout_read_gas_state(file, env%gas_background)
    call inout_read_real(file, "gas_dilute_rate(1/s)", env%gas_dilution_rate)
    call inout_read_aero_dist(file, env%aero_emissions)
    call inout_read_real(file, "aero_emit_rate(1/s)", env%aero_emission_rate)
    call inout_read_aero_dist(file, env%aero_background)
    call inout_read_real(file, "aero_dilute_rat(1/s)", env%aero_dilution_rate)
    call inout_check_comment(file, "end env")

  end subroutine inout_read_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_env(file, env)

    ! Read environment specification from a inout file.

    use pmc_inout

    type(inout_file_t), intent(inout) :: file ! inout file
    type(env_t), intent(out) :: env     ! environment data

    call env_alloc(env)
    call inout_read_real(file, 'rel_humidity', env%rel_humid)
    call inout_read_real(file, 'pressure', env%pressure)
    call inout_read_real(file, 'latitude', env%latitude)
    call inout_read_real(file, 'longitude', env%longitude)
    call inout_read_real(file, 'altitude', env%altitude)
    call inout_read_real(file, 'start_time', env%start_time)
    call inout_read_integer(file, 'start_day', env%start_day)

  end subroutine spec_read_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_average(env_vec, env_avg)
    
    ! Computes the average of an array of env.

    use pmc_util
    use pmc_gas_state
    use pmc_aero_dist
    
    type(env_t), intent(in) :: env_vec(:) ! array of env
    type(env_t), intent(out) :: env_avg   ! average of env_vec

    call average_real(env_vec%temp, env_avg%temp)
    call average_real(env_vec%rel_humid, env_avg%rel_humid)
    call average_real(env_vec%pressure, env_avg%pressure)
    call average_real(env_vec%longitude, env_avg%longitude)
    call average_real(env_vec%latitude, env_avg%latitude)
    call average_real(env_vec%altitude, env_avg%altitude)
    call average_real(env_vec%start_time, env_avg%start_time)
    call average_integer(env_vec%start_day, env_avg%start_day)
    call average_real(env_vec%height, env_avg%height)
    call gas_state_average(env_vec%gas_emissions, env_avg%gas_emissions)
    call average_real(env_vec%gas_emission_rate, env_avg%gas_emission_rate)
    call gas_state_average(env_vec%gas_background, env_avg%gas_background)
    call average_real(env_vec%gas_dilution_rate, env_avg%gas_dilution_rate)
    call aero_dist_average(env_vec%aero_emissions, env_avg%aero_emissions)
    call average_real(env_vec%aero_emission_rate, env_avg%aero_emission_rate)
    call aero_dist_average(env_vec%aero_background, env_avg%aero_background)
    call average_real(env_vec%aero_dilution_rate, env_avg%aero_dilution_rate)
    
  end subroutine env_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_mix(val)

    ! Average val over all processes.

    use pmc_mpi
    
    type(env_t), intent(inout) :: val ! value to average

#ifdef PMC_USE_MPI
    type(env_t) :: val_avg

    call env_alloc(val_avg)
    call pmc_mpi_allreduce_average_real(val%temp, val_avg%temp)
    call pmc_mpi_allreduce_average_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_allreduce_average_real(val%pressure, val_avg%pressure)
    val%temp = val_avg%temp
    val%rel_humid = val_avg%rel_humid
    val%pressure = val_avg%pressure
    call env_free(val_avg)
#endif

  end subroutine env_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_env(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(env_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_size_env = &
         pmc_mpi_pack_size_real(val%temp) &
         + pmc_mpi_pack_size_real(val%rel_humid) &
         + pmc_mpi_pack_size_real(val%pressure) &
         + pmc_mpi_pack_size_real(val%longitude) &
         + pmc_mpi_pack_size_real(val%latitude) &
         + pmc_mpi_pack_size_real(val%altitude) &
         + pmc_mpi_pack_size_real(val%start_time) &
         + pmc_mpi_pack_size_integer(val%start_day) &
         + pmc_mpi_pack_size_real(val%height) &
         + pmc_mpi_pack_size_gas_state(val%gas_emissions) &
         + pmc_mpi_pack_size_real(val%gas_emission_rate) &
         + pmc_mpi_pack_size_gas_state(val%gas_background) &
         + pmc_mpi_pack_size_real(val%gas_dilution_rate) &
         + pmc_mpi_pack_size_aero_dist(val%aero_emissions) &
         + pmc_mpi_pack_size_real(val%aero_emission_rate) &
         + pmc_mpi_pack_size_aero_dist(val%aero_background) &
         + pmc_mpi_pack_size_real(val%aero_dilution_rate)

  end function pmc_mpi_pack_size_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_env(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(env_t), intent(in) :: val ! value to pack

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
    call pmc_mpi_pack_real(buffer, position, val%height)
    call pmc_mpi_pack_gas_state(buffer, position, val%gas_emissions)
    call pmc_mpi_pack_real(buffer, position, val%gas_emission_rate)
    call pmc_mpi_pack_gas_state(buffer, position, val%gas_background)
    call pmc_mpi_pack_real(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_pack_aero_dist(buffer, position, val%aero_emissions)
    call pmc_mpi_pack_real(buffer, position, val%aero_emission_rate)
    call pmc_mpi_pack_aero_dist(buffer, position, val%aero_background)
    call pmc_mpi_pack_real(buffer, position, val%aero_dilution_rate)
    call assert(464101191, position - prev_position == pmc_mpi_pack_size_env(val))
#endif

  end subroutine pmc_mpi_pack_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_env(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(env_t), intent(out) :: val ! value to pack

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
    call pmc_mpi_unpack_real(buffer, position, val%height)
    call pmc_mpi_unpack_gas_state(buffer, position, val%gas_emissions)
    call pmc_mpi_unpack_real(buffer, position, val%gas_emission_rate)
    call pmc_mpi_unpack_gas_state(buffer, position, val%gas_background)
    call pmc_mpi_unpack_real(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_emissions)
    call pmc_mpi_unpack_real(buffer, position, val%aero_emission_rate)
    call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_background)
    call pmc_mpi_unpack_real(buffer, position, val%aero_dilution_rate)
    call assert(205696745, position - prev_position == pmc_mpi_pack_size_env(val))
#endif

  end subroutine pmc_mpi_unpack_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_env(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_mpi

    type(env_t), intent(in) :: val ! value to average
    type(env_t), intent(out) :: val_avg ! result

    call env_alloc(val_avg)
    call env_copy(val, val_avg)
    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)

  end subroutine pmc_mpi_reduce_avg_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env
