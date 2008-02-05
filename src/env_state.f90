! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Environment parameters.
!
! The temperature profile is proscribed as a function of time by
! giving a number of times and the corresponding temperatures. Linear
! interpolation is used between the times, with constant interpolation
! outside of the range of times.

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
#ifdef PMC_USE_MPI
  use mpi
#endif
  
  type env_state_t
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
     real*8 :: aero_dilution_rate       ! aero-background dilute rate (s^{-1})
  end type env_state_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_alloc(env_state)

    ! Allocate an empty environment.

    type(env_state_t), intent(out) :: env_state   ! environment

    env_state%temp = 0d0
    env_state%rel_humid = 0d0
    env_state%pressure = 0d0
    env_state%longitude = 0d0
    env_state%latitude = 0d0
    env_state%altitude = 0d0
    env_state%start_time = 0d0
    env_state%start_day = 0
    env_state%height = 0d0

    call gas_state_alloc(env_state%gas_emissions, 0)
    call gas_state_alloc(env_state%gas_background, 0)
    env_state%gas_emission_rate = 0d0
    env_state%gas_dilution_rate = 0d0
    call aero_dist_alloc(env_state%aero_emissions, 0, 0, 0)
    call aero_dist_alloc(env_state%aero_background, 0, 0, 0)
    env_state%aero_emission_rate = 0d0
    env_state%aero_dilution_rate = 0d0

  end subroutine env_state_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_free(env_state)

    ! Free all storage.

    type(env_state_t), intent(out) :: env_state   ! environment

    call gas_state_free(env_state%gas_emissions)
    call gas_state_free(env_state%gas_background)
    call aero_dist_free(env_state%aero_emissions)
    call aero_dist_free(env_state%aero_background)

  end subroutine env_state_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_add(env_state, env_state_delta)
    
    ! env_state += env_state_delta

    type(env_state_t), intent(inout) :: env_state   ! environment
    type(env_state_t), intent(in) :: env_state_delta ! increment

    env_state%temp = env_state%temp + env_state_delta%temp
    env_state%rel_humid = env_state%rel_humid + env_state_delta%rel_humid
    env_state%pressure = env_state%pressure + env_state_delta%pressure
    env_state%longitude = env_state%longitude + env_state_delta%longitude
    env_state%latitude = env_state%latitude + env_state_delta%latitude
    env_state%altitude = env_state%altitude + env_state_delta%altitude
    env_state%start_time = env_state%start_time + env_state_delta%start_time
    env_state%start_day = env_state%start_day + env_state_delta%start_day
    env_state%height = env_state%height + env_state_delta%height
    call gas_state_add(env_state%gas_emissions, env_state_delta%gas_emissions)
    env_state%gas_emission_rate = env_state%gas_emission_rate &
         + env_state_delta%gas_emission_rate
    call gas_state_add(env_state%gas_background, &
         env_state_delta%gas_background)
    env_state%gas_dilution_rate = env_state%gas_dilution_rate &
         + env_state_delta%gas_dilution_rate
    call aero_dist_add(env_state%aero_emissions, &
         env_state_delta%aero_emissions)
    env_state%aero_emission_rate = env_state%aero_emission_rate &
         + env_state_delta%aero_emission_rate
    call aero_dist_add(env_state%aero_background, &
         env_state_delta%aero_background)
    env_state%aero_dilution_rate = env_state%aero_dilution_rate &
         + env_state_delta%aero_dilution_rate
    
  end subroutine env_state_add
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_scale(env_state, alpha)
    
    ! env_state *= alpha

    type(env_state_t), intent(inout) :: env_state   ! environment
    real*8, intent(in) :: alpha         ! scale factor

    env_state%temp = env_state%temp * alpha
    env_state%rel_humid = env_state%rel_humid * alpha
    env_state%pressure = env_state%pressure * alpha
    env_state%longitude = env_state%longitude * alpha
    env_state%latitude = env_state%latitude * alpha
    env_state%altitude = env_state%altitude * alpha
    env_state%start_time = env_state%start_time * alpha
    env_state%start_day = nint(dble(env_state%start_day) * alpha)
    env_state%height = env_state%height * alpha
    call gas_state_scale(env_state%gas_emissions, alpha)
    env_state%gas_emission_rate = env_state%gas_emission_rate * alpha
    call gas_state_scale(env_state%gas_background, alpha)
    env_state%gas_dilution_rate = env_state%gas_dilution_rate * alpha
    call aero_dist_scale(env_state%aero_emissions, alpha)
    env_state%aero_emission_rate = env_state%aero_emission_rate * alpha
    call aero_dist_scale(env_state%aero_background, alpha)
    env_state%aero_dilution_rate = env_state%aero_dilution_rate * alpha
    
  end subroutine env_state_scale
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_copy(env_from, env_to)
    
    ! env_to = env_from

    type(env_state_t), intent(in) :: env_from ! original
    type(env_state_t), intent(inout) :: env_to ! destination

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
    
  end subroutine env_state_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine env_state_change_water_volume(env_state, aero_data, dv)
    
    ! Adds the given water volume to the water vapor and updates all
    ! environment quantities.
    
    type(env_state_t), intent(inout) :: env_state ! environment state to update
    type(aero_data_t), intent(in) :: aero_data ! aero_data constants
    real*8, intent(in) :: dv            ! conc of water added (m^3/m^3)
    
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
  
  real*8 function env_state_sat_vapor_pressure(env_state) ! (Pa)

    ! Computes the current saturation vapor pressure.
    
    type(env_state_t), intent(in) :: env_state      ! environment state
    
    env_state_sat_vapor_pressure = const%water_eq_vap_press &
         * 10d0**(7.45d0 * (env_state%temp - const%water_freeze_temp) &
         / (env_state%temp - 38d0))
    
  end function env_state_sat_vapor_pressure
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_kappa_rh(aero_particle, aero_data, &
       env_state) ! (1)

    ! Returns the critical relative humidity from the kappa value.

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(env_state_t), intent(in) :: env_state      ! environment state

    real*8 :: kappa, diam, C, A
    
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    A = 4d0 * const%water_surf_eng * const%water_molec_weight &
         / (const%univ_gas_const * env_state%temp * const%water_density)
    C = sqrt(4d0 * A**3 / 27d0)
    diam = vol2diam(aero_particle_volume(aero_particle))
    aero_particle_kappa_rh = C / sqrt(kappa * diam**3) + 1d0

  end function aero_particle_kappa_rh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function env_state_air_den(env_state) ! (kg m^{-3})

    type(env_state_t), intent(in) :: env_state      ! environment state

    env_state_air_den = const%air_molec_weight &
         * env_state_air_molar_den(env_state)

  end function env_state_air_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function env_state_air_molar_den(env_state) ! (mole m^{-3})

    type(env_state_t), intent(in) :: env_state      ! environment state

    env_state_air_molar_den = env_state%pressure &
         / (const%univ_gas_const * env_state%temp)

  end function env_state_air_molar_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_state_mole_dens_to_ppb(gas_state, env_state)
    
    ! Convert (mole m^{-3}) to (ppb).

    type(gas_state_t), intent(inout) :: gas_state ! gas state
    type(env_state_t), intent(in) :: env_state      ! environment state
    
    gas_state%conc = gas_state%conc / env_state_air_molar_den(env_state) * 1d9
    
  end subroutine gas_state_mole_dens_to_ppb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_update_gas_state(env_state, delta_t, old_height, &
       gas_data, gas_state)

    ! Do emissions and background dilution from the environment.

    type(env_state_t), intent(in) :: env_state      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    real*8, intent(in) :: old_height    ! previous height (m)
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(gas_state_t), intent(inout) :: gas_state ! gas state to update

    real*8 :: effective_dilution_rate
    type(gas_state_t) :: emission, dilution

    call gas_state_alloc(emission, gas_data%n_spec)
    call gas_state_alloc(dilution, gas_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env_state%gas_dilution_rate
    if (env_state%height > old_height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env_state%height - old_height) / delta_t / old_height
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

    call gas_state_free(emission)
    call gas_state_free(dilution)

  end subroutine env_state_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_update_aero_state(env_state, delta_t, old_height, &
       bin_grid, aero_data, aero_state, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! particle aerosol distribution.

    type(env_state_t), intent(in) :: env_state      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    real*8, intent(in) :: old_height    ! previous height (m)
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_state_t), intent(inout) :: aero_state ! aero state to update
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    integer :: i
    real*8 :: sample_prop, effective_dilution_rate
    type(aero_state_t) :: aero_state_delta
    type(aero_binned_t) :: aero_binned_delta

    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state_delta)
    call aero_binned_alloc(aero_binned_delta, bin_grid%n_bin, aero_data%n_spec)

    ! account for height changes
    effective_dilution_rate = env_state%aero_dilution_rate
    if (env_state%height > old_height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env_state%height - old_height) / delta_t / old_height
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
    call aero_dist_sample(bin_grid, aero_data, env_state%aero_background, &
         sample_prop, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add_particles(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)
    
    ! emissions
    sample_prop = delta_t * env_state%aero_emission_rate / env_state%height
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_dist_sample(bin_grid, aero_data, env_state%aero_emissions, &
         sample_prop, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add_particles(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)

    call aero_state_free(aero_state_delta)
    call aero_binned_free(aero_binned_delta)

  end subroutine env_state_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_update_aero_binned(env_state, delta_t, old_height, &
       bin_grid, aero_data, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! binned aerosol distribution.

    type(env_state_t), intent(in) :: env_state      ! current environment
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
    effective_dilution_rate = env_state%aero_dilution_rate
    if (env_state%height > old_height) then
       effective_dilution_rate = effective_dilution_rate &
            + (env_state%height - old_height) / delta_t / old_height
    end if

    ! emission = delta_t * aero_emission_rate * aero_emissions
    ! but emissions are #/m^2 so we need to divide by height
    call aero_binned_add_aero_dist(emission, bin_grid, &
         env_state%aero_emissions)
    call aero_binned_scale(emission, &
         delta_t * env_state%aero_emission_rate / env_state%height)

    ! dilution = delta_t * aero_dilution_rate * (aero_background - aero_binned)
    call aero_binned_add_aero_dist(dilution, bin_grid, &
         env_state%aero_background)
    call aero_binned_sub(dilution, aero_binned)
    call aero_binned_scale(dilution, delta_t * effective_dilution_rate)

    call aero_binned_add(aero_binned, emission)
    call aero_binned_add(aero_binned, dilution)

    call aero_binned_free(emission)
    call aero_binned_free(dilution)

  end subroutine env_state_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_env_state(file, env_state)
    
    ! Write full state.
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(env_state_t), intent(in) :: env_state      ! environment to write
    
    call inout_write_comment(file, "begin env_state")
    call inout_write_real(file, "temp(K)", env_state%temp)
    call inout_write_real(file, "rel_humidity(1)", env_state%rel_humid)
    call inout_write_real(file, "pressure(Pa)", env_state%pressure)
    call inout_write_real(file, "longitude(deg)", env_state%longitude)
    call inout_write_real(file, "latitude(deg)", env_state%latitude)
    call inout_write_real(file, "altitude(m)", env_state%altitude)
    call inout_write_real(file, "start_time(s)", env_state%start_time)
    call inout_write_integer(file, "start_day(days)", env_state%start_day)
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

  subroutine inout_read_env_state(file, env_state)
    
    ! Read full state.
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(env_state_t), intent(out) :: env_state      ! environment to read
    
    call inout_check_comment(file, "begin env_state")
    call inout_read_real(file, "temp(K)", env_state%temp)
    call inout_read_real(file, "rel_humidity(1)", env_state%rel_humid)
    call inout_read_real(file, "pressure(Pa)", env_state%pressure)
    call inout_read_real(file, "longitude(deg)", env_state%longitude)
    call inout_read_real(file, "latitude(deg)", env_state%latitude)
    call inout_read_real(file, "altitude(m)", env_state%altitude)
    call inout_read_real(file, "start_time(s)", env_state%start_time)
    call inout_read_integer(file, "start_day(days)", env_state%start_day)
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

  subroutine spec_read_env_state(file, env_state)

    ! Read environment specification from a inout file.

    type(inout_file_t), intent(inout) :: file ! inout file
    type(env_state_t), intent(out) :: env_state     ! environment data

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

  subroutine env_state_average(env_vec, env_avg)
    
    ! Computes the average of an array of env_state.

    type(env_state_t), intent(in) :: env_vec(:) ! array of env_state
    type(env_state_t), intent(out) :: env_avg   ! average of env_vec

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
    
  end subroutine env_state_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_state_mix(val)

    ! Average val over all processes.

    type(env_state_t), intent(inout) :: val ! value to average

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

  integer function pmc_mpi_pack_size_env_state(val)

    ! Determines the number of bytes required to pack the given value.

    type(env_state_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_size_env_state = &
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

  end function pmc_mpi_pack_size_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_env_state(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(env_state_t), intent(in) :: val ! value to pack

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
    call assert(464101191, &
         position - prev_position == pmc_mpi_pack_size_env_state(val))
#endif

  end subroutine pmc_mpi_pack_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_env_state(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(env_state_t), intent(out) :: val ! value to pack

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
    call assert(205696745, &
         position - prev_position == pmc_mpi_pack_size_env_state(val))
#endif

  end subroutine pmc_mpi_unpack_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_env_state(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    type(env_state_t), intent(in) :: val ! value to average
    type(env_state_t), intent(out) :: val_avg ! result

    call env_state_alloc(val_avg)
    call env_state_copy(val, val_avg)
    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)

  end subroutine pmc_mpi_reduce_avg_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env_state
