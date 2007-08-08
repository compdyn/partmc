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
     real*8 :: air_den                  ! air density (kg m^{-3})
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
    env%air_den = 0d0
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
  
  real*8 function env_sat_vapor_pressure(env) ! Pa

    ! Computes the current saturation vapor pressure.
    
    use pmc_constants
    
    type(env_t), intent(in) :: env      ! environment state
    
    env_sat_vapor_pressure = const%p00 * 10d0**(7.45d0 * (env%temp - const%T0) &
         / (env%temp - 38d0)) ! Pa
    
  end function env_sat_vapor_pressure
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_gas_state(env, delta_t, gas_data, gas_state)

    ! Do emissions and background dilution from the environment.

    use pmc_gas_data
    use pmc_gas_state

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(gas_state_t), intent(inout) :: gas_state ! gas state to update

    type(gas_state_t) :: emission, dilution

    call gas_state_alloc(emission, gas_data%n_spec)
    call gas_state_alloc(dilution, gas_data%n_spec)

    ! emission = delta_t * gas_emission_rate * gas_emissions
    call gas_state_copy(env%gas_emissions, emission)
    call gas_state_scale(emission, delta_t * env%gas_emission_rate)

    ! dilution = delta_t * gas_dilution_rate * (gas_background - gas_state)
    call gas_state_copy(env%gas_background, dilution)
    call gas_state_sub(dilution, gas_state)
    call gas_state_scale(dilution, delta_t * env%gas_dilution_rate)

    call gas_state_add(gas_state, emission)
    call gas_state_add(gas_state, dilution)

    call gas_state_free(emission)
    call gas_state_free(dilution)

  end subroutine env_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_aero_state(env, delta_t, bin_grid, &
       aero_data, aero_state, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! particle aerosol distribution.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_aero_binned

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_state_t), intent(inout) :: aero_state ! aero state to update
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    integer :: i
    real*8 :: sample_vol, sample_prop
    type(aero_state_t) :: aero_state_delta
    type(aero_binned_t) :: aero_binned_delta

    call aero_state_alloc(bin_grid%n_bin, aero_data%n_spec, aero_state_delta)
    call aero_binned_alloc(aero_binned_delta, bin_grid%n_bin, aero_data%n_spec)
    aero_state_delta%comp_vol = aero_state%comp_vol

    ! loss to background
    sample_prop = delta_t * env%aero_dilution_rate
    call aero_state_zero(aero_state_delta)
    call aero_state_sample(aero_state, aero_state_delta, sample_prop)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_binned_sub(aero_binned, aero_binned_delta)

    ! addition from background
    sample_vol = delta_t * env%aero_dilution_rate * aero_state%comp_vol
    call aero_state_zero(aero_state_delta)
    call aero_dist_sample(bin_grid, aero_data, env%aero_background, &
         sample_vol, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)
    
    ! emissions
    sample_vol = delta_t * env%aero_emission_rate * aero_state%comp_vol
    call aero_state_zero(aero_state_delta)
    call aero_dist_sample(bin_grid, aero_data, env%aero_emissions, &
         sample_vol, aero_state_delta)
    call aero_state_to_binned(bin_grid, aero_data, aero_state_delta, &
         aero_binned_delta)
    call aero_state_add(aero_state, aero_state_delta)
    call aero_binned_add(aero_binned, aero_binned_delta)

    call aero_state_free(aero_state_delta)
    call aero_binned_free(aero_binned_delta)

  end subroutine env_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_update_aero_binned(env, delta_t, bin_grid, &
       aero_data, aero_binned)

    ! Do emissions and background dilution from the environment for a
    ! binned aerosol distribution.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_binned

    type(env_t), intent(in) :: env      ! current environment
    real*8, intent(in) :: delta_t       ! time increment to update over
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aero data values
    type(aero_binned_t), intent(inout) :: aero_binned ! aero binned to update

    type(aero_binned_t) :: emission, dilution

    call aero_binned_alloc(emission, bin_grid%n_bin, aero_data%n_spec)
    call aero_binned_alloc(dilution, bin_grid%n_bin, aero_data%n_spec)

    ! emission = delta_t * gas_emission_rate * gas_emissions
    call aero_binned_add_aero_dist(emission, bin_grid, env%aero_emissions)
    call aero_binned_scale(emission, delta_t * env%aero_emission_rate)

    ! dilution = delta_t * gas_dilution_rate * (gas_background - aero_binned)
    call aero_binned_add_aero_dist(dilution, bin_grid, env%aero_background)
    call aero_binned_sub(dilution, aero_binned)
    call aero_binned_scale(dilution, delta_t * env%aero_dilution_rate)

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
    
    call inout_write_real(file, "temp(K)", env%temp)
    call inout_write_real(file, "rel_humidity(1)", env%rel_humid)
    call inout_write_real(file, "pressure(Pa)", env%pressure)
    call inout_write_real(file, "air_density(kg/m^3)", env%air_den)
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

  end subroutine inout_write_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_env(file, env)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(env_t), intent(out) :: env      ! environment to read
    
    call inout_read_real(file, "temp(K)", env%temp)
    call inout_read_real(file, "rel_humidity(1)", env%rel_humid)
    call inout_read_real(file, "pressure(Pa)", env%pressure)
    call inout_read_real(file, "air_density(kg/m^3)", env%air_den)
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
    call inout_read_real(file, 'air_density', env%air_den)
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
    call average_real(env_vec%air_den, env_avg%air_den)
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

  integer function pmc_mpi_pack_env_size(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(env_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_env_size = &
         pmc_mpi_pack_real_size(val%temp) &
         + pmc_mpi_pack_real_size(val%rel_humid) &
         + pmc_mpi_pack_real_size(val%pressure) &
         + pmc_mpi_pack_real_size(val%air_den) &
         + pmc_mpi_pack_real_size(val%longitude) &
         + pmc_mpi_pack_real_size(val%latitude) &
         + pmc_mpi_pack_real_size(val%altitude) &
         + pmc_mpi_pack_real_size(val%start_time) &
         + pmc_mpi_pack_integer_size(val%start_day) &
         + pmc_mpi_pack_real_size(val%height) &
         + pmc_mpi_pack_gas_state_size(val%gas_emissions) &
         + pmc_mpi_pack_real_size(val%gas_emission_rate) &
         + pmc_mpi_pack_gas_state_size(val%gas_background) &
         + pmc_mpi_pack_real_size(val%gas_dilution_rate) &
         + pmc_mpi_pack_aero_dist_size(val%aero_emissions) &
         + pmc_mpi_pack_real_size(val%aero_emission_rate) &
         + pmc_mpi_pack_aero_dist_size(val%aero_background) &
         + pmc_mpi_pack_real_size(val%aero_dilution_rate)

  end function pmc_mpi_pack_env_size

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
    call pmc_mpi_pack_real(buffer, position, val%air_den)
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
    call assert(position - prev_position == pmc_mpi_pack_env_size(val))
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
    call pmc_mpi_unpack_real(buffer, position, val%air_den)
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
    call assert(position - prev_position == pmc_mpi_pack_env_size(val))
#endif

  end subroutine pmc_mpi_unpack_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_env(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_mpi

    type(env_t), intent(in) :: val ! value to average
    type(env_t), intent(out) :: val_avg ! result

    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)
    call pmc_mpi_reduce_avg_real(val%air_den, val_avg%air_den)

  end subroutine pmc_mpi_reduce_avg_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env
