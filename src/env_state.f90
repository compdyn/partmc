! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
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
  use pmc_aero_binned
  use pmc_util
  use pmc_gas_data
  use pmc_bin_grid
  use pmc_aero_state
  use pmc_spec_file
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
     real(kind=dp) :: temp
     !> Relative humidity (1).
     real(kind=dp) :: rel_humid
     !> Ambient pressure (Pa).
     real(kind=dp) :: pressure
     !> Longitude (degrees).
     real(kind=dp) :: longitude
     !> Latitude (degrees).
     real(kind=dp) :: latitude
     !> Altitude (m).
     real(kind=dp) :: altitude
     !> Start time (s since 00:00 UTC on \c start_day).
     real(kind=dp) :: start_time
     !> Start day of year (UTC).
     integer :: start_day
     !> Time since \c start_time (s).
     real(kind=dp) :: elapsed_time
     !> Box height (m).
     real(kind=dp) :: height
     !> Gas emissions.
     type(gas_state_t) :: gas_emissions
     !> Gas emisssion rate (s^{-1}).
     real(kind=dp) :: gas_emission_rate
     !> Background gas mixing ratios.
     type(gas_state_t) :: gas_background
     !> Gas-background dilution rate (s^{-1}).
     real(kind=dp) :: gas_dilution_rate
     !> Aerosol emissions.
     type(aero_dist_t) :: aero_emissions
     !> Aerosol emisssion rate (s^{-1}).
     real(kind=dp) :: aero_emission_rate
     !> Aerosol background.
     type(aero_dist_t) :: aero_background
     !> Aero-background dilute rate (s^{-1}).
     real(kind=dp) :: aero_dilution_rate
  end type env_state_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an empty environment.
  subroutine env_state_allocate(env_state)

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

    call gas_state_allocate(env_state%gas_emissions)
    call gas_state_allocate(env_state%gas_background)
    env_state%gas_emission_rate = 0d0
    env_state%gas_dilution_rate = 0d0
    call aero_dist_allocate(env_state%aero_emissions)
    call aero_dist_allocate(env_state%aero_background)
    env_state%aero_emission_rate = 0d0
    env_state%aero_dilution_rate = 0d0

  end subroutine env_state_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine env_state_deallocate(env_state)

    !> Environment.
    type(env_state_t), intent(out) :: env_state

    call gas_state_deallocate(env_state%gas_emissions)
    call gas_state_deallocate(env_state%gas_background)
    call aero_dist_deallocate(env_state%aero_emissions)
    call aero_dist_deallocate(env_state%aero_background)

  end subroutine env_state_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> env_state *= alpha
  subroutine env_state_scale(env_state, alpha)

    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    env_state%temp = env_state%temp * alpha
    env_state%rel_humid = env_state%rel_humid * alpha
    env_state%pressure = env_state%pressure * alpha
    env_state%longitude = env_state%longitude * alpha
    env_state%latitude = env_state%latitude * alpha
    env_state%altitude = env_state%altitude * alpha
    env_state%start_time = env_state%start_time * alpha
    env_state%start_day = nint(real(env_state%start_day, kind=dp) * alpha)
    env_state%elapsed_time = env_state%elapsed_time * alpha
    env_state%height = env_state%height * alpha
    call gas_state_scale(env_state%gas_emissions, alpha)
    env_state%gas_emission_rate = env_state%gas_emission_rate * alpha
    call gas_state_scale(env_state%gas_background, alpha)
    env_state%gas_dilution_rate = env_state%gas_dilution_rate * alpha
    
  end subroutine env_state_scale
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds the given water volume to the water vapor and updates all
  !> environment quantities.
  subroutine env_state_change_water_volume(env_state, aero_data, dv)
    
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Aero_data constants.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration of water added (m^3/m^3).
    real(kind=dp), intent(in) :: dv
    
    real(kind=dp) pmv     ! ambient water vapor pressure (Pa)
    real(kind=dp) mv      ! ambient water vapor density (kg m^{-3})
                   ! pmv and mv are related by the factor molec_weight/(R*T)
    real(kind=dp) dmv     ! change of water density (kg m^{-3})
    
    dmv = dv * aero_data%density(aero_data%i_water)
    pmv = env_state_sat_vapor_pressure(env_state) * env_state%rel_humid
    mv = aero_data%molec_weight(aero_data%i_water) &
         / (const%univ_gas_const*env_state%temp) * pmv
    mv = mv - dmv    
    env_state%rel_humid = const%univ_gas_const * env_state%temp &
         / aero_data%molec_weight(aero_data%i_water) * mv &
         / env_state_sat_vapor_pressure(env_state)
    
  end subroutine env_state_change_water_volume
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the current saturation vapor pressure (Pa).
  real(kind=dp) function env_state_sat_vapor_pressure(env_state)
    
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    
    env_state_sat_vapor_pressure = const%water_eq_vap_press &
         * 10d0**(7.45d0 * (env_state%temp - const%water_freeze_temp) &
         / (env_state%temp - 38d0))
    
  end function env_state_sat_vapor_pressure
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the critical relative humidity from the kappa value (1).
  real(kind=dp) function aero_particle_kappa_rh(aero_particle, aero_data, &
       env_state)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: kappa, diam, C, A
    
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    A = 4d0 * const%water_surf_eng * const%water_molec_weight &
         / (const%univ_gas_const * env_state%temp * const%water_density)
    C = sqrt(4d0 * A**3 / 27d0)
    diam = vol2diam(aero_particle_volume(aero_particle))
    aero_particle_kappa_rh = C / sqrt(kappa * diam**3) + 1d0

  end function aero_particle_kappa_rh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Air density (kg m^{-3}).
  real(kind=dp) function env_state_air_den(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_air_den = const%air_molec_weight &
         * env_state_air_molar_den(env_state)

  end function env_state_air_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Air molar density (mole m^{-3}).
  real(kind=dp) function env_state_air_molar_den(env_state)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    env_state_air_molar_den = env_state%pressure &
         / (const%univ_gas_const * env_state%temp)

  end function env_state_air_molar_den

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert (mole m^{-3}) to (ppb).
  subroutine gas_state_mole_dens_to_ppb(gas_state, env_state)

    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    
    gas_state%mix_rat = gas_state%mix_rat &
         / env_state_air_molar_den(env_state) * 1d9
    
  end subroutine gas_state_mole_dens_to_ppb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment.
  subroutine env_state_update_gas_state(env_state, delta_t, &
       old_env_state, &
       gas_data, gas_state)

    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state to update.
    type(gas_state_t), intent(inout) :: gas_state

    real(kind=dp) :: effective_dilution_rate
    type(gas_state_t) :: emission, dilution

    call gas_state_allocate_size(emission, gas_data%n_spec)
    call gas_state_allocate_size(dilution, gas_data%n_spec)

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

    call gas_state_deallocate(emission)
    call gas_state_deallocate(dilution)

  end subroutine env_state_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment for a
  !> particle aerosol distribution.
  subroutine env_state_update_aero_state(env_state, delta_t, &
       old_env_state, bin_grid, aero_data, aero_state)

    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero state to update.
    type(aero_state_t), intent(inout) :: aero_state

    integer :: i
    real(kind=dp) :: sample_prop, effective_dilution_rate
    type(aero_state_t) :: aero_state_delta

    call aero_state_allocate_size(aero_state_delta, bin_grid%n_bin, &
         aero_data%n_spec)

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
       call die_msg(925788271, &
            'effective dilution rate too high for this timestep')
    end if
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_state_sample(aero_state, aero_state_delta, sample_prop, &
         AERO_INFO_DILUTION)

    ! addition from background
    sample_prop = delta_t * effective_dilution_rate
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_state_add_aero_dist_sample(aero_state_delta, bin_grid, &
         aero_data, env_state%aero_background, sample_prop, &
         env_state%elapsed_time)
    !>DEBUG
    !write(*,*) 'env_state_update_aero_state: calling aero_state_add_particles'
    !<DEBUG
    call aero_state_add_particles(aero_state, aero_state_delta)
    
    ! emissions
    sample_prop = delta_t * env_state%aero_emission_rate / env_state%height
    call aero_state_zero(aero_state_delta)
    aero_state_delta%comp_vol = aero_state%comp_vol
    call aero_state_add_aero_dist_sample(aero_state_delta, bin_grid, &
         aero_data, env_state%aero_emissions, sample_prop, &
         env_state%elapsed_time)
    !>DEBUG
    !write(*,*) 'env_state_update_aero_state: calling aero_state_add_particles'
    !<DEBUG
    call aero_state_add_particles(aero_state, aero_state_delta)

    ! update computational volume
    aero_state%comp_vol = aero_state%comp_vol * env_state%temp &
         / old_env_state%temp

    call aero_state_deallocate(aero_state_delta)

  end subroutine env_state_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment for a
  !> binned aerosol distribution.
  subroutine env_state_update_aero_binned(env_state, delta_t, & 
       old_env_state, &
       bin_grid, aero_data, aero_binned)

    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero binned to update.
    type(aero_binned_t), intent(inout) :: aero_binned

    type(aero_binned_t) :: emission, dilution
    real(kind=dp) :: effective_dilution_rate

    call aero_binned_allocate_size(emission, bin_grid%n_bin, aero_data%n_spec)
    call aero_binned_allocate_size(dilution, bin_grid%n_bin, aero_data%n_spec)

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

    ! dilution = delta_t * aero_dilution_rate
    !            * (aero_background - aero_binned)
    call aero_binned_add_aero_dist(dilution, bin_grid, aero_data, &
         env_state%aero_background)
    call aero_binned_sub(dilution, aero_binned)
    call aero_binned_scale(dilution, delta_t * effective_dilution_rate)

    call aero_binned_add(aero_binned, emission)
    call aero_binned_add(aero_binned, dilution)

    call aero_binned_deallocate(emission)
    call aero_binned_deallocate(dilution)

  end subroutine env_state_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read environment specification from a spec file.
  subroutine spec_file_read_env_state(file, env_state)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Environment data.
    type(env_state_t), intent(out) :: env_state

    call spec_file_read_real(file, 'rel_humidity', env_state%rel_humid)
    call spec_file_read_real(file, 'pressure', env_state%pressure)
    call spec_file_read_real(file, 'latitude', env_state%latitude)
    call spec_file_read_real(file, 'longitude', env_state%longitude)
    call spec_file_read_real(file, 'altitude', env_state%altitude)
    call spec_file_read_real(file, 'start_time', env_state%start_time)
    call spec_file_read_integer(file, 'start_day', env_state%start_day)

  end subroutine spec_file_read_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes.
  subroutine env_state_mix(val)

    !> Value to average.
    type(env_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(env_state_t) :: val_avg

    call env_state_allocate(val_avg)
    call pmc_mpi_allreduce_average_real(val%temp, val_avg%temp)
    call pmc_mpi_allreduce_average_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_allreduce_average_real(val%pressure, val_avg%pressure)
    val%temp = val_avg%temp
    val%rel_humid = val_avg%rel_humid
    val%pressure = val_avg%pressure
    call env_state_deallocate(val_avg)
#endif

  end subroutine env_state_mix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average val over all processes, with the result only on the root
  !> processor.
  subroutine env_state_reduce_avg(val)

    !> Value to average.
    type(env_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    type(env_state_t) :: val_avg

    call env_state_allocate(val_avg)
    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)
    if (pmc_mpi_rank() == 0) then
       val%temp = val_avg%temp
       val%rel_humid = val_avg%rel_humid
       val%pressure = val_avg%pressure
    end if
    call env_state_deallocate(val_avg)
#endif

  end subroutine env_state_reduce_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_env_state(val, val_avg)

    !> Value to average.
    type(env_state_t), intent(in) :: val
    !> Result.
    type(env_state_t), intent(out) :: val_avg

    call env_state_allocate(val_avg)
    call env_state_copy(val, val_avg)
    call pmc_mpi_reduce_avg_real(val%temp, val_avg%temp)
    call pmc_mpi_reduce_avg_real(val%rel_humid, val_avg%rel_humid)
    call pmc_mpi_reduce_avg_real(val%pressure, val_avg%pressure)

  end subroutine pmc_mpi_reduce_avg_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine env_state_output_netcdf(env_state, ncid)
    
    !> Environment state to write.
    type(env_state_t), intent(in) :: env_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_write_real(ncid, env_state%temp, "temperature", "K")
    call pmc_nc_write_real(ncid, env_state%rel_humid, &
         "relative_humidity", "1")
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine env_state_input_netcdf(env_state, ncid)
    
    !> Environment state to read.
    type(env_state_t), intent(inout) :: env_state
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    character(len=1000) :: unit

    call pmc_nc_read_real(ncid, env_state%temp, "temperature", unit)
    call pmc_nc_read_real(ncid, env_state%rel_humid, "relative_humidity", &
         unit)
    call pmc_nc_read_real(ncid, env_state%pressure, "pressure", unit)
    call pmc_nc_read_real(ncid, env_state%longitude, "longitude", unit)
    call pmc_nc_read_real(ncid, env_state%latitude, "latitude", unit)
    call pmc_nc_read_real(ncid, env_state%altitude, "altitude", unit)
    call pmc_nc_read_real(ncid, env_state%start_time, &
         "start_time_of_day", unit)
    call pmc_nc_read_integer(ncid, env_state%start_day, &
         "start_day_of_year", unit)
    call pmc_nc_read_real(ncid, env_state%elapsed_time, "elapsed_time", unit)
    call pmc_nc_read_real(ncid, env_state%height, "height", unit)

  end subroutine env_state_input_netcdf
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env_state
