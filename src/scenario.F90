! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_scenario module.

!> The scenario_t structure and associated subroutines.
module pmc_scenario

  use pmc_gas_state
  use pmc_aero_dist
  use pmc_util
  use pmc_env_state
  use pmc_aero_state
  use pmc_spec_file
  use pmc_aero_data
  use pmc_gas_data
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Type code for an undefined or invalid loss function.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_INVALID  = 0
  !> Type code for a zero loss function.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_ZERO     = 1
  !> Type code for a constant loss function.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_CONSTANT = 2
  !> Type code for a loss rate fuction proportional to volume.
  integer, parameter :: SCENARIO_LOSS_FUNCTION_VOLUME   = 3
  
  !> Scenario data.
  !!
  !! This is everything needed to drive the scenario being simulated.
  !!
  !! The temperature, pressure, emissions and background states are profiles
  !! prescribed as functions of time by giving a number of times and
  !! the corresponding data. Simple data such as temperature and pressure is
  !! linearly interpolated between times, with constant interpolation
  !! outside of the range of times. Gases and aerosols are
  !! interpolated with gas_state_interp_1d() and
  !! aero_dist_interp_1d(), respectively.
  type scenario_t
     !> Temperature set-point times (s).
     real(kind=dp), pointer :: temp_time(:)
     !> Temperatures at set-points (K).
     real(kind=dp), pointer :: temp(:)

     !> Pressure set-point times (s).
     real(kind=dp), pointer :: pressure_time(:)
     !> Pressures at set-points (Pa).
     real(kind=dp), pointer :: pressure(:)

     !> Height set-point times (s).
     real(kind=dp), pointer :: height_time(:)
     !> Heights at set-points (m).
     real(kind=dp), pointer :: height(:)

     !> Gas emission set-point times (s).
     real(kind=dp), pointer :: gas_emission_time(:)
     !> Gas emisssion rate scales at set-points (1).
     real(kind=dp), pointer :: gas_emission_rate_scale(:)
     !> Gas emission rates at set-points (mol m^{-2} s^{-1}).
     type(gas_state_t), pointer :: gas_emission(:)

     !> Gas-background dilution set-point times (s).
     real(kind=dp), pointer :: gas_dilution_time(:)
     !> Gas-background dilution rates at set-points (s^{-1}).
     real(kind=dp), pointer :: gas_dilution_rate(:)
     !> Background gas mixing ratios at set-points (ppb).
     type(gas_state_t), pointer :: gas_background(:)

     !> Aerosol emission set-points times (s).
     real(kind=dp), pointer :: aero_emission_time(:)
     !> Aerosol emission rate scales at set-points (1).
     real(kind=dp), pointer :: aero_emission_rate_scale(:)
     !> Aerosol emissions at set-points (# m^{-2} s^{-1}).
     type(aero_dist_t), pointer :: aero_emission(:)

     !> Aerosol-background dilution set-point times (s).
     real(kind=dp), pointer :: aero_dilution_time(:)
     !> Aerosol-background dilution rates at set-points (s^{-1}).
     real(kind=dp), pointer :: aero_dilution_rate(:)
     !> Aerosol background at set-points (# m^{-3}).
     type(aero_dist_t), pointer :: aero_background(:)
  end type scenario_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an scenario.
  subroutine scenario_allocate(scenario)

    !> Environment data.
    type(scenario_t), intent(out) :: scenario

    allocate(scenario%temp_time(0))
    allocate(scenario%temp(0))

    allocate(scenario%pressure_time(0))
    allocate(scenario%pressure(0))

    allocate(scenario%height_time(0))
    allocate(scenario%height(0))

    allocate(scenario%gas_emission_time(0))
    allocate(scenario%gas_emission_rate_scale(0))
    allocate(scenario%gas_emission(0))

    allocate(scenario%gas_dilution_time(0))
    allocate(scenario%gas_dilution_rate(0))
    allocate(scenario%gas_background(0))

    allocate(scenario%aero_emission_time(0))
    allocate(scenario%aero_emission_rate_scale(0))
    allocate(scenario%aero_emission(0))

    allocate(scenario%aero_dilution_time(0))
    allocate(scenario%aero_dilution_rate(0))
    allocate(scenario%aero_background(0))

  end subroutine scenario_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine scenario_deallocate(scenario)

    !> Environment data.
    type(scenario_t), intent(inout) :: scenario

    integer :: i

    deallocate(scenario%temp_time)
    deallocate(scenario%temp)

    deallocate(scenario%pressure_time)
    deallocate(scenario%pressure)

    deallocate(scenario%height_time)
    deallocate(scenario%height)

    do i = 1,size(scenario%gas_emission)
       call gas_state_deallocate(scenario%gas_emission(i))
    end do
    deallocate(scenario%gas_emission_time)
    deallocate(scenario%gas_emission_rate_scale)
    deallocate(scenario%gas_emission)

    do i = 1,size(scenario%gas_background)
       call gas_state_deallocate(scenario%gas_background(i))
    end do
    deallocate(scenario%gas_dilution_time)
    deallocate(scenario%gas_dilution_rate)
    deallocate(scenario%gas_background)

    do i = 1,size(scenario%aero_emission)
       call aero_dist_deallocate(scenario%aero_emission(i))
    end do
    deallocate(scenario%aero_emission_time)
    deallocate(scenario%aero_emission_rate_scale)
    deallocate(scenario%aero_emission)

    do i = 1,size(scenario%aero_background)
       call aero_dist_deallocate(scenario%aero_background(i))
    end do
    deallocate(scenario%aero_dilution_time)
    deallocate(scenario%aero_dilution_rate)
    deallocate(scenario%aero_background)

  end subroutine scenario_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy structure.
  subroutine scenario_copy(scenario_from, scenario_to)

    !> Source environment data.
    type(scenario_t), intent(in) :: scenario_from
    !> Destination environment data.
    type(scenario_t), intent(inout) :: scenario_to

    integer :: i

    call scenario_deallocate(scenario_to)

    allocate(scenario_to%temp_time( &
         size(scenario_from%temp_time)))
    scenario_to%temp_time = scenario_from%temp_time
    allocate(scenario_to%temp( &
         size(scenario_from%temp)))
    scenario_to%temp = scenario_from%temp

    allocate(scenario_to%pressure_time( &
         size(scenario_from%pressure_time)))
    scenario_to%pressure_time = scenario_from%pressure_time
    allocate(scenario_to%pressure( &
         size(scenario_from%pressure)))
    scenario_to%pressure = scenario_from%pressure

    allocate(scenario_to%height_time( &
         size(scenario_from%height_time)))
    scenario_to%height_time = scenario_from%height_time
    allocate(scenario_to%height( &
         size(scenario_from%height)))
    scenario_to%height = scenario_from%height

    allocate(scenario_to%gas_emission_time( &
         size(scenario_from%gas_emission_time)))
    scenario_to%gas_emission_time = scenario_from%gas_emission_time
    allocate(scenario_to%gas_emission_rate_scale( &
         size(scenario_from%gas_emission_rate_scale)))
    scenario_to%gas_emission_rate_scale = scenario_from%gas_emission_rate_scale
    allocate(scenario_to%gas_emission( &
         size(scenario_from%gas_emission)))
    do i = 1,size(scenario_from%gas_emission)
       call gas_state_allocate(scenario_to%gas_emission(i))
       call gas_state_copy(scenario_from%gas_emission(i), &
            scenario_to%gas_emission(i))
    end do

    allocate(scenario_to%gas_dilution_time( &
         size(scenario_from%gas_dilution_time)))
    scenario_to%gas_dilution_time = scenario_from%gas_dilution_time
    allocate(scenario_to%gas_dilution_rate( &
         size(scenario_from%gas_dilution_rate)))
    scenario_to%gas_dilution_rate = scenario_from%gas_dilution_rate
    allocate(scenario_to%gas_background( &
         size(scenario_from%gas_background)))
    do i = 1,size(scenario_from%gas_background)
       call gas_state_allocate(scenario_to%gas_background(i))
       call gas_state_copy(scenario_from%gas_background(i), &
            scenario_to%gas_background(i))
    end do

    allocate(scenario_to%aero_emission_time( &
         size(scenario_from%aero_emission_time)))
    scenario_to%aero_emission_time = scenario_from%aero_emission_time
    allocate(scenario_to%aero_emission_rate_scale( &
         size(scenario_from%aero_emission_rate_scale)))
    scenario_to%aero_emission_rate_scale &
         = scenario_from%aero_emission_rate_scale
    allocate(scenario_to%aero_emission( &
         size(scenario_from%aero_emission)))
    do i = 1,size(scenario_from%aero_emission)
       call aero_dist_allocate(scenario_to%aero_emission(i))
       call aero_dist_copy(scenario_from%aero_emission(i), &
            scenario_to%aero_emission(i))
    end do

    allocate(scenario_to%aero_dilution_time( &
         size(scenario_from%aero_dilution_time)))
    scenario_to%aero_dilution_time = scenario_from%aero_dilution_time
    allocate(scenario_to%aero_dilution_rate( &
         size(scenario_from%aero_dilution_rate)))
    scenario_to%aero_dilution_rate = scenario_from%aero_dilution_rate
    allocate(scenario_to%aero_background( &
         size(scenario_from%aero_background)))
    do i = 1,size(scenario_from%aero_background)
       call aero_dist_allocate(scenario_to%aero_background(i))
       call aero_dist_copy(scenario_from%aero_background(i), &
            scenario_to%aero_background(i))
    end do

  end subroutine scenario_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the time-dependent contents of the
  !> environment. Thereafter scenario_update_env_state() should be used.
  subroutine scenario_init_env_state(scenario, env_state, time)

    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time

    env_state%temp = interp_1d(scenario%temp_time, scenario%temp, time)
    env_state%pressure = interp_1d(scenario%pressure_time, &
         scenario%pressure, time)
    env_state%height = interp_1d(scenario%height_time, scenario%height, time)
    env_state%elapsed_time = time

  end subroutine scenario_init_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update time-dependent contents of the environment.
  !> scenario_init_env_state() should have been called at the start.
  subroutine scenario_update_env_state(scenario, env_state, time, &
       update_rel_humid)

    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Whether to update the relative humidity.
    logical, intent(in) :: update_rel_humid
    
    !> Ambient water vapor pressure (Pa).
    real(kind=dp) :: pmv_old, pmv_new
    !> Ambient pressure (Pa)
    real(kind=dp) :: pressure_old
    !> Ambient temperature (K)
    real(kind=dp) :: temp_old

    ! Update temperature and pressure and adjust relative humidity to maintain
    ! water mixing ratio.

    pmv_old = env_state_sat_vapor_pressure(env_state) * env_state%rel_humid
    pressure_old = env_state%pressure
    temp_old = env_state%temp

    env_state%temp = interp_1d(scenario%temp_time, scenario%temp, time)
    env_state%pressure = interp_1d(scenario%pressure_time, &
         scenario%pressure, time)

    pmv_new = pmv_old * env_state%pressure / pressure_old

    if (update_rel_humid) then
       env_state%rel_humid = pmv_new / env_state_sat_vapor_pressure(env_state)
    end if

    env_state%height = interp_1d(scenario%height_time, scenario%height, time)
    env_state%elapsed_time = time

  end subroutine scenario_update_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do gas emissions and background dilution.
  !!
  !! Emissions are given as an areal rate \f$e(t)\f$, and then divided
  !! by the current box height \f$h(t)\f$ to obtain a volume
  !! rate. There is also a dimensionless rate scaling \f$r(t)\f$. All
  !! input functions are asusumed constant over the timestep, so the
  !! concentration \f$c(t)\f$ change is given by
  !! \f[
  !!     c(t) = c(0) + \frac{r t}{h} e.
  !! \f]
  !!
  !! We model dilution by considering a gas concentration \f$c(t)\f$
  !! in a box of height \f$h(t)\f$, subject to first-order dilution
  !! with a rate \f$\lambda(t)\f$. Then the effective dilution rate is
  !! \f[
  !!     \lambda_{\rm eff}(t) = \lambda(t) + \frac{\dot{h}(t)}{h(t)}
  !! \f]
  !! and the evolution of \f$c(t)\f$ is given by
  !! \f[
  !!     \dot{c}(t) = - \lambda_{\rm eff}(t) c(t).
  !! \f]
  !! Solving this with separation of variables gives
  !! \f[
  !!     \frac{c(t)}{c(0)} = \frac{h(0)}{h(t)}
  !!                         \exp\left( - \int_0^t \lambda(t)\,dt\right).
  !! \f]
  !! If we define \f$p = c(t)/c(0)\f$ to be the remaining proportion
  !! of the initial concentration, and \f$b\f$ to be the constant
  !! background concentration, then we have
  !! \f[
  !!     c(t) = p(t) c(0) + (1 - p(t)) b.
  !! \f]
  !! We assume constant \f$\lambda\f$ and we only do entrainment with
  !! increasing height \f$h(t)\f$, so we have
  !! \f[
  !!     p(t) = \min\left(1, \frac{h(0)}{h(t)}\right) \exp(-\lambda t).
  !! \f]
  subroutine scenario_update_gas_state(scenario, delta_t, env_state, &
       old_env_state, gas_data, gas_state)

    !> Scenario.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state to update.
    type(gas_state_t), intent(inout) :: gas_state

    real(kind=dp) :: emission_rate_scale, dilution_rate, p
    type(gas_state_t) :: emissions, background

    ! emissions
    call gas_state_allocate_size(emissions, gas_data%n_spec)
    call gas_state_interp_1d(scenario%gas_emission, &
         scenario%gas_emission_time, scenario%gas_emission_rate_scale, &
         env_state%elapsed_time, emissions, emission_rate_scale)
    call gas_state_mole_dens_to_ppb(emissions, env_state)
    p = emission_rate_scale * delta_t / env_state%height
    call gas_state_add_scaled(gas_state, emissions, p)
    call gas_state_deallocate(emissions)

    ! dilution
    call gas_state_allocate_size(background, gas_data%n_spec)
    call gas_state_interp_1d(scenario%gas_background, &
         scenario%gas_dilution_time, scenario%gas_dilution_rate, &
         env_state%elapsed_time, background, dilution_rate)
    p = exp(- dilution_rate * delta_t)
    if (env_state%height > old_env_state%height) then
       p = p * old_env_state%height / env_state%height
    end if
    call gas_state_scale(gas_state, p)
    call gas_state_add_scaled(gas_state, background, 1d0 - p)
    call gas_state_ensure_nonnegative(gas_state)
    call gas_state_deallocate(background)

  end subroutine scenario_update_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution for a particle aerosol
  !> distribution.
  !!
  !! See scenario_update_gas_state() for a description of the
  !! model. We additionally scale the number concentration to account
  !! for temperature changes.
  subroutine scenario_update_aero_state(scenario, delta_t, env_state, &
       old_env_state, aero_data, aero_state, n_emit, n_dil_in, n_dil_out)

    !> Scenario.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero state to update.
    type(aero_state_t), intent(inout) :: aero_state
    !> Number of emitted particles.
    integer, intent(out) :: n_emit
    !> Number of diluted-in particles.
    integer, intent(out) :: n_dil_in
    !> Number of diluted-out particles.
    integer, intent(out) :: n_dil_out

    real(kind=dp) :: emission_rate_scale, dilution_rate, p
    type(aero_dist_t) :: emissions, background
    type(aero_state_t) :: aero_state_delta

    ! emissions
    call aero_dist_allocate(emissions)
    call aero_dist_interp_1d(scenario%aero_emission, &
         scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
         env_state%elapsed_time, emissions, emission_rate_scale)
    p = emission_rate_scale * delta_t / env_state%height
    call aero_state_add_aero_dist_sample(aero_state, aero_data, &
         emissions, p, env_state%elapsed_time, n_emit)
    call aero_dist_deallocate(emissions)

    ! dilution
    call aero_dist_allocate(background)
    call aero_dist_interp_1d(scenario%aero_background, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, &
         env_state%elapsed_time, background, dilution_rate)
    p = exp(- dilution_rate * delta_t)
    if (env_state%height > old_env_state%height) then
       p = p * old_env_state%height / env_state%height
    end if
    ! loss to background
    call aero_state_allocate(aero_state_delta)
    call aero_state_sample_particles(aero_state, aero_state_delta, &
         1d0 - p, AERO_INFO_DILUTION)
    n_dil_out = aero_state_total_particles(aero_state_delta)
    call aero_state_deallocate(aero_state_delta)
    ! addition from background
    call aero_state_add_aero_dist_sample(aero_state, aero_data, &
         background, 1d0 - p, env_state%elapsed_time, n_dil_in)
    call aero_dist_deallocate(background)

    ! update computational volume
    call aero_weight_array_scale(aero_state%awa, &
         old_env_state%temp * env_state%pressure &
         / (env_state%temp * old_env_state%pressure))

  end subroutine scenario_update_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do emissions and background dilution from the environment for a
  !> binned aerosol distribution.
  !!
  !! See scenario_update_gas_state() for a description of the model.
  subroutine scenario_update_aero_binned(scenario, delta_t, env_state, &
       old_env_state, bin_grid, aero_data, aero_binned)

    !> Scenario.
    type(scenario_t), intent(in) :: scenario
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Current environment.
    type(env_state_t), intent(in) :: env_state
    !> Previous environment.
    type(env_state_t), intent(in) :: old_env_state
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero data values.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero binned to update.
    type(aero_binned_t), intent(inout) :: aero_binned

    real(kind=dp) :: emission_rate_scale, dilution_rate, p
    type(aero_dist_t) :: emissions, background
    type(aero_binned_t) :: emissions_binned, background_binned

    ! emissions
    call aero_dist_allocate(emissions)
    call aero_binned_allocate_size(emissions_binned, bin_grid%n_bin, &
         aero_data%n_spec)
    call aero_dist_interp_1d(scenario%aero_emission, &
         scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
         env_state%elapsed_time, emissions, emission_rate_scale)
    call aero_binned_add_aero_dist(emissions_binned, bin_grid, aero_data, &
         emissions)
    p = emission_rate_scale * delta_t / env_state%height
    call aero_binned_add_scaled(aero_binned, emissions_binned, p)
    call aero_dist_deallocate(emissions)
    call aero_binned_deallocate(emissions_binned)

    ! dilution
    call aero_dist_allocate(background)
    call aero_binned_allocate_size(background_binned, bin_grid%n_bin, &
         aero_data%n_spec)
    call aero_dist_interp_1d(scenario%aero_background, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, &
         env_state%elapsed_time, background, dilution_rate)
    call aero_binned_add_aero_dist(background_binned, bin_grid, aero_data, &
         background)
    p = exp(- dilution_rate * delta_t)
    if (env_state%height > old_env_state%height) then
       p = p * old_env_state%height / env_state%height
    end if
    call aero_binned_scale(aero_binned, p)
    call aero_binned_add_scaled(aero_binned, background_binned, 1d0 - p)
    call aero_dist_deallocate(background)
    call aero_binned_deallocate(background_binned)

  end subroutine scenario_update_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evaluate a loss rate function
  real(kind=dp) function scenario_loss_rate(function_id, vol, density, &
       aero_data, temp, press)
       
    !> Id of loss rate function to be used
    integer, intent(in) :: function_id
    !> Volume of particle.
    real(kind=dp), intent(in) :: vol
    !> Density of particle.
    real(kind=dp), intent(in) :: density
    !> Aerosol data. 
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    
    if(function_id == SCENARIO_LOSS_FUNCTION_ZERO) then
      scenario_loss_rate = 0d0
    elseif(function_id == SCENARIO_LOSS_FUNCTION_CONSTANT) then
      scenario_loss_rate = 4d-4
    elseif(function_id == SCENARIO_LOSS_FUNCTION_VOLUME) then
      scenario_loss_rate = vol
    else
       call die_msg(200724934, "Unknown loss function id: " &
            // trim(integer_to_string(function_id)))
    end if
    
  end function scenario_loss_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute and return the max loss rate function for a given volume
  real(kind=dp) function scenario_loss_rate_max(function_id, vol, &
      aero_data, temp, press)

    !> Id of loss rate function to be used
    integer, intent(in) :: function_id
    !> Particle volume.
    real(kind=dp), intent(in) :: vol
    !> Aerosol data. 
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press

    !> Number of density sample points.
    integer, parameter :: n_sample = 3

    real(kind=dp) :: d, d_min, d_max, loss
    integer :: i
    
    d_min = minval(aero_data%density)
    d_max = maxval(aero_data%density)

    scenario_loss_rate_max = 0d0
    do i = 1,n_sample
      d = interp_linear_disc(d_min, d_max, n_sample, i)
      loss = scenario_loss_rate(function_id, vol, d, aero_data, temp, press)
      scenario_loss_rate_max = max(scenario_loss_rate_max, loss)
    end do
  end function scenario_loss_rate_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute an upper bound on the maximum kernel value for each
  !> bin.  Value over_scale is multiplied to the maximum sampled value
  !> to get the upper bound.  A tighter bound may be reached if over_scale
  !> is smaller, but that also risks falling below a kernel value.
  subroutine scenario_loss_rate_bin_max(function_id, bin_grid, &
        aero_data, temp, press, loss_max)
       
    !> Id of loss rate function to be used
    integer, intent(in) :: function_id
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press
    !> Maximum loss vals.
    real(kind=dp), intent(out) :: loss_max(bin_grid%n_bin)
    
    !> Number of sample points per bin.
    integer, parameter :: n_sample = 3
    !> Over-estimation scale factor parameter.
    real(kind=dp), parameter :: over_scale = 2d0
    
    real(kind=dp) :: v_low, v_high, vol, r, r_max
    integer :: b, i
    
    do b = 1,bin_grid%n_bin
      v_low = rad2vol(bin_grid%edge_radius(b))
      v_high = rad2vol(bin_grid%edge_radius(b + 1))
      r_max = 0d0
      do i = 1,n_sample
        vol = interp_linear_disc(v_low, v_high, n_sample, i)
        r = scenario_loss_rate_max(function_id, vol, aero_data, temp, press)
        r_max = max(r_max, r)
      end do
      loss_max(b) = r_max*over_scale
    end do
    
  end subroutine scenario_loss_rate_bin_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Performs stochastic particle loss for one time-step.
  !> If a particle p has a scenario_loss_rate(...) value of rate, then the
  !> probability p will be removed by this function is 1 - exp(-delta_t*rate).
  !> Uses an accept-reject algorithm for efficiency, in which a particle
  !> is first sampled with rate 1 - exp(-delta_t*over_rate)
  !> and then accepted with rate
  !> (1 - exp(-delta_t*rate))/(1 - exp(-delta_t*over_rate)).
  subroutine scenario_particle_loss(function_id, delta_t, aero_data, &
       aero_state, temp, press)

    !> Id of loss rate function to be used
    integer, intent(in) :: function_id
    !> Time increment to update over.
    real(kind=dp), intent(in) :: delta_t
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: press

    type(aero_particle_t), pointer :: aero_particle
    type(aero_info_t) :: aero_info
    integer :: c, b, s, p
    real(kind=dp) :: rate, vol, density, over_rate, over_prob, &
        rand_real, rand_geom
    
    !integer :: init_size, candidates, cand_iter
    
    if(function_id == SCENARIO_LOSS_FUNCTION_ZERO .or. &
        function_id == SCENARIO_LOSS_FUNCTION_INVALID) return
    
    call aero_state_sort(aero_state)
    
    if (.not. aero_state%aero_sorted%removal_rate_bounds_valid) then
      call scenario_loss_rate_bin_max(function_id, &
          aero_state%aero_sorted%bin_grid, aero_data, temp, press, &
          aero_state%aero_sorted%removal_rate_max)
      aero_state%aero_sorted%removal_rate_bounds_valid = .true.
    end if
    
    do c = 1,aero_sorted_n_class(aero_state%aero_sorted)
      do b = 1,aero_sorted_n_bin(aero_state%aero_sorted)
        s = aero_state%aero_sorted%size_class%inverse(b, c)%n_entry + 1
        over_rate = aero_state%aero_sorted%removal_rate_max(b)
        if (delta_t*over_rate <= 0d0) cycle
        over_prob = 1d0 - exp(-delta_t*over_rate)
        do while (.TRUE.)
          rand_real = pmc_random()
          if (rand_real <= 0d0) exit
          rand_geom = -log(rand_real)/(delta_t*over_rate) + 1d0
          if (rand_geom >= s) exit
          s = s - floor(rand_geom)
          
          ! note: floor(rand_geom) is a random geometric variable
          ! with accept probability 1 - exp(-delta*over_rate)
          
          p = aero_state%aero_sorted%size_class%inverse(b, c)%entry(s)
          aero_particle => aero_state%apa%particle(p)
          vol = aero_particle_volume(aero_particle)
          density = aero_particle_density(aero_particle, aero_data)
          rate = scenario_loss_rate(function_id, vol, density, aero_data, &
                  temp, press)
          call warn_assert_msg(295846288, rate <= over_rate, &
              "particle loss upper bound estimation is too tight: " &
              // trim(real_to_string(rate)) // " > " &
              // trim(real_to_string(over_rate)) )
          if (pmc_random()*over_prob > 1d0 - exp(-delta_t*rate)) cycle
          
          call aero_info_allocate(aero_info)
          aero_info%id = aero_particle%id
          aero_info%action = AERO_INFO_DILUTION
          aero_info%other_id = 0
          call aero_state_remove_particle_with_info(aero_state, p, aero_info)
          call aero_info_deallocate(aero_info)
        end do
      end do
    end do
    
    !call aero_state_check_sort(aero_state)

  end subroutine scenario_particle_loss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether any of the contained aerosol modes are of the given type.
  elemental logical function scenario_contains_aero_mode_type(scenario, &
       aero_mode_type)

    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Aerosol mode type to test for.
    integer, intent(in) :: aero_mode_type

    scenario_contains_aero_mode_type &
         = any(aero_dist_contains_aero_mode_type(scenario%aero_emission, &
         aero_mode_type)) &
         .or. any(aero_dist_contains_aero_mode_type(scenario%aero_background, &
         aero_mode_type))

  end function scenario_contains_aero_mode_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read environment data from an spec file.
  subroutine spec_file_read_scenario(file, gas_data, aero_data, scenario)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Environment data.
    type(scenario_t), intent(inout) :: scenario

    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file

    ! note that we have to hard-code the list for doxygen below

    !> \page input_format_scenario Input File Format: Scenario
    !!
    !! The scenario parameters are:
    !! <ul>
    !! <li> \b temp_profile (string): the name of the file from which to
    !!      read the temperature profile --- the file format should be
    !!      \subpage input_format_temp_profile
    !! <li> \b pressure_profile (string): the name of the file from which to
    !!      read the pressure profile --- the file format should be
    !!      \subpage input_format_pressure_profile
    !! <li> \b height_profile (string): the name of the file from which
    !!      to read the mixing layer height profile --- the file format
    !!      should be \subpage input_format_height_profile
    !! <li> \b gas_emissions (string): the name of the file from which to
    !!      read the gas emissions profile --- the file format should be
    !!      \subpage input_format_gas_profile
    !! <li> \b gas_background (string): the name of the file from which
    !!      to read the gas background profile --- the file format should
    !!      be \subpage input_format_gas_profile
    !! <li> \b aero_emissions (string): the name of the file from which
    !!      to read the aerosol emissions profile --- the file format
    !!      should be \subpage input_format_aero_dist_profile
    !! <li> \b aero_background (string): the name of the file from which
    !!      to read the aerosol background profile --- the file format
    !!      should be \subpage input_format_aero_dist_profile
    !! </ul>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    !> \page input_format_temp_profile Input File Format: Temperature Profile
    !!
    !! A temperature profile input file must consist of two lines:
    !! - the first line must begin with \c time and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the times (in
    !!   s after the start of the simulation) of the temperature set
    !!   points --- the times must be in increasing order
    !! - the second line must begin with \c temp and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the
    !!   temperatures (in K) at the corresponding times
    !!
    !! The temperature profile is linearly interpolated between the
    !! specified times, while before the first time it takes the first
    !! temperature value and after the last time it takes the last
    !! temperature value.
    !!
    !! Example:
    !! <pre>
    !! time  0    600  1800  # time (in s) after simulation start
    !! temp  270  290  280   # temperature (in K)
    !! </pre>
    !! Here the temperature starts at 270&nbsp;K at the start of the
    !! simulation, rises to 290&nbsp;K after 10&nbsp;min, and then
    !! falls again to 280&nbsp;K at 30&nbsp;min. Between these times
    !! the temperature is linearly interpolated, while after
    !! 30&nbsp;min it is held constant at 280&nbsp;K.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_scenario --- the environment data
    !!     containing the temperature profile

    !> \page input_format_pressure_profile Input File Format: Pressure Profile
    !!
    !! A pressure profile input file must consist of two lines:
    !! - the first line must begin with \c time and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the times (in
    !!   s after the start of the simulation) of the pressure set
    !!   points --- the times must be in increasing order
    !! - the second line must begin with \c pressure and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the
    !!   pressures (in Pa) at the corresponding times
    !!
    !! The pressure profile is linearly interpolated between the
    !! specified times, while before the first time it takes the first
    !! pressure value and after the last time it takes the last
    !! pressure value.
    !!
    !! Example:
    !! <pre>
    !! time      0    600  1800  # time (in s) after simulation start
    !! pressure  1e5  9e4  7.5e4 # pressure (in Pa)
    !! </pre>
    !! Here the pressure starts at 1e5&nbsp;Pa at the start of the
    !! simulation, decreases to 9e4&nbsp;Pa after 10&nbsp;min, and then
    !! decreases further to 7.5e4&nbsp;Pa at 30&nbsp;min. Between these times
    !! the pressure is linearly interpolated, while after
    !! 30&nbsp;min it is held constant at 7.5e4&nbsp;Pa.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_scenario --- the environment data
    !!     containing the pressure profile

    !> \page input_format_height_profile Input File Format: Mixing Layer Height Profile
    !!
    !! A mixing layer height profile input file must consist of two
    !! lines:
    !! - the first line must begin with \c time and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the times (in
    !!   s after the start of the simulation) of the height set
    !!   points --- the times must be in increasing order
    !! - the second line must begin with \c height and should be
    !!   followed by \f$N\f$ space-separated real scalars, giving the
    !!   mixing layer heights (in m) at the corresponding times
    !!
    !! The mixing layer height profile is linearly interpolated
    !! between the specified times, while before the first time it
    !! takes the first height value and after the last time it takes
    !! the last height value.
    !!
    !! Example:
    !! <pre>
    !! time    0    600   1800  # time (in s) after simulation start
    !! height  500  1000  800   # mixing layer height (in m)
    !! </pre>
    !! Here the mixing layer height starts at 500&nbsp;m at the start
    !! of the simulation, rises to 1000&nbsp;m after 10&nbsp;min, and
    !! then falls again to 800&nbsp;m at 30&nbsp;min. Between these
    !! times the mixing layer height is linearly interpolated, while
    !! after 30&nbsp;min it is held constant at 800&nbsp;m.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_scenario --- the environment data
    !!     containing the mixing layer height profile

    ! temperature profile
    call spec_file_read_string(file, "temp_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "temp", &
         scenario%temp_time, scenario%temp)
    call spec_file_close(sub_file)

    ! pressure profile
    call spec_file_read_string(file, "pressure_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "pressure", &
         scenario%pressure_time, scenario%pressure)
    call spec_file_close(sub_file)

    ! height profile
    call spec_file_read_string(file, "height_profile", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_timed_real_array(sub_file, "height", &
         scenario%height_time, scenario%height)
    call spec_file_close(sub_file)

    ! gas emissions profile
    call spec_file_read_string(file, "gas_emissions", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
         scenario%gas_emission_time, scenario%gas_emission_rate_scale, &
         scenario%gas_emission)
    call spec_file_close(sub_file)

    ! gas background profile
    call spec_file_read_string(file, "gas_background", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_states_times_rates(sub_file, gas_data, &
         scenario%gas_dilution_time, scenario%gas_dilution_rate, &
         scenario%gas_background)
    call spec_file_close(sub_file)

    ! aerosol emissions profile
    call spec_file_read_string(file, "aero_emissions", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dists_times_rates(sub_file, aero_data, &
         scenario%aero_emission_time, scenario%aero_emission_rate_scale, &
         scenario%aero_emission)
    call spec_file_close(sub_file)

    ! aerosol background profile
    call spec_file_read_string(file, "aero_background", sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dists_times_rates(sub_file, aero_data, &
         scenario%aero_dilution_time, scenario%aero_dilution_rate, &
         scenario%aero_background)
    call spec_file_close(sub_file)

  end subroutine spec_file_read_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_scenario(val)

    !> Value to pack.
    type(scenario_t), intent(in) :: val

    integer :: total_size, i, n

    total_size = &
         pmc_mpi_pack_size_real_array(val%temp_time) &
         + pmc_mpi_pack_size_real_array(val%temp) &
         + pmc_mpi_pack_size_real_array(val%pressure_time) &
         + pmc_mpi_pack_size_real_array(val%pressure) &
         + pmc_mpi_pack_size_real_array(val%height_time) &
         + pmc_mpi_pack_size_real_array(val%height) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_time) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_rate_scale) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_time) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_rate) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_time) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_rate_scale) &
         + pmc_mpi_pack_size_real_array(val%aero_dilution_time) &
         + pmc_mpi_pack_size_real_array(val%aero_dilution_rate)
    do i = 1,size(val%gas_emission)
       total_size = total_size &
            + pmc_mpi_pack_size_gas_state(val%gas_emission(i))
    end do
    do i = 1,size(val%gas_background)
       total_size = total_size &
            + pmc_mpi_pack_size_gas_state(val%gas_background(i))
    end do
    do i = 1,size(val%aero_emission)
       total_size = total_size &
            + pmc_mpi_pack_size_aero_dist(val%aero_emission(i))
    end do
    do i = 1,size(val%aero_background)
       total_size = total_size &
            + pmc_mpi_pack_size_aero_dist(val%aero_background(i))
    end do

    pmc_mpi_pack_size_scenario = total_size

  end function pmc_mpi_pack_size_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_scenario(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(scenario_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_pack_real_array(buffer, position, val%temp)
    call pmc_mpi_pack_real_array(buffer, position, val%pressure_time)
    call pmc_mpi_pack_real_array(buffer, position, val%pressure)
    call pmc_mpi_pack_real_array(buffer, position, val%height_time)
    call pmc_mpi_pack_real_array(buffer, position, val%height)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_rate_scale)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, &
         val%aero_emission_rate_scale)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_dilution_time)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_dilution_rate)
    do i = 1,size(val%gas_emission)
       call pmc_mpi_pack_gas_state(buffer, position, val%gas_emission(i))
    end do
    do i = 1,size(val%gas_background)
       call pmc_mpi_pack_gas_state(buffer, position, val%gas_background(i))
    end do
    do i = 1,size(val%aero_emission)
       call pmc_mpi_pack_aero_dist(buffer, position, val%aero_emission(i))
    end do
    do i = 1,size(val%aero_background)
       call pmc_mpi_pack_aero_dist(buffer, position, val%aero_background(i))
    end do
    call assert(639466930, &
         position - prev_position <= pmc_mpi_pack_size_scenario(val))
#endif

  end subroutine pmc_mpi_pack_scenario

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_scenario(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(scenario_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    call scenario_deallocate(val)
    call scenario_allocate(val)
    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%temp)
    call pmc_mpi_unpack_real_array(buffer, position, val%pressure_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%pressure)
    call pmc_mpi_unpack_real_array(buffer, position, val%height_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%height)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, &
         val%gas_emission_rate_scale)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, &
         val%aero_emission_rate_scale)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_rate)
    deallocate(val%gas_emission)
    deallocate(val%gas_background)
    deallocate(val%aero_emission)
    deallocate(val%aero_background)
    allocate(val%gas_emission(size(val%gas_emission_time)))
    allocate(val%gas_background(size(val%gas_dilution_time)))
    allocate(val%aero_emission(size(val%aero_emission_time)))
    allocate(val%aero_background(size(val%aero_dilution_time)))
    do i = 1,size(val%gas_emission)
       call gas_state_allocate(val%gas_emission(i))
       call pmc_mpi_unpack_gas_state(buffer, position, val%gas_emission(i))
    end do
    do i = 1,size(val%gas_background)
       call gas_state_allocate(val%gas_background(i))
       call pmc_mpi_unpack_gas_state(buffer, position, val%gas_background(i))
    end do
    do i = 1,size(val%aero_emission)
       call aero_dist_allocate(val%aero_emission(i))
       call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_emission(i))
    end do
    do i = 1,size(val%aero_background)
       call aero_dist_allocate(val%aero_background(i))
       call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_background(i))
    end do
    call assert(611542570, &
         position - prev_position <= pmc_mpi_pack_size_scenario(val))
#endif

  end subroutine pmc_mpi_unpack_scenario
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a loss function type from a spec file and
  !> generate it.
  subroutine spec_file_read_loss_function_type(file, loss_function_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Function type.
    integer, intent(out) :: loss_function_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: function_name

    !> \page input_format_loss_function Input File Format:
    !!     Loss Rate Function
    !!
    !! The loss rate function is specified by the parameter:
    !!   - \b loss_function (string): the type of loss function ---
    !!     must be one of: \c zero for no particle loss, or \c volume
    !!     for particle loss proportional to particle volume
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'loss_function', function_name)
    if (trim(function_name) == 'zero') then
       loss_function_type = SCENARIO_LOSS_FUNCTION_ZERO
    elseif (trim(function_name) == 'constant') then
       loss_function_type = SCENARIO_LOSS_FUNCTION_CONSTANT
    elseif (trim(function_name) == 'volume') then
       loss_function_type = SCENARIO_LOSS_FUNCTION_VOLUME
    else
       call spec_file_die_msg(518248400, file, &
            "Unknown loss function type: " // trim(function_name))
    end if

  end subroutine spec_file_read_loss_function_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_scenario
