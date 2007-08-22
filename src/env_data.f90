! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Environment data. This is everything needed to support the current
! environment state. The environment data is not time-dependent,
! whereas the environment is time-dependent.
!
! The temperature, emissions and background states are profiles
! proscribed as functions of time by giving a number of times and the
! corresponding data. Linear interpolation is used between the times,
! with constant interpolation outside of the range of times.

module pmc_env_data

  use pmc_gas_state
  use pmc_aero_dist
  
  type env_data_t
     real*8, pointer :: temp_time(:)    ! times at temp set-points (s)
     real*8, pointer :: temp(:)         ! temps at set-points (K)

     real*8, pointer :: height_time(:)  ! times at height set-points (s)
     real*8, pointer :: height(:)       ! heights at set-points (m)

     real*8, pointer :: gas_emission_time(:) ! gas emissions times (s)
     real*8, pointer :: gas_emission_rate(:) ! gas emisssion rates (s^{-1})
     type(gas_state_t), pointer :: gas_emission(:) ! gas emissions

     real*8, pointer :: gas_dilution_time(:) ! gas-backgnd dilute times (s)
     real*8, pointer :: gas_dilution_rate(:) ! gas-backgnd dilute rates (s^{-1})
     type(gas_state_t), pointer :: gas_background(:) ! background gas concs

     real*8, pointer :: aero_emission_time(:) ! aerosol emissions times (s)
     real*8, pointer :: aero_emission_rate(:) ! aerosol emisssion rates (s^{-1})
     type(aero_dist_t), pointer :: aero_emission(:) ! aerosol emissions

     real*8, pointer :: aero_dilution_time(:) ! aero-backgnd dilute times (s)
     real*8, pointer :: aero_dilution_rate(:) ! aero-bkgd dilute rates (s^{-1})
     type(aero_dist_t), pointer :: aero_background(:) ! aerosol backgrounds
  end type env_data_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_data_alloc(env_data)

    ! Allocate an empty environment.

    type(env_data_t), intent(out) :: env_data   ! environment data

    allocate(env_data%temp_time(0))
    allocate(env_data%temp(0))

    allocate(env_data%height_time(0))
    allocate(env_data%height(0))

    allocate(env_data%gas_emission_time(0))
    allocate(env_data%gas_emission_rate(0))
    allocate(env_data%gas_emission(0))

    allocate(env_data%gas_dilution_time(0))
    allocate(env_data%gas_dilution_rate(0))
    allocate(env_data%gas_background(0))

    allocate(env_data%aero_emission_time(0))
    allocate(env_data%aero_emission_rate(0))
    allocate(env_data%aero_emission(0))

    allocate(env_data%aero_dilution_time(0))
    allocate(env_data%aero_dilution_rate(0))
    allocate(env_data%aero_background(0))

  end subroutine env_data_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine env_data_free(env_data)

    ! Free all storage.

    type(env_data_t), intent(out) :: env_data   ! environment data

    integer :: i

    deallocate(env_data%temp_time)
    deallocate(env_data%temp)

    deallocate(env_data%height_time)
    deallocate(env_data%height)

    do i = 1,size(env_data%gas_emission)
       call gas_state_free(env_data%gas_emission(i))
    end do
    deallocate(env_data%gas_emission_time)
    deallocate(env_data%gas_emission_rate)
    deallocate(env_data%gas_emission)

    do i = 1,size(env_data%gas_background)
       call gas_state_free(env_data%gas_background(i))
    end do
    deallocate(env_data%gas_dilution_time)
    deallocate(env_data%gas_dilution_rate)
    deallocate(env_data%gas_background)

    do i = 1,size(env_data%aero_emission)
       call aero_dist_free(env_data%aero_emission(i))
    end do
    deallocate(env_data%aero_emission_time)
    deallocate(env_data%aero_emission_rate)
    deallocate(env_data%aero_emission)

    do i = 1,size(env_data%aero_background)
       call aero_dist_free(env_data%aero_background(i))
    end do
    deallocate(env_data%aero_dilution_time)
    deallocate(env_data%aero_dilution_rate)
    deallocate(env_data%aero_background)

  end subroutine env_data_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine env_data_init_state(env_data, env, time)
    
    ! Initialize the time-dependent contents of the
    ! environment. Thereafter env_data_update_state() should be used.

    use pmc_util
    use pmc_env

    type(env_data_t), intent(in) :: env_data ! environment data
    type(env_t), intent(inout) :: env   ! environment state to update
    real*8, intent(in) :: time          ! current time (s)

    ! init temperature
    env%temp = interp_1d(env_data%temp_time, env_data%temp, time)

    ! init height
    env%height = interp_1d(env_data%height_time, env_data%height, time)

    ! init gas and aerosol emissions and background
    call gas_state_interp_1d(env_data%gas_emission, &
         env_data%gas_emission_time, env_data%gas_emission_rate, &
         time, env%gas_emissions, env%gas_emission_rate)
    call gas_state_interp_1d(env_data%gas_background, &
         env_data%gas_dilution_time, env_data%gas_dilution_rate, &
         time, env%gas_background, env%gas_dilution_rate)
    call aero_dist_interp_1d(env_data%aero_emission, &
         env_data%aero_emission_time, env_data%aero_emission_rate, &
         time, env%aero_emissions, env%aero_emission_rate)
    call aero_dist_interp_1d(env_data%aero_background, &
         env_data%aero_dilution_time, env_data%aero_dilution_rate, &
         time, env%aero_background, env%aero_dilution_rate)
    
  end subroutine env_data_init_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine env_data_update_state(env_data, env, time)
    
    ! Update time-dependent contents of the environment.
    ! env_data_init_state() should have been called at the start.

    use pmc_util
    use pmc_env

    type(env_data_t), intent(in) :: env_data ! environment data
    type(env_t), intent(inout) :: env   ! environment state to update
    real*8, intent(in) :: time          ! current time (s)
    
    real*8 :: pmv ! ambient water vapor pressure (Pa)
    real*8 :: old_height

    ! update temperature and relative humidity
    pmv = env_sat_vapor_pressure(env) * env%rel_humid
    env%temp = interp_1d(env_data%temp_time, env_data%temp, time)
    env%rel_humid = pmv / env_sat_vapor_pressure(env)

    ! update height
    env%height = interp_1d(env_data%height_time, env_data%height, time)

    ! update gas and aerosol emissions and background
    call gas_state_interp_1d(env_data%gas_emission, &
         env_data%gas_emission_time, env_data%gas_emission_rate, &
         time, env%gas_emissions, env%gas_emission_rate)
    call gas_state_interp_1d(env_data%gas_background, &
         env_data%gas_dilution_time, env_data%gas_dilution_rate, &
         time, env%gas_background, env%gas_dilution_rate)
    call aero_dist_interp_1d(env_data%aero_emission, &
         env_data%aero_emission_time, env_data%aero_emission_rate, &
         time, env%aero_emissions, env%aero_emission_rate)
    call aero_dist_interp_1d(env_data%aero_background, &
         env_data%aero_dilution_time, env_data%aero_dilution_rate, &
         time, env%aero_background, env%aero_dilution_rate)

  end subroutine env_data_update_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_env_data(file, env_data)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(env_data_t), intent(in) :: env_data ! environment data to write

    integer :: i
    
    call inout_write_comment(file, "begin env_data")

    call inout_write_real_array(file, "temp_time(s)", env_data%temp_time)
    call inout_write_real_array(file, "temp(K)", env_data%temp)
    
    call inout_write_real_array(file, "height_time(s)", env_data%height_time)
    call inout_write_real_array(file, "height(m)", env_data%height)
    
    call inout_write_integer(file, 'n_gas_emit', &
         size(env_data%gas_emission_time))
    do i = 1,size(env_data%gas_emission_time)
       call inout_write_integer(file, "gas_emit_num", i)
       call inout_write_real(file, "gas_emit_time(s)", &
            env_data%gas_emission_time(i))
       call inout_write_real(file, "gas_emit_rate(1/s)", &
            env_data%gas_emission_rate(i))
       call inout_write_gas_state(file, env_data%gas_emission(i))
    end do

    call inout_write_integer(file, 'n_gas_dilute', &
         size(env_data%gas_dilution_time))
    do i = 1,size(env_data%gas_dilution_time)
       call inout_write_integer(file, "gas_dilute_num", i)
       call inout_write_real(file, "gas_dilute_time(s)", &
            env_data%gas_dilution_time(i))
       call inout_write_real(file, "gas_dilute_rate(1/s)", &
            env_data%gas_dilution_rate(i))
       call inout_write_gas_state(file, env_data%gas_background(i))
    end do

    call inout_write_integer(file, 'n_aero_emit', &
         size(env_data%aero_emission_time))
    do i = 1,size(env_data%aero_emission_time)
       call inout_write_integer(file, "aero_emit_num", i)
       call inout_write_real(file, "aero_emit_time(s)", &
            env_data%aero_emission_time(i))
       call inout_write_real(file, "aero_emit_rate(1/s)", &
            env_data%aero_emission_rate(i))
       call inout_write_aero_dist(file, env_data%aero_emission(i))
    end do

    call inout_write_integer(file, 'n_aero_dilute', &
         size(env_data%aero_dilution_time))
    do i = 1,size(env_data%aero_dilution_time)
       call inout_write_integer(file, "aero_dilute_num", i)
       call inout_write_real(file, "aero_dilute_time(s)", &
            env_data%aero_dilution_time(i))
       call inout_write_real(file, "aero_dilute_rate(1/s)", &
            env_data%aero_dilution_rate(i))
       call inout_write_aero_dist(file, env_data%aero_background(i))
    end do

    call inout_write_comment(file, "end env_data")

  end subroutine inout_write_env_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_env_data(file, env_data)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(env_data_t), intent(out) :: env_data ! environment data to read
    
    integer :: i, n

    call inout_check_comment(file, "begin env_data")

    call inout_read_real_array(file, "temp_time(s)", env_data%temp_time)
    call inout_read_real_array(file, "temp(K)", env_data%temp)
    
    call inout_read_real_array(file, "height_time(s)", env_data%height_time)
    call inout_read_real_array(file, "height(m)", env_data%height)
    
    call inout_read_integer(file, 'n_gas_emit', n)
    allocate(env_data%gas_emission_time(n))
    allocate(env_data%gas_emission_rate(n))
    allocate(env_data%gas_emission(n))
    do i = 1,size(env_data%gas_emission_time)
       call inout_read_integer(file, "gas_emit_num", i)
       call inout_read_real(file, "gas_emit_time(s)", &
            env_data%gas_emission_time(i))
       call inout_read_real(file, "gas_emit_rate(1/s)", &
            env_data%gas_emission_rate(i))
       call inout_read_gas_state(file, env_data%gas_emission(i))
    end do

    call inout_read_integer(file, 'n_gas_dilute', n)
    allocate(env_data%gas_dilution_time(n))
    allocate(env_data%gas_dilution_rate(n))
    allocate(env_data%gas_background(n))
    do i = 1,size(env_data%gas_dilution_time)
       call inout_read_integer(file, "gas_dilute_num", i)
       call inout_read_real(file, "gas_dilute_time(s)", &
            env_data%gas_dilution_time(i))
       call inout_read_real(file, "gas_dilute_rate(1/s)", &
            env_data%gas_dilution_rate(i))
       call inout_read_gas_state(file, env_data%gas_background(i))
    end do

    call inout_read_integer(file, 'n_aero_emit', n)
    allocate(env_data%aero_emission_time(n))
    allocate(env_data%aero_emission_rate(n))
    allocate(env_data%aero_emission(n))
    do i = 1,size(env_data%aero_emission_time)
       call inout_read_integer(file, "aero_emit_num", i)
       call inout_read_real(file, "aero_emit_time(s)", &
            env_data%aero_emission_time(i))
       call inout_read_real(file, "aero_emit_rate(1/s)", &
            env_data%aero_emission_rate(i))
       call inout_read_aero_dist(file, env_data%aero_emission(i))
    end do

    call inout_read_integer(file, 'n_aero_dilute', n)
    allocate(env_data%aero_dilution_time(n))
    allocate(env_data%aero_dilution_rate(n))
    allocate(env_data%aero_background(n))
    do i = 1,size(env_data%aero_dilution_time)
       call inout_read_integer(file, "aero_dilute_num", i)
       call inout_read_real(file, "aero_dilute_time(s)", &
            env_data%aero_dilution_time(i))
       call inout_read_real(file, "aero_dilute_rate(1/s)", &
            env_data%aero_dilution_rate(i))
       call inout_read_aero_dist(file, env_data%aero_background(i))
    end do

    call inout_check_comment(file, "end env_data")

  end subroutine inout_read_env_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_env_data(file, bin_grid, gas_data, &
       aero_data, env_data)

    ! Read environment data from an inout file.

    use pmc_bin_grid
    use pmc_inout
    use pmc_aero_data
    use pmc_gas_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(gas_data_t), intent(in) :: gas_data ! gas data values
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(env_data_t), intent(out) :: env_data ! environment data

    call inout_read_timed_real_array(file, "temp_profile", "temp", &
         env_data%temp_time, env_data%temp)
    call inout_read_timed_real_array(file, "height_profile", "height", &
         env_data%height_time, env_data%height)
    call spec_read_gas_states_times_rates(file, gas_data, &
         'gas_emissions', env_data%gas_emission_time, &
         env_data%gas_emission_rate, env_data%gas_emission)
    call spec_read_gas_states_times_rates(file, gas_data, &
         'gas_background', env_data%gas_dilution_time, &
         env_data%gas_dilution_rate, env_data%gas_background)
    call spec_read_aero_dists_times_rates(file, aero_data, bin_grid, &
         'aero_emissions', env_data%aero_emission_time, &
         env_data%aero_emission_rate, env_data%aero_emission)
    call spec_read_aero_dists_times_rates(file, aero_data, bin_grid, &
         'aero_background', env_data%aero_dilution_time, &
         env_data%aero_dilution_rate, env_data%aero_background)

  end subroutine spec_read_env_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env_data
