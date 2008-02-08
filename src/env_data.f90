! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_env_data module.

!> The env_data_t structure and associated subroutines.
module pmc_env_data

  use pmc_gas_state
  use pmc_aero_dist
  use pmc_util
  use pmc_env_state
  use pmc_inout
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_gas_data
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif
  
  !> Environment data.
  !!
  !! This is everything needed to support the current environment
  !! state. The environment data is not time-dependent, whereas the
  !! environment state in env_state_t is time-dependent.
  !!
  !! The temperature, emissions and background states are profiles
  !! proscribed as functions of time by giving a number of times and
  !! the corresponding data. Simple data such as temperature is
  !! linearly interpoloated between times, with constant interpolation
  !! outside of the range of times. Gases and aerosols are
  !! interpolated with gas_state_interp_1d() and
  !! aero_dist_interp_1d(), respectively.
  type env_data_t
     !> Temperature set-point times (s).
     real*8, pointer :: temp_time(:)
     !> Temperatures at set-points (K).
     real*8, pointer :: temp(:)

     !> Height set-point times (s).
     real*8, pointer :: height_time(:)
     !> Heights at set-points (m).
     real*8, pointer :: height(:)

     !> Gas emission set-point times (s).
     real*8, pointer :: gas_emission_time(:)
     !> Gas emisssion rates at set-points (s^{-1}).
     real*8, pointer :: gas_emission_rate(:)
     !> Gas emissions at set-points.
     type(gas_state_t), pointer :: gas_emission(:)

     !> Gas-background dilution set-point times (s).
     real*8, pointer :: gas_dilution_time(:)
     !> Gas-background dilution rates at set-points (s^{-1}).
     real*8, pointer :: gas_dilution_rate(:)
     !> Background gas concentrations at set-points.
     type(gas_state_t), pointer :: gas_background(:)

     !> Aerosol emission set-points times (s).
     real*8, pointer :: aero_emission_time(:)
     !> Aerosol emission rates at set-points (s^{-1}).
     real*8, pointer :: aero_emission_rate(:)
     !> Aerosol emissions at set-points.
     type(aero_dist_t), pointer :: aero_emission(:)

     !> Aerosol-background dilution set-point times (s).
     real*8, pointer :: aero_dilution_time(:)
     !> Aerosol-background dilution rates at set-points (s^{-1}).
     real*8, pointer :: aero_dilution_rate(:)
     !> Aerosol background at set-points.
     type(aero_dist_t), pointer :: aero_background(:)
  end type env_data_t
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate an empty environment.
  subroutine env_data_alloc(env_data)

    !> Environment data.
    type(env_data_t), intent(out) :: env_data

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

  !> Free all storage.
  subroutine env_data_free(env_data)

    !> Environment data.
    type(env_data_t), intent(out) :: env_data

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

  !> Initialize the time-dependent contents of the
  !> environment. Thereafter env_data_update_state() should be used.
  subroutine env_data_init_state(env_data, env_state, time)

    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real*8, intent(in) :: time

    ! init temperature
    env_state%temp = interp_1d(env_data%temp_time, env_data%temp, time)

    ! init height
    env_state%height = interp_1d(env_data%height_time, env_data%height, time)

    ! init gas and aerosol emissions and background
    call gas_state_interp_1d(env_data%gas_emission, &
         env_data%gas_emission_time, env_data%gas_emission_rate, &
         time, env_state%gas_emissions, env_state%gas_emission_rate)
    call gas_state_interp_1d(env_data%gas_background, &
         env_data%gas_dilution_time, env_data%gas_dilution_rate, &
         time, env_state%gas_background, env_state%gas_dilution_rate)
    call aero_dist_interp_1d(env_data%aero_emission, &
         env_data%aero_emission_time, env_data%aero_emission_rate, &
         time, env_state%aero_emissions, env_state%aero_emission_rate)
    call aero_dist_interp_1d(env_data%aero_background, &
         env_data%aero_dilution_time, env_data%aero_dilution_rate, &
         time, env_state%aero_background, env_state%aero_dilution_rate)
    
  end subroutine env_data_init_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update time-dependent contents of the environment.
  !> env_data_init_state() should have been called at the start.
  subroutine env_data_update_state(env_data, env_state, time)

    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Environment state to update.
    type(env_state_t), intent(inout) :: env_state
    !> Current time (s).
    real*8, intent(in) :: time
    
    !> Ambient water vapor pressure (Pa).
    real*8 :: pmv
    real*8 :: old_height

    ! update temperature and relative humidity
    pmv = env_state_sat_vapor_pressure(env_state) * env_state%rel_humid
    env_state%temp = interp_1d(env_data%temp_time, env_data%temp, time)
    env_state%rel_humid = pmv / env_state_sat_vapor_pressure(env_state)

    ! update height
    env_state%height = interp_1d(env_data%height_time, env_data%height, time)

    ! update gas and aerosol emissions and background
    call gas_state_interp_1d(env_data%gas_emission, &
         env_data%gas_emission_time, env_data%gas_emission_rate, &
         time, env_state%gas_emissions, env_state%gas_emission_rate)
    call gas_state_interp_1d(env_data%gas_background, &
         env_data%gas_dilution_time, env_data%gas_dilution_rate, &
         time, env_state%gas_background, env_state%gas_dilution_rate)
    call aero_dist_interp_1d(env_data%aero_emission, &
         env_data%aero_emission_time, env_data%aero_emission_rate, &
         time, env_state%aero_emissions, env_state%aero_emission_rate)
    call aero_dist_interp_1d(env_data%aero_background, &
         env_data%aero_dilution_time, env_data%aero_dilution_rate, &
         time, env_state%aero_background, env_state%aero_dilution_rate)

  end subroutine env_data_update_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine inout_write_env_data(file, env_data)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Environment data to write.
    type(env_data_t), intent(in) :: env_data

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

  !> Read full state.
  subroutine inout_read_env_data(file, env_data)
    
    !> File to read from.
    type(inout_file_t), intent(inout) :: file
    !> Environment data to read.
    type(env_data_t), intent(out) :: env_data
    
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

  !> Read environment data from an inout file.
  subroutine spec_read_env_data(file, bin_grid, gas_data, &
       aero_data, env_data)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Gas data values.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment data.
    type(env_data_t), intent(out) :: env_data

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

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_env_data(val)

    !> Value to pack.
    type(env_data_t), intent(in) :: val

    integer :: total_size, i, n

    total_size = &
         pmc_mpi_pack_size_real_array(val%temp_time) &
         + pmc_mpi_pack_size_real_array(val%temp) &
         + pmc_mpi_pack_size_real_array(val%height_time) &
         + pmc_mpi_pack_size_real_array(val%height) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_time) &
         + pmc_mpi_pack_size_real_array(val%gas_emission_rate) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_time) &
         + pmc_mpi_pack_size_real_array(val%gas_dilution_rate) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_time) &
         + pmc_mpi_pack_size_real_array(val%aero_emission_rate) &
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

    pmc_mpi_pack_size_env_data = total_size

  end function pmc_mpi_pack_size_env_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_env_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(env_data_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_pack_real_array(buffer, position, val%temp)
    call pmc_mpi_pack_real_array(buffer, position, val%height_time)
    call pmc_mpi_pack_real_array(buffer, position, val%height)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_emission_rate)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_pack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_pack_real_array(buffer, position, val%aero_emission_rate)
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
         position - prev_position == pmc_mpi_pack_size_env_data(val))
#endif

  end subroutine pmc_mpi_pack_env_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_env_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(env_data_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%temp_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%temp)
    call pmc_mpi_unpack_real_array(buffer, position, val%height_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%height)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_emission_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%gas_dilution_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_emission_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_emission_rate)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_time)
    call pmc_mpi_unpack_real_array(buffer, position, val%aero_dilution_rate)
    allocate(val%gas_emission(size(val%gas_emission_time)))
    allocate(val%gas_background(size(val%gas_dilution_time)))
    allocate(val%aero_emission(size(val%aero_emission_time)))
    allocate(val%aero_background(size(val%aero_dilution_time)))
    do i = 1,size(val%gas_emission)
       call pmc_mpi_unpack_gas_state(buffer, position, val%gas_emission(i))
    end do
    do i = 1,size(val%gas_background)
       call pmc_mpi_unpack_gas_state(buffer, position, val%gas_background(i))
    end do
    do i = 1,size(val%aero_emission)
       call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_emission(i))
    end do
    do i = 1,size(val%aero_background)
       call pmc_mpi_unpack_aero_dist(buffer, position, val%aero_background(i))
    end do
    call assert(611542570, &
         position - prev_position == pmc_mpi_pack_size_env_data(val))
#endif

  end subroutine pmc_mpi_unpack_env_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_env_data
