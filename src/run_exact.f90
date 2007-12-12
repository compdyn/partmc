! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Exact solution output.

module pmc_run_exact

  use pmc_inout
  use pmc_aero_dist
  use pmc_process_spec
  use pmc_bin_grid
  use pmc_aero_state
  use pmc_env_data
  use pmc_env
  use pmc_aero_data
  use pmc_output_processed
  use pmc_aero_binned
  use pmc_gas_data
  use pmc_gas_state

  type run_exact_opt_t
     ! FIXME: following few items depend on kernel/soln choice
     real*8 :: num_den                  ! particle number concentration (#/m^3)
     real*8 :: mean_radius              ! mean init radius (m)
     type(aero_dist_t) :: aero_dist_init ! aerosol initial distribution
     real*8 :: rho_p                    ! particle density (kg/m^3)
     real*8 :: t_max                    ! total simulation time
     real*8 :: t_output                 ! interval to output info (s)
     character(len=300) :: prefix       ! output prefix
  end type run_exact_opt_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine run_exact(bin_grid, env_data, env, aero_data, exact_opt, &
       soln, process_spec_list)

    ! FIXME: num_den and mean_radius are really parameters for the
    ! initial value of the particle distribution. They should be
    ! replaced by a n_param, params() pair.

    ! "Run" an exact solution, output data in the same format as
    ! particle-resolved or sectional simulations.
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_data_t), intent(in) :: env_data ! environment data
    type(env_t), intent(inout) :: env   ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(run_exact_opt_t), intent(in) :: exact_opt ! options
    type(process_spec_t), intent(in) :: process_spec_list(:) ! processing spec
    
    integer :: i_time, n_time, ncid
    type(aero_binned_t) :: aero_binned
    real*8 :: time
    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_state
    
    interface
       subroutine soln(bin_grid, aero_data, time, num_den, mean_radius, &
            rho_p, aero_dist_init, env, aero_binned)

         use pmc_bin_grid
         use pmc_env
         use pmc_aero_binned
         use pmc_aero_data

         type(bin_grid_t), intent(in) :: bin_grid ! bin grid
         type(aero_data_t), intent(in) :: aero_data ! aerosol data
         real*8, intent(in) :: time              ! current time
         real*8, intent(in) :: num_den           ! particle number conc (#/m^3)
         real*8, intent(in) :: mean_radius       ! mean init radius (m)
         real*8, intent(in) :: rho_p             ! particle density (kg/m^3)
         type(aero_dist_t), intent(in) :: aero_dist_init ! initial distribution
         type(env_t), intent(in) :: env          ! environment state
         type(aero_binned_t), intent(out) :: aero_binned ! output state
       end subroutine soln
    end interface

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call gas_data_alloc(gas_data, 0)
    call gas_state_alloc(gas_state, 0)
    call env_data_init_state(env_data, env, 0d0)

    call output_processed_open(exact_opt%prefix, 1, ncid)

    n_time = nint(exact_opt%t_max / exact_opt%t_output)
    do i_time = 0,n_time
       time = dble(i_time) / dble(n_time) * exact_opt%t_max
       call env_data_update_state(env_data, env, time)
       call soln(bin_grid, aero_data, time, exact_opt%num_den, &
            exact_opt%mean_radius, exact_opt%rho_p, &
            exact_opt%aero_dist_init, env, aero_binned)
       call output_processed_binned(ncid, process_spec_list, &
            bin_grid, aero_data, aero_binned, gas_data, gas_state, &
            env, i_time + 1, time, exact_opt%t_output)
    end do

    call output_processed_close(ncid)

    call gas_data_free(gas_data)
    call gas_state_free(gas_state)
    call aero_binned_free(aero_binned)
    
  end subroutine run_exact
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_run_exact
