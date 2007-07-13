! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Exact solution output.

module mod_run_exact

  use mod_inout
  use mod_aero_dist

  type run_exact_opt_t
     ! FIXME: following few items depend on kernel/soln choice
     real*8 :: num_conc                 ! particle number concentration (#/m^3)
     real*8 :: mean_vol                 ! mean init volume (m^3)
     type(aero_dist_t) :: aero_dist_init ! aerosol initial distribution
     real*8 :: rho_p                    ! particle density (kg/m^3)
     real*8 :: t_max                    ! total simulation time
     real*8 :: t_output                 ! interval to output info (s)
  end type run_exact_opt_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine run_exact(bin_grid, env, aero_data, exact_opt, soln, &
       summary_file)

    ! FIXME: num_conc and mean_vol are really parameters for the
    ! initial value of the particle distribution. They should be
    ! replaced by a n_param, params() pair.

    ! "Run" an exact solution, output data in the same format as
    ! particle-resolved or sectional simulations.
    
    use mod_bin_grid
    use mod_aero_state
    use mod_environ
    use mod_aero_data
    use mod_output_summary
    use mod_aero_binned
    use mod_gas_data
    use mod_gas_state
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(run_exact_opt_t), intent(in) :: exact_opt ! options
    type(inout_file_t), intent(inout) :: summary_file ! summary output file
    
    integer :: i_time, n_time
    type(aero_binned_t) :: aero_binned
    real*8 :: time
    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_state
    
    interface
       subroutine soln(bin_grid, aero_data, time, num_conc, mean_vol, &
            rho_p, aero_dist_init, env, aero_binned)

         use mod_bin_grid
         use mod_environ
         use mod_aero_binned
         use mod_aero_data

         type(bin_grid_t), intent(in) :: bin_grid ! bin grid
         type(aero_data_t), intent(in) :: aero_data ! aerosol data
         real*8, intent(in) :: time              ! current time
         real*8, intent(in) :: num_conc          ! particle number conc (#/m^3)
         real*8, intent(in) :: mean_vol          ! mean init volume (m^3)
         real*8, intent(in) :: rho_p             ! particle density (kg/m^3)
         type(aero_dist_t), intent(in) :: aero_dist_init ! initial distribution
         type(environ), intent(in) :: env        ! environment state
         type(aero_binned_t), intent(out) :: aero_binned ! output state
       end subroutine soln
    end interface

    call aero_binned_alloc(bin_grid%n_bin, aero_data%n_spec, aero_binned)
    call alloc_gas_data(0, gas_data)
    call gas_state_alloc(0, gas_state)

    n_time = nint(exact_opt%t_max / exact_opt%t_output)
    call init_environ(env, 0d0)
    do i_time = 0,n_time
       time = dble(i_time) / dble(n_time) * exact_opt%t_max
       call update_environ(env, time)
       call soln(bin_grid, aero_data, time, exact_opt%num_conc, &
            exact_opt%mean_vol, exact_opt%rho_p, &
            exact_opt%aero_dist_init, env, aero_binned)

       call output_summary(summary_file, time, bin_grid, aero_data, &
            aero_binned, gas_data, gas_state, env, 1)
    end do
    
  end subroutine run_exact
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_run_exact
