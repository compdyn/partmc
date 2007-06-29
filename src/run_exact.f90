! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Exact solution output.

module mod_run_exact
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine run_exact(n_bin, n_spec, bin_v, bin_g_den, bin_gs_den, &
       bin_n_den, num_conc, mean_vol, rho_p, soln, t_max, t_output, &
       output_unit, env, aero_data)

    ! FIXME: num_conc and mean_vol are really parameters for the
    ! initial value of the particle distribution. They should be
    ! replaced by a n_param, params() pair.

    ! "Run" an exact solution, output data in the same format as
    ! particle-resolved or sectional simulations.
    
    use mod_bin
    use mod_aero_state
    use mod_environ
    use mod_aero_data
    use mod_output
    
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: bin_v(n_bin)  ! volume of bins
    real*8, intent(out) :: bin_g_den(n_bin) ! volume density in bins
    real*8, intent(out) :: bin_n_den(n_bin) ! number density in bins
    real*8, intent(out) :: bin_gs_den(n_bin,n_spec) ! volume density by species
    real*8, intent(in) :: num_conc      ! particle number concentration (#/m^3)
    real*8, intent(in) :: mean_vol      ! mean init volume (m^3)
    real*8, intent(in) :: rho_p         ! particle density (kg/m^3)
    real*8, intent(in) :: t_max         ! total simulation time
    real*8, intent(in) :: t_output      ! interval to output info (seconds)
    integer, intent(in) :: output_unit  ! unit number to output to
    type(environ), intent(inout) :: env ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    
    integer i_time, n_time
    real*8 time
    
    interface
       subroutine soln(n_bin, bin_v, bin_g_den, bin_n_den, &
            time, num_conc, mean_vol, rho_p, env)

         use mod_environ

         integer, intent(in) :: n_bin            ! number of bins
         real*8, intent(in) :: bin_v(n_bin)      ! volume of particles in bins
         real*8, intent(out) :: bin_g_den(n_bin) ! volume density in bins
         real*8, intent(out) :: bin_n_den(n_bin) ! number density in bins
         
         real*8, intent(in) :: time              ! current time
         real*8, intent(in) :: num_conc          ! particle number conc (#/m^3)
         real*8, intent(in) :: mean_vol          ! mean init volume (m^3)
         real*8, intent(in) :: rho_p             ! particle density (kg/m^3)
         type(environ), intent(in) :: env        ! environment state
       end subroutine soln
    end interface
    
    n_time = nint(t_max / t_output)
    call init_environ(env, 0d0)
    do i_time = 0,n_time
       time = dble(i_time) / dble(n_time) * t_max
       call update_environ(env, time)
       call soln(n_bin, bin_v, bin_g_den, bin_n_den, time, &
            num_conc, mean_vol, rho_p, env)
       bin_gs_den(:,1) = bin_g_den

       call output_info_density(output_unit, time, n_bin, n_spec, &
            bin_v, bin_g_den, bin_gs_den, bin_n_den, env, aero_data, 1)
    end do
    
  end subroutine run_exact
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_run_exact
