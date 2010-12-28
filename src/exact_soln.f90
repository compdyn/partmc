! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_exact_soln module.

!> Exact solutions for various simulations.
module pmc_exact_soln

  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_aero_data
  use pmc_aero_dist
  use pmc_kernel
  use pmc_kernel_zero
  use pmc_kernel_additive
  use pmc_kernel_constant
  use pmc_env_state
  use pmc_env_data
  use pmc_aero_binned

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine exact_soln(bin_grid, aero_data, kernel_type, aero_dist_init, &
       env_data, env_state, time, aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
    !> Initial distribution.
    type(aero_dist_t), intent(in) :: aero_dist_init
    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Output state.
    type(aero_binned_t), intent(inout) :: aero_binned

    if (kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       ! FIXME: check env_data has no emissions or dilution
       if (aero_dist_init%n_mode /= 1) then
          call die_msg(827813758, "Exact solution with additive kernel " &
               // "requires exactly 1 initial distribution mode, not: " &
               // integer_to_string(aero_dist_init%n_mode))
       end if
       if (aero_dist_init%mode(1)%type /= AERO_MODE_TYPE_EXP) then
          call die_msg(574495367, "Exact solution with additive kernel " &
               // "requires exactly 1 initial distribution mode of " &
               // "exponential type, not: " &
               // aero_mode_type_to_string(aero_dist_init%mode(1)%type))
       end if
       call soln_additive_exp(bin_grid, aero_data, time, &
            aero_dist_init%mode(1)%num_conc, &
            aero_dist_init%mode(1)%mean_radius, env_state, aero_binned)
    elseif (kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       ! FIXME: check env_data has no emissions or dilution
       if (aero_dist_init%n_mode /= 1) then
          call die_msg(827813758, "Exact solution with constant kernel " &
               // "requires exactly 1 initial distribution mode, not: " &
               // integer_to_string(aero_dist_init%n_mode))
       end if
       if (aero_dist_init%mode(1)%type /= AERO_MODE_TYPE_EXP) then
          call die_msg(574495367, "Exact solution with constant kernel " &
               // "requires exactly 1 initial distribution mode of " &
               // "exponential type, not: " &
               // aero_mode_type_to_string(aero_dist_init%mode(1)%type))
       end if
       call soln_constant_exp(bin_grid, aero_data, time, &
            aero_dist_init%mode(1)%num_conc, &
            aero_dist_init%mode(1)%mean_radius, env_state, aero_binned)
    elseif (kernel_type == COAG_KERNEL_TYPE_ZERO) then
       ! FIXME: check env_data has constant emissions and constant dilution
       call soln_zero(bin_grid, aero_data, time, aero_dist_init, &
            env_state, aero_binned)
    else
       call die_msg(932981721, "No exact solutions with " &
            // "coagulation kernel type " // integer_to_string(kernel_type) &
            // "(" // kernel_type_to_string(kernel_type) // ")")
    end if

  end subroutine exact_soln

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_exact_soln
