! Copyright (C) 2005-2016 Nicole Riemer and Matthew West
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
  use pmc_coag_kernel
  use pmc_coag_kernel_zero
  use pmc_coag_kernel_additive
  use pmc_coag_kernel_constant
  use pmc_env_state
  use pmc_scenario
  use pmc_aero_binned

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact_soln(bin_grid, aero_data, do_coagulation, &
       coag_kernel_type, aero_dist_init, scenario, env_state, time, &
       aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to do coagulation.
    logical, intent(in) :: do_coagulation
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Initial distribution.
    type(aero_dist_t), intent(in) :: aero_dist_init
    !> Environment data.
    type(scenario_t), intent(in) :: scenario
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Output state.
    type(aero_binned_t), intent(inout) :: aero_binned

    if (.not. do_coagulation) then
       call die_msg(287486666, 'Exact solutions require coagulation ' &
            // '(can set coag_kernel to "zero").')
    end if

    if (coag_kernel_type /= COAG_KERNEL_TYPE_ZERO &
       .and. scenario%loss_function_type /= SCENARIO_LOSS_FUNCTION_NONE) then
       call die_msg(189372109, 'Exact solution with particle loss ' &
            // 'requires using the "none" coag_kernel.')
    end if

    if (coag_kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       ! FIXME: check scenario has no emissions or dilution
       if (aero_dist_n_mode(aero_dist_init) /= 1) then
          call die_msg(285407619, "Exact solution with additive kernel " &
               // "requires exactly 1 initial distribution mode, not: " &
               // trim(integer_to_string(aero_dist_n_mode(aero_dist_init))))
       end if
       if (aero_dist_init%mode(1)%type /= AERO_MODE_TYPE_EXP) then
          call die_msg(373749499, "Exact solution with additive kernel " &
               // "requires exactly 1 initial distribution mode of " &
               // "exponential type, not: " &
               // aero_mode_type_to_string(aero_dist_init%mode(1)%type))
       end if
       call soln_additive_exp(bin_grid, aero_data, time, &
            aero_dist_init%mode(1)%num_conc, &
            aero_dist_init%mode(1)%char_radius, env_state, aero_binned)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       ! FIXME: check scenario has no emissions or dilution
       if (aero_dist_n_mode(aero_dist_init) /= 1) then
          call die_msg(827813758, "Exact solution with constant kernel " &
               // "requires exactly 1 initial distribution mode, not: " &
               // trim(integer_to_string(aero_dist_n_mode(aero_dist_init))))
       end if
       if (aero_dist_init%mode(1)%type /= AERO_MODE_TYPE_EXP) then
          call die_msg(574495367, "Exact solution with constant kernel " &
               // "requires exactly 1 initial distribution mode of " &
               // "exponential type, not: " &
               // aero_mode_type_to_string(aero_dist_init%mode(1)%type))
       end if
       call soln_constant_exp(bin_grid, aero_data, time, &
            aero_dist_init%mode(1)%num_conc, &
            aero_dist_init%mode(1)%char_radius, env_state, aero_binned)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_ZERO) then
       ! FIXME: check scenario has constant emissions and constant dilution
       call soln_zero(bin_grid, aero_data, time, aero_dist_init, &
            scenario, env_state, aero_binned)
    else
       call die_msg(932981721, "No exact solutions with " &
            // "coagulation kernel type " &
            // trim(integer_to_string(coag_kernel_type)) &
            // " (" // coag_kernel_type_to_string(coag_kernel_type) // ")")
    end if

  end subroutine exact_soln

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_exact_soln
