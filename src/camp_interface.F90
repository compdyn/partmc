! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_camp_interface module.

!> An interface between PartMC and the CAMP
module pmc_camp_interface

  use pmc_aero_data
  use pmc_aero_particle
  use pmc_aero_state
  use pmc_constants,                  only : dp
  use pmc_gas_data
  use pmc_gas_state
#ifdef PMC_USE_CAMP
  use camp_camp_core
  use camp_camp_state
  use pmc_photolysis
  use camp_rxn_data
  use camp_solver_stats
  use camp_util, only: split_string
#endif
  use pmc_util, only : die_msg, warn_assert_msg, assert_msg

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_CAMP
  !> Run the CAMP module for the current PartMC state
  subroutine pmc_camp_interface_solve(camp_core, camp_state, &
          camp_pre_aero_state, camp_post_aero_state, env_state, aero_data, &
          aero_state, gas_data, gas_state, photolysis, del_t)

    !> CAMP core
    type(camp_core_t), intent(in) :: camp_core
    !> CAMP state
    type(camp_state_t), intent(inout) :: camp_state
    !> Working CAMP state
    type(camp_state_t), intent(inout) :: camp_pre_aero_state
    !> Working CAMP state
    type(camp_state_t), intent(inout) :: camp_post_aero_state
    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state
    type(gas_state_t), intent(inout) :: gas_state
    !> Photolysis calculator
    type(photolysis_t), intent(inout) :: photolysis
    !> Time step (s)
    real(kind=dp), intent(in) :: del_t

    integer :: i_part
    real(kind=dp) :: num_conc
    integer :: camp_state_size
    type(solver_stats_t), target :: solver_stats
    type(camp_state_t) :: camp_state_pre_aero

    ! Set the CAMP environmental state.
    call env_state_set_camp_env_state(env_state, camp_state)

    ! Set the CAMP gas-phase species concentrations
    call gas_state_set_camp_conc(gas_state, camp_state, gas_data)

    ! Recalculate the photolysis rate constants
    call photolysis%update_rate_constants()

    ! Update the number concentrations and composition for all particles
    call pmc_camp_interface_set_camp_aerosol(aero_data, aero_state, &
                                             camp_core, camp_state)

    ! We're modifying particle diameters, so the bin sort is now invalid
    aero_state%valid_sort = .false.

    ! Solve the multi-phase chemical system
    call camp_core%solve(camp_state, del_t, solver_stats = solver_stats)

    call assert_msg(592148911, solver_stats%status_code == 0, &
         "Solver failed for aerosol-phase with code "// &
         integer_to_string(solver_stats%solver_flag))

    ! Update the PartMC aerosol state
    call pmc_camp_interface_set_partmc_aerosol(aero_data, aero_state, &
                                               camp_state)

    ! Update the PartMC gas-phase state
    call gas_state_get_camp_conc(gas_state, camp_state)

  end subroutine pmc_camp_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the CAMP aerosol-phase species and number concentrations
  subroutine pmc_camp_interface_set_camp_aerosol(aero_data, aero_state, &
      camp_core, camp_state)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> CAMP core.
    type(camp_core_t), intent(in) :: camp_core
    !> CAMP state.
    type(camp_state_t), intent(inout) :: camp_state

    real(kind=dp) :: num_conc
    integer :: i_part, i_spec

    associate (aero_rep => aero_data%aero_rep_ptr)
    select type (aero_rep)
       type is (aero_rep_single_particle_t)
          call assert_msg(858496327, aero_state_n_part(aero_state) <= &
               aero_rep%maximum_computational_particles(), &
               "Exceeded CAMP maximum number of particles. Number "// &
               "of computational particles: "// &
               trim(integer_to_string(aero_state_n_part(aero_state)))// &
               ". CAMP maximum: "// &
               trim(integer_to_string( &
                    aero_rep%maximum_computational_particles())))
          do i_part = 1, aero_state_n_part(aero_state)
             associate (part => aero_state%apa%particle(i_part))
             num_conc = aero_weight_array_num_conc(aero_state%awa, part, &
                  aero_data)
             call aero_state%update_number%set_number__n_m3(i_part, num_conc)
             call camp_core%update_data(aero_state%update_number)
             do i_spec = 1, size(aero_data%camp_particle_spec_id)
                camp_state%state_var(aero_data%camp_spec_id(i_part, i_spec)) &
                     = part%vol(i_spec) * aero_data%density(i_spec) ! kg m-3
             end do
             end associate
          end do
          do i_part = aero_state_n_part(aero_state) + 1, &
               aero_rep%maximum_computational_particles()
             call aero_state%update_number%set_number__n_m3(i_part, 0.0d0)
             call camp_core%update_data(aero_state%update_number)
             do i_spec = 1,size(aero_data%camp_particle_spec_id)
                camp_state%state_var(aero_data%camp_spec_id(i_part, i_spec)) &
                     = 0.0d0
             end do
          end do
      class default
         call die_msg(780366884, &
             "Wrong type for PartMC aerosol representation.")
    end select
    end associate

  end subroutine pmc_camp_interface_set_camp_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set PartMC aerosol-phase species concentrations
  subroutine pmc_camp_interface_set_partmc_aerosol(aero_data, aero_state, &
      camp_state)

    !> Aerosol particle
    type(aero_data_t), intent (in) :: aero_data
    !> Aerosol state
    type(aero_state_t), intent(inout) :: aero_state
    !> CAMP state
    type(camp_state_t), intent(in) :: camp_state

    integer :: i_part, i_spec

    associate (aero_rep => aero_data%aero_rep_ptr)
    select type (aero_rep)
      type is (aero_rep_single_particle_t)
        call assert_msg(464490945, aero_state_n_part(aero_state) <= &
             aero_rep%maximum_computational_particles(), &
             "Exceeded CAMP maximum number of particles")
        do i_part = 1, aero_state_n_part(aero_state)
          associate (part => aero_state%apa%particle(i_part))
          do i_spec = 1, size(aero_data%camp_particle_spec_id)
            part%vol(i_spec) = &
                 camp_state%state_var(aero_data%camp_spec_id(i_part, i_spec)) &
                 / aero_data%density(i_spec) ! m3
          end do
          end associate
        end do
      class default
        call die_msg(773649338, &
             "Wrong type for PartMC aerosol representation.")
    end select
    end associate

  end subroutine pmc_camp_interface_set_partmc_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
end module pmc_camp_interface
