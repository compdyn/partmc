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
  use pmc_constants,                  only : i_kind, dp
  use pmc_debug_diff_check,           only : diff_check, &
                                             diff_check_update_only
  use pmc_gas_data
  use pmc_gas_state
  use pmc_camp_core
  use pmc_camp_state
  use pmc_photolysis
  use pmc_rxn_data
  use pmc_solver_stats
  use pmc_util,                       only : die_msg, string_t, &
                                             warn_assert_msg, assert_msg

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the CAMP module for the current PartMC state
  subroutine pmc_camp_interface_solve(camp_core, camp_state, aero_data, &
          aero_state, gas_data, gas_state, photolysis, del_t)

    !> CAMP core
    type(camp_core_t), intent(in) :: camp_core
    !> CAMP state
    type(camp_state_t), intent(inout) :: camp_state
    !> Aerosol data
    type(aero_data_t), intent(inout) :: aero_data
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

    integer(kind=i_kind) :: i_part
    real(kind=dp) :: num_conc
    integer :: camp_state_size
    type(solver_stats_t), target :: solver_stats

    ! Set the camp chem  gas-phase species
    call gas_state%set_camp_conc(camp_state, gas_data)

    call diff_check_update_only( "pre rate update" )

    ! Recalculate the photolysis rate constants
    call photolysis%update_rate_constants()

    call diff_check( "post rate update" )

    ! Solve gas-phase chemistry
    call camp_core%solve(camp_state, del_t, GAS_RXN, &
                         solver_stats = solver_stats)

    call warn_assert_msg(828767454, solver_stats%status_code.eq.0, &
                         "Solver failed for gas-phase with code "// &
                         to_string(solver_stats%solver_flag))

    ! Do phase-transfer and aerosol-phase chemistry for each particle
    ! in the particle array
    do i_part = 1,aero_state%n_part()
      associate (part => aero_state%apa%particle(i_part))

      ! Set the CAMP chem aerosol state
      num_conc = aero_weight_array_num_conc(aero_state%awa, part, aero_data)
      call pmc_camp_interface_set_camp_conc(aero_data, part, camp_state, &
           num_conc)

      ! Update the number for this particle in CAMP
      call aero_data%update_number%set_number__n_m3( num_conc )
      call camp_core%update_data( aero_data%update_number )

      ! Solve the phase-transfer and aerosol-phase chemistry for this particle
      call camp_core%solve(camp_state, del_t, AERO_RXN)

      call assert_msg(592148911, solver_stats%status_code.eq.0, &
                     "Solver failed for aerosol-phase with code "// &
                     to_string(solver_stats%solver_flag))

      ! Update the PartMC aerosol state
      call pmc_camp_interface_get_camp_conc(aero_data, part, camp_state, &
           num_conc)

      end associate
    end do

    ! Update the PartMC gas-phase state
    call gas_state%get_camp_conc(camp_state)

  end subroutine pmc_camp_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the CAMP aerosol-phase species concentrations
  subroutine pmc_camp_interface_set_camp_conc(aero_data, aero_particle, &
            camp_state, num_conc)

    !> Aerosol particle
    type(aero_data_t), intent (in) :: aero_data
    !> Aerosol data
    type(aero_particle_t), intent(in) :: aero_particle
    !> CAMP state
    type(camp_state_t), intent(inout) :: camp_state
    !> Number concentration particle weighting
    real(kind=dp), intent(in) :: num_conc

    integer(kind=i_kind) :: i_spec

    ! Set the aerosol representation state in camp chem
    associate (aero_rep => aero_data%aero_rep_ptr)

    select type (aero_rep)
      type is (aero_rep_single_particle_t)
        do i_spec = 1, size(aero_data%camp_spec_id)
          camp_state%state_var(aero_data%camp_spec_id(i_spec)) = &
                  aero_particle%vol(i_spec) * &
                  aero_data%density(i_spec) ! kg m-3
        end do
      class default
        call die_msg(780366884, &
                "Wrong type for PartMC aerosol representation.")
    end select

    end associate

  end subroutine pmc_camp_interface_set_camp_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the CAMP aerosol-phase species concentrations
  subroutine pmc_camp_interface_get_camp_conc(aero_data, aero_particle, &
            camp_state, num_conc)

    !> Aerosol particle
    type(aero_data_t), intent (in) :: aero_data
    !> Aerosol data
    type(aero_particle_t), intent(inout) :: aero_particle
    !> CAMP state
    type(camp_state_t), intent(inout) :: camp_state
    !> Number concentration particle weighting
    real(kind=dp), intent(in) :: num_conc

    integer(kind=i_kind) :: i_spec

    ! Disassociate the aerosol representation state in camp chem
    associate (aero_rep => aero_data%aero_rep_ptr)

    select type (aero_rep)
      type is (aero_rep_single_particle_t)
        do i_spec = 1, size(aero_data%camp_spec_id)
          aero_particle%vol(i_spec) = &
              camp_state%state_var(aero_data%camp_spec_id(i_spec)) / &
              aero_data%density(i_spec) ! m3
        end do
      class default
        call die_msg(773649338, &
                "Wrong type for PartMC aerosol representation.")
    end select

    end associate

  end subroutine pmc_camp_interface_get_camp_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_camp_interface
