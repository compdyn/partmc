! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_camp_interface module.

!> An interface between PartMC and the CAMP
module pmc_tchem_interface

  use pmc_aero_data
  use pmc_aero_particle
  use pmc_aero_state
  use pmc_constants,                  only : dp
  use pmc_gas_data
  use pmc_gas_state
#ifdef PMC_USE_TCHEM
  use iso_c_binding
  use tchemdriver
#endif
  use pmc_util, only : die_msg, warn_assert_msg, assert_msg

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_TCHEM
  !> Run the CAMP module for the current PartMC state
  subroutine pmc_tchem_interface_solve(env_state, aero_data, aero_state, &
       gas_data, gas_state)

    type(env_state_t), intent(in) :: env_state
    type(aero_data_t), intent(in) :: aero_data
    type(aero_state_t), intent(inout) :: aero_state
    type(gas_data_t), intent(in) :: gas_data
    type(gas_state_t), intent(inout) :: gas_state

    call tchem_from_partmc(gas_data, gas_state)
    call timestep()
    call tchem_to_partmc(gas_data, gas_state)

  end subroutine pmc_tchem_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_tchem_initialize(config_filename, gas_data, gas_state, aero_data)

    use iso_c_binding

    character(len=*), intent(in) :: config_filename
    type(gas_data_t), intent(inout) :: gas_data
    type(gas_state_t), intent(inout) :: gas_state
    type(aero_data_t), intent(inout) :: aero_data
    integer(c_int) :: nSpec
    integer :: i, j

    character(len=5) :: name
    character(kind=c_char), dimension(5) :: name_data
    real(kind=c_double), dimension(:), allocatable :: array 

    ! initialize the model
    call TChem_initialize(trim(config_filename))

    ! feedback to PartMC
    ! what we need:
    !  - number of gas species
    !  - name of gas species to set gas_data
    !  - set the size of gas_state
    nSpec = TChem_getNumberOfSpecies() 
    call ensure_string_array_size(gas_data%name, nSpec)

    call gas_state_set_size(gas_state, nSpec)

    ! name of gas species
    !do i = 1,nSpec
    !   gas_data%name(i) = "H2O"
    !end do   
 
    print*, 'in partmc', nSpec, trim(config_filename), gas_data_n_spec(gas_data)

    call tchem_to_partmc(gas_data, gas_state)
    print*, 'in partmc'
    print*, gas_state%mix_rat 
  
end subroutine

  subroutine pmc_tchem_cleanup()

    call finalize()

  end subroutine pmc_tchem_cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine tchem_to_partmc(gas_data, gas_state)

    !>
    type(gas_data_t), intent(in) :: gas_data
    !>
    type(gas_state_t), intent(inout) :: gas_state

    integer(c_int) :: nSpec
    real(kind=c_double), dimension(:), allocatable :: array

    ! Get gas array
    nSpec = TChem_getNumberOfSpecies()
    allocate(array(nSpec))
    array = 0.0d0
    call TChem_getStateVector(array)
    gas_state%mix_rat = array

  end subroutine tchem_to_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine tchem_from_partmc(gas_data, gas_state)

    !>
    type(gas_data_t), intent(in) :: gas_data
    ! FIXME: Fix intent later
    type(gas_state_t), intent(inout) :: gas_state

    call TChem_setStateVector(gas_state%mix_rat)

  end subroutine tchem_from_partmc
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_tchem_interface
