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

!#ifdef PMC_USE_TCHEM
!    interface
!    end interface
!#endif
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_TCHEM
  !> Run the CAMP module for the current PartMC state
  subroutine pmc_tchem_interface_solve()

       call test()

  end subroutine pmc_tchem_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tchem_initialize(config_filename, gas_data, gas_state, aero_data)

    character(len=*), intent(in) :: config_filename
    type(gas_data_t), intent(inout) :: gas_data
    type(gas_state_t), intent(inout) :: gas_state
    type(aero_data_t), intent(inout) :: aero_data

    integer :: nSpec

    nSpec = 0

    call initialize(trim(config_filename))

    print*, 'in partmc', nSpec, trim(config_filename)

  end subroutine

  subroutine tchem_cleanup()

    call finalize_kokkos()

  end subroutine

#endif
end module pmc_tchem_interface
