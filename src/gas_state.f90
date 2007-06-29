! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Gas state.

module mod_gas_state

  type gas_state_t
     real*8, pointer :: conc(:)          ! length n_spec, concentration (ppb)
  end type gas_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_gas_state(gas_data, gas_state)

    ! Allocate storage for gas species.

    use mod_gas_data

    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(out) :: gas_state ! gas state to be allocated

    allocate(gas_state%conc(gas_data%n_spec))

  end subroutine allocate_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copy_gas_state(from_state, to_state)

    ! Copy to an already allocated to_state.

    type(gas_state_t), intent(in) :: from_state ! existing gas state
    type(gas_state_t), intent(out) :: to_state ! must be allocated already

    integer :: n_spec

    n_spec = size(from_state%conc)
    deallocate(to_state%conc)
    allocate(to_state%conc(n_spec))
    to_state%conc = from_state%conc

  end subroutine copy_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_gas_state
