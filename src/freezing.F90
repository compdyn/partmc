! Copyright (C) 2021 University of Illinois at Urbana-Champaign
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_freezing module.

module pmc_freezing
    use pmc_aero_state
    use pmc_env_state
    use pmc_aero_data
    use pmc_util
    use pmc_aero_particle
    use pmc_constants
    use pmc_rand

contains

    subroutine freeze(aero_state, aero_data, env_state_initial, &
        env_state_final, del_t)

        !> Aerosol state.
        type(aero_state_t), intent(inout) :: aero_state
        !> Aerosol data.
        type(aero_data_t), intent(in) :: aero_data
        !> Environment state at the start of the timestep.
        type(env_state_t), intent(in) :: env_state_initial
        !> Environment state at the end of the timestep. The rel_humid
        !> value will be ignored and overwritten with a new value.
        type(env_state_t), intent(inout) :: env_state_final
        !> Total time to integrate.
        real(kind=dp), intent(in) :: del_t

        integer :: i_part
        real(kind=dp) :: tmp

        do i_part = 1, aero_state_n_part(aero_state)
            tmp = pmc_random()
            !print *, 'i: ', i_part, tmp
        end do

    end subroutine freeze

end module pmc_freezing
