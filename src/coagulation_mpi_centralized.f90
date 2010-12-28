! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coagulation_mpi_centralized module.

!> Aerosol particle coagulation with MPI where each node has its own
!> aero_state, but all particles are transmitted to node 0, which does
!> coagulation, and then sends the particles back to the individual
!> nodes.
module pmc_coagulation_mpi_centralized

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_env_state
  use pmc_aero_state
  use pmc_aero_weight
  use pmc_coagulation
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do coagulation for time del_t in parallel by centralizing on node 0.
  subroutine mc_coag_mpi_centralized(coag_kernel_type, bin_grid, env_state, &
       aero_data, aero_weight, aero_state, del_t, k_max, tot_n_samp, &
       tot_n_coag)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Timestep for coagulation.
    real(kind=dp) :: del_t
    !> Maximum kernel.
    real(kind=dp), intent(in) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    !> Total number of samples tested.
    integer, intent(out) :: tot_n_samp
    !> Number of coagulation events.
    integer, intent(out) :: tot_n_coag

    type(aero_state_t) :: aero_state_total
    
#ifdef PMC_USE_MPI

    call aero_state_allocate(aero_state_total)
    call aero_state_mpi_gather(aero_state, aero_state_total)
    if (pmc_mpi_rank() == 0) then
       call mc_coag(coag_kernel_type, bin_grid, env_state, aero_data, &
            aero_weight, aero_state_total, del_t, k_max, tot_n_samp, &
            tot_n_coag)
    end if
    call aero_state_mpi_scatter(aero_state_total, aero_state, aero_data)
    call aero_state_deallocate(aero_state_total)
    
#endif

  end subroutine mc_coag_mpi_centralized
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coagulation_mpi_centralized
