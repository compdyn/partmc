! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_kernel module.

!> Generic coagulation kernel.
module pmc_kernel

  use pmc_env_state
  use pmc_bin_grid
  use pmc_aero_particle
  use pmc_aero_data
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes an array of kernel values for each bin pair. k(i,j) is
  !> the kernel value at the centers of bins i and j. This assumes the
  !> kernel is only a function of the particle volumes.
  subroutine bin_kernel(n_bin, bin_v, aero_data, kernel, env_state, k)
    
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Volume of particles in bins (m^3).
    real*8, intent(in) :: bin_v(n_bin)
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Kernel values.
    real*8, intent(out) :: k(n_bin,n_bin)
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    interface
       subroutine kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, k)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env_state
         type(aero_particle_t), intent(in) :: aero_particle_1
         type(aero_particle_t), intent(in) :: aero_particle_2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real*8, intent(out) :: k
       end subroutine kernel
    end interface
    
    integer :: i, j
    type(aero_particle_t) :: aero_particle_1, aero_particle_2
    
    call aero_particle_alloc(aero_particle_1, aero_data%n_spec)
    call aero_particle_alloc(aero_particle_2, aero_data%n_spec)
    do i = 1,n_bin
       do j = 1,n_bin
          aero_particle_1%vol(1) = bin_v(i)
          aero_particle_2%vol(1) = bin_v(j)
          call kernel(aero_particle_1, aero_particle_2, aero_data, &
               env_state, k(i,j))
       end do
    end do
    call aero_particle_free(aero_particle_1)
    call aero_particle_free(aero_particle_2)
    
  end subroutine bin_kernel
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Estimate an array of maximum kernel values. Given particles v1 in
  !> bin b1 and v2 in bin b2, it is probably true that kernel(v1,v2)
  !> <= k_max(b1,b2).
  subroutine est_k_max_binned(bin_grid, kernel_max, aero_data, &
       env_state, k_max)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Max kernel vals.
    real*8, intent(out) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    
    interface
       subroutine kernel_max(v1, v2, aero_data, env_state, k_max)
         use pmc_aero_data
         use pmc_env_state
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real*8, intent(out) :: k_max
       end subroutine kernel_max
    end interface
    
    integer i, j
    
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          call est_k_max_for_bin(bin_grid, kernel_max, i, j, aero_data, &
               env_state, k_max(i,j))
       end do
    end do
    
  end subroutine est_k_max_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Samples within bins b1 and b2 to find the maximum value of the
  !> kernel between particles from the two bins.
  subroutine est_k_max_for_bin(bin_grid, kernel_max, b1, b2, aero_data, &
       env_state, k_max)
   
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> First bin.
    integer, intent(in) :: b1
    !> Second bin.
    integer, intent(in) :: b2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Maximum kernel values.
    real*8, intent(out) :: k_max
    
    interface
       subroutine kernel_max(v1, v2, aero_data, env_state, k_max)
         use pmc_aero_data
         use pmc_env_state
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real*8, intent(out) :: k_max
       end subroutine kernel_max
    end interface
    
    real*8 v1, v2, v1_high, v1_low, v2_high, v2_low, k
    integer i, j
    
    !> Number of sample points per bin.
    integer, parameter :: n_sample = 10
    
    ! v1_low < bin_v(b1) < v1_high
    v1_low = bin_grid_edge(bin_grid, b1)
    v1_high = bin_grid_edge(bin_grid, b1 + 1)
    
    ! v2_low < bin_v(b2) < v2_high
    v2_low = bin_grid_edge(bin_grid, b2)
    v2_high = bin_grid_edge(bin_grid, b2 + 1)
    
    k_max = 0d0
    do i = 1,n_sample
       do j = 1,n_sample
          v1 = v1_high * dble(n_sample - i) / dble(n_sample - 1) + &
               v1_low * dble(i - 1) / dble(n_sample - 1)
          v2 = v2_high * dble(n_sample - j) / dble(n_sample - 1) + &
               v2_low * dble(j - 1) / dble(n_sample - 1)
          call kernel_max(v1, v2, aero_data, env_state, k)
          if (k .gt. k_max) k_max = k
       end do
    end do
    
  end subroutine est_k_max_for_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_kernel
