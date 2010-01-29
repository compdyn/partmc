! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
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
  use pmc_aero_weight
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the kernel value with the given weight.
  subroutine weighted_kernel(kernel, aero_particle_1, aero_particle_2, &
       aero_data, aero_weight, env_state, k)

    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel.
    real(kind=dp), intent(out) :: k

#ifndef DOXYGEN_SKIP_DOC
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
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif
    
    real(kind=dp) :: unweighted_k
    real(kind=dp) :: radius_1, radius_2, radius_1_plus_2
    real(kind=dp) :: weight_1, weight_2, weight_1_plus_2

    call kernel(aero_particle_1, aero_particle_2, aero_data, &
            env_state, unweighted_k)
    radius_1 = aero_particle_radius(aero_particle_1)
    radius_2 = aero_particle_radius(aero_particle_2)
    radius_1_plus_2 = vol2rad(rad2vol(radius_1) + rad2vol(radius_2))
    weight_1 = aero_weight_value(aero_weight, radius_1)
    weight_2 = aero_weight_value(aero_weight, radius_2)
    weight_1_plus_2 = aero_weight_value(aero_weight, radius_1_plus_2)
    k = unweighted_k * weight_1 * weight_2 / weight_1_plus_2
    
  end subroutine weighted_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the max kernel value with the given weight.
  subroutine weighted_kernel_max(kernel_max, v1, v2, aero_data, &
       aero_weight, env_state, k_max)

    !> Volume of first particle.
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle.
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel maximum value.
    real(kind=dp), intent(out) :: k_max

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel_max(v1, v2, aero_data, env_state, k_max)
         use pmc_aero_data
         use pmc_env_state
         real(kind=dp), intent(in) :: v1
         real(kind=dp), intent(in) :: v2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k_max
       end subroutine kernel_max
    end interface
#endif

    real(kind=dp) :: unweighted_k_max, weight_1, weight_2, weight_1_plus_2

    call kernel_max(v1, v2, aero_data, env_state, unweighted_k_max)

    weight_1 = aero_weight_value(aero_weight, vol2rad(v1))
    weight_2 = aero_weight_value(aero_weight, vol2rad(v2))
    weight_1_plus_2 = aero_weight_value(aero_weight, vol2rad(v1 + v2))
    k_max = unweighted_k_max * weight_1 * weight_2 / weight_1_plus_2

  end subroutine weighted_kernel_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes an array of kernel values for each bin pair. k(i,j) is
  !> the kernel value at the centers of bins i and j. This assumes the
  !> kernel is only a function of the particle volumes.
  subroutine bin_kernel(n_bin, bin_v, aero_data, kernel, env_state, k)
    
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Volume of particles in bins (m^3).
    real(kind=dp), intent(in) :: bin_v(n_bin)
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Kernel values.
    real(kind=dp), intent(out) :: k(n_bin,n_bin)
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

#ifndef DOXYGEN_SKIP_DOC
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
         real(kind=dp), intent(out) :: k
       end subroutine kernel
    end interface
#endif
    
    integer :: i, j
    type(aero_particle_t) :: aero_particle_1, aero_particle_2
    
    call aero_particle_allocate_size(aero_particle_1, aero_data%n_spec)
    call aero_particle_allocate_size(aero_particle_2, aero_data%n_spec)
    do i = 1,n_bin
       do j = 1,n_bin
          aero_particle_1%vol(1) = bin_v(i)
          aero_particle_2%vol(1) = bin_v(j)
          call kernel(aero_particle_1, aero_particle_2, aero_data, &
               env_state, k(i,j))
       end do
    end do
    call aero_particle_deallocate(aero_particle_1)
    call aero_particle_deallocate(aero_particle_2)
    
  end subroutine bin_kernel
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Estimate an array of maximum kernel values. Given particles v1 in
  !> bin b1 and v2 in bin b2, it is probably true that kernel(v1,v2)
  !> <= k_max(b1,b2).
  subroutine est_k_max_binned(bin_grid, kernel_max, aero_data, &
       aero_weight, env_state, k_max)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Max kernel vals.
    real(kind=dp), intent(out) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    
#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel_max(v1, v2, aero_data, env_state, k_max)
         use pmc_aero_data
         use pmc_env_state
         real(kind=dp), intent(in) :: v1
         real(kind=dp), intent(in) :: v2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k_max
       end subroutine kernel_max
    end interface
#endif
    
    integer i, j
    
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          call est_k_max_for_bin(bin_grid, kernel_max, i, j, aero_data, &
               aero_weight, env_state, k_max(i,j))
       end do
    end do
    
  end subroutine est_k_max_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Samples within bins b1 and b2 to find the maximum value of the
  !> kernel between particles from the two bins.
  subroutine est_k_max_for_bin(bin_grid, kernel_max, b1, b2, aero_data, &
       aero_weight, env_state, k_max)
   
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> First bin.
    integer, intent(in) :: b1
    !> Second bin.
    integer, intent(in) :: b2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Maximum kernel values.
    real(kind=dp), intent(out) :: k_max
    
    !> Number of sample points per bin.
    integer, parameter :: n_sample = 3
    !> Over-estimation scale factor parameter.
    real(kind=dp), parameter :: over_scale = 1.1d0
    
#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine kernel_max(v1, v2, aero_data, env_state, k_max)
         use pmc_aero_data
         use pmc_env_state
         real(kind=dp), intent(in) :: v1
         real(kind=dp), intent(in) :: v2
         type(aero_data_t), intent(in) :: aero_data
         type(env_state_t), intent(in) :: env_state  
         real(kind=dp), intent(out) :: k_max
       end subroutine kernel_max
    end interface
#endif
    
    real(kind=dp) :: v1, v2, v1_high, v1_low, v2_high, v2_low, k
    integer :: i, j
    
    ! v1_low < bin_v(b1) < v1_high
    v1_low = bin_grid_edge(bin_grid, b1)
    v1_high = bin_grid_edge(bin_grid, b1 + 1)
    
    ! v2_low < bin_v(b2) < v2_high
    v2_low = bin_grid_edge(bin_grid, b2)
    v2_high = bin_grid_edge(bin_grid, b2 + 1)
    
    k_max = 0d0
    do i = 1,n_sample
       do j = 1,n_sample
          v1 = interp_linear_disc(v1_low, v1_high, n_sample, i)
          v2 = interp_linear_disc(v2_low, v2_high, n_sample, j)
          call weighted_kernel_max(kernel_max, v1, v2, aero_data, &
               aero_weight, env_state, k)
          if (k .gt. k_max) k_max = k
       end do
    end do
    
    k_max = k_max * over_scale
    
  end subroutine est_k_max_for_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_kernel
