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
  use pmc_kernel_sedi
  use pmc_kernel_additive
  use pmc_kernel_constant
  use pmc_kernel_brown
  use pmc_kernel_zero

  !> Maximum length of a mode type.
  integer, parameter :: COAG_KERNEL_TYPE_LEN = 20

  !> Type code for an undefined or invalid kernel.
  integer, parameter :: COAG_KERNEL_TYPE_INVALID  = 0
  !> Type code for a sedimentation kernel.
  integer, parameter :: COAG_KERNEL_TYPE_SEDI     = 1
  !> Type code for an additive kernel.
  integer, parameter :: COAG_KERNEL_TYPE_ADDITIVE  = 2
  !> Type code for a constant kernel.
  integer, parameter :: COAG_KERNEL_TYPE_CONSTANT = 3
  !> Type code for a Brownian kernel.
  integer, parameter :: COAG_KERNEL_TYPE_BROWN    = 4
  !> Type code for a zero kernel.
  integer, parameter :: COAG_KERNEL_TYPE_ZERO     = 5
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a string representation of a kernel type.
  character(len=COAG_KERNEL_TYPE_LEN) function kernel_type_to_string( &
       kernel_type)

    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
   
    if (kernel_type == COAG_KERNEL_TYPE_INVALID) then
       kernel_type_to_string = "invalid"
    elseif (kernel_type == COAG_KERNEL_TYPE_SEDI) then
       kernel_type_to_string = "sedi"
    elseif (kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       kernel_type_to_string = "additive"
    elseif (kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       kernel_type_to_string = "constant"
    elseif (kernel_type == COAG_KERNEL_TYPE_BROWN) then
       kernel_type_to_string = "brown"
    elseif (kernel_type == COAG_KERNEL_TYPE_ZERO) then
       kernel_type_to_string = "zero"
    else
       kernel_type_to_string = "unknown"
    end if

  end function kernel_type_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evalulate a coagulation kernel function.
  subroutine kernel(kernel_type, aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)

    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Kernel k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: k

    if (kernel_type == COAG_KERNEL_TYPE_SEDI) then
       call kernel_sedi(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       call kernel_additive(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       call kernel_constant(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (kernel_type == COAG_KERNEL_TYPE_BROWN) then
       call kernel_brown(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (kernel_type == COAG_KERNEL_TYPE_ZERO) then
       call kernel_zero(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    else
       call die_msg(200724934, "Unknown kernel type: " &
            // integer_to_string(kernel_type))
    end if

  end subroutine kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the maximum coagulation kernel.
  subroutine kernel_max(kernel_type, v1, v2, aero_data, env_state, k_max)

    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Maximum kernel value (m^3/s).
    real(kind=dp), intent(out) :: k_max

    if (kernel_type == COAG_KERNEL_TYPE_SEDI) then
       call kernel_sedi_max(v1, v2, aero_data, env_state, k_max)
    elseif (kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       call kernel_additive_max(v1, v2, aero_data, env_state, k_max)
    elseif (kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       call kernel_constant_max(v1, v2, aero_data, env_state, k_max)
    elseif (kernel_type == COAG_KERNEL_TYPE_BROWN) then
       call kernel_brown_max(v1, v2, aero_data, env_state, k_max)
    elseif (kernel_type == COAG_KERNEL_TYPE_ZERO) then
       call kernel_zero_max(v1, v2, aero_data, env_state, k_max)
    else
       call die_msg(330498208, "Unknown kernel type: " &
            // integer_to_string(kernel_type))
    end if

  end subroutine kernel_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the kernel value with the given weight.
  subroutine weighted_kernel(kernel_type, aero_particle_1, aero_particle_2, &
       aero_data, aero_weight, env_state, k)

    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
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

    real(kind=dp) :: unweighted_k
    real(kind=dp) :: radius_1, radius_2, radius_1_plus_2
    real(kind=dp) :: weight_1, weight_2, weight_1_plus_2, weight_min

    call kernel(kernel_type, aero_particle_1, aero_particle_2, aero_data, &
            env_state, unweighted_k)
    radius_1 = aero_particle_radius(aero_particle_1)
    radius_2 = aero_particle_radius(aero_particle_2)
    radius_1_plus_2 = vol2rad(rad2vol(radius_1) + rad2vol(radius_2))
    weight_1 = aero_weight_value(aero_weight, radius_1)
    weight_2 = aero_weight_value(aero_weight, radius_2)
    weight_1_plus_2 = aero_weight_value(aero_weight, radius_1_plus_2)
    weight_min = min(weight_1, weight_2, weight_1_plus_2)
    k = unweighted_k * weight_1 * weight_2 / weight_min
    
  end subroutine weighted_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the max kernel value with the given weight.
  subroutine weighted_kernel_max(kernel_type, v1, v2, aero_data, &
       aero_weight, env_state, k_max)

    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
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

    real(kind=dp) :: unweighted_k_max, weight_1, weight_2, weight_1_plus_2
    real(kind=dp) :: weight_min

    call kernel_max(kernel_type, v1, v2, aero_data, env_state, &
         unweighted_k_max)

    weight_1 = aero_weight_value(aero_weight, vol2rad(v1))
    weight_2 = aero_weight_value(aero_weight, vol2rad(v2))
    weight_1_plus_2 = aero_weight_value(aero_weight, vol2rad(v1 + v2))
    weight_min = min(weight_1, weight_2, weight_1_plus_2)
    k_max = unweighted_k_max * weight_1 * weight_2 / weight_min

  end subroutine weighted_kernel_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes an array of kernel values for each bin pair. k(i,j) is
  !> the kernel value at the centers of bins i and j. This assumes the
  !> kernel is only a function of the particle volumes.
  subroutine bin_kernel(n_bin, bin_r, aero_data, kernel_type, env_state, k)
    
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Radii of particles in bins (m).
    real(kind=dp), intent(in) :: bin_r(n_bin)
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Kernel values.
    real(kind=dp), intent(out) :: k(n_bin,n_bin)

    integer :: i, j
    type(aero_particle_t) :: aero_particle_1, aero_particle_2
    
    call aero_particle_allocate_size(aero_particle_1, aero_data%n_spec)
    call aero_particle_allocate_size(aero_particle_2, aero_data%n_spec)
    do i = 1,n_bin
       do j = 1,n_bin
          aero_particle_1%vol(1) = rad2vol(bin_r(i))
          aero_particle_2%vol(1) = rad2vol(bin_r(j))
          call kernel(kernel_type, aero_particle_1, aero_particle_2, &
               aero_data, env_state, k(i,j))
       end do
    end do
    call aero_particle_deallocate(aero_particle_1)
    call aero_particle_deallocate(aero_particle_2)
    
  end subroutine bin_kernel
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Estimate an array of maximum kernel values. Given particles v1 in
  !> bin b1 and v2 in bin b2, it is probably true that kernel(v1,v2)
  !> <= k_max(b1,b2).
  subroutine est_k_max_binned(bin_grid, kernel_type, aero_data, &
       aero_weight, env_state, k_max)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Max kernel vals.
    real(kind=dp), intent(out) :: k_max(bin_grid%n_bin,bin_grid%n_bin)
    
    integer i, j
    
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          call est_k_max_for_bin(bin_grid, kernel_type, i, j, aero_data, &
               aero_weight, env_state, k_max(i,j))
       end do
    end do
    
  end subroutine est_k_max_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Samples within bins b1 and b2 to find the maximum value of the
  !> kernel between particles from the two bins.
  subroutine est_k_max_for_bin(bin_grid, kernel_type, b1, b2, aero_data, &
       aero_weight, env_state, k_max)
   
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Coagulation kernel type.
    integer, intent(in) :: kernel_type
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
    
    real(kind=dp) :: v1, v2, v1_high, v1_low, v2_high, v2_low, k
    integer :: i, j
    
    ! v1_low < bin_v(b1) < v1_high
    v1_low = rad2vol(bin_grid%edge_radius(b1))
    v1_high = rad2vol(bin_grid%edge_radius(b1 + 1))
    
    ! v2_low < bin_v(b2) < v2_high
    v2_low = rad2vol(bin_grid%edge_radius(b2))
    v2_high = rad2vol(bin_grid%edge_radius(b2 + 1))
    
    k_max = 0d0
    do i = 1,n_sample
       do j = 1,n_sample
          v1 = interp_linear_disc(v1_low, v1_high, n_sample, i)
          v2 = interp_linear_disc(v2_low, v2_high, n_sample, j)
          call weighted_kernel_max(kernel_type, v1, v2, aero_data, &
               aero_weight, env_state, k)
          if (k .gt. k_max) k_max = k
       end do
    end do
    
    k_max = k_max * over_scale
    
  end subroutine est_k_max_for_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a kernel type from a spec file and
  !> generate it.
  subroutine spec_file_read_kernel_type(file, kernel_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Kernel type.
    integer, intent(out) :: kernel_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: kernel_name

    !> \page input_format_coag_kernel Input File Format: Coagulation Kernel
    !!
    !! The coagulation kernel is specified by the parameter:
    !!   - \b kernel (string): the type of coagulation kernel --- must
    !!     be one of: \c sedi for the gravitational sedimentation
    !!     kernel; \c additive for the additive kernel;
    !!     \c constant for the constant kernel; \c brown for the
    !!     Brownian kernel, or \c zero for no coagulation
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'kernel', kernel_name)
    if (trim(kernel_name) == 'sedi') then
       kernel_type = COAG_KERNEL_TYPE_SEDI
    elseif (trim(kernel_name) == 'additive') then
       kernel_type = COAG_KERNEL_TYPE_ADDITIVE
    elseif (trim(kernel_name) == 'constant') then
       kernel_type = COAG_KERNEL_TYPE_CONSTANT
    elseif (trim(kernel_name) == 'brown') then
       kernel_type = COAG_KERNEL_TYPE_BROWN
    elseif (trim(kernel_name) == 'zero') then
       kernel_type = COAG_KERNEL_TYPE_ZERO
    else
       call spec_file_die_msg(494684716, file, &
            "Unknown kernel type: " // trim(kernel_name))
    end if

  end subroutine spec_file_read_kernel_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_kernel
