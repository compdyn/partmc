! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coag_kernel module.

!> Generic coagulation kernel.
module pmc_coag_kernel

  use pmc_env_state
  use pmc_bin_grid
  use pmc_aero_particle
  use pmc_aero_data
  use pmc_aero_weight
  use pmc_aero_weight_array
  use pmc_coag_kernel_sedi
  use pmc_coag_kernel_additive
  use pmc_coag_kernel_constant
  use pmc_coag_kernel_brown
  use pmc_coag_kernel_zero
  use pmc_coag_kernel_brown_free
  use pmc_coag_kernel_brown_cont

  !> Maximum length of a mode type.
  integer, parameter :: COAG_KERNEL_TYPE_LEN = 20

  !> Type code for an undefined or invalid kernel.
  integer, parameter :: COAG_KERNEL_TYPE_INVALID      = 0
  !> Type code for a sedimentation kernel.
  integer, parameter :: COAG_KERNEL_TYPE_SEDI         = 1
  !> Type code for an additive kernel.
  integer, parameter :: COAG_KERNEL_TYPE_ADDITIVE     = 2
  !> Type code for a constant kernel.
  integer, parameter :: COAG_KERNEL_TYPE_CONSTANT     = 3
  !> Type code for a Brownian kernel.
  integer, parameter :: COAG_KERNEL_TYPE_BROWN        = 4
  !> Type code for a zero kernel.
  integer, parameter :: COAG_KERNEL_TYPE_ZERO         = 5
  !> Type code for a Brownian kernel in free molecular regime from Vemury
  !> and Pratsinis [1995].
  integer, parameter :: COAG_KERNEL_TYPE_BROWN_FREE   = 6
  !> Type code for a Brownian kernel in continuum regime from Vemury and
  !> Pratsinis [1995].
  integer, parameter :: COAG_KERNEL_TYPE_BROWN_CONT   = 7

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a string representation of a kernel type.
  character(len=COAG_KERNEL_TYPE_LEN) function coag_kernel_type_to_string( &
       coag_kernel_type)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type

    if (coag_kernel_type == COAG_KERNEL_TYPE_INVALID) then
       coag_kernel_type_to_string = "invalid"
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_SEDI) then
       coag_kernel_type_to_string = "sedi"
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       coag_kernel_type_to_string = "additive"
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       coag_kernel_type_to_string = "constant"
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN) then
       coag_kernel_type_to_string = "brown"
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_ZERO) then
       coag_kernel_type_to_string = "zero"
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN_FREE) then
       coag_kernel_type_to_string = "brown_free"
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN_CONT) then
       coag_kernel_type_to_string = "brown_cont"
    else
       coag_kernel_type_to_string = "unknown"
    end if

  end function coag_kernel_type_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evalulate a coagulation kernel function.
  subroutine kernel(coag_kernel_type, aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
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

    if (coag_kernel_type == COAG_KERNEL_TYPE_SEDI) then
       call kernel_sedi(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       call kernel_additive(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       call kernel_constant(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN) then
       call kernel_brown(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_ZERO) then
       call kernel_zero(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN_FREE) then
       call kernel_brown_free(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN_CONT) then
       call kernel_brown_cont(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)
    else
       call die_msg(200724934, "Unknown kernel type: " &
            // trim(integer_to_string(coag_kernel_type)))
    end if

  end subroutine kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the minimum and maximum coagulation kernel.
  subroutine kernel_minmax(coag_kernel_type, v1, v2, aero_data, env_state, &
       k_min, k_max)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Minimum kernel value (m^3/s).
    real(kind=dp), intent(out) :: k_min
    !> Maximum kernel value (m^3/s).
    real(kind=dp), intent(out) :: k_max

    if (coag_kernel_type == COAG_KERNEL_TYPE_SEDI) then
       call kernel_sedi_minmax(v1, v2, aero_data, env_state, k_min, k_max)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
       call kernel_additive_minmax(v1, v2, aero_data, env_state, k_min, k_max)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_CONSTANT) then
       call kernel_constant_minmax(v1, v2, aero_data, env_state, k_min, k_max)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN) then
       call kernel_brown_minmax(v1, v2, aero_data, env_state, k_min, k_max)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_ZERO) then
       call kernel_zero_minmax(v1, v2, aero_data, env_state, k_min, k_max)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN_FREE) then
       call kernel_brown_free_minmax(v1, v2, aero_data, env_state, &
            k_min, k_max)
    elseif (coag_kernel_type == COAG_KERNEL_TYPE_BROWN_CONT) then
       call kernel_brown_cont_minmax(v1, v2, aero_data, env_state, &
            k_min, k_max)
    else
       call die_msg(330498208, "Unknown kernel type: " &
            // trim(integer_to_string(coag_kernel_type)))
    end if

  end subroutine kernel_minmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the kernel value with the given number concentration
  !> weighting.
  subroutine num_conc_weighted_kernel(coag_kernel_type, aero_particle_1, &
       aero_particle_2, i_class, j_class, ij_class, aero_data, &
       aero_weight_array, env_state, k)

    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Weight class of first particle.
    integer, intent(in) :: i_class
    !> Weight class of second particle.
    integer, intent(in) :: j_class
    !> Weight class of combined particle.
    integer, intent(in) :: ij_class
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel.
    real(kind=dp), intent(out) :: k

    real(kind=dp) :: unweighted_k, i_r, j_r

    call kernel(coag_kernel_type, aero_particle_1, aero_particle_2, &
         aero_data, env_state, unweighted_k)
    i_r = aero_particle_radius(aero_particle_1, aero_data)
    j_r = aero_particle_radius(aero_particle_2, aero_data)
    k = unweighted_k * coag_num_conc_factor(aero_weight_array, aero_data, &
         i_r, j_r, i_class, j_class, ij_class)

  end subroutine num_conc_weighted_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes an array of kernel values for each bin pair. k(i,j) is
  !> the kernel value at the centers of bins i and j. This assumes the
  !> kernel is only a function of the particle volumes.
  subroutine bin_kernel(n_bin, bin_r, aero_data, coag_kernel_type, &
       env_state, k)

    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Radii of particles in bins (m).
    real(kind=dp), intent(in) :: bin_r(n_bin)
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Kernel values.
    real(kind=dp), intent(out) :: k(n_bin,n_bin)

    integer :: i, j
    type(aero_particle_t) :: aero_particle_1, aero_particle_2

    do i = 1,n_bin
       do j = 1,n_bin
          call aero_particle_zero(aero_particle_1, aero_data)
          call aero_particle_zero(aero_particle_2, aero_data)
          aero_particle_1%vol(1) = aero_data_rad2vol(aero_data, bin_r(i))
          aero_particle_2%vol(1) = aero_data_rad2vol(aero_data, bin_r(j))
          call kernel(coag_kernel_type, aero_particle_1, aero_particle_2, &
               aero_data, env_state, k(i,j))
       end do
    end do

  end subroutine bin_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Estimate an array of minimum and maximum kernel values. Given
  !> particles v1 in bin b1 and v2 in bin b2, it is probably true that
  !> <tt>k_min(b1,b2) <= kernel(v1,v2) <= k_max(b1,b2)</tt>.
  subroutine est_k_minmax_binned_unweighted(bin_grid, coag_kernel_type, &
       aero_data, env_state, k_min, k_max)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Minimum kernel vals.
    real(kind=dp), intent(out) :: k_min(bin_grid_size(bin_grid), &
         bin_grid_size(bin_grid))
    !> Maximum kernel vals.
    real(kind=dp), intent(out) :: k_max(bin_grid_size(bin_grid), &
         bin_grid_size(bin_grid))

    integer i, j

    do i = 1,bin_grid_size(bin_grid)
       do j = 1,bin_grid_size(bin_grid)
          call est_k_minmax_for_bin_unweighted(bin_grid, coag_kernel_type, &
               i, j, aero_data, env_state, k_min(i,j), k_max(i,j))
       end do
    end do

  end subroutine est_k_minmax_binned_unweighted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Samples within bins b1 and b2 to find the minimum and maximum
  !> value of the kernel between particles from the two bins.
  subroutine est_k_minmax_for_bin_unweighted(bin_grid, coag_kernel_type, &
       b1, b2, aero_data, env_state, k_min, k_max)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Coagulation kernel type.
    integer, intent(in) :: coag_kernel_type
    !> First bin.
    integer, intent(in) :: b1
    !> Second bin.
    integer, intent(in) :: b2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Minimum kernel value.
    real(kind=dp), intent(out) :: k_min
    !> Maximum kernel value.
    real(kind=dp), intent(out) :: k_max

    !> Number of sample points per bin.
    integer, parameter :: n_sample = 3
    !> Over-estimation scale factor parameter.
    real(kind=dp), parameter :: over_scale = 2d0

    real(kind=dp) :: v1, v2, v1_high, v1_low, v2_high, v2_low
    real(kind=dp) :: new_k_min, new_k_max
    integer :: i, j

    ! v1_low < bin_v(b1) < v1_high
    v1_low = aero_data_rad2vol(aero_data, bin_grid%edges(b1))
    v1_high = aero_data_rad2vol(aero_data, bin_grid%edges(b1 + 1))

    ! v2_low < bin_v(b2) < v2_high
    v2_low = aero_data_rad2vol(aero_data, bin_grid%edges(b2))
    v2_high = aero_data_rad2vol(aero_data, bin_grid%edges(b2 + 1))

    do i = 1,n_sample
       do j = 1,n_sample
          v1 = interp_linear_disc(v1_low, v1_high, n_sample, i)
          v2 = interp_linear_disc(v2_low, v2_high, n_sample, j)
          call kernel_minmax(coag_kernel_type, v1, v2, aero_data, &
               env_state, new_k_min, new_k_max)
          if ((i == 1) .and. (j == 1)) then
             k_min = new_k_min
             k_max = new_k_max
          else
             k_min = min(k_min, new_k_min)
             k_max = max(k_max, new_k_max)
          end if
       end do
    end do

    k_max = k_max * over_scale

  end subroutine est_k_minmax_for_bin_unweighted

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Coagulation scale factor due to number concentrations.
  real(kind=dp) function coag_num_conc_factor(aero_weight_array, aero_data, &
       i_r, j_r, i_class, j_class, ij_class)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Radius of first particle.
    real(kind=dp), intent(in) :: i_r
    !> Radius of second particle.
    real(kind=dp), intent(in) :: j_r
    !> Weight class of first particle.
    integer, intent(in) :: i_class
    !> Weight class of second particle.
    integer, intent(in) :: j_class
    !> Weight class of combined particle.
    integer, intent(in) :: ij_class

    real(kind=dp) :: ij_r, i_nc, j_nc, ij_nc, nc_min

    ij_r = aero_data_vol2rad(aero_data, aero_data_rad2vol(aero_data, i_r) &
         + aero_data_rad2vol(aero_data, j_r))
    i_nc = aero_weight_array_num_conc_at_radius(aero_weight_array, i_class, &
         i_r)
    j_nc = aero_weight_array_num_conc_at_radius(aero_weight_array, j_class, &
         j_r)
    ij_nc = aero_weight_array_num_conc_at_radius(aero_weight_array, ij_class, &
         ij_r)
    nc_min = min(i_nc, j_nc, ij_nc)
    coag_num_conc_factor = i_nc * j_nc / nc_min

  end function coag_num_conc_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the weight class in which coagulated particles will be placed.
  integer function coag_dest_class(aero_weight_array, aero_data, &
       bin_grid, i_bin, j_bin, i_class, j_class)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> First bin number.
    integer, intent(in) :: i_bin
    !> Second bin number.
    integer, intent(in) :: j_bin
    !> Weight class of first particle.
    integer, intent(in) :: i_class
    !> Weight class of second particle.
    integer, intent(in) :: j_class

    real(kind=dp) :: i_r, j_r, ij_r, ij_nc_i, ij_nc_j

    i_r = bin_grid%centers(i_bin)
    j_r = bin_grid%centers(i_bin)
    ij_r = aero_data_vol2rad(aero_data, aero_data_rad2vol(aero_data, i_r) &
         + aero_data_rad2vol(aero_data, j_r))
    ij_nc_i = aero_weight_array_num_conc_at_radius(aero_weight_array, &
         i_class, ij_r)
    ij_nc_j = aero_weight_array_num_conc_at_radius(aero_weight_array, &
         j_class, ij_r)
    if (ij_nc_i < ij_nc_j) then
       coag_dest_class = i_class
    else
       coag_dest_class = j_class
    end if

  end function coag_dest_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the minimum and maximum number concentration factors
  !> for coagulation.
  subroutine max_coag_num_conc_factor(aero_weight_array, aero_data, &
       bin_grid, i_bin, j_bin, i_class, j_class, ij_class, f_max)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> First bin number.
    integer, intent(in) :: i_bin
    !> Second bin number.
    integer, intent(in) :: j_bin
    !> Weight class of first particle.
    integer, intent(in) :: i_class
    !> Weight class of second particle.
    integer, intent(in) :: j_class
    !> Weight class of coagulated particle.
    integer, intent(in) :: ij_class
    !> Maximum coagulation factor.
    real(kind=dp), intent(out) :: f_max

    integer, parameter :: n_sample = 5

    real(kind=dp) :: i_r_min, i_r_max, j_r_min, j_r_max, i_r, j_r, f
    integer :: i_sample, j_sample

    i_r_min = bin_grid%edges(i_bin)
    i_r_max = bin_grid%edges(i_bin + 1)
    j_r_min = bin_grid%edges(j_bin)
    j_r_max = bin_grid%edges(j_bin + 1)

    f_max = 0d0
    do i_sample = 1,n_sample
       do j_sample = 1,n_sample
          i_r = interp_linear_disc(i_r_min, i_r_max, n_sample, i_sample)
          j_r = interp_linear_disc(j_r_min, j_r_max, n_sample, j_sample)
          f = coag_num_conc_factor(aero_weight_array, aero_data, i_r, j_r, &
               i_class, j_class, ij_class)
          f_max = max(f_max, f)
       end do
    end do

  end subroutine max_coag_num_conc_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a kernel type from a spec file and
  !> generate it.
  subroutine spec_file_read_coag_kernel_type(file, coag_kernel_type)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Kernel type.
    integer, intent(out) :: coag_kernel_type

    character(len=SPEC_LINE_MAX_VAR_LEN) :: kernel_name

    !> \page input_format_coag_kernel Input File Format: Coagulation Kernel
    !!
    !! The coagulation kernel is specified by the parameter:
    !!   - \b coag_kernel (string): the type of coagulation kernel ---
    !!     must be one of: \c sedi for the gravitational sedimentation
    !!     kernel; \c additive for the additive kernel; \c constant
    !!     for the constant kernel; \c brown for the Brownian kernel,
    !!     or \c zero for no coagulation
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'coag_kernel', kernel_name)
    if (trim(kernel_name) == 'sedi') then
       coag_kernel_type = COAG_KERNEL_TYPE_SEDI
    elseif (trim(kernel_name) == 'additive') then
       coag_kernel_type = COAG_KERNEL_TYPE_ADDITIVE
    elseif (trim(kernel_name) == 'constant') then
       coag_kernel_type = COAG_KERNEL_TYPE_CONSTANT
    elseif (trim(kernel_name) == 'brown') then
       coag_kernel_type = COAG_KERNEL_TYPE_BROWN
    elseif (trim(kernel_name) == 'zero') then
       coag_kernel_type = COAG_KERNEL_TYPE_ZERO
    elseif (trim(kernel_name) == 'brown_free') then
       coag_kernel_type = COAG_KERNEL_TYPE_BROWN_FREE
    elseif (trim(kernel_name) == 'brown_cont') then
       coag_kernel_type = COAG_KERNEL_TYPE_BROWN_CONT
    else
       call spec_file_die_msg(920761229, file, &
            "Unknown coagulation kernel type: " // trim(kernel_name))
    end if

  end subroutine spec_file_read_coag_kernel_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coag_kernel
