! Copyright (C) 2012 Jeff Curtis
! Copyright (C) 2012-2015 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_weight_array module.

!> The aero_weight_array_t structure and associated subroutines.
module pmc_aero_weight_array

  use pmc_util
  use pmc_constants
  use pmc_rand
  use pmc_spec_file
  use pmc_aero_particle
  use pmc_netcdf
  use pmc_mpi
  use pmc_aero_weight
  use pmc_aero_data
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> An array of aerosol size distribution weighting functions.
  !!
  !! Each particle has a weight group and a weight class, which
  !! defines the weighting function associated with it. To determine
  !! the actual weight for a particle, the harmonic mean of the weight
  !! from each group is computed (for the given class). All particles
  !! of a given size within a weight class will thus have the same
  !! weight, irrespective of which group they are in. The group is
  !! only important when doubling or halving occurs, which takes place
  !! for a specific group and class.
  type aero_weight_array_t
     !> Aero weight array, <tt>aero_weight_array_n_group(aero_weight_array) x
     !> aero_weight_array_n_class(aero_weight_array)</tt>.
     type(aero_weight_t), allocatable :: weight(:, :)
  end type aero_weight_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the number of weight groups and classes.
  subroutine aero_weight_array_set_sizes(aero_weight_array, n_group, &
       n_class)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array
    !> Number of weight groups.
    integer,intent(in) :: n_group
    !> Number of weight classes.
    integer,intent(in) :: n_class

    if (allocated(aero_weight_array%weight)) then
       deallocate(aero_weight_array%weight)
    end if
    allocate(aero_weight_array%weight(n_group, n_class))

  end subroutine aero_weight_array_set_sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an \c aero_weight_array as flat weightings.
  subroutine aero_weight_array_set_flat(aero_weight_array, n_class)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(out) :: aero_weight_array
    !> Number of weight classes.
    integer,intent(in) :: n_class

    call aero_weight_array_set_sizes(aero_weight_array, 1, n_class)
    aero_weight_array%weight%type = AERO_WEIGHT_TYPE_NONE
    aero_weight_array%weight%magnitude = 1d0
    aero_weight_array%weight%exponent = 0d0

  end subroutine aero_weight_array_set_flat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an \c aero_weight_array as power weightings.
  subroutine aero_weight_array_set_power(aero_weight_array, n_class, &
       exponent)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(out) :: aero_weight_array
    !> Number of weight classes.
    integer,intent(in) :: n_class
    !> Exponent for power-law.
    real(kind=dp), intent(in) :: exponent

    call aero_weight_array_set_sizes(aero_weight_array, 1, n_class)
    aero_weight_array%weight%type = AERO_WEIGHT_TYPE_POWER
    aero_weight_array%weight%magnitude = 1d0
    aero_weight_array%weight%exponent = exponent

  end subroutine aero_weight_array_set_power

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an \c aero_weight_array as joint flat/power-3 weightings..
  subroutine aero_weight_array_set_nummass(aero_weight_array, n_class)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(out) :: aero_weight_array
    !> Number of weight classes.
    integer,intent(in) :: n_class

    call aero_weight_array_set_sizes(aero_weight_array, 2, n_class)
    aero_weight_array%weight(1, :)%type = AERO_WEIGHT_TYPE_NONE
    aero_weight_array%weight(1, :)%magnitude = 1d0
    aero_weight_array%weight(1, :)%exponent = 0d0
    aero_weight_array%weight(2, :)%type = AERO_WEIGHT_TYPE_POWER
    aero_weight_array%weight(2, :)%magnitude = 1d0
    aero_weight_array%weight(2, :)%exponent = -3d0

  end subroutine aero_weight_array_set_nummass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Normalizes the \c aero_weight_array to a non-zero value.
  elemental subroutine aero_weight_array_normalize(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array

    call aero_weight_normalize(aero_weight_array%weight)

  end subroutine aero_weight_array_normalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the number of weight groups.
  integer function aero_weight_array_n_group(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array

    aero_weight_array_n_group = size(aero_weight_array%weight, 1)

  end function aero_weight_array_n_group

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the number of weight classes.
  integer function aero_weight_array_n_class(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array

    aero_weight_array_n_class = size(aero_weight_array%weight, 2)

  end function aero_weight_array_n_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale the weights by the given factor, so
  !> <tt>new_weight = old_weight * factor</tt>.
  subroutine aero_weight_array_scale(aero_weight_array, factor)

    !> Aerosol weight array to scale.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array
    !> Factor to scale by.
    real(kind=dp), intent(in) :: factor

    call aero_weight_scale(aero_weight_array%weight, factor)

  end subroutine aero_weight_array_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Combine \c aero_weight_array_delta into \c aero_weight_array
  !> with a harmonic mean.
  subroutine aero_weight_array_combine(aero_weight_array, &
      aero_weight_array_delta)

    !> Aerosol weight array to combine into.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array
    !> Aerosol weight array to combine from.
    type(aero_weight_array_t), intent(in) :: aero_weight_array_delta

    call aero_weight_combine(aero_weight_array%weight, &
            aero_weight_array_delta%weight)

  end subroutine aero_weight_array_combine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adjust source and destination weights to reflect moving \c
  !> sample_prop proportion of particles from
  !> \c aero_weight_array_from to \c aero_weight_array_to.
  subroutine aero_weight_array_shift(aero_weight_array_from, &
       aero_weight_array_to, sample_prop, overwrite_to)

    !> Aerosol weight array to shift from.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array_from
    !> Aerosol weight array to shift to.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array_to
    !> Proportion of particles being transfered.
    real(kind=dp), intent(in) :: sample_prop
    !> Whether to overwrite the destination weight (default: no).
    logical, intent(in), optional :: overwrite_to

    call aero_weight_shift(aero_weight_array_from%weight, &
         aero_weight_array_to%weight, sample_prop, overwrite_to)

  end subroutine aero_weight_array_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_array_single_num_conc(aero_weight_array, &
       aero_particle, aero_data)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_weight_array_single_num_conc = aero_weight_num_conc( &
         aero_weight_array%weight(aero_particle%weight_group, &
         aero_particle%weight_class), aero_particle, aero_data)

  end function aero_weight_array_single_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the total number concentration at a given radius (m^3).
  real(kind=dp) function aero_weight_array_num_conc_at_radius( &
       aero_weight_array, i_class, radius)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Weight class number.
    integer, intent(in) :: i_class
    !> Radius to compute number concentration at (m).
    real(kind=dp), intent(in) :: radius

    integer :: i_group
    real(kind=dp) :: num_conc(size(aero_weight_array%weight, 1))

    do i_group = 1,size(aero_weight_array%weight, 1)
       num_conc(i_group) = aero_weight_num_conc_at_radius( &
            aero_weight_array%weight(i_group, i_class), radius)
    end do
    ! harmonic mean (same as summing the computational volumes)
    aero_weight_array_num_conc_at_radius = 1d0 / sum(1d0 / num_conc)

  end function aero_weight_array_num_conc_at_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_array_num_conc(aero_weight_array, &
       aero_particle, aero_data)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_weight_array_num_conc = aero_weight_array_num_conc_at_radius( &
         aero_weight_array, aero_particle%weight_class, &
         aero_particle_radius(aero_particle, aero_data))

  end function aero_weight_array_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether a given aero_weight array is flat in total.
  logical function aero_weight_array_check_flat(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array

    integer :: i_group, i_class

    ! check that all weights have correctly set exponents
    do i_group = 1,size(aero_weight_array%weight, 1)
       do i_class = 1,size(aero_weight_array%weight, 2)
          call aero_weight_check_valid_exponent( &
               aero_weight_array%weight(i_group, i_class))
       end do
    end do

    if (all(abs(sum(aero_weight_array%weight%exponent, 1)) &
         < 1d-20 * sum(abs(aero_weight_array%weight%exponent), 1))) then
       aero_weight_array_check_flat = .true.
    else
       aero_weight_array_check_flat = .false.
    end if

  end function aero_weight_array_check_flat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine whether all weight functions in an array are monotone
  !> increasing, monotone decreasing, or neither.
  subroutine aero_weight_array_check_monotonicity(aero_weight_array, &
       monotone_increasing, monotone_decreasing)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Whether all weights are monotone increasing.
    logical, intent(out) :: monotone_increasing
    !> Whether all weights are monotone decreasing.
    logical, intent(out) :: monotone_decreasing

    integer :: i_group, i_class
    logical :: mono_increasing_array(size(aero_weight_array%weight, 1), &
         size(aero_weight_array%weight, 2))
    logical :: mono_decreasing_array(size(aero_weight_array%weight, 1), &
         size(aero_weight_array%weight, 2))

    do i_group = 1,size(aero_weight_array%weight, 1)
       do i_class = 1,size(aero_weight_array%weight, 2)
          call aero_weight_check_monotonicity( &
               aero_weight_array%weight(i_group, i_class), &
               mono_increasing_array(i_group, i_class), &
               mono_decreasing_array(i_group, i_class))
       end do
    end do

    monotone_increasing = all(mono_increasing_array)
    monotone_decreasing = all(mono_decreasing_array)

  end subroutine aero_weight_array_check_monotonicity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the maximum and minimum number concentrations between the
  !> given radii.
  subroutine aero_weight_array_minmax_num_conc(aero_weight_array, i_class, &
       radius_1, radius_2, num_conc_min, num_conc_max)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Weight class.
    integer, intent(in) :: i_class
    !> First radius.
    real(kind=dp), intent(in) :: radius_1
    !> Second radius.
    real(kind=dp), intent(in) :: radius_2
    !> Minimum number concentration.
    real(kind=dp), intent(out) :: num_conc_min
    !> Maximum number concentration.
    real(kind=dp), intent(out) :: num_conc_max

    real(kind=dp) :: num_conc_1, num_conc_2
    logical :: monotone_increasing, monotone_decreasing

    call aero_weight_array_check_monotonicity(aero_weight_array, &
         monotone_increasing, monotone_decreasing)
    call assert(857727714, monotone_increasing .or. monotone_decreasing)

    num_conc_1 = aero_weight_array_num_conc_at_radius(aero_weight_array, &
         i_class, radius_1)
    num_conc_2 = aero_weight_array_num_conc_at_radius(aero_weight_array, &
         i_class, radius_2)
    num_conc_min = min(num_conc_1, num_conc_2)
    num_conc_max = max(num_conc_1, num_conc_2)

  end subroutine aero_weight_array_minmax_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random group at the given radius, with probability
  !> inversely proportional to group weight at that radius.
  integer function aero_weight_array_rand_group(aero_weight_array, i_class, &
       radius)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Weight class to select from.
    integer, intent(in) :: i_class
    !> Radius to sample group at (m).
    real(kind=dp), intent(in) :: radius

    real(kind=dp) :: comp_vols(size(aero_weight_array%weight, 1))
    integer :: i_group

    do i_group = 1,size(aero_weight_array%weight, 1)
       comp_vols(i_group) = 1d0 / aero_weight_num_conc_at_radius( &
            aero_weight_array%weight(i_group, i_class), radius)
    end do
    aero_weight_array_rand_group = sample_cts_pdf(comp_vols)

  end function aero_weight_array_rand_group

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an aero_weight_array from a spec file.
  subroutine spec_file_read_aero_weight_array(file, aero_weight_array)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol weight array.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array

    integer :: i_group, i_class

    do i_group = 1,size(aero_weight_array%weight, 1)
       do i_class = 1,size(aero_weight_array%weight, 2)
          call spec_file_read_aero_weight(file, &
               aero_weight_array%weight(i_group, i_class))
       end do
    end do

  end subroutine spec_file_read_aero_weight_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_weight_array(val)

    !> Value to pack.
    type(aero_weight_array_t), intent(in) :: val

    integer :: i_group, i_class, total_size
    logical :: is_allocated

    total_size = 0
    is_allocated = allocated(val%weight)
    total_size = total_size + pmc_mpi_pack_size_logical(is_allocated)
    if (is_allocated) then
       total_size = total_size &
            + pmc_mpi_pack_size_integer(size(val%weight, 1)) &
            + pmc_mpi_pack_size_integer(size(val%weight, 2))
       do i_group = 1,size(val%weight, 1)
          do i_class = 1,size(val%weight, 2)
             total_size = total_size &
                  + pmc_mpi_pack_size_aero_weight(val%weight(i_group, i_class))
          end do
       end do
    end if
    pmc_mpi_pack_size_aero_weight_array = total_size

  end function pmc_mpi_pack_size_aero_weight_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_weight_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_weight_array_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i_group, i_class
    logical :: is_allocated

    prev_position = position
    is_allocated = allocated(val%weight)
    call pmc_mpi_pack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       call pmc_mpi_pack_integer(buffer, position, size(val%weight, 1))
       call pmc_mpi_pack_integer(buffer, position, size(val%weight, 2))
       do i_group = 1,size(val%weight, 1)
          do i_class = 1,size(val%weight, 2)
             call pmc_mpi_pack_aero_weight(buffer, position, &
                  val%weight(i_group, i_class))
          end do
       end do
    end if
    call assert(84068036, &
         position - prev_position <= pmc_mpi_pack_size_aero_weight_array(val))
#endif

  end subroutine pmc_mpi_pack_aero_weight_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_weight_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_weight_array_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, n_group, n_class, i_group, i_class
    logical :: is_allocated

    prev_position = position
    call pmc_mpi_unpack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       call pmc_mpi_unpack_integer(buffer, position, n_group)
       call pmc_mpi_unpack_integer(buffer, position, n_class)
       call aero_weight_array_set_sizes(val, n_group, n_class)
       do i_group = 1,size(val%weight, 1)
          do i_class = 1,size(val%weight, 2)
             call pmc_mpi_unpack_aero_weight(buffer, position, &
                  val%weight(i_group, i_class))
          end do
       end do
    else
       if (allocated(val%weight)) then
          deallocate(val%weight)
       end if
    end if
    call assert(321022868, &
         position - prev_position <= pmc_mpi_pack_size_aero_weight_array(val))
#endif

  end subroutine pmc_mpi_unpack_aero_weight_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the \c aero_weight_group dimension to the given NetCDF file
  !> if it is not already present and in any case return the
  !> associated dimid.
  subroutine aero_weight_netcdf_dim_aero_weight_group(aero_weight_array, &
       ncid, dimid_aero_weight_group)

    !> Aero_weight structure array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_aero_weight_group

    integer :: status, i_group, n_group
    integer :: varid_aero_weight_group
    integer :: aero_weight_group_centers(size(aero_weight_array%weight, 1))

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_weight_group", dimid_aero_weight_group)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define it
    call pmc_nc_check(nf90_redef(ncid))

    n_group = size(aero_weight_array%weight, 1)

    call pmc_nc_check(nf90_def_dim(ncid, "aero_weight_group", n_group, &
         dimid_aero_weight_group))
    call pmc_nc_check(nf90_def_var(ncid, "aero_weight_group", NF90_INT, &
         dimid_aero_weight_group, varid_aero_weight_group))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_weight_group, &
         "description", "dummy dimension variable (no useful value)"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_group = 1,n_group
       aero_weight_group_centers(i_group) = i_group
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_weight_group, &
         aero_weight_group_centers))

  end subroutine aero_weight_netcdf_dim_aero_weight_group

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the \c aero_weight_class dimension to the given NetCDF file
  !> if it is not already present and in any case return the
  !> associated dimid.
  subroutine aero_weight_netcdf_dim_aero_weight_class(aero_weight_array, &
       ncid, dimid_aero_weight_class)

    !> Aero_weight structure array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_aero_weight_class

    integer :: status, i_class, n_class
    integer :: varid_aero_weight_class
    integer :: aero_weight_class_centers(size(aero_weight_array%weight, 2))

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_weight_class", dimid_aero_weight_class)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define it
    call pmc_nc_check(nf90_redef(ncid))

    n_class = size(aero_weight_array%weight, 2)

    call pmc_nc_check(nf90_def_dim(ncid, "aero_weight_class", n_class, &
         dimid_aero_weight_class))
    call pmc_nc_check(nf90_def_var(ncid, "aero_weight_class", NF90_INT, &
         dimid_aero_weight_class, varid_aero_weight_class))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_weight_class, &
         "description", "dummy dimension variable (no useful value)"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_class = 1,n_class
       aero_weight_class_centers(i_class) = i_class
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_weight_class, &
         aero_weight_class_centers))

  end subroutine aero_weight_netcdf_dim_aero_weight_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full aero_weight_array.
  subroutine aero_weight_array_output_netcdf(aero_weight_array, ncid)

    !> Aero weight array to write.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_weight_group, dimid_aero_weight_class

    !> \page output_format_aero_weight_array Output File Format: Aerosol Weighting Functions
    !!
    !! The aerosol weighting function NetCDF dimensions are:
    !!   - \b aero_weight_group: number of aerosol weighting groups
    !!   - \b aero_weight_class: number of aerosol weighting classes
    !!
    !! The aerosol weighting function NetCDF variables are:
    !!   - \b weight_type (no unit, dim \c aero_weight): the type of
    !!     each weighting function, with 0 = invalid weight, 1 = no
    !!     weight (\f$w(D) = 1\f$), 2 = power weight (\f$w(D) =
    !!     (D/D_0)^\alpha\f$), 3 = MFA weight (\f$w(D) =
    !!     (D/D_0)^{-3}\f$)
    !!   - \b weight_magnitude (unit m^{-3}, dim \c aero_weight): the
    !!     number concentration magnitude associated with each
    !!     weighting function
    !!   - \b weight_exponent (no unit, dim \c aero_weight): for each
    !!     weighting function, specifies the exponent \f$\alpha\f$ for
    !!     the power \c weight_type, the value -3 for the MFA \c
    !!     weight_type, and zero for any other \c weight_type

    call aero_weight_netcdf_dim_aero_weight_group(aero_weight_array, ncid, &
         dimid_aero_weight_group)
    call aero_weight_netcdf_dim_aero_weight_class(aero_weight_array, ncid, &
         dimid_aero_weight_class)

    call pmc_nc_write_integer_2d(ncid, aero_weight_array%weight%type, &
         "weight_type", &
         (/ dimid_aero_weight_group, dimid_aero_weight_class /), &
         description="type of each aerosol weighting function: 0 = invalid, " &
         // "1 = none (w(D) = 1), 2 = power (w(D) = (D/D_0)^alpha), " &
         // "3 = MFA (mass flow) (w(D) = (D/D_0)^(-3))")
    call pmc_nc_write_real_2d(ncid, aero_weight_array%weight%magnitude, &
         "weight_magnitude", &
         (/ dimid_aero_weight_group, dimid_aero_weight_class /), &
         unit="m^{-3}", &
         description="magnitude for each weighting function")
    call pmc_nc_write_real_2d(ncid, aero_weight_array%weight%exponent, &
         "weight_exponent", &
         (/ dimid_aero_weight_group, dimid_aero_weight_class /), unit="1", &
         description="exponent alpha for the power weight_type, " &
         // "set to -3 for MFA, and zero otherwise")

 end subroutine aero_weight_array_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full aero_weight_array.
  subroutine aero_weight_array_input_netcdf(aero_weight_array, ncid)

    !> Aero weight array to read.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_weight_group, dimid_aero_weight_class, n_group
    integer :: n_class
    character(len=1000) :: name
    integer, allocatable :: type(:, :)
    real(kind=dp), allocatable :: magnitude(:, :), exponent(:, :)

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_weight_group", &
         dimid_aero_weight_group))
    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_weight_class", &
         dimid_aero_weight_class))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, &
         dimid_aero_weight_group, name, n_group))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, &
         dimid_aero_weight_class, name, n_class))
    call assert(719221386, n_group < 1000)
    call assert(520105999, n_class < 1000)

    call pmc_nc_read_integer_2d(ncid, type, "weight_type")
    call pmc_nc_read_real_2d(ncid, magnitude, "weight_magnitude")
    call pmc_nc_read_real_2d(ncid, exponent, "weight_exponent")

    call assert(309191498, size(magnitude) == size(type))
    call assert(588649520, size(magnitude) == size(exponent))

    call aero_weight_array_set_sizes(aero_weight_array, n_group, n_class)

    aero_weight_array%weight%type = type
    aero_weight_array%weight%magnitude = magnitude
    aero_weight_array%weight%exponent = exponent

  end subroutine aero_weight_array_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_weight_array
