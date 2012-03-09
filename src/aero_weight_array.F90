! Copyright (C) 2012 Jeff Curtis
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
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> An array of aerosol size distribution weighting functions.
  type aero_weight_array_t
     !> Aero weight array.
     type(aero_weight_t), pointer :: weight(:)
  end type aero_weight_array_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_weight_array.
  subroutine aero_weight_array_allocate(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(out) :: aero_weight_array

    allocate(aero_weight_array%weight(0))

  end subroutine aero_weight_array_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_weight_array to the given size.
  subroutine aero_weight_array_allocate_size(aero_weight_array, n_weight)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(out) :: aero_weight_array
    !> Number of weights.
    integer,intent(in) :: n_weight

    allocate(aero_weight_array%weight(n_weight))

    call aero_weight_allocate(aero_weight_array%weight)

  end subroutine aero_weight_array_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aero_weight_array_deallocate(aero_weight_array)

    !> Structure to deallocate.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array

    call aero_weight_deallocate(aero_weight_array%weight)

    deallocate(aero_weight_array%weight)

  end subroutine aero_weight_array_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the contents of the \c aero_weight_array.
  subroutine aero_weight_array_zero(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array

    call aero_weight_zero(aero_weight_array%weight)

  end subroutine aero_weight_array_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the \c aero_weight_array computational volume.
  elemental subroutine aero_weight_array_zero_comp_vol(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array

    call aero_weight_zero_comp_vol(aero_weight_array%weight)

  end subroutine aero_weight_array_zero_comp_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an \c aero_weight_array.
  subroutine aero_weight_array_copy(aero_weight_array_from, &
       aero_weight_array_to)

    !> Origin aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array_from
    !> Destination aerosol weight array.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array_to

    if (size(aero_weight_array_to%weight)  &
         /= size(aero_weight_array_from%weight)) then
       deallocate(aero_weight_array_to%weight)
       call aero_weight_array_allocate_size(aero_weight_array_to, &
            size(aero_weight_array_from%weight))
    end if
    call aero_weight_copy(aero_weight_array_from%weight, &
         aero_weight_array_to%weight)

  end subroutine aero_weight_array_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale the computational volume by the given fraction, so
  !> <tt>new_comp_vol = old_comp_vol * fraction</tt>.
  subroutine aero_weight_array_scale_comp_vol(aero_weight_array, fraction)

    !> Aerosol weight array to halve.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array
    !> Fraction to scale computational volume by.
    real(kind=dp), intent(in) :: fraction

    integer :: i

    call aero_weight_scale_comp_vol(aero_weight_array%weight, fraction)

  end subroutine aero_weight_array_scale_comp_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add the computational volume of \c aero_weight_array_delta to
  !> \c aero_weight_array.
  subroutine aero_weight_array_add_comp_vol(aero_weight_array, &
      aero_weight_array_delta)

    !> Aerosol weight array to add volume to.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array
    !> Aerosol weight array to add volume from.
    type(aero_weight_array_t), intent(in) :: aero_weight_array_delta

    call aero_weight_add_comp_vol(aero_weight_array%weight, &
            aero_weight_array_delta%weight)

  end subroutine aero_weight_array_add_comp_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Transfer the computational volume from \c aero_weight_from to \c
  !> aero_weight_to, weighted by \c sample_prop.
  subroutine aero_weight_array_transfer_comp_vol( &
       aero_weight_array_from,aero_weight_array_to, sample_prop)

    !> Aerosol weight array to take volume from.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array_from
    !> Aerosol weight array to add volume to.
    type(aero_weight_array_t), intent(inout) :: aero_weight_array_to
    !> Proportion of from volume to transfer.
    real(kind=dp), intent(in) :: sample_prop

    real(kind=dp) :: transfer_comp_vol

    call aero_weight_transfer_comp_vol(aero_weight_array_from%weight, &
         aero_weight_array_to%weight,sample_prop)

  end subroutine aero_weight_array_transfer_comp_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_array_single_num_conc(aero_weight_array, &
       aero_particle)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_weight_array_single_num_conc = aero_weight_num_conc( &
         aero_weight_array%weight(aero_particle%weight_group), aero_particle)

  end function aero_weight_array_single_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the total number concentration at a given radius (m^3).
  real(kind=dp) function aero_weight_array_num_conc_at_radius( &
       aero_weight_array, radius)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Radius to compute number concentration at (m).
    real(kind=dp), intent(in) :: radius

    integer :: i_group
    real(kind=dp) :: num_conc(size(aero_weight_array%weight))

    do i_group = 1,size(aero_weight_array%weight)
       num_conc(i_group) = aero_weight_num_conc_at_radius( &
            aero_weight_array%weight(i_group), radius)
    end do
    ! harmonic mean (same as summing the computational volumes)
    aero_weight_array_num_conc_at_radius = 1d0 / sum(1d0 / num_conc)

  end function aero_weight_array_num_conc_at_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the number concentration for a particle (m^{-3}).
  real(kind=dp) function aero_weight_array_num_conc(aero_weight_array, &
       aero_particle)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Aerosol particle to compute number concentration for.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_weight_array_num_conc = aero_weight_array_num_conc_at_radius( &
         aero_weight_array, aero_particle_radius(aero_particle))

  end function aero_weight_array_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether a given aero_weight array is flat in total.
  logical function aero_weight_array_check_flat(aero_weight_array)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array

    integer :: i_group

    ! check we know about all the weight types
    do i_group = 1,size(aero_weight_array%weight)
       call aero_weight_check_flat(aero_weight_array%weight(i_group))
    end do

    if (abs(sum(aero_weight_array%weight%exponent)) &
         < 1d-20 * sum(abs(aero_weight_array%weight%exponent))) then
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

    integer :: i_group

    monotone_increasing = .true.
    monotone_decreasing = .true.
    do i_group = 1,size(aero_weight_array%weight)
       call aero_weight_check_monotonicity(aero_weight_array%weight(i_group), &
            monotone_increasing, monotone_decreasing)
    end do

  end subroutine aero_weight_array_check_monotonicity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the maximum and minimum number concentrations between the
  !> given radii.
  subroutine aero_weight_array_minmax_num_conc(aero_weight_array, radius_1, &
       radius_2, num_conc_min, num_conc_max)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
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
         radius_1)
    num_conc_2 = aero_weight_array_num_conc_at_radius(aero_weight_array, &
         radius_2)
    num_conc_min = min(num_conc_1, num_conc_2)
    num_conc_max = max(num_conc_1, num_conc_2)

  end subroutine aero_weight_array_minmax_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Choose a random group at the given radius, with probability
  !> proportional to group volume at that radius.
  integer function aero_weight_array_rand_group(aero_weight_array, radius)

    !> Aerosol weight array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> Radius to sample group at (m).
    real(kind=dp), intent(in) :: radius

    real(kind=dp) :: comp_vols(size(aero_weight_array%weight))
    integer :: i

    do i = 1,size(aero_weight_array%weight)
       comp_vols(i) = 1d0 / aero_weight_num_conc_at_radius( &
            aero_weight_array%weight(i), radius)
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

    integer :: i

    do i = 1,size(aero_weight_array%weight)
       call spec_file_read_aero_weight(file,aero_weight_array%weight(i))
    end do

  end subroutine spec_file_read_aero_weight_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_weight_array(val)

    !> Value to pack.
    type(aero_weight_array_t), intent(in) :: val

    integer :: i

    pmc_mpi_pack_size_aero_weight_array = &
         pmc_mpi_pack_size_integer(size(val%weight))

    do i = 1, size(val%weight)
       pmc_mpi_pack_size_aero_weight_array = &
            pmc_mpi_pack_size_aero_weight_array  &
            + pmc_mpi_pack_size_aero_weight(val%weight(i))
    end do

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
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, size(val%weight))
    do i = 1, size(val%weight)
       call pmc_mpi_pack_aero_weight(buffer,position,val%weight(i))
    end do
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
    integer :: prev_position, array_size

    call aero_weight_array_deallocate(val)
    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, array_size)
    call aero_weight_array_allocate_size(val, array_size)
    do i = 1, size(val%weight)
       call pmc_mpi_unpack_aero_weight(buffer, position, val%weight(i))
    end do
    call assert(321022868, &
         position - prev_position <= pmc_mpi_pack_size_aero_weight_array(val))
#endif

  end subroutine pmc_mpi_unpack_aero_weight_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the \c aero_weight dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine aero_weight_netcdf_dim_aero_weight(aero_weight_array, ncid, &
       dimid_aero_weight)

    !> Aero_weight structure array.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_aero_weight

    integer :: status, i_weight, n_weight
    integer :: varid_aero_weight
    integer :: aero_weight_centers(size(aero_weight_array%weight))

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_weight", dimid_aero_weight)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define it
    call pmc_nc_check(nf90_redef(ncid))

    n_weight = size(aero_weight_array%weight)

    call pmc_nc_check(nf90_def_dim(ncid, "aero_weight", n_weight, &
         dimid_aero_weight))
    call pmc_nc_check(nf90_def_var(ncid, "aero_weight", NF90_INT, &
         dimid_aero_weight, varid_aero_weight))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_weight, "description", &
         "dummy dimension variable (no useful value)"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_weight = 1,n_weight
       aero_weight_centers(i_weight) = i_weight
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_weight, &
         aero_weight_centers))

  end subroutine aero_weight_netcdf_dim_aero_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full aero_weight_array.
  subroutine aero_weight_array_output_netcdf(aero_weight_array, ncid)

    !> Aero weight array to write.
    type(aero_weight_array_t), intent(in) :: aero_weight_array
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_weight

    !> \page output_format_aero_weight_array Output File Format: Aerosol Weighting Functions
    !!
    !! The aerosol weighting function NetCDF dimensions are:
    !!   - \b aero_weight: number of aerosol weighting functions
    !!
    !! The aerosol weighting function NetCDF variables are:
    !!   - \b weight_comp_vol (unit m^3, dim \c aero_weight): the
    !!     computational volume associated with each weighting
    !!     function
    !!   - \b weight_type (no unit, dim \c aero_weight): the type of
    !!     each weighting function, with 0 = invalid weight, 1 = no
    !!     weight (\f$w(D) = 1\f$), 2 = power weight (\f$w(D) =
    !!     (D/D_0)^\alpha\f$), 3 = MFA weight (\f$w(D) =
    !!     (D/D_0)^{-3}\f$)
    !!   - \b weight_exponent (no unit, dim \c aero_weight): for each
    !!     weighting function, specifies the exponent \f$\alpha\f$ for
    !!     the power \c weight_type, the value -3 for the MFA \c
    !!     weight_type, and zero for any other \c weight_type

    call aero_weight_netcdf_dim_aero_weight(aero_weight_array, ncid, &
         dimid_aero_weight)

    call pmc_nc_write_real_1d(ncid, aero_weight_array%weight%comp_vol, &
         "weight_comp_vol", (/ dimid_aero_weight /), unit="m^3", &
         description="computational volume for each weighting function")
    call pmc_nc_write_integer_1d(ncid, aero_weight_array%weight%type, &
         "weight_type", (/ dimid_aero_weight /), &
         description="type of each aerosol weighting function: 0 = invalid, " &
         // "1 = none (w(D) = 1), 2 = power (w(D) = (D/D_0)^alpha), " &
         // "3 = MFA (mass flow) (w(D) = (D/D_0)^(-3))")
    call pmc_nc_write_real_1d(ncid, aero_weight_array%weight%exponent, &
         "weight_exponent", (/ dimid_aero_weight /), unit="1", &
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

    integer :: dimid_aero_weight, n_weight
    character(len=1000) :: name
    real(kind=dp), allocatable :: comp_vol(:), exponent(:)
    integer, allocatable :: type(:)

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_weight", dimid_aero_weight))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, &
         dimid_aero_weight, name, n_weight))
    call assert(719221386, n_weight < 1000)

    allocate(comp_vol(n_weight))
    allocate(type(n_weight))
    allocate(exponent(n_weight))

    call pmc_nc_read_real_1d(ncid, comp_vol, "weight_comp_vol")
    call pmc_nc_read_integer_1d(ncid, type, "weight_type")
    call pmc_nc_read_real_1d(ncid, exponent, "weight_exponent")

    call assert(309191498, size(comp_vol) == size(type))
    call assert(588649520, size(comp_vol) == size(exponent))

    call aero_weight_array_deallocate(aero_weight_array)
    call aero_weight_array_allocate_size(aero_weight_array, size(comp_vol))

    do i = 1,size(aero_weight_array%weight)
       aero_weight_array%weight(i)%comp_vol = comp_vol(i)
       aero_weight_array%weight(i)%type = type(i)
       aero_weight_array%weight(i)%ref_radius = 1d0
       aero_weight_array%weight(i)%exponent = exponent(i)
    end do

    deallocate(comp_vol)
    deallocate(type)
    deallocate(exponent)

  end subroutine aero_weight_array_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_weight_array
