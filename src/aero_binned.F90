! Copyright (C) 2005-2015 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for the.

!> \file
!> The pmc_aero_binned module.

!> The aero_binned_t structure and associated subroutines.
module pmc_aero_binned

  use pmc_bin_grid
  use pmc_aero_particle
  use pmc_spec_file
  use pmc_util
  use pmc_bin_grid
  use pmc_aero_dist
  use pmc_mpi
  use pmc_aero_data
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Aerosol number and volume distributions stored per bin.
  !!
  !! These quantities are densities both in volume (per m^3) and in
  !! radius (per log_width). The total concentration per volume is computed as
  !! sum(aero_binned%%num_conc * bin_grid%%log_width).
  !!
  !! An aero_binned_t is similar to an aero_dist_t in that they both
  !! store binned aerosol distributions. The difference is that an
  !! aero_dist_t has the same composition in every bin, whereas an
  !! aero_binned_t can have aerosol composition that varies per bin.
  !!
  !! By convention, if aero_binned_is_allocated() return \c .false.,
  !! then the aero_binned is treated as zero for all operations on
  !! it. This will be the case for new \c aero_binned_t structures.
  type aero_binned_t
     !> Number concentration per bin (#/m^3/log_width).
     !! Array length is \c bin_grid_size(bin_grid).
     real(kind=dp), allocatable :: num_conc(:)
     !> Volume concentration per bin and per species (m^3/m^3/log_width).
     !! Array size is <tt>bin_grid_size(bin_grid) x
     !! aero_data_n_spec(aero_data)</tt>.
     real(kind=dp), allocatable :: vol_conc(:,:)
  end type aero_binned_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine whether the \c aero_binned is correctly allocated.
  logical function aero_binned_is_allocated(aero_binned)

    !> Aerosol binned to check.
    type(aero_binned_t), intent(in) :: aero_binned

    logical :: valid

    valid = .true.
    valid = valid .and. allocated(aero_binned%num_conc)
    valid = valid .and. allocated(aero_binned%vol_conc)
    valid = valid &
         .and. (size(aero_binned%num_conc) == size(aero_binned%num_conc, 1))
    aero_binned_is_allocated = valid

  end function aero_binned_is_allocated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the number of bins and species in an aero_binned_t.
  subroutine aero_binned_set_sizes(aero_binned, n_bin, n_spec)

    !> Structure to be allocated.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Number of aerosol bins to allocate (typically \c bin_grid%%n_bin).
    integer, intent(in) :: n_bin
    !> Number of aerosol species to allocate (typically
    !> \c aero_data%%n_spec).
    integer, intent(in) :: n_spec

    if (allocated(aero_binned%num_conc)) deallocate(aero_binned%num_conc)
    if (allocated(aero_binned%vol_conc)) deallocate(aero_binned%vol_conc)
    allocate(aero_binned%num_conc(n_bin))
    allocate(aero_binned%vol_conc(n_bin, n_spec))
    call aero_binned_zero(aero_binned)

  end subroutine aero_binned_set_sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set all internal data in an aero_binned_t structure to zero.
  subroutine aero_binned_zero(aero_binned)

    !> Structure to zero.
    type(aero_binned_t), intent(inout) :: aero_binned

    if (aero_binned_is_allocated(aero_binned)) then
       aero_binned%num_conc = 0d0
       aero_binned%vol_conc = 0d0
    end if

  end subroutine aero_binned_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add two aero_binned_t structures together.
  !!
  !! Symbolically does aero_binned = aero_binned + aero_binned_delta.
  subroutine aero_binned_add(aero_binned, aero_binned_delta)

    !> Base aero_binned_t structure that will be added to.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to add to aero_binned.
    type(aero_binned_t), intent(in) :: aero_binned_delta

    if (aero_binned_is_allocated(aero_binned_delta)) then
       if (aero_binned_is_allocated(aero_binned)) then
          aero_binned%num_conc = aero_binned%num_conc &
               + aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned%vol_conc &
               + aero_binned_delta%vol_conc
       else
          aero_binned%num_conc = aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned_delta%vol_conc
       end if
    end if

  end subroutine aero_binned_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a scaled \c aero_binned_t structure to an existing one.
  !!
  !! Symbolically does aero_binned = aero_binned + alpha * aero_binned_delta.
  subroutine aero_binned_add_scaled(aero_binned, aero_binned_delta, alpha)

    !> Base aero_binned_t structure that will be added to.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to add to aero_binned.
    type(aero_binned_t), intent(in) :: aero_binned_delta
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    if (aero_binned_is_allocated(aero_binned_delta)) then
       if (aero_binned_is_allocated(aero_binned)) then
          aero_binned%num_conc = aero_binned%num_conc &
               + alpha * aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned%vol_conc &
               + alpha * aero_binned_delta%vol_conc
       else
          aero_binned%num_conc = aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned_delta%vol_conc
       end if
    end if

  end subroutine aero_binned_add_scaled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Subtract one aero_binned_t structure from another.
  !!
  !! Symbolically does aero_binned = aero_binned - aero_binned_delta.
  subroutine aero_binned_sub(aero_binned, aero_binned_delta)

    !> Base aero_binned_t structure that will be subtracted from.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to subtract from aero_binned.
    type(aero_binned_t), intent(in) :: aero_binned_delta

    if (aero_binned_is_allocated(aero_binned_delta)) then
       if (aero_binned_is_allocated(aero_binned)) then
          aero_binned%num_conc = aero_binned%num_conc &
               - aero_binned_delta%num_conc
          aero_binned%vol_conc = aero_binned%vol_conc &
               - aero_binned_delta%vol_conc
       else
          aero_binned%num_conc = - aero_binned_delta%num_conc
          aero_binned%vol_conc = - aero_binned_delta%vol_conc
       end if
    end if

  end subroutine aero_binned_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scale an aero_binned_t by a real number.
  !!
  !! Symbolically does aero_binned = aero_binned * alpha.
  subroutine aero_binned_scale(aero_binned, alpha)

    !> Base aero_binned to scale.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Scale factor.
    real(kind=dp), intent(in) :: alpha

    if (aero_binned_is_allocated(aero_binned)) then
       aero_binned%num_conc = aero_binned%num_conc * alpha
       aero_binned%vol_conc = aero_binned%vol_conc * alpha
    end if

  end subroutine aero_binned_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scales an aero_binned_t element-wise by an array of reals.
  subroutine aero_binned_scale_by_array(aero_binned, alpha_array)

    !> Base aero_binned_t structure that will be scaled.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Structure to scale aero_binned.
    real(kind=dp), allocatable, intent(in) :: alpha_array(:)

    integer :: i

    do i = 1, size(aero_binned%num_conc)
       aero_binned%num_conc(i) = alpha_array(i)*aero_binned%num_conc(i)
       aero_binned%vol_conc(i,:) = alpha_array(i)*aero_binned%vol_conc(i,:)
    end do

  end subroutine aero_binned_scale_by_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add an aero_dist_t to an aero_binned_t.
  !!
  !! Symbolically does aero_binned = aero_binned + aero_dist.
  subroutine aero_binned_add_aero_dist(aero_binned, bin_grid, &
       aero_data, aero_dist)

    !> Base aero_binned_t structure to add to.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> The aero_dist_t structure to add.
    type(aero_dist_t), intent(in) :: aero_dist

    real(kind=dp) :: dist_num_conc(bin_grid_size(bin_grid))
    real(kind=dp) :: dist_vol_conc(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    call aero_dist_num_conc(aero_dist, bin_grid, aero_data, &
         dist_num_conc)
    call aero_dist_vol_conc(aero_dist, bin_grid, aero_data, &
         dist_vol_conc)
    if (aero_binned_is_allocated(aero_binned)) then
       aero_binned%num_conc = aero_binned%num_conc + dist_num_conc
       aero_binned%vol_conc = aero_binned%vol_conc + dist_vol_conc
    else
       aero_binned%num_conc = dist_num_conc
       aero_binned%vol_conc = dist_vol_conc
    end if

  end subroutine aero_binned_add_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the number of bytes required to pack the structure.
  !!
  !! See pmc_mpi for usage details.
  integer function pmc_mpi_pack_size_aero_binned(val)

    !> Structure to pack.
    type(aero_binned_t), intent(in) :: val

    pmc_mpi_pack_size_aero_binned = &
         pmc_mpi_pack_size_real_array(val%num_conc) &
         + pmc_mpi_pack_size_real_array_2d(val%vol_conc)

  end function pmc_mpi_pack_size_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the structure into the buffer and advance position.
  !!
  !! See pmc_mpi for usage details.
  subroutine pmc_mpi_pack_aero_binned(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Structure to pack.
    type(aero_binned_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%num_conc)
    call pmc_mpi_pack_real_array_2d(buffer, position, val%vol_conc)
    call assert(348207873, &
         position - prev_position <= pmc_mpi_pack_size_aero_binned(val))
#endif

  end subroutine pmc_mpi_pack_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the structure from the buffer and advance position.
  !!
  !! See pmc_mpi for usage details.
  subroutine pmc_mpi_unpack_aero_binned(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Structure to unpack into (must not be allocated).
    type(aero_binned_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%num_conc)
    call pmc_mpi_unpack_real_array_2d(buffer, position, val%vol_conc)
    call assert(878267066, &
         position - prev_position <= pmc_mpi_pack_size_aero_binned(val))
#endif

  end subroutine pmc_mpi_unpack_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of the structure across all processes,
  !> storing the result on the root process.
  subroutine pmc_mpi_reduce_avg_aero_binned(val, val_avg)

    !> Per-process value to average.
    type(aero_binned_t), intent(in) :: val
    !> Averaged result (only valid on root process).
    type(aero_binned_t), intent(inout) :: val_avg

    call pmc_mpi_reduce_avg_real_array(val%num_conc, val_avg%num_conc)
    call pmc_mpi_reduce_avg_real_array_2d(val%vol_conc, val_avg%vol_conc)

  end subroutine pmc_mpi_reduce_avg_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_binned_output_netcdf(aero_binned, ncid, bin_grid, &
       aero_data)

    !> Aero_binned to write.
    type(aero_binned_t), intent(in) :: aero_binned
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data

    integer :: dimid_aero_diam, dimid_aero_species
    real(kind=dp) :: mass_conc(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))
    integer :: i_bin

    !> \page output_format_aero_binned Output File Format: Aerosol Binned Sectional State
    !!
    !! The aerosol size distributions (number and mass) are stored on
    !! a logarmithmic grid (see the \ref output_format_diam_bin_grid
    !! section). To compute the total number or mass concentration,
    !! compute the sum over \c i of <tt>aero_number_concentration(i) *
    !! aero_diam_widths(i)</tt>, for example.
    !!
    !! The aerosol binned sectional state uses the \c aero_species
    !! NetCDF dimension as specified in the \ref
    !! output_format_aero_data section, as well as the \c aero_diam
    !! NetCDF dimension specified in the \ref
    !! output_format_diam_bin_grid section.
    !!
    !! The aerosol binned sectional state NetCDF variables are:
    !!   - \b aero_number_concentration (unit 1/m^3, dim \c aero_diam): the
    !!     number size distribution for the aerosol population,
    !!     \f$ dN(r)/d\ln r \f$, per bin
    !!   - \b aero_mass_concentration (unit kg/m^3, dim
    !!     <tt>dimid_aero_diam x dimid_aero_species</tt>): the mass size
    !!     distribution for the aerosol population,
    !!     \f$ dM(r,s)/d\ln r \f$, per bin and per species

    ! output_format_diam_bin_grid is here, as this is the only place it's used

    !> \page output_format_diam_bin_grid Output File Format: Diameter Bin Grid Data
    !!
    !! The aerosol diameter bin grid data NetCDF dimensions are:
    !!   - \b aero_diam: number of bins (grid cells) on the diameter axis
    !!   - \b aero_diam_edges: number of bin edges (grid cell edges) on
    !!     the diameter axis --- always equal to <tt>aero_diam + 1</tt>
    !!
    !! The aerosol diameter bin grid data NetCDF variables are:
    !!   - \b aero_diam (unit m, dim \c aero_diam): aerosol diameter axis
    !!     bin centers --- centered on a logarithmic scale from the edges, so
    !!     that <tt>aero_diam(i) / aero_diam_edges(i) =
    !!     sqrt(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>
    !!   - \b aero_diam_edges (unit m, dim \c aero_diam_edges): aersol
    !!     diameter axis bin edges (there is one more edge than center)
    !!   - \b aero_diam_widths (dimensionless, dim \c aero_diam):
    !!     the base-e logarithmic bin widths --- <tt>aero_diam_widths(i)
    !!     = ln(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>, so
    !!     all bins have the same width
    !!
    !! See also:
    !!   - \ref input_format_diam_bin_grid --- the corresponding input format

    do i_bin = 1,bin_grid_size(bin_grid)
       mass_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) * aero_data%density
    end do

    call bin_grid_netcdf_dim(bin_grid, ncid, "aero_diam", "m", &
         dimid_aero_diam, "aerosol diameter", scale=2d0)
    call aero_data_netcdf_dim_aero_species(aero_data, ncid, dimid_aero_species)

    call pmc_nc_write_real_1d(ncid, aero_binned%num_conc, &
         "aero_number_concentration", (/ dimid_aero_diam /), &
         unit="1/m^3", &
         long_name="aerosol number size concentration distribution", &
         description="logarithmic number size concentration, " &
         // "d N(r)/d ln D --- multiply by aero_diam_widths(i) " &
         // "and sum over i to compute the total number concentration")
    call pmc_nc_write_real_2d(ncid, mass_conc, &
         "aero_mass_concentration", &
         (/ dimid_aero_diam, dimid_aero_species /), unit="kg/m^3", &
         long_name="aerosol mass size concentration distribution", &
         description="logarithmic mass size concentration per species, " &
         // "d M(r,s)/d ln D --- multiply by aero_diam_widths(i) " &
         // "and sum over i to compute the total mass concentration of " &
         // "species s")

  end subroutine aero_binned_output_netcdf

  ! output_format_diam_bin_grid is here, as this is the only place it's used

  ! this belongs in the subroutine above, but is outside because
  ! Doxygen 1.8.7 doesn't resolve references when multiple \page
  ! blocks are in one subroutine

  !> \page output_format_diam_bin_grid Output File Format: Diameter Bin Grid Data
  !!
  !! The aerosol diameter bin grid data NetCDF dimensions are:
  !!   - \b aero_diam: number of bins (grid cells) on the diameter axis
  !!   - \b aero_diam_edges: number of bin edges (grid cell edges) on
  !!     the diameter axis --- always equal to <tt>aero_diam + 1</tt>
  !!
  !! The aerosol diameter bin grid data NetCDF variables are:
  !!   - \b aero_diam (unit m, dim \c aero_diam): aerosol diameter axis
  !!     bin centers --- centered on a logarithmic scale from the edges, so
  !!     that <tt>aero_diam(i) / aero_diam_edges(i) =
  !!     sqrt(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>
  !!   - \b aero_diam_edges (unit m, dim \c aero_diam_edges): aersol
  !!     diameter axis bin edges (there is one more edge than center)
  !!   - \b aero_diam_widths (dimensionless, dim \c aero_diam):
  !!     the base-e logarithmic bin widths --- <tt>aero_diam_widths(i)
  !!     = ln(aero_diam_edges(i+1) / aero_diam_edges(i))</tt>, so
  !!     all bins have the same width
  !!
  !! See also:
  !!   - \ref input_format_diam_bin_grid --- the corresponding input format

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_binned_input_netcdf(aero_binned, ncid, bin_grid, &
       aero_data)

    !> Aero_binned to write.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> aero_data structure.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_bin

    call pmc_nc_read_real_1d(ncid, aero_binned%num_conc, &
         "aero_number_concentration")
    call pmc_nc_read_real_2d(ncid, aero_binned%vol_conc, &
         "aero_mass_concentration")
    ! convert mass concentation to volume concentration
    do i_bin = 1,bin_grid_size(bin_grid)
       aero_binned%vol_conc(i_bin,:) = aero_binned%vol_conc(i_bin,:) &
            / aero_data%density
    end do

  end subroutine aero_binned_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_binned
