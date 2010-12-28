! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_bin_grid module.

!> The bin_grid_t structure and associated subroutines.
module pmc_bin_grid

  use pmc_constants
  use pmc_util
  use pmc_spec_file
  use pmc_mpi
  use pmc_netcdf
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> 1D grid of size bins.
  !!
  !! The grid of bins is logarithmically spaced in volume, an
  !! assumption that is quite heavily incorporated into the code. At
  !! some point in the future it would be nice to relax this
  !! assumption.
  type bin_grid_t
     !> Number of bins.
     integer :: n_bin
     !> Len n_bin, bin center radii (m^3).
     real(kind=dp), pointer :: center_radius(:)
     !> Len (n_bin + 1), bin edge radii (m^3).
     real(kind=dp), pointer :: edge_radius(:)
     !> Bin logarithmic width, equal to <tt>log(edge_radius(i+1)) -
     !> log(edge_radius(i))</tt> for any \c i (dimensionless).
     real(kind=dp) :: log_width
  end type bin_grid_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a bin_grid.
  subroutine bin_grid_allocate(bin_grid)

    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid

    bin_grid%n_bin = 0
    allocate(bin_grid%center_radius(0))
    allocate(bin_grid%edge_radius(0))

  end subroutine bin_grid_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a bin_grid of the given size.
  subroutine bin_grid_allocate_size(bin_grid, n_bin)

    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Number of bins.
    integer, intent(in) :: n_bin

    bin_grid%n_bin = n_bin
    allocate(bin_grid%center_radius(n_bin))
    allocate(bin_grid%edge_radius(n_bin + 1))

  end subroutine bin_grid_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all memory.
  subroutine bin_grid_deallocate(bin_grid)

    !> Bin_grid to free.
    type(bin_grid_t), intent(inout) :: bin_grid

    deallocate(bin_grid%center_radius)
    deallocate(bin_grid%edge_radius)

  end subroutine bin_grid_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a concentration f(vol)d(vol) to f(ln(r))d(ln(r))
  !> where vol = 4/3 pi r^3.
  subroutine vol_to_lnr(r, f_vol, f_lnr)
    
    !> Radius (m).
    real(kind=dp), intent(in) :: r
    !> Concentration as a function of volume.
    real(kind=dp), intent(in) :: f_vol
    !> Concentration as a function of ln(r).
    real(kind=dp), intent(out) :: f_lnr
    
    f_lnr = f_vol * 4d0 * const%pi * r**3
    
  end subroutine vol_to_lnr
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates the bin grid given the range and number of bins.
  subroutine bin_grid_make(bin_grid, n_bin, r_min, r_max)
    
    !> New bin grid.
    type(bin_grid_t), intent(inout) :: bin_grid
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Minimum radius (m^3).
    real(kind=dp), intent(in) :: r_min
    !> Minimum radius (m^3).
    real(kind=dp), intent(in) :: r_max

    integer :: i_bin

    call assert_msg(538534122, n_bin > 0, &
         "bin_grid requires a positive n_bin, not: " &
         // trim(integer_to_string(n_bin)))
    call assert_msg(966541762, r_min > 0d0, &
         "bin_grid requires a positive r_min, not: " &
         // trim(real_to_string(r_min)))
    call assert_msg(966541762, r_min < r_max, &
         "bin_grid requires r_min < r_max, not: " &
         // trim(real_to_string(r_min)) // " and " &
         // trim(real_to_string(r_max)))
    call bin_grid_deallocate(bin_grid)
    call bin_grid_allocate_size(bin_grid, n_bin)
    call logspace(r_min, r_max, bin_grid%edge_radius)
    do i_bin = 1,n_bin
       bin_grid%center_radius(i_bin) &
            = exp(0.5d0 * log(bin_grid%edge_radius(i_bin)) &
            + 0.5d0 * log(bin_grid%edge_radius(i_bin + 1)))
    end do
    bin_grid%log_width = (log(r_max) - log(r_min)) / real(n_bin, kind=dp)

  end subroutine bin_grid_make

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the bin number that contains a given particle.
  !!
  !! This assumes logarithmically spaced bins. If a particle is below
  !! the smallest bin or above the largest bin, then it is returned as
  !! being in the smallest or largest bin, respectively.
  integer function bin_grid_particle_in_bin(bin_grid, radius)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Radius of particle.
    real(kind=dp), intent(in) :: radius

    call assert(448215689, bin_grid%n_bin >= 1)
    bin_grid_particle_in_bin = logspace_find(bin_grid%edge_radius(1), &
         bin_grid%edge_radius(bin_grid%n_bin + 1), bin_grid%n_bin + 1, &
         radius)

  end function bin_grid_particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a bin_grid from a spec file and
  !> generate it.
  subroutine spec_file_read_bin_grid(file, bin_grid)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Bin grid.
    type(bin_grid_t), intent(inout) :: bin_grid

    integer :: n_bin
    real(kind=dp) :: d_min, d_max

    !> \page input_format_bin_grid Input File Format: Diameter Axis Bin Grid
    !!
    !! The bin grid is logarithmic in diameter, consisting of
    !! \f$n_{\rm bin}\f$ bins with centers \f$c_i\f$ (\f$i =
    !! 1,\ldots,n_{\rm bin}\f$) and edges \f$e_i\f$ (\f$i =
    !! 1,\ldots,(n_{\rm bin} + 1)\f$) such that \f$e_{i+1}/e_i\f$ is a
    !! constant and \f$c_i/e_i = \sqrt{e_{i+1}/e_i}\f$. That is,
    !! \f$\ln(e_i)\f$ are uniformly spaced and \f$\ln(c_i)\f$ are the
    !! arithmetic centers.
    !!
    !! The diameter axis bin grid is specified by the parameters:
    !!   - \b n_bin (integer): The number of bins \f$n_{\rm bin}\f$.
    !!   - \b d_min (real, unit m): The left edge of the left-most bin,
    !!     \f$e_1\f$.
    !!   - \b d_max (real, unit m): The right edge of the right-most bin,
    !!     \f$e_{n_{\rm bin} + 1}\f$.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_bin_grid --- the corresponding output format

    call spec_file_read_integer(file, 'n_bin', n_bin)
    call spec_file_read_real(file, 'd_min', d_min)
    call spec_file_read_real(file, 'd_max', d_max)
    call bin_grid_make(bin_grid, n_bin, diam2rad(d_min), diam2rad(d_max))

  end subroutine spec_file_read_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_bin_grid(val)

    !> Value to pack.
    type(bin_grid_t), intent(in) :: val

    pmc_mpi_pack_size_bin_grid = &
         pmc_mpi_pack_size_integer(val%n_bin) &
         + pmc_mpi_pack_size_real_array(val%center_radius) &
         + pmc_mpi_pack_size_real_array(val%edge_radius) &
         + pmc_mpi_pack_size_real(val%log_width)

  end function pmc_mpi_pack_size_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_bin_grid(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(bin_grid_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_bin)
    call pmc_mpi_pack_real_array(buffer, position, val%center_radius)
    call pmc_mpi_pack_real_array(buffer, position, val%edge_radius)
    call pmc_mpi_pack_real(buffer, position, val%log_width)
    call assert(385455586, &
         position - prev_position == pmc_mpi_pack_size_bin_grid(val))
#endif

  end subroutine pmc_mpi_pack_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_bin_grid(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(bin_grid_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_bin)
    call pmc_mpi_unpack_real_array(buffer, position, val%center_radius)
    call pmc_mpi_unpack_real_array(buffer, position, val%edge_radius)
    call pmc_mpi_unpack_real(buffer, position, val%log_width)
    call assert(741838730, &
         position - prev_position == pmc_mpi_pack_size_bin_grid(val))
#endif

  end subroutine pmc_mpi_unpack_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the aero_diam dimension to the given NetCDF file if it is
  !> not already present and in any case return the associated dimid.
  subroutine bin_grid_netcdf_dim_aero_diam(bin_grid, ncid, &
       dimid_aero_diam)

    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the aero_diam dimension.
    integer, intent(out) :: dimid_aero_diam

    integer :: status, varid_aero_diam
    integer :: dimid_aero_diam_edges, varid_aero_diam_edges, &
         varid_aero_diam_widths
    real(kind=dp) :: aero_diam_centers(bin_grid%n_bin), &
         aero_diam_edges(bin_grid%n_bin + 1)
    real(kind=dp) :: aero_diam_widths(bin_grid%n_bin)

    status = nf90_inq_dimid(ncid, "aero_diam", dimid_aero_diam)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_diam", &
         bin_grid%n_bin, dimid_aero_diam))
    call pmc_nc_check(nf90_def_dim(ncid, "aero_diam_edges", &
         bin_grid%n_bin + 1, dimid_aero_diam_edges))

    call pmc_nc_check(nf90_enddef(ncid))

    aero_diam_centers = rad2diam(bin_grid%center_radius)
    aero_diam_widths = bin_grid%log_width
    aero_diam_edges = rad2diam(bin_grid%edge_radius)

    call pmc_nc_write_real_1d(ncid, aero_diam_centers, &
         "aero_diam", (/ dimid_aero_diam /), unit="m", &
         long_name="aerosol radius axis bin centers", &
         description="logarithmically spaced centers of radius axis grid, " &
         // "so that aero_diam(i) / aero_diam_edges(i) = " &
         // "0.5 * aero_diam_edges(i+1) / aero_diam_edges(i)")
    call pmc_nc_write_real_1d(ncid, aero_diam_edges, &
         "aero_diam_edges", (/ dimid_aero_diam_edges /), unit="m", &
         long_name="aerosol radius axis bin edges", &
         description="logarithmically spaced edges of radius axis grid, " &
         // "with one more edge than center")
    call pmc_nc_write_real_1d(ncid, aero_diam_widths, &
         "aero_diam_widths", (/ dimid_aero_diam /), unit="m", &
         long_name="aerosol radius axis bin widths", &
         description="base-e logarithmic widths of radius axis grid, " &
         // "so that aero_diam_widths(i)" &
         // "= ln(aero_diam_edges(i+1) / aero_diam_edges(i)) and " &
         // "all bins have the same width")

  end subroutine bin_grid_netcdf_dim_aero_diam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine bin_grid_output_netcdf(bin_grid, ncid)
    
    !> bin_grid to write.
    type(bin_grid_t), intent(in) :: bin_grid
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_diam

    !> \page output_format_bin_grid Output File Format: Bin Grid Data
    !!
    !! The bin grid data NetCDF dimensions are:
    !!   - \b aero_diam: number of bins (grid cells) on the diameter axis
    !!   - \b aero_diam_edges: number of bin edges (grid cell edges) on
    !!     the diameter axis --- always equal to <tt>aero_diam + 1</tt>
    !!
    !! The bin grid data NetCDF variables are:
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
    !!   - \ref input_format_bin_grid --- the corresponding input format

    call bin_grid_netcdf_dim_aero_diam(bin_grid, ncid, &
         dimid_aero_diam)

    ! no need to write any more data as it's all contained in the
    ! dimension and associated variables

  end subroutine bin_grid_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine bin_grid_input_netcdf(bin_grid, ncid)
    
    !> bin_grid to read.
    type(bin_grid_t), intent(inout) :: bin_grid
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_diam, n_bin
    character(len=1000) :: name
    real(kind=dp), allocatable :: aero_diam_centers(:)
    real(kind=dp), allocatable :: aero_diam_edges(:)
    real(kind=dp), allocatable :: aero_diam_widths(:)

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_diam", dimid_aero_diam))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_diam, name, &
         n_bin))

    call bin_grid_deallocate(bin_grid)
    call bin_grid_allocate_size(bin_grid, n_bin)

    allocate(aero_diam_centers(n_bin))
    allocate(aero_diam_edges(n_bin + 1))
    allocate(aero_diam_widths(n_bin))

    call pmc_nc_read_real_1d(ncid, aero_diam_centers, &
         "aero_diam")
    call pmc_nc_read_real_1d(ncid, aero_diam_edges, &
         "aero_diam_edges")
    call pmc_nc_read_real_1d(ncid, aero_diam_widths, &
         "aero_diam_widths")

    bin_grid%center_radius = diam2rad(aero_diam_centers)
    bin_grid%edge_radius = diam2rad(aero_diam_edges)
    bin_grid%log_width = aero_diam_widths(1)

    deallocate(aero_diam_centers)
    deallocate(aero_diam_edges)
    deallocate(aero_diam_widths)

  end subroutine bin_grid_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_bin_grid
