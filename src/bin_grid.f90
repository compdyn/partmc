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
     !> Len n_bin, bin center volumes (m^3).
     real(kind=dp), pointer :: v(:)
     !> Bin scale factor (1).
     real(kind=dp) :: dlnr
  end type bin_grid_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a bin_grid.
  subroutine bin_grid_allocate(bin_grid)

    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid

    bin_grid%n_bin = 0
    allocate(bin_grid%v(0))

  end subroutine bin_grid_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a bin_grid of the given size.
  subroutine bin_grid_allocate_size(bin_grid, n_bin)

    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Number of bins.
    integer, intent(in) :: n_bin

    bin_grid%n_bin = n_bin
    allocate(bin_grid%v(n_bin))

  end subroutine bin_grid_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all memory.
  subroutine bin_grid_deallocate(bin_grid)

    !> Bin_grid to free.
    type(bin_grid_t), intent(inout) :: bin_grid

    deallocate(bin_grid%v)

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
  subroutine bin_grid_make(bin_grid, n_bin, v_min, v_max)
    
    !> New bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Minimum volume (m^3).
    real(kind=dp), intent(in) :: v_min
    !> Minimum volume (m^3).
    real(kind=dp), intent(in) :: v_max

    call bin_grid_deallocate(bin_grid)
    call bin_grid_allocate_size(bin_grid, n_bin)
    call logspace(v_min, v_max, n_bin, bin_grid%v)
    ! dlnr = ln(r(i) / r(i-1))
    bin_grid%dlnr = log(vol2rad(v_max) / vol2rad(v_min)) / real(n_bin - 1, kind=dp)

  end subroutine bin_grid_make

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Given a bin_grid (which stores the center points of the bins),
  !> find the given edge volume (m^3).
  !!
  !! With n_bin bin centers there are (n_bin + 1) bin edges, so bin
  !! center bin_grid%v(i) is between bin edges i and (i + 1). This
  !! code currently assumes a logarithmically spaced bin grid and
  !! returns logarithmically spaced edges.
  real(kind=dp) function bin_grid_edge(bin_grid, i)
    
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Edge number (1 <= i <= n_bin + 1).
    integer, intent(in) :: i

    real(kind=dp) :: log_v_min, log_v_max, log_delta

    call assert(440393735, bin_grid%n_bin > 1)
    log_v_min = log(bin_grid%v(1))
    log_v_max = log(bin_grid%v(bin_grid%n_bin))
    log_delta = (log_v_max - log_v_min) / real(bin_grid%n_bin - 1, kind=dp)
    bin_grid_edge = exp(log_v_min + (real(i, kind=dp) - 1.5d0) * log_delta)
    
  end function bin_grid_edge
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the bin number that contains a given particle. This assumes
  !> logarithmically spaced bins.
  integer function bin_grid_particle_in_bin(bin_grid, v)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume of particle.
    real(kind=dp), intent(in) :: v

    real(kind=dp) :: log_v_min, log_v_max, log_edge_min, log_edge_max
    real(kind=dp) :: half_log_delta
    integer :: k

    call assert(448215689, bin_grid%n_bin > 2)
    log_v_min = log(bin_grid%v(1))
    log_v_max = log(bin_grid%v(bin_grid%n_bin))
    half_log_delta = (log_v_max - log_v_min) / real(2 * (bin_grid%n_bin - 1), kind=dp)
    log_edge_min = log_v_min + half_log_delta
    log_edge_max = log_v_max - half_log_delta
    k = ceiling((log(v) - log_edge_min) / (log_edge_max - log_edge_min) &
         * real(bin_grid%n_bin - 2, kind=dp)) + 1
    k = max(k, 1)
    k = min(k, bin_grid%n_bin)
    bin_grid_particle_in_bin = k
    !FIXME: above should be equivalent to:
    !    i_bin = ceiling((log(radius) - log(r_min)) &
    !         / (log(r_max) - log(r_min)) * real(n_bin - 1, kind=dp) + 0.5d0)
    
  end function bin_grid_particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a bin_grid from a spec file and
  !> generate it.
  subroutine spec_file_read_bin_grid(file, bin_grid)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid

    integer :: n_bin
    real(kind=dp) :: r_min, r_max

    call spec_file_read_integer(file, 'n_bin', n_bin)
    call spec_file_read_real(file, 'r_min', r_min)
    call spec_file_read_real(file, 'r_max', r_max)
    call bin_grid_make(bin_grid, n_bin, rad2vol(r_min), rad2vol(r_max))

  end subroutine spec_file_read_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_bin_grid(val)

    !> Value to pack.
    type(bin_grid_t), intent(in) :: val

    pmc_mpi_pack_size_bin_grid = &
         pmc_mpi_pack_size_integer(val%n_bin) &
         + pmc_mpi_pack_size_real_array(val%v) &
         + pmc_mpi_pack_size_real(val%dlnr)

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
    call pmc_mpi_pack_real_array(buffer, position, val%v)
    call pmc_mpi_pack_real(buffer, position, val%dlnr)
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
    type(bin_grid_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_bin)
    call pmc_mpi_unpack_real_array(buffer, position, val%v)
    call pmc_mpi_unpack_real(buffer, position, val%dlnr)
    call assert(741838730, &
         position - prev_position == pmc_mpi_pack_size_bin_grid(val))
#endif

  end subroutine pmc_mpi_unpack_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the aero_radius dimension to the given NetCDF file if it is
  !> not already present and in any case return the associated dimid.
  subroutine bin_grid_netcdf_dim_aero_radius(bin_grid, ncid, &
       dimid_aero_radius)

    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the aero_radius dimension.
    integer, intent(out) :: dimid_aero_radius

    integer :: status, i_bin, varid_aero_radius
    integer :: dimid_aero_radius_edges, varid_aero_radius_edges, &
         varid_aero_radius_widths
    real(kind=dp) :: aero_radius_centers(bin_grid%n_bin), &
         aero_radius_edges(bin_grid%n_bin + 1)
    real(kind=dp) :: aero_radius_widths(bin_grid%n_bin)

    status = nf90_inq_dimid(ncid, "aero_radius", dimid_aero_radius)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_radius", &
         bin_grid%n_bin, dimid_aero_radius))
    call pmc_nc_check(nf90_def_dim(ncid, "aero_radius_edges", &
         bin_grid%n_bin + 1, dimid_aero_radius_edges))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_bin = 1,bin_grid%n_bin
       aero_radius_centers(i_bin) = vol2rad(bin_grid%v(i_bin))
       aero_radius_widths(i_bin) = bin_grid%dlnr
    end do
    do i_bin = 1,(bin_grid%n_bin + 1)
       aero_radius_edges(i_bin) = vol2rad(bin_grid_edge(bin_grid, i_bin))
    end do

    call pmc_nc_write_real_1d(ncid, aero_radius_centers, &
         "aero_radius", (/ dimid_aero_radius /), unit="m", &
         long_name="aerosol radius axis bin centers", &
         description="logarithmically spaced centers of radius axis grid, " &
         // "so that aero_radius(i) / aero_radius_edges(i) = " &
         // "0.5 * aero_radius_edges(i+1) / aero_radius_edges(i)")
    call pmc_nc_write_real_1d(ncid, aero_radius_edges, &
         "aero_radius_edges", (/ dimid_aero_radius_edges /), unit="m", &
         long_name="aerosol radius axis bin edges", &
         description="logarithmically spaced edges of radius axis grid, " &
         // "with one more edge than center")
    call pmc_nc_write_real_1d(ncid, aero_radius_widths, &
         "aero_radius_widths", (/ dimid_aero_radius /), unit="m", &
         long_name="aerosol radius axis bin widths", &
         description="base-e logarithmic widths of radius axis grid, " &
         // "so that aero_radius_widths(i)" &
         // "= ln(aero_radius_edges(i+1) / aero_radius_edges(i)) and " &
         // "all bins have the same width")

  end subroutine bin_grid_netcdf_dim_aero_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine bin_grid_output_netcdf(bin_grid, ncid)
    
    !> bin_grid to write.
    type(bin_grid_t), intent(in) :: bin_grid
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_radius

    !> \page output_format_bin_grid Output File Format: Bin Grid Data
    !!
    !! The bin grid data NetCDF dimensions are:
    !!   - \b aero_radius: number of bins (grid cells) on the radius axis
    !!   - \b aero_radius_edges: number of bin edges (grid cell edges) on
    !!     the radius axis --- always equal to <tt>aero_radius + 1</tt>
    !!
    !! The bin grid data NetCDF variables are:
    !!   - \b aero_radius (unit m, dim \c aero_radius): aerosol radius axis
    !!     bin centers --- centered on a logarithmic scale from the edges, so
    !!     that <tt>aero_radius(i) / aero_radius_edges(i) =
    !!     0.5 * aero_radius_edges(i+1) / aero_radius_edges(i)</tt>
    !!   - \b aero_radius_edges (unit m, dim \c aero_radius_edges): aersol
    !!     radius axis bin edges (there is one more edge than center)
    !!   - \b aero_radius_widths (dimensionless, dim \c aero_radius):
    !!     the base-e logarithmic bin widths --- <tt>aero_radius_widths(i)
    !!     = ln(aero_radius_edges(i+1) / aero_radius_edges(i))</tt>, so
    !!     all bins have the same width

    call bin_grid_netcdf_dim_aero_radius(bin_grid, ncid, &
         dimid_aero_radius)

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

    integer :: dimid_aero_radius
    character(len=1000) :: unit, name
    integer :: n_bin, i_bin
    real(kind=dp), allocatable :: aero_radius_centers(:)
    real(kind=dp), allocatable :: aero_radius_widths(:)

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_radius", dimid_aero_radius))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_radius, name, &
         n_bin))

    call bin_grid_deallocate(bin_grid)
    call bin_grid_allocate_size(bin_grid, n_bin)

    allocate(aero_radius_centers(n_bin))
    allocate(aero_radius_widths(n_bin))

    call pmc_nc_read_real_1d(ncid, aero_radius_centers, &
         "aero_radius", unit)
    call pmc_nc_read_real_1d(ncid, aero_radius_widths, &
         "aero_radius_widths", unit)

    do i_bin = 1,n_bin
       bin_grid%v(i_bin) = rad2vol(aero_radius_centers(i_bin))
    end do
    bin_grid%dlnr = aero_radius_widths(1)

    deallocate(aero_radius_centers)
    deallocate(aero_radius_widths)

  end subroutine bin_grid_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_bin_grid
