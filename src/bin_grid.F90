! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
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

  !> Invalid type of bin grid.
  integer, parameter :: BIN_GRID_TYPE_INVALID = 0
  !> Logarithmically spaced bin grid.
  integer, parameter :: BIN_GRID_TYPE_LOG = 1
  !> Linearly spaced bin grid.
  integer, parameter :: BIN_GRID_TYPE_LINEAR = 2

  !> 1D grid, either logarithmic or linear.
  !!
  !! The grid of bins is logarithmically spaced in volume, an
  !! assumption that is quite heavily incorporated into the code. At
  !! some point in the future it would be nice to relax this
  !! assumption.
  type bin_grid_t
     !> Type of grid spacing (BIN_GRID_TYPE_LOG, etc).
     integer :: type
     !> Number of bins.
     integer :: n_bin
     !> Minimum bin edge.
     real(kind=dp) :: min
     !> Maximum bin edge.
     real(kind=dp) :: max
  end type bin_grid_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a bin_grid.
  subroutine bin_grid_allocate(bin_grid)

    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid

    call bin_grid_zero(bin_grid)

  end subroutine bin_grid_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all memory.
  subroutine bin_grid_deallocate(bin_grid)

    !> Bin_grid to free.
    type(bin_grid_t), intent(inout) :: bin_grid

  end subroutine bin_grid_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets a bin grid to a zero state.
  subroutine bin_grid_zero(bin_grid)

    !> Bin_grid to zero.
    type(bin_grid_t), intent(inout) :: bin_grid

    bin_grid%type = BIN_GRID_TYPE_INVALID
    bin_grid%n_bin = 0
    bin_grid%min = 0d0
    bin_grid%max = 0d0

  end subroutine bin_grid_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copies a bin grid.
  subroutine bin_grid_copy(bin_grid_from, bin_grid_to)

    !> Bin_grid to copy from.
    type(bin_grid_t), intent(in) :: bin_grid_from
    !> Bin_grid to copy to.
    type(bin_grid_t), intent(inout) :: bin_grid_to

    bin_grid_to%type = bin_grid_from%type
    bin_grid_to%n_bin = bin_grid_from%n_bin
    bin_grid_to%min = bin_grid_from%min
    bin_grid_to%max = bin_grid_from%max

  end subroutine bin_grid_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the center of a bin in a grid.
  real(kind=dp) function bin_grid_center(bin_grid, i_bin)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Bin number to return center of.
    integer, intent(in) :: i_bin

    call assert_msg(550416095, bin_grid%n_bin > 0, &
         "cannot locate bin center with non-positive n_bin")
    call assert_msg(614416072, &
         (i_bin >= 1) .and. (i_bin <= bin_grid%n_bin), &
         "invalid bin number " // integer_to_string(i_bin) &
         // " in grid with " // integer_to_string(bin_grid%n_bin) // " bins")

    if (bin_grid%type == BIN_GRID_TYPE_LOG) then
       bin_grid_center = exp(interp_linear_disc(log(bin_grid%min), &
            log(bin_grid%max), 2 * bin_grid%n_bin + 1, 2 * i_bin))
    elseif (bin_grid%type == BIN_GRID_TYPE_LOG) then
       bin_grid_center = interp_linear_disc(bin_grid%min, bin_grid%max, &
            2 * bin_grid%n_bin + 1, 2 * i_bin)
    else
       call die_msg(801471019, "unknown bin_grid type: " &
            // integer_to_string(bin_grid%type))
    end if

  end function bin_grid_center

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the edge of a bin in a grid.
  !!
  !! The edge number \c i_edge must be in the range 1,...,<tt>(n_bin + 1)</tt>.
  !! Bin number \c i has edges numbers \c i and <tt>i + 1</tt>.
  real(kind=dp) function bin_grid_edge(bin_grid, i_edge)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Edge number to return (1 to <tt>n_bin + 1</tt>).
    integer, intent(in) :: i_edge

    call assert_msg(677240492, bin_grid%n_bin > 0, &
         "cannot locate bin edge with non-positive n_bin")
    call assert_msg(905667628, &
         (i_edge >= 1) .and. (i_edge <= bin_grid%n_bin + 1), &
         "invalid edge number " // integer_to_string(i_edge) &
         // " in grid with " // integer_to_string(bin_grid%n_bin) // " bins")
    if (i_edge == 1) then
       bin_grid_edge = bin_grid%min
       return
    elseif (i_edge == bin_grid%n_bin + 1) then
       bin_grid_edge = bin_grid%max
       return
    end if

    if (bin_grid%type == BIN_GRID_TYPE_LOG) then
       bin_grid_edge = exp(interp_linear_disc(log(bin_grid%min), &
            log(bin_grid%max), bin_grid%n_bin + 1, i_edge))
    elseif (bin_grid%type == BIN_GRID_TYPE_LOG) then
       bin_grid_edge = interp_linear_disc(bin_grid%min, bin_grid%max, &
            bin_grid%n_bin + 1, i_edge)
    else
       call die_msg(749216544, "unknown bin_grid type: " &
            // integer_to_string(bin_grid%type))
    end if

  end function bin_grid_edge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the width of a bin in a grid.
  real(kind=dp) function bin_grid_width(bin_grid, i_bin)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Bin number to return size of.
    integer, intent(in) :: i_bin

    call assert_msg(604053727, bin_grid%n_bin > 0, &
         "cannot locate bin width with non-positive n_bin")
    call assert_msg(683341327, &
         (i_bin >= 1) .and. (i_bin <= bin_grid%n_bin), &
         "invalid bin number " // integer_to_string(i_bin) &
         // " in grid with " // integer_to_string(bin_grid%n_bin) // " bins")

    if (bin_grid%type == BIN_GRID_TYPE_LOG) then
       bin_grid_width = (log(bin_grid%max) - log(bin_grid%min)) &
            / bin_grid%n_bin
    elseif (bin_grid%type == BIN_GRID_TYPE_LOG) then
       bin_grid_width = (bin_grid%max - bin_grid%min) / bin_grid%n_bin
    else
       call die_msg(801471019, "unknown bin_grid type: " &
            // integer_to_string(bin_grid%type))
    end if

  end function bin_grid_width

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an array of the centers of a bin in a grid.
  function bin_grid_center_array(bin_grid)

    real(kind=dp), allocatable :: bin_grid_center_array(:)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

    integer :: i

    allocate(bin_grid_center_array(bin_grid%n_bin))
    do i = 1,bin_grid%n_bin
       bin_grid_center_array(i) = bin_grid_center(bin_grid, i)
    end do

  end function bin_grid_center_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an array of the widths of a bin in a grid.
  function bin_grid_edge_array(bin_grid)

    real(kind=dp), allocatable :: bin_grid_edge_array(:)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

    integer :: i

    allocate(bin_grid_edge_array(bin_grid%n_bin + 1))
    do i = 1,bin_grid%n_bin + 1
       bin_grid_edge_array(i) = bin_grid_edge(bin_grid, i)
    end do

  end function bin_grid_edge_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns an array of the widths of a bin in a grid.
  function bin_grid_width_array(bin_grid)

    real(kind=dp), allocatable :: bin_grid_width_array(:)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid

    integer :: i

    allocate(bin_grid_width_array(bin_grid%n_bin))
    do i = 1,bin_grid%n_bin
       bin_grid_width_array(i) = bin_grid_width(bin_grid, i)
    end do

  end function bin_grid_width_array

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
  subroutine bin_grid_make(bin_grid, type, n_bin, min, max)
    
    !> New bin grid.
    type(bin_grid_t), intent(inout) :: bin_grid
    !> Type of bin grid.
    integer, intent(in) :: type
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Minimum edge value.
    real(kind=dp), intent(in) :: min
    !> Minimum edge value.
    real(kind=dp), intent(in) :: max

    call assert_msg(538534122, n_bin >= 0, &
         "bin_grid requires a non-negative n_bin, not: " &
         // trim(integer_to_string(n_bin)))
    if (n_bin > 0) then
       if (type == BIN_GRID_TYPE_LOG) then
          call assert_msg(966541762, min > 0d0, &
               "log bin_grid requires a positive min value, not: " &
               // trim(real_to_string(min)))
       end if
       call assert_msg(711537859, min < max, &
            "bin_grid requires min < max, not: " &
            // trim(real_to_string(min)) // " and " &
            // trim(real_to_string(max)))
    end if
    bin_grid%type = type
    bin_grid%n_bin = n_bin
    bin_grid%min = min
    bin_grid%max = max

  end subroutine bin_grid_make

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the bin number that contains a given value.
  !!
  !! If a particle is below the smallest bin then its bin number is
  !! 0. If a particle is above the largest bin then its bin number is
  !! <tt>n_bin + 1</tt>.
  integer function bin_grid_find(bin_grid, val)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Value to locate bin for.
    real(kind=dp), intent(in) :: val

    call assert(448215689, bin_grid%n_bin >= 0)
    if (bin_grid%n_bin == 0) then
       bin_grid_find = 0
       return
    end if
    if (bin_grid%type == BIN_GRID_TYPE_LOG) then
       bin_grid_find = logspace_find(bin_grid%min, bin_grid%max, &
            bin_grid%n_bin + 1, val)
    elseif (bin_grid%type == BIN_GRID_TYPE_LINEAR) then
       bin_grid_find = linspace_find(bin_grid%min, bin_grid%max, &
            bin_grid%n_bin + 1, val)
    else
       call die_msg(348908641, "unknown bin_grid type: " &
            // integer_to_string(bin_grid%type))
    end if

  end function bin_grid_find
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Make a histogram with of the given weighted data, scaled by the
  !> bin sizes.
  subroutine bin_grid_histogram_1d(x_bin_grid, x_data, weight_data, hist)

    !> x-axis bin grid.
    type(bin_grid_t), intent(in) :: x_bin_grid
    !> Data values on the x-axis.
    real(kind=dp), intent(in) :: x_data(:)
    !> Data value weights.
    real(kind=dp), intent(in) :: weight_data(size(x_data))
    !> Histogram to compute.
    real(kind=dp), intent(out) :: hist(x_bin_grid%n_bin)

    integer :: i_data, i_bin

    hist = 0d0
    do i_data = 1,size(x_data)
       i_bin = bin_grid_find(x_bin_grid, x_data(i_data))
       if ((i_bin >= 1) .and. (i_bin <= x_bin_grid%n_bin)) then
          hist(i_bin) = hist(i_bin) &
               + weight_data(i_data) / bin_grid_width(x_bin_grid, i_bin)
       end if
    end do

  end subroutine bin_grid_histogram_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a radius bin_grid from a spec file.
  subroutine spec_file_read_radius_bin_grid(file, bin_grid)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Radius bin grid.
    type(bin_grid_t), intent(inout) :: bin_grid

    integer :: n_bin
    real(kind=dp) :: d_min, d_max

    !> \page input_format_diam_bin_grid Input File Format: Diameter Axis Bin Grid
    !!
    !! The diameter bin grid is logarithmic, consisting of
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
    !!   - \ref output_format_diam_bin_grid --- the corresponding output format

    call spec_file_read_integer(file, 'n_bin', n_bin)
    call spec_file_read_real(file, 'd_min', d_min)
    call spec_file_read_real(file, 'd_max', d_max)
    call bin_grid_make(bin_grid, BIN_GRID_TYPE_LOG, n_bin, diam2rad(d_min), &
         diam2rad(d_max))

  end subroutine spec_file_read_radius_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_bin_grid(val)

    !> Value to pack.
    type(bin_grid_t), intent(in) :: val

    pmc_mpi_pack_size_bin_grid = &
         pmc_mpi_pack_size_integer(val%type) &
         + pmc_mpi_pack_size_integer(val%n_bin) &
         + pmc_mpi_pack_size_real(val%min) &
         + pmc_mpi_pack_size_real(val%max)

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
    call pmc_mpi_pack_integer(buffer, position, val%type)
    call pmc_mpi_pack_integer(buffer, position, val%n_bin)
    call pmc_mpi_pack_real(buffer, position, val%min)
    call pmc_mpi_pack_real(buffer, position, val%max)
    call assert(385455586, &
         position - prev_position <= pmc_mpi_pack_size_bin_grid(val))
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
    call pmc_mpi_unpack_integer(buffer, position, val%type)
    call pmc_mpi_unpack_integer(buffer, position, val%n_bin)
    call pmc_mpi_unpack_real(buffer, position, val%min)
    call pmc_mpi_unpack_real(buffer, position, val%max)
    call assert(741838730, &
         position - prev_position <= pmc_mpi_pack_size_bin_grid(val))
#endif

  end subroutine pmc_mpi_unpack_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether all processors have the same value.
  logical function pmc_mpi_allequal_bin_grid(val)

    !> Value to compare.
    type(bin_grid_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    if (.not. pmc_mpi_allequal_integer(val%type)) then
       pmc_mpi_allequal_bin_grid = .false.
       return
    end if

    if (.not. pmc_mpi_allequal_integer(val%n_bin)) then
       pmc_mpi_allequal_bin_grid = .false.
       return
    end if

    if (val%n_bin == 0) then
       pmc_mpi_allequal_bin_grid = .true.
       return
    end if

    if (pmc_mpi_allequal_real(val%min) &
         .and. pmc_mpi_allequal_real(val%max)) then
       pmc_mpi_allequal_bin_grid = .true.
    else
       pmc_mpi_allequal_bin_grid = .false.
    end if
#else
    pmc_mpi_allequal_bin_grid = .true.
#endif

  end function pmc_mpi_allequal_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a bin grid dimension to the given NetCDF file if it is
  !> not already present and in any case return the associated dimid.
  subroutine bin_grid_netcdf_dim(bin_grid, ncid, dim_name, long_name, &
       unit_name, dimid, scale)

    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimension name.
    character(len=*), intent(in) :: dim_name
    !> Long dimension name to use.
    character(len=*), intent(in) :: long_name
    !> Units for the grid.
    character(len=*), intent(in) :: unit_name
    !> Dimid of the grid dimension.
    integer, intent(out) :: dimid
    !> Factor to scale grid by before output.
    real(kind=dp), intent(in), optional :: scale

    integer :: status, varid, dimid_edges, varid_edges, varid_widths, i
    real(kind=dp) :: centers(bin_grid%n_bin), edges(bin_grid%n_bin + 1)
    real(kind=dp) :: widths(bin_grid%n_bin)
    character(len=(len_trim(dim_name)+10)) :: dim_name_edges

    status = nf90_inq_dimid(ncid, dim_name, dimid)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    dim_name_edges = trim(dim_name) // "_edges"

    call pmc_nc_check(nf90_def_dim(ncid, dim_name, bin_grid%n_bin, dimid))
    call pmc_nc_check(nf90_def_dim(ncid, dim_name_edges, bin_grid%n_bin + 1, &
         dimid_edges))

    call pmc_nc_check(nf90_enddef(ncid))

    do i = 1,bin_grid%n_bin
       centers(i) = bin_grid_center(bin_grid, i)
       edges(i) = bin_grid_edge(bin_grid, i)
       widths(i) = bin_grid_width(bin_grid, i)
    end do
    edges(bin_grid%n_bin + 1) = bin_grid_edge(bin_grid, bin_grid%n_bin + 1)

    if (bin_grid%type == BIN_GRID_TYPE_LOG) then
       if (present(scale)) then
          centers = centers * scale
          edges = edges * scale
       end if
       call pmc_nc_write_real_1d(ncid, centers, dim_name, (/ dimid /), &
            unit=unit_name, long_name=(trim(long_name) // " grid centers"), &
            description=("logarithmically spaced centers of " &
            // trim(long_name) // " grid, so that " // trim(dim_name) &
            // "(i) is the geometric mean of " // trim(dim_name_edges) &
            // "(i) and " // trim(dim_name_edges) // "(i + 1)"))
       call pmc_nc_write_real_1d(ncid, edges, dim_name_edges, &
            (/ dimid_edges /), unit=unit_name, &
            long_name=(trim(long_name) // " grid edges"), &
            description=("logarithmically spaced edges of " &
            // trim(long_name) // " grid, with one more edge than center"))
       call pmc_nc_write_real_1d(ncid, widths, trim(dim_name) // "_widths", &
            (/ dimid /), unit="1", &
            long_name=(trim(long_name) // " grid widths"), &
            description=("base-e logarithmic widths of " &
            // trim(long_name) // " grid, with " // trim(dim_name) &
            // "_widths(i) = ln(" // trim(dim_name_edges) // "(i + 1) / " &
            // trim(dim_name_edges) // "(i))"))
    elseif (bin_grid%type == BIN_GRID_TYPE_LINEAR) then
       if (present(scale)) then
          centers = centers * scale
          edges = edges * scale
          widths = widths * scale
       end if
       call pmc_nc_write_real_1d(ncid, centers, dim_name, (/ dimid /), &
            unit=unit_name, long_name=(trim(long_name) // " grid centers"), &
            description=("linearly spaced centers of " // trim(long_name) &
            // " grid, so that " // trim(dim_name) // "(i) is the mean of " &
            // trim(dim_name_edges) // "(i) and " // trim(dim_name_edges) &
            // "(i + 1)"))
       call pmc_nc_write_real_1d(ncid, edges, dim_name_edges, &
            (/ dimid_edges /), unit=unit_name, &
            long_name=(trim(long_name) // " grid edges"), &
            description=("linearly spaced edges of " &
            // trim(long_name) // " grid, with one more edge than center"))
       call pmc_nc_write_real_1d(ncid, widths, trim(dim_name) // "_widths", &
            (/ dimid /), unit=unit_name, &
            long_name=(trim(long_name) // " grid widths"), &
            description=("widths of " // trim(long_name) // " grid, with " &
            // trim(dim_name) // "_widths(i) = " // trim(dim_name_edges) &
            // "(i + 1) - " // trim(dim_name_edges) // "(i)"))
    else
       call die_msg(942560572, "unknown bin_grid type: " &
            // integer_to_string(bin_grid%type))
    end if

  end subroutine bin_grid_netcdf_dim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine bin_grid_input_netcdf(bin_grid, ncid, dim_name, scale)
    
    !> bin_grid to read.
    type(bin_grid_t), intent(inout) :: bin_grid
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimension name.
    character(len=*), intent(in) :: dim_name
    !> Factor to scale grid by after input.
    real(kind=dp), intent(in), optional :: scale

    integer :: dimid, varid, n_bin, type
    character(len=1000) :: name, description
    real(kind=dp), allocatable :: edges(:)

    call pmc_nc_check(nf90_inq_dimid(ncid, dim_name, dimid))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid, name, n_bin))
    call pmc_nc_check(nf90_inq_varid(ncid, dim_name, varid))
    call pmc_nc_check(nf90_get_att(ncid, varid, "description", description))

    allocate(edges(n_bin + 1))
    call pmc_nc_read_real_1d(ncid, edges, dim_name // "_edges")

    if (starts_with(description, "logarithmically")) then
       type = BIN_GRID_TYPE_LOG
    elseif (starts_with(description, "logarithmically")) then
       type = BIN_GRID_TYPE_LINEAR
    else
       call die_msg(792158584, "cannot identify grid type for NetCDF " &
            // "dimension: " // trim(dim_name))
    end if

    if (present(scale)) then
       call bin_grid_make(bin_grid, type, n_bin, scale * edges(1), &
            scale * edges(n_bin + 1))
    else
       call bin_grid_make(bin_grid, type, n_bin, edges(1), edges(n_bin + 1))
    end if

    deallocate(edges)

  end subroutine bin_grid_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_bin_grid
