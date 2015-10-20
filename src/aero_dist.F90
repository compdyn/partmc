! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_dist module.

!> The aero_dist_t structure and associated subroutines.
!!
!! The initial size distributions are computed as number densities, so
!! they can be used for both sectional and particle-resolved
!! simulations. The routine dist_to_n() converts a number concentration
!! distribution to an actual number of particles ready for a
!! particle-resolved simulation.
!!
!! Initial distributions should be normalized so that <tt>sum(n_den) =
!! 1/log_width</tt>.
module pmc_aero_dist

  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_spec_file
  use pmc_aero_data
  use pmc_aero_mode
  use pmc_mpi
  use pmc_rand
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> A complete aerosol distribution, consisting of several modes.
  type aero_dist_t
     !> Internally mixed modes.
     type(aero_mode_t), allocatable :: mode(:)
  end type aero_dist_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the number of modes.
  elemental integer function aero_dist_n_mode(aero_dist)

    !> Aero data structure.
    type(aero_dist_t), intent(in) :: aero_dist

    if (allocated(aero_dist%mode)) then
       aero_dist_n_mode = size(aero_dist%mode)
    else
       aero_dist_n_mode = 0
    end if

  end function aero_dist_n_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number concentration of a distribution. (#/m^3)
  real(kind=dp) function aero_dist_total_num_conc(aero_dist)

    !> Aerosol distribution.
    type(aero_dist_t), intent(in) :: aero_dist

    integer :: i_mode

    aero_dist_total_num_conc = 0d0
    do i_mode = 1,aero_dist_n_mode(aero_dist)
       aero_dist_total_num_conc = aero_dist_total_num_conc &
            + aero_mode_total_num_conc(aero_dist%mode(i_mode))
    end do

  end function aero_dist_total_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of particles of a distribution.
  real(kind=dp) function aero_dist_number(aero_dist, aero_weight)

    !> Aerosol distribution.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Aerosol weight.
    type(aero_weight_t), intent(in) :: aero_weight

    integer :: i_mode

    aero_dist_number = 0d0
    do i_mode = 1,aero_dist_n_mode(aero_dist)
       aero_dist_number = aero_dist_number &
            + aero_mode_number(aero_dist%mode(i_mode), aero_weight)
    end do

  end function aero_dist_number

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the binned number concentration for an aero_dist.
  subroutine aero_dist_num_conc(aero_dist, bin_grid, aero_data, &
       num_conc)

    !> Aero dist for which to compute number concentration.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number concentration (#(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: num_conc(bin_grid_size(bin_grid))

    integer :: i_mode
    real(kind=dp) :: mode_num_conc(size(num_conc, 1))

    num_conc = 0d0
    do i_mode = 1,aero_dist_n_mode(aero_dist)
       call aero_mode_num_conc(aero_dist%mode(i_mode), bin_grid, &
            aero_data, mode_num_conc)
       num_conc = num_conc + mode_num_conc
    end do

  end subroutine aero_dist_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the binned per-species volume concentration for an
  !> aero_dist.
  subroutine aero_dist_vol_conc(aero_dist, bin_grid, aero_data, &
       vol_conc)

    !> Aero dist for which to compute volume concentration.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration (V(ln(r))d(ln(r))).
    real(kind=dp), intent(out) :: vol_conc(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    integer :: i_mode
    real(kind=dp) :: mode_vol_conc(size(vol_conc, 1), size(vol_conc, 2))

    vol_conc = 0d0
    do i_mode = 1,aero_dist_n_mode(aero_dist)
       call aero_mode_vol_conc(aero_dist%mode(i_mode), bin_grid, &
            aero_data, mode_vol_conc)
       vol_conc = vol_conc + mode_vol_conc
    end do

  end subroutine aero_dist_vol_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether any of the modes are of the given type.
  elemental logical function aero_dist_contains_aero_mode_type(aero_dist, &
       aero_mode_type)

    !> Aerosol distribution.
    type(aero_dist_t), intent(in) :: aero_dist
    !> Aerosol mode type to test for.
    integer, intent(in) :: aero_mode_type

    integer :: i_mode

    aero_dist_contains_aero_mode_type = .false.
    do i_mode = 1,aero_dist_n_mode(aero_dist)
       aero_dist_contains_aero_mode_type = aero_dist_contains_aero_mode_type &
            .or. (aero_dist%mode(i_mode)%type == aero_mode_type)
    end do

  end function aero_dist_contains_aero_mode_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the current aero_dist and rate by interpolating at the
  !> current time with the lists of aero_dists and rates.
  subroutine aero_dist_interp_1d(aero_dist_list, time_list, &
         rate_list, time, aero_dist, rate)

    !> Gas states.
    type(aero_dist_t), intent(in) :: aero_dist_list(:)
    !> Times (s).
    real(kind=dp), intent(in) :: time_list(size(aero_dist_list))
    !> Rates (s^{-1}).
    real(kind=dp), intent(in) :: rate_list(size(aero_dist_list))
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current gas state.
    type(aero_dist_t), intent(inout) :: aero_dist
    !> Current rate (s^{-1}).
    real(kind=dp), intent(out) :: rate

    integer :: n, p, n_bin, n_spec, i, i_new
    real(kind=dp) :: y, alpha

    n = size(aero_dist_list)
    p = find_1d(n, time_list, time)
    if (p == 0) then
       ! before the start, just use the first state and rate
       aero_dist = aero_dist_list(1)
       rate = rate_list(1)
    elseif (p == n) then
       ! after the end, just use the last state and rate
       aero_dist = aero_dist_list(n)
       rate = rate_list(n)
    else
       ! in the middle, use the previous dist
       aero_dist = aero_dist_list(p)
       rate = rate_list(p)
    end if

  end subroutine aero_dist_interp_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read continuous aerosol distribution composed of several modes.
  subroutine spec_file_read_aero_dist(file, aero_data, aero_dist)

    !> Spec file to read data from.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Aerosol dist.
    type(aero_dist_t), intent(inout) :: aero_dist

    type(aero_mode_t) :: aero_mode
    type(aero_mode_t), allocatable :: new_modes(:)
    integer :: n, i
    logical :: eof

    ! note that the <p> is needed below to force correct paragraph
    ! breaking by doxygen

    !> \page input_format_aero_dist Input File Format: Aerosol Distribution
    !!
    !! <p>An aerosol distribution file consists of zero or more modes,
    !! each in the format described by \subpage input_format_aero_mode
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aero_mode --- the format for each mode
    !!     of an aerosol distribution

    allocate(aero_dist%mode(0))
    call spec_file_read_aero_mode(file, aero_data, aero_mode, eof)
    do while (.not. eof)
       n = size(aero_dist%mode)
       allocate(new_modes(n + 1))
       do i = 1,n
          new_modes(i) = aero_dist%mode(i)
       end do
       new_modes(n + 1) = aero_mode
       ! Next line is commented out due to a reported memory leak with valgrind
       ! using both gcc 4.6.3 and 4.9.2.
       !call move_alloc(new_modes, aero_dist%mode)
       ! Next two lines are here to avoid issues with reallocation of lhs by
       ! assignment.
       deallocate(aero_dist%mode)
       allocate(aero_dist%mode(n+1))
       aero_dist%mode = new_modes
       deallocate(new_modes)
       ! Next line appears to work in gcc 4.9.2 but led to problems with older
       ! versions of gcc (4.6.3). This will be ideal in the future as it avoids
       ! the copying of the aero_dist%mode array each loop.
       !aero_dist%mode = [aero_dist%mode, aero_mode]
       call spec_file_read_aero_mode(file, aero_data, aero_mode, eof)
    end do

  end subroutine spec_file_read_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read an array of aero_dists with associated times and rates from
  !> the given file.
  subroutine spec_file_read_aero_dists_times_rates(file, aero_data, &
       times, rates, aero_dists)

    !> Spec file to read data from.
    type(spec_file_t), intent(inout) :: file
    !> Aero data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Times (s).
    real(kind=dp), allocatable :: times(:)
    !> Rates (s^{-1}).
    real(kind=dp), allocatable :: rates(:)
    !> Aero dists.
    type(aero_dist_t), allocatable :: aero_dists(:)

    type(spec_line_t) :: aero_dist_line
    type(spec_file_t) :: aero_dist_file
    integer :: n_time, i_time
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: names(:)
    real(kind=dp), allocatable :: data(:,:)

    !> \page input_format_aero_dist_profile Input File Format: Aerosol Distribution Profile
    !!
    !! An aerosol distribution profile input file must consist of
    !! three lines:
    !! - the first line must begin with \c time and should be followed
    !!   by \f$N\f$ space-separated real scalars, giving the times (in
    !!   s after the start of the simulation) of the aerosol
    !!   distrbution set points --- the times must be in increasing
    !!   order
    !! - the second line must begin with \c rate and should be
    !!   followed by \f$N\f$ space-separated real scalars, giving the
    !!   values at the corresponding times
    !! - the third line must begin with \c dist and should be followed
    !!   by \f$N\f$ space-separated filenames, each specifying an
    !!   aerosol distribution in the format \ref input_format_aero_dist
    !!   at the corresponding time
    !!
    !! The units of the \c rate line depend on the type of aerosol
    !! distribution profile:
    !! - Emissions aerosol profiles have rates with units m/s ---
    !!   the aerosol distribution number concentrations are multiplied
    !!   by the rate to give an emission rate with unit #/(m^2 s)
    !!   which is then divided by the current mixing layer height
    !!   to give a per-volume emission rate.
    !! - Background aerosol profiles have rates with units \f$\rm
    !!   s^{-1}\f$, which is the dilution rate between the background
    !!   and the simulated air parcel. That is, if the simulated
    !!   number concentration is \f$N\f$ and the background number
    !!   concentration is \f$N_{\rm back}\f$, then dilution is modeled
    !!   as \f$\dot{N} = r N_{\rm back} - r N\f$, where \f$r\f$ is the
    !!   rate.
    !!
    !! Between the specified times the aerosol profile is interpolated
    !! step-wise and kept constant at its last value. That is, if the
    !! times are \f$t_i\f$, the rates are \f$r_i\f$, and the aerosol
    !! distributions are \f$a_i\f$ (all with \f$i = 1,\ldots,n\f$),
    !! then between times \f$t_i\f$ and \f$t_{i+1}\f$ the aerosol
    !! state is constant at \f$r_i a_i\f$. Before time \f$t_1\f$ the
    !! aerosol state is \f$r_1 a_1\f$, while after time \f$t_n\f$ it
    !! is \f$r_n a_n\f$.
    !!
    !! Example: an emissions aerosol profile could be:
    !! <pre>
    !! time  0          600        1800       # time (in s) after sim start
    !! rate  1          0.5        1          # scaling factor in m/s
    !! dist  dist1.dat  dist2.dat  dist3.dat  # aerosol distribution files
    !! </pre>
    !! Here the emissions between 0&nbsp;min and 10&nbsp;min are given
    !! by <tt>dist1.dat</tt> (with the number concentration
    !! interpreted as having units 1/(m^2 s)), the emissions between
    !! 10&nbsp;min and 30&nbsp;min are given by <tt>dist2.dat</tt>
    !! (scaled by 0.5), while the emissions after 30&nbsp;min are
    !! given by <tt>dist3.dat</tt>.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aero_data --- the aerosol species list
    !!     and material data
    !!   - \ref input_format_aero_dist --- the format of the
    !!     instantaneous aerosol distribution files

    ! read the data from the file
    call spec_file_read_real_named_array(file, 2, names, data)
    call spec_file_read_line_no_eof(file, aero_dist_line)
    call spec_file_check_line_name(file, aero_dist_line, "dist")
    call spec_file_check_line_length(file, aero_dist_line, size(data, 2))

    ! check the data size
    if (size(names) /= 2) then
       call die_msg(530745700, &
            trim(file%name) // ' must contain exactly two data lines')
    end if
    if (trim(names(1)) /= 'time') then
       call die_msg(570205795, 'row 1 in ' // trim(file%name) &
            // ' must start with: time not: ' // trim(names(1)))
    end if
    if (trim(names(2)) /= 'rate') then
       call die_msg(221270915, 'row 2 in ' // trim(file%name) &
            // ' must start with: rate not: ' // trim(names(1)))
    end if
    n_time = size(data, 2)
    if (n_time < 1) then
       call die_msg(457229710, 'each line in ' // trim(file%name) &
            // ' must contain at least one data value')
    end if

    ! copy over the data
    times = data(1,:)
    rates = data(2,:)
    if (allocated(aero_dists)) deallocate(aero_dists)
    allocate(aero_dists(n_time))
    do i_time = 1,n_time
       call spec_file_open(aero_dist_line%data(i_time), aero_dist_file)
       call spec_file_read_aero_dist(aero_dist_file, &
            aero_data, aero_dists(i_time))
       call spec_file_close(aero_dist_file)
    end do

  end subroutine spec_file_read_aero_dists_times_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_dist(val)

    !> Value to pack.
    type(aero_dist_t), intent(in) :: val

    integer :: i, total_size

    if (allocated(val%mode)) then
       total_size = pmc_mpi_pack_size_integer(aero_dist_n_mode(val))
       do i = 1,size(val%mode)
          total_size = total_size + pmc_mpi_pack_size_aero_mode(val%mode(i))
       end do
    else
       total_size = pmc_mpi_pack_size_integer(-1)
    end if
    pmc_mpi_pack_size_aero_dist = total_size

  end function pmc_mpi_pack_size_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_dist(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_dist_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    if (allocated(val%mode)) then
       call pmc_mpi_pack_integer(buffer, position, aero_dist_n_mode(val))
       do i = 1,size(val%mode)
          call pmc_mpi_pack_aero_mode(buffer, position, val%mode(i))
       end do
    else
       call pmc_mpi_pack_integer(buffer, position, -1)
    end if
    call assert(440557910, &
         position - prev_position <= pmc_mpi_pack_size_aero_dist(val))
#endif

  end subroutine pmc_mpi_pack_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_dist(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_dist_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    if (allocated(val%mode)) deallocate(val%mode)
    if (n >= 0) then
       allocate(val%mode(n))
       do i = 1,n
          call pmc_mpi_unpack_aero_mode(buffer, position, val%mode(i))
       end do
    end if
    call assert(742535268, &
         position - prev_position <= pmc_mpi_pack_size_aero_dist(val))
#endif

  end subroutine pmc_mpi_unpack_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_dist
