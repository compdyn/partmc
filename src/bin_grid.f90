! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Functions that deal with the bin grid.
!
! The grid of bins is logarithmically spaced in volume, an assumption
! that is quite heavily incorporated into the code. At some point in
! the future it would be nice to relax this assumption.

module mod_bin_grid

  type bin_grid_t
     integer :: n_bin                   ! number of bins
     real*8, pointer :: v(:)            ! len n_bin, bin center volumes (m^3)
     real*8 :: dlnr                     ! bin scale factor (1)
  end type bin_grid_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_bin_grid(n_bin, bin_grid)

    ! Allocates a bin_grid.

    integer, intent(in) :: n_bin        ! number of bins
    type(bin_grid_t), intent(out) :: bin_grid ! bin grid

    bin_grid%n_bin = n_bin
    allocate(bin_grid%v(n_bin))

  end subroutine alloc_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_bin_grid(bin_grid)

    ! Frees all memory.

    type(bin_grid_t), intent(inout) :: bin_grid ! bin_grid to free

    deallocate(bin_grid%v)

  end subroutine free_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine vol_to_lnr(r, f_vol, f_lnr)
    
    ! Convert a density f(vol)d(vol) to f(ln(r))d(ln(r))
    ! where vol = 4/3 pi r^3.
    
    use mod_constants
    
    real*8, intent(in) :: r             ! radius (m)
    real*8, intent(in) :: f_vol         ! density as a function of volume
    real*8, intent(out) :: f_lnr        ! density as a function of ln(r)
    
    f_lnr = f_vol * 4d0 * const%pi * r**3
    
  end subroutine vol_to_lnr
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bin_grid_make(n_bin, v_min, v_max, bin_grid)

    ! Generates the bin grid given the range and number of bins.
    
    use mod_util

    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: v_min         ! minimum volume (m^3)
    real*8, intent(in) :: v_max         ! minimum volume (m^3)
    type(bin_grid_t), intent(out) :: bin_grid ! new bin grid, will be allocated

    call alloc_bin_grid(n_bin, bin_grid)
    call logspace(v_min, v_max, n_bin, bin_grid%v)
    ! dlnr = ln(r(i) / r(i-1))
    bin_grid%dlnr = log(vol2rad(v_max) / vol2rad(v_min)) / dble(n_bin - 1)

  end subroutine bin_grid_make

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bin_edge(bin_grid, i, v_edge)

    ! Given a bin_grid (which stores the center points of the bins),
    ! find the given edge volume. With n_bin bin centers there are
    ! (n_bin + 1) bin edges, so bin center bin_grid%v(i) is between
    ! bin edges i and (i + 1).
    
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    integer, intent(in) :: i            ! edge number (1 <= i <= n_bin + 1)
    real*8, intent(out) :: v_edge       ! volume at edge
    
    if (i .eq. 1) then
       v_edge = bin_grid%v(1) - (bin_grid%v(2) - bin_grid%v(1)) / 2d0
    elseif (i .eq. (bin_grid%n_bin + 1)) then
       v_edge = bin_grid%v(bin_grid%n_bin) &
            + (bin_grid%v(bin_grid%n_bin) &
            - bin_grid%v(bin_grid%n_bin - 1)) / 2d0
    else
       v_edge = (bin_grid%v(i - 1) + bin_grid%v(i)) / 2d0
    end if
    
  end subroutine bin_edge
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine particle_in_bin_fast(v, bin_grid, k)
    
    ! Find the bin number that contains a given particle. This assumes
    ! logarithmically spaced bins.

    use mod_util
    
    real*8, intent(in) :: v             ! volume of particle
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    integer, intent(out) :: k           ! bin number containing particle

    real*8 :: log_v_min, log_v_max, log_edge_min, log_edge_max
    real*8 :: half_log_delta

    call assert(bin_grid%n_bin > 2)
    log_v_min = log(bin_grid%v(1))
    log_v_max = log(bin_grid%v(bin_grid%n_bin))
    half_log_delta = (log_v_max - log_v_min) / dble(2 * (bin_grid%n_bin - 1))
    log_edge_min = log_v_min + half_log_delta
    log_edge_max = log_v_max - half_log_delta
    k = ceiling((log(v) - log_edge_min) / (log_edge_max - log_edge_min) &
         * dble(bin_grid%n_bin - 2)) + 1
    k = max(k, 1)
    k = min(k, bin_grid%n_bin)
    
  end subroutine particle_in_bin_fast
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine particle_in_bin(v, bin_grid, k)
    
    ! Find the bin number that contains a given particle.
    
    real*8, intent(in) :: v             ! volume of particle
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    integer, intent(out) :: k           ! bin number containing particle
    
    real*8 :: edge

    k = 0
300 k = k + 1
    if (k .lt. bin_grid%n_bin) then
       call bin_edge(bin_grid, k + 1, edge)
       if (v .gt. edge) then
          goto 300
       end if
    end if
    
  end subroutine particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_bin_grid(file, bin_grid)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid to write

    call inout_write_integer(file, "n_bin", bin_grid%n_bin)
    call inout_write_real_array(file, "center_volumes(m^3)", bin_grid%v)
    call inout_write_real(file, "dlnr", bin_grid%dlnr)
    
  end subroutine inout_write_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_bin_grid(file, bin_grid)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(bin_grid_t), intent(out) :: bin_grid ! bin_grid to read

    call inout_read_integer(file, "n_bin", bin_grid%n_bin)
    call inout_read_real_array(file, "center_volumes(m^3)", bin_grid%v)
    call inout_read_real(file, "dlnr", bin_grid%dlnr)
    
  end subroutine inout_read_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_bin_grid(file, bin_grid)

    ! Read the specification for a bin_grid from a inout file and
    ! generate it.

    use mod_inout
    use mod_util

    type(inout_file_t), intent(inout) :: file ! inout file
    type(bin_grid_t), intent(out) :: bin_grid ! bin grid

    integer :: n_bin
    real*8 :: r_min, r_max

    call inout_read_integer(file, 'n_bin', n_bin)
    call inout_read_real(file, 'r_min', r_min)
    call inout_read_real(file, 'r_max', r_max)
    call bin_grid_make(n_bin, rad2vol(r_min), rad2vol(r_max), bin_grid)

  end subroutine spec_read_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_bin_grid
