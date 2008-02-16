! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_bin_grid module.

!> The bin_grid_t structure and associated subroutines.
module pmc_bin_grid

  use pmc_constants
  use pmc_util
  use pmc_inout
  use pmc_mpi
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
     real*8, pointer :: v(:)
     !> Bin scale factor (1).
     real*8 :: dlnr
  end type bin_grid_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates a bin_grid.
  subroutine bin_grid_alloc(bin_grid, n_bin)

    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Number of bins.
    integer, intent(in) :: n_bin

    bin_grid%n_bin = n_bin
    allocate(bin_grid%v(n_bin))

  end subroutine bin_grid_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all memory.
  subroutine bin_grid_free(bin_grid)

    !> Bin_grid to free.
    type(bin_grid_t), intent(inout) :: bin_grid

    deallocate(bin_grid%v)

  end subroutine bin_grid_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a density f(vol)d(vol) to f(ln(r))d(ln(r))
  !> where vol = 4/3 pi r^3.
  subroutine vol_to_lnr(r, f_vol, f_lnr)
    
    !> Radius (m).
    real*8, intent(in) :: r
    !> Density as a function of volume.
    real*8, intent(in) :: f_vol
    !> Density as a function of ln(r).
    real*8, intent(out) :: f_lnr
    
    f_lnr = f_vol * 4d0 * const%pi * r**3
    
  end subroutine vol_to_lnr
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Generates the bin grid given the range and number of bins.
  subroutine bin_grid_make(bin_grid, n_bin, v_min, v_max)
    
    !> New bin grid, will be allocated.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Minimum volume (m^3).
    real*8, intent(in) :: v_min
    !> Minimum volume (m^3).
    real*8, intent(in) :: v_max

    call bin_grid_alloc(bin_grid, n_bin)
    call logspace(v_min, v_max, n_bin, bin_grid%v)
    ! dlnr = ln(r(i) / r(i-1))
    bin_grid%dlnr = log(vol2rad(v_max) / vol2rad(v_min)) / dble(n_bin - 1)

  end subroutine bin_grid_make

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Given a bin_grid (which stores the center points of the bins),
  !> find the given edge volume (m^3).
  !!
  !! With n_bin bin centers there are (n_bin + 1) bin edges, so bin
  !! center bin_grid%v(i) is between bin edges i and (i + 1). This
  !! code currently assumes a logarithmically spaced bin grid and
  !! returns logarithmically spaced edges.
  real*8 function bin_edge(bin_grid, i)
    
    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Edge number (1 <= i <= n_bin + 1).
    integer, intent(in) :: i

    real*8 :: log_v_min, log_v_max, log_delta

    call assert(440393735, bin_grid%n_bin > 1)
    log_v_min = log(bin_grid%v(1))
    log_v_max = log(bin_grid%v(bin_grid%n_bin))
    log_delta = (log_v_max - log_v_min) / dble(bin_grid%n_bin - 1)
    bin_edge = exp(log_v_min + (dble(i) - 1.5d0) * log_delta)
    
  end function bin_edge
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the bin number that contains a given particle. This assumes
  !> logarithmically spaced bins.
  integer function bin_grid_particle_in_bin(bin_grid, v)

    !> Bin_grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume of particle.
    real*8, intent(in) :: v

    real*8 :: log_v_min, log_v_max, log_edge_min, log_edge_max
    real*8 :: half_log_delta
    integer :: k

    call assert(448215689, bin_grid%n_bin > 2)
    log_v_min = log(bin_grid%v(1))
    log_v_max = log(bin_grid%v(bin_grid%n_bin))
    half_log_delta = (log_v_max - log_v_min) / dble(2 * (bin_grid%n_bin - 1))
    log_edge_min = log_v_min + half_log_delta
    log_edge_max = log_v_max - half_log_delta
!DEBUG
    write(*,*) 'bin_grid_particle_in_bin() **********************************'
    write(*,*) 'v = ', v
    write(*,*) 'log_v_min = ', log_v_min
    write(*,*) 'log_v_max = ', log_v_max
    write(*,*) 'half_log_delta = ', half_log_delta
    write(*,*) 'log_edge_min = ', log_edge_min
    write(*,*) 'log_edge_max = ', log_edge_max
    write(*,*) '1 = ', (log(v) - log_edge_min)
    write(*,*) '2 = ', (log_edge_max - log_edge_min)
    write(*,*) '3 = ', dble(bin_grid%n_bin - 2)
    write(*,*) '4 = ', (log(v) - log_edge_min) / (log_edge_max - log_edge_min)
    write(*,*) '5 = ', (log(v) - log_edge_min) / (log_edge_max - log_edge_min) &
         * dble(bin_grid%n_bin - 2)
    write(*,*) '6 = ', ceiling((log(v) - log_edge_min) &
         / (log_edge_max - log_edge_min) &
         * dble(bin_grid%n_bin - 2))
!DEBUG
    k = ceiling((log(v) - log_edge_min) / (log_edge_max - log_edge_min) &
         * dble(bin_grid%n_bin - 2)) + 1
!DEBUG
    write(*,*) 'k = ', k
!DEBUG
    k = max(k, 1)
    k = min(k, bin_grid%n_bin)
    bin_grid_particle_in_bin = k
    
  end function bin_grid_particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine inout_write_bin_grid(file, bin_grid)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Bin_grid to write.
    type(bin_grid_t), intent(in) :: bin_grid

    call inout_write_comment(file, "begin bin_grid")
    call inout_write_integer(file, "n_bin", bin_grid%n_bin)
    call inout_write_real_array(file, "center_volumes(m^3)", bin_grid%v)
    call inout_write_real(file, "dlnr", bin_grid%dlnr)
    call inout_write_comment(file, "end bin_grid")
    
  end subroutine inout_write_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine inout_read_bin_grid(file, bin_grid)
    
    !> File to read from.
    type(inout_file_t), intent(inout) :: file
    !> Bin_grid to read.
    type(bin_grid_t), intent(out) :: bin_grid

    call inout_check_comment(file, "begin bin_grid")
    call inout_read_integer(file, "n_bin", bin_grid%n_bin)
    call inout_read_real_array(file, "center_volumes(m^3)", bin_grid%v)
    call inout_read_real(file, "dlnr", bin_grid%dlnr)
    call inout_check_comment(file, "end bin_grid")
    
  end subroutine inout_read_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a bin_grid from a inout file and
  !> generate it.
  subroutine spec_read_bin_grid(file, bin_grid)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid

    integer :: n_bin
    real*8 :: r_min, r_max

    call inout_read_integer(file, 'n_bin', n_bin)
    call inout_read_real(file, 'r_min', r_min)
    call inout_read_real(file, 'r_max', r_max)
    call bin_grid_make(bin_grid, n_bin, rad2vol(r_min), rad2vol(r_max))

  end subroutine spec_read_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_bin_grid(val)

    !> Value to pack.
    type(bin_grid_t), intent(in) :: val

    pmc_mpi_pack_size_bin_grid = &
         pmc_mpi_pack_size_integer(val%n_bin) &
         + pmc_mpi_pack_size_real_array(val%v) &
         + pmc_mpi_pack_size_real(val%dlnr)

  end function pmc_mpi_pack_size_bin_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_bin_grid
