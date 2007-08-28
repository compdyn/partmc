! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Functions that deal with the binned aerosol distributions.

module pmc_aero_binned

  type aero_binned_t
     real*8, pointer :: num_den(:)    ! len n_bin, number density (#/m^3)
     real*8, pointer :: vol_den(:,:)  ! n_bin x n_spec, volume density (m^3/m^3)
  end type aero_binned_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_alloc(aero_binned, n_bin, n_spec)

    ! Allocates an aero_binned.

    type(aero_binned_t), intent(out) :: aero_binned ! bin distribution
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species

    allocate(aero_binned%num_den(n_bin))
    allocate(aero_binned%vol_den(n_bin, n_spec))
    call aero_binned_zero(aero_binned)

  end subroutine aero_binned_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_free(aero_binned)

    ! Frees all memory.

    type(aero_binned_t), intent(inout) :: aero_binned ! aero_binned to free

    deallocate(aero_binned%num_den)
    deallocate(aero_binned%vol_den)

  end subroutine aero_binned_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_zero(aero_binned)

    ! Zeros an aero_binned.

    type(aero_binned_t), intent(inout) :: aero_binned ! bin distribution

    aero_binned%num_den = 0d0
    aero_binned%vol_den = 0d0

  end subroutine aero_binned_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_add_particle_in_bin(aero_binned, bin_grid, &
       bin, comp_vol, aero_particle)

    ! Updates binned data structures for the addition of the given
    ! particle that must be in the given bin.

    use pmc_bin_grid
    use pmc_aero_particle

    type(aero_binned_t), intent(inout) :: aero_binned ! binned distributions
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    integer, intent(in) :: bin          ! bin number
    real*8, intent(in) :: comp_vol      ! computational volume (m^3)
    type(aero_particle_t), intent(in) :: aero_particle ! particle to add

    aero_binned%num_den(bin) = aero_binned%num_den(bin) &
         + 1d0 / comp_vol / bin_grid%dlnr
    aero_binned%vol_den(bin,:) = aero_binned%vol_den(bin,:) &
         + aero_particle%vol / comp_vol / bin_grid%dlnr

  end subroutine aero_binned_add_particle_in_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_add_particle(aero_binned, bin_grid, &
       comp_vol, aero_particle)

    ! Updates binned data structures for the addition of the given
    ! particle.

    use pmc_bin_grid
    use pmc_aero_particle

    type(aero_binned_t), intent(inout) :: aero_binned ! binned distributions
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(in) :: comp_vol      ! computational volume (m^3)
    type(aero_particle_t), intent(in) :: aero_particle ! particle to add

    call aero_binned_add_particle_in_bin(aero_binned, bin_grid, &
         aero_particle_in_bin(aero_particle, bin_grid), comp_vol, aero_particle)

  end subroutine aero_binned_add_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_remove_particle_in_bin(aero_binned, bin_grid, &
       bin, comp_vol, aero_particle)

    ! Updates binned data structures for the removal of the given
    ! particle that must be in the given bin.

    use pmc_bin_grid
    use pmc_aero_particle

    type(aero_binned_t), intent(inout) :: aero_binned ! binned distributions
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    integer, intent(in) :: bin          ! bin number
    real*8, intent(in) :: comp_vol      ! computational volume (m^3)
    type(aero_particle_t), intent(in) :: aero_particle ! particle to remove

    aero_binned%num_den(bin) = aero_binned%num_den(bin) &
         - 1d0 / comp_vol / bin_grid%dlnr
    aero_binned%vol_den(bin,:) = aero_binned%vol_den(bin,:) &
         - aero_particle%vol / comp_vol / bin_grid%dlnr

  end subroutine aero_binned_remove_particle_in_bin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_remove_particle(aero_binned, bin_grid, &
       comp_vol, aero_particle)

    ! Updates binned data structures for the removal of the given
    ! particle.

    use pmc_bin_grid
    use pmc_aero_particle

    type(aero_binned_t), intent(inout) :: aero_binned ! binned distributions
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    real*8, intent(in) :: comp_vol      ! computational volume (m^3)
    type(aero_particle_t), intent(in) :: aero_particle ! particle to remove

    call aero_binned_remove_particle_in_bin(aero_binned, bin_grid, &
         aero_particle_in_bin(aero_particle, bin_grid), comp_vol, aero_particle)

  end subroutine aero_binned_remove_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_binned(file, aero_binned)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_binned_t), intent(in) :: aero_binned ! aero_binned to write

    call inout_write_comment(file, "begin aero_binned")
    call inout_write_real_array(file, "num_dens(num/m^3)", &
         aero_binned%num_den)
    call inout_write_real_array_2d(file, "vol_dens(num/m^3)", &
         aero_binned%vol_den)
    call inout_write_comment(file, "end aero_binned")
    
  end subroutine inout_write_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_binned(file, aero_binned)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_binned_t), intent(out) :: aero_binned ! aero_binned to read

    call inout_check_comment(file, "begin aero_binned")
    call inout_read_real_array(file, "num_dens(num/m^3)", &
         aero_binned%num_den)
    call inout_read_real_array_2d(file, "vol_dens(num/m^3)", &
         aero_binned%vol_den)
    call inout_check_comment(file, "end aero_binned")
    
  end subroutine inout_read_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_average(aero_binned_vec, aero_binned_avg)
    
    ! Computes the average of an array of aero_binned.

    use pmc_util

    type(aero_binned_t), intent(in) :: aero_binned_vec(:) ! array of aero_binned
    type(aero_binned_t), intent(out) :: aero_binned_avg   ! average of vec

    integer :: n_bin, n_spec, i_bin, i_spec, n, i

    n_bin = size(aero_binned_vec(1)%vol_den, 1)
    n_spec = size(aero_binned_vec(1)%vol_den, 2)
    call aero_binned_alloc(aero_binned_avg, n_bin, n_spec)
    n = size(aero_binned_vec)
    do i_bin = 1,n_bin
       call average_real((/(aero_binned_vec(i)%num_den(i_bin),i=1,n)/), &
            aero_binned_avg%num_den(i_bin))
       do i_spec = 1,n_spec
          call average_real((/(aero_binned_vec(i)%vol_den(i_bin,i_spec),&
               i=1,n)/), &
               aero_binned_avg%vol_den(i_bin,i_spec))
       end do
    end do
    
  end subroutine aero_binned_average
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_add(aero_binned, aero_binned_delta)

    ! Adds aero_binned_delta to aero_binned.

    type(aero_binned_t), intent(inout) :: aero_binned ! base aero_binned
    type(aero_binned_t), intent(in) :: aero_binned_delta ! aero_binned to add

    aero_binned%num_den = aero_binned%num_den + aero_binned_delta%num_den
    aero_binned%vol_den = aero_binned%vol_den + aero_binned_delta%vol_den

  end subroutine aero_binned_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_sub(aero_binned, aero_binned_delta)

    ! Subtracts aero_binned_delta from aero_binned.

    type(aero_binned_t), intent(inout) :: aero_binned ! base aero_binned
    type(aero_binned_t), intent(in) :: aero_binned_delta ! aero_binned to sub

    aero_binned%num_den = aero_binned%num_den - aero_binned_delta%num_den
    aero_binned%vol_den = aero_binned%vol_den - aero_binned_delta%vol_den

  end subroutine aero_binned_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_scale(aero_binned, alpha)

    ! Scales by a real number.

    type(aero_binned_t), intent(inout) :: aero_binned ! base aero_binned
    real*8, intent(in) :: alpha         ! scale factor

    aero_binned%num_den = aero_binned%num_den * alpha
    aero_binned%vol_den = aero_binned%vol_den * alpha

  end subroutine aero_binned_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_copy(aero_binned_from, aero_binned_to)

    ! Copies all data.

    type(aero_binned_t), intent(in) :: aero_binned_from ! base aero_binned
    type(aero_binned_t), intent(out) :: aero_binned_to ! already alloced

    integer :: n_bin, n_spec

    call aero_binned_free(aero_binned_to)
    n_bin = size(aero_binned_from%vol_den, 1)
    n_spec = size(aero_binned_from%vol_den, 2)
    call aero_binned_alloc(aero_binned_to, n_bin, n_spec)
    aero_binned_to%num_den = aero_binned_from%num_den
    aero_binned_to%vol_den = aero_binned_from%vol_den

  end subroutine aero_binned_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_add_aero_dist(aero_binned, bin_grid, aero_dist)

    ! Converts an aero_dist to an aero_binned.

    use pmc_bin_grid
    use pmc_aero_dist

    type(aero_binned_t), intent(out) :: aero_binned ! must be alloced
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_dist_t), intent(in) :: aero_dist ! source aero_dist

    integer :: i_mode, i_bin

    do i_mode = 1,aero_dist%n_mode
       do i_bin = 1,bin_grid%n_bin
          aero_binned%num_den(i_bin) = aero_binned%num_den(i_bin) &
               + aero_dist%mode(i_mode)%num_den(i_bin)
          aero_binned%vol_den(i_bin,:) = aero_binned%vol_den(i_bin,:) &
               + bin_grid%v(i_bin) * aero_dist%mode(i_mode)%num_den(i_bin) &
               * aero_dist%mode(i_mode)%vol_frac
       end do
    end do

  end subroutine aero_binned_add_aero_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_aero_binned(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(aero_binned_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_size_aero_binned = &
         pmc_mpi_pack_size_real_array(val%num_den) &
         + pmc_mpi_pack_size_real_array_2d(val%vol_den)

  end function pmc_mpi_pack_size_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_aero_binned(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_binned_t), intent(in) :: val          ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%num_den)
    call pmc_mpi_pack_real_array_2d(buffer, position, val%vol_den)
    call assert(position - prev_position == pmc_mpi_pack_size_aero_binned(val))
#endif

  end subroutine pmc_mpi_pack_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_aero_binned(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_binned_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%num_den)
    call pmc_mpi_unpack_real_array_2d(buffer, position, val%vol_den)
    call assert(position - prev_position == pmc_mpi_pack_size_aero_binned(val))
#endif

  end subroutine pmc_mpi_unpack_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_aero_binned(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_mpi

    type(aero_binned_t), intent(in) :: val ! value to average
    type(aero_binned_t), intent(out) :: val_avg ! result

    call pmc_mpi_reduce_avg_real_array(val%num_den, val_avg%num_den)
    call pmc_mpi_reduce_avg_real_array_2d(val%vol_den, val_avg%vol_den)

  end subroutine pmc_mpi_reduce_avg_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_write_summary(aero_binned, aero_data, &
       bin_grid, time, index, out_unit)

    use pmc_aero_data
    use pmc_util
    use pmc_bin_grid
    
    type(aero_binned_t), intent(in) :: aero_binned ! aero_binned
    type(aero_data_t), intent(in) :: aero_data ! aero_data
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    integer, intent(in) :: out_unit     ! unit number to write to
    
    integer :: col_num, i_bin, i_spec
    
    write(out_unit, '(a,e20.10)') '# time(s) = ', time
    write(out_unit, '(a,i10)') '# index = ', index
    write(out_unit, '(a)') &
         '# VL species are volume density (m^3/m^3)'
    write(out_unit, '(a)') &
         '# MS species are mass density (kg/m^3)'
    write(out_unit, '(a)') &
         '# ML species are mole density (mole/m^3)'
    write(out_unit, '(a1)', advance='no') '#'
    write(out_unit, '(a24)', advance='no') 'radius(m)'
    write(out_unit, '(a25)', advance='no') 'num_dens(#/m^3)'
    col_num = 2
    do i_spec = 1,aero_data%n_spec
       col_num = col_num + 1
       write(out_unit, '(i6,a4,a15)', advance='no') &
            col_num, '-VL/', aero_data%name(i_spec)
    end do
    do i_spec = 1,aero_data%n_spec
       col_num = col_num + 1
       write(out_unit, '(i6,a4,a15)', advance='no') &
            col_num, '-MS/', aero_data%name(i_spec)
    end do
    do i_spec = 1,aero_data%n_spec
       col_num = col_num + 1
       write(out_unit, '(i6,a4,a15)', advance='no') &
            col_num, '-ML/', aero_data%name(i_spec)
    end do
    write(out_unit, *) ''
    do i_bin = 1,bin_grid%n_bin
       write(out_unit, '(e25.15,e25.15)', advance='no') &
            vol2rad(bin_grid%v(i_bin)), &
            aero_binned%num_den(i_bin)
       do i_spec = 1,aero_data%n_spec
          write(out_unit, '(e25.15)', advance='no') &
               aero_binned%vol_den(i_bin, i_spec)
       end do
       do i_spec = 1,aero_data%n_spec
          write(out_unit, '(e25.15)', advance='no') &
               aero_binned%vol_den(i_bin, i_spec) &
               * aero_data%density(i_spec)
       end do
       do i_spec = 1,aero_data%n_spec
          write(out_unit, '(e25.15)', advance='no') &
               aero_binned%vol_den(i_bin, i_spec) &
               * aero_data%density(i_spec) &
               / aero_data%molec_weight(i_spec)
       end do
       write(out_unit, *) ''
    end do

  end subroutine aero_binned_write_summary
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_binned
