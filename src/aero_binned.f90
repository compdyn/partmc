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

    call inout_write_real_array(file, "num_dens(num/m^3)", &
         aero_binned%num_den)
    call inout_write_real_array_2d(file, "vol_dens(num/m^3)", &
         aero_binned%vol_den)
    
  end subroutine inout_write_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_binned(file, aero_binned)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_binned_t), intent(out) :: aero_binned ! aero_binned to read

    call inout_read_real_array(file, "num_dens(num/m^3)", &
         aero_binned%num_den)
    call inout_read_real_array_2d(file, "vol_dens(num/m^3)", &
         aero_binned%vol_den)
    
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
  
end module pmc_aero_binned
