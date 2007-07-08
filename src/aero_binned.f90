! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Functions that deal with the binned aerosol distributions.

module mod_aero_binned

  type aero_binned_t
     real*8, pointer :: num_den(:)    ! len n_bin, number density (#/m^3)
     real*8, pointer :: vol_den(:,:)  ! n_bin x n_spec, volume density (m^3/m^3)
  end type aero_binned_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_alloc(n_bin, n_spec, aero_binned)

    ! Allocates an aero_binned.

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_binned_t), intent(out) :: aero_binned ! bin distribution

    allocate(aero_binned%num_den(n_bin))
    allocate(aero_binned%vol_den(n_bin, n_spec))

  end subroutine aero_binned_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_binned_free(aero_binned)

    ! Frees all memory.

    type(aero_binned_t), intent(inout) :: aero_binned ! aero_binned to free

    deallocate(aero_binned%num_den)
    deallocate(aero_binned%vol_den)

  end subroutine aero_binned_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_binned(file, aero_binned)
    
    ! Write full state.
    
    use mod_inout
    
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
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_binned_t), intent(out) :: aero_binned ! aero_binned to read

    call inout_read_real_array(file, "num_dens(num/m^3)", &
         aero_binned%num_den)
    call inout_read_real_array_2d(file, "vol_dens(num/m^3)", &
         aero_binned%vol_den)
    
  end subroutine inout_read_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_aero_binned(aero_binned_vec, aero_binned_avg)
    
    ! Computes the average of an array of aero_binned.

    use mod_util

    type(aero_binned_t), intent(in) :: aero_binned_vec(:) ! array of aero_binned
    type(aero_binned_t), intent(out) :: aero_binned_avg   ! average of vec

    integer :: n_bin, n_spec, i_bin, i_spec, n, i

    n_bin = size(aero_binned_vec(1)%vol_den, 1)
    n_spec = size(aero_binned_vec(1)%vol_den, 2)
    call aero_binned_alloc(n_bin, n_spec, aero_binned_avg)
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
    
  end subroutine average_aero_binned
  
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
    call aero_binned_alloc(n_bin, n_spec, aero_binned_to)
    aero_binned_to%num_den = aero_binned_from%num_den
    aero_binned_to%vol_den = aero_binned_from%vol_den

  end subroutine aero_binned_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_dist_to_binned(bin_grid, aero_dist, aero_binned)

    ! Converts an aero_dist to an aero_binned.

    use mod_bin_grid
    use mod_aero_dist

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_dist_t), intent(in) :: aero_dist ! source aero_dist
    type(aero_binned_t), intent(out) :: aero_binned ! must be alloced

    integer :: i_mode, i_bin

    aero_binned%num_den = 0d0
    aero_binned%vol_den = 0d0
    do i_mode = 1,aero_dist%n_modes
       do i_bin = 1,bin_grid%n_bin
          aero_binned%num_den(i_bin) = aero_binned%num_den(i_bin) &
               + aero_dist%modes(i_mode)%n_den(i_bin)
          aero_binned%vol_den(i_bin,:) = aero_binned%vol_den(i_bin,:) &
               + bin_grid%v(i_bin) * aero_dist%modes(i_mode)%n_den(i_bin) &
               * aero_dist%modes(i_mode)%vol_frac
       end do
    end do

  end subroutine aero_dist_to_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_aero_binned
