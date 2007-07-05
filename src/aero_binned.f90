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

  subroutine alloc_aero_binned(n_bin, n_spec, aero_binned)

    ! Allocates an aero_binned.

    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    type(aero_binned_t), intent(out) :: aero_binned ! bin distribution

    allocate(aero_binned%num_den(n_bin))
    allocate(aero_binned%vol_den(n_bin, n_spec))

  end subroutine alloc_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine free_aero_binned(aero_binned)

    ! Frees all memory.

    type(aero_binned_t), intent(inout) :: aero_binned ! aero_binned to free

    deallocate(aero_binned%num_den)
    deallocate(aero_binned%vol_den)

  end subroutine free_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_aero_binned(file, aero_binned)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_binned_t), intent(in) :: aero_binned ! aero_binned to write

    call inout_write_real_array(file, "num_density(#/m^3)", &
         aero_binned%num_den)
    call inout_write_real_array_2d(file, "vol_density(#/m^3)", &
         aero_binned%vol_den)
    
  end subroutine inout_write_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_binned(file, aero_binned)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_binned_t), intent(out) :: aero_binned ! aero_binned to read

    call inout_read_real_array(file, "num_density(#/m^3)", &
         aero_binned%num_den)
    call inout_read_real_array_2d(file, "vol_density(#/m^3)", &
         aero_binned%vol_den)
    
  end subroutine inout_read_aero_binned

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_aero_binned
