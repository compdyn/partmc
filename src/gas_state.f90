! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Gas state.

module mod_gas_state

  type gas_state_t
     real*8, pointer :: conc(:)          ! length n_spec, concentration (ppb)
  end type gas_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_gas_state(n_spec, gas_state)

    ! Allocate storage for gas species.

    integer, intent(in) :: n_spec       ! number of species
    type(gas_state_t), intent(out) :: gas_state ! gas state to be allocated

    allocate(gas_state%conc(n_spec))

  end subroutine alloc_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine copy_gas_state(from_state, to_state)

    ! Copy to an already allocated to_state.

    type(gas_state_t), intent(in) :: from_state ! existing gas state
    type(gas_state_t), intent(out) :: to_state ! must be allocated already

    integer :: n_spec

    n_spec = size(from_state%conc)
    deallocate(to_state%conc)
    allocate(to_state%conc(n_spec))
    to_state%conc = from_state%conc

  end subroutine copy_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_gas_state(file, gas_state)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(gas_state_t), intent(in) :: gas_state ! gas_state to write

    call inout_write_real_array(file, "conc(ppb)", gas_state%conc)
    
  end subroutine inout_write_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_gas_state(file, gas_state)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(gas_state_t), intent(out) :: gas_state ! gas_state to read

    call inout_read_real_array(file, "conc(ppb)", gas_state%conc)
    
  end subroutine inout_read_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_gas_state(file, gas_data, name, gas_state)

    ! Read gas state from the file named on the line read from file.

    use mod_inout
    use mod_gas_data

    type(inout_file_t), intent(inout) :: file ! inout file
    type(gas_data_t), intent(out) :: gas_data ! gas data
    character(len=*), intent(in) :: name ! name of data line for filename
    type(gas_state_t), intent(out) :: gas_state ! gas data

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)
    integer :: species_data_shape(2)

    ! read the filename then read the data from that file
    call inout_read_string(file, name, read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 0, species_name, species_data)
    call inout_close(read_file)

    ! check the data size
    species_data_shape = shape(species_data)
    n_species = species_data_shape(1)
    if (.not. ((species_data_shape(2) == 1) .or. (n_species == 0))) then
       write(0,*) 'ERROR: each line in ', trim(read_name), &
            ' should contain exactly one data value'
       call exit(1)
    end if

    ! copy over the data
    call alloc_gas_state(gas_data%n_spec, gas_state)
    gas_state%conc = 0d0
    do i = 1,n_species
       species = gas_spec_by_name(gas_data, species_name(i))
       if (species == 0) then
          write(0,*) 'ERROR: unknown species ', trim(species_name(i)), &
               ' in file ', trim(read_name)
          call exit(1)
       end if
       gas_state%conc(species) = species_data(i,1)
    end do

  end subroutine spec_read_gas_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine average_gas_state(gas_state_vec, gas_state_avg)
    
    ! Computes the average of an array of gas_state.

    use mod_util

    type(gas_state_t), intent(in) :: gas_state_vec(:) ! array of gas_state
    type(gas_state_t), intent(out) :: gas_state_avg   ! average of gas_state_vec

    integer :: n_spec, i_spec, i, n

    n_spec = size(gas_state_vec(1)%conc)
    call alloc_gas_state(n_spec, gas_state_avg)
    n = size(gas_state_vec)
    do i_spec = 1,n_spec
       call average_real((/(gas_state_vec(i)%conc(i_spec),i=1,n)/), &
            gas_state_avg%conc(i_spec))
    end do
    
  end subroutine average_gas_state
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_gas_state
