! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process the saved state files to obtain summary data.

program process_state

  use pmc_inout
  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env
  use pmc_output_state

  character(len=100) :: filename        ! input filename
  character(len=100) :: basename        ! basename of the input filename
  type(bin_grid_t) :: bin_grid          ! bin_grid structure
  type(aero_data_t) :: aero_data        ! aero_data structure
  type(aero_state_t) :: aero_state      ! aero_state structure
  type(gas_data_t) :: gas_data          ! gas_data structure
  type(gas_state_t) :: gas_state        ! gas_state structure
  type(env_t) :: env                    ! env structure
  real*8 :: time                        ! current time

  call get_filename(filename, basename)

  call inout_read_state(filename, bin_grid, aero_data, aero_state, &
       gas_data, gas_state, env, time)

  call process_n_orig_part(basename, bin_grid, aero_data, aero_state)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_filename(filename, basename)
    
    character(len=*), intent(out) :: filename ! input filename
    character(len=*), intent(out) :: basename ! basename of the input filename

    integer :: i
    
    ! check there is exactly one commandline argument
    if (iargc() .ne. 1) then
       write(0,*) 'Usage: process_state <filename.d>'
       call exit(1)
    end if
    
    ! get and check first commandline argument (must be "filename.d")
    call getarg(1, filename)
    i = len_trim(filename)
    if (i > len(filename)) then
       write(0,*) 'ERROR: filename too long'
       call exit(1)
    end if
    if ((filename(i:i) /= 'd') .or. &
         (filename((i-1):(i-1)) /= '.')) then
       write(0,*) 'ERROR: Filename must end in .d'
       call exit(1)
    end if
    
    ! chop .d off the end of the filename to get the basename
    basename = filename
    basename((i-1):i) = '  '
    
  end subroutine get_filename
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_n_orig_part(basename, bin_grid, aero_data, aero_state)

    ! Compute a histogram of the number of coagulation events per
    ! particle.

    use pmc_util

    character(len=*), intent(in) :: basename ! basename of the input filename
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    integer, allocatable :: n_orig_part(:,:)
    integer :: i_bin, i_part, n, n_orig_part_max, f_out
    character(len=len(basename)+50) :: filename

    ! determine the max number of coag events
    n_orig_part_max = 0
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          n = aero_state%bins(i_bin)%particle(i_part)%n_orig_part
          if (n > n_orig_part_max) then
             n_orig_part_max = n
          end if
       end do
    end do

    ! compute the histogram
    allocate(n_orig_part(bin_grid%n_bin, n_orig_part_max))
    n_orig_part = 0
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          n = aero_state%bins(i_bin)%particle(i_part)%n_orig_part
          call assert(n > 0)
          n_orig_part(i_bin, n) = n_orig_part(i_bin, n) + 1
       end do
    end do

    ! write output
    filename = basename
    filename((len_trim(filename)+1):) = '_n_orig_part.d'
    f_out = get_unit()
    open(f_out, file=filename)
    do i_bin = 1,bin_grid%n_bin
       do n = 1,n_orig_part_max
          write(f_out, '(i20)', advance='no') n_orig_part(i_bin, n)
       end do
       write(f_out, *) ''
    end do
    close(unit=f_out)
    
    deallocate(n_orig_part)
    
  end subroutine process_n_orig_part
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program process_state
