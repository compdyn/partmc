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

  write(*,'(a,e20.10)') 'time (s) = ', time
  call process_env(env)
  call process_info(bin_grid, aero_data, aero_state)
  call process_moments(basename, bin_grid, aero_data, aero_state)
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
  
  subroutine open_output(basename, suffix, out_unit)
    
    ! Allocate a new unit and open it with a filename given by
    ! basename + suffix.

    use pmc_util

    character(len=*), intent(in) :: basename ! basename of the output file
    character(len=*), intent(in) :: suffix ! suffix of the output file
    integer, intent(out) :: out_unit    ! unit for the file

    character(len=len(basename)+len(suffix)) :: filename
    
    filename = basename
    filename((len_trim(filename)+1):) = suffix
    out_unit = get_unit()
    open(out_unit, file=filename)
    write(*,'(a,a)') 'Writing ', trim(filename)
    
  end subroutine open_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_env(env)

    type(env_t), intent(in) :: env      ! environment state

    write(*,'(a,e20.10)') 'temp (K) = ', env%temp
    write(*,'(a,e20.10)') 'rel_humid (1) = ', env%rel_humid
    write(*,'(a,e20.10)') 'pressure (Pa) = ', env%pressure
    write(*,'(a,e20.10)') 'air_den (kg/m^3) = ', env%air_den
    write(*,'(a,e20.10)') 'height (m) = ', env%height

  end subroutine process_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_info(bin_grid, aero_data, aero_state)

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    integer :: n_part

    n_part = total_particles(aero_state)
    write(*,'(a,i20)') 'total particles = ', n_part
    write(*,'(a,e20.10)') 'comp_vol (m^3) = ', aero_state%comp_vol
    write(*,'(a,e20.10)') 'num_dens (#/m^3) = ', &
         dble(n_part) / aero_state%comp_vol

  end subroutine process_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_moments(basename, bin_grid, aero_data, aero_state)

    use pmc_util
    use pmc_aero_binned

    character(len=*), intent(in) :: basename ! basename of the input filename
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    type(aero_binned_t) :: aero_binned
    integer :: f_out, i_bin, i_spec

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)
    call open_output(basename, "_aero_binned.d", f_out)
    write(f_out, '(a1)', advance='no') '#'
    write(f_out, '(a19)', advance='no') 'radius(m)'
    write(f_out, '(a20)', advance='no') 'num_dens(#/m^3)'
    do i_spec = 1,aero_data%n_spec
       write(f_out, '(i4,a1,a15)', advance='no') (i_spec + 2), '/', &
            aero_data%name(i_spec)
    end do
    write(f_out, *) ''
    do i_bin = 1,bin_grid%n_bin
       write(f_out, '(e20.10,e20.10)', advance='no') &
            vol2rad(bin_grid%v(i_bin)), &
            aero_binned%num_den(i_bin)
       do i_spec = 1,aero_data%n_spec
          write(f_out, '(e20.10)', advance='no') &
               aero_binned%vol_den(i_bin, i_spec)
       end do
       write(f_out, *) ''
    end do
    close(unit=f_out)
    call aero_binned_free(aero_binned)

  end subroutine process_moments

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
    call open_output(basename, "_n_orig_part.d", f_out)
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
