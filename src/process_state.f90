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

  if (iargc() < 1) then
     call print_usage()
     call exit(1)
  end if
  
  call get_filename(filename, basename)

  call inout_read_state(filename, bin_grid, aero_data, aero_state, &
       gas_data, gas_state, env, time)

  write(*,'(a,e20.10)') 'time (s) = ', time
  if (iargc() == 1) then
     call process_env(env)
     call process_info(bin_grid, aero_data, aero_state)
     call process_moments(basename, bin_grid, aero_data, aero_state)
     call process_n_orig_part(basename, bin_grid, aero_data, aero_state)
  else
     call process_comp(basename, bin_grid, aero_data, aero_state)
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_usage()

    write(0,*) 'Usage: process_state <filename.d>'
    write(0,*) '       process_state <filename.d> <comp_suffix>' &
         // ' <n_steps> -a <A species> -b <B species>'
    
  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_filename(filename, basename)
    
    character(len=*), intent(out) :: filename ! input filename
    character(len=*), intent(out) :: basename ! basename of the input filename

    integer :: i
    
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
    call aero_binned_write_summary(aero_binned, aero_data, bin_grid, f_out)
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

  subroutine process_comp(basename, bin_grid, aero_data, aero_state)

    ! Compute a histogram of the composition of the particles.

    use pmc_util

    character(len=*), intent(in) :: basename ! basename of the input filename
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    character(len=100) :: comp_suffix, total_comp_suffix, tmp_str
    character(len=100) :: gnuplot_comp_suffix
    integer, allocatable :: comp(:,:), total_comp(:)
    real*8, allocatable :: comp_dens(:,:), total_comp_dens(:)
    integer :: n_steps, n_a, n_b, i, n_zero_vol, i_bin, i_step, i_part
    integer :: f_out
    integer, allocatable :: a_species(:), b_species(:)
    logical :: error
    real*8, allocatable :: props(:), left_props(:)
    type(aero_particle_t), pointer :: particle
    real*8 :: a_vol, b_vol, prop, bin_width

    ! process commandline
    call getarg(2, tmp_str)
    comp_suffix(1:1) = "_"
    comp_suffix(2:) = tmp_str
    comp_suffix((len_trim(comp_suffix)+1):) = '.d'
    total_comp_suffix(1:1) = "_"
    total_comp_suffix(2:) = tmp_str
    total_comp_suffix((len_trim(total_comp_suffix)+1):) = '_total.d'
    gnuplot_comp_suffix(1:1) = "_"
    gnuplot_comp_suffix(2:) = tmp_str
    gnuplot_comp_suffix((len_trim(gnuplot_comp_suffix)+1):) = '_gnuplot.d'
    call getarg(3, tmp_str)
    n_steps = string_to_integer(tmp_str)
    call getarg(4, tmp_str)
    if (tmp_str /= "-a") then
       write(0,*) 'ERROR: argument 4 must be "-a"'
       call exit(1)
    end if
    ! figure out how many A and B species we have
    n_a = 0
    n_b = 0
    do i = 5,iargc()
       call getarg(i, tmp_str)
       if (tmp_str == "-b") then
          n_a = i - 5
          n_b = iargc() - i
          exit
       end if
    end do
    if (n_b == 0) then
       write(0,*) 'ERROR: no B species specified'
       call exit(1)
    end if
    if (n_a == 0) then
       write(0,*) 'ERROR: no A species specified'
       call exit(1)
    end if
    ! read the A and B species
    allocate(a_species(n_a))
    allocate(b_species(n_b))
    error = .false.
    do i = 1,n_a
       call getarg(4 + i, tmp_str)
       a_species(i) = aero_data_spec_by_name(aero_data, tmp_str)
       if (a_species(i) == 0) then
          write(0,'(a,a)') 'ERROR: unknown species: ', trim(tmp_str)
          error = .true.
       end if
    end do
    do i = 1,n_b
       call getarg(5 + n_a + i, tmp_str)
       b_species(i) = aero_data_spec_by_name(aero_data, tmp_str)
       if (b_species(i) == 0) then
          write(0,'(a,a)') 'ERROR: unknown species: ', trim(tmp_str)
          error = .true.
       end if
    end do
    if (error) then
       call exit(1)
    end if
    
    ! compute compositions
    allocate(props(n_steps),left_props(n_steps))
    allocate(comp(bin_grid%n_bin, n_steps), total_comp(n_steps))
    allocate(comp_dens(bin_grid%n_bin, n_steps), total_comp_dens(n_steps))
    do i_step = 1,n_steps
       props(i_step) = (dble(i_step) - 0.5d0) / dble(n_steps)
       left_props(i_step) = dble(i_step - 1) / dble(n_steps)
    end do
    comp = 0
    total_comp = 0
    n_zero_vol = 0
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          particle => aero_state%bins(i_bin)%particle(i_part)
          a_vol = 0d0
          do i = 1,n_a
             a_vol = a_vol + particle%vol(a_species(i))
          end do
          b_vol = 0d0
          do i = 1,n_b
             b_vol = b_vol + particle%vol(b_species(i))
          end do
          call assert(a_vol >= 0d0)
          call assert(b_vol >= 0d0)
          if ((a_vol == 0d0) .and. (b_vol == 0d0)) then
             n_zero_vol = n_zero_vol + 1
          else
             prop = b_vol / (a_vol + b_vol)
             i_step = floor(dble(n_steps) * prop)
             if (i_step == n_steps) then
                i_step = n_steps - 1
             end if
             comp(i_bin, i_step) = comp(i_bin, i_step) + 1
             total_comp(i_step) = total_comp(i_step) + 1
          end if
       end do
    end do
    if (n_zero_vol > 0) then
       write(0,*) 'WARNING: number of particles without any A or B volume: ', &
            n_zero_vol
    end if
    bin_width = 1d0 / dble(n_steps)
    do i_step = 1,n_steps
       total_comp_dens(i_step) = dble(total_comp(i_step)) &
            / aero_state%comp_vol / bin_width
       do i_bin = 1,bin_grid%n_bin
          comp_dens(i_bin, i_step) = dble(comp(i_bin, i_step)) &
               / aero_state%comp_vol / bin_grid%dlnr / bin_width
       end do
    end do

    ! write comp output
    call open_output(basename, comp_suffix, f_out)
    write(f_out, '(a)') '# rows are bins'
    write(f_out, '(a)') '# columns are proportions'
    write(f_out, '(a)') '# entries are number densities'
    write(f_out, '(a1,a24)', advance='no') '#', 'radius(m)'
    do i_step = 1,n_steps
       write(f_out, '(f24.1,a1)', advance='no') (props(i_step) * 100d0), '%'
    end do
    write(f_out, *) ''
    do i_bin = 1,bin_grid%n_bin
       write(f_out, '(e25.15)', advance='no') vol2rad(bin_grid%v(i_bin))
       do i_step = 1,n_steps
          write(f_out, '(e25.15)', advance='no') comp_dens(i_bin, i_step)
       end do
       write(f_out, *) ''
    end do
    close(unit=f_out)

    ! write comp output in gnuplot pm3d format
    call open_output(basename, gnuplot_comp_suffix, f_out)
    write(f_out, '(a)') '# number densities are scaled for bin widths in' &
         // ' both bins and proportions'
    write(f_out, '(a1,a24,a25,a25)') '#', 'radius(m)', 'proportion(0-1)', &
         'num_dens(#/m^3)'
    do i_bin = 1,bin_grid%n_bin
       do i_step = 1,n_steps
          write(f_out, '(e25.15,e25.15,e25.15)') vol2rad(bin_grid%v(i_bin)), &
               left_props(i_step), comp_dens(i_bin, i_step)
       end do
       ! extra dummy value at 1 to terminate the plot
       write(f_out, '(e25.15,e25.15,e25.15)') vol2rad(bin_grid%v(i_bin)), &
            1d0, 0d0
       write(f_out, *) ''
    end do
    ! FIXME: the resulting plot will discard the largest bin, which
    ! could be fixed by using bin edges rather than centers
    close(unit=f_out)
    
    ! write total_comp output
    call open_output(basename, total_comp_suffix, f_out)
    write(f_out, '(a1,a24,a25)') '#', 'prop', 'num_dens'
    do i_step = 1,n_steps
       write(f_out, '(e25.15,e25.15)') props(i_step), total_comp_dens(i_step)
    end do
    close(unit=f_out)
    
    deallocate(comp, total_comp, comp_dens, total_comp_dens)
    deallocate(a_species, b_species, props, left_props)
    
  end subroutine process_comp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program process_state
