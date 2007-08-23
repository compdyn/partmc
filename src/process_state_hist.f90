! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process the saved state files to obtain summary data.

module pmc_process_state_hist
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine process_hist(basename, type, bin_grid, env, aero_data, &
       aero_state, step_comp_grid, step_comp)

    ! Compute a histogram by calling the step_comp() function on each
    ! particle and write out various representations of it.
    
    use pmc_util
    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env
    
    character(len=*), intent(in) :: basename ! basename of the input filename
    character(len=*), intent(in) :: type ! type of the input filename
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    interface
       subroutine step_comp_grid(bin_grid, env, aero_data, aero_state, &
            n_step, step_grid)
         use pmc_bin_grid
         use pmc_aero_data
         use pmc_aero_state
         use pmc_env
         type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
         type(env_t), intent(in) :: env ! environment state
         type(aero_data_t), intent(in) :: aero_data ! aero_data structure
         type(aero_state_t), intent(in) :: aero_state ! aero_state structure
         integer, intent(out) :: n_step  ! number of histogram steps
         real*8, pointer :: step_grid(:) ! len n_step+1, step grid edges
       end subroutine step_comp_grid

       integer function step_comp(bin_grid, env, aero_data, n_step, &
            aero_particle)
         use pmc_bin_grid
         use pmc_aero_data
         use pmc_aero_particle
         use pmc_env
         type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
         type(env_t), intent(in) :: env ! environment state
         type(aero_data_t), intent(in) :: aero_data ! aero_data structure
         integer, intent(in) :: n_step  ! number of histogram steps
         type(aero_particle_t), intent(in) :: aero_particle ! particle
       end function step_comp
    end interface

    real*8, allocatable :: num_den(:,:), num_den_tot(:)
    real*8, allocatable :: vol_den(:,:), vol_den_tot(:)
    real*8, allocatable :: mass_den(:,:), mass_den_tot(:)
    real*8, allocatable :: mole_den(:,:), mole_den_tot(:)
    integer :: n_step, i_step, i_bin, i_part, n_invalid
    type(aero_particle_t), pointer :: particle
    real*8 :: scale, scale_bin, num, vol, mass, mole
    real*8, pointer :: step_grid(:) ! length n_step + 1

    call step_comp_grid(bin_grid, env, aero_data, aero_state, n_step, step_grid)
    allocate(num_den(bin_grid%n_bin, n_step), num_den_tot(n_step))
    allocate(vol_den(bin_grid%n_bin, n_step), vol_den_tot(n_step))
    allocate(mass_den(bin_grid%n_bin, n_step), mass_den_tot(n_step))
    allocate(mole_den(bin_grid%n_bin, n_step), mole_den_tot(n_step))

    num_den = 0d0
    num_den_tot = 0d0
    vol_den = 0d0
    vol_den_tot = 0d0
    mass_den = 0d0
    mass_den_tot = 0d0
    mole_den = 0d0
    mole_den_tot = 0d0

    scale = 1d0 / bin_grid%dlnr / dble(n_step)
    scale_bin = 1d0 / bin_grid%dlnr
    n_invalid = 0
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          particle => aero_state%bins(i_bin)%particle(i_part)
          i_step = step_comp(bin_grid, env, aero_data, n_step, particle)
          if (i_step == 0) then
             n_invalid = n_invalid + 1
             continue
          end if
          num = 1d0
          vol = aero_particle_volume(particle)
          mass = aero_particle_mass(particle, aero_data)
          mole = aero_particle_moles(particle, aero_data)
          num_den(i_bin, i_step) = num_den(i_bin, i_step) + num * scale
          num_den_tot(i_step) = num_den_tot(i_step) + num * scale_bin
          vol_den(i_bin, i_step) = vol_den(i_bin, i_step) + vol * scale
          vol_den_tot(i_step) = vol_den_tot(i_step) + vol * scale_bin
          mass_den(i_bin, i_step) = mass_den(i_bin, i_step) + mass * scale
          mass_den_tot(i_step) = mass_den_tot(i_step) + mass * scale_bin
          mole_den(i_bin, i_step) = mole_den(i_bin, i_step) + mole * scale
          mole_den_tot(i_step) = mole_den_tot(i_step) + mole * scale_bin
       end do
    end do
    if (n_invalid > 0) then
       write(0,*) 'WARNING: number of particles without bin: ', n_invalid
    end if

    call write_hist_matrix(basename, type, "_num", "number density", &
         n_step, bin_grid, step_grid, num_den, num_den_tot)
    call write_hist_matrix(basename, type, "_vol", "volume density", &
         n_step, bin_grid, step_grid, vol_den, vol_den_tot)
    call write_hist_matrix(basename, type, "_mass", "mass density", &
         n_step, bin_grid, step_grid, mass_den, mass_den_tot)
    call write_hist_matrix(basename, type, "_mole", "molar density", &
         n_step, bin_grid, step_grid, mole_den, mole_den_tot)

    deallocate(step_grid)
    deallocate(num_den, num_den_tot)
    deallocate(vol_den, vol_den_tot)
    deallocate(mass_den, mass_den_tot)
    deallocate(mole_den, mole_den_tot)

  end subroutine process_hist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_hist_matrix(basename, type, unit, unit_descript, &
       n_step, bin_grid, step_grid, den, den_tot)

    ! Helper function for process_hist() to write out the histograms.

    use pmc_util
    use pmc_bin_grid

    character(len=*), intent(in) :: basename ! basename of the filename
    character(len=*), intent(in) :: type ! type of the filename
    character(len=*), intent(in) :: unit ! unit of the data
    character(len=*), intent(in) :: unit_descript ! description of the unit
    integer, intent(in) :: n_step       ! number of histogram steps
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    real*8, intent(in) :: step_grid(n_step + 1) ! step grid edges
    real*8, intent(in) :: den(bin_grid%n_bin,n_step) ! density per bin per step
    real*8, intent(in) :: den_tot(n_step) ! density per step

    character(len=200) :: outname
    integer :: f_out, i_bin, i_step
    real*8 :: d

    outname = basename
    outname((len_trim(outname)+1):) = type
    outname((len_trim(outname)+1):) = unit

    ! write full output as a matrix
    call open_output(outname, "_matrix.d", f_out)
    write(f_out, '(a)') '# histogram matrix'
    write(f_out, '(a)') '# rows are size bins'
    write(f_out, '(a,a,a)') '# columns are ', type, ' bins'
    write(f_out, '(a,a)') '# entries are ', unit_descript
    write(f_out, '(a)') '# first row (from 2nd column) is step edges'
    write(f_out, '(a)') '# last row (from 2nd column) is junk'
    write(f_out, '(a)') '# first column (from 2nd row) is radius (m)'
    write(f_out, '(a)') '# last column (from 2nd row) is junk'
    write(f_out, '(e25.15)', advance='no') 0d0
    do i_step = 1,(n_step + 1)
       write(f_out, '(e25.15)', advance='no') step_grid(i_step)
    end do
    write(f_out, *) ''
    do i_bin = 1,(bin_grid%n_bin + 1)
       write(f_out, '(e25.15)', advance='no') vol2rad(bin_edge(bin_grid, i_bin))
       do i_step = 1,(n_step + 1)
          if ((i_bin <= bin_grid%n_bin) .and. (i_step <= n_step)) then
             d = den(i_bin, i_step)
          else
             d = 0d0
          end if
          write(f_out, '(e25.15)', advance='no') d
       end do
       write(f_out, *) ''
    end do
    close(unit=f_out)

    ! write full output in gnuplot pm3d format
    call open_output(outname, "_gnuplot.d", f_out)
    write(f_out, '(a,a)') '# quantity is ', unit_descript
    write(f_out, '(a1,a24,a25,a25)') '#', 'radius(m)', 'step', 'quantity'
    do i_bin = 1,(bin_grid%n_bin + 1)
       do i_step = 1,(n_step + 1)
          if ((i_bin <= bin_grid%n_bin) .and. (i_step <= n_step)) then
             d = den(i_bin, i_step)
          else
             d = 0d0
          end if
          write(f_out, '(e25.15,e25.15,e25.15)') &
               vol2rad(bin_edge(bin_grid, i_bin)), step_grid(i_step), d
       end do
       write(f_out, *) ''
    end do
    close(unit=f_out)
    
    ! write totaled output
    call open_output(outname, "_total.d", f_out)
    write(f_out, '(a)') '# histogram totals'
    write(f_out, '(a,a)') '# quantities are ', unit_descript
    write(f_out, '(a,a,a)') '# steps are ', type
    write(f_out, '(a)') '# last row quantity is junk'
    write(f_out, '(a1,a24,a25)') '#', 'step', 'quantity'
    do i_step = 1,n_step
       write(f_out, '(e25.15,e25.15)') step_grid(i_step), den_tot(i_step)
    end do
    write(f_out, '(e25.15,e25.15)') step_grid(n_step + 1), 0d0
    close(unit=f_out)

  end subroutine write_hist_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine orig_part_step_comp_grid(bin_grid, env, aero_data, aero_state, &
       n_step, step_grid)

    ! Histogram helper for n_orig_part.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure
    integer, intent(out) :: n_step  ! number of histogram steps
    real*8, pointer :: step_grid(:) ! step grid edges

    integer :: i_bin, i_part, n, n_orig_part_max

    ! determine the max of n_orig_part
    n_orig_part_max = 0
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          n = aero_state%bins(i_bin)%particle(i_part)%n_orig_part
          if (n > n_orig_part_max) then
             n_orig_part_max = n
          end if
       end do
    end do
    
    n_step = n_orig_part_max
    allocate(step_grid(n_step + 1))
    do n = 1,(n_step + 1)
       step_grid(n) = dble(n - 1)
    end do

  end subroutine orig_part_step_comp_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer function orig_part_step_comp(bin_grid, env, aero_data, n_step, &
       aero_particle)

    ! Histogram helper for n_orig_part.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_particle
    use pmc_env

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(in) :: n_step       ! number of histogram steps
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    orig_part_step_comp = aero_particle%n_orig_part

  end function orig_part_step_comp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine comp_step_comp_grid(bin_grid, env, aero_data, aero_state, &
       n_step, step_grid)

    ! Histogram helper for composition.

    use pmc_util
    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure
    integer, intent(out) :: n_step  ! number of histogram steps
    real*8, pointer :: step_grid(:) ! step grid edges

    integer :: n_a, n_b, i
    integer, pointer :: a_species(:), b_species(:)
    logical :: error
    character(len=100) :: tmp_str

    common/comp_step_comp_c/ n_a, n_b, a_species, b_species

    ! process commandline
    call getarg(3, tmp_str)
    n_step = string_to_integer(tmp_str)
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

    ! allocate the step grid
    allocate(step_grid(n_step + 1))
    do i = 1,(n_step + 1)
       step_grid(i) = dble(i - 1) / dble(n_step)
    end do

  end subroutine comp_step_comp_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer function comp_step_comp(bin_grid, env, aero_data, n_step, &
       aero_particle)

    ! Histogram helper for composition.

    use pmc_util
    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_particle
    use pmc_env

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(in) :: n_step       ! number of histogram steps
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    integer :: n_a, n_b
    integer, pointer :: a_species(:), b_species(:)
    logical :: error
    integer :: i, i_step
    real*8 :: a_vol, b_vol, prop

    common/comp_step_comp_c/ n_a, n_b, a_species, b_species

    a_vol = 0d0
    do i = 1,n_a
       a_vol = a_vol + aero_particle%vol(a_species(i))
    end do
    b_vol = 0d0
    do i = 1,n_b
       b_vol = b_vol + aero_particle%vol(b_species(i))
    end do
    call assert(a_vol >= 0d0)
    call assert(b_vol >= 0d0)
    if ((a_vol == 0d0) .and. (b_vol == 0d0)) then
       comp_step_comp = 0
    else
       prop = b_vol / (a_vol + b_vol)
       i_step = floor(dble(n_step) * prop) + 1
       if (i_step > n_step) then
          i_step = n_step
       end if
       comp_step_comp = i_step
    end if
    
  end function comp_step_comp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine kappa_step_comp_grid(bin_grid, env, aero_data, aero_state, &
       n_step, step_grid)

    ! Histogram helper for kappa.

    use pmc_util
    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure
    integer, intent(out) :: n_step  ! number of histogram steps
    real*8, pointer :: step_grid(:) ! step grid edges

    character(len=100) :: tmp_str
    integer :: i_step, i_bin, i_part
    real*8 :: rh, min_rh, max_rh
    logical :: first_time

    common/kappa_step_comp_c/ min_rh, max_rh

    ! process commandline
    call getarg(3, tmp_str)
    n_step = string_to_integer(tmp_str)

    ! find max and min kappas
    first_time = .true.
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          rh = aero_particle_kappa_rh(aero_state%bins(i_bin)%particle(i_part), &
               aero_data, env)
          if (first_time) then
             min_rh = rh
             max_rh = rh
          else
             if (rh < min_rh) min_rh = rh
             if (rh > max_rh) max_rh = rh
          end if
          first_time = .false.
       end do
    end do

    ! make step grid
    allocate(step_grid(n_step + 1))
    call linspace(min_rh, max_rh, n_step + 1, step_grid)

  end subroutine kappa_step_comp_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function kappa_step_comp(bin_grid, env, aero_data, n_step, &
       aero_particle)

    ! Histogram helper for kappa.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_particle
    use pmc_env

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(in) :: n_step  ! number of histogram steps
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    real*8 :: rh, min_rh, max_rh

    common/kappa_step_comp_c/ min_rh, max_rh

    rh = aero_particle_kappa_rh(aero_particle, aero_data, env)
    kappa_step_comp = floor((rh - min_rh) / (max_rh - min_rh) &
         * dble(n_step)) + 1
    if (kappa_step_comp > n_step) then
       kappa_step_comp = n_step
    end if

  end function kappa_step_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_process_state_hist
