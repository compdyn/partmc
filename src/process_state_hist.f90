! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process data with a histogram over bins and a second scalar quantity.

module pmc_process_state_hist
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine process_hist(basename, type, bin_grid, env, aero_data, &
       aero_state, step_comp_grid, step_comp, particle_func, do_sum, &
       time, index)

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
    logical, intent(in) :: do_sum       ! if quantity can be summed
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index

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

       real*8 function particle_func(aero_particle, aero_data, env)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env
         type(aero_particle_t), intent(in) :: aero_particle ! particle
         type(aero_data_t), intent(in) :: aero_data ! aero_data structure
         type(env_t), intent(in) :: env ! environment state
       end function particle_func
    end interface

    real*8, allocatable :: num_den(:,:), num_den_tot(:), num_den_spec(:,:)
    real*8, allocatable :: vol_den(:,:), vol_den_tot(:), vol_den_spec(:,:)
    real*8, allocatable :: mass_den(:,:), mass_den_tot(:), mass_den_spec(:,:)
    real*8, allocatable :: mole_den(:,:), mole_den_tot(:), mole_den_spec(:,:)
    real*8, allocatable :: sum_den(:)
    real*8, allocatable :: num_cum_tot(:), num_cum_spec(:,:)
    real*8, allocatable :: vol_cum_tot(:), vol_cum_spec(:,:)
    real*8, allocatable :: mass_cum_tot(:), mass_cum_spec(:,:)
    real*8, allocatable :: mole_cum_tot(:), mole_cum_spec(:,:)
    integer :: n_step, i_step, i_bin, i_part, n_invalid
    type(aero_particle_t), pointer :: aero_particle
    real*8 :: scale, scale_bin, scale_step, num, vol, mass, mole, val
    real*8, pointer :: step_grid(:) ! length n_step + 1

    call step_comp_grid(bin_grid, env, aero_data, aero_state, &
         n_step, step_grid)
    allocate(num_den(bin_grid%n_bin, n_step), num_den_tot(n_step))
    allocate(num_den_spec(n_step, aero_data%n_spec))
    allocate(vol_den(bin_grid%n_bin, n_step), vol_den_tot(n_step))
    allocate(vol_den_spec(n_step, aero_data%n_spec))
    allocate(mass_den(bin_grid%n_bin, n_step), mass_den_tot(n_step))
    allocate(mass_den_spec(n_step, aero_data%n_spec))
    allocate(mole_den(bin_grid%n_bin, n_step), mole_den_tot(n_step))
    allocate(mole_den_spec(n_step, aero_data%n_spec))
    allocate(sum_den(bin_grid%n_bin))
    allocate(num_cum_tot(n_step), num_cum_spec(n_step,aero_data%n_spec))
    allocate(vol_cum_tot(n_step), vol_cum_spec(n_step,aero_data%n_spec))
    allocate(mass_cum_tot(n_step), mass_cum_spec(n_step,aero_data%n_spec))
    allocate(mole_cum_tot(n_step), mole_cum_spec(n_step,aero_data%n_spec))

    num_den = 0d0
    num_den_tot = 0d0
    num_den_spec = 0d00
    vol_den = 0d0
    vol_den_tot = 0d0
    vol_den_spec = 0d0
    mass_den = 0d0
    mass_den_tot = 0d0
    mass_den_spec = 0d0
    mole_den = 0d0
    mole_den_tot = 0d0
    mole_den_spec = 0d0
    sum_den = 0d0

    scale_bin = 1d0 / bin_grid%dlnr / aero_state%comp_vol
    scale_step = 1d0 / dble(n_step) / aero_state%comp_vol
    scale = scale_bin * scale_step

    n_invalid = 0
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          aero_particle => aero_state%bins(i_bin)%particle(i_part)
          i_step = step_comp(bin_grid, env, aero_data, n_step, aero_particle)
          if (i_step == 0) then
             n_invalid = n_invalid + 1
             continue
          end if
          num = 1d0
          vol = aero_particle_volume(aero_particle)
          mass = aero_particle_mass(aero_particle, aero_data)
          mole = aero_particle_moles(aero_particle, aero_data)
          num_den(i_bin, i_step) = num_den(i_bin, i_step) + num * scale
          num_den_tot(i_step) = num_den_tot(i_step) + num * scale_step
          vol_den(i_bin, i_step) = vol_den(i_bin, i_step) + vol * scale
          vol_den_tot(i_step) = vol_den_tot(i_step) + vol * scale_step
          vol_den_spec(i_step, :) = vol_den_spec(i_step, :) &
               + aero_particle%vol * scale_step
          mass_den(i_bin, i_step) = mass_den(i_bin, i_step) + mass * scale
          mass_den_tot(i_step) = mass_den_tot(i_step) + mass * scale_step
          mass_den_spec(i_step, :) = mass_den_spec(i_step, :) &
               + aero_particle%vol * aero_data%density * scale_step
          mole_den(i_bin, i_step) = mole_den(i_bin, i_step) + mole * scale
          mole_den_tot(i_step) = mole_den_tot(i_step) + mole * scale_step
          mole_den_spec(i_step, :) = mole_den_spec(i_step, :) &
               + aero_particle%vol * aero_data%density &
               / aero_data%molec_weight * scale_step
          if (do_sum) then
             val = particle_func(aero_particle, aero_data, env)
             sum_den(i_bin) = sum_den(i_bin) + val * scale_bin
          end if
       end do
    end do

    num_cum_tot(1) = num_den_tot(1)
    num_cum_spec(1,:) = num_den_spec(1,:)
    vol_cum_tot(1) = vol_den_tot(1)
    vol_cum_spec(1,:) = vol_den_spec(1,:)
    mass_cum_tot(1) = mass_den_tot(1)
    mass_cum_spec(1,:) = mass_den_spec(1,:)
    mole_cum_tot(1) = mole_den_tot(1)
    mole_cum_spec(1,:) = mole_den_spec(1,:)
    do i_step = 2,n_step
        num_cum_tot(i_step) = num_cum_tot(i_step-1) + num_den_tot(i_step)
        num_cum_spec(i_step,:) = num_cum_spec(i_step-1,:) &
             + num_den_spec(i_step,:)
        vol_cum_tot(i_step) = vol_cum_tot(i_step-1) + vol_den_tot(i_step)
        vol_cum_spec(i_step,:) = vol_cum_spec(i_step-1,:) & 
             + vol_den_spec(i_step,:)
        mass_cum_tot(i_step) = mass_cum_tot(i_step-1) + mass_den_tot(i_step)
        mass_cum_spec(i_step,:) = mass_cum_spec(i_step-1,:) & 
             + mass_den_spec(i_step,:)
        mole_cum_tot(i_step) = mole_cum_tot(i_step-1) + mole_den_tot(i_step)
        mole_cum_spec(i_step,:) = mole_cum_spec(i_step-1,:) & 
             + mole_den_spec(i_step,:)
    enddo

    if (n_invalid > 0) then
       write(0,*) 'WARNING: number of particles without bin: ', n_invalid
    end if

    call write_hist_matrix(basename, time, index, type, "_num", &
         "number density", n_step, bin_grid, aero_data, step_grid, &
         num_den, num_den_tot, num_den_spec, sum_den, do_sum, &
         num_cum_tot, num_cum_spec)
    call write_hist_matrix(basename, time, index, type, "_vol", &
         "volume density", n_step, bin_grid, aero_data, step_grid, &
         vol_den, vol_den_tot, vol_den_spec, sum_den, do_sum, &
         vol_cum_tot, vol_cum_spec)
    call write_hist_matrix(basename, time, index, type, "_mass", &
         "mass density", n_step, bin_grid, aero_data, step_grid, &
         mass_den, mass_den_tot, mass_den_spec, sum_den, do_sum, &
         mass_cum_tot, mass_cum_spec)
    call write_hist_matrix(basename, time, index, type, "_mole", &
         "molar density", n_step, bin_grid, aero_data, step_grid, &
         mole_den, mole_den_tot, mole_den_spec, sum_den, do_sum, &
         mole_cum_tot, mole_cum_spec)

    deallocate(step_grid)
    deallocate(num_den, num_den_tot)
    deallocate(num_den_spec)
    deallocate(vol_den, vol_den_tot)
    deallocate(vol_den_spec)
    deallocate(mass_den, mass_den_tot)
    deallocate(mass_den_spec)
    deallocate(mole_den, mole_den_tot)
    deallocate(mole_den_spec)
    deallocate(sum_den)
    deallocate(num_cum_tot, num_cum_spec)
    deallocate(vol_cum_tot, vol_cum_spec)
    deallocate(mass_cum_tot, mass_cum_spec)
    deallocate(mole_cum_tot, mole_cum_spec)

  end subroutine process_hist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_hist_matrix(basename, time, index, type, unit, &
       unit_descript, n_step, bin_grid, aero_data, step_grid, den, &
       den_tot, den_spec, sum_den, do_sum, cum_tot, cum_spec)

    ! Helper function for process_hist() to write out the histograms.

    use pmc_util
    use pmc_bin_grid
    use pmc_aero_data

    character(len=*), intent(in) :: basename ! basename of the filename
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    character(len=*), intent(in) :: type ! type of the filename
    character(len=*), intent(in) :: unit ! unit of the data
    character(len=*), intent(in) :: unit_descript ! description of the unit
    integer, intent(in) :: n_step       ! number of histogram steps
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    real*8, intent(in) :: step_grid(n_step + 1) ! step grid edges
    real*8, intent(in) :: den(bin_grid%n_bin,n_step) ! density per bin per step
    real*8, intent(in) :: den_tot(n_step) ! density per step
    real*8, intent(in) :: den_spec(n_step,aero_data%n_spec) ! species density
                                                            ! per step
    real*8, intent(in) :: sum_den(bin_grid%n_bin) ! summed density
    logical, intent(in) :: do_sum       ! whether to write sum_den
    real*8, intent(in) :: cum_tot(n_step) ! cumulative density per step
    real*8, intent(in) :: cum_spec(n_step,aero_data%n_spec) ! cumul. 
                                                            ! species density
                                                            ! per step

    character(len=200) :: outname
    integer :: f_out, i_bin, i_step, i_spec
    real*8 :: d

    outname = basename
    outname((len_trim(outname)+1):) = type
    outname((len_trim(outname)+1):) = unit

    ! write full output as a matrix
    call open_output(outname, "_matrix.d", f_out)
    write(f_out, '(a)') '# histogram matrix'
    write(f_out, '(a,e25.15)') '# time is ', time
    write(f_out, '(a,i10)') '# index is ', index
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
       write(f_out, '(e25.15)', advance='no') &
            vol2rad(bin_edge(bin_grid, i_bin))
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
    write(f_out, '(a,e25.15)') '# time is ', time
    write(f_out, '(a,i10)') '# index is ', index
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
    write(f_out, '(a,e25.15)') '# time is ', time
    write(f_out, '(a,i10)') '# index is ', index
    write(f_out, '(a,a)') '# quantities are ', unit_descript
    write(f_out, '(a,a,a)') '# steps are ', type
    write(f_out, '(a)') '# last row quantity is junk'
    write(f_out, '(a1,a24,a25)') '#', 'step', 'quantity'
    do i_step = 1,n_step
       write(f_out, '(e25.15,e25.15)') step_grid(i_step), den_tot(i_step)
    end do
    write(f_out, '(e25.15,e25.15)') step_grid(n_step + 1), 0d0
    close(unit=f_out)

    ! write histogram species output
    call open_output(outname, "_spec.d", f_out)
    write(f_out, '(a)') '# histogram species'
    write(f_out, '(a,e25.15)') '# time is ', time
    write(f_out, '(a,i10)') '# index is ', index
    write(f_out, '(a,a)') '# quantities are ', unit_descript
    write(f_out, '(a,a,a)') '# steps are ', type
    write(f_out, '(a)') '# last row quantity is junk'
    write(f_out, '(a1,a25)', advance='no') '#',' step'
    do i_spec = 1,aero_data%n_spec
       write(f_out, '(i9,a1,a15)', advance='no') &
            i_spec+1, '/', aero_data%name(i_spec)
    end do
    write(f_out, *) ''
    do i_step = 1,n_step
       write(f_out, '(e25.15)', advance='no') step_grid(i_step)
       do i_spec = 1,aero_data%n_spec
          write(f_out, '(e25.15)', advance='no') &
               den_spec(i_step,i_spec)
       end do
       write(f_out, *) ''
    end do
    write(f_out, '(e25.15)', advance='no') step_grid(n_step + 1)
    do i_spec = 1,aero_data%n_spec
          write(f_out, '(e25.15)', advance='no') 0d0
    end do
    write(f_out, *) ''
    close(unit=f_out)

    ! write summed output
    if (do_sum) then
       call open_output(outname, "_sum.d", f_out)
       write(f_out, '(a)') '# histogram sum'
       write(f_out, '(a,e25.15)') '# time is ', time
       write(f_out, '(a,i10)') '# index is ', index
       write(f_out, '(a,a)') '# quantities are ', type
       write(f_out, '(a)') '# last row quantity is junk'
       write(f_out, '(a1,a24,a25)') '#', 'bin', 'quantity'
       do i_bin = 1,bin_grid%n_bin
          write(f_out, '(e25.15,e25.15)') bin_grid%v(i_bin), sum_den(i_bin)
       end do
       close(unit=f_out)
    end if

    ! write cumulative total output
    call open_output(outname, "_cum_total.d", f_out)
    write(f_out, '(a)') '# cumulative totals'
    write(f_out, '(a,e25.15)') '# time is ', time
    write(f_out, '(a,i10)') '# index is ', index
    write(f_out, '(a,a)') '# quantities are ', unit_descript
    write(f_out, '(a,a,a)') '# steps are ', type
    write(f_out, '(a)') '# last row quantity is junk'
    write(f_out, '(a1,a24,a25)') '#', 'step', 'quantity'
    do i_step = 1,n_step
       write(f_out, '(e25.15,e25.15)') step_grid(i_step), cum_tot(i_step)
    end do
    write(f_out, '(e25.15,e25.15)') step_grid(n_step + 1), 0d0
    close(unit=f_out)

    ! write cumulative species output
    call open_output(outname, "_cum_spec.d", f_out)
    write(f_out, '(a)') '# cumulative species'
    write(f_out, '(a,e25.15)') '# time is ', time
    write(f_out, '(a,i10)') '# index is ', index
    write(f_out, '(a,a)') '# quantities are ', unit_descript
    write(f_out, '(a,a,a)') '# steps are ', type
    write(f_out, '(a)') '# last row quantity is junk'
    write(f_out, '(a1,a25)', advance='no') '#','step'
    do i_spec = 1,aero_data%n_spec
       write(f_out, '(i9,a1,a15)', advance='no') &
            i_spec+1, '/', aero_data%name(i_spec)
    end do
    write(f_out, *) ''
    do i_step = 1,n_step
       write(f_out, '(e25.15)', advance='no') step_grid(i_step)
       do i_spec = 1,aero_data%n_spec
          write(f_out, '(e25.15)', advance='no') &
               cum_spec(i_step,i_spec)
       end do
       write(f_out, *) ''
    end do
    write(f_out, '(e25.15)', advance='no') step_grid(n_step + 1)
    do i_spec = 1,aero_data%n_spec
          write(f_out, '(e25.15)', advance='no') 0d0
    end do
    write(f_out, *) ''
    close(unit=f_out)

  end subroutine write_hist_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine per_particle_step_comp_grid(bin_grid, env, aero_data, &
       aero_state, n_step, step_grid, particle_func, min_val, max_val, &
       logscale)

    ! Generic histogram helper for per-particle scalar quantity
    ! determined by particle_func.

    use pmc_util
    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_state
    use pmc_env

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure
    integer, intent(out) :: n_step      ! number of histogram steps
    real*8, pointer :: step_grid(:)     ! step grid edges
    real*8, intent(out) :: min_val      ! minimum value of per-particle quantity
    real*8, intent(out) :: max_val      ! maximum value of per-particle quantity
    logical, intent(in) :: logscale     ! whether to use a log scale
    
    interface
       real*8 function particle_func(aero_particle, aero_data, env)
         use pmc_aero_particle
         use pmc_aero_data
         use pmc_env
         type(aero_particle_t), intent(in) :: aero_particle ! particle
         type(aero_data_t), intent(in) :: aero_data ! aero_data structure
         type(env_t), intent(in) :: env ! environment state
       end function particle_func
    end interface

    character(len=100) :: tmp_str
    integer :: i_step, i_bin, i_part
    real*8 :: rh, val
    logical :: first_time

    ! process commandline
    call getarg(3, tmp_str)
    n_step = string_to_integer(tmp_str)

    ! find max and min
    first_time = .true.
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          val = particle_func(aero_state%bins(i_bin)%particle(i_part), &
               aero_data, env)
          if (first_time) then
             min_val = val
             max_val = val
          else
             if (val < min_val) min_val = val
             if (val > max_val) max_val = val
          end if
          first_time = .false.
       end do
    end do

    ! make step grid
    allocate(step_grid(n_step + 1))
    if (max_val <= min_val) then
       write(0,*) 'ERROR: min_val is not greater than max_val'
       call exit(1)
    end if
    if (logscale) then
       call logspace(min_val, max_val, n_step + 1, step_grid)
    else
       call linspace(min_val, max_val, n_step + 1, step_grid)
    end if

  end subroutine per_particle_step_comp_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function orig_part_particle_func(aero_particle, aero_data, env)

    use pmc_aero_particle
    use pmc_aero_data
    use pmc_env

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(env_t), intent(in) :: env      ! environment state

    orig_part_particle_func = dble(aero_particle%n_orig_part)

  end function orig_part_particle_func

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

  real*8 function comp_particle_func(aero_particle, aero_data, env)

    use pmc_aero_particle
    use pmc_aero_data
    use pmc_env
    use pmc_util

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(env_t), intent(in) :: env      ! environment state

    integer :: n_a, n_b, i
    integer, pointer :: a_species(:), b_species(:)
    real*8 :: a_vol, b_vol

    common/comp_step_comp_c/ n_a, n_b, a_species, b_species

    a_vol = 0d0
    do i = 1,n_a
       a_vol = a_vol + aero_particle%vol(a_species(i))
    end do
    b_vol = 0d0
    do i = 1,n_b
       b_vol = b_vol + aero_particle%vol(b_species(i))
    end do
    call assert(880038232, a_vol >= 0d0)
    call assert(715496111, b_vol >= 0d0)
    if ((a_vol == 0d0) .and. (b_vol == 0d0)) then
       comp_particle_func = 0d0
    else
       comp_particle_func = b_vol / (a_vol + b_vol)
    end if

  end function comp_particle_func

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
    call linspace(0d0, 1d0, n_step + 1, step_grid)

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

    comp_step_comp = linspace_find(0d0, 1d0, n_step + 1, &
         comp_particle_func(aero_particle, aero_data, env))
    
  end function comp_step_comp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function kappa_particle_func(aero_particle, aero_data, env)

    use pmc_aero_particle
    use pmc_aero_data
    use pmc_env

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(env_t), intent(in) :: env      ! environment state

    real*8 :: rh, supersat

    rh = aero_particle_kappa_rh(aero_particle, aero_data, env)
    supersat = rh - 1d0
    kappa_particle_func = supersat

  end function kappa_particle_func

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

    logical :: logscale
    real*8 :: min_val, max_val

    common/kappa_step_comp_c/ min_val, max_val

    logscale = .true.
    call per_particle_step_comp_grid(bin_grid, env, aero_data, aero_state, &
         n_step, step_grid, kappa_particle_func, min_val, max_val, logscale)

  end subroutine kappa_step_comp_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function kappa_step_comp(bin_grid, env, aero_data, n_step, &
       aero_particle)

    ! Histogram helper for kappa.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_particle
    use pmc_env
    use pmc_util

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(in) :: n_step  ! number of histogram steps
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    real*8 :: min_val, max_val

    common/kappa_step_comp_c/ min_val, max_val

    kappa_step_comp = logspace_find(min_val, max_val, n_step + 1, &
         kappa_particle_func(aero_particle, aero_data, env))

  end function kappa_step_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function absorb_particle_func(aero_particle, aero_data, env)

    use pmc_aero_particle
    use pmc_aero_data
    use pmc_env

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(env_t), intent(in) :: env      ! environment state

    real*8 :: rh, supersat

    absorb_particle_func = aero_particle%absorb_cross_sect

  end function absorb_particle_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine absorb_step_comp_grid(bin_grid, env, aero_data, aero_state, &
       n_step, step_grid)

    ! Histogram helper for absorb.

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

    logical :: logscale
    real*8 :: min_val, max_val

    common/absorb_step_comp_c/ min_val, max_val

    logscale = .false.
    call per_particle_step_comp_grid(bin_grid, env, aero_data, aero_state, &
         n_step, step_grid, absorb_particle_func, min_val, max_val, logscale)

  end subroutine absorb_step_comp_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function absorb_step_comp(bin_grid, env, aero_data, n_step, &
       aero_particle)

    ! Histogram helper for absorb.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_particle
    use pmc_env
    use pmc_util

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(in) :: n_step  ! number of histogram steps
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    real*8 :: min_val, max_val

    common/absorb_step_comp_c/ min_val, max_val

    absorb_step_comp = linspace_find(min_val, max_val, n_step + 1, &
         absorb_particle_func(aero_particle, aero_data, env))

  end function absorb_step_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function scatter_particle_func(aero_particle, aero_data, env)

    use pmc_aero_particle
    use pmc_aero_data
    use pmc_env

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(env_t), intent(in) :: env      ! environment state

    real*8 :: rh, supersat

    scatter_particle_func = aero_particle%scatter_cross_sect

  end function scatter_particle_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine scatter_step_comp_grid(bin_grid, env, aero_data, aero_state, &
       n_step, step_grid)

    ! Histogram helper for scatter.

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

    logical :: logscale
    real*8 :: min_val, max_val

    common/scatter_step_comp_c/ min_val, max_val

    logscale = .false.
    call per_particle_step_comp_grid(bin_grid, env, aero_data, aero_state, &
         n_step, step_grid, scatter_particle_func, min_val, max_val, logscale)

  end subroutine scatter_step_comp_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function scatter_step_comp(bin_grid, env, aero_data, n_step, &
       aero_particle)

    ! Histogram helper for scatter.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_particle
    use pmc_env
    use pmc_util

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(in) :: n_step  ! number of histogram steps
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    real*8 :: min_val, max_val

    common/scatter_step_comp_c/ min_val, max_val

    scatter_step_comp = linspace_find(min_val, max_val, n_step + 1, &
         scatter_particle_func(aero_particle, aero_data, env))

  end function scatter_step_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function extinct_particle_func(aero_particle, aero_data, env)

    use pmc_aero_particle
    use pmc_aero_data
    use pmc_env

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(env_t), intent(in) :: env      ! environment state

    real*8 :: rh, supersat

    extinct_particle_func = aero_particle%absorb_cross_sect &
         + aero_particle%scatter_cross_sect

  end function extinct_particle_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine extinct_step_comp_grid(bin_grid, env, aero_data, aero_state, &
       n_step, step_grid)

    ! Histogram helper for extinct.

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

    logical :: logscale
    real*8 :: min_val, max_val

    common/extinct_step_comp_c/ min_val, max_val

    logscale = .false.
    call per_particle_step_comp_grid(bin_grid, env, aero_data, aero_state, &
         n_step, step_grid, extinct_particle_func, min_val, max_val, logscale)

  end subroutine extinct_step_comp_grid
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function extinct_step_comp(bin_grid, env, aero_data, n_step, &
       aero_particle)

    ! Histogram helper for extinct.

    use pmc_bin_grid
    use pmc_aero_data
    use pmc_aero_particle
    use pmc_env
    use pmc_util

    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(env_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(in) :: n_step  ! number of histogram steps
    type(aero_particle_t), intent(in) :: aero_particle ! particle

    real*8 :: min_val, max_val

    common/extinct_step_comp_c/ min_val, max_val

    extinct_step_comp = linspace_find(min_val, max_val, n_step + 1, &
         extinct_particle_func(aero_particle, aero_data, env))

  end function extinct_step_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_process_state_hist
