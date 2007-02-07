! -*- mode: f90; -*-
! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Top level driver that reads .spec file and calls simulation routines.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program partbox

  use mod_read_spec

  integer, parameter :: in_unit = 32
  type(spec_file) :: spec
  character(len=300) :: in_name
  character(len=100) :: run_type
  integer :: i
  
  ! check there is exactly one commandline argument
  if (iargc() .ne. 1) then
     write(6,*) 'Usage: run_mc <filename.d>'
     call exit(2)
  end if
  
  ! get and check first commandline argument (must be "filename.spec")
  call getarg(1, in_name)
  i = len_trim(in_name)
  if (in_name((i-4):i) /= '.spec') then
     write(6,*) 'ERROR: input filename must end in .spec'
     call exit(2)
  end if

  call open_spec(spec, in_name, in_unit)

  call read_string(spec, 'run_type', run_type)

  if (trim(run_type) == 'mc') then
     call partbox_mc(spec)
  elseif (trim(run_type) == 'exact') then
     call partbox_exact(spec)
  elseif (trim(run_type) == 'sect') then
     call partbox_sect(spec)
  else
     write(*,*) 'ERROR: unknown run_type: ', trim(run_type)
     call exit(1)
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partbox_mc(spec)

    use mod_bin
    use mod_array
    use mod_init_dist
    use mod_condensation
    use mod_kernel_sedi
    use mod_kernel_golovin
    use mod_material
    use mod_environ
    use mod_run_mc
    use mod_read_spec

    type(spec_file), intent(out) :: spec     ! spec file

    integer, parameter :: MAX_DIST_ARGS = 10
    
    integer :: MM, M, M_new, i_loop, i, i_bin
    character(len=300) :: out_file_name
    integer, allocatable :: MH(:), bin_n(:)
    type(bin_p), allocatable ::  VH(:)
    real*8, allocatable :: bin_v(:), n_den(:), bin_g(:), bin_gs(:,:)
    real*8 :: dlnr, t_wall_start
    
    character(len=300) :: output_name ! name of output files
    integer :: n_loop             ! number of Monte Carlo loops
    real*8 :: N_0                 ! particle concentration (#/m^3)
    character(len=100) :: kernel_name ! coagulation kernel name
    
    real*8 :: t_max               ! total simulation time (s)
    real*8 :: del_t               ! timestep (s)
    real*8 :: t_output            ! output interval (0 disables) (s)
    real*8 :: t_state             ! state output interval (0 disables) (s)
    real*8 :: t_progress          ! progress printing interval (0 disables) (s)
    
    type(material) :: mat         ! material data
    type(environ) :: env          ! environment data
    
    integer :: n_init_dist        ! number of initial distributions
    integer, allocatable :: dist_n_part(:) ! distribution particle numbers
    character(len=100), allocatable :: dist_types(:) ! distribution names
    real*8, allocatable :: dist_args(:,:) ! distribution arguments
    real*8, allocatable :: dist_vol_frac(:,:) ! distribution composition
    
    integer :: n_bin              ! number of bins
    real*8 :: v_min               ! volume of smallest bin (m^3)
    integer :: scal               ! scale factor (integer)
    
    integer :: rand_init          ! random initialization (0 to auto-init)
    logical :: do_coagulation     ! do coagulation? (yes/no)
    logical :: do_condensation    ! do condensation? (yes/no)
    logical :: do_restart         ! restart from stored state? (yes/no)
    character(len=300) :: restart_name ! filename to restart from
    
    call read_string(spec, 'output_name', output_name)
    call read_integer(spec, 'n_loop', n_loop)
    call read_real(spec, 'N_0', N_0)
    call read_string(spec, 'kernel', kernel_name)
    
    call read_real(spec, 't_max', t_max)
    call read_real(spec, 'del_t', del_t)
    call read_real(spec, 't_output', t_output)
    call read_real(spec, 't_state', t_state)
    call read_real(spec, 't_progress', t_progress)
    
    call read_material(spec, mat)
    call read_environ(spec, env)
    
    call read_integer(spec, 'n_init_dist', n_init_dist)
    allocate(dist_n_part(n_init_dist))
    allocate(dist_vol_frac(n_init_dist, mat%n_spec))
    allocate(dist_types(n_init_dist))
    allocate(dist_args(n_init_dist, MAX_DIST_ARGS))
    
    do i = 1,n_init_dist
       call read_integer(spec, 'n_p', dist_n_part(i))
       call read_real_array(spec, mat%n_spec, 'vol_frac', dist_vol_frac(i,:))
       call read_init_dist(spec, dist_types(i), dist_args(i,:))
    end do
    
    call read_integer(spec, 'n_bin', n_bin)
    call read_real(spec, 'v_min', v_min)
    call read_integer(spec, 'scal', scal)
    
    call read_integer(spec, 'rand_init', rand_init)
    call read_logical(spec, 'do_coagulation', do_coagulation)
    call read_logical(spec, 'do_condensation', do_condensation)
    call read_logical(spec, 'do_restart', do_restart)
    call read_string(spec, 'restart_name', restart_name)
    
    call close_spec(spec)

    ! finished reading .spec data, now do the run
    
    MM = sum(dist_n_part)
    allocate(MH(n_bin), VH(n_bin), bin_v(n_bin), n_den(n_bin))
    allocate(bin_g(n_bin), bin_gs(n_bin,mat%n_spec), bin_n(n_bin))
    call init_array(mat%n_spec, MH, VH)
    call make_bin_grid(n_bin, scal, v_min, bin_v, dlnr)
    
    write(out_file_name, '(a,a,a)') 'out_', trim(output_name), '.d'
    open(30,file=out_file_name)
    call print_header(n_loop, n_bin, mat%n_spec, nint(t_max / t_output) + 1)
    
    if (rand_init /= 0) then
       call srand(rand_init)
    else
       call srand(time())
    end if

    call cpu_time(t_wall_start)
    do i_loop = 1,n_loop
       
       call zero_array(mat%n_spec, MH, VH)
       
       do i = 1,n_init_dist
          call init_dist(dist_types(i), dist_args(i,:), n_bin, bin_v, n_den)
          call dist_to_n(dist_n_part(i), dlnr, n_bin, bin_v, n_den, bin_n)
          call add_particles(n_bin, mat%n_spec, dist_vol_frac(i,:), &
               bin_v, bin_n, MH, VH)
       end do

       M = sum(MH)
       env%V_comp = dble(M) / N_0
       
       ! equlibriate all particles if condensation is active
       if (do_condensation) then
          do i_bin = 1,n_bin
             do i = 1,MH(i_bin)
                call equilibriate_particle(mat%n_spec, VH(i_bin)%p(i,:), &
                     env, mat)
             end do
          end do
       end if
       
       if (trim(kernel_name) == 'sedi') then
          call run_mc(MM, M, mat%n_spec, n_bin, MH, VH, &
               bin_v, bin_g, bin_gs, bin_n, dlnr, &
               kernel_sedi, t_max, t_output, t_state, t_progress, del_t, &
               do_coagulation, do_condensation, do_restart, restart_name, &
               i_loop, n_loop, t_wall_start, env, mat)
       elseif (trim(kernel_name) == 'golovin') then
          call run_mc(MM, M, mat%n_spec, n_bin, MH, VH, &
               bin_v, bin_g, bin_gs, bin_n, dlnr, &
               kernel_golovin, t_max, t_output, t_state, t_progress, del_t, &
               do_coagulation, do_condensation, do_restart, restart_name, &
               i_loop, n_loop, t_wall_start, env, mat)
       else
          write(*,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
          call exit(1)
       end if
       
    end do

  end subroutine partbox_mc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partbox_exact(spec)

    use mod_bin
    use mod_array
    use mod_init_dist
    use mod_condensation
    use mod_kernel_sedi
    use mod_kernel_golovin
    use mod_material
    use mod_environ
    use mod_run_exact
    use mod_read_spec

    type(spec_file), intent(out) :: spec     ! spec file

    integer, allocatable :: bin_n(:)
    real*8, allocatable :: bin_v(:), bin_g(:), bin_gs(:,:)
    real*8 :: dlnr
    character(len=300) :: out_file_name
    
    character(len=300) :: output_name ! name of output files
    real*8 :: N_0                 ! particle concentration (#/m^3)

    character(len=100) :: soln_name ! exact solution name
    real*8 :: mean_vol              ! mean volume of initial distribution
    
    real*8 :: t_max               ! total simulation time (s)
    real*8 :: t_output            ! output interval (0 disables) (s)
    
    type(material) :: mat         ! material data
    type(environ) :: env          ! environment data

    integer :: n_bin              ! number of bins
    real*8 :: v_min               ! volume of smallest bin (m^3)
    integer :: scal               ! scale factor (integer)
    
    call read_string(spec, 'output_name', output_name)
    call read_real(spec, 'N_0', N_0)

    call read_string(spec, 'soln', soln_name)

    if (trim(soln_name) == 'golovin_exp') then
       call read_real(spec, 'mean_vol', mean_vol)
    else
       write(*,*) 'ERROR: unknown solution type: ', trim(soln_name)
       call exit(1)
    end if
    
    call read_real(spec, 't_max', t_max)
    call read_real(spec, 't_output', t_output)
    
    call read_material(spec, mat)
    call read_environ(spec, env)

    call read_integer(spec, 'n_bin', n_bin)
    call read_real(spec, 'v_min', v_min)
    call read_integer(spec, 'scal', scal)

    call close_spec(spec)

    ! finished reading .spec data, now do the run
    
    write(out_file_name, '(a,a,a)') 'out_', trim(output_name), '.d'
    open(30,file=out_file_name)
    call print_header(1, n_bin, mat%n_spec, nint(t_max / t_output) + 1)

    allocate(bin_n(n_bin), bin_v(n_bin), bin_g(n_bin))
    allocate(bin_gs(n_bin,mat%n_spec))
    
    call make_bin_grid(n_bin, scal, v_min, bin_v, dlnr)
    
    if (trim(soln_name) == 'golovin_exp') then
       call run_exact(n_bin, mat%n_spec, bin_v, bin_g, bin_gs, &
            bin_n, dlnr, N_0, mean_vol, mat%rho(1), soln_golovin_exp, t_max, &
            t_output, env, mat)
    else
       write(*,*) 'ERROR: unknown solution type: ', trim(soln_name)
       call exit(1)
    end if
    
  end subroutine partbox_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine partbox_sect(spec)

    use mod_material
    use mod_environ
    use mod_run_sect
    use mod_kernel_sedi
    use mod_kernel_golovin
    use mod_bin
    use mod_array
    use mod_init_dist

    type(spec_file), intent(out) :: spec     ! spec file

    integer, parameter :: MAX_DIST_ARGS = 10

    character(len=300) :: out_file_name
    integer, allocatable :: bin_n(:)
    real*8, allocatable :: bin_v(:), n_den(:), bin_g(:), bin_gs(:,:)
    real*8 :: dlnr
    
    character(len=300) :: output_name ! name of output files
    real*8 :: N_0                 ! particle concentration (#/m^3)
    character(len=100) :: kernel_name ! coagulation kernel name
    
    real*8 :: t_max               ! total simulation time (s)
    real*8 :: del_t               ! timestep (s)
    real*8 :: t_output            ! output interval (0 disables) (s)
    real*8 :: t_progress          ! progress printing interval (0 disables) (s)
    
    type(material) :: mat         ! material data
    type(environ) :: env          ! environment data
    
    character(len=100) :: dist_type ! initial distribution
    real*8 :: dist_args(MAX_DIST_ARGS) ! distribution arguments
    
    integer :: n_bin              ! number of bins
    real*8 :: v_min               ! volume of smallest bin (m^3)
    integer :: scal               ! scale factor (integer)
    
    call read_string(spec, 'output_name', output_name)
    call read_real(spec, 'N_0', N_0)

    call read_string(spec, 'kernel', kernel_name)

    call read_init_dist(spec, dist_type, dist_args)

    call read_real(spec, 't_max', t_max)
    call read_real(spec, 'del_t', del_t)
    call read_real(spec, 't_output', t_output)
    call read_real(spec, 't_progress', t_progress)

    call read_material(spec, mat)
    call read_environ(spec, env)

    call read_integer(spec, 'n_bin', n_bin)
    call read_real(spec, 'v_min', v_min)
    call read_integer(spec, 'scal', scal)

    call close_spec(spec)

    ! finished reading .spec data, now do the run

    allocate(bin_v(n_bin), n_den(n_bin))
    allocate(bin_g(n_bin), bin_gs(n_bin,mat%n_spec), bin_n(n_bin))
    
    write(out_file_name, '(a,a,a)') 'out_', trim(output_name), '.d'
    open(30,file=out_file_name)
    call print_header(1, n_bin, mat%n_spec, nint(t_max / t_output) + 1)
    
    call make_bin_grid(n_bin, scal, v_min, bin_v, dlnr)
    
    call init_dist(dist_type, dist_args, n_bin, bin_v, n_den)

    if (trim(kernel_name) == 'sedi') then
       call run_sect(n_bin, bin_v, dlnr, n_den, N_0, kernel_sedi, &
            t_max, del_t, t_output, t_progress, mat, env)
    elseif (trim(kernel_name) == 'golovin') then
       call run_sect(n_bin, bin_v, dlnr, n_den, N_0, kernel_golovin, &
            t_max, del_t, t_output, t_progress, mat, env)
    else
       write(*,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
       call exit(1)
    end if
    
  end subroutine partbox_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program partbox

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
