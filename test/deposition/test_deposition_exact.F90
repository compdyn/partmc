program test_deposition_exact

  use pmc_deposition
  use pmc_constants
  use pmc_run_part
  use pmc_spec_file
  use pmc_util
  use pmc_aero_state
  use pmc_aero_data
  use pmc_gas_state
  use pmc_gas_data
  use pmc_scenario
  use pmc_env_state

  implicit none

  !> Output unit number.
  integer, parameter :: out_unit = 34
  !> Output filename.
  character(len=*), parameter :: out_name = "out/deposition_data.txt"

  character(len=300) :: spec_name
  real(kind=dp),allocatable, dimension(:) :: diameter
  real(kind=dp),allocatable, dimension(:) :: num_conc
  real(kind=dp),allocatable, dimension(:) :: particle_mass
  real(kind=dp),allocatable, dimension(:) :: total_mass
  real(kind=dp),allocatable, dimension(:) :: density
  real(kind=dp),allocatable, dimension(:) :: total_initial_mass
  integer :: i_spec
  integer :: i_mode
  integer :: optind
  real(kind=dp) :: time
  real(kind=dp) :: del_t
  real(kind=dp) :: c
  real(kind=dp) :: vd
  real(kind=dp) :: vs
  real(kind=dp) :: rs
  real(kind=dp), allocatable, dimension(:) :: vd_rate
  integer :: i
  integer :: ntime
  type(run_part_opt_t) :: run_part_opt
  type(spec_file_t) :: file
  character(len=100) :: run_type,restart_filename
  logical :: do_restart
  real(kind=dp) :: n_part
  type(scenario_t) :: scenario
  type(aero_data_t) :: aero_data
  type(gas_data_t) :: gas_data
  type(aero_dist_t) :: aero_dist_init
  character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
  type(env_state_t) :: env_state
  type(gas_state_t) :: gas_state_init
  type(spec_file_t) :: sub_file

  real(kind=dp) :: A
  real(kind=dp) :: gamma
  real(kind=dp) :: alpha
  real(kind=dp) :: lambda
  real(kind=dp) :: viscosk
  real(kind=dp) :: viscosd
  real(kind=dp) :: rho_air

  real(kind=dp) :: ustar
  real(kind=dp) :: aer_res_a

  ! Allocates
  call gas_data_allocate(gas_data)
  call gas_state_allocate(gas_state_init)
  call aero_data_allocate(aero_data)
  call aero_dist_allocate(aero_dist_init)
  call env_state_allocate(env_state)
  call scenario_allocate(scenario)

  optind = 1
  call get_command_argument(optind, spec_name)
  print*,spec_name
  call spec_file_open(spec_name, file)

  ! Read in spec file stuff
  call spec_file_read_string(file, 'run_type', run_type)
  call spec_file_read_string(file, 'output_prefix', &
       run_part_opt%output_prefix)
  call spec_file_read_integer(file, 'n_repeat', run_part_opt%n_repeat)
  call spec_file_read_real(file, 'n_part', n_part)
  call spec_file_read_logical(file, 'restart', do_restart)
  if (do_restart) then
     call spec_file_read_string(file, 'restart_file', restart_filename)
  end if
  call spec_file_read_real(file, 't_max', run_part_opt%t_max)
  call spec_file_read_real(file, 'del_t', run_part_opt%del_t)
  call spec_file_read_real(file, 't_output', run_part_opt%t_output)
  call spec_file_read_real(file, 't_progress', run_part_opt%t_progress)

  print*, 'dt = ', run_part_opt%del_t,'t output = ', run_part_opt%t_output

  ! Dummy gas info
  call spec_file_read_string(file, 'gas_data', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_gas_data(sub_file, gas_data)
  call spec_file_close(sub_file)

  call spec_file_read_string(file, 'gas_init', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_gas_state(sub_file, gas_data, &
       gas_state_init)
  call spec_file_close(sub_file)
  ! Aerosol info
  call spec_file_read_string(file, 'aerosol_data', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_data(sub_file, aero_data)
  call spec_file_close(sub_file)
  call spec_file_read_string(file, 'aerosol_init', sub_filename)
  call spec_file_open(sub_filename, sub_file)
  call spec_file_read_aero_dist(sub_file, aero_data, aero_dist_init)
  call spec_file_close(sub_file)
  ! Other info
  call spec_file_read_scenario(file, gas_data, aero_data, scenario)
  call spec_file_read_env_state(file, env_state)
  call scenario_init_env_state(scenario, env_state, &
       env_state%elapsed_time)
  ! Various error checking?

  allocate(density(aero_data%n_spec))
  allocate(diameter(aero_data%n_spec))
  allocate(num_conc(aero_data%n_spec))

  density = 0.0d0
  diameter = 0.0d0
  num_conc = 0.0d0

  do i_mode = 1, size(aero_dist_init%mode)
     ! Check if not monodisperse
     if (aero_dist_init%mode(i_mode)%type .ne. 3) then
        print*, 'Must be monodisperse'
        stop
     end if

     do i_spec = 1, aero_data%n_spec
        if ((aero_dist_init%mode(i_mode)%vol_frac(i_spec) == 1.0d0)) then
           if (density(i_spec) == 0.0d0) then
              diameter(i_spec) = 2.0* aero_dist_init%mode(i_mode)%char_radius
              num_conc(i_spec) = aero_dist_init%mode(i_mode)%num_conc
              density(i_spec) = density(i_spec) &
                 + aero_dist_init%mode(i_mode)%vol_frac(i_spec) &
                 * aero_data%density(i_spec)
           else
              print*, 'Each species must have at most 1 mode'
              stop
           end if
        end if
     end do
  end do

  print*, 'Mean diameter = ', diameter
  print*, 'Mean density = ', density
  print*, 'Number concentration = ', num_conc

  call spec_file_close(file)

  allocate(particle_mass(aero_data%n_spec))
  allocate(total_initial_mass(aero_data%n_spec))

  do i_spec = 1, aero_data%n_spec
     particle_mass(i_spec) = density(i_spec) * (1.0d0/6.0d0) * const%pi &
        * diameter(i_spec)**3.0d0
     total_initial_mass(i_spec) = num_conc(i_spec)*particle_mass(i_spec)
  end do

  print*,'Total mass concentration = ', total_initial_mass

  ! Other things that are still hard coded
  aer_res_a = .0d0
  ustar = 1.0
  gamma = .6d0
  A = 2.0/1000.0
  alpha =  .8d0
  ! Environmental values
  rho_air = env_state_air_den(env_state)
  viscosd = env_state_air_dynamic_viscosity(env_state)
  viscosk = viscosd / rho_air
  lambda = env_state_air_mean_free_path(env_state)
  allocate(vd_rate(aero_data%n_spec))
  do i_spec = 1, aero_data%n_spec
     if (diameter(i_spec) > 0.0d0) then
     vs = calculate_vs(diameter(i_spec), density(i_spec), viscosd, lambda)
     rs = calculate_rs(diameter(i_spec), vs, env_state%temp, ustar, gamma, A, &
          alpha, lambda)
     vd = calculate_vd(aer_res_a, rs, vs)
     vd_rate(i_spec) = vd/env_state%height
     else
       vd_rate(i_spec) = 0.0d0
     end if
  end do
  ! Compute analytical solution

  print*,vd_rate

  open(unit=out_unit, file=out_name)

  time = 0d0
  ntime = int(run_part_opt%t_max / run_part_opt%t_output)
  print*, 'time steps = ', ntime
  write(out_unit,*) time, 0.0, 0.0, total_initial_mass

  allocate(total_mass(aero_data%n_spec))

  do i = 1,ntime
     time =  time + run_part_opt%t_output

     do i_spec = 1, aero_data%n_spec
        total_mass(i_spec) = total_initial_mass(i_spec) * exp(-vd_rate(i_spec) * time)
     end do

     write(out_unit,*) time, 0.0, 0.0, total_mass

  end do

  close(out_unit)

  ! Deallocates
  call gas_data_deallocate(gas_data)
  call gas_state_deallocate(gas_state_init)
  call aero_data_deallocate(aero_data)
  call aero_dist_deallocate(aero_dist_init)
  call env_state_deallocate(env_state)
  call scenario_deallocate(scenario)

  deallocate(diameter)
  deallocate(vd_rate)
  deallocate(particle_mass)
  deallocate(total_mass)
  deallocate(density)
  deallocate(num_conc)
  deallocate(total_initial_mass)

end program test_deposition_exact
