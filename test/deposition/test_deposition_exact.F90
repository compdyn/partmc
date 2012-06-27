program test_deposition_exact

  !use pmc_deposition
  !use pmc_constants

  !> Output unit number.
  integer, parameter :: out_unit = 34
  !> Output filename.
  character(len=*), parameter :: out_name = "out/deposition_data.txt"
  character(len=100) :: in_prefix
!  real(kind=dp) :: diameter
!  real(kind=dp) :: density
  integer :: optind

  ! All the PartMC structures we want values from
!  aero_data
!  aero_dist_init
!  scenario
!  env_state

  ! We need the specfile
  optind = 1
  call get_command_argument(optind, in_prefix)
  print*,in_prefix

!  call spec_file_read_string(file, 'aerosol_data', sub_filename)
!  call spec_file_open(sub_filename, sub_file)
!  call spec_file_read_aero_data(sub_file, aero_data)
!  call spec_file_close(sub_file)
!
!  call spec_file_read_string(file, 'aerosol_init', sub_filename)
!  call spec_file_open(sub_filename, sub_file)
!  call spec_file_read_aero_dist(sub_file, aero_data, aero_dist_init)
!  call spec_file_close(sub_file)
!
!!    call spec_file_read_string(file, "temp_profile", sub_filename)
!    call spec_file_open(sub_filename, sub_file)
!!    call spec_file_read_timed_real_array(sub_file, "temp", &
!         scenario%temp_time, scenario%temp)
!    call spec_file_close(sub_file)
!    call spec_file_read_string(file, "height_profile", sub_filename)
!    call spec_file_open(sub_filename, sub_file)
!    call spec_file_read_timed_real_array(sub_file, "height", &
!         scenario%height_time, scenario%height)
!    call spec_file_close(sub_file)
!
!  ! Read in spec file stuff
!!  ! del_t
!  !call spec_file_read_real(file, 'del_t', run_part_opt%del_t)
!  !call spec_file_read_real(file, 't_output', run_part_opt%t_output)
!  !call spec_file_read_real(file, 't_max', run_part_opt%t_max)
!  ! diameter
!  ! species so we have density
!!  ! Ideal number concentration so we can find mass
!
!  diameter = 1.0d0
!  density = 1.0d0
!  ! Other things
!  aer_res_a = 0.0d0
!  ustar = 1.0d0
!  gamma = 1.0d0
!  A = 1.0d0
!  alpha = 1.0d0
!  call env_state_allocate(env_state)
!  call spec_file_read_env_state(file, env_state)
!
!  rho_air = env_state_air_den(env_state)
!!  viscosd = env_state_air_dynamic_viscosity(env_state)
!  viscosk = viscosd / rho_air
!  lambda = env_state_air_mean_free_path(env_state)
!  vs = calculate_vs(diameter, density, viscosd, lambda)
!!  rs = calculate_rs(diameter, vs, env_state%temp, ustar, gamma, A, &
!       alpha, lambda)
!  vd = calculate_vd(aer_res_a, rs, vs)
!
!  ! Compute analytical solution
!   
!  ntime = 10
!  open(unit=out_unit, file=out_name)
!  t = 0d0
!  do i = 0,ntime-1
!     t = i*del_t
!     c = c*exp(-vd*t)
!     write(30,*) t, c
!!  end do
!  close(out_unit)
!
end program test_deposition_exact
