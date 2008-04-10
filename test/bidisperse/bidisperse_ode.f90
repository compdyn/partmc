! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Compute the evolution of a bidisperse distribution with the
! sedimentation kernel.
!
! This is a complete hack and is not supposed to be general or
! re-usable.
! 
! The initial distribution consists of n_small small particles of
! size v_small, and one big particle of size v_big. The
! sedimentation kernel is zero between same sized particles, so all
! that happens is the number of small particles decreases (but the
! remaining ones keep the initial volume) while the big particle
! remains just one particle but grows in volume. This is thus really
! a one-dimensional ODE which we treat as being defined in terms of
! the current number of small particles.

program bidisperse_ode
  
  use pmc_kernel_sedi
  use pmc_env_state
  use pmc_util
  use pmc_bin_grid
  
  ! Volume of one small particle (m^3).
  real*8, parameter :: v_small = 0.3967039840458671239139981d-14
  !> Initial volume of big particle (m^3).
  real*8, parameter :: v_big_init = 0.4188790204786393454605195d-11
  !> Init number of small particles.
  real*8, parameter :: n_small_init = 10000d0
  !> Total simulation time.
  real*8, parameter :: t_max = 600d0
  !> Timestep.
  real*8, parameter :: del_t = 0.001d0
  !> How often to print progress.
  real*8, parameter :: t_progress = 10d0
  !> Particle number conc (#/m^3).
  real*8, parameter :: num_den_small = 1d9
  !> Number of bins.
  integer, parameter :: n_bin = 250
  !> Minimum bin radius (m).
  real*8, parameter :: bin_r_min = 1d-8
  !> Minimum bin radius (m).
  real*8, parameter :: bin_r_max = 1d0
  !> Scale factor for bins.
  integer, parameter :: scal = 3
  !> Output unit number.
  integer, parameter :: out_unit = 33
  !> Output filename.
  character(len=*), parameter :: out_name = "out/bidisperse_ode_counts.d"
  
  type(env_state_t) :: env_state
  integer :: i_step, n_step
  real*8 :: comp_vol, n_small, time, v_big, num_den
  type(bin_grid_t) :: bin_grid

  num_den = num_den_small * (n_small_init + 1d0) / n_small_init
  comp_vol = (n_small_init + 1d0) / num_den
  call bin_grid_make(bin_grid, n_bin, rad2vol(bin_r_min), rad2vol(bin_r_max))

  open(unit=out_unit, file=out_name)
  time = 0d0
  n_small = n_small_init
  n_step = nint(t_max / del_t) + 1
  v_big = v_big_init + (n_small_init - n_small) * v_small
  write(*,'(a8,a14,a14,a9)') &
       't', 'n_small', 'v_big', 'n_coag'
  write(*,'(f8.1,e14.5,e14.5,f9.2)') &
       time, n_small / comp_vol, v_big / comp_vol, &
       n_small_init - n_small
  write(out_unit,'(e20.10,e20.10,e20.10)') &
       time, n_small / comp_vol, v_big / comp_vol
  do i_step = 2,n_step
     time = dble(i_step - 1) * del_t
     call bidisperse_step(v_small, v_big_init, n_small_init, &
          env_state, comp_vol, del_t, n_small)
     v_big = v_big_init + (n_small_init - n_small) * v_small
     if (mod(i_step - 1, nint(t_progress / del_t)) .eq. 0) then
        write(*,'(a8,a14,a14,a9)') &
             't', 'n_small', 'v_big', 'n_coag'
        write(*,'(f8.1,e14.5,e14.5,f9.2)') &
             time, n_small / comp_vol, v_big / comp_vol, &
             n_small_init - n_small
        write(out_unit,'(e20.10,e20.10,e20.10)') &
             time, n_small / comp_vol, v_big / comp_vol
     end if
  end do

  close(out_unit)
  call bin_grid_free(bin_grid)
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bidisperse_f(n_small, v_small, v_big_init, &
       n_small_init, env_state, comp_vol, n_small_dot)
    
    use pmc_env_state

    !> Current number of small particles.
    real*8, intent(in) :: n_small
    !> Volume of one small particle.
    real*8, intent(in) :: v_small
    !> Initial volume of the big particle.
    real*8, intent(in) :: v_big_init
    !> Initial number of small particles.
    real*8, intent(in) :: n_small_init
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Computational volume (m^3).
    real*8, intent(in) :: comp_vol
    !> Derivative of n_small.
    real*8, intent(out) :: n_small_dot

    integer :: n_spec
    real*8 :: v_big, k
    type(aero_particle_t) :: aero_particle_1, aero_particle_2
    type(aero_data_t) :: aero_data
    
    v_big = v_big_init + (n_small_init - n_small) * v_small
    n_spec = 1
    call aero_data_alloc(aero_data, n_spec)
    call aero_particle_alloc(aero_particle_1, n_spec)
    call aero_particle_alloc(aero_particle_2, n_spec)
    aero_particle_1%vol(1) = v_small
    aero_particle_2%vol(1) = v_big
    call kernel_sedi(aero_particle_1, aero_particle_2, aero_data, env_state, k)
    call aero_particle_free(aero_particle_1)
    call aero_particle_free(aero_particle_2)
    call aero_data_free(aero_data)
    n_small_dot = - k / comp_vol * n_small
    
  end subroutine bidisperse_f
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bidisperse_step(v_small, v_big_init, n_small_init, &
       env_state, comp_vol, del_t, n_small)
    
    use pmc_env_state

    !> Volume of one small particle.
    real*8, intent(in) :: v_small
    !> Initial volume of the big particle.
    real*8, intent(in) :: v_big_init
    !> Initial number of small particles.
    real*8, intent(in) :: n_small_init
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Computational volume (m^3).
    real*8, intent(in) :: comp_vol
    !> Timestep.
    real*8, intent(in) :: del_t
    !> Current number of small particles.
    real*8, intent(inout) :: n_small
    
    real*8 n_small_dot, k1, k2, k3, k4
    
    ! integrate ODE with Runge-Kutta-4
    
    call bidisperse_f(n_small, &
         v_small, v_big_init, n_small_init, env_state, comp_vol, n_small_dot)
    k1 = del_t * n_small_dot
    
    call bidisperse_f(n_small + k1/2d0, &
         v_small, v_big_init, n_small_init, env_state, comp_vol, n_small_dot)
    k2 = del_t * n_small_dot
    
    call bidisperse_f(n_small + k2/2d0, &
         v_small, v_big_init, n_small_init, env_state, comp_vol, n_small_dot)
    k3 = del_t * n_small_dot
    
    call bidisperse_f(n_small + k3, &
         v_small, v_big_init, n_small_init, env_state, comp_vol, n_small_dot)
    k4 = del_t * n_small_dot
    
    n_small = n_small + k1/6d0 + k2/3d0 + k3/3d0 + k4/6d0
    
  end subroutine bidisperse_step
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program bidisperse_ode
