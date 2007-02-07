! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Compute the evolution of a bidisperse distribution with the
! sedimentation kernel.
! 
! The initial distribution consists of n_small small particles of
! size v_small, and one big particle of size v_big. The
! sedimentation kernel is zero between same sized particles, so all
! that happens is the number of small particles decreases (but the
! remaining ones keep the initial volume) while the big particle
! remains just one particle but grows in volume. This is thus really
! a one-dimensional ODE which we treat as being defined in terms of
! the current number of small particles.

program run_sedi_ode
  
  use mod_kernel_sedi
  use mod_environ
  
  real*8 v_small, v_big_init, n_small_init, del_t, t_max, N_0
  integer scal
  type(environ) :: env                    ! environment state

  parameter (v_small = 4.8531435d-15)     ! volume of one small particle
  parameter (v_big_init = 3.94438917d-12) ! init volume of the big particle
  parameter (n_small_init = 9999d0)       ! init number of small particles
  parameter (del_t = 0.001d0)             ! timestep
  parameter (t_max = 600d0)               ! total simulation time
  parameter (N_0 = 1d9)                   ! particle number concentration
  parameter (scal = 3)                    ! scale factor for bins
  
  integer i_step, n_step
  real*8 n_small, time, v_big, dlnr
  
  env%V_comp = dble(n_small_init + 1d0) / N_0
  dlnr = dlog(2d0) / (3d0 * dble(scal))
  
  time = 0d0
  n_small = n_small_init
  n_step = nint(t_max / del_t) + 1
  v_big = v_big_init + (n_small_init - n_small) * v_small
  write(*,'(a8,a14,a14,a9)') &
       't', 'n_small', 'v_big', 'n_coag'
  write(*,'(f8.1,e14.5,e14.5,f9.2)') &
       time, n_small / env%V_comp / dlnr, v_big / env%V_comp / dlnr, &
       n_small_init - n_small
  do i_step = 1,n_step
     time = dble(i_step - 1) * del_t
     call bidisperse_step(v_small, v_big_init, n_small_init, &
          env, del_t, n_small)
     v_big = v_big_init + (n_small_init - n_small) * v_small
     if (mod(i_step - 1, nint(1d0 / del_t)) .eq. 0) then
        write(*,'(a8,a14,a14,a9)') &
             't', 'n_small', 'v_big', 'n_coag'
        write(*,'(f8.1,e14.5,e14.5,f9.2)') &
             time, n_small / env%V_comp / dlnr, v_big / env%V_comp / dlnr, &
             n_small_init - n_small
     end if
  end do
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bidisperse_f(n_small, v_small, v_big_init, &
       n_small_init, env, n_small_dot)
    
    use mod_environ

    real*8, intent(in) :: n_small        !  current number of small particles
    real*8, intent(in) :: v_small        !  volume of one small particle
    real*8, intent(in) :: v_big_init     !  initial volume of the big particle
    real*8, intent(in) :: n_small_init   !  initial number of small particles
    type(environ), intent(in) :: env        ! environment state
    real*8, intent(out) :: n_small_dot    !  derivative of n_small
    
    real*8 v_big, k
    
    v_big = v_big_init + (n_small_init - n_small) * v_small
    call kernel_sedi(v_small, v_big, env, k)
    n_small_dot = - (k * 1d0/env%V_comp * n_small)
    
  end subroutine bidisperse_f
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine bidisperse_step(v_small, v_big_init, n_small_init, &
       env, del_t, n_small)
    
    use mod_environ

    real*8, intent(in) :: v_small        !  volume of one small particle
    real*8, intent(in) :: v_big_init     !  initial volume of the big particle
    real*8, intent(in) :: n_small_init   !  initial number of small particles
    real*8, intent(in) :: del_t          !  timestep
    type(environ), intent(in) :: env        ! environment state
    real*8, intent(inout) :: n_small        !  current number of small particles
    
    real*8 n_small_dot, k1, k2, k3, k4
    
    ! integrate ODE with Runge-Kutta-4
    
    call bidisperse_f(n_small, &
         v_small, v_big_init, n_small_init, env, n_small_dot)
    k1 = del_t * n_small_dot
    
    call bidisperse_f(n_small + k1/2d0, &
         v_small, v_big_init, n_small_init, env, n_small_dot)
    k2 = del_t * n_small_dot
    
    call bidisperse_f(n_small + k2/2d0, &
         v_small, v_big_init, n_small_init, env, n_small_dot)
    k3 = del_t * n_small_dot
    
    call bidisperse_f(n_small + k3, &
         v_small, v_big_init, n_small_init, env, n_small_dot)
    k4 = del_t * n_small_dot
    
    n_small = n_small + k1/6d0 + k2/3d0 + k3/3d0 + k4/6d0
    
  end subroutine bidisperse_step
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program run_sedi_ode
