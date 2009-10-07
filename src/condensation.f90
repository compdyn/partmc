! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Copyright (C) 2009 Joseph Ching
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_condensation module.

!> Water condensation onto aerosol particles.
module pmc_condensation

  use pmc_aero_state
  use pmc_bin_grid
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
  use dvode_f90_m
  
  logical, parameter :: PMC_COND_CHECK_DERIVS = .true.
  real(kind=dp), parameter :: PMC_COND_CHECK_EPS = 1d-8
  real(kind=dp), parameter :: PMC_COND_CHECK_REL_TOL = 5d-2

  integer, parameter :: PMC_COND_N_REAL_PARAMS = 11
  integer, parameter :: PMC_COND_N_INTEGER_PARAMS = 0

  integer, parameter :: PMC_COND_U = 1
  integer, parameter :: PMC_COND_V = 2
  integer, parameter :: PMC_COND_W = 3
  integer, parameter :: PMC_COND_X = 4
  integer, parameter :: PMC_COND_Y = 5
  integer, parameter :: PMC_COND_Z = 6
  integer, parameter :: PMC_COND_K_A = 7
  integer, parameter :: PMC_COND_D_V = 8
  integer, parameter :: PMC_COND_KAPPA = 9
  integer, parameter :: PMC_COND_V_DRY = 10
  integer, parameter :: PMC_COND_S = 11

  !> Diameter at which the pmc_cond_saved_delta_star was evaluated.
  real(kind=dp), save :: pmc_cond_saved_diameter
  !> Saved value of delta_star evaluated at pmc_cond_saved_diameter.
  real(kind=dp), save :: pmc_cond_saved_delta_star

  real(kind=dp), save :: save_real_params(PMC_COND_N_REAL_PARAMS)
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> vapor.
  subroutine condense_particles(bin_grid, env_state, aero_data, &
       aero_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    
    integer :: i_bin, j, new_bin, k
    real(kind=dp) :: pv, pre_water, post_water
    type(aero_particle_t), pointer :: particle
    real(kind=dp) :: kappa

    pre_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          particle => aero_state%bin(i_bin)%particle(j)
          pre_water = pre_water + particle%vol(aero_data%i_water)
       end do
    end do

    !write(*,*) 'RH (before condense) =', env_state%rel_humid

    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          !write(*,*) 'i_bin= ', i_bin, 'j =',  j
          call condense_particle_vode(del_t, env_state, aero_data, &
               aero_state%bin(i_bin)%particle(j))
!          write(*,*) 'stop after the first particle'
!          STOP
       end do
    end do

    ! We resort the particles in the bins after only all particles
    ! have condensation done, otherwise we will lose track of which
    ! ones have had condensation and which have not.
    call aero_state_resort(bin_grid, aero_state)

    post_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          particle => aero_state%bin(i_bin)%particle(j)
          post_water = post_water + particle%vol(aero_data%i_water)
       end do
    end do

    !write(*,*) 'pre_water (BF cond) =', pre_water
    !write(*,*) 'post_water (AF cond) =', post_water

    !write(*,*) 'RH (before update due to cond. ) =', env_state%rel_humid

    ! update the environment due to condensation of water
    call env_state_change_water_volume(env_state, aero_data, &
         (post_water - pre_water) / aero_state%comp_vol)

    !write(*,*) 'RH (after update due to cond. ) =', env_state%rel_humid 

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the Kappa value of individual particle based on volume
  !> fractions
  subroutine compute_kappa(aero_data, aero_particle, tot_kappa)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Kappa values for individual particle. 
    real(kind=dp), intent(out) :: tot_kappa
    !> Quantity to total. 
!    real(kind=dp), intent(in) :: quantity(:)

    ! local variables
    real(kind=dp) :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp) :: g_water, g_solute
    real(kind=dp) :: indiv_kappa
    real(kind=dp) :: tot_vol, vol_frac, tot_vol_frac
    real(kind=dp) :: water_vol,other_vol
    integer :: j, k
    
       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       tot_kappa=0d0
       tot_vol=0d0
       tot_vol_frac=0d0

       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
          ! total volume of solute (non-water components)
          tot_vol=tot_vol+aero_particle%vol(k)
!          write(*,*) 'k=',k ,'vol= ', aero_particle%vol(k)
          endif
       enddo
       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
              if(aero_data%kappa(k) > 500d0) then
          indiv_kappa=(M_water/rho_water)*(aero_data%density(k)/aero_data%molec_weight(k))    !volume fraction
!          write(*,*) 'k=', k, 'mole weight =', M_water, 'den =', rho_water
!          write(*,*) 'k=', k, 'aero den =', aero_data%density(k), 'aero weight=', aero_data%molec_weight(k)   
!          write(*,*) 'inorganic', k, indiv_kappa
              else
          indiv_kappa=aero_data%kappa(k)
!          write(*,*) 'organic', k, indiv_kappa
              endif

          vol_frac=aero_particle%vol(k)/tot_vol
!          write(*,*) 'vol= ', vol_frac, 'kappa= ', indiv_kappa
          tot_kappa=tot_kappa+vol_frac*indiv_kappa
          tot_vol_frac=tot_vol_frac+vol_frac
          endif
       enddo
     
!       write(*,*) 'total vol frac=', tot_vol_frac
 
  end subroutine compute_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the condensation growth or decay ODE for total time
  !> del_t for a single particle.
  subroutine condense_particle(del_t, env_state, aero_data, &
       aero_particle)

    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp) :: time_step, time, pre_water, post_water
    logical :: done

    !write(*,*) 'in condense_particle '

    time = 0d0
    done = .false.
    do while (.not. done)
       !write(*,*) 'in do loop in condense_particle', 'del_t=', del_t , 'time=', time 
       !write(*,*)  'time_step=', time_step, 'done=', done
       call condense_step_euler(del_t - time, time_step, done, env_state, &
            aero_data, aero_particle)
       time = time + time_step
    end do
    
  end subroutine condense_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the condensation growth or decay ODE for total time
  !> del_t for a single particle.
  subroutine condense_particle_vode(del_t, env_state, aero_data, &
       aero_particle)

    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    integer, parameter :: real_work_len = 33 ! for scalar equation
    integer, parameter :: integer_work_len = 31 ! for scalar equation

    integer :: n_eqn, tol_types, itask, istate, iopt, method_flag
    real(kind=dp) :: val(1), init_time, final_time, rel_tol(1)
    real(kind=dp) :: abs_tol(1), real_work(real_work_len)
    integer :: integer_work(integer_work_len)
    real(kind=dp) :: real_params(PMC_COND_N_REAL_PARAMS)
    integer :: integer_params(PMC_COND_N_INTEGER_PARAMS)
    real(kind=dp) :: rho_water, M_water
    real(kind=dp) :: thermal_conductivity, molecular_diffusion
    type(vode_opts) :: options

!    interface
!       ! interfaces copied by hand from vode.f
!       subroutine dvode(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, &
!            ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, &
!            MF, RPAR, IPAR)
!         integer :: LRW, LIW
!         double precision :: Y(:), T, TOUT, RTOL(:), ATOL(:), &
!              RWORK(LRW), RPAR(:)
!         integer :: NEQ, ITOL, ITASK, ISTATE, IOPT, IWORK(LIW), &
!              MF, IPAR(:)
!         interface
!            subroutine F(NEQ, T, Y, YDOT, RPAR, IPAR)
!              integer :: NEQ
!              double precision :: T, Y(NEQ), YDOT(NEQ), RPAR(:)
!              integer :: IPAR(:)
!            end subroutine F
!            subroutine JAC(NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
!              integer :: NEQ, NROWPD
!              double precision :: T, Y(NEQ), PD(NROWPD,NEQ), RPAR(:)
!              integer :: ML, MU, IPAR(:)
!            end subroutine JAC
!         end interface
!       end subroutine dvode
!    end interface

    thermal_conductivity = 1d-3 * (4.39d0 + 0.071d0 &
         * env_state%temp) ! FIXME: supposedly in J m^{-1} s^{-1} K^{-1}
    molecular_diffusion = 0.211d-4 &
         / (env_state%pressure / const%air_std_press) &
         * (env_state%temp / 273d0)**1.94d0 ! FIXME: supposedly in m^2 s^{-1}
    rho_water = aero_particle_water_density(aero_data)
    M_water = aero_particle_water_molec_weight(aero_data)

    real_params(PMC_COND_U) = const%water_latent_heat * rho_water &
         / (4d0 * env_state%temp)
    real_params(PMC_COND_V) = 4d0 * M_water &
         * env_state_sat_vapor_pressure(env_state) &
         / (rho_water * const%univ_gas_const * env_state%temp)
    real_params(PMC_COND_W) = const%water_latent_heat * M_water &
         / (const%univ_gas_const * env_state%temp)
    real_params(PMC_COND_X) = 4d0 * M_water * const%water_surf_eng &
         / (const%univ_gas_const * env_state%temp * rho_water) 
    real_params(PMC_COND_Y) = 2d0 * thermal_conductivity &
         / (const%accom_coeff * env_state_air_den(env_state) &
         * const%water_spec_heat) &
         * sqrt(2d0 * const%pi * const%air_molec_weight &
         / (const%univ_gas_const * env_state%temp))
    real_params(PMC_COND_Z) = 2d0 * molecular_diffusion &
         / const%accom_coeff * sqrt(2d0 * const%pi * M_water &
         / (const%univ_gas_const * env_state%temp))
    real_params(PMC_COND_K_A) = thermal_conductivity
    real_params(PMC_COND_D_V) = molecular_diffusion
    real_params(PMC_COND_KAPPA) &
         = aero_particle_solute_kappa(aero_particle, aero_data)
    real_params(PMC_COND_V_DRY) &
         = aero_particle_solute_volume(aero_particle, aero_data)
    real_params(PMC_COND_S) = env_state%rel_humid - 1d0

    ! set VODE inputs
    n_eqn = 1
    val(1) = vol2diam(aero_particle_volume(aero_particle))
    init_time = 0d0
    final_time = del_t
    tol_types = 1 ! both rel_tol and abs_tol are scalars
    rel_tol(1) = 1d-6
    abs_tol(1) = val(1) * 1d-6
    itask = 1 ! just output val at final_time
    istate = 1 ! first call for this ODE
    iopt = 0 ! no optional inputs
    real_work = 0d0
    integer_work = 0
    method_flag = 21 ! stiff (BDF) method, user-supplied full Jacobian

    ! call ODE integrator
    !call dvode(condense_vode_f, n_eqn, val, init_time, &
    !     final_time, tol_types, rel_tol, abs_tol, itask, istate, &
    !     iopt, real_work, real_work_len, integer_work, integer_work_len, &
    !     condense_vode_jac, method_flag, real_params, integer_params)
    !if (istate /= 2) then
    !   call die_msg(982335370, "DVODE error code: " &
    !        // trim(integer_to_string(istate)))
    !end if

    options = set_opts(dense_j = .true., abserr = abs_tol(1), relerr = rel_tol(1), &
         user_supplied_jacobian = .true.)
    save_real_params = real_params
    call dvode_f90(condense_vode_f, n_eqn, val, init_time, final_time, &
         itask, istate, options, j_fcn = condense_vode_jac)

    ! translate output back to particle
    aero_particle%vol(aero_data%i_water) = diam2vol(val(1)) &
         - aero_particle_solute_volume(aero_particle, aero_data)

    ! ensure volumes stay positive
    aero_particle%vol(aero_data%i_water) = max(0d0, &
         aero_particle%vol(aero_data%i_water))
    
  end subroutine condense_particle_vode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function \$ \dot{D} (D) \$.
  subroutine condense_vode_f(n_eqn, time, state, state_dot)
    !, real_params, &
    !   integer_params)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Time derivative of state vector.
    real(kind=dp), intent(out) :: state_dot(n_eqn)
    !> Real parameters.
    !real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    !integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: real_params(PMC_COND_N_REAL_PARAMS)
    integer :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: delta, diameter, k_ap

    real_params = save_real_params

    call assert(383442283, n_eqn == 1)
    delta = 0d0
    diameter = state(1)
    write(*,*) 'condense_vode_f: t,D = ', time, diameter
    call condense_vode_delta_star_newton(delta, diameter, real_params, &
         integer_params)
    k_ap = corrected_thermal_conductivity(diameter, real_params, &
         integer_params)
    state_dot(1) = k_ap * delta / (real_params(PMC_COND_U) * diameter)

    ! save value of delta to be used by condense_vode_jac
    pmc_cond_saved_diameter = diameter
    pmc_cond_saved_delta_star = delta

  end subroutine condense_vode_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Derivative \$ \partial \dot{D}(D) / \partial D \$.
  subroutine condense_vode_jac(n_eqn, time, state, lower_band_width, &
       upper_band_width, state_jac, n_row_jac)
    !, real_params, integer_params)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Lower band width for banded Jacobian (unused).
    integer, intent(in) :: lower_band_width
    !> Upper band width for banded Jacobian (unused).
    integer, intent(in) :: upper_band_width
    !> Number of rows in the Jacobian.
    integer, intent(in) :: n_row_jac
    !> Jacobian of time derivative of state vector.
    real(kind=dp), intent(out) :: state_jac(n_row_jac, n_eqn)
    !> Real parameters.
    !real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    !integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: real_params(PMC_COND_N_REAL_PARAMS)
    integer :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: delta_star, D, U, k_ap, dkap_dD
    real(kind=dp) :: dh_dD, dh_ddelta, ddeltastar_dD

    real_params = save_real_params

    call assert(973429886, n_eqn == 1)
    call assert(230741766, n_row_jac == 1)

    U = real_params(PMC_COND_U)
    D = state(1)
    call assert_msg(232625737, D == pmc_cond_saved_diameter, &
         "state diameter does not match pmc_cond_saved_diameter")
    delta_star = pmc_cond_saved_delta_star
    k_ap = corrected_thermal_conductivity(D, real_params, &
         integer_params)
    dkap_dD = corrected_thermal_conductivity_deriv(D, real_params, &
         integer_params)
    dh_dD = condense_vode_implicit_dh_dD(delta_star, D, real_params, &
         integer_params)
    dh_ddelta = condense_vode_implicit_dh_ddelta(delta_star, D, &
         real_params, integer_params)
    ddeltastar_dD = - dh_dD / dh_ddelta
    state_jac(1,1) = dkap_dD * delta_star / (U * D) &
         + k_ap * ddeltastar_dD / (U * D) &
         - k_ap * delta_star / (U * D**2)

  end subroutine condense_vode_jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scalar Newton's method for solving the implicit delta_star
  !> function.
  subroutine condense_vode_delta_star_newton(delta, diameter, real_params, &
       integer_params)
    
    !> Variable that we are solving for.
    real(kind=dp), intent(inout) :: delta
    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    !> Delta convergence tolerance.
    real(kind=dp), parameter :: delta_tol = 1d-16
    !> Function convergence tolerance.
    real(kind=dp), parameter :: h_tol = 1d-10
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100

    integer :: iter
    real(kind=dp) :: h, dh_ddelta, old_h, delta_step, h_step
    real(kind=dp) :: check_delta, check_h, finite_diff_dh, rel_error

    h = condense_vode_implicit_h(delta, diameter, real_params, &
         integer_params)
    dh_ddelta = condense_vode_implicit_dh_ddelta(delta, diameter, &
         real_params, integer_params)
    old_h = h

    iter = 0
    do
       iter = iter + 1
       delta_step = - h / dh_ddelta
       write(*,'(a5,a15,a15,a15,a15)') 'iter', 'delta', 'h', 'dh_ddelta', 'delta_step'
       write(*,'(i5,e15.5,e15.5,e15.5,e15.5)') iter, delta, h, dh_ddelta, delta_step
       
       delta = delta + delta_step
       h = condense_vode_implicit_h(delta, diameter, real_params, &
            integer_params)
       dh_ddelta = condense_vode_implicit_dh_ddelta(delta, diameter, &
            real_params, integer_params)
       h_step = h - old_h
       old_h = h

       if (PMC_COND_CHECK_DERIVS) then
          check_delta = delta * (1d0 + PMC_COND_CHECK_EPS)
          check_h = condense_vode_implicit_h(check_delta, diameter, real_params, &
               integer_params)
          finite_diff_dh = (check_h - h) / (check_delta - delta)
          rel_error = abs(finite_diff_dh - dh_ddelta) &
               / (abs(finite_diff_dh) + abs(dh_ddelta))
          if (rel_error > PMC_COND_CHECK_REL_TOL) then
             write(0,*) "ERROR: cond_newt: incorrect derivative"
             write(0,*) "delta ", delta
             write(0,*) "check_delta ", check_delta
             write(0,*) "delta_delta ", (check_delta - delta)
             write(0,*) "h ", h
             write(0,*) "check_h ", check_h
             write(0,*) "delta_h ", (check_h - h)
             write(0,*) "dh_ddelta ", dh_ddelta
             write(0,*) "finite_diff_dh ", finite_diff_dh
             write(0,*) "rel_error ", rel_error
          end if
       end if

       if (iter .ge. iter_max) then
          call die_msg(136296873, 'Newton iteration failed to terminate')
       end if
       
       ! FIXME: gfortran 4.1.1 requires the "then" in the following
       ! statement, rather than using a single-line "if" statement.
       if ((abs(delta_step) .lt. delta_tol) &
            .and. (abs(h_step) .lt. h_tol)) then
          exit
       end if
    end do
 
  end subroutine condense_vode_delta_star_newton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function \$ h(\delta,D) \$.
  real(kind=dp) function condense_vode_implicit_h(delta, diameter, &
       real_params, integer_params)

    !> Growth parameter (units???).
    real(kind=dp), intent(in) :: delta
    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: D, U, V, W, X, S, k_ap, D_vp, a_w

    D = diameter
    U = real_params(PMC_COND_U)
    V = real_params(PMC_COND_V)
    W = real_params(PMC_COND_W)
    X = real_params(PMC_COND_X)
    S = real_params(PMC_COND_S)
    k_ap = corrected_thermal_conductivity(diameter, real_params, &
         integer_params)
    D_vp = corrected_molecular_diffusion(diameter, real_params, &
         integer_params)
    a_w = water_activity(diameter, real_params, integer_params)

    condense_vode_implicit_h = k_ap * delta - U * V * D_vp &
         * (S - a_w / (1d0 + delta) * exp(W * delta / (1d0 + delta) &
         + (X / D) / (1d0 + delta)))

  end function condense_vode_implicit_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Derivative \$ \partial h(\delta,D) / \partial \delta \$.
  real(kind=dp) function condense_vode_implicit_dh_ddelta(delta, diameter, &
       real_params, integer_params)

    !> Growth parameter (units???).
    real(kind=dp), intent(in) :: delta
    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: D, U, V, W, X, S, k_ap, D_vp, a_w
    real(kind=dp) :: dkap_dD, DDvp_dD, daw_dD

    D = diameter
    U = real_params(PMC_COND_U)
    V = real_params(PMC_COND_V)
    W = real_params(PMC_COND_W)
    X = real_params(PMC_COND_X)
    S = real_params(PMC_COND_S)
    k_ap = corrected_thermal_conductivity(diameter, real_params, &
         integer_params)
    D_vp = corrected_molecular_diffusion(diameter, real_params, &
         integer_params)
    a_w = water_activity(diameter, real_params, integer_params)
    dkap_dD = corrected_thermal_conductivity_deriv(diameter, real_params, &
         integer_params)
    dDvp_dD = corrected_molecular_diffusion_deriv(diameter, real_params, &
         integer_params)
    daw_dD = water_activity_deriv(diameter, real_params, integer_params)

    condense_vode_implicit_dh_ddelta = &
         k_ap - U * V * D_vp * a_w / (1d0 + delta)**2 &
         * (1d0 - W / (1d0 + delta) + (X / D) / (1d0 + delta)) &
         * exp(W * delta / (1d0 + delta) + (X / D) / (1d0 + delta))
    
  end function condense_vode_implicit_dh_ddelta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Derivative \$ \partial h(\delta,D) / \partial D \$.
  real(kind=dp) function condense_vode_implicit_dh_dD(delta, diameter, &
       real_params, integer_params)

    !> Growth parameter (units???).
    real(kind=dp), intent(in) :: delta
    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: D, U, V, W, X, S, k_ap, D_vp, a_w
    real(kind=dp) :: dkap_dD, DDvp_dD, daw_dD

    D = diameter
    U = real_params(PMC_COND_U)
    V = real_params(PMC_COND_V)
    W = real_params(PMC_COND_W)
    X = real_params(PMC_COND_X)
    S = real_params(PMC_COND_S)
    k_ap = corrected_thermal_conductivity(diameter, real_params, &
         integer_params)
    D_vp = corrected_molecular_diffusion(diameter, real_params, &
         integer_params)
    a_w = water_activity(diameter, real_params, integer_params)
    dkap_dD = corrected_thermal_conductivity_deriv(diameter, real_params, &
         integer_params)
    dDvp_dD = corrected_molecular_diffusion_deriv(diameter, real_params, &
         integer_params)
    daw_dD = water_activity_deriv(diameter, real_params, integer_params)

    condense_vode_implicit_dh_dD = dkap_dD * delta &
         - U * V * dDvp_dD * S + U * V * (a_w * dDvp_dD + D_vp * daw_dD &
         - D_vp * a_w * (X / D**2) / (1d0 + delta)) * (1d0 / (1d0 + delta)) &
         * exp((W * delta) / (1d0 + delta) + (X / D) / (1d0 + delta))
    
  end function condense_vode_implicit_dh_dD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function \$ k_{\rm a}'(D) \$.
  real(kind=dp) function corrected_thermal_conductivity(diameter, &
       real_params, integer_params)

    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: k_a, Y, D

    k_a = real_params(PMC_COND_K_A)
    Y = real_params(PMC_COND_Y)
    D = diameter

    corrected_thermal_conductivity = k_a / (1d0 + Y / D)

  end function corrected_thermal_conductivity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Derivative \$ \partial k_{\rm a}'(D) / \partial D \$.
  real(kind=dp) function corrected_thermal_conductivity_deriv(diameter, &
       real_params, integer_params)

    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: k_a, Y, D

    k_a = real_params(PMC_COND_K_A)
    Y = real_params(PMC_COND_Y)
    D = diameter

    corrected_thermal_conductivity_deriv = k_a * Y / (D + Y)**2

  end function corrected_thermal_conductivity_deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function \$ D_{\rm v}'(D) \$.
  real(kind=dp) function corrected_molecular_diffusion(diameter, &
       real_params, integer_params)

    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: D_v, Z, D

    D_v = real_params(PMC_COND_D_V)
    Z = real_params(PMC_COND_Z)
    D = diameter

    corrected_molecular_diffusion = D_v / (1d0 + Z / D)

  end function corrected_molecular_diffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Derivative \$ \partial D_{\rm v}'(D) / \partial D \$.
  real(kind=dp) function corrected_molecular_diffusion_deriv(diameter, &
       real_params, integer_params)

    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: D_v, Z, D

    D_v = real_params(PMC_COND_D_V)
    Z = real_params(PMC_COND_Z)
    D = diameter

    corrected_molecular_diffusion_deriv = D_v * Z / (D + Z)**2

  end function corrected_molecular_diffusion_deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function \$ a_{\rm w}(D) \$.
  real(kind=dp) function water_activity(diameter, real_params, &
       integer_params)

    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: D, kappa, V_dry

    D = diameter
    kappa = real_params(PMC_COND_KAPPA)
    V_dry = real_params(PMC_COND_V_DRY)

    water_activity = (const%pi / 6d0 * D**3 - V_dry) &
         / (const%pi / 6d0 * D**3 + (kappa - 1d0) * V_dry)

  end function water_activity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Derivative \$ \partial a_{\rm w}(D) / \partial D \$.
  real(kind=dp) function water_activity_deriv(diameter, real_params, &
       integer_params)

    !> Current diameter (m)
    real(kind=dp), intent(in) :: diameter
    !> Real parameters.
    real(kind=dp), intent(in) :: real_params(PMC_COND_N_REAL_PARAMS)
    !> Integer parameters.
    integer, intent(in) :: integer_params(PMC_COND_N_INTEGER_PARAMS)

    real(kind=dp) :: D, kappa, V_dry

    D = diameter
    kappa = real_params(PMC_COND_KAPPA)
    V_dry = real_params(PMC_COND_V_DRY)

    water_activity_deriv = const%pi / 2d0 * D**2 * kappa * V_dry &
         / (const%pi / 6d0 * D**3 + (kappa - 1d0) * V_dry)**2

  end function water_activity_deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does one timestep (determined by this subroutine) of the
  !> condensation ODE.
  !!
  !! The timestep will not exceed max_dt, but might be less. If we in
  !! fact step all the way to max_dt then done will be true. This uses
  !! the explicit (forward) Euler integrator.
  subroutine condense_step_euler(max_dt, dt, done, env_state, aero_data, &
       aero_particle)
    
    !> Maximum timestep to integrate.
    real(kind=dp), intent(in) :: max_dt
    !> Actual timestep used.
    real(kind=dp), intent(out) :: dt
    !> Did we reach the maximum timestep?.
    logical, intent(out) :: done
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp) dvdt

    !write(*,*) 'do euler integration' 

    ! get timestep
    done = .false.

    call find_condense_timestep_variable(dt, env_state, aero_data, &
         aero_particle)
    !write(*,*) 'after find_condense_timestep_variable', 'dt=',dt,'max_dt=',max_dt
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if
    !write(*,*) 'usig ', ' dt=',dt

    ! do condensation
      
    call cond_growth_rate(dvdt, env_state, aero_data, aero_particle)
    aero_particle%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * dvdt

    ! ensure volumes stay positive
    aero_particle%vol(aero_data%i_water) = max(0d0, &
         aero_particle%vol(aero_data%i_water))
   
  end subroutine condense_step_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does one timestep (determined by this subroutine) of the
  !> condensation ODE.
  !!
  !! The timestep will not exceed max_dt, but might be less. If we in
  !! fact step all the way to max_dt then done will be true. This uses
  !! the explicit 4th-order Runge-Kutta integrator.
  subroutine condense_step_rk_fixed(max_dt, &
       dt, done, env_state, aero_data, aero_particle)
    
    !> Maximum timestep to integrate.
    real(kind=dp), intent(in) :: max_dt
    !> Actual timestep used.
    real(kind=dp), intent(out) :: dt
    !> Did we reach the maximum timestep?.
    logical, intent(out) :: done
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    ! get timestep
    done = .false.
    call find_condense_timestep_variable(dt, env_state, aero_data, &
         aero_particle)
    if (dt .ge. max_dt) then
       dt = max_dt
       done = .true.
    end if

    ! do condensation
    call condense_step_rk(dt, env_state, aero_data, aero_particle)
   
  end subroutine condense_step_rk_fixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Does one fixed timestep of Runge-Kutta-4.
  subroutine condense_step_rk(dt, env_state, aero_data, aero_particle)

    !> Timestep.
    real(kind=dp), intent(out) :: dt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    ! local variables
    real(kind=dp) k1, k2, k3, k4
    type(aero_particle_t) :: aero_particle_tmp

    call aero_particle_allocate_size(aero_particle_tmp, aero_data%n_spec)
    call aero_particle_copy(aero_particle, aero_particle_tmp)

    ! step 1
    call cond_growth_rate(k1, env_state, aero_data, aero_particle_tmp)

    ! step 2
    aero_particle_tmp%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * k1 / 2d0
    aero_particle_tmp%vol(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vol(aero_data%i_water))
    call cond_growth_rate(k2, env_state, aero_data, aero_particle_tmp)

    ! step 3
    aero_particle_tmp%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * k2 / 2d0
    aero_particle_tmp%vol(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vol(aero_data%i_water))
    call cond_growth_rate(k3, env_state, aero_data, aero_particle_tmp)

    ! step 4
    aero_particle_tmp%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) + dt * k3
    aero_particle_tmp%vol(aero_data%i_water) = &
         max(0d0, aero_particle_tmp%vol(aero_data%i_water))
    call cond_growth_rate(k4, env_state, aero_data, aero_particle_tmp)

    aero_particle%vol(aero_data%i_water) = &
         aero_particle%vol(aero_data%i_water) &
         + dt * (k1 / 6d0 + k2 / 3d0 + k3 / 3d0 + k4 / 6d0)
    aero_particle%vol(aero_data%i_water) = &
         max(0d0, aero_particle%vol(aero_data%i_water))

    call aero_particle_deallocate(aero_particle_tmp)
   
  end subroutine condense_step_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Just returns a constant timestep.
  subroutine find_condense_timestep_constant(dt, env_state, aero_data, &
       aero_particle)

    !> Timestep to use.
    real(kind=dp), intent(out) :: dt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

    dt = 5d-3

  end subroutine find_condense_timestep_constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes a timestep proportional to V / (dV/dt).
  subroutine find_condense_timestep_variable(dt, env_state, aero_data, &
       aero_particle)

    !> Timestep to use.
    real(kind=dp), intent(out) :: dt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

    !> Scale factor for timestep.
    real(kind=dp), parameter :: scale = 0.1d0

    real(kind=dp) pv, dvdt, di

    !> pv = total particle volume
    pv = aero_particle_volume(aero_particle)
    call cond_growth_rate(dvdt, env_state, aero_data, aero_particle)
    dt = abs(scale * pv / dvdt)
    !write(*,*) 'in find condense timestep variable pv =' , pv
!
    di=vol2diam(pv)
    !write(*,*) 'di =' , di
!
  end subroutine find_condense_timestep_variable
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the water volume growth rate due to condensation.
  subroutine cond_growth_rate(dvdt, env_state, aero_data, aero_particle)

    !> Dv/dt (m^3 s^{-1}).
    real(kind=dp), intent(out) :: dvdt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

!    call cond_growth_rate_one(dvdt, env_state, aero_data, aero_particle)
!    call cond_growth_rate_two(dvdt, env_state, aero_data, aero_particle)
!    call cond_growth_rate_three(dvdt, env_state, aero_data, aero_particle)
    call cond_growth_rate_four(dvdt, env_state, aero_data, aero_particle)

  end subroutine cond_growth_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the growth rate according to formula 15.73 (SP) using Newton method. 
  subroutine cond_growth_rate_one(dvdt, env_state, aero_data, aero_particle)

    !> Dv/dt (m^3 s^{-1}).
    real(kind=dp), intent(out) :: dvdt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

    !> Relative dm/dt convergence tol.
    real(kind=dp), parameter :: dmdt_rel_tol = 1d-8
    !> Function convergence tolerance.
    real(kind=dp), parameter :: f_tol = 1d-15
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100

    real(kind=dp) dmdt, pm, dmdt_tol

    !write(6,*) 'i am in scheme 1'

    pm = aero_particle_mass(aero_particle, aero_data)
    dmdt_tol = pm * dmdt_rel_tol

    dmdt = 0d0
    call cond_newt(dmdt, env_state, aero_data, cond_growth_rate_func, &
         dmdt_tol, f_tol, iter_max, aero_particle)
    
    dvdt = dmdt / aero_data%density(aero_data%i_water)

  end subroutine cond_growth_rate_one

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the growth rate according to formula 15.73 (SP) using Taylor series expansion (2nd order). 
  subroutine cond_growth_rate_two(dvdt, env_state, aero_data, aero_particle)

    !> Dv/dt (m^3 s^{-1}).
    real(kind=dp), intent(out) :: dvdt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle


    real(kind=dp) dmdt

    !write(6,*) 'i am in scheme 2'

    dmdt = 0d0
    call cond_taylor_quad(env_state, aero_data, dmdt, &
         aero_particle)
    
    dvdt = dmdt / aero_data%density(aero_data%i_water)

  end subroutine cond_growth_rate_two

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the growth rate according to formula 15.74 (SP). 
  subroutine cond_growth_rate_three(dvdt, env_state, aero_data, aero_particle)

    !> Dv/dt (m^3 s^{-1}).
    real(kind=dp), intent(out) :: dvdt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle


    real(kind=dp) dmdt

    !write(6,*) 'i am in scheme 3'

    dmdt = 0d0
    call cond_seinfeld(env_state, aero_data, dmdt, &
         aero_particle)
    
    dvdt = dmdt / aero_data%density(aero_data%i_water)

  end subroutine cond_growth_rate_three

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find the growth rate according to formula in terms of aw (eqn 74). 
  subroutine cond_growth_rate_four(dvdt, env_state, aero_data, aero_particle)

    !> Dv/dt (m^3 s^{-1}).
    real(kind=dp), intent(out) :: dvdt
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Relative dm/dt convergence tol.
    real(kind=dp), parameter :: pv_rel_tol = 1d-3
    !> Function convergence tolerance.
    real(kind=dp), parameter :: f_tol = 1d-15
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100
     
    real(kind=dp) :: dw, dwdt_init, dwdt, dw_tol
    real(kind=dp) :: A, k_a, k_ap_div, k_ap, pv, rho_water   
    integer :: k

!    write(*,*) 'i am in scheme 4 (aw)'

!   need to calculate diameter of the droplet dw

    pv = aero_particle_volume(aero_particle)
    dw = vol2diam(pv)

    dwdt_init = 0d0
    !write(*,*) 'id =', aero_particle%id, 'part. growth rate (before iteration) =', dwdt_init
    dwdt = dwdt_init
!    write(*,*) 'in cond growth rate four', 'pv =', pv , 'dw =', dw
  
    ! convert volume relative tolerance to diameter absolute tolerance
    dw_tol = vol2diam(pv * (1d0 + pv_rel_tol)) - vol2diam(pv)
    
    !write(*,*) 'dw_tol=',dw_tol
    call cond_newt(dwdt, env_state, aero_data, cond_growth_rate_aw, &
         dw_tol, f_tol, iter_max, aero_particle)
   
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 &
            * env_state%temp) ! (J m^{-1} s^{-1} K^{-1})
       k_ap_div = 1d0 + 2d0 &  ! dimensionless
            * k_a / (const%accom_coeff * dw * env_state_air_den(env_state) &
            * const%water_spec_heat) &
            * (2d0 * const%pi * const%air_molec_weight &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! (J m^{-1} s^{-1} K^{-1})

       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       A = const%water_latent_heat * rho_water &
            / (4d0 * k_ap * env_state%temp)
 
    dwdt= dwdt / A/ dw      ! dwdt = delta

    dvdt = dwdt * (dw **2d0) * (const%pi /2d0)
    !write(*,*) 'end of cond four, dvdt= ', dvdt
  end subroutine cond_growth_rate_four

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Scalar Newton's method for solving the implicit condensation
  !> functions.
  subroutine cond_newt(x, env_state, aero_data, func, x_tol, f_tol, &
       iter_max, aero_particle)

    !> Variable (set to init value on call).
    real(kind=dp), intent(inout) :: x
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> X convergence tolerance.
    real(kind=dp), intent(in) :: x_tol
    !> F convergence tolerance.
    real(kind=dp), intent(in) :: f_tol
    !> Maximum number of iterations.
    integer, intent(in) :: iter_max
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

#ifndef DOXYGEN_SKIP_DOC
    interface
       subroutine func(env_state, aero_data, init, x, f, df, aero_particle)
         use pmc_env_state
         use pmc_aero_data
         use pmc_aero_particle
         !> Environment state.
         type(env_state_t), intent(in) :: env_state
         !> Aerosol data.
         type(aero_data_t), intent(in) :: aero_data
         !> True if first Newton loop.
         logical, intent(in) :: init
         !> Independent variable to solve for.
         real(kind=dp), intent(in) :: x
         !> Function to solve.
         real(kind=dp), intent(out) :: f
         !> Derivative df/dx.
         real(kind=dp), intent(out) :: df
         !> Particle.
         type(aero_particle_t), intent(in) :: aero_particle
       end subroutine func
    end interface
#endif
    
    integer iter, k
    real(kind=dp) delta_f, delta_x, f, old_f, df
    real(kind=dp) :: check_x, check_f, check_df, finite_diff_df, rel_error

    call func(env_state, aero_data, .true., x, f, df, aero_particle)
    old_f = f

    iter = 0
    do
       !write(*,*) 'no. of iteration', iter
       iter = iter + 1
       delta_x = f / df
       !write(*,*) 'old x =' , x , 'old f=' ,old_f
       
       if (PMC_COND_CHECK_DERIVS) then
          check_x = x * (1d0 + PMC_COND_CHECK_EPS)
          call func(env_state, aero_data, .false., check_x, check_f, &
               check_df, aero_particle)
          finite_diff_df = (check_f - f) / (check_x - x)
          rel_error = abs(finite_diff_df - df) / (abs(finite_diff_df) + abs(df))
          if (rel_error > PMC_COND_CHECK_REL_TOL) then
             write(0,*) "ERROR: cond_newt: incorrect derivative"
             write(0,*) "x ", x
             write(0,*) "check_x ", check_x
             write(0,*) "delta_x ", (check_x - x)
             write(0,*) "f ", f
             write(0,*) "check_f ", check_f
             write(0,*) "delta_f ", (check_f - f)
             write(0,*) "df ", df
             write(0,*) "finite_diff_df ", finite_diff_df
             write(0,*) "rel_error ", rel_error
          end if
       end if

       x = x - delta_x
       call func(env_state, aero_data, .false., x, f, df, aero_particle)
       delta_f = f - old_f
       old_f = f

       !write(*,*) 'new x =' , x , 'new f=' , f
       !write(*,*) 'delta x =', delta_x , 'delta_f=', delta_f
       if (iter .ge. iter_max) then
          call die_msg(449343787, 'Newton iteration failed to terminate')
       end if
       
       ! FIXME: gfortran 4.1.1 requires the "then" in the following
       ! statement, rather than using a single-line "if" statement.
       if ((abs(delta_x) .lt. x_tol) &
            .and. (abs(delta_f) .lt. f_tol)) then
          exit
       end if
    end do
 
    !write(*,*) 'id =', aero_particle%id, 'finish using', iter ,'iteration(s)'
  end subroutine cond_newt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit growth rate function.
  subroutine cond_growth_rate_func(env_state, aero_data, init, dmdt, &
       f, df, aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Mass growth rate dm/dt (kg s^{-1}).
    real(kind=dp), intent(in) :: dmdt
    !> Error.
    real(kind=dp), intent(out) :: f
    !> Derivative of error with respect to x.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    
    ! local variables
    real(kind=dp), save :: k_a, k_ap, k_ap_div, D_v, D_v_div, D_vp, d_p, pv
    real(kind=dp), save :: rat, fact1, fact2, c1, c2, c3, c4, c5
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute

    real(kind=dp) T_a ! droplet temperature (K), determined as part of solve

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)
       d_p = vol2diam(pv) ! (m)
       
       ! molecular diffusion coefficient uncorrected
       D_v = 0.211d-4 / (env_state%pressure / const%air_std_press) &
            * (env_state%temp / 273d0)**1.94d0 ! (m^2 s^{-1})
       ! molecular diffusion coefficient corrected for non-continuum effects
       D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (const%accom_coeff * d_p)) &
            * (2d0 * const%pi * M_water &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       D_vp = D_v / D_v_div
       
       ! TEST: use the basic expression for D_vp
       ! D_vp = D_v                ! (m^2 s^{-1})
       ! FIXME: check whether we can reinstate the correction
       
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 &
            * env_state%temp) ! (J m^{-1} s^{-1} K^{-1})
       k_ap_div = 1d0 + 2d0 &  ! dimensionless
            * k_a / (const%accom_coeff * d_p * env_state_air_den(env_state) &
            * const%water_spec_heat) &
            * (2d0 * const%pi * const%air_molec_weight &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! (J m^{-1} s^{-1} K^{-1})

       rat = env_state_sat_vapor_pressure(env_state) &
            / (const%univ_gas_const * env_state%temp)
       fact1 = const%water_latent_heat * M_water &
            / (const%univ_gas_const * env_state%temp)
       fact2 = const%water_latent_heat &
            / (2d0 * const%pi * d_p * k_ap * env_state%temp)

       c1 = 2d0 * const%pi * d_p * D_vp * M_water * rat
       c2 = 4d0 * M_water &
            * const%water_surf_eng / (const%univ_gas_const * rho_water * d_p)
       c3 = c1 * fact1 * fact2
       c4 = const%water_latent_heat / (2d0 * const%pi * d_p * k_ap)
       ! incorrect expression from Majeed and Wexler:
       !     c5 = nu * eps * M_water * rho_solute * r_n**3d0 &
       !         / (M_solute * rho_water * ((d_p / 2)**3d0 - r_n**3))
       ! c5 = nu * eps * M_water / M_solute * g_solute / g_water
       ! corrected according to Jim's note:
       !    c5 = nu * eps * M_water / M_solute * g_solute / &
       !         (g_water + (rho_water / rho_solute) * eps * g_solute)
           c5 = M_water / M_solute * nu * g_solute / &
                (g_water - rho_water / rho_solute * (1d0 - eps) * g_solute ) 

    end if

    T_a = env_state%temp + c4 * dmdt ! (K)
    
    f = dmdt - c1 * (env_state%rel_humid - exp(c2 / T_a - c5)) &
         / (1d0 + c3 * exp(c2 / T_a - c5))
    
    df = 1d0 + c1 * env_state%rel_humid &
         * (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * &
         exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 + c1 * &
         (exp(c2 / T_a - c5) * (-1d0) * c2 * c4 / T_a**2d0 * (1d0 + c3 &
         * exp(c2 / T_a -c5))**(-1d0) + exp(c2 / T_a - c5) * (-1d0) * &
         (1d0 + c3 * exp(c2 / T_a -c5))**(-2d0) * c3 * exp(c2 / T_a - &
         c5) * (-1d0) * c2 * c4 / T_a**2d0)

    ! NOTE: we could return T_a (the droplet temperature) if we have
    ! any need for it.

  end subroutine cond_growth_rate_func
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cond_taylor_quad(env_state, aero_data, dmdt, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Mass growth rate dm/dt (kg s^{-1}).
    real(kind=dp), intent(out) :: dmdt
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    
    ! local variables
    real(kind=dp), save :: k_a, k_ap, k_ap_div, D_v, D_v_div, D_vp, d_p, pv
    real(kind=dp), save :: rat, fact1, fact2, fact3, c0, c1, c2, c3, c4, c5
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute

    real(kind=dp) T_a ! droplet temperature (K), determined as part of solve

   
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)
       d_p = vol2diam(pv) ! (m)
       
       ! molecular diffusion coefficient uncorrected
       D_v = 0.211d-4 / (env_state%pressure / const%air_std_press) &
            * (env_state%temp / 273d0)**1.94d0 ! (m^2 s^{-1})
       ! molecular diffusion coefficient corrected for non-continuum effects
       D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (const%accom_coeff * d_p)) &
            * (2d0 * const%pi * M_water &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       D_vp = D_v / D_v_div
       
       ! TEST: use the basic expression for D_vp
       ! D_vp = D_v                ! (m^2 s^{-1})
       ! FIXME: check whether we can reinstate the correction
       
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 &
            * env_state%temp) ! (J m^{-1} s^{-1} K^{-1})
       k_ap_div = 1d0 + 2d0 &  ! dimensionless
            * k_a / (const%accom_coeff * d_p * env_state_air_den(env_state) &
            * const%water_spec_heat) &
            * (2d0 * const%pi * const%air_molec_weight &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! (J m^{-1} s^{-1} K^{-1})

       rat = env_state_sat_vapor_pressure(env_state) &
            / (const%univ_gas_const * env_state%temp)

       c0 = const%water_latent_heat &
            / (2d0 * const%pi * d_p * k_ap)

       c1 = const%water_latent_heat * rho_water &
            / (4d0 * k_ap * env_state%temp)

       c2 = 4d0 * D_vp * M_water * rat &
            / rho_water

       c3 = const%water_latent_heat * M_water &
            / (const%univ_gas_const * env_state%temp)

       c4 = 4d0 * M_water * const%water_surf_eng &
            / (const%univ_gas_const * rho_water * d_p * env_state%temp)

       ! incorrect expression from Majeed and Wexler:
       !     c5 = nu * eps * M_water * rho_solute * r_n**3d0 &
       !         / (M_solute * rho_water * ((d_p / 2)**3d0 - r_n**3))
       ! c5 = nu * eps * M_water / M_solute * g_solute / g_water
       ! corrected according to Jim's note:
       !    c5 = nu * eps * M_water / M_solute * g_solute / &
       !         (g_water + (rho_water / rho_solute) * eps * g_solute)
           c5 = M_water / M_solute * nu * g_solute / &
                (g_water - rho_water / rho_solute * (1d0 - eps) * g_solute ) 


    T_a = env_state%temp + c0 * dmdt ! (K)

    ! forming the quadratic equation coefficients
    
      fact1 = -1d0/2d0 * &
              c2*(exp(c4 - c5)*(-2d0+4d0*(c3 - c4) - (c3 - c4)*(c3 - c4))) &
              * (c1 * 2d0 / (const%pi * rho_water * d_p))**2d0 
 
      fact2 = (1d0/c1 - c2 * exp(c4 - c5)*(1d0 - c3 + c4)) &
              * (c1 * 2d0 / (const%pi * rho_water * d_p)) 
 
      fact3 = -1d0*c2*(env_state%rel_humid - exp(c4 - c5))  
 
   ! solve for dmdt using analytic quadratic equation and taking the larger root

      dmdt = (-1d0*fact2 + sqrt (fact2*fact2 - 4d0*fact1*fact3)) &
              /(2d0*fact1)


    ! NOTE: we could return T_a (the droplet temperature) if we have
    ! any need for it.

  end subroutine cond_taylor_quad
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cond_seinfeld(env_state, aero_data, dmdt, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Mass growth rate dm/dt (kg s^{-1}).
    real(kind=dp), intent(out) :: dmdt
    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    
    ! local variables
    real(kind=dp), save :: k_a, k_ap, k_ap_div, D_v, D_v_div, D_vp, d_p, pv
    real(kind=dp), save :: rat, c0, c1, c2, c3, c4, c5
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute

    real(kind=dp) T_a ! droplet temperature (K), determined as part of solve


       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)
       d_p = vol2diam(pv) ! (m)
       
       ! molecular diffusion coefficient uncorrected
       D_v = 0.211d-4 / (env_state%pressure / const%air_std_press) &
            * (env_state%temp / 273d0)**1.94d0 ! (m^2 s^{-1})
       ! molecular diffusion coefficient corrected for non-continuum effects
       D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (const%accom_coeff * d_p)) &
            * (2d0 * const%pi * M_water &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       D_vp = D_v / D_v_div

       ! TEST: use the basic expression for D_vp
       ! D_vp = D_v                ! (m^2 s^{-1})
       ! FIXME: check whether we can reinstate the correction
       
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 &
            * env_state%temp) ! (J m^{-1} s^{-1} K^{-1})
       k_ap_div = 1d0 + 2d0 &  ! dimensionless
            * k_a / (const%accom_coeff * d_p * env_state_air_den(env_state) &
            * const%water_spec_heat) &
            * (2d0 * const%pi * const%air_molec_weight &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! (J m^{-1} s^{-1} K^{-1})

       rat = env_state_sat_vapor_pressure(env_state) &
            / (const%univ_gas_const * env_state%temp)
       c0 = const%water_latent_heat &
            / (2d0 * const%pi * d_p * k_ap)

       c1 = rho_water / (4d0 * D_vp * M_water * rat)

       c2 = 4d0 * M_water * const%water_surf_eng & 
            / (const%univ_gas_const * rho_water * d_p * env_state%temp)

       c3 = const%water_latent_heat * M_water &
            / (const%univ_gas_const * env_state%temp)

       c4 = const%water_latent_heat * rho_water / (4d0 * k_ap * env_state%temp) 
       
       ! incorrect expression from Majeed and Wexler:
       !     c5 = nu * eps * M_water * rho_solute * r_n**3d0 &
       !         / (M_solute * rho_water * ((d_p / 2)**3d0 - r_n**3))
       ! c5 = nu * eps * M_water / M_solute * g_solute / g_water
       ! corrected according to Jim's note:
       !    c5 = nu * eps * M_water / M_solute * g_solute / &
       !         (g_water + (rho_water / rho_solute) * eps * g_solute)
           c5 = M_water / M_solute * nu * g_solute / &
                (g_water - rho_water / rho_solute * (1d0 - eps) * g_solute ) 


    T_a = env_state%temp + c0 * dmdt ! (K)
    
    dmdt = (env_state%rel_humid - exp (c2 - c5)) * &
            (1d0 / 2d0) * const%pi * rho_water * d_p / &
             (c1 + c4 * (c3 - 1d0))  
    
    ! NOTE: we could return T_a (the droplet temperature) if we have
    ! any need for it.

  end subroutine cond_seinfeld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve the particle growth equation in terms of aw using newton solver
  !> This subroutine calcualte the growth equation and then input to the 
  !> solver cond_newt
  subroutine cond_growth_rate_aw(env_state, aero_data, init, d_diamdt, & 
       f, df, aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol particle.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Mass growth rate dm/dt (kg s^{-1}).
    real(kind=dp), intent(in) :: d_diamdt
    !> Error.
    real(kind=dp), intent(out) :: f
    !> Derivative of error with respect to x.
    real(kind=dp), intent(out) :: df

    ! local variables
    real(kind=dp), save :: k_a, k_ap, k_ap_div, D_v, D_v_div, D_vp, d_p, pv
    real(kind=dp), save :: A, B, C, delta, dw 
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute, df_ddelta

    real(kind=dp) T_a ! droplet temperature (K), determined as part of solve
    real(kind=dp) sat ! environmental saturation 

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)
       d_p = vol2diam(pv) ! (m)
       
       ! molecular diffusion coefficient uncorrected
       D_v = 0.211d-4 / (env_state%pressure / const%air_std_press) &
            * (env_state%temp / 273d0)**1.94d0 ! (m^2 s^{-1})
       ! molecular diffusion coefficient corrected for non-continuum effects
       D_v_div = 1d0 + (2d0 * D_v * 1d-4 / (const%accom_coeff * d_p)) &
            * (2d0 * const%pi * M_water &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       D_vp = D_v / D_v_div

       ! TEST: use the basic expression for D_vp
       ! D_vp = D_v                ! (m^2 s^{-1})
       ! FIXME: check whether we can reinstate the correction
       
       ! thermal conductivity uncorrected
       k_a = 1d-3 * (4.39d0 + 0.071d0 &
            * env_state%temp) ! (J m^{-1} s^{-1} K^{-1})
       k_ap_div = 1d0 + 2d0 &  ! dimensionless
            * k_a / (const%accom_coeff * d_p * env_state_air_den(env_state) &
            * const%water_spec_heat) &
            * (2d0 * const%pi * const%air_molec_weight &
            / (const%univ_gas_const * env_state%temp))**0.5d0
       ! thermal conductivity corrected
       k_ap = k_a / k_ap_div     ! (J m^{-1} s^{-1} K^{-1})

       sat=env_state%rel_humid

       A = const%water_latent_heat * rho_water &
            / (4d0 * k_a * env_state%temp)

       B = 4d0 * D_v * M_water * env_state_sat_vapor_pressure(env_state) &
            / (rho_water * const%univ_gas_const * env_state%temp)

       C = const%water_latent_heat * M_water &
            / (const%univ_gas_const * env_state%temp)

    end if

!   get the wet diameter to calculate the diameter growth rate

!       write(*,*) 'we are in cond_growth_rate_aw' 
       pv = aero_particle_volume(aero_particle)
       dw = vol2diam(pv)

!       write(*,*) 'in cond growth rate aw' ,'pv =', pv , 'dw =', dw
 
       sat=env_state%rel_humid
       delta= A * dw * d_diamdt

!       write(*,*) 'const%water_latent_heat = ', const%water_latent_heat
!       write(*,*) 'density of water = ', rho_water
!       write(*,*) 'k_a =', k_a
!       write(*,*) 'k_ap =', k_ap
!       write(*,*) 'k_ap_div =', k_ap_div
!       write(*,*) 'M_water =', M_water 
!       write(*,*) 'temp =' , env_state%temp
!       write(*,*) 'sat=', sat
!       write(*,*) 'A =', A
!       write(*,*) 'latent heat=', const%water_latent_heat
!       write(*,*) 'density of water=', rho_water
!       write(*,*) 'heat diff=', k_a
!       write(*,*) 'temp=', env_state%temp
!       write(*,*) 'B =', B
!       write(*,*) 'D_vp =', D_vp
!       write(*,*) 'D_v =', D_v
!       write(*,*) 'D_v_div =', D_v_div
!       write(*,*) 'sat pressure=', env_state_sat_vapor_pressure(env_state)
!       write(*,*) 'C =', C
!       write(*,*) 'delta= ', delta
!       write(*,*) 'd_diamdt =',  d_diamdt

       T_a = env_state%temp + (const%water_latent_heat * rho_water / (4d0 * k_ap)) & 
              * dw * d_diamdt  ! (K)
!       write(*,*) 'aw ', aw(aero_data, aero_particle, dw)
!       write(*,*)'kelvin_argu ', kelvin_argu(aero_particle, aero_data, env_state, dw) 
       f = delta / A - B * (sat - (1d0 / (1d0 + delta)) & 
            * aw(aero_data, aero_particle, dw) * & 
            exp(kelvin_argu(aero_particle, aero_data, env_state, dw) &
            * (1d0 / (1d0+delta)) + & 
            C * (delta / (1d0+delta))))

       df_ddelta = 1d0 / A - B * (aw(aero_data, aero_particle, dw) &
            * (1d0 / ((1d0 + delta)**2d0))   &
            * exp( C * (delta/(1d0+delta)) + kelvin_argu(aero_particle, aero_data, env_state, dw) &
            *(1d0 / (1d0 + delta))) * (1d0 - (1d0/(1d0 + delta))*(C & 
            - kelvin_argu (aero_particle, aero_data, env_state, dw))))  
       df = df_ddelta * A * dw
       
!       write(*,*) 'f= ', f
!       write(*,*) 'df= ', df 

  end subroutine cond_growth_rate_aw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add water to the particle until it is in equilibrium.
  subroutine equilibriate_particle(env_state, aero_data, aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    ! parameters
    !> Maximum iterations.
    integer, parameter :: it_max = 400
    !> Volume relative convergence tolerance.
!    real(kind=dp), parameter :: pv_rel_tol = 1d-6
    real(kind=dp), parameter :: pv_rel_tol = 1d-3
    !> Function convergence tolerance.
    real(kind=dp), parameter :: f_tol = 1d-15
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100
    !> Initial value.
    !> want to get variable initial guess as the dry diameter
!    real(kind=dp), parameter :: dw_init = 1d-6
    
    real(kind=dp) dw ! wet diameter of particle
    real(kind=dp) dw_tol, pv, di
    real(kind=dp) dw_init     

    !> total volume of the particle before equilibriation
    pv = aero_particle_volume(aero_particle)
    di = vol2diam(pv)
    !write(*,*) 'id =', aero_particle%id, 'wet diam (before equilibriation)', di  
    !> dry diameter of the particle
    dw_init = vol2diam(aero_particle_solute_volume(aero_particle, aero_data))
    !write(*,*) 'id =', aero_particle%id, 'dry diam (before iteration) =', dw_init
    dw = dw_init

    ! convert volume relative tolerance to diameter absolute tolerance
    dw_tol = vol2diam(pv * (1d0 + pv_rel_tol)) - vol2diam(pv)
    !write(*,*) 'sat (in equilibriate_particle) BF newton =' , env_state%rel_humid
!    write(*,*) 'dw_tol = ', dw_tol
!    write(*,*) 'id =', aero_particle%id 

    call cond_newt(dw, env_state, aero_data, equilibriate_func, &
         dw_tol, f_tol, iter_max, aero_particle)

    !write(*,*) 'sat (in equilibriate_particle) AF newton =' , env_state%rel_humid
!    aero_particle%vol(aero_data%i_water) = diam2vol(dw) - pv
    !> the amount of water after equilibriation
    !> dw is the wet diameter after equilibriation 
    !> dw_init is the dry diameter before equilibriation
    aero_particle%vol(aero_data%i_water) = diam2vol(dw) - diam2vol(dw_init)
    di=vol2diam(pv)

    !write(*,*)  'id = ', aero_particle%id, 'wet diam (after iteration) = ', dw
  end subroutine equilibriate_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

!    call equilibriate_func_old(env_state, aero_data, init, dw, f, df, &
!       aero_particle)

!    call equilibriate_func_new(env_state, aero_data, init, dw, f, df, &
!       aero_particle)

!     write(*,*) 'sat (in equilibriate_func) 1 =', env_state%rel_humid
    call equilibriate_func_aw(env_state, aero_data, init, dw, f, df, &
       aero_particle)
    
!     write(*,*) 'sat (in equiliraiate_func) 2 =', env_state%rel_humid
!    write(*,*) 'dw = ', dw
!    write(*,*) 'f = ', f
!    write(*,*) 'df = ', df
!    write(*,*) 'id = ', aero_particle%id

    end subroutine equilibriate_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func_old(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp), save :: c0, c1, c3, c4, dc0, dc2, dc3
    real(kind=dp), save :: A, B
    real(kind=dp), save ::  pv
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)
    
       A = 4d0 * M_water * const%water_surf_eng &
            / (const%univ_gas_const * env_state%temp * rho_water)
       
       B = nu * eps * M_water * rho_solute &
            * vol2rad(pv)**3.d0 / (M_solute * rho_water)
       
       c4 = log(env_state%rel_humid) / 8d0
       c3 = A / 8d0
       
       dc3 = log(env_state%rel_humid) / 2d0
       dc2 = 3d0 * A / 8d0
       
       c1 = B - log(env_state%rel_humid) * vol2rad(pv)**3d0
       c0 = A * vol2rad(pv)**3d0
       dc0 = c1

    endif
    
    f = c4 * dw**4d0 - c3 * dw**3d0 + c1 * dw + c0
    df = dc3 * dw**3d0 -dc2 * dw**2d0 + dc0

  end subroutine equilibriate_func_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func_new(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp), save :: c1, c2, c3, c4, c5, c6
    real(kind=dp), save :: A, B
    real(kind=dp), save ::  pv
    real(kind=dp), save :: M_water, M_solute, rho_water, rho_solute
    real(kind=dp), save :: eps, nu, g_water, g_solute

    if (init) then
       ! Start of new Newton loop, compute all constants

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       M_solute = aero_particle_solute_molec_weight(aero_particle, aero_data)
       ! (kg/mole)
       nu = aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)
       eps = aero_particle_solute_solubility(aero_particle, aero_data) ! (1)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)
       rho_solute = aero_particle_solute_density(aero_particle, aero_data)
       ! (kg/m^3)
       g_water = aero_particle_water_mass(aero_particle, aero_data) ! (kg)
       g_solute = aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

       pv = aero_particle_volume(aero_particle) ! (m^3)

       c1=log(env_state%rel_humid)
       c2 = 4d0 * M_water &
            * const%water_surf_eng / (env_state%temp * const%univ_gas_const * rho_water)
       c3=(6d0*nu*M_water) / (const%pi*rho_water*M_solute) * g_solute 
       
!       c4=1d-9                   ! hardcode insoluble diameter
        c4=0d0
       c5=c2*(c4**3d0)
       c6=c1*(c4**3d0) 
      
      endif
    
    f = c1*dw**4d0-c2*dw**3d0+(c3-c6)*dw+c5     
    df= 4d0*c1*dw**3d0 - 3d0*c2*dw**2d0 + (c3-c6) 

    ! NOTE: we could return T_a (the droplet temperature) if we have
    ! any need for it.


  end subroutine equilibriate_func_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the error function value and its derivative for the
  !> implicit function that determines the equilibrium state of a
  !> particle.
  subroutine equilibriate_func_aw(env_state, aero_data, init, dw, f, df, &
       aero_particle)

    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> True if first Newton loop.
    logical, intent(in) :: init
    !> Wet diameter (m).
    real(kind=dp), intent(in) :: dw
    !> Function value.
    real(kind=dp), intent(out) :: f
    !> Function derivative df/dx.
    real(kind=dp), intent(out) :: df
    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle

    real(kind=dp) :: sat, tot_kappa, kappa
    
    integer :: k
 
!       write(*,*) 'sat (in equilibriate_func_aw) 0 = ', env_state%rel_humid
!       write(*,*) 'init =', init
    if (init) then
          
       ! Start of new Newton loop, compute all constants

       call compute_kappa(aero_data, aero_particle, tot_kappa)   
       kappa=tot_kappa
!       write(*,*) ' id = ', aero_particle%id, 'kappa = ' , kappa       

       sat = env_state%rel_humid
!       write(*,*) 'sat (in equilibriate_func_aw) 1 = ', sat
   
    endif
         
!       write(*,*) 'now start to calculate f and df in equilibration'
!       write(*,*) ' id = ', aero_particle%id, 'kappa = ', kappa

       sat = env_state%rel_humid
       
        f= sat - aw(aero_data, aero_particle, dw) & 
              * exp(kelvin_argu(aero_particle, aero_data, env_state, dw))
       df= -  aw_deri(aero_data, aero_particle, dw) & 
              * exp(kelvin_argu(aero_particle, aero_data, env_state, dw)) &
              + kelvin_deri(aero_particle, aero_data, env_state, dw) &
              * aw(aero_data, aero_particle, dw)

!       write(*,*) 'f = ', f 
!       write(*,*) 'df = ', df
!       write(*,*)  'wet diam (after iteration) = ', dw
!       write(*,*) 'sat (in equilibriate_func_aw) 2 = ', sat
    
   end subroutine equilibriate_func_aw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the argument of the Kelvin term.
  real(kind=dp) function kelvin_argu(aero_particle, aero_data, env_state, dw)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol Particle.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment.
    type(env_state_t), intent(in) :: env_state
    !> Wet diameter
    real(kind=dp), intent(in) :: dw

    real(kind=dp), save :: M_water, rho_water

!      write(*,*) 'now in function kelvin_term'
!      write(*,*) 'input dw to the funtion kelvin_argu = ', dw

       M_water = aero_particle_water_molec_weight(aero_data) ! (kg/mole)
       rho_water = aero_particle_water_density(aero_data) ! (kg/m^3)

    kelvin_argu= 4d0 * M_water * const%water_surf_eng &
         / (env_state%temp * const%univ_gas_const * rho_water * dw) 
    
!    write(*,*) 'input dw to the funtion aw = ', dw
!    write(*,*) 'mole weight of h20 =', M_water
!    write(*,*) 'temp =', env_state%temp
!    write(*,*) 'gas const =', const%univ_gas_const
!    write(*,*) 'rho water =', rho_water 
!    write(*,*) 'surface tension =', const%water_surf_eng
!    write(*,*) aero_particle%id, 'kelvin argu= ', kelvin_argu

  end function kelvin_argu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the derivative of the Kelvin term.
  real(kind=dp) function kelvin_deri(aero_particle, aero_data, env_state, dw)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol Particle.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment.
    type(env_state_t), intent(in) :: env_state
    !> Wet diameter.
    real(kind=dp), intent(in) :: dw

!      write(*,*) 'now in function kelvin_deri'      
!      write(*,*) 'input dw to the funtion aw = ', dw

    kelvin_deri=exp(kelvin_argu(aero_particle, aero_Data, env_state, dw)) & 
             *kelvin_argu(aero_particle, aero_data, env_state, dw) / dw

!    write(*,*) aero_particle%id, 'kelvin deri= ', kelvin_deri
  end function kelvin_deri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the water activity.
  real(kind=dp) function aw(aero_data, aero_particle, dw)

    !> Particle.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol data
    type(aero_particle_t), intent(in) :: aero_particle

    real(kind=dp), intent(in) :: dw
    
    real(kind=dp) :: vs, vt, kappa, tot_kappa

    integer :: k

!       write(*,*) 'now in function aw' 

       vs=0d0
       vt=0d0
!   vs is the dry particle volume (non-water components)
       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
          vs=vs+aero_particle%vol(k)
          endif
          vt=vt+aero_particle%vol(k)
       enddo
    
       call compute_kappa(aero_data, aero_particle, tot_kappa)   
       kappa=tot_kappa

!      write(*,*) aero_particle%id, 'kappa = ', kappa
!      write(*,*) 'input dw to the funtion aw = ', dw
!      write(*,*) 'dry salt volume in the function aw = ', vs
!      write(*,*) 'wet salt volume (dw^3) ', const%pi/6d0*(dw**3d0)     
       
      aw=(const%pi/6d0*(dw**3d0)-vs)/(const%pi/6d0*(dw**3d0)+(kappa-1d0)*vs)
  
!      write(*,*) 'kappa =', kappa   
!      write(*,*) aero_particle%id, 'aw= ', aw
  end function aw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> calculate the derivative of the water activity.
  real(kind=dp) function aw_deri(aero_data, aero_particle, dw)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data
    type(aero_data_t), intent(in) :: aero_data

    real(kind=dp), intent(in) :: dw
    
    real(kind=dp) :: vs, vt, kappa, tot_kappa
    
    integer :: k

!       write(*,*) 'now in function aw_deri'

       vs=0d0
       vt=0d0
!   vs is the dry particle volume (non-water components)
       do k=1,aero_data%n_spec 
          if (k /= aero_data%i_water) then
          vs=vs+aero_particle%vol(k)
          endif
          vt=vt+aero_particle%vol(k)
       enddo
    
       call compute_kappa(aero_data, aero_particle, tot_kappa)   
       kappa=tot_kappa
      
!      write(*,*) 'input dw to the funtion aw = ', dw
!      write(*,*) 'dry salt volume in the function aw = ', vs

       aw_deri= (kappa*vs*const%pi/2d0*(dw**2d0))   &
       /((const%pi/6d0*(dw**3d0)+(kappa-1d0)*vs)**2d0)

!    write(*,*) aero_particle%id, 'aw_deri= ', aw_deri
  end function aw_deri

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Call equilibriate_particle() on each particle in the aerosol.
  subroutine aero_state_equilibriate(bin_grid, env_state, aero_data, &
       aero_state)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
 
    integer :: i_bin, i
    
    do i_bin = 1,bin_grid%n_bin
       do i = 1,aero_state%bin(i_bin)%n_part
         !write(*,*) 'hello id =', aero_particle%id
         env_state%rel_humid=9.50d-1
         !write(*,*) 'RH (in aero_state_equi) =', env_state%rel_humid
          call equilibriate_particle(env_state, aero_data, &
               aero_state%bin(i_bin)%particle(i))
!     temporary stop statement
          !write(*,*) 'stop at particle i', i
!        STOP
       end do
    end do

  end subroutine aero_state_equilibriate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_condensation
