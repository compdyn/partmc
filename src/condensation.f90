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
  use pmc_env_data
  use pmc_env_state
  use pmc_aero_data
  use pmc_util
  use pmc_aero_particle
  use pmc_constants
  use dvode_f90_m
  
  logical, parameter :: PMC_COND_CHECK_DERIVS = .false.
  real(kind=dp), parameter :: PMC_COND_CHECK_EPS = 1d-8
  real(kind=dp), parameter :: PMC_COND_CHECK_REL_TOL = 5d-2

  integer, parameter :: PMC_COND_N_REAL_PARAMS = 16
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
  integer, parameter :: PMC_COND_V_COMP = 12
  integer, parameter :: PMC_COND_TIME_0 = 13
  integer, parameter :: PMC_COND_TIME_1 = 14
  integer, parameter :: PMC_COND_TEMP_0 = 15
  integer, parameter :: PMC_COND_TEMP_1 = 16

  !> Diameter at which the pmc_cond_saved_delta_star was evaluated.
  real(kind=dp), save :: pmc_cond_saved_diameter
  !> Saved value of delta_star evaluated at pmc_cond_saved_diameter.
  real(kind=dp), save :: pmc_cond_saved_delta_star

  real(kind=dp), save :: save_real_params(PMC_COND_N_REAL_PARAMS)

  logical, save :: rh_first = .true.
  integer, save :: d_offset = 1 ! 1 if rh_first, otherwise 0

  !>DEBUG
  !logical, save :: kill_flag = .false.
  !integer, save :: i_call_f
  !integer, pointer, save :: n_call_f(:)
  !<DEBUG

  type(aero_data_t) :: condense_saved_aero_data
  type(env_data_t) :: condense_saved_env_data
  type(env_state_t) :: condense_saved_env_state
  real(kind=dp) :: condense_saved_V_comp
  real(kind=dp), pointer, save :: cond_kappa(:) ! FIXME: rename to saved_kappa, make allocatable, no save
  real(kind=dp), pointer, save :: cond_D_dry(:) ! FIXME: rename to saved_D_dry, make allocatable, no save

  type condense_params_t
     real(kind=dp) :: T
     real(kind=dp) :: H
     real(kind=dp) :: p
     real(kind=dp) :: V_comp
     real(kind=dp) :: rho_w
     real(kind=dp) :: M_w
     real(kind=dp) :: P_0
     real(kind=dp) :: rho_air
     real(kind=dp) :: k_a
     real(kind=dp) :: D_v
     real(kind=dp) :: U
     real(kind=dp) :: V
     real(kind=dp) :: W
     real(kind=dp) :: X
     real(kind=dp) :: Y
     real(kind=dp) :: Z
     real(kind=dp) :: D
     real(kind=dp) :: kappa
     real(kind=dp) :: D_dry
     real(kind=dp) :: k_ap
     real(kind=dp) :: dkap_dD
     real(kind=dp) :: D_vp
     real(kind=dp) :: dDvp_dD
     real(kind=dp) :: a_w
     real(kind=dp) :: daw_dD
     real(kind=dp) :: delta_star
     real(kind=dp) :: Ddot
     real(kind=dp) :: Hdot
     real(kind=dp) :: dh_ddelta
     real(kind=dp) :: dh_dD
     real(kind=dp) :: dh_dH
     real(kind=dp) :: ddeltastar_dD
     real(kind=dp) :: ddeltastar_dH
     real(kind=dp) :: dDdot_dD
     real(kind=dp) :: dDdot_dH
     real(kind=dp) :: dHdot_dD
     real(kind=dp) :: dHdot_dH
  end type condense_params_t

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> vapor.
  subroutine condense_particles(bin_grid, env_state, env_data, aero_data, &
       aero_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t

    call condense_particles_new_unified(bin_grid, env_state, env_data, &
         aero_data, aero_state, del_t)

  end subroutine condense_particles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> vapor.
  subroutine condense_particles_new_unified(bin_grid, env_state, env_data, &
       aero_data, aero_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    
    integer :: i_bin, j, new_bin, k, i_part
    real(kind=dp) :: pre_water, post_water
    type(aero_particle_t), pointer :: aero_particle
    integer :: n_eqn, itask, istate
    real(kind=dp) :: state(aero_state%n_part + 1), init_time, final_time
    real(kind=dp) :: abs_tol_vector(aero_state%n_part + 1)
    type(vode_opts) :: options

    pre_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          pre_water = pre_water + aero_particle%vol(aero_data%i_water)
       end do
    end do

    call aero_data_allocate(condense_saved_aero_data)
    call env_data_allocate(condense_saved_env_data)
    call env_state_allocate(condense_saved_env_state)

    call aero_data_copy(aero_data, condense_saved_aero_data)
    call env_data_copy(env_data, condense_saved_env_data)
    call env_state_copy(env_state, condense_saved_env_state)

    allocate(cond_kappa(aero_state%n_part))
    allocate(cond_D_dry(aero_state%n_part))
    i_part = 0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          i_part = i_part + 1
          aero_particle => aero_state%bin(i_bin)%particle(j)
          cond_kappa(i_part) = aero_particle_solute_kappa(aero_particle, aero_data)
          cond_D_dry(i_part) = vol2diam(aero_particle_solute_volume(aero_particle, aero_data))
          state(i_part + d_offset) = aero_particle_diameter(aero_particle)
          abs_tol_vector(i_part + d_offset) = max(1d-12, &
               1d-8 * (state(i_part + d_offset) - cond_D_dry(i_part)))
       end do
    end do
    if (rh_first) then
       state(1) = env_state%rel_humid
       abs_tol_vector(1) = 1d-10
    else
       state(aero_state%n_part + 1) = env_state%rel_humid
       abs_tol_vector(aero_state%n_part + 1) = 1d-10
    end if

    ! set VODE inputs
    n_eqn = aero_state%n_part + 1
    init_time = 0d0
    final_time = del_t
    itask = 1 ! just output val at final_time
    istate = 1 ! first call for this ODE

    options = set_opts(sparse_j = .true., &
         user_supplied_jacobian = .true., &
         abserr_vector = abs_tol_vector, relerr = 1d-8)
    call dvode_f90(condense_vode_new_unified_f, n_eqn, state, init_time, final_time, &
         itask, istate, options, j_fcn = condense_vode_new_unified_jac)

    if (rh_first) then
       env_state%rel_humid = state(1)
    else
       env_state%rel_humid = state(aero_state%n_part + 1)
    end if

    i_part = 0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          i_part = i_part + 1
          aero_particle => aero_state%bin(i_bin)%particle(j)

          ! translate output back to particle
          aero_particle%vol(aero_data%i_water) = diam2vol(state(i_part + d_offset)) &
               - aero_particle_solute_volume(aero_particle, aero_data)

          ! ensure volumes stay positive
          aero_particle%vol(aero_data%i_water) = max(0d0, &
               aero_particle%vol(aero_data%i_water))
       end do
    end do

    ! We resort the particles in the bins after only all particles
    ! have condensation done, otherwise we will lose track of which
    ! ones have had condensation and which have not.
    call aero_state_resort(bin_grid, aero_state)

    post_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          post_water = post_water + aero_particle%vol(aero_data%i_water)
       end do
    end do

    !>DEBUG
    write(*,*) 'RH post cond = ', env_state%rel_humid
    !write(*,*) 'time,min_call_f,max_call_f,mean_call_f = ', env_state%elapsed_time, &
    !     minval(n_call_f), maxval(n_call_f), (real(sum(n_call_f),kind=dp) / real(aero_state%n_part,kind=dp))
    !<DEBUG

    deallocate(cond_kappa)
    deallocate(cond_D_dry)

    call aero_data_deallocate(condense_saved_aero_data)
    call env_data_deallocate(condense_saved_env_data)
    call env_state_deallocate(condense_saved_env_state)

  end subroutine condense_particles_new_unified

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> vapor.
  subroutine condense_particles_unified(bin_grid, env_state, env_data, &
       aero_data, aero_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Environment data.
    type(env_data_t), intent(in) :: env_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t
    
    integer :: i_bin, j, new_bin, k, i_part
    real(kind=dp) :: pv, pre_water, post_water
    type(aero_particle_t), pointer :: aero_particle
    integer :: n_eqn, tol_types, itask, istate, iopt, method_flag
    real(kind=dp) :: state(aero_state%n_part + 1), init_time, final_time, rel_tol(1)
    real(kind=dp) :: abs_tol(1), abs_tol_vector(aero_state%n_part + 1)
    real(kind=dp) :: real_params(PMC_COND_N_REAL_PARAMS)
    integer :: integer_params(PMC_COND_N_INTEGER_PARAMS)
    real(kind=dp) :: rho_water, M_water
    real(kind=dp) :: thermal_conductivity, molecular_diffusion
    type(vode_opts) :: options

    pre_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          pre_water = pre_water + aero_particle%vol(aero_data%i_water)
       end do
    end do

    !>DEBUG
    !write(*,*) '*******************************************************'
    !write(*,*) 'time ', env_state%elapsed_time
    !write(*,*) 'id ', aero_particle%id
    !<DEBUG
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
         * const%air_spec_heat) &
         * sqrt(2d0 * const%pi * const%air_molec_weight &
         / (const%univ_gas_const * env_state%temp))
    real_params(PMC_COND_Z) = 2d0 * molecular_diffusion &
         / const%accom_coeff * sqrt(2d0 * const%pi * M_water &
         / (const%univ_gas_const * env_state%temp))
    real_params(PMC_COND_K_A) = thermal_conductivity
    real_params(PMC_COND_D_V) = molecular_diffusion
    real_params(PMC_COND_V_COMP) = aero_state%comp_vol

    allocate(cond_kappa(aero_state%n_part))
    allocate(cond_D_dry(aero_state%n_part))
    i_part = 0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          i_part = i_part + 1
          aero_particle => aero_state%bin(i_bin)%particle(j)
          cond_kappa(i_part) = aero_particle_solute_kappa(aero_particle, aero_data)
          cond_D_dry(i_part) = vol2diam(aero_particle_solute_volume(aero_particle, aero_data))
          state(i_part + d_offset) = aero_particle_diameter(aero_particle)
          abs_tol_vector(i_part + d_offset) = max(1d-12, &
               1d-8 * (state(i_part + d_offset) - cond_D_dry(i_part)))
       end do
    end do
    if (rh_first) then
       state(1) = env_state%rel_humid
       abs_tol_vector(1) = 1d-10
    else
       state(aero_state%n_part + 1) = env_state%rel_humid
       abs_tol_vector(aero_state%n_part + 1) = 1d-10
    end if

    !>DEBUG
    !allocate(n_call_f(aero_state%n_part))
    !n_call_f = 0
    !i_call_f = 0
    !<DEBUG

    ! set VODE inputs
    n_eqn = aero_state%n_part + 1
    init_time = 0d0
    final_time = del_t
    rel_tol(1) = 1d-8
    !abs_tol(1) = 1d-10
    itask = 1 ! just output val at final_time
    istate = 1 ! first call for this ODE
    iopt = 0 ! no optional inputs
    method_flag = 26 ! stiff (BDF) method, user-supplied sparse Jacobian

    options = set_opts(sparse_j = .true., &
         user_supplied_jacobian = .true., &
         abserr_vector = abs_tol_vector, relerr = rel_tol(1))
    save_real_params = real_params
    call dvode_f90(condense_vode_unified_f, n_eqn, state, init_time, final_time, &
         itask, istate, options, j_fcn = condense_vode_unified_jac)

    if (rh_first) then
       env_state%rel_humid = state(1)
    else
       env_state%rel_humid = state(aero_state%n_part + 1)
    end if

    i_part = 0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          i_part = i_part + 1
          aero_particle => aero_state%bin(i_bin)%particle(j)

          ! translate output back to particle
          aero_particle%vol(aero_data%i_water) = diam2vol(state(i_part + d_offset)) &
               - aero_particle_solute_volume(aero_particle, aero_data)

          ! ensure volumes stay positive
          aero_particle%vol(aero_data%i_water) = max(0d0, &
               aero_particle%vol(aero_data%i_water))
       end do
    end do

    ! We resort the particles in the bins after only all particles
    ! have condensation done, otherwise we will lose track of which
    ! ones have had condensation and which have not.
    call aero_state_resort(bin_grid, aero_state)

    post_water = 0d0
    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(j)
          post_water = post_water + aero_particle%vol(aero_data%i_water)
       end do
    end do

    !>DEBUG
    write(*,*) 'RH post cond = ', env_state%rel_humid
    !write(*,*) 'time,min_call_f,max_call_f,mean_call_f = ', env_state%elapsed_time, &
    !     minval(n_call_f), maxval(n_call_f), (real(sum(n_call_f),kind=dp) / real(aero_state%n_part,kind=dp))
    !<DEBUG

    deallocate(cond_kappa)
    deallocate(cond_D_dry)

  end subroutine condense_particles_unified

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do condensation to all the particles for a given time interval,
  !> including updating the environment to account for the lost
  !> vapor.
  subroutine condense_particles_individual(bin_grid, env_state, env_data, aero_data, &
       aero_state, del_t)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Environment data.
    type(env_data_t), intent(inout) :: env_data
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
    !>DEBUG
    write(*,*) 'condense_particles A: ', aero_data%num_ions(1)
    !<DEBUG

    do i_bin = 1,bin_grid%n_bin
       do j = 1,aero_state%bin(i_bin)%n_part
          !write(*,*) 'i_bin= ', i_bin, 'j =',  j
          !>DEBUG
          !i_call_f = i_call_f + 1
          !<DEBUG
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
    !>DEBUG
    write(*,*) 'RH post cond = ', env_state%rel_humid
    !write(*,*) 'time,min_call_f,max_call_f,mean_call_f = ', env_state%elapsed_time, &
    !     minval(n_call_f), maxval(n_call_f), (real(sum(n_call_f),kind=dp) / real(aero_state%n_part,kind=dp))
    !<DEBUG

    !write(*,*) 'RH (after update due to cond. ) =', env_state%rel_humid 

  end subroutine condense_particles_individual

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
    integer :: constraint_indices(1)
    real(kind=dp) :: constraint_lower(1), constraint_upper(1)

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

    !>DEBUG
    !write(*,*) '*******************************************************'
    !write(*,*) 'time ', env_state%elapsed_time
    !write(*,*) 'id ', aero_particle%id
    !<DEBUG
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
         * const%air_spec_heat) &
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
    real_params(PMC_COND_S) = env_state%rel_humid

    !>DEBUG
    !write(*,*) '************************************************************'
    !write(*,*) 'id ', aero_particle%id
    !write(*,*) 'D_dry ', vol2diam(real_params(PMC_COND_V_DRY))
    !write(*,*) 'time ', env_state%elapsed_time
    !write(*,*) 'RH ', env_state%rel_humid
    !write(*,*) 'temp ', env_state%temp
    !write(*,*) 'PMC_COND_U ', real_params(PMC_COND_U)
    !write(*,*) 'PMC_COND_V ', real_params(PMC_COND_V)
    !write(*,*) 'PMC_COND_W ', real_params(PMC_COND_W)
    !write(*,*) 'PMC_COND_X ', real_params(PMC_COND_X)
    !write(*,*) 'PMC_COND_Y ', real_params(PMC_COND_Y)
    !write(*,*) 'PMC_COND_Z ', real_params(PMC_COND_Z)
    !write(*,*) 'PMC_COND_K_A ', real_params(PMC_COND_K_A)
    !write(*,*) 'PMC_COND_D_V ', real_params(PMC_COND_D_V)
    !write(*,*) 'PMC_COND_KAPPA ', real_params(PMC_COND_KAPPA)
    !write(*,*) 'PMC_COND_V_DRY ', real_params(PMC_COND_V_DRY)
    !write(*,*) 'PMC_COND_S ', real_params(PMC_COND_S)
    !<DEBUG

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
    constraint_indices(1) = 1 ! constrain element 1 of the state
    constraint_lower(1) = vol2diam(real_params(PMC_COND_V_DRY)) ! lower bound
    constraint_upper(1) = 1d0 ! upper bound

    ! call ODE integrator
    !call dvode(condense_vode_f, n_eqn, val, init_time, &
    !     final_time, tol_types, rel_tol, abs_tol, itask, istate, &
    !     iopt, real_work, real_work_len, integer_work, integer_work_len, &
    !     condense_vode_jac, method_flag, real_params, integer_params)
    !if (istate /= 2) then
    !   call die_msg(982335370, "DVODE error code: " &
    !        // trim(integer_to_string(istate)))
    !end if

    !dense_j = .true., &
    !>DEBUG
    !write(*,*) 'clower ', constraint_lower(1)
    !write(*,*) 'cupper ', constraint_upper(1)
    !n_call_f = 0
    !<DEBUG
    options = set_opts(method_flag = method_flag, &
         abserr = abs_tol(1), relerr = rel_tol(1), &
         user_supplied_jacobian = .true.)
    !, constrained = constraint_indices, &
    !     clower = constraint_lower, cupper = constraint_upper)
    save_real_params = real_params
    call dvode_f90(condense_vode_f, n_eqn, val, init_time, final_time, &
         itask, istate, options, j_fcn = condense_vode_jac)

    ! translate output back to particle
    aero_particle%vol(aero_data%i_water) = diam2vol(val(1)) &
         - aero_particle_solute_volume(aero_particle, aero_data)

    ! ensure volumes stay positive
    aero_particle%vol(aero_data%i_water) = max(0d0, &
         aero_particle%vol(aero_data%i_water))

    !>DEBUG
    !write(*,*) 'condense_particle_vode: time, id, n_call_f = ', &
    !     env_state%elapsed_time, aero_particle%id, n_call_f
    !if (kill_flag) then
    !   stop
    !end if
    !<DEBUG
    
  end subroutine condense_particle_vode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_vode_unified_f(n_eqn, time, state, state_dot)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Time derivative of state vector.
    real(kind=dp), intent(out) :: state_dot(n_eqn)

    real(kind=dp) :: diameter(1), diameter_dot(1)
    real(kind=dp) :: rel_humid, rel_humid_dot
    integer :: i_part

    if (rh_first) then
       rel_humid = state(1)
    else
       rel_humid = state(n_eqn)
    end if

    rel_humid_dot = 0d0
    do i_part = 1,(n_eqn - 1)
       diameter(1) = state(i_part + d_offset)
       save_real_params(PMC_COND_KAPPA) = cond_kappa(i_part)
       save_real_params(PMC_COND_V_DRY) = diam2vol(cond_D_dry(i_part))
       save_real_params(PMC_COND_S) = rel_humid
       call condense_vode_f(1, time, diameter, diameter_dot)
       state_dot(i_part + d_offset) = diameter_dot(1)
       rel_humid_dot = rel_humid_dot - 2d0 * const%pi * diameter(1)**2 &
            / (save_real_params(PMC_COND_V) * save_real_params(PMC_COND_V_COMP)) &
            * diameter_dot(1)
    end do

    if (rh_first) then
       state_dot(1) = rel_humid_dot
    else
       state_dot(n_eqn) = rel_humid_dot
    end if

  end subroutine condense_vode_unified_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fills in the \c env_state with the current environment state,
  !> taken from the \c state vector and from global variables.
  subroutine condense_current_env_state(n_eqn, time, state, env_state)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Current environment state.
    type(env_state_t), intent(inout) :: env_state

    call env_state_copy(condense_saved_env_state, env_state)
    call env_data_update_state(condense_saved_env_data, &
         env_state, env_state%elapsed_time + time)

    if (rh_first) then
       env_state%rel_humid = state(1)
    else
       env_state%rel_humid = state(n_eqn)
    end if

  end subroutine condense_current_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_params_global(p, env_state, aero_data, V_comp)

    !> Condensation parameters.
    type(condense_params_t), intent(out) :: p
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Computational volume (m^3).
    real(kind=dp), intent(in) :: V_comp

    p%T = env_state%temp
    p%H = env_state%rel_humid
    p%p = env_state%pressure
    p%V_comp = V_comp

    p%rho_w = aero_particle_water_density(aero_data)
    p%M_w = aero_particle_water_molec_weight(aero_data)
    p%P_0 = env_state_sat_vapor_pressure(env_state)
    p%rho_air = env_state_air_den(env_state)

    p%k_a = 1d-3 * (4.39d0 + 0.071d0 * p%T)
    p%D_v = 0.211d-4 / (p%p / const%air_std_press) * (p%T / 273d0)**1.94d0
    p%U = const%water_latent_heat * p%rho_w / (4d0 * p%T)
    p%V = 4d0 * p%M_w * p%P_0 / (p%rho_w * const%univ_gas_const * p%T)
    p%W = const%water_latent_heat * p%M_w / (const%univ_gas_const * p%T)
    p%X = 4d0 * p%M_w * const%water_surf_eng &
         / (const%univ_gas_const * p%T * p%rho_w) 
    p%Y = 2d0 * p%k_a / (const%accom_coeff * p%rho_air &
         * const%air_spec_heat) &
         * sqrt(2d0 * const%pi * const%air_molec_weight &
         / (const%univ_gas_const * p%T))
    p%Z = 2d0 * p%D_v / const%accom_coeff * sqrt(2d0 * const%pi * p%M_w &
         / (const%univ_gas_const * p%T))

  end subroutine condense_params_global

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_params_per_particle(D, kappa, D_dry, p)

    !> Particle diameter (m).
    real(kind=dp), intent(in) :: D
    !> Kappa (1).
    real(kind=dp), intent(in) :: kappa
    !> Dry particle diameter (m).
    real(kind=dp), intent(in) :: D_dry
    !> Condensation parameters.
    type(condense_params_t), intent(inout) :: p

    real(kind=dp) :: delta

    p%D = D
    p%kappa = kappa
    p%D_dry = D_dry

    p%k_ap = p%k_a / (1d0 + p%Y / p%D)
    p%dkap_dD = p%k_a * p%Y / (p%D + p%Y)**2
    p%D_vp = p%D_v / (1d0 + p%Z / p%D)
    p%dDvp_dD = p%D_v * p%Z / (p%D + p%Z)**2
    if (p%D > p%D_dry) then
       p%a_w = (p%D**3 - p%D_dry**3) &
            / (p%D**3 + (p%kappa - 1d0) * p%D_dry**3)
    else
       p%a_w = 0d0
    end if
    p%daw_dD = 3d0 * p%D**2 * p%kappa * p%D_dry &
         / (p%D**3 + (p%kappa - 1d0) * p%D_dry**3)**2

    delta = 0d0
    call condense_vode_delta_star_newton_new(delta, p)
    p%delta_star = delta
    
    p%Ddot = p%k_ap * p%delta_star / (p%U * p%D)
    p%Hdot = - 2d0 * const%pi * p%D**2 / (p%V * p%V_comp) * p%Ddot

    p%dh_ddelta = condense_vode_implicit_dh_ddelta_new(p%delta_star, p)
    p%dh_dD = p%dkap_dD * p%delta_star &
         - p%U * p%V * p%dDvp_dD * p%H + p%U * p%V &
         * (p%a_w * p%dDvp_dD + p%D_vp * p%daw_dD &
         - p%D_vp * p%a_w * (p%X / p%D**2) / (1d0 + p%delta_star)) &
         * (1d0 / (1d0 + p%delta_star)) &
         * exp((p%W * p%delta_star) / (1d0 + p%delta_star) &
         + (p%X / p%D) / (1d0 + p%delta_star))
    p%dh_dH = - p%U * p%V * p%D_vp

    p%ddeltastar_dD = - p%dh_dD / p%dh_ddelta
    p%ddeltastar_dH = - p%dh_dH / p%dh_ddelta

    p%dDdot_dD = p%dkap_dD * p%delta_star / (p%U * p%D) &
         + p%k_ap * p%ddeltastar_dD / (p%U * p%D) &
         - p%k_ap * p%delta_star / (p%U * p%D**2)
    p%dDdot_dH = p%k_ap / (p%U * p%D) * p%ddeltastar_dH
    p%dHdot_dD = - 2d0 * const%pi / (p%V * p%V_comp) &
         * (2d0 * p%D * p%Ddot + p%D**2 * p%dDdot_dD)
    p%dHdot_dH = - 2d0 * const%pi / (p%V * p%V_comp) * p%D**2 * p%dDdot_dH

  end subroutine condense_params_per_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Function \$ h(\delta,D) \$.
  real(kind=dp) function condense_vode_implicit_h_new(delta, p)

    !> Growth parameter (units???).
    real(kind=dp), intent(in) :: delta
    !> Condensation parameters.
    type(condense_params_t), intent(in) :: p

    condense_vode_implicit_h_new = p%k_ap * delta - p%U * p%V * p%D_vp &
         * (p%H - p%a_w / (1d0 + delta) * exp(p%W * delta / (1d0 + delta) &
         + (p%X / p%D) / (1d0 + delta)))

  end function condense_vode_implicit_h_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Derivative \$ \partial h(\delta,D) / \partial \delta \$.
  real(kind=dp) function condense_vode_implicit_dh_ddelta_new(delta, p)

    !> Growth parameter (units???).
    real(kind=dp), intent(in) :: delta
    !> Condensation parameters.
    type(condense_params_t), intent(in) :: p

    condense_vode_implicit_dh_ddelta_new = &
         p%k_ap - p%U * p%V * p%D_vp * p%a_w / (1d0 + delta)**2 &
         * (1d0 - p%W / (1d0 + delta) + (p%X / p%D) / (1d0 + delta)) &
         * exp(p%W * delta / (1d0 + delta) + (p%X / p%D) / (1d0 + delta))
    
  end function condense_vode_implicit_dh_ddelta_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_vode_new_unified_f(n_eqn, time, state, state_dot)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Time derivative of state vector.
    real(kind=dp), intent(out) :: state_dot(n_eqn)

    real(kind=dp) :: D, kappa, D_dry, Hdot
    real(kind=dp) :: rel_humid, rel_humid_dot
    integer :: i_part
    type(env_state_t) :: env_state
    type(condense_params_t) :: p

    call env_state_allocate(env_state)
    call condense_current_env_state(n_eqn, time, state, env_state)
    call condense_params_global(p, env_state, condense_saved_aero_data, &
         condense_saved_V_comp)
    call env_state_deallocate(env_state)

    Hdot = 0d0
    do i_part = 1,(n_eqn - 1)
       D = state(i_part + d_offset)
       kappa = cond_kappa(i_part)
       D_dry = cond_D_dry(i_part)
       call condense_params_per_particle(D, kappa, D_dry, p)
       state_dot(i_part + d_offset) = p%Ddot
       Hdot = Hdot + p%Hdot
    end do

    if (rh_first) then
       state_dot(1) = Hdot
    else
       state_dot(n_eqn) = Hdot
    end if

  end subroutine condense_vode_new_unified_f

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
    !>DEBUG
    real(kind=dp) :: D, kappa, V_dry
    !<DEBUG

    real_params = save_real_params

    call assert(383442283, n_eqn == 1)
    delta = 0d0
    diameter = state(1)
    !>DEBUG
    !n_call_f(i_call_f) = n_call_f(i_call_f) + 1
    !write(*,*) 'condense_vode_f: t,D,De =    ', time, diameter, &
    !     (diameter - vol2diam(real_params(PMC_COND_V_DRY)))
    !D = diameter
    !kappa = real_params(PMC_COND_KAPPA)
    !V_dry = real_params(PMC_COND_V_DRY)
    !write(*,*) 'const%pi / 6d0 * D**3 = ', (const%pi / 6d0 * D**3)
    !write(*,*) 'V_dry = ', (V_dry)
    !write(*,*) '(kappa - 1d0) = ', ((kappa - 1d0))
    !write(*,*) '(kappa - 1d0) * V_dry = ', ((kappa - 1d0) * V_dry)
    !write(*,*) 'const%pi / 6d0 * D**3 - V_dry = ', (const%pi / 6d0 * D**3 - V_dry)
    !write(*,*) 'const%pi / 6d0 * D**3 + (kappa - 1d0) * V_dry = ', (const%pi / 6d0 * D**3 + (kappa - 1d0) * V_dry)
    !write(*,*) 'a_w ', water_activity(diameter, real_params, integer_params)
    !write(*,*) 'k_ap ', corrected_thermal_conductivity(diameter, real_params, integer_params)
    !write(*,*) 'D_vp ', corrected_molecular_diffusion(diameter, real_params, integer_params)
    !<DEBUG
    call condense_vode_delta_star_newton(delta, diameter, real_params, &
         integer_params)
    k_ap = corrected_thermal_conductivity(diameter, real_params, &
         integer_params)
    state_dot(1) = k_ap * delta / (real_params(PMC_COND_U) * diameter)
    !>DEBUG
    !write(*,*) 'condense_vode_f: t,D,De,Dd = ', time, diameter, &
    !     (diameter - vol2diam(real_params(PMC_COND_V_DRY))), state_dot(1)
    !write(*,*) 'condense_vode_f: t,D = ', time, diameter
    !write(*,*) 'condense_vode_f: delta, Ddot = ', delta, state_dot(1)
    !<DEBUG

    ! save value of delta to be used by condense_vode_jac
    pmc_cond_saved_diameter = diameter
    pmc_cond_saved_delta_star = delta

  end subroutine condense_vode_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_vode_new_unified_jac(n_eqn, time, state, ia, ja, &
       n_non_zero, state_jac)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Non-zero column offsets.
    integer, intent(out) :: ia(*)
    !> Non-zero row locations.
    integer, intent(out) :: ja(*)
    !> Number of non-zero elements in the Jacobian.
    integer, intent(inout) :: n_non_zero
    !> Sparse Jacobian of time derivative of state vector.
    real(kind=dp), intent(out) :: state_jac(*)

    real(kind=dp) :: dDdot_dD(n_eqn - 1), dDdot_dH(n_eqn - 1)
    real(kind=dp) :: dHdot_dD(n_eqn - 1), dHdot_dH
    integer :: i_nz, i_part
    real(kind=dp) :: D, kappa, D_dry
    type(env_state_t) :: env_state
    type(condense_params_t) :: p

    ! signal from vode to initialize number of non-zeros
    if (n_non_zero == 0) then
       n_non_zero = 3 * n_eqn - 2
       return
    end if

    ! if not initializing, this should be correct
    call assert(395158320, n_non_zero == 3 * n_eqn - 2)

    call env_state_allocate(env_state)
    call condense_current_env_state(n_eqn, time, state, env_state)
    call condense_params_global(p, env_state, condense_saved_aero_data, &
         condense_saved_V_comp)
    call env_state_deallocate(env_state)

    dHdot_dH = 0d0
    do i_part = 1,(n_eqn - 1)
       D = state(i_part + d_offset)
       kappa = cond_kappa(i_part)
       D_dry = cond_D_dry(i_part)
       call condense_params_per_particle(D, kappa, D_dry, p)
       dDdot_dD(i_part) = p%dDdot_dD
       dDdot_dH(i_part) = p%dDdot_dH
       dHdot_dD(i_part) = p%dHdot_dD
       dHdot_dH = dHdot_dH + p%dHdot_dH
    end do

    ! Copied from dvode_f90_m documentation:
    !
    ! IA defines the number of nonzeros including the
    ! diagonal in each column of the Jacobian. Define
    ! IA(1) = 1 and for J = 1,..., N,
    ! IA(J+1) = IA(J) + number of nonzeros in column J.
    ! Diagonal elements must be included even if they are
    ! zero. You should check to ensure that IA(N+1)-1 = NZ.
    ! JA defines the rows in which the nonzeros occur.
    ! For I = 1,...,NZ, JA(I) is the row in which the Ith
    ! element of the Jacobian occurs. JA must also include
    ! the diagonal elements of the Jacobian.
    ! PD defines the numerical value of the Jacobian
    ! elements. For I = 1,...,NZ, PD(I) is the numerical
    ! value of the Ith element in the Jacobian. PD must
    ! also include the diagonal elements of the Jacobian.
    if (rh_first) then
       ia(1) = 1
       ia(2) = ia(1) + n_eqn
       do i_part = 1,(n_eqn - 1)
          ia(2 + i_part) = ia(1 + i_part) + 2
       end do
       call assert(351522940, ia(n_eqn + 1) - 1 == n_non_zero)

       i_nz = 0
       i_nz = i_nz + 1
       ja(i_nz) = 1
       state_jac(i_nz) = dHdot_dH
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = 1 + i_part
          state_jac(i_nz) = dDdot_dH(i_part)
       end do
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = 1
          state_jac(i_nz) = dHdot_dD(i_part)
          i_nz = i_nz + 1
          ja(i_nz) = i_part + 1
          state_jac(i_nz) = dDdot_dD(i_part)
       end do
       call assert(260941580, i_nz == n_non_zero)
    else
       ia(1) = 1
       do i_part = 1,(n_eqn - 1)
          ia(1 + i_part) = ia(i_part) + 2
       end do
       ia(n_eqn + 1) = ia(n_eqn) + n_eqn
       call assert(824345254, ia(n_eqn + 1) - 1 == n_non_zero)

       i_nz = 0
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = i_part
          state_jac(i_nz) = dDdot_dD(i_part)
          i_nz = i_nz + 1
          ja(i_nz) = n_eqn
          state_jac(i_nz) = dHdot_dD(i_part)
       end do
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = i_part
          state_jac(i_nz) = dDdot_dH(i_part)
       end do
       i_nz = i_nz + 1
       ja(i_nz) = n_eqn
       state_jac(i_nz) = dHdot_dH
       call assert(901219708, i_nz == n_non_zero)
    end if

  end subroutine condense_vode_new_unified_jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine condense_vode_unified_jac(n_eqn, time, state, ia, ja, &
       n_non_zero, state_jac)

    !> Length of state vector.
    integer, intent(in) :: n_eqn
    !> Current time (s).
    real(kind=dp), intent(in) :: time
    !> Current state vector.
    real(kind=dp), intent(in) :: state(n_eqn)
    !> Non-zero column offsets.
    integer, intent(out) :: ia(*)
    !> Non-zero row locations.
    integer, intent(out) :: ja(*)
    !> Number of non-zero elements in the Jacobian.
    integer, intent(inout) :: n_non_zero
    !> Sparse Jacobian of time derivative of state vector.
    real(kind=dp), intent(out) :: state_jac(*)

    real(kind=dp) :: dDdot_dD(n_eqn - 1), dDdot_dH(n_eqn - 1)
    real(kind=dp) :: dHdot_dD(n_eqn - 1), dHdot_dH
    real(kind=dp) :: diameter(1), diameter_dot(1), dDdot_dD_single(1)
    real(kind=dp) :: delta_star, D, U, k_ap, dkap_dD, D_vp
    real(kind=dp) :: dh_dH, dh_ddelta, ddeltastar_dH
    real(kind=dp) :: rel_humid, V, V_comp
    real(kind=dp) :: real_params(PMC_COND_N_REAL_PARAMS)
    integer :: integer_params(PMC_COND_N_INTEGER_PARAMS)
    integer :: i_nz, i_part

    if (n_non_zero == 0) then
       ! initialization
       n_non_zero = 3 * n_eqn - 2
       return
    end if
    
    call assert(746653549, n_non_zero == 3 * n_eqn - 2)
    !call assert(118353281, size(ia) == n_eqn + 1)
    !call assert(796050890, size(ja) == n_non_zero)
    !call assert(622939472, size(state_jac) == n_non_zero)

    if (rh_first) then
       rel_humid = state(1)
    else
       rel_humid = state(n_eqn)
    end if

    U = save_real_params(PMC_COND_U)
    V = save_real_params(PMC_COND_V)
    V_comp = save_real_params(PMC_COND_V_COMP)

    dHdot_dH = 0d0
    do i_part = 1,(n_eqn - 1)
       diameter(1) = state(i_part + d_offset)
       save_real_params(PMC_COND_KAPPA) = cond_kappa(i_part)
       save_real_params(PMC_COND_V_DRY) = diam2vol(cond_D_dry(i_part))
       save_real_params(PMC_COND_S) = rel_humid
       ! call f to set up saved paramters
       call condense_vode_f(1, time, diameter, diameter_dot)
       call condense_vode_jac(1, time, diameter, 0, 0, dDdot_dD_single, 1)
       dDdot_dD(i_part) = dDdot_dD_single(1)

       real_params = save_real_params
       D = diameter(1)
       call assert_msg(168527146, D == pmc_cond_saved_diameter, &
            "state diameter does not match pmc_cond_saved_diameter")
       delta_star = pmc_cond_saved_delta_star
       k_ap = corrected_thermal_conductivity(D, real_params, &
            integer_params)
       D_vp = corrected_molecular_diffusion(D, real_params, &
            integer_params)
       dkap_dD = corrected_thermal_conductivity_deriv(D, real_params, &
            integer_params)
       dh_dH = - U * V * D_vp
       dh_ddelta = condense_vode_implicit_dh_ddelta(delta_star, D, &
            real_params, integer_params)
       ddeltastar_dH = - dh_dH / dh_ddelta

       dDdot_dH(i_part) = k_ap / (U * D) * ddeltastar_dH
       dHdot_dD(i_part) = - 2d0 * const%pi / (V * V_comp) &
            * (2d0 * diameter(1) * diameter_dot(1) &
            + diameter(1)**2 * dDdot_dD(i_part))
       dHdot_dH = dHdot_dH - 2d0 * const%pi / (V * V_comp) &
            * diameter(1)**2 * dDdot_dH(i_part)
    end do

    ! Copied from dvode_f90_m documentation:
    !
    ! IA defines the number of nonzeros including the
    ! diagonal in each column of the Jacobian. Define
    ! IA(1) = 1 and for J = 1,..., N,
    ! IA(J+1) = IA(J) + number of nonzeros in column J.
    ! Diagonal elements must be included even if they are
    ! zero. You should check to ensure that IA(N+1)-1 = NZ.
    ! JA defines the rows in which the nonzeros occur.
    ! For I = 1,...,NZ, JA(I) is the row in which the Ith
    ! element of the Jacobian occurs. JA must also include
    ! the diagonal elements of the Jacobian.
    ! PD defines the numerical value of the Jacobian
    ! elements. For I = 1,...,NZ, PD(I) is the numerical
    ! value of the Ith element in the Jacobian. PD must
    ! also include the diagonal elements of the Jacobian.
    if (rh_first) then
       ia(1) = 1
       ia(2) = ia(1) + n_eqn
       do i_part = 1,(n_eqn - 1)
          ia(2 + i_part) = ia(1 + i_part) + 2
       end do
       call assert(742852055, ia(n_eqn + 1) - 1 == n_non_zero)

       i_nz = 0
       i_nz = i_nz + 1
       ja(i_nz) = 1
       state_jac(i_nz) = dHdot_dH
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = 1 + i_part
          state_jac(i_nz) = dDdot_dH(i_part)
       end do
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = 1
          state_jac(i_nz) = dHdot_dD(i_part)
          i_nz = i_nz + 1
          ja(i_nz) = i_part + 1
          state_jac(i_nz) = dDdot_dD(i_part)
       end do
       call assert(885711184, i_nz == n_non_zero)
    else
       ia(1) = 1
       do i_part = 1,(n_eqn - 1)
          ia(1 + i_part) = ia(i_part) + 2
       end do
       ia(n_eqn + 1) = ia(n_eqn) + n_eqn
       call assert(915113988, ia(n_eqn + 1) - 1 == n_non_zero)

       i_nz = 0
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = i_part
          state_jac(i_nz) = dDdot_dD(i_part)
          i_nz = i_nz + 1
          ja(i_nz) = n_eqn
          state_jac(i_nz) = dHdot_dD(i_part)
       end do
       do i_part = 1,(n_eqn - 1)
          i_nz = i_nz + 1
          ja(i_nz) = i_part
          state_jac(i_nz) = dDdot_dH(i_part)
       end do
       i_nz = i_nz + 1
       ja(i_nz) = n_eqn
       state_jac(i_nz) = dHdot_dH
       call assert(723317905, i_nz == n_non_zero)
    end if

  end subroutine condense_vode_unified_jac

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
  subroutine condense_vode_delta_star_newton_new(delta, p)
    
    !> Variable that we are solving for.
    real(kind=dp), intent(inout) :: delta
    !> Condensation parameters.
    type(condense_params_t), intent(in) :: p

    !> Delta convergence tolerance.
    real(kind=dp), parameter :: delta_tol = 1d-16
    !> Function convergence tolerance.
    real(kind=dp), parameter :: h_tol = 1d-10
    !> Maximum number of iterations.
    integer, parameter :: iter_max = 100

    integer :: iter
    real(kind=dp) :: h, dh_ddelta, old_h, delta_step, h_step
    real(kind=dp) :: check_delta, check_h, finite_diff_dh, rel_error

    h = condense_vode_implicit_h_new(delta, p)
    dh_ddelta = condense_vode_implicit_dh_ddelta_new(delta, p)
    old_h = h

    iter = 0
    do
       iter = iter + 1
       delta_step = - h / dh_ddelta
       
       delta = delta + delta_step
       h = condense_vode_implicit_h_new(delta, p)
       dh_ddelta = condense_vode_implicit_dh_ddelta_new(delta, p)
       h_step = h - old_h
       old_h = h

       if (PMC_COND_CHECK_DERIVS) then
          check_delta = delta * (1d0 + PMC_COND_CHECK_EPS)
          check_h = condense_vode_implicit_h_new(check_delta, p)
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

             write(0,*) 'ERROR: deriv rel_error = ', rel_error
          end if
       end if

       if (iter .ge. iter_max) then
          call die_msg(720859756, 'Newton iteration failed to terminate')
       end if
       
       ! FIXME: gfortran 4.1.1 requires the "then" in the following
       ! statement, rather than using a single-line "if" statement.
       if ((abs(delta_step) .lt. delta_tol) &
            .and. (abs(h_step) .lt. h_tol)) then
          exit
       end if
    end do
 
  end subroutine condense_vode_delta_star_newton_new

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

    !>DEBUG
    !write(*,*) 'PMC_COND_U ', real_params(PMC_COND_U)
    !write(*,*) 'PMC_COND_V ', real_params(PMC_COND_V)
    !write(*,*) 'PMC_COND_W ', real_params(PMC_COND_W)
    !write(*,*) 'PMC_COND_X ', real_params(PMC_COND_X)
    !write(*,*) 'PMC_COND_Y ', real_params(PMC_COND_Y)
    !write(*,*) 'PMC_COND_Z ', real_params(PMC_COND_Z)
    !write(*,*) 'PMC_COND_K_A ', real_params(PMC_COND_K_A)
    !write(*,*) 'PMC_COND_D_V ', real_params(PMC_COND_D_V)
    !write(*,*) 'PMC_COND_KAPPA ', real_params(PMC_COND_KAPPA)
    !write(*,*) 'PMC_COND_V_DRY ', real_params(PMC_COND_V_DRY)
    !write(*,*) 'PMC_COND_S ', real_params(PMC_COND_S)
    !write(*,*) 'a_w ', water_activity(diameter, real_params, integer_params)
    !write(*,*) 'k_ap ', corrected_thermal_conductivity(diameter, real_params, integer_params)
    !write(*,*) 'D_vp ', corrected_molecular_diffusion(diameter, real_params, integer_params)
    !write(*,*) 'D_dry ', vol2diam(real_params(PMC_COND_V_DRY))
    !write(*,'(a5,a15,a15,a15)') 'iter', 'delta', 'h', 'dh_ddelta'
    !do iter = -100,100
    !   check_delta = real(iter, kind=dp) / 100d0
    !   h = condense_vode_implicit_h(check_delta, diameter, real_params, &
    !        integer_params)
    !   dh_ddelta = condense_vode_implicit_dh_ddelta(check_delta, diameter, &
    !        real_params, integer_params)
    !   write(*,'(i5,e15.5,e15.5,e15.5)') iter, check_delta, h, dh_ddelta
    !end do
    !<DEBUG

    h = condense_vode_implicit_h(delta, diameter, real_params, &
         integer_params)
    dh_ddelta = condense_vode_implicit_dh_ddelta(delta, diameter, &
         real_params, integer_params)
    old_h = h

    iter = 0
    do
       iter = iter + 1
       delta_step = - h / dh_ddelta
       !>DEBUG
       !write(*,'(a5,a15,a15,a15,a15)') 'iter', 'delta', 'h', 'dh_ddelta', 'delta_step'
       !write(*,'(i5,e15.5,e15.5,e15.5,e15.5)') iter, delta, h, dh_ddelta, delta_step
       !<DEBUG
       
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

             write(0,*) 'ERROR: deriv rel_error = ', rel_error
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

    if (const%pi / 6d0 * D**3 <= V_dry) then
       water_activity = 0d0
       !>DEBUG
       !kill_flag = .true.
       !<DEBUG
    end if

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_condensation
