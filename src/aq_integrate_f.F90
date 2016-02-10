! Copyright (C) 2015 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_integrate_f module.

!> Interface module for c ODE solver

module pmc_aq_integrate_f

  use pmc_aq_mech_data
  use pmc_aq_spec_data
  use pmc_aq_state
  use pmc_env_state

  use iso_c_binding

  implicit none

  !> Relative tolerance for integration
  real(kind=dp), parameter :: AQ_INTEGRATE_RELTOL = 1.0D-5
  !> Maximum number of iterations for solver
  integer, parameter :: AQ_INTEGRATE_MAX_STEPS = 10000000

  !> Result code indicating successful completion.
  integer, parameter :: PMC_AQ_INTEGRATE_SUCCESS         = 0
  !> Result code indicating no available integration routine
  integer, parameter :: PMC_AQ_INTEGRATE_NO_AVAIL_SOLVER = 12

#ifdef PMC_USE_SUNDIALS
  !> Result code indicating failure to allocate \c y vector.
  integer, parameter :: PMC_AQ_INTEGRATE_INIT_Y          = 1
  !> Result code indicating failure to allocate \c abstol vector.
  integer, parameter :: PMC_AQ_INTEGRATE_INIT_ABSTOL     = 2
  !> Result code indicating failure to create the solver.
  integer, parameter :: PMC_AQ_INTEGRATE_INIT_CVODE_MEM  = 3
  !> Result code indicating failure to initialize the solver.
  integer, parameter :: PMC_AQ_INTEGRATE_INIT_CVODE      = 4
  !> Result code indicating failure to set tolerances.
  integer, parameter :: PMC_AQ_INTEGRATE_SVTOL           = 5
  !> Result code indicating failure to set maximum steps.
  integer, parameter :: PMC_AQ_INTEGRATE_SET_MAX_STEPS   = 6
  !> Result code indicating failure of the solver.
  integer, parameter :: PMC_AQ_INTEGRATE_FAIL            = 7
  !> Result code indicating failure to set dense Jacobian solver
  integer, parameter :: PMC_AQ_INTEGRATE_DENSE_JAC       = 8
  !> Result code indicating failure to set Jacobian function
  integer, parameter :: PMC_AQ_INTEGRATE_JAC_FUNC        = 9
  !> Result code indicating failure to set user data
  integer, parameter :: PMC_AQ_INTEGRATE_SET_USER_DATA   = 10
  !> Result code indicating SUNDIALS realtype is not set to double precision
  integer, parameter :: PMC_AQ_INTEGRATE_WRONG_PRECISION = 11
#endif

  !> Structure containing required system data for use during calls 
  !> to the ODE solver
  type aq_integrate_data_t
    !> Initial temperature
    real(kind=dp) :: temp_initial
    !> Final temperature
    real(kind=dp) :: temp_final
    !> Initial solar zenith angle (radians)
    real(kind=dp) :: solar_zenith_angle_initial
    !> Final solar zenith angle (radians)
    real(kind=dp) :: solar_zenith_angle_final
    !> Interation time
    real(kind=dp) :: del_t
    !> Particle radius (m)
    real(kind=dp) :: radius
    !> Particle number concentration (#/cc)
    real(kind=dp) :: n_particle
    !> Pointer to mechanism
    type(aq_mech_data_t), pointer :: aq_mech_data
    !> Pointer to species data
    type(aq_spec_data_t), pointer :: aq_spec_data
  end type aq_integrate_data_t

#ifndef DOXYGEN_SKIP_DOC
  ! Interface to c ODE solver function
  interface
    integer(kind=c_int) function aq_integrate_solver(neq, x, abstol, reltol, &
                                 max_steps, t_initial, t_final, sysdata) bind(c)
        use iso_c_binding
        !> number of equations (i.e. species)
        integer(kind=c_int), value :: neq
        !> aq. chemistry state (i.e. concentrations)
        type(c_ptr), value :: x
        !> absolute tolerance for each species
        type(c_ptr), value :: abstol
        !> relative tolerance for all species
        real(kind=c_double), value :: reltol
        !> maximum number of iterations for solver
        integer(kind=c_int), value :: max_steps
        !> Initial time (s)
        real(kind=c_double), value :: t_initial
        !> Final time (s)
        real(kind=c_double), value :: t_final
        !> Pointer to system data required by integration subfunctions
        type(c_ptr), value :: sysdata
    end function aq_integrate_solver
  end interface
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do aq. phase chemistry and phase transfer over a given time interval
  subroutine aq_integrate(aq_state, aq_mech_data, aq_spec_data, env_state_initial, &
       env_state_final, del_t)

    !> Aq. Phase Species State
    type(aq_state_t), intent(inout) :: aq_state
    !> Aq. Phase Mechanism data.
    type(aq_mech_data_t), intent(in), target :: aq_mech_data
    !> Aq. Phase Species data.
    type(aq_spec_data_t), intent(in), target :: aq_spec_data
    !> Environment state at the start of the timestep.
    type(env_state_t), intent(in) :: env_state_initial
    !> Environment state at the end of the timestep.
    type(env_state_t), intent(in) :: env_state_final
    !> Total time to integrate.
    real(kind=dp), intent(in) :: del_t

    ! single variable to reference system data in ODE solver
    type(aq_integrate_data_t),target :: aq_integrate_data

    ! Variables requried by c ODE solver function
    ! as c_* types
    ! number of equations (i.e. species)
    integer(kind=c_int) :: neq_c
    ! aq. chemistry state (i.e. concentrations)
    type(c_ptr) :: state_c_p
    ! absolute tolerance for each species
    type(c_ptr) :: abstol_c_p
    ! relative tolerance for all species
    real(kind=c_double) :: reltol_c
    ! Initial time (s)
    real(kind=c_double) :: t_initial_c
    ! Final time (s)
    real(kind=c_double) :: t_final_c
    ! Pointer to system data
    type(c_ptr) :: sysdata_c_p

    ! Mixing ratio vector
    real(kind=c_double), target :: state_c(aq_spec_data%n_spec)
    ! Absolute tolerance vector
    real(kind=c_double), target :: abstol_c(aq_spec_data%n_spec)

    ! Current temperature, used in case of integration failure
    real(kind=dp) :: curr_temp

    ! Return value from c ODE solver function
    integer(kind=c_int) :: solver_stat

    ! Set up aq_integrate variable with system data
    aq_integrate_data%temp_initial = env_state_initial%temp
    aq_integrate_data%temp_final   = env_state_final%temp
    aq_integrate_data%solar_zenith_angle_initial = env_state_initial%solar_zenith_angle
    aq_integrate_data%solar_zenith_angle_final = env_state_final%solar_zenith_angle
    aq_integrate_data%del_t = del_t
    aq_integrate_data%radius = aq_state%radius
    aq_integrate_data%n_particle = aq_state%n_particle
    aq_integrate_data%aq_mech_data => aq_mech_data
    aq_integrate_data%aq_spec_data => aq_spec_data

    ! Set up c ODE solver variables
    neq_c       = int(aq_spec_data%n_spec, kind=c_int)
    state_c(:)  = real(aq_state%mix_rat(:), c_double)
    state_c_p   = c_loc(state_c)
    abstol_c(:) = real(aq_spec_data%abstol(:), c_double)
    abstol_c_p  = c_loc(abstol_c)
    reltol_c    = real(AQ_INTEGRATE_RELTOL, c_double)
    t_initial_c = real(0, c_double)
    t_final_c   = real(del_t, c_double)
    sysdata_c_p = c_loc(aq_integrate_data)

    ! Call c ODE solver function
    solver_stat = aq_integrate_solver(neq_c, state_c_p, abstol_c_p, reltol_c, &
                    AQ_INTEGRATE_MAX_STEPS, t_initial_c, t_final_c, sysdata_c_p)

    ! DIAGNOSTIC - Print system info on integration failure
    if (solver_stat.ne.PMC_AQ_INTEGRATE_SUCCESS) then
        write(*,*)  " "
        write(*,*)  "deltaT (s) = ", aq_integrate_data%del_t
        write(*,*)  "Temp init = ", aq_integrate_data%temp_initial
        write(*,*)  "Temp final = ", aq_integrate_data%temp_final
        write(*,*)  "SZA init = ", aq_integrate_data%solar_zenith_angle_initial
        write(*,*)  "SZA final = ", aq_integrate_data%solar_zenith_angle_final
        write(*,*)  "[Particles] (#/cc) = ", aq_integrate_data%n_particle
        write(*,*)  "Particle Radius (m) = ", aq_integrate_data%radius
        curr_temp = (aq_integrate_data%temp_initial + aq_integrate_data%temp_final)/2.0
        write(*,*)  "HL conversion factor (atm*L/mol; approx.) = ", &
            aq_integrate_get_gas_conversion_factor(aq_integrate_data, curr_temp)
        write(*,*)  " "
        call pmc_aq_state_print(aq_state, aq_spec_data)
        call aq_mech_data_print(aq_mech_data, aq_spec_data)
    endif

    ! Check for errors
    call aq_chem_check_solve(solver_stat)

    ! Update aq. chemistry state
    aq_state%mix_rat(:) = real(state_c(:), dp)

  end subroutine aq_integrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the return code from the condense_solver() function.
  subroutine aq_chem_check_solve(value)

    !> Return code to check.
    integer(kind=c_int), intent(in) :: value


    if (value == PMC_AQ_INTEGRATE_SUCCESS) then
       return
    elseif (value == PMC_AQ_INTEGRATE_NO_AVAIL_SOLVER) then
       call die_msg(592214019, "aq_integrate: " &
            // "no available solver")

#ifdef PMC_USE_SUNDIALS
    elseif (value == PMC_AQ_INTEGRATE_INIT_Y) then
       call die_msg(592432510, "aq_integrate: " &
            // "failed to allocate y vector")
    elseif (value == PMC_AQ_INTEGRATE_INIT_ABSTOL) then
       call die_msg(592617387, "aq_integrate: " &
            // "failed to allocate abstol vector")
    elseif (value == PMC_AQ_INTEGRATE_INIT_CVODE_MEM) then
       call die_msg(592785457, "aq_integrate: " &
            // "failed to create the solver")
    elseif (value == PMC_AQ_INTEGRATE_INIT_CVODE) then
       call die_msg(592970334, "aq_integrate: " &
            // "failure to initialize the solver")
    elseif (value == PMC_AQ_INTEGRATE_SVTOL) then
       call die_msg(593172018, "aq_integrate: " &
            // "failed to set tolerances")
    elseif (value == PMC_AQ_INTEGRATE_SET_MAX_STEPS) then
       call die_msg(593810684, "aq_integrate: " &
            // "failed to set maximum steps")
    elseif (value == PMC_AQ_INTEGRATE_FAIL) then
       call die_msg(594113210, "aq_integrate: solver failed")
    elseif (value==PMC_AQ_INTEGRATE_DENSE_JAC) then
       call die_msg(594281280, "aq_integrate: " &
            // "failed to set dense jacobian solver")
    elseif (value== PMC_AQ_INTEGRATE_JAC_FUNC) then
       call die_msg(594415736, "aq_integrate: " &
            // "failed to set jacobian function")
    elseif (value==PMC_AQ_INTEGRATE_SET_USER_DATA) then
       call die_msg(594533385, "aq_integrate: " &
            // "failed to set user data")
    elseif (value==PMC_AQ_INTEGRATE_WRONG_PRECISION) then
       call die_msg(594651034, "aq_integrate: " &
            // "SUNDIALS was not compiled for double precision variables")
#endif

    else
       call die_msg(594768683, "aq_integrate: unknown return code: " &
            // trim(integer_to_string(value)))
    end if

    return

  end subroutine aq_chem_check_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether a solver is available for the aqueous chemistry mechanism
  logical function aq_integrate_is_solver_available()

#ifdef PMC_USE_SUNDIALS
    aq_integrate_is_solver_available = .true.
#else
    aq_integrate_is_solver_available = .false.
#endif

  end function aq_integrate_is_solver_available

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calcualte f(t,y) for the aqueous-phase chemical mechanism
  subroutine aq_integrate_f(n_eqn_c, curr_time_c, state_c_p, f_c_p, &
                              sysdata_c_p) bind(c)

    use iso_c_binding

    !> number of equations (i.e. species)
    integer(kind=c_int), intent(in), value :: n_eqn_c
    !> current integration time (s)
    real(kind=c_double), intent(in), value :: curr_time_c
    !> current aq. chemistry state (i.e. mixing ratios)
    type(c_ptr), intent(in), value :: state_c_p
    !> f(t,y) to be calculated
    type(c_ptr), intent(in), value :: f_c_p
    !> pointer to system data in an aq_integrate_t variable
    type(c_ptr), intent(in), value :: sysdata_c_p

    ! Equivalent fortran variables
    ! number of equations (i.e. species)
    integer :: n_eqn
    ! current integration time (s)
    real(kind=dp) :: curr_time
    ! current aq. chemistry state (i.e. mixing ratios)
    real(kind=dp), pointer :: state_p(:)
    ! f(t,y) to be calculated
    real(kind=dp), pointer :: f_p(:)
    ! pointer to system data in an aq_integrate_t variable
    type(aq_integrate_data_t), pointer :: sysdata_p

    ! forward and backward rate constants
    real(kind=dp) :: rc_forward, rc_backward
    ! current temperature (K)
    real(kind=dp) :: curr_temp
    ! current solar zenith angle (radians)
    real(kind=dp) :: curr_solar_zenith_angle
    ! conversion factor for gas-phase species in phase-transfer reactions
    real(kind=dp) :: conv_factor, curr_conv_factor

    integer :: i, j, species_index
    real(kind=dp) :: forward_rate, backward_rate
    real(kind=dp) :: ret_val

    ! Map c variables to fortran variables
    n_eqn = int(n_eqn_c)
    curr_time = real(curr_time_c, dp)
    call c_f_pointer(state_c_p, state_p, (/ n_eqn /))
    call c_f_pointer(f_c_p, f_p, (/ n_eqn /))
    call c_f_pointer(sysdata_c_p, sysdata_p)

    ! Reset f(t,y)
    f_p(:) = 0.0

    ! Interpolate current temperature and solar zenith angle
    curr_temp = sysdata_p%temp_initial + (sysdata_p%temp_final-sysdata_p%temp_initial) &
                    * curr_time / sysdata_p%del_t
    curr_solar_zenith_angle = sysdata_p%solar_zenith_angle_initial + &
                    (sysdata_p%solar_zenith_angle_final-sysdata_p%solar_zenith_angle_initial) &
                    * curr_time / sysdata_p%del_t

    ! Get conversion factor for gas-phase species in phase-transfer reactions
    conv_factor = aq_integrate_get_gas_conversion_factor(sysdata_p, curr_temp)

    ! Calculate f(t,y)
    do i=1, sysdata_p%aq_mech_data%n_rxn

        ! Get rate constant for reaction i
        ret_val = aq_rxn_data_get_rate_constant(rc_forward, rc_backward, &
                    sysdata_p%aq_mech_data%rxn(i), sysdata_p%aq_spec_data, &
                    curr_temp, curr_solar_zenith_angle, sysdata_p%radius)
        if (ret_val.ne.0) then
            call die_msg(595138437, "error getting rate constant for rxn: " // &
                         trim(integer_to_string(i)))
        endif

        ! Calculate the forward rate for rxn i
        forward_rate = rc_forward
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%reactant)
            forward_rate = forward_rate * state_p(sysdata_p%aq_mech_data%rxn(i)%reactant(j))
        enddo

        ! Calculate the backward rate for rxn i
        backward_rate = rc_backward
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%product)
            backward_rate = backward_rate * state_p(sysdata_p%aq_mech_data%rxn(i)%product(j))
        enddo

        ! If this is a phase-transfer reaction, the change in gas-phase
        ! reactant concentration must be converted to units of atm/s
        if (trim(aq_rxn_data_get_class_name(sysdata_p%aq_mech_data%rxn(i)%class_index)) &
            .eq. "HENRY") then
            curr_conv_factor = conv_factor
        else
            curr_conv_factor = 1.0
        endif

        ! Set production and destruction terms for reactant species
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%reactant)
            species_index = sysdata_p%aq_mech_data%rxn(i)%reactant(j)
            if (.not.sysdata_p%aq_spec_data%const_conc(species_index)) then
                f_p(species_index) = f_p(species_index) + (-forward_rate + &
                        backward_rate) * curr_conv_factor
            endif
        enddo

        ! Set production and destruction terms for product species
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%product)
            species_index = sysdata_p%aq_mech_data%rxn(i)%product(j)
            if (.not.sysdata_p%aq_spec_data%const_conc(species_index)) then
                f_p(species_index) = f_p(species_index) &
                    + forward_rate * sysdata_p%aq_mech_data%rxn(i)%prod_yield(j) &
                    - backward_rate
            endif
        enddo

    enddo

  end subroutine aq_integrate_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate J(t,y) = df/dy for the aqueous-phase chemical mechanism
  subroutine aq_integrate_jac(n_eqn_c, curr_time_c, state_c_p, jac_c_p, &
                              sysdata_c_p) bind(c)

    use iso_c_binding

    !> number of equations (i.e. species)
    integer(kind=c_int), intent(in), value :: n_eqn_c
    !> current integration time (s)
    real(kind=c_double), intent(in), value :: curr_time_c
    !> current aq. chemistry state (i.e. mixing ratios)
    type(c_ptr), intent(in), value :: state_c_p
    !> jacobian matrix to be calculated
    type(c_ptr), intent(in), value :: jac_c_p
    !> pointer to system data in an aq_integrate_t variable
    type(c_ptr), intent(in), value :: sysdata_c_p

    ! Equivalent fortran variables
    ! number of equations (i.e. species)
    integer :: n_eqn
    ! current integration time (s)
    real(kind=dp) :: curr_time
    ! current aq. chemistry state (i.e. mixing ratios)
    real(kind=dp), pointer :: state_p(:)
    ! jacobian matrix to be calculated
    real(kind=dp), pointer :: jac_p(:,:)
    ! pointer to system data in an aq_integrate_t variable
    type(aq_integrate_data_t), pointer :: sysdata_p

    ! forward and backward rate constants
    real(kind=dp) :: rc_forward, rc_backward
    ! current temperature (K)
    real(kind=dp) :: curr_temp
    ! current solar zenith angle (radians)
    real(kind=dp) :: curr_solar_zenith_angle
    ! conversion factor for gas-phase species in phase-transfer reactions
    real(kind=dp) :: conv_factor, curr_conv_factor

    integer :: i, j, k, species_1_index, species_2_index
    real(kind=dp) :: forward_rate, backward_rate
    real(kind=dp) :: ret_val

    ! Map c variables to fortran variables
    n_eqn = int(n_eqn_c)
    curr_time = real(curr_time_c, dp)
    call c_f_pointer(state_c_p, state_p, (/ n_eqn /))
    call c_f_pointer(jac_c_p, jac_p, (/ n_eqn, n_eqn /))
    call c_f_pointer(sysdata_c_p, sysdata_p)

    ! Reset J(t,y)
    jac_p(:,:) = 0.0

    ! Interpolate current temperature and solar zenith angle
    curr_temp = sysdata_p%temp_initial + (sysdata_p%temp_final-sysdata_p%temp_initial) &
                    * curr_time / sysdata_p%del_t
    curr_solar_zenith_angle = sysdata_p%solar_zenith_angle_initial + &
                    (sysdata_p%solar_zenith_angle_final-sysdata_p%solar_zenith_angle_initial) &
                    * curr_time / sysdata_p%del_t

    ! Get conversion factor for gas-phase species in phase-transfer reactions
    conv_factor = aq_integrate_get_gas_conversion_factor(sysdata_p, curr_temp)

    ! Calculate J(t,y)
    do i=1, sysdata_p%aq_mech_data%n_rxn

        ! Get rate constant for reaction i
        ret_val = aq_rxn_data_get_rate_constant(rc_forward, rc_backward, &
                    sysdata_p%aq_mech_data%rxn(i), sysdata_p%aq_spec_data, &
                    curr_temp, curr_solar_zenith_angle, sysdata_p%radius)
        if (ret_val.ne.0) then
            call die_msg(595709875, "error getting rate constant for rxn: " // &
                         trim(integer_to_string(i)))
        endif

        ! Calculate the forward rate for rxn i
        forward_rate = rc_forward
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%reactant)
            forward_rate = forward_rate * state_p(sysdata_p%aq_mech_data%rxn(i)%reactant(j))
        enddo

        ! Calculate the backward rate for rxn i
        backward_rate = rc_backward
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%product)
            backward_rate = backward_rate * state_p(sysdata_p%aq_mech_data%rxn(i)%product(j))
        enddo

        ! If this is a phase-transfer reaction, the change in gas-phase
        ! reactant concentration must be converted to units of atm/s
        if (trim(aq_rxn_data_get_class_name(sysdata_p%aq_mech_data%rxn(i)%class_index)) &
            .eq. "HENRY") then
            curr_conv_factor = conv_factor
        else
            curr_conv_factor = 1.0
        endif

        ! Set partial production and destruction terms for reactant species by ...
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%reactant)
            species_1_index = sysdata_p%aq_mech_data%rxn(i)%reactant(j)
            ! ... each reactant
            if (forward_rate.ne.0.0 .and. &
                .not.sysdata_p%aq_spec_data%const_conc(species_1_index)) then
                do k=1, size(sysdata_p%aq_mech_data%rxn(i)%reactant)
                    species_2_index = sysdata_p%aq_mech_data%rxn(i)%reactant(k)
                    jac_p(species_1_index,species_2_index) = &
                        jac_p(species_1_index,species_2_index) - &
                        forward_rate / state_p(species_2_index) * curr_conv_factor
                enddo
            endif
            ! ... each product
            if (backward_rate.ne.0.0 .and. &
                .not.sysdata_p%aq_spec_data%const_conc(species_1_index)) then
                do k=1, size(sysdata_p%aq_mech_data%rxn(i)%product)
                    species_2_index = sysdata_p%aq_mech_data%rxn(i)%product(k)
                    jac_p(species_1_index,species_2_index) = &
                        jac_p(species_1_index,species_2_index) + &
                        backward_rate / state_p(species_2_index) * curr_conv_factor
                enddo
            endif
        enddo

        ! Set partial production and destruction terms for product species by ...
        do j=1, size(sysdata_p%aq_mech_data%rxn(i)%product)
            species_1_index = sysdata_p%aq_mech_data%rxn(i)%product(j)
            !> ... each reactant
            if (forward_rate.ne.0.0 .and. &
                .not.sysdata_p%aq_spec_data%const_conc(species_1_index)) then
                do k=1, size(sysdata_p%aq_mech_data%rxn(i)%reactant)
                    species_2_index = sysdata_p%aq_mech_data%rxn(i)%reactant(k)
                    jac_p(species_1_index,species_2_index) = &
                        jac_p(species_1_index,species_2_index) + &
                        forward_rate / state_p(species_2_index) &
                        * sysdata_p%aq_mech_data%rxn(i)%prod_yield(j)
                enddo
            endif
            ! ... each product
            if (backward_rate.ne.0.0 .and. &
                .not.sysdata_p%aq_spec_data%const_conc(species_1_index)) then
                do k=1, size(sysdata_p%aq_mech_data%rxn(i)%product)
                    species_2_index = sysdata_p%aq_mech_data%rxn(i)%product(k)
                    jac_p(species_1_index,species_2_index) = &
                        jac_p(species_1_index,species_2_index) - &
                        backward_rate / state_p(species_2_index)
                enddo
            endif
        enddo

    enddo

  end subroutine aq_integrate_jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate conversion factor for gas phase species in phase-transfer reactions
  real(kind=dp) function aq_integrate_get_gas_conversion_factor(sysdata_p, temp)

    !> System data
    type(aq_integrate_data_t), pointer, intent(in) :: sysdata_p
    !> current temperature (K)
    real(kind=dp), intent(in) :: temp
    !> Gas constant (cc*atm/K/mol)
    real(kind=dp), parameter :: R_gas = 82.05736

    !> Henry's Law phase transfer reaction rates are in units of M/s. In
    !! order to calculate the change in gas-phase concentrations, this 
    !! must be converted to units of atm/s. 
    !! (see Schwartz, "Mass transport considerations pertinent to aqueous 
    !!  phase reactions of gases in liquid water clouds", in Chemistry of 
    !!  Atmospheric Systems, NATO ASI vol. 6, 1986 pg 451)
    !!
    !! d[G]/dt(atm/s) = d[G]/dt(mol/L/s) * conv_factor
    !! conv_factor = 1000(L/m^3) * 4/3pi*r^3(m^3/particle)
    !!       * [particle](particle/cc) * R_gas(cc*atm/K/mol) * T(K)
    !!
    aq_integrate_get_gas_conversion_factor = 1000.0 * 4.0/3.0 * const%pi &
            * sysdata_p%radius**3.0 * sysdata_p%n_particle * R_gas * temp

  end function aq_integrate_get_gas_conversion_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aq_integrate_f






