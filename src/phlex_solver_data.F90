! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_phlex_solver_data module.

!> The phlex_solver_data_t structure and associated subroutines.
module pmc_phlex_solver_data
#define PHLEX_SOLVER_SUCCESS 0
#define PHLEX_SOLVER_FAIL 1

  use pmc_constants,                   only : i_kind, dp
  use pmc_util,                        only : assert_msg, to_string, die_msg
  use pmc_phlex_state
  use pmc_mechanism_data
  use pmc_rxn_data
  use pmc_rxn_factory

  use iso_c_binding

  implicit none
  private

  public :: phlex_solver_data_t

  !> Default relative tolerance for integration
  real(kind=dp), parameter :: PMC_SOLVER_DEFAULT_REL_TOL = 1.0D-8
  !> Default max number of integration steps
  integer(kind=i_kind), parameter :: PMC_SOLVER_DEFAULT_MAX_STEPS = 10000000
  !> Default maximum number of integration convergence failures
  integer(kind=i_kind), parameter :: PMC_SOLVER_DEFAULT_MAX_CONV_FAILS = 10

  !> Result code indicating successful completion
  integer, parameter :: PMC_SOLVER_SUCCESS = 0

  !> Interface to c ODE solver functions
  interface
    !> Get a new solver 
    type(c_ptr) function solver_new(n_state_var, var_type, &
                    n_rxn, n_int_param, n_float_param) bind (c)
      use iso_c_binding
      !> Number of variables on the state array (including const, PSSA, etc.)
      integer(kind=c_int), value :: n_state_var
      !> Pointer to array of state variable types (solver, constant, PSSA)
      type(c_ptr), value :: var_type
      !> Number of reactions to solve
      integer(kind=c_int), value :: n_rxn
      !> Total number of integer parameters for all reactions
      integer(kind=c_int), value :: n_int_param
      !> Total number of floating-point parameters for all reactions
      integer(kind=c_int), value :: n_float_param
    end function solver_new

    !> Solver initialization
    subroutine solver_initialize(solver_data, abs_tol, rel_tol, max_steps, &
                    max_conv_fails) bind (c)
      use iso_c_binding
      !> Pointer to a SolverData object
      type(c_ptr), value :: solver_data
      !> Pointer to array of absolute toleracnes
      type(c_ptr), value :: abs_tol
      !> Relative integration tolerance
      real(kind=c_double), value :: rel_tol
      !> Maximum number of internal integration steps
      integer(kind=c_int), value :: max_steps
      !> Maximum number of convergence failures
      integer(kind=c_int), value :: max_conv_fails
    end subroutine solver_initialize
  
    !> Run the solver
    integer(kind=c_int) function solver_run(solver_data, state, env, &
                    t_initial, t_final) bind (c)
      use iso_c_binding
      !> Pointer to the initialized solver data
      type(c_ptr), value :: solver_data
      !> Pointer to the state array
      type(c_ptr), value :: state
      !> Pointer to the environmental state array
      type(c_ptr), value :: env
      !> Initial time (s)
      real(kind=c_double), value :: t_initial
      !> Final time (s)
      real(kind=c_double), value :: t_final
    end function solver_run

    !> Add condensed data to the solver data block
    subroutine rxn_add_condensed_data(rxn_type, n_int_param, &
                    n_float_param, int_param, float_param, solver_data) &
                    bind (c)
      use iso_c_binding
      !> Reaction type
      integer(kind=c_int), value :: rxn_type
      !> Number of integer parameters to add
      integer(kind=c_int), value :: n_int_param
      !> Number of floating-point parameters to add
      integer(kind=c_int), value :: n_float_param
      !> Pointer to the integer parameter array
      type(c_ptr), value :: int_param
      !> Pointer to the floating-point parameter array
      type(c_ptr), value :: float_param
      !> Pointer to the solver data
      type(c_ptr), value :: solver_data
    end subroutine rxn_add_condensed_data

    !> Update a photolysis rate constant
    subroutine rxn_set_photo_rate(photo_id, base_rate, solver_data) bind(c)
      use iso_c_binding
      !> Id used to identify reactions to update
      integer(kind=c_int), value :: photo_id
      !> Base (unscaled) photolysis rate
      real(kind=c_double), value :: base_rate
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine rxn_set_photo_rate

    !> Print the solver data
    subroutine rxn_print_data(solver_data) bind(c)
      use iso_c_binding
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine rxn_print_data

    !> Get the current reaction rates
    function rxn_get_rates(solver_data, state, env, n_rxn) bind(c)
      use iso_c_binding
      !> Reaction rates
      type(c_ptr) :: rxn_get_rates
      !> Solver data
      type(c_ptr), value :: solver_data
      !> Pointer to the state array
      type(c_ptr), value :: state
      !> Pointer to the environmental state array
      type(c_ptr), value :: env
      !> Number of reactions
      integer(kind=c_int) :: n_rxn
    end function rxn_get_rates

  end interface

  !> Solver data
  !!
  !! Acts as the interface between the phlex-chem module and the solver. 
  !! Instances of the type hold a pointer to a solver c object and provide
  !! functions to initialize the solver with model data and run the solver
  !! over a specified time step.
  type :: phlex_solver_data_t
    private
    !> C Solver object
    type(c_ptr), public :: solver_c_ptr
    !> Relative tolerance for the integration
    real(kind=dp), public :: rel_tol = PMC_SOLVER_DEFAULT_REL_TOL
    !> Maximum number of timesteps
    integer(kind=i_kind), public :: max_steps = PMC_SOLVER_DEFAULT_MAX_STEPS
    !> Maximum number of convergence failures
    integer(kind=i_kind), public :: max_conv_fails = PMC_SOLVER_DEFAULT_MAX_CONV_FAILS
  contains
    !> Initialize the solver
    procedure :: initialize
    !> Set a photolysis rate constant
    procedure :: set_photo_rate
    !> Integrate over a given time step
    procedure :: solve
    !> Checks whether a solver is available
    procedure :: is_solver_available
    !> Get the current reaction rates
    procedure :: get_rates
    !> Print the solver data
    procedure :: print => do_print
  end type phlex_solver_data_t

  ! Constructor for phlex_solver_data_t
  interface phlex_solver_data_t
    procedure :: constructor
  end interface phlex_solver_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for phlex_solver_data_t
  function constructor() result (new_obj)
    
    !> New solver variable
    type(phlex_solver_data_t), pointer :: new_obj

    !> Allocate space for the new object
    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the solver
  subroutine initialize(this, var_type, abs_tol, mechanisms, rxn_phase)

    !> Solver data
    class(phlex_solver_data_t), intent(inout) :: this
    !> Variable type for each species in the state array. This array must be
    !! of the same length as the state array.
    integer(kind=i_kind), pointer, intent(in) :: var_type(:)
    !> Absolute tolerance for each species in the state array. This array must
    !! be of the same length as the state array. Values for CONST and PSSA
    !! species will be ignored by the solver.
    real(kind=dp), allocatable, intent(in) :: abs_tol(:)
    !> Mechanisms to include in solver
    type(mechanism_data_t), pointer, intent(in) :: mechanisms(:)
    !> Reactions phase to solve -- gas, aerosol, or both (default)
    !! Use parameters in pmc_rxn_data to specify phase:
    !! GAS_RXN, AERO_RXN, GAS_AERO_RXN
    integer(kind=i_kind), intent(in) :: rxn_phase

    ! Variable types
    integer(kind=c_int), pointer :: var_type_c(:)
    ! Absolute tolerances
    real(kind=c_double), pointer :: abs_tol_c(:)
    ! Indices for iteration
    integer(kind=i_kind) :: i_mech, i_rxn
    ! Reaction pointer
    class(rxn_data_t), pointer :: rxn
    ! Reaction factory object for getting reaction type
    type(rxn_factory_t) :: rxn_factory
    ! Integer parameters being transfered
    integer(kind=c_int), pointer :: int_param(:)
    ! Floating point parameters being transfered
    real(kind=c_double), pointer :: float_param(:)
    ! Number of reactions
    integer(kind=c_int) :: n_rxn
    ! Number of integer reaction parameters
    integer(kind=c_int) :: n_int_param
    ! Number of floating-point reaction parameters
    integer(kind=c_int) :: n_float_param

    ! Make sure the variable type and absolute tolerance arrays are of 
    ! equal length
    call assert_msg(825843466, size(abs_tol).eq.size(var_type), &
            "Mismatched absolute tolerance and variable type arrays: "// &
            "abs_tol size: "//trim(to_string(size(abs_tol)))// &
            "; var_type: "//trim(to_string(size(var_type))))

    ! Set the absolute tolerance and variable type arrays
    allocate(var_type_c(size(var_type)))
    allocate(abs_tol_c(size(abs_tol)))
    var_type_c(:) = int(var_type(:), kind=c_int)
    abs_tol_c(:) = real(abs_tol(:), kind=c_double)

    ! Initialize the counters
    n_rxn = 0
    n_int_param = 0
    n_float_param = 0

    ! Calculate the number of reactions and the size of the condensed data
    do i_mech=1, size(mechanisms)
      do i_rxn=1, mechanisms(i_mech)%size()
        associate (rxn => mechanisms(i_mech)%get_rxn(i_rxn))
        select case (rxn%rxn_phase)
          case (GAS_RXN)
            if (rxn_phase.eq.AERO_RXN) cycle
          case (AERO_RXN)
            if (rxn_phase.eq.GAS_RXN) cycle
          case (GAS_AERO_RXN)
            if (rxn_phase.eq.GAS_RXN) cycle
        end select
        n_rxn = n_rxn + 1
        n_int_param = n_int_param + size(rxn%condensed_data_int)
        n_float_param = n_float_param + size(rxn%condensed_data_real)
        end associate
      end do
    end do

    ! Get a new solver object
    this%solver_c_ptr = solver_new( &
            int(size(var_type_c), kind=c_int),          & ! Size of the state variable
            c_loc(var_type_c),                          & ! Variable types
            n_rxn,                                      & ! Number of reactions
            n_int_param,                                & ! Number of integer parameters
            n_float_param                               & ! Number of floating-point parameters
            )

    ! Add all the condensed reaction data to the solver data block for reactions
    ! of the specified phase
    do i_mech=1, size(mechanisms)
      do i_rxn=1, mechanisms(i_mech)%size()
        
        ! Assign rxn to the current reaction
        associate (rxn => mechanisms(i_mech)%get_rxn(i_rxn))
        
        ! Check reaction phase
        select case (rxn%rxn_phase)
          case (GAS_RXN)
            if (rxn_phase.eq.AERO_RXN) cycle
          case (AERO_RXN)
            if (rxn_phase.eq.GAS_RXN) cycle
          case (GAS_AERO_RXN)
            if (rxn_phase.eq.GAS_RXN) cycle
        end select

        ! Load temporary data arrays
        allocate(int_param(size(rxn%condensed_data_int)))
        allocate(float_param(size(rxn%condensed_data_real)))
        int_param(:) = int(rxn%condensed_data_int(:), kind=c_int)
        float_param(:) = real(rxn%condensed_data_real(:), kind=c_double)
      
        ! Send the condensed data to the solver
        call rxn_add_condensed_data ( &
                rxn_factory%get_type(rxn),      & ! Reaction type
                size(int_param),                & ! Size of integer parameter array
                size(float_param),              & ! Size of floating-point parameter array
                c_loc(int_param),               & ! Pointer to integer parameter array
                c_loc(float_param),             & ! Pointer to floating-point parameter array
                this%solver_c_ptr               & ! Pointer to solver data
                )

        ! Deallocate temporary arrays
        deallocate(int_param)
        deallocate(float_param)

        end associate
      end do
    end do

    ! Initialize the solver
    call solver_initialize( &
            this%solver_c_ptr,                          & ! Pointer to solver data
            c_loc(abs_tol_c),                           & ! Absolute tolerances
            real(this%rel_tol, kind=c_double),          & ! Relative tolerance
            int(this%max_steps, kind=c_int),            & ! Maximum number of integration steps
            int(this%max_conv_fails, kind=c_int)        & ! Maximum number of convergence failures
            )

    ! Free memory allocated for solver initialization
    deallocate(abs_tol_c)
    deallocate(var_type_c)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set a base (unscaled) photolysis rate constant based on the index
  !! provided by an external phototlysis module during intialization.
  subroutine set_photo_rate(this, photo_id, base_rate)

    !> Solver data
    class(phlex_solver_data_t), intent(inout) :: this
    !> Id to use to identify reactions to update
    integer(kind=i_kind), intent(in) :: photo_id
    !> New base (unscaled) photolysis rate
    real(kind=dp), intent(in) :: base_rate

    call rxn_set_photo_rate( &
            int(photo_id, kind=c_int),          & ! Id of photo rate to update
            real(base_rate, kind=c_double),     & ! New value for the base photo rate constant
            this%solver_c_ptr                   & ! Pointer to solver data
            )

  end subroutine set_photo_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solve the mechanism(s) for a specified timestep
  subroutine solve(this, phlex_state, t_initial, t_final)

    !> Solver data
    class(phlex_solver_data_t), intent(inout) :: this
    !> Model state
    type(phlex_state_t), target, intent(inout) :: phlex_state
    !> Start time (s)
    real(kind=dp), intent(in) :: t_initial
    !> End time (s)
    real(kind=dp), intent(in) :: t_final

    integer(kind=c_int) :: solver_status
    
    ! Run the solver
    solver_status = solver_run( &
            this%solver_c_ptr,                  & ! Pointer to intialized solver
            c_loc(phlex_state%state_var),       & ! Pointer to state array
            c_loc(phlex_state%env_var),         & ! Pointer to environmental variables
            real(t_initial, kind=c_double),     & ! Start time (s)
            real(t_final, kind=c_double)        & ! Final time (s)
            )

    call assert_msg(997420005, solver_status.eq.0, "Solver failed");

  end subroutine solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether a solver is available for the integration 
  logical function is_solver_available(this)

    !> Solver data
    class(phlex_solver_data_t), intent(in) :: this

#ifdef PMC_USE_SUNDIALS
    is_solver_available = .true.
#else
    is_solver_available = .false.
#endif

  end function is_solver_available

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the current reaction rates
  function get_rates(this, phlex_state)

    !> Reaction rates
    real(kind=dp), allocatable :: get_rates(:)
    !> Solver data
    class(phlex_solver_data_t), intent(in) :: this
    !> Model state
    type(phlex_state_t), pointer, intent(in) :: phlex_state

    real(kind=c_double), pointer :: rates(:)
    type(c_ptr) :: rates_c_ptr
    integer(kind=c_int) :: n_rxn

    rates_c_ptr = rxn_get_rates(                &
            this%solver_c_ptr,                  & ! Pointer to solver data
            c_loc(phlex_state%state_var),       & ! Pointer to state array
            c_loc(phlex_state%env_var),         & ! Pointer to environmental variables
            n_rxn                               & ! Number of reactions
            )
    call c_f_pointer(rates_c_ptr, rates, [n_rxn])
    allocate(get_rates(n_rxn))
    get_rates(:) = real(rates(:), kind=dp)
    deallocate(rates)

  end function get_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the solver data
  subroutine do_print(this)

    !> Solver data
    class(phlex_solver_data_t), intent(in) :: this

    call rxn_print_data(this%solver_c_ptr)

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef PHLEX_SOLVER_SUCCESS
#undef PHLEX_SOLVER_FAIL
end module pmc_phlex_solver_data
