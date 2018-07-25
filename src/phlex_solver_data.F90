! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_phlex_solver_data module.

!> The phlex_solver_data_t structure and associated subroutines.
module pmc_phlex_solver_data
#define PHLEX_SOLVER_SUCCESS 0
#define PHLEX_SOLVER_FAIL 1

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_constants,                   only : i_kind, dp
  use pmc_mechanism_data
  use pmc_phlex_state
  use pmc_rxn_data
  use pmc_rxn_factory
  use pmc_sub_model_data
  use pmc_sub_model_factory
  use pmc_util,                        only : assert_msg, to_string, die_msg

  use iso_c_binding

  implicit none
  private

  public :: phlex_solver_data_t

  !> Default relative tolerance for integration
  real(kind=dp), parameter :: PMC_SOLVER_DEFAULT_REL_TOL = 1.0D-8
  !> Default max number of integration steps
  integer(kind=i_kind), parameter :: PMC_SOLVER_DEFAULT_MAX_STEPS = 10000000
  !> Default maximum number of integration convergence failures
  integer(kind=i_kind), parameter :: PMC_SOLVER_DEFAULT_MAX_CONV_FAILS = 1000000

  !> Result code indicating successful completion
  integer, parameter :: PMC_SOLVER_SUCCESS = 0

  !> Interface to c ODE solver functions
  interface
    !> Get a new solver 
    type(c_ptr) function solver_new(n_state_var, var_type, &
                    n_rxn, n_rxn_int_param, n_rxn_float_param, &
                    n_aero_phase, n_aero_phase_int_param, &
                    n_aero_phase_float_param, n_aero_rep, &
                    n_aero_rep_int_param, n_aero_rep_float_param, &
                    n_sub_model, n_sub_model_int_param, &
                    n_sub_model_float_param) bind (c)
      use iso_c_binding
      !> Number of variables on the state array (including const, PSSA, etc.)
      integer(kind=c_int), value :: n_state_var
      !> Pointer to array of state variable types (solver, constant, PSSA)
      type(c_ptr), value :: var_type
      !> Number of reactions to solve
      integer(kind=c_int), value :: n_rxn
      !> Total number of integer parameters for all reactions
      integer(kind=c_int), value :: n_rxn_int_param
      !> Total number of floating-point parameters for all reactions
      integer(kind=c_int), value :: n_rxn_float_param
      !> Number of aerosol phases 
      integer(kind=c_int), value :: n_aero_phase
      !> Total number of integer parameters for all aerosol phases
      integer(kind=c_int), value :: n_aero_phase_int_param
      !> Total number of floating-point parameters for all aerosol phases
      integer(kind=c_int), value :: n_aero_phase_float_param
      !> Number of aerosol representations
      integer(kind=c_int), value :: n_aero_rep
      !> Total number of integer parameters for all aerosol representations
      integer(kind=c_int), value :: n_aero_rep_int_param
      !> Total number of floating-point parameters for all aerosol
      !! representations
      integer(kind=c_int), value :: n_aero_rep_float_param
      !> Number of sub models
      integer(kind=c_int), value :: n_sub_model
      !> Total number of integer parameters for all sub models
      integer(kind=c_int), value :: n_sub_model_int_param
      !> Total number of floating-point parameters for all sub models
      integer(kind=c_int), value :: n_sub_model_float_param
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

    !> Add condensed reaction data to the solver data block
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

    !> Update reaction data
    subroutine rxn_update_data(rxn_type, update_data, solver_data) bind(c)
      use iso_c_binding
      !> Reaction type to updateto update
      integer(kind=c_int), value :: rxn_type
      !> Data required by reaction for updates
      type(c_ptr), value :: update_data
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine rxn_update_data

    !> Print the solver data
    subroutine rxn_print_data(solver_data) bind(c)
      use iso_c_binding
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine rxn_print_data

    !> Add condensed sub model data to the solver data block
    subroutine sub_model_add_condensed_data(sub_model_type, n_int_param, &
                  n_float_param, int_param, float_param, solver_data) bind(c)
      use iso_c_binding
      !> Sub model  type
      integer(kind=c_int), value :: sub_model_type
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
    end subroutine sub_model_add_condensed_data

    !> Update reaction data
    subroutine sub_model_update_data(sub_model_type, update_data, &
              solver_data) bind(c)
      use iso_c_binding
      !> Reaction type to updateto update
      integer(kind=c_int), value :: sub_model_type
      !> Data required by reaction for updates
      type(c_ptr), value :: update_data
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine sub_model_update_data

    !> Print the solver data
    subroutine sub_model_print_data(solver_data) bind(c)
      use iso_c_binding
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine sub_model_print_data

    !> Get a sub model parameter id
    function sub_model_get_parameter_id_sd(solver_data, sub_model_type, &
                  identifiers) bind (c)
      use iso_c_binding
      !> Parameter id
      integer(kind=c_int) :: sub_model_get_parameter_id_sd
      !> Solver data
      type(c_ptr), value :: solver_data
      !> Sub model type
      integer(kind=c_int), value :: sub_model_type
      !> Sub model identifiers
      type(c_ptr), value :: identifiers
    end function sub_model_get_parameter_id_sd

    !> Get a sub model parameter value
    function sub_model_get_parameter_value_sd(solver_data, parameter_id) &
                  bind (c)
      use iso_c_binding
      !> Parameter value
      real(kind=c_double) :: sub_model_get_parameter_value_sd
      !> Solver data
      type(c_ptr), value :: solver_data
      !> Parameter id
      integer(kind=c_int), value :: parameter_id
    end function sub_model_get_parameter_value_sd

    !> Add condensed aerosol phase data to the solver data block
    subroutine aero_phase_add_condensed_data(n_int_param, n_float_param, &
                  int_param, float_param, solver_data) bind(c)
      use iso_c_binding
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
    end subroutine aero_phase_add_condensed_data

    !> Print the solver data
    subroutine aero_phase_print_data(solver_data) bind(c)
      use iso_c_binding
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine aero_phase_print_data

    !> Add condensed aerosol representation data to the solver data block
    subroutine aero_rep_add_condensed_data(aero_rep_type, n_int_param, &
                  n_float_param, int_param, float_param, solver_data) bind(c)
      use iso_c_binding
      !> Aerosol representation type
      integer(kind=c_int), value :: aero_rep_type
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
    end subroutine aero_rep_add_condensed_data

    !> Update aerosol representation data
    subroutine aero_rep_update_data(aero_rep_type, update_data, solver_data) &
              bind(c)
      use iso_c_binding
      !> Aerosol representation type to updateto update
      integer(kind=c_int), value :: aero_rep_type
      !> Data required by aerosol representation for updates
      type(c_ptr), value :: update_data
      !> Solver data
      type(c_ptr), value :: solver_data
    end subroutine aero_rep_update_data

    !> Print the aerosol representation data
    subroutine aero_rep_print_data(solver_data) bind(c)
      use iso_c_binding
      !> Pointer to the solver data
      type(c_ptr), value :: solver_data
    end subroutine aero_rep_print_data

    !> Free the memory associated with a solver
    pure subroutine solver_free(solver_data) bind(c)
      use iso_c_binding
      !> Pointer to the solver data
      type(c_ptr), value, intent(in) :: solver_data
    end subroutine solver_free

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
    integer(kind=i_kind), public :: max_conv_fails = &
            PMC_SOLVER_DEFAULT_MAX_CONV_FAILS
    !> Flag indicating whether the solver was intialized
    logical :: initialized = .false.
  contains
    !> Initialize the solver
    procedure :: initialize
    !> Update sub-model data
    procedure :: update_sub_model_data
    !> Update reactions data
    procedure :: update_rxn_data
    !> Update aerosol representation data
    procedure :: update_aero_rep_data
    !> Integrate over a given time step
    procedure :: solve
    !> Checks whether a solver is available
    procedure :: is_solver_available
    !> Get a parameter id
    procedure :: get_sub_model_parameter_id
    !> Get a current parameter value
    procedure :: get_sub_model_parameter_value
    !> Print the solver data
    procedure :: print => do_print
    !> Finalize the solver data
    final :: finalize
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
  subroutine initialize(this, var_type, abs_tol, mechanisms, aero_phases, &
                  aero_reps, sub_models, rxn_phase)

    !> Solver data
    class(phlex_solver_data_t), intent(inout) :: this
    !> Variable type for each species in the state array. This array must be
    !! of the same length as the state array.
    integer(kind=i_kind), allocatable, intent(in) :: var_type(:)
    !> Absolute tolerance for each species in the state array. This array must
    !! be of the same length as the state array. Values for CONST and PSSA
    !! species will be ignored by the solver.
    real(kind=dp), allocatable, intent(in) :: abs_tol(:)
    !> Mechanisms to include in solver
    type(mechanism_data_ptr), pointer, intent(in) :: mechanisms(:)
    !> Aerosol phases to include
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phases(:)
    !> Aerosol representations to include
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_reps(:)
    !> Sub models to include
    type(sub_model_data_ptr), pointer, intent(in) :: sub_models(:)
    !> Reactions phase to solve -- gas, aerosol, or both (default)
    !! Use parameters in pmc_rxn_data to specify phase:
    !! GAS_RXN, AERO_RXN, GAS_AERO_RXN
    integer(kind=i_kind), intent(in) :: rxn_phase

    ! Variable types
    integer(kind=c_int), pointer :: var_type_c(:)
    ! Absolute tolerances
    real(kind=c_double), pointer :: abs_tol_c(:)
    ! Indices for iteration
    integer(kind=i_kind) :: i_mech, i_rxn, i_aero_phase, i_aero_rep, &
            i_sub_model
    ! Reaction pointer
    class(rxn_data_t), pointer :: rxn
    ! Reaction factory object for getting reaction type
    type(rxn_factory_t) :: rxn_factory
    ! Aerosol phase pointer
    type(aero_phase_data_t), pointer :: aero_phase
    ! Aerosol representation pointer
    class(aero_rep_data_t), pointer :: aero_rep
    ! Aerosol representation factory object for getting aerosol 
    ! representation type
    type(aero_rep_factory_t) :: aero_rep_factory
    ! Sub model pointer
    class(sub_model_data_t), pointer :: sub_model
    ! Sub model factory for getting sub model type
    type(sub_model_factory_t) :: sub_model_factory
    ! Integer parameters being transfered
    integer(kind=c_int), pointer :: int_param(:)
    ! Floating point parameters being transfered
    real(kind=c_double), pointer :: float_param(:)
    ! Number of reactions
    integer(kind=c_int) :: n_rxn
    ! Number of integer reaction parameters
    integer(kind=c_int) :: n_rxn_int_param
    ! Number of floating-point reaction parameters
    integer(kind=c_int) :: n_rxn_float_param
    ! Number of aerosol phases
    integer(kind=c_int) :: n_aero_phase
    ! Number of integer aerosol phase  parameters
    integer(kind=c_int) :: n_aero_phase_int_param
    ! Number of floating-point aerosol phase parameters
    integer(kind=c_int) :: n_aero_phase_float_param
    ! Number of aerosol representations
    integer(kind=c_int) :: n_aero_rep
    ! Number of integer aerosol representation parameters
    integer(kind=c_int) :: n_aero_rep_int_param
    ! Number of floating-point aerosol representation parameters
    integer(kind=c_int) :: n_aero_rep_float_param
    ! Number of sub models
    integer(kind=c_int) :: n_sub_model
    ! Number of integer sub model parameters
    integer(kind=c_int) :: n_sub_model_int_param
    ! Number of floating-point sub model parameters
    integer(kind=c_int) :: n_sub_model_float_param

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
    n_rxn_int_param = 0
    n_rxn_float_param = 0

    ! Calculate the number of reactions and the size of the condensed data
    do i_mech=1, size(mechanisms)
      do i_rxn=1, mechanisms(i_mech)%val%size()
        rxn => mechanisms(i_mech)%val%get_rxn(i_rxn)
        select case (rxn%rxn_phase)
          case (GAS_RXN)
            if (rxn_phase.eq.AERO_RXN) cycle
          case (AERO_RXN)
            if (rxn_phase.eq.GAS_RXN) cycle
          case (GAS_AERO_RXN)
            if (rxn_phase.eq.GAS_RXN) cycle
        end select
        n_rxn = n_rxn + 1
        n_rxn_int_param = n_rxn_int_param + size(rxn%condensed_data_int)
        n_rxn_float_param = n_rxn_float_param + size(rxn%condensed_data_real)
      end do
    end do
    rxn => null()

    ! Initialize the aerosol phase counters
    n_aero_phase = size(aero_phases)
    n_aero_phase_int_param = 0
    n_aero_phase_float_param = 0

    ! Calculate the size of the aerosol phases condensed data
    do i_aero_phase=1, n_aero_phase
      aero_phase => aero_phases(i_aero_phase)%val
      n_aero_phase_int_param = n_aero_phase_int_param + &
              size(aero_phase%condensed_data_int)
      n_aero_phase_float_param = n_aero_phase_float_param + &
              size(aero_phase%condensed_data_real)
    end do
    aero_phase => null()

    ! Initialize the aerosol representation counters
    n_aero_rep = size(aero_reps)
    n_aero_rep_int_param = 0
    n_aero_rep_float_param = 0

    ! Calculate the size of the aerosol representations condensed data
    do i_aero_rep=1, n_aero_rep
      aero_rep => aero_reps(i_aero_rep)%val
      n_aero_rep_int_param = n_aero_rep_int_param + &
              size(aero_rep%condensed_data_int)
      n_aero_rep_float_param = n_aero_rep_float_param + &
              size(aero_rep%condensed_data_real)
    end do
    aero_rep => null()

    ! Initialize the sub model counters
    n_sub_model = size(sub_models)
    n_sub_model_int_param = 0
    n_sub_model_float_param = 0

    ! Calculate the size of the sub model condensed data
    do i_sub_model=1, n_sub_model
      sub_model => sub_models(i_sub_model)%val
      n_sub_model_int_param = n_sub_model_int_param + &
              size(sub_model%condensed_data_int)
      n_sub_model_float_param = n_sub_model_float_param + &
              size(sub_model%condensed_data_real)
    end do
    sub_model => null()

    ! Get a new solver object
    this%solver_c_ptr = solver_new( &
            int(size(var_type_c), kind=c_int), & ! Size of the state variable
            c_loc(var_type_c),                 & ! Variable types
            n_rxn,                             & ! # of reactions
            n_rxn_int_param,                   & ! # of rxn data int params
            n_rxn_float_param,                 & ! # of rxn data real params
            n_aero_phase,                      & ! # of aero phases
            n_aero_phase_int_param,            & ! # of aero phase int params
            n_aero_phase_float_param,          & ! # of aero phase real params
            n_aero_rep,                        & ! # of aero reps
            n_aero_rep_int_param,              & ! # of aero rep int params
            n_aero_rep_float_param,            & ! # of aero rep real params
            n_sub_model,                       & ! # of sub models
            n_sub_model_int_param,             & ! # of sub model int params
            n_sub_model_float_param            & ! # of sub model real params
            )

    ! Add all the condensed reaction data to the solver data block for
    ! reactions of the specified phase
    do i_mech=1, size(mechanisms)
      do i_rxn=1, mechanisms(i_mech)%val%size()
       
        ! FIXME Put ZSR aerosol water first, so water is available for other
        ! reactions - and find a better way to account for inter-dependence 
        ! of reactions/sub-models

        ! Assign rxn to the current reaction
        rxn => mechanisms(i_mech)%val%get_rxn(i_rxn)
        
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
                int(rxn_factory%get_type(rxn), kind=c_int),& ! Rxn type
                int(size(int_param), kind=c_int),          & ! Int array size
                int(size(float_param), kind=c_int),        & ! Real array size
                c_loc(int_param),                          & ! Int array ptr
                c_loc(float_param),                        & ! Real array ptr
                this%solver_c_ptr                          & ! Solver data ptr
                )

        ! Deallocate temporary arrays
        deallocate(int_param)
        deallocate(float_param)

      end do
    end do
    rxn => null()

    ! Add all the condensed aerosol phase data to the solver data block
    do i_aero_phase=1, size(aero_phases)

      ! Assign aero_phase to the current aerosol phase
      aero_phase => aero_phases(i_aero_phase)%val

      ! Load temporary data arrays
      allocate(int_param(size(aero_phase%condensed_data_int)))
      allocate(float_param(size(aero_phase%condensed_data_real)))
      int_param(:) = int(aero_phase%condensed_data_int(:), kind=c_int)
      float_param(:) = real(aero_phase%condensed_data_real(:), kind=c_double)

      ! Send the condensed data to the solver
      call aero_phase_add_condensed_data ( &
              int(size(int_param), kind=c_int),   & ! Int array size
              int(size(float_param), kind=c_int), & ! Real array size
              c_loc(int_param),                   & ! Int array ptr
              c_loc(float_param),                 & ! Real array ptr
              this%solver_c_ptr                   & ! Solver data ptr
              )

      ! Deallocate temporary arrays
      deallocate(int_param)
      deallocate(float_param)

    end do
    aero_phase => null()

    ! Add all the condensed aerosol representation data to the solver data
    ! block
    do i_aero_rep=1, size(aero_reps)

      ! Assign aero_rep to the current aerosol representation
      aero_rep => aero_reps(i_aero_rep)%val

      ! Load temporary data arrays
      allocate(int_param(size(aero_rep%condensed_data_int)))
      allocate(float_param(size(aero_rep%condensed_data_real)))
      int_param(:) = int(aero_rep%condensed_data_int(:), kind=c_int)
      float_param(:) = real(aero_rep%condensed_data_real(:), kind=c_double)

      ! Send the condensed data to the solver
      call aero_rep_add_condensed_data ( &
              int(aero_rep_factory%get_type(aero_rep), kind=c_int), & 
                                                    ! Aero rep type
              int(size(int_param), kind=c_int),   & ! Int array size
              int(size(float_param), kind=c_int), & ! Real array size
              c_loc(int_param),                   & ! Int array ptr
              c_loc(float_param),                 & ! Real array ptr
              this%solver_c_ptr                   & ! Solver data ptr
              )

      ! Deallocate temporary arrays
      deallocate(int_param)
      deallocate(float_param)

    end do
    aero_rep => null()

    ! Add all the condensed sub model data to the solver data block
    do i_sub_model=1, size(sub_models)

      ! Assign sub_model to the current sub model
      sub_model => sub_models(i_sub_model)%val

      ! Load temporary data arrays
      allocate(int_param(size(sub_model%condensed_data_int)))
      allocate(float_param(size(sub_model%condensed_data_real)))
      int_param(:) = int(sub_model%condensed_data_int(:), kind=c_int)
      float_param(:) = real(sub_model%condensed_data_real(:), kind=c_double)

      ! Send the condensed data to the solver
      call sub_model_add_condensed_data ( &
              int(sub_model_factory%get_type(sub_model), kind=c_int), & 
                                                    ! Sub model type
              int(size(int_param), kind=c_int),   & ! Int array size
              int(size(float_param), kind=c_int), & ! Real array size
              c_loc(int_param),                   & ! Int array ptr
              c_loc(float_param),                 & ! Real array ptr
              this%solver_c_ptr                   & ! Solver data ptr
              )

      ! Deallocate temporary arrays
      deallocate(int_param)
      deallocate(float_param)

    end do
    sub_model => null()

    ! Initialize the solver
    call solver_initialize( &
            this%solver_c_ptr,                  & ! Pointer to solver data
            c_loc(abs_tol_c),                   & ! Absolute tolerances
            real(this%rel_tol, kind=c_double),  & ! Relative tolerance
            int(this%max_steps, kind=c_int),    & ! Max # of integration steps
            int(this%max_conv_fails, kind=c_int)& ! Max # of convergence fails
            )

    ! Flag the solver as initialized
    this%initialized = .true.

    ! Free memory allocated for solver initialization
    deallocate(abs_tol_c)
    deallocate(var_type_c)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update sub-model data
  subroutine update_sub_model_data(this, update_data)

    !> Solver data
    class(phlex_solver_data_t), intent(inout) :: this
    !> Update data
    class(sub_model_update_data_t), intent(in) :: update_data

    call sub_model_update_data( &
            update_data%get_type(),     & ! Sub-model type to update
            update_data%get_data(),     & ! Data needed to perform update
            this%solver_c_ptr           & ! Pointer to solver data
            )

  end subroutine update_sub_model_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update reaction data
  subroutine update_rxn_data(this, update_data)

    !> Solver data
    class(phlex_solver_data_t), intent(inout) :: this
    !> Update data
    class(rxn_update_data_t), intent(in) :: update_data

    call rxn_update_data( &
            update_data%get_type(),     & ! Reaction type to update
            update_data%get_data(),     & ! Data needed to perform update 
            this%solver_c_ptr           & ! Pointer to solver data
            )

  end subroutine update_rxn_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update aerosol representation data based on data passed from the host
  !! model related to aerosol properties
  subroutine update_aero_rep_data(this, update_data)

    !> Solver data
    class(phlex_solver_data_t), intent(inout) :: this
    !> Update data
    class(aero_rep_update_data_t), intent(in) :: update_data

    call aero_rep_update_data( &
            update_data%get_type(),     & ! Aerosol representation type
            update_data%get_data(),     & ! Data needed to perform update
            this%solver_c_ptr           & ! Pointer to solver data
            )

  end subroutine update_aero_rep_data

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
            this%solver_c_ptr,              & ! Pointer to intialized solver
            c_loc(phlex_state%state_var),   & ! Pointer to state array
            c_loc(phlex_state%env_var),     & ! Pointer to environmental vars
            real(t_initial, kind=c_double), & ! Start time (s)
            real(t_final, kind=c_double)    & ! Final time (s)
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

  !> Get a sub model parameter id
  !!
  !! Returns the id for use with sub_model_get_parameter_value to get a
  !! current sub-model parameter from the solver data
  function get_sub_model_parameter_id(this, sub_model_type, identifiers) &
      result (parameter_id)

    !> The parameter id
    integer(kind=c_int) :: parameter_id
    !> Solver data
    class(phlex_solver_data_t), intent(in) :: this
    !> Sub-model type
    integer(kind=i_kind), intent(in) :: sub_model_type
    !> Identifiers required by the sub-model
    type(c_ptr), intent(in) :: identifiers

    integer(kind=c_int) :: sub_model_type_c

    sub_model_type_c = int(sub_model_type, kind=c_int)

    parameter_id = sub_model_get_parameter_id_sd( &
            this%solver_c_ptr,          & ! Pointer to solver data
            sub_model_type_c,           & ! Sub model type
            identifiers                 & ! Indentifiers needed
            )

  end function get_sub_model_parameter_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a sub model parameter
  !!
  !! Returns the value associate with a sub model output, based on the current
  !! solver state.
  function get_sub_model_parameter_value(this, parameter_id) &
      result (parameter_value)

    !> Value of the parameter
    real(kind=dp) :: parameter_value
    !> Solver data
    class(phlex_solver_data_t), intent(in) :: this
    !> Parameter id
    integer(kind=c_int), intent(in) :: parameter_id

    real(kind=c_double) :: parameter_value_c

    parameter_value_c = sub_model_get_parameter_value_sd( &
            this%solver_c_ptr,          & ! Pointer to solver data
            parameter_id                & ! Id of the parameter
            )

    parameter_value = real(parameter_value_c, kind=dp)

  end function get_sub_model_parameter_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the solver data
  subroutine do_print(this)

    !> Solver data
    class(phlex_solver_data_t), intent(in) :: this

    call rxn_print_data(this%solver_c_ptr)
    call aero_phase_print_data(this%solver_c_ptr)
    call aero_rep_print_data(this%solver_c_ptr)
    call sub_model_print_data(this%solver_c_ptr)

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the solver data
  elemental subroutine finalize(this)

    !> Solver data
    type(phlex_solver_data_t), intent(inout) :: this

    if (this%initialized) call solver_free(this%solver_c_ptr)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef PHLEX_SOLVER_SUCCESS
#undef PHLEX_SOLVER_FAIL
end module pmc_phlex_solver_data
