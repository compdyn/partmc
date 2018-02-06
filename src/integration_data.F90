! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_integration_data module.

!> The integration_data_t structure and associated subroutines.
module pmc_integration_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_util,                       only : assert_msg, to_string, die_msg

  use iso_c_binding

  implicit none
  private

  !> Default relative tolerance for integration
  real(kind=dp), parameter :: PMC_INTEGRATION_DEFAULT_REL_TOL = 1.0D-5
  !> Default max number of integration steps
  integer(kind=i_kind), parameter :: PMC_INTEGRATION_DEFAULT_MAX_STEPS = 10000000

  !> Result code indicating successful completion.
  integer, parameter :: PMC_INTEGRATION_SUCCESS             = 0
  !> Result code indicating no available integration routine
  integer, parameter :: PMC_INTEGRATION_NO_AVAIL_SOLVER     = 12
#ifdef PMC_USE_SUNDIALS
  !> Result code indicating failure to allocate \c y vector.
  integer, parameter :: PMC_INTEGRATION_INIT_Y              = 1
  !> Result code indicating failure to allocate \c abstol vector.
  integer, parameter :: PMC_INTEGRATION_INIT_ABSTOL         = 2
  !> Result code indicating failure to create the solver.
  integer, parameter :: PMC_INTEGRATION_INIT_CVODE_MEM      = 3
  !> Result code indicating failure to initialize the solver.
  integer, parameter :: PMC_INTEGRATION_INIT_CVODE          = 4
  !> Result code indicating failure to set tolerances.
  integer, parameter :: PMC_INTEGRATION_SVTOL               = 5
  !> Result code indicating failure to set maximum steps.
  integer, parameter :: PMC_INTEGRATION_SET_MAX_STEPS       = 6
  !> Result code indicating failure of the solver.
  integer, parameter :: PMC_INTEGRATION_FAIL                = 7
  !> Result code indicating failure to set dense Jacobian solver
  integer, parameter :: PMC_INTEGRATION_DENSE_JAC           = 8
  !> Result code indicating failure to set Jacobian function
  integer, parameter :: PMC_INTEGRATION_JAC_FUNC            = 9
  !> Result code indicating failure to set user data
  integer, parameter :: PMC_INTEGRATION_SET_USER_DATA       = 10
  !> Result code indicating SUNDIALS realtype is not set to double precision
  integer, parameter :: PMC_INTEGRATION_WRONG_PRECISION     = 11
  !> Result code indicating failure to get the dense linear solver
  integer, parameter :: PMC_INTEGRATION_DENSE_LINEAR_SOLVER = 13
  !> Result code indicating failure to set the dense linear solver
  integer, parameter :: PMC_INTEGRATION_SET_LINEAR_SOLVER   = 14
  !> Result code indicating failure to set the maximum number of convergence failures
  integer, parameter :: PMC_INTEGRATION_SET_MAX_CONV_FAILS  = 15
  !> Result code indicating failure to set the SPGMR solver
  integer, parameter :: PMC_INTEGRATION_SPGMR_LINEAR_SOLVER = 16
  !> Result code indicating failure to set the preconditioner functions
  integer, parameter :: PMC_INTEGRATION_SET_PRECONDITIONER  = 17

  !> Interface to c ODE solver function
  interface
    integer(kind=c_int) function integrate_solver(neq, x, abs_tol, rel_tol, &
                                 max_steps, t_initial, t_final, sysdata) bind(c)
      use iso_c_binding
      !> Number of equations (i.e. species)
      integer(kind=c_int), value :: neq
      !> State array
      type(c_ptr), value :: x
      !> Absolute tolerance for each species
      type(c_ptr), value :: abs_tol
      !> Relative tolerance for all species
      real(kind=c_double), value :: rel_tol
      !> Maximum number of iterations for the solver
      integer(kind=c_int), value :: max_steps
      !> Initial time (s)
      real(kind=c_double), value :: t_initial
      !> Final time (s)
      real(kind=c_double), value :: t_final
      !> Pointer to model data required by integration functions
      type(c_ptr), value :: sysdata
    end function integrate_solver
  end interface
#endif

  !> Integration data
  !!
  !! Contains pointers to model data, state, and derivative functions, and 
  !! provides an interface to the solver(s) for integration. 
  type :: integration_data_t
    private
    !> PartMC model data
    type(c_ptr) :: phlex_core_c_ptr
    !> PartMC model state
    type(c_ptr) :: phlex_state_c_ptr
    !> Time derivative function pointer
    procedure(integration_data_deriv_func), pointer, nopass :: deriv_func_ptr => null()
    !> Jacobian matrix function pointer
    procedure(integration_data_jac_func), pointer, nopass :: jac_func_ptr => null()
    !> Absolute tolerance for each species
    real(kind=c_double), pointer :: abs_tol_c(:)
    !> Relative tolerance for the integration
    real(kind=dp), public :: rel_tol = PMC_INTEGRATION_DEFAULT_REL_TOL
    !> Maximum number of timesteps
    integer(kind=i_kind), public :: max_steps = PMC_INTEGRATION_DEFAULT_MAX_STEPS
  contains
    !> Integrate over a given time step
    procedure :: solve
    !> Checks the result status code from an integration
    procedure :: check_status
    !> Checks whether a solver is available
    procedure :: is_solver_available
    !> Gets the absolute integration tolerance for a species by index
    procedure :: get_abs_tol
    !> Print the integration data
    procedure :: print => do_print
  end type integration_data_t

  ! Constructor for integration_data_t
  interface integration_data_t
    procedure :: constructor, constructor_empty
  end interface integration_data_t

  abstract interface 
    !> Interface for the time derivative function
    subroutine integration_data_deriv_func(curr_time, deriv, &
                    phlex_core_c_ptr, phlex_state_c_ptr)
      use pmc_constants,              only: dp
      use iso_c_binding
      
      !> Current solver time
      real(kind=dp), intent(in) :: curr_time
      !> Time derivative to calculate
      real(kind=dp), pointer, intent(inout) :: deriv(:)
      !> Pointer to model data
      type(c_ptr), intent(in) :: phlex_core_c_ptr
      !> Pointer to model state
      type(c_ptr), intent(in) :: phlex_state_c_ptr
    end subroutine integration_data_deriv_func

    !> Interface for the Jacobian function
    subroutine integration_data_jac_func(curr_time, jac, &
                    phlex_core_c_ptr, phlex_state_c_ptr)
      use pmc_constants,              only: dp
      use iso_c_binding

      !> Current solver time
      real(kind=dp), intent(in) :: curr_time
      !> Jacobian matrix to calculate
      real(kind=dp), pointer, intent(inout) :: jac(:,:)
      !> Pointer to model data
      type(c_ptr), intent(in) :: phlex_core_c_ptr
      !> Pointer to model state
      type(c_ptr), intent(in) :: phlex_state_c_ptr
    end subroutine integration_data_jac_func
  end interface 

  public :: integration_data_t, integration_data_deriv_func, &
          integration_data_jac_func

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for integration_data_t
  function constructor(phlex_core_c_ptr, deriv_func_ptr, &
                  jac_func_ptr, abs_tol) result(new_obj)

    !> New integration variable
    type(integration_data_t), pointer :: new_obj
    !> Model data pointer
    type(c_ptr), intent(in) :: phlex_core_c_ptr
    !> Time derivative function pointer
    procedure(integration_data_deriv_func), pointer, intent(in) :: deriv_func_ptr
    !> Jacobian matrix function pointer
    procedure(integration_data_jac_func), pointer, intent(in) :: jac_func_ptr
    !> Absolute tolerance for each species
    real(kind=dp), pointer, intent(in) :: abs_tol(:)

    allocate(new_obj)
    call assert_msg(668626727, new_obj%is_solver_available(), &
            "No solver available for integration.")
    new_obj%phlex_core_c_ptr = phlex_core_c_ptr
    new_obj%deriv_func_ptr => deriv_func_ptr
    new_obj%jac_func_ptr => jac_func_ptr
    allocate(new_obj%abs_tol_c(size(abs_tol)))
    new_obj%abs_tol_c(:) = real(abs_tol(:), kind=c_double)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> General constructor for integration_data_t
  function constructor_empty() result(new_obj)

    !> New integration variable
    type(integration_data_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor_empty

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the integration over a specified time step
  function solve(this, state_array, phlex_state_c_ptr, &
                  time_step) result(solver_status)

    !> Solver status
    integer(kind=i_kind) :: solver_status
    !> Integration data
    class(integration_data_t), target, intent(inout) :: this
    !> Model state array. This should be a pointer to a 1-D array of species 
    !! concentrations (y) for use by the solver. It should hold the initial 
    !! (time t) concentrations for all species. After the integration, it will
    !! hold the final (t+dt) concentrations. It should be located within the
    !! phlex_state_c_ptr, so that during derivative function calls, the 
    !! current model state is correct.
    real(kind=dp), pointer, intent(inout) :: state_array(:)
    !> Model state pointer. This is a void pointer for use by the derivative
    !! functions to allow access to state variables including those that are 
    !! not solved for during the integration (e.g., temperature, pressure)
    type(c_ptr), intent(in) :: phlex_state_c_ptr
    !> Time step (s)
    real(kind=dp), intent(in) :: time_step

    ! Variables required by the c ODE solver
    !
    ! Number of equations
    integer(kind=c_int) :: neq_c
    ! Relative tolerance for each species
    real(kind=c_double) :: rel_tol_c
    ! Maximum number of integration steps
    integer(kind=c_int) :: max_steps_c
    ! Initial time (s)
    real(kind=c_double) :: t_initial_c
    ! Final time (s)
    real(kind=c_double) :: t_final_c
    ! Non-polymorphic pointer for fortran
    type(integration_data_t), pointer :: this_ptr => null()

    this_ptr => this

    ! Set up solver variables
    neq_c = int(size(state_array), kind=c_int)
    rel_tol_c = real(this%rel_tol, kind=c_double)
    max_steps_c = int(this%max_steps, kind=c_int)
    t_initial_c = real(0.0, kind=c_double)
    t_final_c = real(time_step, kind=c_double)

    ! Set up the state pointer
    this%phlex_state_c_ptr = phlex_state_c_ptr

#ifdef PMC_USE_SUNDIALS
    ! Run the solver
    solver_status = int(integrate_solver(neq_c, c_loc(state_array), &
            c_loc(this%abs_tol_c), rel_tol_c, max_steps_c, t_initial_c, &
            t_final_c, c_loc(this_ptr)), kind=i_kind)
#else
    call die_msg(764363956, "No available solver.")
#endif
  end function solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the time derivative f(t,y) when requested by the solver
  subroutine deriv_func(num_spec_c, curr_time_c, &
                  state_array_c_p, deriv_c_p, integration_data_c_p) bind (c)

    use iso_c_binding

    !> Number of species
    integer(kind=c_int), intent(in), value :: num_spec_c
    !> Current integration time (s)
    real(kind=c_double), intent(in), value :: curr_time_c
    !> State array
    type(c_ptr), intent(in), value :: state_array_c_p
    !> Time derivative to be calculated f(t,y)
    type(c_ptr), intent(in), value :: deriv_c_p
    !> Integration data
    type(c_ptr), intent(in), value :: integration_data_c_p

    ! Equivalent fortran variables
    ! Number of species
    integer(kind=i_kind) :: num_spec
    ! Current integration time (s)
    real(kind=dp) :: curr_time
    ! Time derivative to be calculated f(t,y)
    real(kind=dp), pointer :: deriv(:)
    ! Integration data
    type(integration_data_t), pointer :: integration_data

    ! Map c variables to fortran variables
    num_spec = int(num_spec_c, kind=i_kind)
    curr_time = real(curr_time_c, kind=dp)
    call c_f_pointer(deriv_c_p, deriv, (/ num_spec /))
    call c_f_pointer(integration_data_c_p, integration_data)

    ! Reset the time derivative
    deriv(:) = 0.0

    ! Calculate the time derivative
    call integration_data%deriv_func_ptr(curr_time, deriv, &
            integration_data%phlex_core_c_ptr, &
            integration_data%phlex_state_c_ptr)

  end subroutine deriv_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the Jacobian matrix J(t,y) when requested by the solver
  subroutine jac_func(num_spec_c, curr_time_c, &
                  state_array_c_p, jac_c_p, integration_data_c_p) bind (c)

    use iso_c_binding

    !> Number of species
    integer(kind=c_int), intent(in), value :: num_spec_c
    !> Current integration time (s)
    real(kind=c_double), intent(in), value :: curr_time_c
    !> State array
    type(c_ptr), intent(in), value :: state_array_c_p
    !> Jacobian matrix to be calculated J(t,y)
    type(c_ptr), intent(in), value :: jac_c_p
    !> Integration data
    type(c_ptr), intent(in), value :: integration_data_c_p

    ! Equivalent fortran variables
    ! Number of species
    integer(kind=i_kind) :: num_spec
    ! Current integration time (s)
    real(kind=dp) :: curr_time
    ! Jacobian matrix to be calculated J(t,y)
    real(kind=dp), pointer :: jac(:,:)
    ! Integration data
    type(integration_data_t), pointer :: integration_data

    ! Map c variables to fortran variables
    num_spec = int(num_spec_c, kind=i_kind)
    curr_time = real(curr_time_c, kind=dp)
    call c_f_pointer(jac_c_p, jac, (/ num_spec, num_spec /))
    call c_f_pointer(integration_data_c_p, integration_data)

    ! Reset the Jacobian matrix
    jac(:,:) = 0.0

    ! Calculate the Jacobian matrix
    call integration_data%jac_func_ptr(curr_time, jac, &
            integration_data%phlex_core_c_ptr, &
            integration_data%phlex_state_c_ptr)

  end subroutine jac_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the return code from the solve() function.
  subroutine check_status(this, value)

    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> Return code to check.
    integer(kind=i_kind), intent(in) :: value

    if (value == PMC_INTEGRATION_SUCCESS) then
       return
    elseif (value == PMC_INTEGRATION_NO_AVAIL_SOLVER) then
       call die_msg(592214019, "integrate_data: " &
            // "no available solver")

#ifdef PMC_USE_SUNDIALS
    elseif (value == PMC_INTEGRATION_INIT_Y) then
       call die_msg(592432510, "integrate_data: " &
            // "failed to allocate y vector")
    elseif (value == PMC_INTEGRATION_INIT_ABSTOL) then
       call die_msg(592617387, "integrate_data: " &
            // "failed to allocate abstol vector")
    elseif (value == PMC_INTEGRATION_INIT_CVODE_MEM) then
       call die_msg(592785457, "integrate_data: " &
            // "failed to create the solver")
    elseif (value == PMC_INTEGRATION_INIT_CVODE) then
       call die_msg(592970334, "integrate_data: " &
            // "failure to initialize the solver")
    elseif (value == PMC_INTEGRATION_SVTOL) then
       call die_msg(593172018, "integrate_data: " &
            // "failed to set tolerances")
    elseif (value == PMC_INTEGRATION_SET_MAX_STEPS) then
       call die_msg(593810684, "integrate_data: " &
            // "failed to set maximum steps")
    elseif (value == PMC_INTEGRATION_SET_MAX_CONV_FAILS) then
       call die_msg(763861174, "integrate_data: " &
            // "failed to set maximum convergence failures")
    elseif (value == PMC_INTEGRATION_FAIL) then
       call die_msg(594113210, "integrate_data: solver failed")
    elseif (value==PMC_INTEGRATION_DENSE_JAC) then
       call die_msg(594281280, "integrate_data: " &
            // "failed to set dense jacobian solver")
    elseif (value== PMC_INTEGRATION_JAC_FUNC) then
       call die_msg(594415736, "integrate_data: " &
            // "failed to set jacobian function")
    elseif (value==PMC_INTEGRATION_SET_USER_DATA) then
       call die_msg(594533385, "integrate_data: " &
            // "failed to set user data")
    elseif (value==PMC_INTEGRATION_SPGMR_LINEAR_SOLVER ) then
       call die_msg(786820319, "integrate_data: " &
            // "failed to set SPGMR linear solver")
    elseif (value==PMC_INTEGRATION_SET_PRECONDITIONER ) then
       call die_msg(281613914, "integrate_data: " &
            // "failed to set preconditioner functions")
    elseif (value==PMC_INTEGRATION_WRONG_PRECISION) then
       call die_msg(594651034, "integrate_data: " &
            // "SUNDIALS was not compiled for double precision variables")
#endif

    else
       call die_msg(594768683, "integrate_data: unknown return code: " &
            // trim(to_string(value)))
    end if

    return

  end subroutine check_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check whether a solver is available for the integration 
  logical function is_solver_available(this)

    !> Integration data
    class(integration_data_t), intent(in) :: this

#ifdef PMC_USE_SUNDIALS
    is_solver_available = .true.
#else
    is_solver_available = .false.
#endif

  end function is_solver_available

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets the absolute integration tolerance for species by index
  real(kind=dp) function get_abs_tol(this, i_spec)

    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> Species index
    integer(kind=i_kind), intent(in) :: i_spec

    call assert_msg(279436430, i_spec.gt.0 .and. i_spec.le.size(this%abs_tol_c), &
            "Request for absolute integration tolerance out-of-bounds: "// &
            to_string(i_spec))

    get_abs_tol = this%abs_tol_c(i_spec)

  end function get_abs_tol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the integration data
  subroutine do_print(this, file_unit)

    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), optional :: file_unit

    integer(kind=i_kind) :: i_spec
    integer(kind=i_kind) :: f_unit = 6

    if (present(file_unit)) f_unit = file_unit

    write(f_unit,*) "*** Integration data ***"
    write(f_unit,*) "  Number of integration species: ", size(this%abs_tol_c)
    write(f_unit,*) "  Absolute tolerances:"
    do i_spec = 1, size(this%abs_tol_c)
      write(f_unit,*) "    index: ", i_spec, " tolerance: ", this%abs_tol_c(i_spec)
    end do
    write(f_unit,*) "  Relative tolerance: ", this%rel_tol
    write(f_unit,*) "  Maximum number of integration steps: ", this%max_steps

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_integration_data
