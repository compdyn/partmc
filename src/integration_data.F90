! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_integration_data module.

!> The integration_data_t structure and associated subroutines.
module pmc_integration_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_util,                       only : assert_msg, to_string, die_msg
  use pmc_phlex_state

  use iso_c_binding

  implicit none
  private

  !> Default relative tolerance for integration
  real(kind=dp), parameter :: PMC_INTEGRATION_DEFAULT_REL_TOL = 1.0D-8
  !> Default max number of integration steps
  integer(kind=i_kind), parameter :: PMC_INTEGRATION_DEFAULT_MAX_STEPS = 10000000

  !> Result code indicating successful completion.
  integer, parameter :: PMC_INTEGRATION_SUCCESS              = 0
  !> Result code indicating no available integration routine
  integer, parameter :: PMC_INTEGRATION_NO_AVAIL_SOLVER      = 12
#ifdef PMC_USE_SUNDIALS
  !> Result code indicating failure to allocate \c y vector.
  integer, parameter :: PMC_INTEGRATION_INIT_Y               = 1
  !> Result code indicating failure to allocate \c abstol vector.
  integer, parameter :: PMC_INTEGRATION_INIT_ABSTOL          = 2
  !> Result code indicating failure to create the solver.
  integer, parameter :: PMC_INTEGRATION_INIT_CVODE_MEM       = 3
  !> Result code indicating failure to initialize the solver.
  integer, parameter :: PMC_INTEGRATION_INIT_CVODE           = 4
  !> Result code indicating failure to set tolerances.
  integer, parameter :: PMC_INTEGRATION_SVTOL                = 5
  !> Result code indicating failure to set maximum steps.
  integer, parameter :: PMC_INTEGRATION_SET_MAX_STEPS        = 6
  !> Result code indicating failure of the solver.
  integer, parameter :: PMC_INTEGRATION_FAIL                 = 7
  !> Result code indicating failure to set sparse Jacobian solver
  integer, parameter :: PMC_INTEGRATION_SPARSE_JAC           = 8
  !> Result code indicating failure to set Jacobian function
  integer, parameter :: PMC_INTEGRATION_JAC_FUNC             = 9
  !> Result code indicating failure to set user data
  integer, parameter :: PMC_INTEGRATION_SET_USER_DATA        = 10
  !> Result code indicating SUNDIALS realtype is not set to double precision
  integer, parameter :: PMC_INTEGRATION_WRONG_PRECISION      = 11
  !> Result code indicating failure to get the sparse linear solver
  integer, parameter :: PMC_INTEGRATION_SPARSE_LINEAR_SOLVER = 13
  !> Result code indicating failure to set the linear solver
  integer, parameter :: PMC_INTEGRATION_SET_LINEAR_SOLVER    = 14
  !> Result code indicating failure to set the maximum number of convergence failures
  integer, parameter :: PMC_INTEGRATION_SET_MAX_CONV_FAILS   = 15
  !> Result code indicating failure to set the SPGMR solver
  integer, parameter :: PMC_INTEGRATION_SPGMR_LINEAR_SOLVER  = 16
  !> Result code indicating failure to set the preconditioner functions
  integer, parameter :: PMC_INTEGRATION_SET_PRECONDITIONER   = 17

  !> Interface to c ODE solver function
  interface
    integer(kind=c_int) function integrate_solver(neq, x, abs_tol, rel_tol, &
                                 max_steps, t_initial, t_final, sysdata, &
                                 n_jac_elem, jac_col_ptrs, jac_row_ids) bind(c)
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
      !> Number of non-zero Jacobian elements
      integer(kind=c_int), value :: n_jac_elem
      !> Array containing the index of the first element in each column of
      !! the Jacobian matrix with an additional element at the end of the 
      !! array that equals n_jac_elem. The size of the array should be 
      !! neq + 1.
      type(c_ptr), value :: jac_col_ptrs
      !> Array containing the row id of each element in the Jacobian matrix
      !! data arranged as a flat array. The array should contain n_jac_elem
      !! elements.
      type(c_ptr), value :: jac_row_ids
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
    type(c_ptr), public :: phlex_core_c_ptr
    !> PartMC model state
    type(phlex_state_t), pointer, public :: phlex_state
    !> Current time
    real(kind=dp) :: curr_time
    ! Time derivative being calculated f(t,y)
    real(kind=dp), pointer :: deriv(:)
    ! Sparse Jacobian matrix being calculated J(t,y)
    real(kind=dp), pointer :: jac(:)
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
    !> Number of time derivative elements
    integer(kind=i_kind) :: n_deriv_elem
    !> Number of non-zero Jacobian elements  
    integer(kind=i_kind) :: n_jac_elem
    !> Array containing the index of the first element in each column of
    !! the Jacobian matrix with an additional element at the end of the 
    !! array that equals n_jac_elem. The size of the array should be 
    !! neq + 1.
    integer(kind=c_int), pointer :: jac_col_ptrs_c(:)
    !> Array containing the row id of each element in the Jacobian matrix
    !! data arranged as a flat array. The array should contain n_jac_elem
    !! elements.
    integer(kind=c_int), pointer :: jac_row_ids_c(:)
  contains
    !> Integrate over a given time step
    procedure :: solve
    !> Checks the result status code from an integration
    procedure :: check_status
    !> Checks whether a solver is available
    procedure :: is_solver_available
    !> Gets the absolute integration tolerance for a species by index
    procedure :: get_abs_tol
    !> Get the current time during integration
    procedure :: get_curr_time
    !> Add a contribution to the time derivative
    procedure :: add_deriv_contrib
    !> Add a contribution to the Jacobian matrix
    procedure :: add_jac_contrib
    !> Get the id of a Jacobian element in the sparse matrix
    procedure :: get_jac_elem_id
    !> Get a time derivative element value
    procedure :: get_deriv_elem
    !> Get a Jacobian element value
    procedure :: get_jac_elem
    !> Get a pointer to the time derivative function
    procedure :: get_deriv_func
    !> Get a pointer to the Jacobian function
    procedure :: get_jac_func
    !> Reset the time derivative vector
    procedure :: reset_deriv
    !> Reset the Jacobian matrix
    procedure :: reset_jac
    !> Allocate space for the solver arrays. For use in tests
    procedure :: init_solver_arrays
    !> Deallocate space for the solver arrays. For use in tests
    procedure :: kill_solver_arrays
    !> Print the integration data
    procedure :: print => do_print
  end type integration_data_t

  ! Constructor for integration_data_t
  interface integration_data_t
    procedure :: constructor, constructor_empty
  end interface integration_data_t

  abstract interface 
    !> Interface for the time derivative function
    subroutine integration_data_deriv_func(integration_data)
      import integration_data_t
      
      !> Pointer to integration data
      type(integration_data_t), pointer, intent(inout) :: integration_data  
    end subroutine integration_data_deriv_func

    !> Interface for the Jacobian function
    subroutine integration_data_jac_func(integration_data)
      import integration_data_t

      !> Pointer to integration data
      type(integration_data_t), pointer, intent(inout) :: integration_data  
    end subroutine integration_data_jac_func
  end interface 

  public :: integration_data_t, integration_data_deriv_func, &
          integration_data_jac_func, deriv_func, jac_func

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for integration_data_t
  function constructor(phlex_core_c_ptr, deriv_func_ptr, &
                  jac_func_ptr, abs_tol, n_jac_elem, n_jac_col_elem, &
                  jac_row_ids) result(new_obj)
    
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
    !> Number of non-zero Jacobian elements  
    integer(kind=i_kind), intent(in) :: n_jac_elem
    !> Array containing the number of non-zero elements in each column
    !! of the Jacobian matrix. Their sum should equal n_jac_elem
    integer(kind=i_kind), pointer, intent(in) :: n_jac_col_elem(:)
    !> Array containing the row id of each element in the Jacobian matrix
    !! data arranged as a flat array. The array should contain n_jac_elem
    !! elements.
    integer(kind=i_kind), pointer, intent(in)  :: jac_row_ids(:)

    integer(kind=i_kind) :: i_col, col_id

    allocate(new_obj)
    new_obj%phlex_core_c_ptr = phlex_core_c_ptr
    new_obj%deriv_func_ptr => deriv_func_ptr
    new_obj%jac_func_ptr => jac_func_ptr
    allocate(new_obj%abs_tol_c(size(abs_tol)))
    new_obj%abs_tol_c(:) = real(abs_tol(:), kind=c_double)
    new_obj%n_deriv_elem = size(n_jac_col_elem)
    new_obj%n_jac_elem = n_jac_elem
    allocate(new_obj%jac_col_ptrs_c(new_obj%n_deriv_elem+1))
    ! Column ids in c SUNMatrix start at 0
    col_id = 0
    do i_col = 1, new_obj%n_deriv_elem
      new_obj%jac_col_ptrs_c(i_col) = int(col_id, kind=c_int)
      col_id = col_id + n_jac_col_elem(i_col)
    end do
    new_obj%jac_col_ptrs_c(new_obj%n_deriv_elem+1) = &
            new_obj%jac_col_ptrs_c(new_obj%n_deriv_elem)
    allocate(new_obj%jac_row_ids_c(size(jac_row_ids)))
    ! Row ids in c SUNMatrix start at 0
    new_obj%jac_row_ids_c(:) = int(jac_row_ids(:)-1, kind=c_int)

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
  function solve(this, phlex_state, time_step) result(solver_status)

    !> Solver status
    integer(kind=i_kind) :: solver_status
    !> Integration data
    class(integration_data_t), target, intent(inout) :: this
    !> Model state pointer. This is a pointer for use by the derivative
    !! functions to allow access to state variables including those that are 
    !! not solved for during the integration (e.g., temperature, pressure)
    type(phlex_state_t), pointer, intent(in) :: phlex_state
    !> Time step (s)
    real(kind=dp), intent(in) :: time_step

    ! Variables required by the c ODE solver
    !
    ! State variable array (Species to solve for)
    real(kind=dp), pointer :: state_array(:)
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
    ! Number of non-zero Jacobian elements
    integer(kind=c_int) :: n_jac_elem_c

    this_ptr => this

    ! Set up solver variables
    allocate(state_array(size(phlex_state%state_var)))
    state_array(:) = phlex_state%state_var(:)
    neq_c = int(size(state_array), kind=c_int)
    rel_tol_c = real(this%rel_tol, kind=c_double)
    max_steps_c = int(this%max_steps, kind=c_int)
    t_initial_c = real(0.0, kind=c_double)
    t_final_c = real(time_step, kind=c_double)
    n_jac_elem_c = int(this%n_jac_elem, kind=c_int)

    ! Set up the state pointer
    this%phlex_state => phlex_state

#ifdef PMC_USE_SUNDIALS
    ! Run the solver
    solver_status = int(integrate_solver(neq_c, c_loc(state_array), &
            c_loc(this%abs_tol_c), rel_tol_c, max_steps_c, t_initial_c, &
            t_final_c, c_loc(this_ptr), n_jac_elem_c, &
            c_loc(this%jac_col_ptrs_c), c_loc(this%jac_row_ids_c)), kind=i_kind)
#else
    call die_msg(764363956, "No available solver.")
#endif

    ! Copy solved concentrations
    phlex_state%state_var(:) = state_array(:)

    ! Deallocate state array
    ! TODO move allocation / deallocation to initialization routine
    deallocate(state_array)

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
    ! Integration data
    type(integration_data_t), pointer :: integration_data
    ! State array
    real(kind=c_double), pointer :: state_array(:)

    ! Attach to the integration data object
    call c_f_pointer(integration_data_c_p, integration_data)

    ! Map c variables to fortran variables
    num_spec = int(num_spec_c, kind=i_kind)
    integration_data%curr_time = real(curr_time_c, kind=dp)
    call c_f_pointer(deriv_c_p, integration_data%deriv, (/ num_spec /))
    call c_f_pointer(state_array_c_p, state_array, (/ num_spec /))

    ! Update the species concentrations
    integration_data%phlex_state%state_var(:) = real(state_array(:), kind=dp)

    ! Reset the time derivative
    integration_data%deriv(:) = 0.0

    ! Calculate the time derivative
    call integration_data%deriv_func_ptr(integration_data)

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
    ! State array
    real(kind=c_double), pointer :: state_array(:)

    ! Equivalent fortran variables
    ! Number of species
    integer(kind=i_kind) :: num_spec
    ! Integration data
    type(integration_data_t), pointer :: integration_data

    integer :: i_col, i_row, i_elem

    ! Attach to the integration data object
    call c_f_pointer(integration_data_c_p, integration_data)

    ! Map c variables to fortran variables
    num_spec = int(num_spec_c, kind=i_kind)
    integration_data%curr_time = real(curr_time_c, kind=dp)
    call c_f_pointer(jac_c_p, integration_data%jac, &
            (/ integration_data%n_jac_elem /))
    call c_f_pointer(state_array_c_p, state_array, (/ num_spec /))

    ! Update the species concentrations
    integration_data%phlex_state%state_var(:) = real(state_array(:), kind=dp)

    ! Reset the Jacobian matrix
    integration_data%jac(:) = 0.0

    ! Calculate the Jacobian matrix
    call integration_data%jac_func_ptr(integration_data)

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
    elseif (value==PMC_INTEGRATION_SPARSE_JAC) then
       call die_msg(594281280, "integrate_data: " &
            // "failed to initialize sparse jacobian")
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

  !> Get the curret time during integration
  function get_curr_time(this) result (curr_time)

    !> Current model time (s)
    real(kind=dp) :: curr_time
    !> Integration data
    class(integration_data_t), intent(in) :: this

    curr_time = this%curr_time

  end function get_curr_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a contribution to a time derivative vector
  subroutine add_deriv_contrib(this, spec_id, contrib_value)

    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> Species ID to add contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    !> Contribution value to add
    real(kind=dp), intent(in) :: contrib_value

    this%deriv(spec_id) = this%deriv(spec_id) + contrib_value

  end subroutine add_deriv_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a contribution to a Jacobian matrix
  subroutine add_jac_contrib(this, jac_id, contrib_value)

    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> Index in the sparse Jacobian data array to add contribution to
    integer(kind=i_kind), intent(in) :: jac_id
    !> Contribution value to add
    real(kind=dp), intent(in) :: contrib_value

    this%jac(jac_id) = this%jac(jac_id) + contrib_value

  end subroutine add_jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index of an element in the sparse Jacobian matrix
  function get_jac_elem_id(this, dep_spec_id, indep_spec_id) result (jac_id)

    !> Sparse Jacobian id of requested element
    integer(kind=i_kind) :: jac_id
    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> Index of dependent species
    integer(kind=i_kind), intent(in) :: dep_spec_id
    !> Index of independent species
    integer(kind=i_kind), intent(in) :: indep_spec_id

    integer(kind=i_kind) :: i_elem

    jac_id = 0
    ! indexes for c SUNMatrix start at 0
    do i_elem = this%jac_col_ptrs_c(indep_spec_id) + 1, &
            this%jac_col_ptrs_c(indep_spec_id+1)
      if (this%jac_row_ids_c(i_elem).eq.dep_spec_id-1) jac_id = i_elem
    end do

    if (jac_id.eq.0) call this%print()

    call assert_msg(247151806, jac_id.gt.0, "Missing Jacobian element for "// &
            "i="//trim(to_string(dep_spec_id))//"; j="// &
            trim(to_string(indep_spec_id)))

  end function get_jac_elem_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a time derivative element value
  function get_deriv_elem(this, spec_id) result (val)

    !> Value of element in derivative vector
    real(kind=dp) :: val
    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> Species id
    integer(kind=i_kind), intent(in) :: spec_id

    val = this%deriv(spec_id)

  end function get_deriv_elem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a Jacobian element value
  function get_jac_elem(this, dep_spec_id, indep_spec_id) result (val)

    !> Value of element in Jacobian matrix
    real(kind=dp) :: val
    !> Integration data
    class(integration_data_t), intent(in) :: this
    !> Dependent species id
    integer(kind=i_kind), intent(in) :: dep_spec_id
    !> Independent species id
    integer(kind=i_kind), intent(in) :: indep_spec_id

    integer(kind=i_kind) :: i_elem

    val = real(0.0, kind=dp)

    ! Indexes for c SUNMatrix start at 0
    do i_elem = this%jac_col_ptrs_c(indep_spec_id) + 1, &
            this%jac_col_ptrs_c(indep_spec_id+1)
      if (this%jac_row_ids_c(i_elem).eq.dep_spec_id-1) then
        val = this%jac(i_elem)
      end if
    end do

  end function get_jac_elem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a pointer to the time derivative function
  function get_deriv_func(this) result (func_ptr)

    !> Pointer to the time derivative function
    procedure(integration_data_deriv_func), pointer :: func_ptr
    !> Integration data
    class(integration_data_t), intent(in) :: this

    func_ptr => this%deriv_func_ptr

  end function get_deriv_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a pointer to the Jacobian function
  function get_jac_func(this) result (func_ptr)

    !> Pointer to the Jacobian function
    procedure(integration_data_jac_func), pointer :: func_ptr
    !> Integration data
    class(integration_data_t), intent(in) :: this

    func_ptr => this%jac_func_ptr

  end function get_jac_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reset the time derivative vector
  subroutine reset_deriv(this)

    !> Integration data
    class(integration_data_t), intent(inout) :: this

    this%deriv(:) = 0.0

  end subroutine reset_deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reset the Jacobian matrix 
  subroutine reset_jac(this)

    !> Integration data
    class(integration_data_t), intent(inout) :: this

    this%jac(:) = 0.0

  end subroutine reset_jac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate space for time derivative and Jacobian arrays. For use when
  !! calculating time derivative and Jacobian for testing
  subroutine init_solver_arrays(this, phlex_state)

    !> Integration data
    class(integration_data_t), intent(inout) :: this
    !> Phlex state
    type(phlex_state_t), pointer, intent(in) :: phlex_state

    this%phlex_state => phlex_state
    allocate(this%deriv(this%n_deriv_elem))
    allocate(this%jac(this%n_jac_elem))

  end subroutine init_solver_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deallocate space for time derivative and Jacobian arrays. For use when
  !! calculating time derivative and Jacobian for testing
  subroutine kill_solver_arrays(this)

    !> Integration data
    class(integration_data_t), intent(inout) :: this

    this%phlex_state => NULL()
    deallocate(this%deriv)
    deallocate(this%jac)

  end subroutine kill_solver_arrays

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
    write(f_unit,*) "  Number of non-zero Jacobian elements: ", this%n_jac_elem
    write(f_unit,*) "  Index of each Jacobian column: ", this%jac_col_ptrs_c
    write(f_unit,*) "  Row ids for each Jacobian element: ", this%jac_row_ids_c

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_integration_data
