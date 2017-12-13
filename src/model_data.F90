! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_model_data module.

!> The model_data_t structure and associated subroutines.
module pmc_model_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_model_state
  use pmc_chem_spec_state
  use pmc_chem_spec_data
  use pmc_mechanism_data
  use pmc_integration_data

  implicit none
  private

  public :: model_data_t

  !> Part-MC model data
  !!
  !! Contains all time-invariant data for a Part-MC model run.
  type :: model_data_t
    !> Chemical mechanisms
    type(mechanism_data_t), pointer :: mechanism(:)
    !> Chemical species data
    type(chem_spec_data_t), pointer :: chem_spec_data
    !> Integration data
    type(integration_data_t), pointer, private :: integration_data => null()
  contains
    !> Load model data
    procedure :: load => pmc_model_data_load
    !> Initialize the model
    procedure :: initialize => pmc_model_data_initialize
    !> Find a mechanism by name
    procedure :: find_mechanism => pmc_model_data_find_mechanism
    !> Add a mechanism to the model
    procedure :: add_mechanism => pmc_model_data_add_mechanism
    !> Get a new model state variable
    procedure :: new_state => pmc_model_data_new_state
    !> Run the chemical mechanisms
    procedure :: do_chemistry => pmc_model_data_do_chemistry
    !> Determine the number of bytes required to pack the variable
    procedure :: pack_size => pmc_model_data_pack_size
    !> Pack the given variable into a buffer, advancing position
    procedure :: bin_pack => pmc_model_data_bin_pack
    !> Unpack the given variable from a buffer, advancing position
    procedure :: bin_unpack => pmc_model_data_bin_unpack
  end type model_data_t

  !> Constructor for chem_spec_data_t
  interface model_data_t
    procedure :: pmc_model_data_constructor
  end interface model_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for model_data_t
  function pmc_model_data_constructor(input_file_path) result(new_obj)

    !> A new set of model parameters
    type(model_data_t), pointer :: new_obj
    !> Part-MC input file paths
    type(string_t), allocatable, optional :: input_file_path(:)

    allocate(new_obj)
    allocate(new_obj%mechanism(0))
    new_obj%chem_spec_data => chem_spec_data_t()

    if (present(input_file_path)) then
      call new_obj%load(input_file_path)
    end if

  end function pmc_model_data_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load model data from input files
#ifdef PMC_USE_JSON
  !! Reads json files containing model object data. The general format of the
  !! json files should be the following:
  !!
  !! {"pmc-data" : [
  !!   {
  !!     "type" : "OBJECT_TYPE",
  !!     ...
  !!   },
  !!   {
  !!     "type" : "OBJECT_TYPE",
  !!     ...
  !!   },
  !!   ...
  !! ]}
  !!
  !! Each json input file should contain exactly one json object with a 
  !! single key-value pair "pmc-data" whose value is an array of json objects.
  !! Additional top-level key-value pairs will be ignored.
  !! Each of the json objects in the pmc-data array must contain a key-value
  !! pair "type" whose value is a string referenceing a valid Part-MC object.
  !! The valid values for type are:
  !!     MECHANISM
  !!     GAS_SPEC
  !!     AERO_SPEC
  !!     AERO_REP
  !! Refer to specific Part-MC data type documentation for the required format
  !! for input objects.
#endif
  subroutine pmc_model_data_load(this, input_file_path)

    !> Model data
    class(model_data_t), intent(inout) :: this
    !> Part-MC input file paths
    type(string_t), allocatable :: input_file_path(:)

    integer(kind=i_kind) :: i_file, i_mech
#ifdef PMC_USE_JSON
    type(json_core), target :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next

    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val
    character(len=:), allocatable :: str_val

    logical :: found

    j_obj => null()
    j_next => null()
    do i_file = 1, size(input_file_path)
      call j_file%initialize()
      call j_file%get_core(json)
      call j_file%load_file(filename = input_file_path(i_file)%string)
      call j_file%get('pmc-data(1)', j_obj)
      do while (associated(j_obj))
        call json%get(j_obj, 'type', unicode_str_val, found)
        call assert_msg(689470331, found, "Missing type in json input file "//&
                input_file_path(i_file)%string)
        str_val = unicode_str_val
        if (str_val.eq.'MECHANISM') then
          call json%get(j_obj, 'name', unicode_str_val, found)
          call assert_msg(822680732, found, "Missing mechanism name in file "//&
                  input_file_path(i_file)%string)
          str_val = unicode_str_val
          if (.not.this%find_mechanism(str_val, i_mech)) then
            call this%add_mechanism(str_val)
            i_mech = size(this%mechanism)
          end if
          call this%mechanism(i_mech)%load(json, j_obj)
        else if (str_val.eq.'GAS_SPEC') then
          call this%chem_spec_data%load(json, j_obj)
        else if (str_val.eq.'AERO_SPEC') then
          call this%chem_spec_data%load(json, j_obj)
        else
          call die_msg(448039776, "Received invalid json input object type: "//&
                  str_val)
        end if
        j_next => j_obj
        call json%get_next(j_next, j_obj)
      end do
      call j_file%destroy()
    end do
#else
    call warn_msg(350136328, "No support for input files.");
#endif

  end subroutine pmc_model_data_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the model data
  subroutine pmc_model_data_initialize(this)

    use iso_c_binding
          
    !> Model data
    class(model_data_t), target, intent(inout) :: this

    integer(kind=i_kind) :: i_mech
    procedure(integration_data_deriv_func), pointer :: deriv_func
    procedure(integration_data_jac_func), pointer :: jac_func
    real(kind=dp), pointer :: abs_tol(:)
    type(model_data_t), pointer :: this_ptr

    this_ptr => this
    deriv_func => pmc_model_data_calc_derivative
    jac_func => pmc_model_data_calc_jacobian

    do i_mech = 1, size(this%mechanism)
      call this%mechanism(i_mech)%initialize(this%chem_spec_data)
    end do

    ! Set up the integrator
    abs_tol => this%chem_spec_data%get_abs_tolerances()
    this%integration_data => integration_data_t(c_loc(this_ptr), deriv_func, &
            jac_func, abs_tol)

    deallocate(abs_tol)

  end subroutine pmc_model_data_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find a mechanism by name in the model data
  logical function pmc_model_data_find_mechanism(this, mech_name, mech_id) &
                  result(found)

    !> Model data
    class(model_data_t), intent(in) :: this
    !> Mechanism name to search for
    character(len=:), allocatable :: mech_name
    !> Index of mechanism in the array
    integer(kind=i_kind) :: mech_id

    found = .false.
    
    do mech_id = 1, size(this%mechanism)
      if (this%mechanism(mech_id)%name().eq.mech_name) then
        found = .true.
        return
      end if
    end do
    mech_id = 0

  end function pmc_model_data_find_mechanism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a chemical mechanism to the model data
  subroutine pmc_model_data_add_mechanism(this, mech_name)

    !> Model data
    class(model_data_t), intent(inout) :: this
    !> Mechanism name
    character(len=:), allocatable :: mech_name

    type(mechanism_data_t), pointer :: new_mechanism(:)

    allocate(new_mechanism(size(this%mechanism)+1))

    new_mechanism(1:size(this%mechanism)) = &
            this%mechanism(1:size(this%mechanism))

    new_mechanism(size(new_mechanism)) = mechanism_data_t(mech_name)

    deallocate(this%mechanism)
    this%mechanism => new_mechanism

  end subroutine pmc_model_data_add_mechanism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a model state variable based on the this set of model data
  function pmc_model_data_new_state(this) result(new_state)

    !> New model state
    type(model_state_t) :: new_state
    !> Chemical model
    class(model_data_t), intent(in) :: this

    new_state%chem_spec_state = this%chem_spec_data%new_state()

  end function pmc_model_data_new_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the chemical mechanism
  subroutine pmc_model_data_do_chemistry(this, model_state, time_step)
 
    use iso_c_binding

    !> Chemical model
    class(model_data_t), intent(in) :: this
    !> Current model state
    type(model_state_t), intent(inout), target :: model_state
    !> Time step over which to integrate (s)
    real(kind=dp), intent(in) :: time_step

    integer(kind=i_kind) :: solver_status
    real(kind=dp), pointer :: state_array(:)

    state_array => model_state%chem_spec_state%conc

    ! Run integration
    solver_status = this%integration_data%solve(state_array, &
            c_loc(model_state), time_step)

    ! Evaluate the solver status
    call this%integration_data%check_status(solver_status)

  end subroutine pmc_model_data_do_chemistry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the mechanism
  integer(kind=i_kind) function pmc_model_data_pack_size(this) &
                  result (pack_size)

    !> Chemical model
    class(model_data_t), intent(in) :: this
    
    integer(kind=i_kind) :: i_mech

    pack_size =  pmc_mpi_pack_size_integer(size(this%mechanism))
    do i_mech = 1, size(this%mechanism)
      pack_size = pack_size + this%mechanism(i_mech)%pack_size()
    end do

  end function pmc_model_data_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_model_data_bin_pack(this, buffer, pos)

    !> Chemical model
    class(model_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: i_mech, prev_position

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, size(this%mechanism)
    do i_mech = 1, size(this%mechanism)
      call this%mechanism(i_mech)%bin_pack(buffer, pos)
    end do
    call assert(184050835, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_model_data_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine pmc_model_data_bin_unpack(this, buffer, pos)

    !> Chemical model
    class(model_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: i_mech, prev_position, num_mech

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, num_mech)
    allocate(this%mechanism(num_mech))
    do i_mech = 1, size(this%mechanism)
      call this%mechanism(i_mech)%bin_unpack(buffer, pos)
    end do
    call assert(291557168, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_model_data_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the time derivative f(t,y)
  subroutine pmc_model_data_calc_derivative(curr_time, deriv, &
                  model_data_c_ptr, model_state_c_ptr)

    use iso_c_binding

    !> Current solver time (s)
    real(kind=dp), intent(in) :: curr_time
    !> Time derivative to calculate
    real(kind=dp), intent(inout), pointer :: deriv(:)
    !> Pointer to model data
    type(c_ptr), intent(in) :: model_data_c_ptr
    !> Pointer to model state
    type(c_ptr), intent(in) :: model_state_c_ptr

    type(model_data_t), pointer :: model_data
    type(model_state_t), pointer :: model_state
    integer(kind=i_kind) :: i_mech

    call c_f_pointer(model_data_c_ptr, model_data)
    call c_f_pointer(model_state_c_ptr, model_state)

    ! Calculate f(t,y)
    do i_mech=1, size(model_data%mechanism)
      call model_data%mechanism(i_mech)%get_func_contrib(model_state, deriv)
    end do

  end subroutine pmc_model_data_calc_derivative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the Jacobian matrix J(t,y)
  subroutine pmc_model_data_calc_jacobian(curr_time, jac, &
                  model_data_c_ptr, model_state_c_ptr)

    use iso_c_binding

    !> Current solver time (s)
    real(kind=dp), intent(in) :: curr_time
    !> Jacobian matrix to calculate
    real(kind=dp), intent(inout), pointer :: jac(:,:)
    !> Pointer to model data
    type(c_ptr), intent(in) :: model_data_c_ptr
    !> Pointer to model state
    type(c_ptr), intent(in) :: model_state_c_ptr

    type(model_data_t), pointer :: model_data
    type(model_state_t), pointer :: model_state
    integer(kind=i_kind) :: i_mech

    call c_f_pointer(model_data_c_ptr, model_data)
    call c_f_pointer(model_state_c_ptr, model_state)

    ! Calculate J(t,y)
    do i_mech=1, size(model_data%mechanism)
      call model_data%mechanism(i_mech)%get_jac_contrib(model_state, jac)
    end do

  end subroutine pmc_model_data_calc_jacobian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_model_data
