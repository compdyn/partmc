! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_phlex_core module.

!> The phlex_core_t structure and associated subroutines.
module pmc_phlex_core

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_phlex_state
  use pmc_chem_spec_data
  use pmc_mechanism_data
  use pmc_integration_data
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_phase_data

  implicit none
  private

  public :: phlex_core_t

  !> Part-MC model data
  !!
  !! Contains all time-invariant data for a Part-MC model run.
  type :: phlex_core_t
    !> Chemical mechanisms
    type(mechanism_data_t), pointer :: mechanism(:)
    !> Chemical species data
    type(chem_spec_data_t), pointer :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer :: aero_rep(:)
    !> Aerosol phases
    type(aero_phase_data_t), pointer :: aero_phase(:)
    !> Integration data
    type(integration_data_t), pointer, private :: integration_data => null()
  contains
    !> Load a set of configuration files
    procedure :: load_files => pmc_phlex_core_load_files
    !> Load model data from a configuration file
    procedure :: load => pmc_phlex_core_load
    !> Initialize the model
    procedure :: initialize => pmc_phlex_core_initialize
    !> Find a mechanism by name
    procedure :: find_mechanism => pmc_phlex_core_find_mechanism
    !> Add a mechanism to the model
    procedure :: add_mechanism => pmc_phlex_core_add_mechanism
    !> Find an aerosol representation index by name
    procedure :: find_aero_rep_id_by_name => pmc_phlex_core_find_aero_rep_id_by_name
    !> Find an aerosol representation by name
    procedure :: find_aero_rep_by_name => pmc_phlex_core_find_aero_rep_by_name
    !> Find an aerosol representation or it's id
    generic :: find_aero_rep => find_aero_rep_id_by_name, find_aero_rep_by_name
    !> Add an aerosol representation to the model
    procedure :: add_aero_rep => pmc_phlex_core_add_aero_rep
    !> Find an aerosol phase by name
    procedure :: find_aero_phase => pmc_phlex_core_find_aero_phase
    !> Add an aerosol phase to the model
    procedure :: add_aero_phase => pmc_phlex_core_add_aero_phase
    !> Get a new model state variable
    procedure :: new_state => pmc_phlex_core_new_state
    !> Run the chemical mechanisms
    procedure :: solve => pmc_phlex_core_solve
    !> Determine the number of bytes required to pack the variable
    procedure :: pack_size => pmc_phlex_core_pack_size
    !> Pack the given variable into a buffer, advancing position
    procedure :: bin_pack => pmc_phlex_core_bin_pack
    !> Unpack the given variable from a buffer, advancing position
    procedure :: bin_unpack => pmc_phlex_core_bin_unpack
  end type phlex_core_t

  !> Constructor for phlex_core_t
  interface phlex_core_t
    procedure :: pmc_phlex_core_constructor
  end interface phlex_core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for phlex_core_t
  function pmc_phlex_core_constructor(input_file_path) result(new_obj)

    !> A new set of model parameters
    type(phlex_core_t), pointer :: new_obj
    !> Part-MC input file paths
    character(len=:), allocatable, optional :: input_file_path

    allocate(new_obj)
    allocate(new_obj%mechanism(0))
    new_obj%chem_spec_data => chem_spec_data_t()

    if (present(input_file_path)) then
      call new_obj%load_files(input_file_path)
    end if

  end function pmc_phlex_core_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a set of model data files
#ifdef PMC_USE_JSON
  !! Reads a list of files containing model data. The format of the json file
  !! should be the following:
  !!
  !! { "pmc-files" : [
  !!   "file_one.json",
  !!   "some_dir/file_two.json",
  !!   ...
  !! ]}
  !!
  !! The input file should be in json format and contain a single key-value
  !! pair named pmc-files whose value is an array of string with paths to the
  !! set of configuration files to load. Input files should be json format.
#endif
  subroutine pmc_phlex_core_load_files(this, input_file_path)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Part-MC input file paths
    character(len=:), allocatable :: input_file_path

#ifdef PMC_USE_JSON
    type(json_core), target :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next

    logical(kind=json_lk) :: found
    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val

    integer(kind=json_ik) :: i_file, num_files
    type(string_t), allocatable :: file_list(:)

    call j_file%initialize()
    call j_file%get_core(json)
    call j_file%load_file(filename = input_file_path)
    call j_file%get('pmc-files(1)', j_obj, found)
    call assert_msg(405149265, found, &
            "Could not find pmc-files object in input file: "//input_file_path)
    call j_file%info('pmc-files', n_children = num_files)
    call assert_msg(411804027, num_files.gt.0, &
            "No file names were found in "//input_file_path)
    allocate(file_list(num_files))
    j_next => null()
    i_file = 1
    do while (associated(j_obj))
      call json%get(j_obj, unicode_str_val)
      file_list(i_file)%string = unicode_str_val
      i_file = i_file + 1
      call json%get_next(j_obj, j_next)
      j_obj => j_next
    end do
    call j_file%destroy()

    call this%load(file_list)
#else
    call warn_msg(171627969, "No support for input files.");
#endif

  end subroutine pmc_phlex_core_load_files
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load model data from input files
#ifdef PMC_USE_JSON
  !! Reads json files containing model object data. The general format of the
  !! json files should be the following:
  !!
  !! { "pmc-data" : [
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
  subroutine pmc_phlex_core_load(this, input_file_path)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Part-MC input file paths
    type(string_t), allocatable :: input_file_path(:)

    integer(kind=i_kind) :: i_file, i_mech
#ifdef PMC_USE_JSON
    type(json_core), target :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next

    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val
    character(len=:), allocatable :: str_val

    type(aero_phase_data_t), pointer :: new_aero_phase(:)
    type(aero_rep_data_ptr), pointer :: new_aero_rep(:)
    type(aero_rep_factory_t), pointer :: aero_rep_factory

    type(aero_phase_data_t), pointer :: aero_phase
    type(aero_rep_data_ptr), pointer :: aero_rep_ptr

    integer(kind=i_kind) :: i_rep, i_phase
    logical :: found

    aero_rep_factory = aero_rep_factory_t()

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
        else if (str_val(1:8).eq.'AERO_REP') then
          aero_rep_ptr%val => aero_rep_factory%load(json, j_obj)
          if (this%find_aero_rep(aero_rep_ptr%val%name(), i_rep)) then
            deallocate(aero_rep_ptr)
            call this%aero_rep(i_rep)%val%load(json, j_obj)
          else
            allocate(new_aero_rep(size(this%aero_rep)+1))
            new_aero_rep(1:size(this%aero_rep)) = this%aero_rep(1:size(this%aero_rep))
            new_aero_rep(size(new_aero_rep))%val => aero_rep_ptr%val
            deallocate(aero_rep_ptr)
            deallocate(this%aero_rep)
            this%aero_rep => new_aero_rep
          end if
        else if (str_val.eq.'AERO_PHASE') then
          aero_phase = aero_phase_data_t()
          call aero_phase%load(json, j_obj)
          if (this%find_aero_phase(aero_phase%name(), i_phase)) then
            deallocate(aero_phase)
            call this%aero_phase(i_phase)%load(json, j_obj)
          else
            allocate(new_aero_phase(size(this%aero_phase)+1))
            new_aero_phase(1:size(this%aero_phase)) = this%aero_phase(1:size(this%aero_phase))
            new_aero_phase(size(new_aero_phase)) = aero_phase
            deallocate(this%aero_phase)
            this%aero_phase = new_aero_phase
          end if
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

  end subroutine pmc_phlex_core_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the model data
  subroutine pmc_phlex_core_initialize(this)

    use iso_c_binding
          
    !> Model data
    class(phlex_core_t), target, intent(inout) :: this

    integer(kind=i_kind) :: i_mech, i_phase, i_aero_rep, i_state_var
    procedure(integration_data_deriv_func), pointer :: deriv_func
    procedure(integration_data_jac_func), pointer :: jac_func
    real(kind=dp), pointer :: abs_tol(:)
    type(phlex_core_t), pointer :: this_ptr

    this_ptr => this
    deriv_func => pmc_phlex_core_calc_derivative
    jac_func => pmc_phlex_core_calc_jacobian

    ! Get the size of the gas-phase species on the state array
    i_state_var = this%chem_spec_data%size() + 1

    ! Initialize the aerosol phases
    do i_phase = 1, size(this%aero_phase)
      call this%aero_phase(i_phase)%initialize(this%chem_spec_data)
    end do

    ! Initialize the aerosol representations
    do i_aero_rep = 1, size(this%aero_rep)
      call this%aero_rep(i_aero_rep)%val%initialize(this%aero_phase, &
              i_state_var, i_aero_rep, this%chem_spec_data)
      i_state_var = i_state_var + this%aero_rep(i_aero_rep)%val%size()
    end do

    ! Initialize the mechanisms
    do i_mech = 1, size(this%mechanism)
      call this%mechanism(i_mech)%initialize(this%chem_spec_data)
    end do

    ! Set up the integrator
    abs_tol => this%chem_spec_data%get_abs_tolerances()
    this%integration_data => integration_data_t(c_loc(this_ptr), deriv_func, &
            jac_func, abs_tol)

    deallocate(abs_tol)

  end subroutine pmc_phlex_core_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find a mechanism by name in the model data
  logical function pmc_phlex_core_find_mechanism(this, mech_name, mech_id) &
                  result(found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
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

  end function pmc_phlex_core_find_mechanism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a chemical mechanism to the model data
  subroutine pmc_phlex_core_add_mechanism(this, mech_name)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Mechanism name
    character(len=:), allocatable :: mech_name

    type(mechanism_data_t), pointer :: new_mechanism(:)

    allocate(new_mechanism(size(this%mechanism)+1))

    new_mechanism(1:size(this%mechanism)) = &
            this%mechanism(1:size(this%mechanism))

    new_mechanism(size(new_mechanism)) = mechanism_data_t(mech_name)

    deallocate(this%mechanism)
    this%mechanism => new_mechanism

  end subroutine pmc_phlex_core_add_mechanism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find an aerosol representation id by name in the model data
  logical function pmc_phlex_core_find_aero_rep_id_by_name(this, rep_name, &
                  rep_id) result(found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Aerosol representation name to search for
    character(len=:), allocatable :: rep_name
    !> Index of the representation in the array
    integer(kind=i_kind) :: rep_id

    found = .false.
    
    do rep_id = 1, size(this%aero_rep)
      if (this%aero_rep(rep_id)%val%name().eq.rep_name) then
        found = .true.
        return
      end if
    end do
    rep_id = 0

  end function pmc_phlex_core_find_aero_rep_id_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find an aerosol representation by name in the model data
  logical function pmc_phlex_core_find_aero_rep_by_name(this, rep_name, rep) &
                  result(found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Aerosol representation name to search for
    character(len=:), allocatable :: rep_name
    !> Aerosol representation
    class(aero_rep_data_t), pointer :: rep

    integer(kind=i_kind) :: rep_id

    found = this%find_aero_rep_id_by_name(rep_name, rep_id)
    if (.not.found) return
    
    rep => this%aero_rep(rep_id)%val

  end function pmc_phlex_core_find_aero_rep_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a aerosol representation to the model data
  subroutine pmc_phlex_core_add_aero_rep(this, rep_name)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Aerosol representation name
    character(len=:), allocatable :: rep_name

    type(aero_rep_data_ptr), pointer :: new_aero_rep(:)
    type(aero_rep_factory_t), pointer :: aero_rep_factory

    aero_rep_factory = aero_rep_factory_t()
    allocate(new_aero_rep(size(this%aero_rep)+1))

    new_aero_rep(1:size(this%aero_rep)) = &
            this%aero_rep(1:size(this%aero_rep))
    new_aero_rep(size(new_aero_rep))%val => aero_rep_factory%create(rep_name)

    deallocate(this%aero_rep)
    this%aero_rep => new_aero_rep

  end subroutine pmc_phlex_core_add_aero_rep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find an aerosol phase by name in the model data
  logical function pmc_phlex_core_find_aero_phase(this, phase_name, rep_id) &
                  result(found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Aerosol phase name to search for
    character(len=:), allocatable :: phase_name
    !> Index of the phase in the array
    integer(kind=i_kind) :: rep_id

    found = .false.
    
    do rep_id = 1, size(this%aero_phase)
      if (this%aero_phase(rep_id)%name().eq.phase_name) then
        found = .true.
        return
      end if
    end do
    rep_id = 0

  end function pmc_phlex_core_find_aero_phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a aerosol phase to the model data
  subroutine pmc_phlex_core_add_aero_phase(this, phase_name)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Aerosol phase name
    character(len=:), allocatable :: phase_name

    type(aero_phase_data_t), pointer :: new_aero_phase(:)

    allocate(new_aero_phase(size(this%aero_phase)+1))

    new_aero_phase(1:size(this%aero_phase)) = &
            this%aero_phase(1:size(this%aero_phase))

    new_aero_phase(size(new_aero_phase)) = aero_phase_data_t(phase_name)

    deallocate(this%aero_phase)
    this%aero_phase => new_aero_phase

  end subroutine pmc_phlex_core_add_aero_phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a model state variable based on the this set of model data
  function pmc_phlex_core_new_state(this) result(new_state)

    !> New model state
    type(phlex_state_t), pointer :: new_state
    !> Chemical model
    class(phlex_core_t), intent(in) :: this

    integer(kind=i_kind) :: state_size, i_rep

    new_state => phlex_state_t()

    ! Set up the state variable array
    state_size = this%chem_spec_data%size()
    do i_rep = 1, size(this%aero_rep)
      state_size = state_size + this%aero_rep(i_rep)%val%size()
    end do
    allocate(new_state%state_var(state_size))

    ! Create the aerosol representation states
    allocate(new_state%aero_rep_state(size(this%aero_rep)))
    do i_rep = 1, size(this%aero_rep)
      new_state%aero_rep_state(i_rep)%val => this%aero_rep(i_rep)%val%new_state()
    end do

  end function pmc_phlex_core_new_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the chemical mechanism
  subroutine pmc_phlex_core_solve(this, phlex_state, time_step, rxn_phase)

    use pmc_rxn_data 
    use iso_c_binding

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(inout), target :: phlex_state
    !> Time step over which to integrate (s)
    real(kind=dp), intent(in) :: time_step
    !> Phase to solve - gas, aerosol, or both (default)
    !! Use parameters in pmc_rxn_data to specify phase:
    !! GAS_RXN, AERO_RXN, GAS_AERO_RXN
    integer(kind=i_kind), intent(in), optional :: rxn_phase

    integer(kind=i_kind) :: solver_status, phase
    real(kind=dp), pointer :: state_array(:)

    if (present(rxn_phase)) then
      phase = rxn_phase
    else
      phase = GAS_AERO_RXN
    end if

    if (phase.ne.GAS_RXN .and. &
        phase.ne.AERO_RXN .and. &
        phase.ne.GAS_AERO_RXN) then
      call die_msg(704896254, "Invalid rxn phase specified for chemistry "// &
              "solver: "//to_string(phase))
    end if

    phlex_state%rxn_phase = phase

    state_array => phlex_state%state_var

    ! Run integration
    solver_status = this%integration_data%solve(state_array, &
            c_loc(phlex_state), time_step)

    ! Evaluate the solver status
    call this%integration_data%check_status(solver_status)

  end subroutine pmc_phlex_core_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the mechanism
  integer(kind=i_kind) function pmc_phlex_core_pack_size(this) &
                  result (pack_size)

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    
    integer(kind=i_kind) :: i_mech

    pack_size =  pmc_mpi_pack_size_integer(size(this%mechanism))
    do i_mech = 1, size(this%mechanism)
      pack_size = pack_size + this%mechanism(i_mech)%pack_size()
    end do

  end function pmc_phlex_core_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_phlex_core_bin_pack(this, buffer, pos)

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: i_mech, prev_position

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, size(this%mechanism))
    do i_mech = 1, size(this%mechanism)
      call this%mechanism(i_mech)%bin_pack(buffer, pos)
    end do
    call assert(184050835, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_phlex_core_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine pmc_phlex_core_bin_unpack(this, buffer, pos)

    !> Chemical model
    class(phlex_core_t), intent(inout) :: this
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

  end subroutine pmc_phlex_core_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the time derivative f(t,y)
  subroutine pmc_phlex_core_calc_derivative(curr_time, deriv, &
                  phlex_core_c_ptr, phlex_state_c_ptr)

    use iso_c_binding

    !> Current solver time (s)
    real(kind=dp), intent(in) :: curr_time
    !> Time derivative to calculate
    real(kind=dp), intent(inout), pointer :: deriv(:)
    !> Pointer to model data
    type(c_ptr), intent(in) :: phlex_core_c_ptr
    !> Pointer to model state
    type(c_ptr), intent(in) :: phlex_state_c_ptr

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    integer(kind=i_kind) :: i_mech

    call c_f_pointer(phlex_core_c_ptr, phlex_core)
    call c_f_pointer(phlex_state_c_ptr, phlex_state)

    ! Generate a unique state id
    call phlex_state%reset_id()

    ! Calculate f(t,y)
    do i_mech=1, size(phlex_core%mechanism)
      call phlex_core%mechanism(i_mech)%get_func_contrib(phlex_state, deriv)
    end do

  end subroutine pmc_phlex_core_calc_derivative

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the Jacobian matrix J(t,y)
  subroutine pmc_phlex_core_calc_jacobian(curr_time, jac, &
                  phlex_core_c_ptr, phlex_state_c_ptr)

    use iso_c_binding

    !> Current solver time (s)
    real(kind=dp), intent(in) :: curr_time
    !> Jacobian matrix to calculate
    real(kind=dp), intent(inout), pointer :: jac(:,:)
    !> Pointer to model data
    type(c_ptr), intent(in) :: phlex_core_c_ptr
    !> Pointer to model state
    type(c_ptr), intent(in) :: phlex_state_c_ptr

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    integer(kind=i_kind) :: i_mech

    call c_f_pointer(phlex_core_c_ptr, phlex_core)
    call c_f_pointer(phlex_state_c_ptr, phlex_state)

    ! Generate a unique state id
    call phlex_state%reset_id()

    ! Calculate J(t,y)
    do i_mech=1, size(phlex_core%mechanism)
      call phlex_core%mechanism(i_mech)%get_jac_contrib(phlex_state, jac)
    end do

  end subroutine pmc_phlex_core_calc_jacobian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_phlex_core
