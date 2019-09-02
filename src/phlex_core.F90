! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_phlex_core module.

!> \page phlex_chem Phlexible Module for Chemistry
!!
!! The Phlexible Module for Chemistry (\ref phlex_chem "phlex-chem") is a
!! module within PartMC designed to provide a flexible framework for
!! incorporating chemical mechanisms into PartMC or another host model. In
!! general, the \ref phlex_chem "phlex-chem" module solves one or more
!! \ref phlex_mechanism "mechanisms" composed of a set of \ref phlex_rxn
!! "reactions" over a time-step specified by the host model. \ref phlex_rxn
!! "Reactions" can take place in the gas phase, in one of several \ref
!! phlex_aero_phase "aerosol phases", or across an interface between phases
!! (gas or aerosol). The \ref phlex_chem "phlex-chem" module is designed to
!! work with any \ref phlex_aero_rep "aerosol representation" used by the
!! host model (e.g., binned, modal, single particle) by abstracting the
!! chemistry from the \ref phlex_aero_rep "aerosol representation" through the
!! use of custom extending types of the abstract
!! \c pmc_aero_rep_data::aero_rep_data_t type that implement a set of
!! \ref phlex_aero_phase "aerosol phases" based on the configuration of the
!! host model. A set of \ref phlex_sub_model "sub-models" may also be included
!! to calculate parameters needed by \ref phlex_rxn "reactions" during solving.
!!
!! The \ref phlex_chem "phlex-chem" module uses \ref ss_json
!! "json input files" to load \ref input_format_species "chemical species",
!! \ref input_format_mechanism "mechanisms", \ref input_format_aero_phase
!! "aerosol phases", \ref input_format_aero_rep "aerosol representations",
!! and \ref input_format_sub_model "sub-models" at runtime. This allows a user
!! to modify any of these data without recompiling the model, permits host
!! models to choose which mechanisms to solve based on model conditions, and
!! allows multiple mechanisms to be solved simultaneously.
!!
!! # Phlex-Chem Input Classes #
!!
!!  - \subpage phlex_aero_phase "Aerosol Phases"
!!  - \subpage phlex_aero_rep "Aerosol Representations"
!!  - \subpage phlex_species "Chemical Species"
!!  - \subpage phlex_mechanism "Mechanisms"
!!  - \subpage phlex_rxn "Reactions"
!!  - \subpage phlex_sub_model "Sub-Models"
!!
!! # Usage #
!!
!! ## Compiling ##
!!
!! To include the \ref phlex_chem "phlex-chem" module in the PartMC library,
!! set the ccmake flags \c ENABLE_JSON and \c ENABLE_SUNDIALS to \c ON.
!! (<a href="http://www.llnl.gov/casc/sundials/">SUNDIALS</a> and
!! <a href="https://github.com/jacobwilliams/json-fortran">json-fortran</a>
!! must be installed).
!!
!! ## Input files ##
!!
!! The \ref phlex_chem "phlex-chem" module uses two types of input files:
!!
!!  - \subpage input_format_phlex_file_list "File List" A \c json file
!!             containing a list of \ref phlex_chem "phlex-chem" configuration
!!             file names.
!!  - \subpage input_format_phlex_config "Configuration File" One or more
!!             \c json files containing all the \ref phlex_chem "phlex-chem"
!!             configuration data.
!!
!! To initialize the \ref phlex_chem "phlex-chem" module, the path to the
!! \ref input_format_phlex_file_list "file list" must be passed to the
!! \ref pmc_phlex_core::phlex_core_t constructor. The method by which this is
!! done depends on the host model configuration.
!!
!! ## PartMC scenarios ##
!!
!! Using \ref phlex_chem "phlex-chem" in a PartMC scenario requires modifying
!! the \ref input_format "spec file" and providing a \ref
!! input_format_phlex_file_list "phlex-chem file list" file and one or more
!! \ref input_format_phlex_config "phlex-chem configuration" files that
!! describe the \ref phlex_species "chemical species", \ref phlex_mechanism
!! "mechanism(s)", \ref phlex_aero_phase "aerosol phase(s)", \ref
!! phlex_aero_rep "aerosol representation", and \ref phlex_sub_model
!! "sub-model(s)". A description of the input files required for a PartMC run
!! can be found \ref input_format "here".
!!
!! ## Phlex-chem in another host model ##
!!
!! Incorporating the \ref phlex_chem "phlex-chem" module into another host
!! model can be done in the following steps:
!!
!! TODO: Finish
!!

!> The phlex_core_t structure and associated subroutines.
module pmc_phlex_core

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_chem_spec_data
  use pmc_constants,                  only : i_kind, dp
  use pmc_env_state
  use pmc_mechanism_data
  use pmc_mpi
  use pmc_phlex_solver_data
  use pmc_phlex_state
  use pmc_rxn_data
  use pmc_solver_stats
  use pmc_sub_model_data
  use pmc_sub_model_factory
  use pmc_util,                       only : die_msg, string_t

  implicit none
  private

  public :: phlex_core_t

  !> Part-MC model data
  !!
  !! Contains all time-invariant data for a Part-MC model run.
  type :: phlex_core_t
  private
    !> Chemical mechanisms
    !! FIXME set up an iterator for external modules to use and
    !! make mechanisms private
    type(mechanism_data_ptr), pointer, public :: mechanism(:) => null()
    !> Chemical species data
    type(chem_spec_data_t), pointer :: chem_spec_data => null()
    !> Sub models
    type(sub_model_data_ptr), pointer :: sub_model(:) => null()
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer :: aero_rep(:) => null()
    !> Aerosol phases
    type(aero_phase_data_ptr), pointer :: aero_phase(:) => null()
    !> Size of the state array per grid cell
    integer(kind=i_kind) :: size_state_per_cell
    !> Number of cells to compute
    integer(kind=i_kind) :: n_cells = 1
    !> Initial state values
    real(kind=dp), allocatable :: init_state(:)
    !> Flag to split gas- and aerosol-phase reactions
    !! (for large aerosol representations, like single-particle)
    logical :: split_gas_aero = .false.
    !> Relative integration tolerance
    real(kind=dp) :: rel_tol = 0.0
    ! Absolute integration tolerances
    ! (Values for non-solver species will be ignored)
    real(kind=dp), allocatable :: abs_tol(:)
    ! Variable types
    integer(kind=i_kind), allocatable :: var_type(:)
    !> Solver data (gas-phase reactions)
    type(phlex_solver_data_t), pointer, public :: solver_data_gas => null()
    !> Solver data (aerosol-phase reactions)
    type(phlex_solver_data_t), pointer, public :: solver_data_aero => null()
    !> Solver data (mixed gas- and aerosol-phase reactions)
    type(phlex_solver_data_t), pointer, public :: solver_data_gas_aero => null()
    !> Flag indicating the model data has been initialized
    logical :: core_is_initialized = .false.
    !> Flag indicating the solver has been initialized
    logical :: solver_is_initialized = .false.
  contains
    !> Load a set of configuration files
    procedure :: load_files
    !> Load model data from a configuration file
    procedure :: load
    !> Initialize the model
    procedure :: initialize
    !> Indicate whether the core has been initialized
    procedure :: is_initialized
    !> Get a pointer to an aerosol phase by name
    procedure :: get_aero_phase
    !> Get a pointer to an aerosol representation by name
    procedure :: get_aero_rep
    !> Get a pointer to the set of chemical species
    procedure :: get_chem_spec_data
    !> Get a pointer to a mechanism by name
    procedure :: get_mechanism
    !> Get a pointer to a sub-model by name
    procedure :: get_sub_model
    !> Get a new model state variable
    procedure :: new_state
    !> Get the size of the state array
    procedure :: state_size
    !> Get the size of the state array per grid cell
    procedure :: state_size_per_cell
    !> Initialize the solver
    procedure :: solver_initialize
    !> Free the solver
    procedure :: free_solver
    !> Update aerosol representation data
    procedure :: update_aero_rep_data
    !> Update reaction data
    procedure :: update_rxn_data
    !> Update sub-model data
    procedure :: update_sub_model_data
    !> Run the chemical mechanisms
    procedure :: solve
    !> Determine the number of bytes required to pack the variable
    procedure :: pack_size
    !> Pack the given variable into a buffer, advancing position
    procedure :: bin_pack
    !> Unpack the given variable from a buffer, advancing position
    procedure :: bin_unpack
    !> Print the core data
    procedure :: print => do_print
    !> Finalize the core
    final :: finalize

    ! Private functions
    !> Add an aerosol phase to the model
    procedure, private :: add_aero_phase
    !> Add an aerosol representation to the model
    procedure, private :: add_aero_rep
    !> Add a mechanism to the model
    procedure, private :: add_mechanism
    !> Add a sub-model to the model
    procedure, private :: add_sub_model
  end type phlex_core_t

  !> Constructor for phlex_core_t
  interface phlex_core_t
    procedure :: constructor
  end interface phlex_core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for phlex_core_t
  function constructor(input_file_path, n_cells) result(new_obj)

    !> A new set of model parameters
    type(phlex_core_t), pointer :: new_obj
    !> Part-MC input file paths
    character(len=*), intent(in), optional :: input_file_path
    !> Num cells to compute simulatenously
    integer(kind=i_kind), optional :: n_cells

    allocate(new_obj)
    allocate(new_obj%mechanism(0))
    new_obj%chem_spec_data => chem_spec_data_t()
    allocate(new_obj%aero_phase(0))
    allocate(new_obj%aero_rep(0))
    allocate(new_obj%sub_model(0))

    if (present(n_cells)) then
      new_obj%n_cells=n_cells
    end if

    if (present(input_file_path)) then
      call new_obj%load_files(trim(input_file_path))
    end if

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_phlex_file_list Input File Format: Phlex-Chem Configuration File List
  !!
  !! A list of files containing configuration data for the \ref phlex_chem
  !! "Phlexible Module for Chemistry". The file should be in \c json format
  !! and the general structure should be the following:
  !! \code{.json}
  !! { "pmc-files" : [
  !!   "file_one.json",
  !!   "some_dir/file_two.json",
  !!   ...
  !! ]}
  !! \endcode
  !! The file should contain a single key-value pair named \b pmc-files whose
  !! value is an array of \b strings with paths to the set of \ref
  !! input_format_phlex_config "configuration" files to load. Input files
  !! should be in \c json format.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a set of model data files.
  !!
  !! See \ref input_format_phlex_file_list for the input file format.
  subroutine load_files(this, input_file_path)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Part-MC input file paths
    character(len=*), intent(in) :: input_file_path

#ifdef PMC_USE_JSON
    type(json_core), target :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next

    logical(kind=json_lk) :: found, valid
    character(kind=json_ck, len=:), allocatable :: unicode_str_val
    character(kind=json_ck, len=:), allocatable :: json_err_msg
    integer(kind=json_ik) :: i_file, num_files
    type(string_t), allocatable :: file_list(:)
    logical :: file_exists

    ! load the file containing the paths to the configuration files
    call j_file%initialize()
    call j_file%get_core(json)
    call assert_msg(600888426, trim(input_file_path).ne."", &
            "Received empty string for file path")
    inquire( file=trim(input_file_path), exist=file_exists )
    call assert_msg(433777575, file_exists, "Cannot find file: "//&
            trim(input_file_path))
    call j_file%load_file(filename = trim(input_file_path))

    ! get the set of configuration file names
    call j_file%get('pmc-files(1)', j_obj, found)
    call assert_msg(405149265, found, &
            "Could not find pmc-files object in input file: "// &
            input_file_path)
    call json%validate(j_obj, valid, json_err_msg)
    if(.not.valid) then
      call die_msg(959537834, "Bad JSON format in file '"// &
                   trim(input_file_path)//"': "//trim(json_err_msg))
    end if
    call j_file%info('pmc-files', n_children = num_files)
    call assert_msg(411804027, num_files.gt.0, &
            "No file names were found in "//input_file_path)

    ! allocate space for the configurtaion file names
    allocate(file_list(num_files))

    ! cycle through the list of file names, adding each to the list
    j_next => null()
    i_file = 1
    do while (associated(j_obj))
      call json%get(j_obj, unicode_str_val)
      file_list(i_file)%string = unicode_str_val
      i_file = i_file + 1
      call json%get_next(j_obj, j_next)
      j_obj => j_next
    end do

    ! free the json file
    call j_file%destroy()

    ! load all the configuration files
    call this%load(file_list)

#else
    call warn_msg(171627969, "No support for input files.");
#endif

  end subroutine load_files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_phlex_config Input File Format: Phlex-Chem Configuration Data
  !!
  !! Configuration data for the
  !! \ref phlex_chem "Phlexible Module for Chemistry". The files are in
  !! \c json format and their general structure should be the following:
  !! \code{.json}
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
  !! \endcode
  !! Each input file should contain exactly one \c json object with a single
  !! key-value pair \b pmc-data whose value is an array of \c json objects.
  !! Additional top-level key-value pairs will be ignored. Each of the \c json
  !! objects in the \b pmc-data array must contain a key-value pair \b type
  !! whose value is a string referencing a valid PartMC object.
  !!
  !! The valid values for \b type are:
  !!
  !!   - \subpage input_format_mechanism "MECHANISM"
  !!   - \subpage input_format_species "CHEM_SPEC"
  !!   - \subpage input_format_aero_phase "AERO_PHASE"
  !!   - \subpage input_format_aero_rep "AERO_REP_*"
  !!   - \subpage input_format_sub_model "SUB_MODEL_*"
  !!
  !! The arrangement of objects within the \b pmc-data array and between input
  !! files is arbitrary. Additionally, some objects, such as \ref
  !! input_format_species "chemical species" and \ref input_format_mechanism
  !! "mechanisms" may be split into multiple objects within the \b pmc-data
  !! array and/or between files, and will be combined based on their unique
  !! name. This flexibility is provided so that the chemical mechanism data
  !! can be organized in a way that makes sense to the designer of the
  !! mechanism. For example, files could be split based on species source
  !! (biogenic, fossil fuel, etc.) or based on properties (molecular weight,
  !! density, etc.) or any combination of criteria. However, if a single
  !! property of an object (e.g., the molecular weight of a chemical species)
  !! is set in more than one location, this will cause an error.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load model data from input files
  !!
  !! See \ref input_format_phlex_config for the input file format.
  subroutine load(this, input_file_path)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Part-MC input file paths
    type(string_t), allocatable, intent(in) :: input_file_path(:)

    integer(kind=i_kind) :: i_file, cell
#ifdef PMC_USE_JSON
    type(json_core), pointer :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next

    logical(kind=json_lk) :: valid
    character(kind=json_ck, len=:), allocatable :: unicode_str_val
    character(kind=json_ck, len=:), allocatable :: json_err_msg
    character(len=:), allocatable :: str_val
    real(kind=json_rk) :: real_val
    logical :: file_exists, found

    ! mechansim
    type(mechanism_data_t), pointer :: mech_ptr

    ! sub models
    type(sub_model_data_ptr), pointer :: new_sub_model(:)
    type(sub_model_factory_t) :: sub_model_factory
    type(sub_model_data_ptr) :: sub_model_ptr
    class(sub_model_data_t), pointer :: existing_sub_model_ptr
    logical :: sub_model_placed
    integer(kind=i_kind) :: i_sub_model, j_sub_model

    ! aerosol representations
    type(aero_rep_data_ptr), pointer :: new_aero_rep(:)
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_data_ptr) :: aero_rep_ptr
    class(aero_rep_data_t), pointer :: existing_aero_rep_ptr

    ! aerosol phases
    type(aero_phase_data_ptr), pointer :: new_aero_phase(:)
    class(aero_phase_data_t), pointer :: aero_phase, existing_aero_phase

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! cycle through the list of configuration files !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    j_obj => null()
    j_next => null()
    allocate(json)
    do i_file = 1, size(input_file_path)

      ! load the configuration file
      call j_file%initialize()
      call j_file%get_core(json)
      call assert_msg(366175417, allocated(input_file_path(i_file)%string), &
              "Received non-allocated string for file path")
      call assert_msg(936390222, trim(input_file_path(i_file)%string).ne."", &
              "Received empty string for file path")
      inquire( file=input_file_path(i_file)%string, exist=file_exists )
      call assert_msg(910660557, file_exists, "Cannot file file: "// &
              input_file_path(i_file)%string)
      call j_file%load_file(filename = input_file_path(i_file)%string)

      ! get the phlex-chem objects
      call j_file%get('pmc-data(1)', j_obj)
      call json%validate(j_obj, valid, json_err_msg)
      if (.not.valid) then
        call die_msg(560270545, "Bad JSON format in file '"// &
                     trim(input_file_path(i_file)%string)//"': "// &
                     trim(json_err_msg))
      end if
      do while (associated(j_obj))

        ! get the object type and load data into the appropriate phlex-chem
        ! derived type
        call json%get(j_obj, 'type', unicode_str_val, found)
        call assert_msg(689470331, found, &
                "Missing type in json input file "// &
                input_file_path(i_file)%string)
        str_val = unicode_str_val

        !!!!!!!!!!!!!!!!!!!!!!!!
        !!! load a mechanism !!!
        !!!!!!!!!!!!!!!!!!!!!!!!
        if (str_val.eq.'MECHANISM') then
          call json%get(j_obj, 'name', unicode_str_val, found)
          call assert_msg(822680732, found, &
                  "Missing mechanism name in file "// &
                  input_file_path(i_file)%string)
          str_val = unicode_str_val

          ! if a mechanism with the same name already exists, add data to it
          ! otherwise, add a new mechanism
          if (.not.this%get_mechanism(str_val, mech_ptr)) then
            call this%add_mechanism(str_val)
            call assert(105816325, this%get_mechanism(str_val, mech_ptr))
          end if
          call mech_ptr%load(json, j_obj)

        ! load a chemical species
        else if (str_val.eq.'CHEM_SPEC') then
          call this%chem_spec_data%load(json, j_obj)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! load an aerosol representation !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (str_val(1:8).eq.'AERO_REP') then
          aero_rep_ptr%val => aero_rep_factory%load(json, j_obj)
          str_val = aero_rep_ptr%val%name()

          ! if an aerosol representation with the same name already exists,
          ! add data to it. otherwise, add a new aerosol representation
          if (this%get_aero_rep(str_val, existing_aero_rep_ptr)) then
            deallocate(aero_rep_ptr%val)
            call existing_aero_rep_ptr%load(json, j_obj)
          else
            allocate(new_aero_rep(size(this%aero_rep)+1))
            new_aero_rep(1:size(this%aero_rep)) = &
                    this%aero_rep(1:size(this%aero_rep))
            new_aero_rep(size(new_aero_rep))%val => aero_rep_ptr%val
            call this%aero_rep(:)%dereference()
            deallocate(this%aero_rep)
            this%aero_rep => new_aero_rep
            call aero_rep_ptr%dereference()
          end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! load an aerosol phase !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (str_val.eq.'AERO_PHASE') then
          aero_phase => aero_phase_data_t()
          call aero_phase%load(json, j_obj)
          str_val = aero_phase%name()

          ! if an aerosol phase with the same name already exists, add data to
          ! it. otherwise, add a new aerosol phase
          if (this%get_aero_phase(str_val, existing_aero_phase)) then
            deallocate(aero_phase)
            call existing_aero_phase%load(json, j_obj)
          else
            allocate(new_aero_phase(size(this%aero_phase)+1))
            new_aero_phase(1:size(this%aero_phase)) = &
                    this%aero_phase(1:size(this%aero_phase))
            new_aero_phase(size(new_aero_phase))%val => aero_phase
            call this%aero_phase(:)%dereference()
            deallocate(this%aero_phase)
            this%aero_phase => new_aero_phase
          end if

        !!!!!!!!!!!!!!!!!!!!!!!!
        !!! load a sub-model !!!
        !!!!!!!!!!!!!!!!!!!!!!!!
        else if (str_val(1:9).eq.'SUB_MODEL') then
          sub_model_ptr%val => sub_model_factory%load(json, j_obj)
          str_val = sub_model_ptr%val%name()

          ! if an sub-model with the same name already exists, add data to it.
          ! otherwise, add a new sub-model
          if (this%get_sub_model(str_val, existing_sub_model_ptr)) then
            deallocate(sub_model_ptr%val)
            call existing_sub_model_ptr%load(json, j_obj)
          else
            sub_model_placed = .false.
            allocate(new_sub_model(size(this%sub_model)+1))
            j_sub_model = 1
            do i_sub_model = 1, size(this%sub_model)
              if (.not.sub_model_placed .and. &
                  sub_model_ptr%val%priority() < &
                  this%sub_model(i_sub_model)%val%priority()) then
                    sub_model_placed = .true.
                    new_sub_model(j_sub_model)%val => sub_model_ptr%val
                    j_sub_model = j_sub_model + 1
              end if
              new_sub_model(j_sub_model) = &
                      this%sub_model(i_sub_model)
              j_sub_model = j_sub_model + 1
            end do
            if (.not.sub_model_placed) then
              new_sub_model(j_sub_model)%val => sub_model_ptr%val
            end if
            call this%sub_model(:)%dereference()
            deallocate(this%sub_model)
            this%sub_model => new_sub_model
            call sub_model_ptr%dereference()
          end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! set the relative tolerance for the model !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (str_val.eq.'RELATIVE_TOLERANCE') then
          call json%get(j_obj, 'value', real_val, found)
          call assert_msg(761842352, found, &
                  "Missing value for relative tolerance")
          call assert_msg(162564706, real_val.gt.0.0.and.real_val.lt.1.0, &
                  "Invalid relative tolerance: "// &
                  trim(to_string(real(real_val, kind=dp))))
          this%rel_tol = real(real_val, kind=dp)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! set whether to solve gas and aerosol phases separately !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else if (str_val.eq.'SPLIT_GAS_AERO') then
          this%split_gas_aero = .true.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! fail on invalid object type !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        else
          call die_msg(448039776, &
                  "Received invalid json input object type: "//str_val)
        end if

        ! get the next object
        j_next => j_obj
        call json%get_next(j_next, j_obj)
      end do

      ! reset the json objects
      call j_file%destroy()
      call json%destroy()
    end do

    ! free the json core
    deallocate(json)
#else
    call warn_msg(350136328, "No support for input files.");
#endif

  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the model data
  subroutine initialize(this)

    use iso_c_binding

    !> Model data
    class(phlex_core_t), target, intent(inout) :: this

    ! Indices for iteration
    integer(kind=i_kind) :: i_mech, i_phase, i_aero_rep, i_sub_model
    integer(kind=i_kind) :: i_state_var, i_spec, i_cell

    ! Variables for setting initial state values
    class(aero_rep_data_t), pointer :: rep
    integer(kind=i_kind) :: i_state_elem, i_name

    ! Species name for looking up properties
    character(len=:), allocatable :: spec_name
    ! Gas-phase species names
    type(string_t), allocatable :: gas_spec_names(:)
    ! Aerosol species
    type(string_t), allocatable :: unique_names(:)

    ! make sure the core has not already been initialized
    call assert_msg(157261665, .not.this%core_is_initialized, &
            "Attempting to initialize a phlex_core_t object twice.")

    ! Initialize the species database
    call this%chem_spec_data%initialize()

    ! Get the next index on the state array after the gas-phase species
    i_state_var = this%chem_spec_data%size(spec_phase=CHEM_SPEC_GAS_PHASE) + 1

    ! Initialize the aerosol phases
    do i_phase = 1, size(this%aero_phase)
      call assert(254948966, associated(this%aero_phase(i_phase)%val))
      call this%aero_phase(i_phase)%val%initialize(this%chem_spec_data)
    end do

    ! Initialize the aerosol representations
    do i_aero_rep = 1, size(this%aero_rep)
      call assert(251590193, associated(this%aero_rep(i_aero_rep)%val))
      call this%aero_rep(i_aero_rep)%val%initialize(this%aero_phase, &
              i_state_var)
      i_state_var = i_state_var + this%aero_rep(i_aero_rep)%val%size()
    end do

    ! Initialize the sub-models
    do i_sub_model = 1, size(this%sub_model)
      call assert(565644925, associated(this%sub_model(i_sub_model)%val))
      call this%sub_model(i_sub_model)%val%initialize(this%aero_rep, &
                this%aero_phase, this%chem_spec_data)
    end do

    ! Set the size of the state array per grid cell
    this%size_state_per_cell = i_state_var - 1

    ! Initialize the mechanisms
    do i_mech = 1, size(this%mechanism)
      call this%mechanism(i_mech)%val%initialize(this%chem_spec_data, &
              this%aero_rep, this%n_cells)
    end do

    ! Allocate space for the variable types and absolute tolerances
    allocate(this%abs_tol(this%size_state_per_cell))
    allocate(this%var_type(this%size_state_per_cell))

    ! Start at the first state array element
    i_state_var = 0

    ! Add gas-phase species variable types and absolute tolerances
    gas_spec_names = &
            this%chem_spec_data%get_spec_names(spec_phase = &
            CHEM_SPEC_GAS_PHASE)
    do i_spec = 1, size(gas_spec_names)
      i_state_var = i_state_var + 1
      call assert(716433999, &
              this%chem_spec_data%get_abs_tol(gas_spec_names(i_spec)%string, &
              this%abs_tol(i_state_var)))
      call assert(888496437, &
              this%chem_spec_data%get_type(gas_spec_names(i_spec)%string, &
              this%var_type(i_state_var)))
    end do

    ! Loop through aerosol representations
    do i_aero_rep = 1, size(this%aero_rep)
      ! Set aerosol-phase species variable types and absolute tolerances
      ! TODO Move this to the aerosol representations, so they have control
      ! of their portion on the state array and what is stored there
      call assert(666823548, associated(this%aero_rep(i_aero_rep)%val))
      unique_names = this%aero_rep(i_aero_rep)%val%unique_names()
      do i_spec = 1, this%aero_rep(i_aero_rep)%val%size()
        i_state_var = i_state_var + 1
        spec_name = this%aero_rep(i_aero_rep)%val%spec_name( &
                  unique_names(i_spec)%string)
        call assert(709716453, &
                this%chem_spec_data%get_abs_tol(spec_name, &
                this%abs_tol(i_state_var)))
        call assert(257084300, &
                this%chem_spec_data%get_type(spec_name, &
                this%var_type(i_state_var)))
      end do
      deallocate(unique_names)
    end do

    ! Make sure absolute tolerance and variable type arrays are completely
    ! filled
    call assert_msg(501609702, i_state_var.eq.this%size_state_per_cell, &
            "Internal error. Filled "//trim(to_string(i_state_var))// &
            " of "//trim(to_string(this%size_state_per_cell))// &
            " elements of absolute tolerance and variable type arrays")

    this%core_is_initialized = .true.

    ! Set the initial state values
    allocate(this%init_state(this%size_state_per_cell * this%n_cells))

    ! Set species concentrations to zero
    this%init_state(:) = 0.0

    ! Set activity coefficients to 1.0
      do i_cell = 0, this%n_cells - 1
        do i_aero_rep = 1, size(this%aero_rep)

          rep => this%aero_rep(i_aero_rep)%val

          ! Get the ion pairs for which activity coefficients can be calculated
          unique_names = rep%unique_names(tracer_type = CHEM_SPEC_ACTIVITY_COEFF)

          ! Set the activity coefficients to 1.0 as default
          do i_name = 1, size(unique_names)
            i_state_elem = rep%spec_state_id(unique_names(i_name)%string)
            this%init_state(i_state_elem + i_cell * this%size_state_per_cell) = &
                    real(1.0d0, kind=dp)
          end do

          deallocate(unique_names)

        end do
      end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Inidicate whether the core has been initialized
  logical function is_initialized(this)

    !> Model data
    class(phlex_core_t), intent(in) :: this

    is_initialized = this%core_is_initialized

  end function is_initialized

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a pointer to an aerosol phase by name
  logical function get_aero_phase(this, aero_phase_name, aero_phase) &
            result (found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Aerosol phase name to search for
    character(len=:), allocatable, intent(in) :: aero_phase_name
    !> Pointer to the aerosol phase
    class(aero_phase_data_t), pointer, intent(out) :: aero_phase

    integer(kind=i_kind) :: i_aero_phase

    found = .false.
    aero_phase => null()
    if (.not.associated(this%aero_phase)) return
    do i_aero_phase = 1, size(this%aero_phase)
      if (this%aero_phase(i_aero_phase)%val%name().eq.aero_phase_name) then
        found = .true.
        aero_phase => this%aero_phase(i_aero_phase)%val
        return
      end if
    end do

  end function get_aero_phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a pointer to an aerosol representation by name
  logical function get_aero_rep(this, aero_rep_name, aero_rep) result (found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Aerosol representation name to search for
    character(len=:), allocatable, intent(in) :: aero_rep_name
    !> Aerosol representation
    class(aero_rep_data_t), pointer, intent(out) :: aero_rep

    integer(kind=i_kind) :: i_aero_rep

    found = .false.
    aero_rep => null()
    if (.not.associated(this%aero_rep)) return
    do i_aero_rep = 1, size(this%aero_rep)
      if (this%aero_rep(i_aero_rep)%val%name().eq.aero_rep_name) then
        aero_rep => this%aero_rep(i_aero_rep)%val
        found = .true.
        return
      end if
    end do

  end function get_aero_rep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a pointer to the chemical species data
  logical function get_chem_spec_data(this, chem_spec_data) result (found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Pointer to the chemical species data
    type(chem_spec_data_t), pointer :: chem_spec_data

    found = .false.
    chem_spec_data => this%chem_spec_data
    if (associated(chem_spec_data)) found = .true.

  end function get_chem_spec_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a pointer to a mechanism by name
  logical function get_mechanism(this, mech_name, mechanism) result (found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Mechanism name to search for
    character(len=:), allocatable, intent(in) :: mech_name
    !> Pointer to the mechanism
    type(mechanism_data_t), pointer, intent(out) :: mechanism

    integer(kind=i_kind) :: i_mech

    found = .false.
    mechanism => null()
    if (.not.associated(this%mechanism)) return
    do i_mech = 1, size(this%mechanism)
      if (this%mechanism(i_mech)%val%name().eq.mech_name) then
        found = .true.
        mechanism => this%mechanism(i_mech)%val
        return
      end if
    end do

  end function get_mechanism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find an sub-model by name
  logical function get_sub_model(this, sub_model_name, sub_model) &
            result (found)

    !> Model data
    class(phlex_core_t), intent(in) :: this
    !> Sub model name to search for
    character(len=:), allocatable, intent(in) :: sub_model_name
    !> Sub model
    class(sub_model_data_t), pointer, intent(out) :: sub_model

    integer(kind=i_kind) :: i_sub_model

    found = .false.
    sub_model => null()
    if (.not.associated(this%sub_model)) return
    do i_sub_model = 1, size(this%sub_model)
      if (this%sub_model(i_sub_model)%val%name().eq.sub_model_name) then
        sub_model => this%sub_model(i_sub_model)%val
        found = .true.
        return
      end if
    end do

  end function get_sub_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a model state variable based on the this set of model data
  function new_state(this, env_state)

    !> New model state
    type(phlex_state_t), pointer :: new_state
    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> Environmental state array
    !! (one element per grid cell to solve simultaneously)
    type(env_state_t), optional, target, intent(in) :: env_state

    ! Initialize phlex_state
    new_state => phlex_state_t(env_state,this%n_cells)

    ! Set up the state variable array
    allocate(new_state%state_var, source=this%init_state)

  end function new_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the state array
  function state_size(this)

    !> State size
    integer(kind=i_kind) :: state_size
    !> Chemical model
    class(phlex_core_t), intent(in) :: this

    call assert_msg(629102639, allocated(this%init_state), &
                    "Trying to get the size of the state array before "// &
                    "initializing the phlex_core")

    state_size = size(this%init_state)

  end function state_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the state array for each grid cell
  function state_size_per_cell(this) result(state_size)

    !> State size
    integer(kind=i_kind) :: state_size
    !> Chemical model
    class(phlex_core_t), intent(in) :: this

    call assert_msg(175845182, allocated(this%init_state), &
                    "Trying to get the size of the state array before "// &
                    "initializing the phlex_core")

    state_size = this%size_state_per_cell

  end function state_size_per_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the solver
  subroutine solver_initialize(this)

    !> Chemical model
    class(phlex_core_t), intent(inout) :: this

    call assert_msg(662920365, .not.this%solver_is_initialized, &
            "Attempting to initialize the solver twice.")

    ! Set up either two solvers (gas and aerosol) or one solver (combined)
    if (this%split_gas_aero) then

      ! Create the new solver data objects
      this%solver_data_gas => phlex_solver_data_t()
      this%solver_data_aero => phlex_solver_data_t()

      ! Set custom relative integration tolerance, if present
      if (this%rel_tol.ne.real(0.0, kind=dp)) then
        this%solver_data_gas%rel_tol = this%rel_tol
        this%solver_data_aero%rel_tol = this%rel_tol
      end if

      ! Initialize the solvers
      call this%solver_data_gas%initialize( &
                this%var_type,   & ! State array variable types
                this%abs_tol,    & ! Absolute tolerances for each state var
                this%mechanism,  & ! Pointer to the mechanisms
                this%aero_phase, & ! Pointer to the aerosol phases
                this%aero_rep,   & ! Pointer to the aerosol representations
                this%sub_model,  & ! Pointer to the sub-models
                GAS_RXN,         & ! Reaction phase
                this%n_cells   & ! # of cells computed simultaneosly
                )
      call this%solver_data_aero%initialize( &
                this%var_type,   & ! State array variable types
                this%abs_tol,    & ! Absolute tolerances for each state var
                this%mechanism,  & ! Pointer to the mechanisms
                this%aero_phase, & ! Pointer to the aerosol phases
                this%aero_rep,   & ! Pointer to the aerosol representations
                this%sub_model,  & ! Pointer to the sub-models
                AERO_RXN,        & ! Reaction phase
                this%n_cells   & ! # of cells computed simultaneosly
                )
    else

      ! Create a new solver data object
      this%solver_data_gas_aero => phlex_solver_data_t()

      ! Set custom relative integration tolerance, if present
      if (this%rel_tol.ne.0.0) then
        this%solver_data_gas_aero%rel_tol = this%rel_tol
      end if

      ! Initialize the solver
      call this%solver_data_gas_aero%initialize( &
                this%var_type,   & ! State array variable types
                this%abs_tol,    & ! Absolute tolerances for each state var
                this%mechanism,  & ! Pointer to the mechanisms
                this%aero_phase, & ! Pointer to the aerosol phases
                this%aero_rep,   & ! Pointer to the aerosol representations
                this%sub_model,  & ! Pointer to the sub-models
                GAS_AERO_RXN,    & ! Reaction phase
                this%n_cells   & ! # of cells computed simultaneosly
                )

    end if

    this%solver_is_initialized = .true.

  end subroutine solver_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free the solver memory
  subroutine free_solver(this)

    !> Phlex-core
    class(phlex_core_t), intent(inout) :: this

    if( associated( this%solver_data_gas ) )  &
        deallocate( this%solver_data_gas )
    if( associated( this%solver_data_aero ) ) &
        deallocate( this%solver_data_aero )
    if( associated( this%solver_data_gas_aero ) ) &
        deallocate( this%solver_data_gas_aero )

  end subroutine free_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update data associated with an aerosol representation. This function
  !! should be called by an external aerosol microphysics model whenever
  !! the aerosol condensed data needs updated based on changes in, e.g.,
  !! particle size or number concentration. The update types are aerosol-
  !! representation specific.
  subroutine update_aero_rep_data(this, update_data)

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> Update data
    class(aero_rep_update_data_t), intent(in) :: update_data

    if (associated(this%solver_data_gas)) &
            call this%solver_data_gas%update_aero_rep_data(update_data)
    if (associated(this%solver_data_aero)) &
            call this%solver_data_aero%update_aero_rep_data(update_data)
    if (associated(this%solver_data_gas_aero)) &
            call this%solver_data_gas_aero%update_aero_rep_data(update_data)

  end subroutine update_aero_rep_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update data associated with a reaction. This function should be called
  !! when reaction parameters need updated from the host model. For example,
  !! this function can be called to update photolysis rates from a host
  !! model's photolysis module.
  subroutine update_rxn_data(this, update_data)

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> Update data
    class(rxn_update_data_t), intent(in) :: update_data

    if (associated(this%solver_data_gas)) &
            call this%solver_data_gas%update_rxn_data(update_data)
    if (associated(this%solver_data_aero)) &
            call this%solver_data_aero%update_rxn_data(update_data)
    if (associated(this%solver_data_gas_aero)) &
            call this%solver_data_gas_aero%update_rxn_data(update_data)

  end subroutine update_rxn_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update data associated with a sub-model. This function should be called
  !! when sub-model parameters need updated from the host model.
  subroutine update_sub_model_data(this, update_data)

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> Update data
    class(sub_model_update_data_t), intent(in) :: update_data

    if (associated(this%solver_data_gas)) &
            call this%solver_data_gas%update_sub_model_data(update_data)
    if (associated(this%solver_data_aero)) &
            call this%solver_data_aero%update_sub_model_data(update_data)
    if (associated(this%solver_data_gas_aero)) &
            call this%solver_data_gas_aero%update_sub_model_data(update_data)

  end subroutine update_sub_model_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Integrate the chemical mechanism
  subroutine solve(this, phlex_state, time_step, rxn_phase, solver_stats)

    use pmc_rxn_data
    use pmc_solver_stats
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
    !> Return solver statistics to the host model
    type(solver_stats_t), intent(inout), optional, target :: solver_stats

    ! Phase to solve
    integer(kind=i_kind) :: phase
    ! Pointer to solver data
    type(phlex_solver_data_t), pointer :: solver

    call assert_msg(593328365, this%solver_is_initialized,                   &
                    "Trying to solve system with uninitialized solver" )

    ! Get the phase(s) to solve for
    if (present(rxn_phase)) then
      phase = rxn_phase
    else
      phase = GAS_AERO_RXN
    end if

    ! Determine the solver to use
    if (phase.eq.GAS_RXN) then
        solver => this%solver_data_gas
    else if (phase.eq.AERO_RXN) then
        solver => this%solver_data_aero
    else if (phase.eq.GAS_AERO_RXN) then
        solver => this%solver_data_gas_aero
    else
      call die_msg(704896254, "Invalid rxn phase specified for chemistry "// &
              "solver: "//to_string(phase))
    end if

    ! Update the environmental state array
    ! TODO May move this into the solver functions to allow user to vary
    ! environmental parameters with time during the chemistry time step
    if (this%n_cells.eq.1) then ! Make this not necessary -> discuss with matt
      call phlex_state%update_env_state(0)
    end if

    ! Make sure the requested solver was loaded
    call assert_msg(730097030, associated(solver), "Invalid solver requested")

    ! Run the integration
    if (present(solver_stats)) then
      call solver%solve(phlex_state, real(0.0, kind=dp), time_step,          &
                        solver_stats)
    else
      call solver%solve(phlex_state, real(0.0, kind=dp), time_step)
    end if
  end subroutine solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the mechanism
  integer(kind=i_kind) function pack_size(this, comm)

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> MPI communicator
    integer, intent(in), optional :: comm

    type(aero_rep_factory_t) :: aero_rep_factory
    type(sub_model_factory_t) :: sub_model_factory
    class(aero_rep_data_t), pointer :: aero_rep
    class(sub_model_data_t), pointer :: sub_model
    integer(kind=i_kind) :: i_mech, i_phase, i_rep, i_sub_model, l_comm

#ifdef PMC_USE_MPI
    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    call assert_msg(143374295, this%core_is_initialized, &
            "Trying to get the buffer size of an uninitialized core.")

    pack_size =  pmc_mpi_pack_size_integer(size(this%mechanism),  l_comm) + &
                 pmc_mpi_pack_size_integer(size(this%aero_phase), l_comm) + &
                 pmc_mpi_pack_size_integer(size(this%aero_rep),   l_comm) + &
                 pmc_mpi_pack_size_integer(size(this%sub_model),  l_comm)
    do i_mech = 1, size(this%mechanism)
      pack_size = pack_size + this%mechanism(i_mech)%val%pack_size(l_comm)
    end do
    do i_phase = 1, size(this%aero_phase)
      pack_size = pack_size + this%aero_phase(i_phase)%val%pack_size(l_comm)
    end do
    do i_rep = 1, size(this%aero_rep)
      aero_rep => this%aero_rep(i_rep)%val
      pack_size = pack_size + aero_rep_factory%pack_size(aero_rep, l_comm)
      aero_rep => null()
    end do
    do i_sub_model = 1, size(this%sub_model)
      sub_model => this%sub_model(i_sub_model)%val
      pack_size = pacK_size + sub_model_factory%pack_size(sub_model, l_comm)
      sub_model => null()
    end do
    pack_size = pack_size + &
                pmc_mpi_pack_size_integer(this%size_state_per_cell, l_comm) + &
                pmc_mpi_pack_size_logical(this%split_gas_aero, l_comm) + &
                pmc_mpi_pack_size_real(this%rel_tol, l_comm) + &
                pmc_mpi_pack_size_real_array(this%abs_tol, l_comm) + &
                pmc_mpi_pack_size_integer_array(this%var_type, l_comm) + &
                pmc_mpi_pack_size_real_array(this%init_state, l_comm)
#else
    pack_size = 0
#endif
  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos, comm)

    !> Chemical model
    class(phlex_core_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    type(aero_rep_factory_t) :: aero_rep_factory
    type(sub_model_factory_t) :: sub_model_factory
    class(aero_rep_data_t), pointer :: aero_rep
    class(sub_model_data_t), pointer :: sub_model
    integer(kind=i_kind) :: i_mech, i_phase, i_rep, i_sub_model, &
            prev_position, l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    call assert_msg(143374295, this%core_is_initialized, &
            "Trying to pack an uninitialized core.")

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, size(this%mechanism),  l_comm)
    call pmc_mpi_pack_integer(buffer, pos, size(this%aero_phase), l_comm)
    call pmc_mpi_pack_integer(buffer, pos, size(this%aero_rep),   l_comm)
    call pmc_mpi_pack_integer(buffer, pos, size(this%sub_model),  l_comm)
    do i_mech = 1, size(this%mechanism)
      call this%mechanism(i_mech)%val%bin_pack(buffer, pos, l_comm)
    end do
    do i_phase = 1, size(this%aero_phase)
      call this%aero_phase(i_phase)%val%bin_pack(buffer, pos, l_comm)
    end do
    do i_rep = 1, size(this%aero_rep)
      aero_rep => this%aero_rep(i_rep)%val
      call aero_rep_factory%bin_pack(aero_rep, buffer, pos, l_comm)
      aero_rep => null()
    end do
    do i_sub_model = 1, size(this%sub_model)
      sub_model => this%sub_model(i_sub_model)%val
      call sub_model_factory%bin_pack(sub_model, buffer, pos, l_comm)
      sub_model => null()
    end do
    call pmc_mpi_pack_integer(buffer, pos, this%size_state_per_cell, l_comm)
    call pmc_mpi_pack_logical(buffer, pos, this%split_gas_aero, l_comm)
    call pmc_mpi_pack_real(buffer, pos, this%rel_tol, l_comm)
    call pmc_mpi_pack_real_array(buffer, pos, this%abs_tol, l_comm)
    call pmc_mpi_pack_integer_array(buffer, pos, this%var_type, l_comm)
    call pmc_mpi_pack_real_array(buffer, pos, this%init_state, l_comm)
    call assert(184050835, &
         pos - prev_position <= this%pack_size(l_comm))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value from the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos, comm)

    !> Chemical model
    class(phlex_core_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in), optional :: comm

#ifdef PMC_USE_MPI
    type(aero_rep_factory_t) :: aero_rep_factory
    type(sub_model_factory_t) :: sub_model_factory
    integer(kind=i_kind) :: i_mech, i_phase, i_rep, i_sub_model, &
            prev_position, num_mech, num_phase, num_rep, num_sub_model, &
            l_comm

    if (present(comm)) then
      l_comm = comm
    else
      l_comm = MPI_COMM_WORLD
    endif

    call finalize(this)
    this%chem_spec_data => chem_spec_data_t()

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, num_mech,      l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, num_phase,     l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, num_rep,       l_comm)
    call pmc_mpi_unpack_integer(buffer, pos, num_sub_model, l_comm)
    allocate(this%mechanism(num_mech))
    allocate(this%aero_phase(num_phase))
    allocate(this%aero_rep(num_rep))
    allocate(this%sub_model(num_sub_model))
    do i_mech = 1, num_mech
      this%mechanism(i_mech)%val => mechanism_data_t()
      call this%mechanism(i_mech)%val%bin_unpack(buffer, pos, l_comm)
    end do
    do i_phase = 1, num_phase
      this%aero_phase(i_phase)%val => aero_phase_data_t()
      call this%aero_phase(i_phase)%val%bin_unpack(buffer, pos, l_comm)
    end do
    do i_rep = 1, num_rep
      this%aero_rep(i_rep)%val => aero_rep_factory%bin_unpack(buffer, pos, l_comm)
    end do
    do i_sub_model = 1, num_sub_model
      this%sub_model(i_sub_model)%val => &
              sub_model_factory%bin_unpack(buffer, pos, l_comm)
    end do
    call pmc_mpi_unpack_integer(buffer, pos, this%size_state_per_cell, l_comm)
    call pmc_mpi_unpack_logical(buffer, pos, this%split_gas_aero, l_comm)
    call pmc_mpi_unpack_real(buffer, pos, this%rel_tol, l_comm)
    call pmc_mpi_unpack_real_array(buffer, pos, this%abs_tol, l_comm)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%var_type, l_comm)
    call pmc_mpi_unpack_real_array(buffer, pos, this%init_state, l_comm)
    this%core_is_initialized = .true.
    call assert(291557168, &
         pos - prev_position <= this%pack_size(l_comm))
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the core data
  subroutine do_print(this, file_unit, solver_data_only)

    !> Core data
    class(phlex_core_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), intent(in), optional :: file_unit
    !> Print only the solver data (can be used during runtime for debugging)
    logical, intent(in), optional :: solver_data_only

    integer(kind=i_kind) :: i_gas_spec, i_spec, j_spec, i_phase, i_aero_rep, i_mech
    integer(kind=i_kind) :: i_sub_model, i_solver_spec
    integer(kind=i_kind) :: f_unit=6
    type(string_t), allocatable :: state_names(:), rep_spec_names(:)
    logical :: sd_only = .false.

    if (present(file_unit)) f_unit = file_unit
    if (present(solver_data_only)) sd_only = solver_data_only

    write(f_unit,*) "*********************"
    write(f_unit,*) "** Phlex core data **"
    write(f_unit,*) "*********************"
    if (.not.sd_only ) then
      write(f_unit,*) "Number of grid cells to solve simultaneously: ", &
                      this%n_cells
      write(f_unit,*) "Relative integration tolerance: ", this%rel_tol
      call this%chem_spec_data%print(f_unit)
      write(f_unit,*) "*** Aerosol Phases ***"
      do i_phase=1, size(this%aero_phase)
        call this%aero_phase(i_phase)%val%print(f_unit)
      end do
      write(f_unit,*) "*** Aerosol Representations ***"
      do i_aero_rep=1, size(this%aero_rep)
        write(f_unit,*) "Aerosol representation ", i_aero_rep
        call this%aero_rep(i_aero_rep)%val%print(f_unit)
      end do
      write(f_unit,*) "*** Sub Models ***"
      do i_sub_model=1, size(this%sub_model)
        write(f_unit,*) "Sub model: ", i_sub_model
        call this%sub_model(i_sub_model)%val%print(f_unit)
      end do
      write(f_unit,*) "*** Mechanisms ***"
      write(f_unit,*) "Number of mechanisms: ", size(this%mechanism)
      do i_mech=1, size(this%mechanism)
        call this%mechanism(i_mech)%val%print(f_unit)
      end do
      write(f_unit,*) "*** State Array ***"
      write(f_unit,*) "Number of species on the state array per grid cell: ", &
                      this%size_state_per_cell
      allocate(state_names(this%size_state_per_cell))
      i_spec = 0
      do i_gas_spec = 1, &
              this%chem_spec_data%size(spec_phase=CHEM_SPEC_GAS_PHASE)
        i_spec = i_gas_spec
        state_names(i_spec)%string = &
                this%chem_spec_data%gas_state_name(i_gas_spec)
      end do
      write(f_unit,*) "Gas-phase species: ", i_spec
      do i_aero_rep = 1, size(this%aero_rep)
        rep_spec_names = this%aero_rep(i_aero_rep)%val%unique_names()
        call assert(620697091, allocated(rep_spec_names))
        call assert(787495222, size(rep_spec_names).gt.0)
        forall (j_spec=1:size(rep_spec_names)) &
            state_names(i_spec+j_spec)%string = &
                rep_spec_names(j_spec)%string
        i_spec = i_spec + size(rep_spec_names)
        write(f_unit,*) "Aerosol rep ", &
                this%aero_rep(i_aero_rep)%val%rep_name, &
                " species: ", size(rep_spec_names)
        deallocate(rep_spec_names)
      end do
      do i_spec = 1, size(state_names)
        write(f_unit,*) i_spec-1, state_names(i_spec)%string
      end do

      write(f_unit,*) "*** Solver Data ***"
      write(f_unit,*) "Solver variables:"
      i_solver_spec = 0
      do i_spec = 1, size(state_names)
        if (this%var_type(i_spec).eq.CHEM_SPEC_VARIABLE) then
          write(f_unit,*) "solver spec", i_solver_spec, state_names(i_spec)%string
          i_solver_spec = i_solver_spec + 1
        end if
      end do
      write(f_unit,*) ""
      deallocate(state_names)
    end if

    flush(f_unit)

    if (associated(this%solver_data_gas)) &
            call this%solver_data_gas%print()
    if (associated(this%solver_data_gas_aero)) &
            call this%solver_data_gas_aero%print()

    flush(f_unit)

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the core
  elemental subroutine finalize(this)

    !> Phlex-core data
    type(phlex_core_t), intent(inout) :: this

    if (associated(this%mechanism)) &
            deallocate(this%mechanism)
    if (associated(this%chem_spec_data)) &
            deallocate(this%chem_spec_data)
    if (associated(this%sub_model)) &
            deallocate(this%sub_model)
    if (associated(this%aero_rep)) &
            deallocate(this%aero_rep)
    if (associated(this%aero_phase)) &
            deallocate(this%aero_phase)
    if (allocated(this%abs_tol)) &
            deallocate(this%abs_tol)
    if (allocated(this%var_type)) &
            deallocate(this%var_type)
    if (associated(this%solver_data_gas)) &
            deallocate(this%solver_data_gas)
    if (associated(this%solver_data_aero)) &
            deallocate(this%solver_data_aero)
    if (associated(this%solver_data_gas_aero)) &
            deallocate(this%solver_data_gas_aero)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a aerosol phase to the model data
  subroutine add_aero_phase(this, phase_name)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Aerosol phase name
    character(len=:), allocatable, intent(in) :: phase_name

    type(aero_phase_data_ptr), pointer :: new_aero_phase(:)

    allocate(new_aero_phase(size(this%aero_phase)+1))

    new_aero_phase(1:size(this%aero_phase)) = &
            this%aero_phase(1:size(this%aero_phase))

    new_aero_phase(size(new_aero_phase))%val => aero_phase_data_t(phase_name)

    deallocate(this%aero_phase)
    this%aero_phase => new_aero_phase

  end subroutine add_aero_phase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a aerosol representation to the model data
  subroutine add_aero_rep(this, rep_name)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Aerosol representation name
    character(len=:), allocatable, intent(in) :: rep_name

    type(aero_rep_data_ptr), pointer :: new_aero_rep(:)
    type(aero_rep_factory_t) :: aero_rep_factory

    !TODO: Improve this multiple reallocation
    allocate(new_aero_rep(size(this%aero_rep)+1))

    new_aero_rep(1:size(this%aero_rep)) = &
            this%aero_rep(1:size(this%aero_rep))
    new_aero_rep(size(new_aero_rep))%val => aero_rep_factory%create(rep_name)

    call this%aero_rep(:)%dereference()
    deallocate(this%aero_rep)
    this%aero_rep => new_aero_rep

  end subroutine add_aero_rep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a chemical mechanism to the model data
  subroutine add_mechanism(this, mech_name)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Mechanism name
    character(len=:), allocatable, intent(in) :: mech_name

    type(mechanism_data_ptr), pointer :: new_mechanism(:)

    allocate(new_mechanism(size(this%mechanism)+1))

    new_mechanism(1:size(this%mechanism)) = &
            this%mechanism(1:size(this%mechanism))

    new_mechanism(size(new_mechanism))%val => mechanism_data_t(mech_name)

    call this%mechanism(:)%dereference()
    deallocate(this%mechanism)
    this%mechanism => new_mechanism

  end subroutine add_mechanism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a sub-model to the model data
  subroutine add_sub_model(this, sub_model_name)

    !> Model data
    class(phlex_core_t), intent(inout) :: this
    !> Sub model name
    character(len=:), allocatable, intent(in) :: sub_model_name

    type(sub_model_data_ptr), pointer :: new_sub_model(:)
    type(sub_model_factory_t) :: sub_model_factory

    allocate(new_sub_model(size(this%sub_model)+1))

    new_sub_model(1:size(this%sub_model)) = &
            this%sub_model(1:size(this%sub_model))
    new_sub_model(size(new_sub_model))%val => &
            sub_model_factory%create(sub_model_name)

    deallocate(this%sub_model)
    this%sub_model => new_sub_model

  end subroutine add_sub_model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_phlex_core
