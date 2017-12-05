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
  contains
    !> Load model data
    procedure :: load => pmc_model_data_load
    !> Initialize the model
    procedure :: initialize => pmc_model_data_initialize
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

    !> Private functions
    !> Find mechanism by name
    procedure :: find_mechanism => pmc_model_data_find_mechanism
    !> Add a mechanism to the model data
    procedure :: add_mechanism => pmc_model_data_add_mechanism
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
    character(len=:), allocatable, optional :: input_file_path(:)

    allocate(new_obj)
    allocate(new_obj%mechanism(0))
    new_obj%chem_spec_data => chem_spec_data_t()

    if (present(input_file_path)) then
      call this%load(input_file_path)
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
    character(len=:), allocatable :: input_file_path(:)

    integer(kind=i_kind) :: i_file, i_mech
#ifdef PMC_USE_JSON
    type(json_core), pointer :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next

    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val
    character(len=:), allocatable :: str_val

    do i_file = 1, size(input_file_path)
      call j_file%initialize()
      call j_file%load_file(filename = input_file_path(i_file))      
      call j_file%get('pmc-data(1)', j_obj)
      do while (associated(j_obj))
        call json%get(j_obj, 'type', unicode_str_val, found)
        call assert_msg(689470331, found, "Missing type in json input file "//&
                input_file_path(i_file))
        str_val = unicode_str_val
        if (str_val.eq.'MECHANISM') then
          call json%get(j_obj, 'name', unicode_str_val, found)
          call assert_msg(822680732, "Missing mechanism name in file "//&
                  input_file_path(i_file))
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

  end subroutine pmc_model_data_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the model data
  subroutine pmc_model_data_initialize(this)

    !> Model data
    class(model_data_t), intent(inout) :: this

    forall this%mechanism
      call this%mechanism%initialize(this%chem_spec_data)
    end forall

  end subroutine pmc_model_data_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find a mechanism by name in the model data
  logical function pmc_model_data_find_mechanism(this, mech_name, mech_id) &
                  result(found)

    !> Model data
    type(model_data_t), intent(in) :: this
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
    type(model_data_t), intent(inout) :: this
    !> Mechanism name
    character(len=:), allocatable :: mech_name

    type(mechanism_data_t), pointer :: new_mechanism(:)

    allocate(new_mechanism(size(this%mechanism)+1))

    new_mechanism(1:size(this%mechanism)) = &
            this%mechanism(1:size(this%mechanism))

    deallocate(this%mechanism)
    this%mechanism => new_mechanism

  end subroutine pmc_model_data_add_mechanism

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end model pmc_model_data
