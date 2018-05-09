! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_mechanism_data module.

!> \page phlex_mechanism Phlexible Module for Chemistry: Chemical Mechanism
!!
!! A mechanism in the \ref phlex_chem "phlex-chem" module is a set of
!! \ref phlex_rxn "reactions" that occur in the gas-phase or within one of
!! several \ref phlex_aero_phase "aerosol phases" or across an interface
!! between two phases (gas or aerosol). One or several mechanisms may be
!! included in a \ref phlex_chem "phlex-chem" model run. 
!!
!! Every mechanism in a \ref phlex_chem "phlex-chem" run will have access to
!! the same set of \ref phlex_species "chemical species" and \ref
!! phlex_aero_phase "aerosol phases", so phase and species names must be 
!! consistent across all concurrently loaded mechanisms. The division of \ref
!! phlex_rxn "reactions" into distinct mechanisms permits a host model to
!! specificy which mechanisms should be solved during a call to 
!! \c pmc_phlex_core::phlex_core_t::solve().
!!
!! The input format for mechanism data can be found \ref
!! input_format_mechanism "here".

!> The mechanism_data_t structure and associated subroutines.
module pmc_mechanism_data

  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t
#ifdef PMC_USE_MPI
  use mpi
#endif
#ifdef PMC_USE_JSON
  use json_module
#endif
  ! Reaction modules
  use pmc_rxn_data
  use pmc_rxn_factory
  use pmc_chem_spec_data
  use pmc_phlex_state
  use pmc_aero_rep_data

  implicit none
  private

  public :: mechanism_data_t

  !> Reallocation increment
  integer(kind=i_kind), parameter :: REALLOC_INC = 50
  !> Fixed module file unit
  integer(kind=i_kind), parameter :: MECH_FILE_UNIT = 16

  !> A chemical mechanism
  !!
  !! Instances of mechanism_data_t represent complete \ref phlex_mechanism 
  !! chemical mechanisms (e.g., CACM/MPMPO, CB-5, EQSAM). Multiple mechanisms
  !! may be used during one model run and will be solved simultaneously.
  type :: mechanism_data_t
    private
    !> Number of reactions
    integer(kind=i_kind) :: num_rxn = 0
    !> Mechanism name
    character(len=:), allocatable :: mech_name
    !> Path and prefix for fixed module output
    character(len=:), allocatable :: fixed_file_prefix
    !> Reactions
    type(rxn_data_ptr), pointer :: rxn_ptr(:) => null()
  contains
    !> Load reactions from an input file
    procedure :: load
    !> Initialize the mechanism
    procedure :: initialize
    !> Get the mechanism name
    procedure :: name => get_name
    !> Get the size of the reaction database
    procedure :: size => get_size
    !> Get a reaction by its index
    procedure :: get_rxn
    !> Build a fixed mechanism module
    procedure :: build_fixed_module
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Print the mechanism data
    procedure :: print => do_print

    ! Private functions
    !> Ensure there is enough room in the reaction dataset to add a
    !! specified number of reactions
    procedure, private :: ensure_size
  end type mechanism_data_t

  ! Constructor for mechanism_data_t
  interface mechanism_data_t
    procedure :: constructor
  end interface mechanism_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for mechanism_data_t
  function constructor(mech_name, init_size) result(new_obj)

    !> Chemical mechanism
    type(mechanism_data_t), pointer :: new_obj
    !> Name of the mechanism
    character(len=:), allocatable, intent(in) :: mech_name
    !> Number of reactions to allocate space for initially
    integer(i_kind), intent(in), optional :: init_size

    integer(i_kind) :: alloc_size = REALLOC_INC

    allocate(new_obj)
    if (present(init_size)) alloc_size = init_size
    new_obj%mech_name = mech_name
    allocate(new_obj%rxn_ptr(alloc_size))

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Ensure there is enough room in the reaction dataset to add a specified
  !! number of reactions
  subroutine ensure_size(this, num_rxn)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> Number of new reactions to ensure space for
    integer(i_kind), intent(in) :: num_rxn

    integer(kind=i_kind) :: new_size
    integer(kind=i_kind), allocatable :: new_rxn_type(:)
    type(rxn_data_ptr), pointer :: new_rxn_ptr(:)

    if (size(this%rxn_ptr) .ge. this%num_rxn + num_rxn) return
    new_size = this%num_rxn + num_rxn + REALLOC_INC
    allocate(new_rxn_ptr(new_size))
    new_rxn_ptr(1:this%num_rxn) = this%rxn_ptr(1:this%num_rxn)
    deallocate(this%rxn_ptr)
    this%rxn_ptr => new_rxn_ptr

  end subroutine ensure_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> \page input_format_mechanism Input JSON Object Format: Mechanism
  !!
  !! A \c json object containing information about a \ref phlex_mechanism
  !! "chemical mechanism" of the form:
  !! \code{.json}
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my mechanism",
  !!     "type" : "MECHANISM",
  !!     "reactions" : [
  !!       ...
  !!     ]
  !!   }
  !! ]}
  !! \endcode
  !! A \ref phlex_mechanism "mechanism" object must have a unique \b name, 
  !! a \b type of \b MECHANISM and an array of \ref input_format_rxn
  !! "reaction objects" labelled \b reactions. Mechanism data may be split
  !! into multiple mechanism objects - they will be combined based on the
  !! mechanism name.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a chemical mechanism from an input file
#ifdef PMC_USE_JSON
  subroutine load(this, json, j_obj)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    type(json_value), pointer :: child, next
    character(kind=json_ck, len=:), allocatable :: unicode_str_val
    type(rxn_factory_t) :: rxn_factory
    logical :: found

    ! Cycle through the set of reactions in the json file
    next => null()
   
    ! Get the reaction set
    call json%get(j_obj, 'reactions(1)', child, found)
    do while (associated(child) .and. found)

      ! Increase the size of the mechanism
      call this%ensure_size(1)
      this%num_rxn = this%num_rxn + 1
      
      ! Load the reaction into the mechanism
      this%rxn_ptr(this%num_rxn)%val => rxn_factory%load(json, child)
      
      ! Get the next reaction in the json file
      call json%get_next(child, next)
      child => next
    end do

    ! Determine whether and where to build fixed module code
    call json%get(j_obj, 'build fixed module', unicode_str_val, found)
    if (found) then
      call assert_msg(410823202, .not.allocated(this%fixed_file_prefix), &
              "Received multiple file prefixes for fixed mechanism module.")
      this%fixed_file_prefix = trim(unicode_str_val)
    end if

#else
  subroutine load(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this

    call warn_msg(384838139, "No support for input files")
#endif

  end subroutine load
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the mechanism
  subroutine initialize(this, chem_spec_data, aero_rep_data)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representation data
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep_data(:)

    integer(kind=i_kind) :: i_rxn

    do i_rxn = 1, this%num_rxn
      call assert(340397127, associated(this%rxn_ptr(i_rxn)%val))
      call this%rxn_ptr(i_rxn)%val%initialize(chem_spec_data, aero_rep_data)
    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the current size of the chemical mechanism
  integer(kind=i_kind) function get_size(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    get_size = this%num_rxn

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a reaction by its index
  function get_rxn(this, rxn_id) result (rxn_ptr)

    !> Pointer to the reaction
    class(rxn_data_t), pointer :: rxn_ptr
    !> Mechanism data
    class(mechanism_data_t), intent(in) :: this
    !> Reaction index
    integer(kind=i_kind), intent(in) :: rxn_id

    call assert_msg(129484547, rxn_id.gt.0 .and. rxn_id .le. this%num_rxn, &
            "Invalid reaction id: "//trim(to_string(rxn_id))//&
            "exptected a value between 1 and "//trim(to_string(this%num_rxn)))

    rxn_ptr => this%rxn_ptr(rxn_id)%val

  end function get_rxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of the mechanism
  function get_name(this) result(mech_name)

    !> Name of the mechanism
    character(len=:), allocatable :: mech_name
    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    mech_name = this%mech_name

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build a fixed mechanism module and save it with the specified file prefix
  subroutine build_fixed_module(this, state_array_size)

    !> Mechanism data
    class(mechanism_data_t), intent(in) :: this
    !> Size of the state array
    integer(kind=i_kind), intent(in) :: state_array_size

    ! Output file name
    character(len=:), allocatable :: file_name
    ! File unit
    integer(kind=i_kind) :: U = MECH_FILE_UNIT
    ! Counters
    integer(kind=i_kind) :: i_rxn
    ! Temporary string for building file lines
    character(len=:), allocatable :: str_temp

    ! Determine whether to create fixed module
    if (.not.allocated(this%fixed_file_prefix)) return

    ! Open the output file
    file_name = trim(this%fixed_file_prefix)//"_mechanism.F90"
    open(unit=U, file=file_name, status="replace", action="write")

    write(*,*) "Building output file: "//file_name

    write(U,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(U,*) "!! ",trim(this%mech_name)," fixed mechanism module"
    write(U,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write(U,*) ""
    write(U,*) "module mech_"//trim(this%mech_name)
    write(U,*) ""
    write(U,*) "  use pmc_util,                                   only : i_kind, dp"
    write(U,*) "  implicit none"
    write(U,*) ""
    write(U,*) "  integer(kind=i_kind), parameter :: NUM_RXN = ", this%num_rxn
    write(U,*) "  integer(kind=i_kind), parameter :: NUM_VAR = ", state_array_size
    write(U,*) "  real(kind=dp) :: RATE_CONST(NUM_RXN)"
    write(U,*) "  real(kind=dp) :: DERIV(NUM_VAR)"
    write(U,*) ""
    write(U,*) "contains"
    write(U,*) ""
    write(U,*) "  subroutine calc_rate_const(temp_K, press_atm)"
    write(U,*) ""
    write(U,*) "    real(kind=dp), intent(in) :: temp_K"
    write(U,*) "    real(kind=dp), intent(in) :: press_atm"
    write(U,*) ""
    write(U,*) "    real(kind=dp) :: press_Pa = press_atm * 101325"
    write(U,*) ""
    do i_rxn = 1, this%num_rxn
      str_temp = this%rxn_ptr(i_rxn)%val%build_rate_const_expr(i_rxn)
    write(U,*) "    RATE_CONST("//trim(to_string(i_rxn))//") = "//trim(str_temp)
    end do
    write(U,*) ""
    write(U,*) "  end subroutine calc_rate_const"

    write(U,*) ""
    write(U,*) "end module mech_"//trim(this%mech_name)

    close(U)

  end subroutine build_fixed_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the mechanism
  integer(kind=i_kind) function pack_size(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
   
    type(rxn_factory_t) :: rxn_factory 
    integer(kind=i_kind) :: i_rxn

    pack_size =  pmc_mpi_pack_size_integer(this%num_rxn) + &
                 pmc_mpi_pack_size_string(this%mech_name)
    do i_rxn = 1, this%num_rxn
      associate (rxn => this%rxn_ptr(i_rxn)%val)
      pack_size = pack_size + rxn_factory%pack_size(rxn)
      end associate
    end do

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    type(rxn_factory_t) :: rxn_factory 
    integer :: i_rxn, prev_position

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, this%num_rxn)
    call pmc_mpi_pack_string(buffer, pos, this%mech_name)
    do i_rxn = 1, this%num_rxn
      associate (rxn => this%rxn_ptr(i_rxn)%val)
      call rxn_factory%bin_pack(rxn, buffer, pos)
      end associate
    end do
    call assert(669506045, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    type(rxn_factory_t) :: rxn_factory 
    integer :: i_rxn, prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, this%num_rxn)
    call pmc_mpi_unpack_string(buffer, pos, this%mech_name)
    call this%ensure_size(this%num_rxn)
    do i_rxn = 1, this%num_rxn
      this%rxn_ptr(i_rxn)%val => rxn_factory%bin_unpack(buffer, pos)
    end do
    call assert(360900030, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the mechanism data
  subroutine do_print(this, file_unit)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), optional :: file_unit

    integer :: i_rxn
    integer(kind=i_kind) :: f_unit = 6

    if (present(file_unit)) f_unit = file_unit

    write(f_unit,*) "Mechanism: "//trim(this%name())
    do i_rxn = 1, this%num_rxn
      call this%rxn_ptr(i_rxn)%val%print(f_unit)
    end do
    write(f_unit,*) "End mechanism: "//trim(this%name())

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mechanism_data
