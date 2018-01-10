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

  implicit none
  private

  public :: mechanism_data_t

  !> Reallocation increment
  integer(kind=i_kind), parameter :: REALLOC_INC = 50

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
    !> Reactions
    type(rxn_data_ptr), pointer, public :: rxn_ptr(:) => null()
  contains
    !> Load reactions from an input file
    procedure :: load
    !> Initialize the mechanism
    procedure :: initialize
    !> Get the mechanism name
    procedure :: name => get_name
    !> Get the size of the species database
    procedure :: size => get_size
    !> Get constributions of mechanism reactions to the time derivative
    !! vector
    procedure :: get_func_contrib
    !> Get contributions of mechanism reactions to the Jaobian matrix
    procedure :: get_jac_contrib
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
    type(rxn_factory_t) :: rxn_factory

    ! Cycle through the set of reactions in the json file
    next => null()
    call json%get(j_obj, 'reactions(1)', child)
    do while (associated(child))

      ! Increase the size of the mechanism
      call this%ensure_size(1)
      this%num_rxn = this%num_rxn + 1
      
      ! Load the reaction into the mechanism
      this%rxn_ptr(this%num_rxn)%val => rxn_factory%load(json, child)
      
      ! Get the next reaction in the json file
      call json%get_next(child, next)
      child => next
    end do

#else
  subroutine load(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this

    call warn_msg(384838139, "No support for input files")
#endif

  end subroutine load
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the mechanism
  subroutine initialize(this, chem_spec_data)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    integer(kind=i_kind) :: i_rxn

    do i_rxn = 1, this%num_rxn
      call this%rxn_ptr(i_rxn)%val%initialize(chem_spec_data)
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

  !> Get the name of the mechanism
  function get_name(this) result(mech_name)

    !> Name of the mechanism
    character(len=:), allocatable :: mech_name
    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    mech_name = this%mech_name

  end function get_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get contributions of the mechanism reactions to the time derivative
  !! vector
  subroutine get_func_contrib(this, phlex_state, func)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Time derivative vector
    real(kind=dp), pointer, intent(inout) :: func(:)

    integer(kind=i_kind) :: i_rxn

    do i_rxn = 1, this%num_rxn
      if (this%rxn_ptr(i_rxn)%val%check_phase(phlex_state%rxn_phase)) then
        call this%rxn_ptr(i_rxn)%val%func_contrib(phlex_state, func)
      end if
    end do

  end subroutine get_func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get contributions of the mechanism reactions to the Jacobian matrix
  subroutine get_jac_contrib(this, phlex_state, jac_matrix)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Time derivative vector
    real(kind=dp), pointer, intent(inout) :: jac_matrix(:,:)

    integer(kind=i_kind) :: i_rxn

    do i_rxn = 1, this%num_rxn
      if (this%rxn_ptr(i_rxn)%val%check_phase(phlex_state%rxn_phase)) then
        call this%rxn_ptr(i_rxn)%val%jac_contrib(phlex_state, jac_matrix)
      end if
    end do

  end subroutine get_jac_contrib

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
  subroutine do_print(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    integer :: i_rxn

    write(*,*) "Mechanism: "//trim(this%name())
    do i_rxn = 1, this%num_rxn
      call this%rxn_ptr(i_rxn)%val%print()
    end do
    write(*,*) "End mechanism: "//trim(this%name())

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mechanism_data
