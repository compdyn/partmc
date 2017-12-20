! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_mechanism_data module.

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
  use pmc_rxn_arrhenius
  use pmc_chem_spec_data
  use pmc_model_state

  implicit none
  private

  public :: mechanism_data_t

  !> Reallocation increment
  integer(kind=i_kind), parameter :: REALLOC_INC = 50

  !> Arrhenius equation
  integer(kind=i_kind), parameter :: RXN_ARRHENIUS = 1

  !> A chemical mechanism
  !!
  !! Instances of mechanism_data_t represent complete chemical
  !! mechanisms (e.g., CACM/MPMPO, CB-5, EQSAM). Multiple mechanisms
  !! may be used during one model run and will be integrated
  !! together.
  type :: mechanism_data_t
    private
    !> Number of reactions
    integer(kind=i_kind) :: num_rxn = 0
    !> Mechanism name
    character(len=:), allocatable :: mech_name
    !> Reaction type
    integer(kind=i_kind), allocatable :: rxn_type(:)
    !> Reactions
    type(rxn_data_ptr), pointer :: rxn_ptr(:) => null()
  contains
    !> Load reactions from an input file
    procedure :: load => pmc_mechanism_data_load
    !> Initialize the mechanism
    procedure :: initialize => pmc_mechanism_data_initialize
    !> Get the mechanism name
    procedure :: name => pmc_mechanism_data_name
    !> Get the size of the species database
    procedure :: size => pmc_mechanism_data_size
    !> Get constributions of mechanism reactions to the time derivative
    !! vector
    procedure :: get_func_contrib => pmc_mechanism_data_get_func_contrib
    !> Get contributions of mechanism reactions to the Jaobian matrix
    procedure :: get_jac_contrib => pmc_mechanism_data_get_jac_contrib
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size => pmc_mechanism_data_pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack => pmc_mechanism_data_bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack => pmc_mechanism_data_bin_unpack
    !> Print the mechanism data
    procedure :: print => pmc_mechanism_data_print

    !> Private functions
    !> Ensure there is enough room in the reaction dataset to add a
    !! specified number of reactions
    procedure, private :: ensure_size => pmc_mechanism_data_ensure_size
  end type mechanism_data_t

  !> Constructor for mechanism_data_t
  interface mechanism_data_t
    procedure :: pmc_mechanism_data_constructor
  end interface mechanism_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for mechanism_data_t
  function pmc_mechanism_data_constructor(mech_name, init_size) result(new_obj)

    !> Chemical mechanism
    type(mechanism_data_t), target :: new_obj
    !> Name of the mechanism
    character(len=:), allocatable :: mech_name
    !> Number of reactions to allocate space for initially
    integer(i_kind), intent(in), optional :: init_size

    integer(i_kind) :: alloc_size = REALLOC_INC

    if (present(init_size)) alloc_size = init_size
    new_obj%mech_name = mech_name
    allocate(new_obj%rxn_type(alloc_size))
    allocate(new_obj%rxn_ptr(alloc_size))

  end function pmc_mechanism_data_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Ensure there is enough room in the reaction dataset to add a specified
  !! number of reactions
  subroutine pmc_mechanism_data_ensure_size(this, num_rxn)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> Number of new reactions to ensure space for
    integer(i_kind) :: num_rxn

    integer(kind=i_kind) :: new_size
    integer(kind=i_kind), allocatable :: new_rxn_type(:)
    type(rxn_data_ptr), pointer :: new_rxn_ptr(:)

    if (size(this%rxn_ptr) .ge. this%num_rxn + num_rxn) return
    new_size = this%num_rxn + num_rxn + REALLOC_INC
    allocate(new_rxn_type(new_size))
    allocate(new_rxn_ptr(new_size))
    new_rxn_type(1:this%num_rxn) = this%rxn_type(1:this%num_rxn)
    new_rxn_ptr(1:this%num_rxn) = this%rxn_ptr(1:this%num_rxn)
    deallocate(this%rxn_type)
    deallocate(this%rxn_ptr)
    this%rxn_type = new_rxn_type
    this%rxn_ptr => new_rxn_ptr

  end subroutine pmc_mechanism_data_ensure_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a chemical mechanism from an input file
#ifdef PMC_USE_JSON
  !! j_obj is expected to be a JSON object containing data related to a
  !! chemical mechanism. It should be of the form:
  !! 
  !! { "pmc-data" : [
  !!   {
  !!     "name" : "my mechanism",
  !!     "type" : "MECHANISM",
  !!     "reactions" : [
  !!       ...
  !!     ]
  !!   }
  !! ]}
  !!
  !! For the structure of the reaction objects, see the pmc_rxn_data module.
  !! Mechanism data may be split into multiple mechanism objects - they will
  !! be combined based on the mechanism name. All mechanisms objects must have
  !! a name, a type = "MECHANISM" and a reaction array containing at least
  !! one reaction.
  subroutine pmc_mechanism_data_load(this, json, j_obj)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    logical :: found
    type(json_value), pointer :: child, next
    character(kind=json_ck, len=:), allocatable :: key, unicode_str_val
    character(len=:), allocatable :: str_val

    next => null()
    call json%get(j_obj, 'reactions(1)', child)
    do while (associated(child))
      call this%ensure_size(1)
      this%num_rxn = this%num_rxn + 1
      call json%get(child, "type", unicode_str_val, found)
      call assert_msg(303914867, found, "Missing reaction type") 
      str_val = unicode_str_val
      if (str_val .eq. "ARRHENIUS") then
        this%rxn_type = RXN_ARRHENIUS
        this%rxn_ptr(this%num_rxn)%val => rxn_arrhenius_t()
      else 
        call die_msg(359134071, "Invalid reaction type: "//trim(str_val))
      end if
      call this%rxn_ptr(this%num_rxn)%val%load(json, child)
      call json%get_next(child, next)
      child => next
    end do

#else
  subroutine pmc_mechanism_data_load(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this

    call warn_msg(384838139, "No support for input files")
#endif

  end subroutine pmc_mechanism_data_load
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the mechanism
  subroutine pmc_mechanism_data_initialize(this, chem_spec_data)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    integer(kind=i_kind) :: i_rxn

    do i_rxn = 1, this%num_rxn
      call this%rxn_ptr(i_rxn)%val%initialize(chem_spec_data)
    end do

  end subroutine pmc_mechanism_data_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the current size of the chemical mechanism
  integer(kind=i_kind) function pmc_mechanism_data_size(this) &
                  result(mech_size)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    mech_size = this%num_rxn

  end function pmc_mechanism_data_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of the mechanism
  function pmc_mechanism_data_name(this) result(mech_name)

    !> Name of the mechanism
    character(len=:), allocatable :: mech_name
    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    mech_name = this%mech_name

  end function pmc_mechanism_data_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get contributions of the mechanism reactions to the time derivative
  !! vector
  subroutine pmc_mechanism_data_get_func_contrib(this, model_state, func)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Current model state
    type(model_state_t), intent(in) :: model_state
    !> Time derivative vector
    real(kind=dp), pointer, intent(inout) :: func(:)

    integer(kind=i_kind) :: i_rxn

    do i_rxn = 1, this%num_rxn
      if (this%rxn_ptr(i_rxn)%val%check_phase(model_state%rxn_phase)) then
        call this%rxn_ptr(i_rxn)%val%func_contrib(model_state, func)
      end if
    end do

  end subroutine pmc_mechanism_data_get_func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get contributions of the mechanism reactions to the Jacobian matrix
  subroutine pmc_mechanism_data_get_jac_contrib(this, model_state, jac_matrix)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Current model state
    type(model_state_t), intent(in) :: model_state
    !> Time derivative vector
    real(kind=dp), pointer, intent(inout) :: jac_matrix(:,:)

    integer(kind=i_kind) :: i_rxn

    do i_rxn = 1, this%num_rxn
      if (this%rxn_ptr(i_rxn)%val%check_phase(model_state%rxn_phase)) then
        call this%rxn_ptr(i_rxn)%val%jac_contrib(model_state, jac_matrix)
      end if
    end do

  end subroutine pmc_mechanism_data_get_jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the mechanism
  integer(kind=i_kind) function pmc_mechanism_data_pack_size(this) &
                  result (pack_size)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    
    integer(kind=i_kind) :: i_rxn

    pack_size =  pmc_mpi_pack_size_integer(this%num_rxn) + &
                 pmc_mpi_pack_size_string(this%mech_name) + &
                 pmc_mpi_pack_size_integer_array(this%rxn_type)
    do i_rxn = 1, this%num_rxn
      pack_size = pack_size + this%rxn_ptr(i_rxn)%val%pack_size()
    end do

  end function pmc_mechanism_data_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine pmc_mechanism_data_bin_pack(this, buffer, pos)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: i_rxn, prev_position

    prev_position = pos
    call pmc_mpi_pack_integer(buffer, pos, this%num_rxn)
    call pmc_mpi_pack_string(buffer, pos, this%mech_name)
    call pmc_mpi_pack_integer_array(buffer, pos, this%rxn_type)
    do i_rxn = 1, this%num_rxn
      call this%rxn_ptr(i_rxn)%val%bin_pack(buffer, pos)
    end do
    call assert(669506045, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_mechanism_data_bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine pmc_mechanism_data_bin_unpack(this, buffer, pos)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: i_rxn, prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, this%num_rxn)
    call pmc_mpi_unpack_string(buffer, pos, this%mech_name)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%rxn_type)
    do i_rxn = 1, this%num_rxn
      select (this%rxn_type(i_rxn))
        case (RXN_ARRHENIUS)
          this%rxn_ptr(i_rxn)%val => rxn_arrhenius_t()
        case default
          call die_msg(655357132, "Trying to unpack unknown reaction type: "&
                  //to_string(this%rxn_type(i_rxn)))
      end select
      call this%rxn_ptr(i_rxn)%val%bin_unpack(buffer, pos)
    end do
    call assert(360900030, &
         pos - prev_position <= this%pack_size())
#endif

  end subroutine pmc_mechanism_data_bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the mechanism data
  subroutine pmc_mechanism_data_print(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    integer :: i_rxn

    write(*,*) "Mechanism: "//trim(this%name())
    do i_rxn = 1, this%num_rxn
      call this%rxn_ptr(i_rxn)%val%print()
    end do
    write(*,*) "End mechanism: "//trim(this%name())

  end subroutine pmc_mechanism_data_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mechanism_data
