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
  use rxn_arrhenius


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
    character(len=:), allocatable :: name
    !> Reaction type
    integer(kind=i_kind), pointer :: rxn_type(:) => null()
    !> Reactions
    type(rxn_data_t), pointer :: reaction(:) => null()
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
    new_obj%name = mech_name
    allocate(new_obj%rxn_type(alloc_size))
    allocate(new_obj%reactions(alloc_size))

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
    integer(kind=i_kind) :: new_rxn_type(:)
    type(rxn_data_t), pointer :: new_reaction(:)

    if (size(this%reaction) .ge. this%num_rxn + num_rxn) return
    new_size = this%num_rxn + num_rxn + REALLOC_INC
    allocate(new_rxn_type(new_size))
    allocate(new_reaction(new_size))
    new_rxn_type(1:this%num_spec) = this%rxn_type(1:this%num_spec)
    call this%reaction(1:this%num_spec)%move(new_reaction(1:this%num_spec))
    deallocate(this%rxn_type)
    deallocate(this%reaction)
    this%rxn_type => new_rxn_type
    this%reaction => new_reaction

  end subroutine mechanism_data_ensure_size

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
        this%reaction(this%num_rxn) = rxn_arrhenius_t()
      else 
        die_msg(359134071, "Invalid reaction type: "//trim(str_val))
      end if
      this%reaction(this%num_rxn)%load(json, child)
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
  subroutine pmc_mechanism_data_initialize(this)

    !> Chemical mechanism
    class(mechanism_data_t), intent(inout) :: this

    forall this%reaction
      call this%reaction%initialize()
    end forall

  end subroutine pmc_mechanism_data_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the current size of the chemical mechanism
  integer(kind=i_kind) function pmc_mechanism_data_size(this) &
                  result(mech_size)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    mech_size = this%num_rxn

  end function pmc_mechanism_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the name of the mechanism
  function pmc_mechanism_name(this) result(mech_name)

    !> Name of the mechanism
    character(len=:), allocatable :: mech_name
    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this

    mech_name = this%name

  end function pmc_mechanism_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get contributions of the mechanism reactions to the time derivative
  !! vector
  subroutine pmc_mechanism_data_get_func_contrib(this, model_state, func)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Current model state
    type(model_state_t), intent(in) :: model_state
    !> Time derivative vector
    real(kind=dp), allocatable, intent(inout) :: func(:)

    forall this%reactions
      call this%reaction%func_contrib(model_state, func)
    end forall

  end subroutine pmc_mechanism_data_get_func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get contributions of the mechanism reactions to the Jacobian matrix
  subroutine pmc_mechanism_data_get_jac_contrib(this, model_state, jac_matrix)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    !> Current model state
    type(model_state_t), intent(in) :: model_state
    !> Time derivative vector
    real(kind=dp), allocatable, intent(inout) :: jac_matrix(:,:)

    forall this%reactions
      call this%reaction%jac_contrib(model_state, jac_matrix)
    end forall

  end subroutine pmc_mechanism_data_get_jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the mechanism
  integer(kind=i_kind) function pmc_mechanism_data_pack_size(this) &
                  result (pack_size)

    !> Chemical mechanism
    class(mechanism_data_t), intent(in) :: this
    
    integer(kind=i_kind) :: i_rxn

    pack_size =  pmc_mpi_pack_size_integer(this%num_rxn) + &
                 pmc_pack_size_string(this%name) + &
                 pmc_mpi_pack_size_integer_array(this%rxn_type)
    do i_rxn = 1, this%num_rxn
      pack_size = pack_size + this%reaction(i_rxn)%pack_size()
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
    call pmc_mpi_pack_string(buffer, pos, this%name)
    call pmc_mpi_pack_integer_array(buffer, pos, this%rxn_type)
    do i_rxn = 1, this%num_rxn
      this%reactions(i_rxn)%bin_pack(buffer, pos)
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
    call pmc_mpi_unpack_string(buffer, pos, this%name)
    call pmc_mpi_unpack_integer_array(buffer, pos, this%rxn_type)
    do i_rxn = 1, this%num_rxn
      select (this%rxn_type(i_rxn))
        case (RXN_ARRHENIUS)
          this%reaction(i_rxn) => rxn_arrhenius_t()
        case default
          call die_msg(655357132, "Trying to unpack unknown reaction type: "&
                  //to_string(this%rxn_type(i_rxn)))
      end select
      this%reactions(i_rxn)%bin_unpack(buffer, pos)
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

    write(*,*) "Mechanism: "//trim(this%name)
    do i_rxn = 1, this%num_rxn
      this%reaction(i_rxn)%print()
    end do
    write(*,*) "End mechanism: "//trim(this%name)

  end subroutine pmc_mechanism_data_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mechanism_data
