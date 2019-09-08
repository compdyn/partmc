! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_phlex_state module.

!> The phlex_state_t structure and associated subroutines.
module pmc_phlex_state

! Define array size for contain temperature and pressure
#define PHLEX_STATE_NUM_ENV_PARAM 2

#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_env_state
  use pmc_mpi
  use pmc_rand
  use pmc_util,                       only : die_msg, string_t

  implicit none
  private

  public :: phlex_state_t, phlex_state_ptr

  !> Model state
  !!
  !! Temporal state of the model
  type phlex_state_t
    !> Environmental state array. This array will include one entry
    !! for every environmental variable requried to solve the
    !! chemical mechanism(s)
    real(kind=dp), allocatable :: env_var(:)
    !> State variable array. This array includes one entry for each
    !! variable whose state will be solved for during the mechanism
    !! integration.
    real(kind=dp), allocatable :: state_var(:)
    !> Environmental conditions
    type(env_state_t), pointer :: env_state => null()
    !> Flag indicating whether the env_state object is owned by the
    !! state object
    logical, private :: owns_env_state = .false.
  contains
    !> Update the environmental state array
    procedure :: update_env_state
    !> Determine the size of a binary required to pack a given variable
    procedure :: pack_size
    !> Pack the given value to the buffer, advancing position
    procedure :: bin_pack
    !> Unpack the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Finalize the state
    final :: finalize
  end type phlex_state_t

  ! Constructor for phlex_state_t
  interface phlex_state_t
    procedure :: constructor
  end interface phlex_state_t

  !> Pointer type for building arrays
  type phlex_state_ptr
    type(phlex_state_t), pointer :: val => null()
  contains
    !> Dereference the pointer
    procedure :: dereference
    !> Finalize the pointer
    final :: ptr_finalize
  end type phlex_state_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for phlex_state_t
  function constructor(env_state, n_cells) result (new_obj)

    !> New model state
    type(phlex_state_t), pointer :: new_obj
    !> Num cells to compute simulatenously
    integer(kind=i_kind), optional :: n_cells
    !> Environmental state
    type(env_state_t), target, intent(in), optional :: env_state

    ! Allocate space for the new object
    allocate(new_obj)

    ! Set the environmental state (if present)
    if (present(env_state)) then
      new_obj%env_state => env_state
    else
      allocate(new_obj%env_state)
      new_obj%owns_env_state = .true.
    end if

    ! Set up the environmental state array
    allocate(new_obj%env_var(PHLEX_STATE_NUM_ENV_PARAM*n_cells))

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the environmental state array
  subroutine update_env_state(this, grid_cell)

    !> Model state
    class(phlex_state_t), intent(inout) :: this
    !> Grid cell to update
    integer, optional :: grid_cell

    integer :: grid_offset = 0

    if (present(grid_cell)) &
      grid_offset = (grid_cell-1)*PHLEX_STATE_NUM_ENV_PARAM

    call assert_msg(618562571, grid_offset >= 0 .and. &
                               grid_offset <= size(this%env_var)-2, &
                    "Invalid env state offset: "//trim(to_string(grid_offset)))

    this%env_var(grid_offset+1) = this%env_state%temp     ! Temperature (K)
    this%env_var(grid_offset+2) = this%env_state%pressure ! Pressure (Pa)

  end subroutine update_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack a given variable
  integer(kind=i_kind) function pack_size(this)

    !> Chemical species states
    class(phlex_state_t), intent(in) :: this

    pack_size = &
            this%env_state%pack_size()

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

    !> Chemical species states
    class(phlex_state_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
     integer :: prev_position

     prev_position = pos
     call this%env_state%bin_pack(buffer, pos)
     call assert(204737650, &
             pos - prev_position <= this%pack_size())
#endif

   end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

    !> Chemical species states
    class(phlex_state_t), intent(inout) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
     integer :: prev_position

     prev_position = pos
     call this%env_state%bin_pack(buffer, pos)
     call assert(820809161, &
             pos - prev_position <= this%pack_size())
#endif

   end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the state
  elemental subroutine finalize(this)

    !> Phlex model state
    type(phlex_state_t), intent(inout) :: this

    if (allocated(this%env_var)) deallocate(this%env_var)
    if (allocated(this%state_var)) deallocate(this%state_var)
    if (associated(this%env_state) .and. this%owns_env_state) &
            deallocate(this%env_state)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deference a pointer to a phlex state
  elemental subroutine dereference(this)

    !> Pointer to the phlex state
    class(phlex_state_ptr), intent(inout) :: this

    this%val => null()

  end subroutine dereference

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a pointer to a phlex state
  elemental subroutine ptr_finalize(this)

    !> Pointer to the phlex state
    type(phlex_state_ptr), intent(inout) :: this

    if (associated(this%val)) deallocate(this%val)

  end subroutine ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef PHLEX_STATE_NUM_ENV_PARAM
end module pmc_phlex_state
