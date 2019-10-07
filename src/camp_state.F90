! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_camp_state module.

!> The camp_state_t structure and associated subroutines.
module pmc_camp_state

! Define array size for contain temperature and pressure
#define CAMP_STATE_NUM_ENV_PARAM 2

#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_env_state
  use pmc_mpi
  use pmc_rand
  use pmc_util,                       only : die_msg, string_t

  implicit none
  private

  public :: camp_state_t, camp_state_ptr

  !> Model state
  !!
  !! Temporal state of the model
  type camp_state_t
    !> Environmental state array. This array will include one entry
    !! for every environmental variable requried to solve the
    !! chemical mechanism(s)
    real(kind=dp), allocatable :: env_var(:)
    !> State variable array. This array includes one entry for each
    !! variable whose state will be solved for during the mechanism
    !! integration.
    real(kind=dp), allocatable :: state_var(:)
    !> Environmental conditions
    type(env_state_ptr), pointer :: env_states(:)
    !> Flag indicating whether the env_state object is owned by the
    !! state object
    logical, private :: owns_env_states = .false.
  contains
    !> Update the environmental state array
    procedure :: update_env_state
    !> Finalize the state
    final :: finalize
  end type camp_state_t

  ! Constructor for camp_state_t
  interface camp_state_t
    procedure :: constructor_one_cell, constructor_multi_cell
  end interface camp_state_t

  !> Pointer type for building arrays
  type camp_state_ptr
    type(camp_state_t), pointer :: val => null()
  contains
    !> Dereference the pointer
    procedure :: dereference
    !> Finalize the pointer
    final :: ptr_finalize
  end type camp_state_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for camp_state_t
  function constructor_one_cell(env_state) result (new_obj)

    !> New model state
    type(camp_state_t), pointer :: new_obj
    !> Environmental state
    type(env_state_t), target, intent(in), optional :: env_state

    ! Allocate space for the new object
    allocate(new_obj)
    allocate(new_obj%env_states(1))

    ! Set the environmental state (if present)
    if (present(env_state)) then
      new_obj%env_states(1)%val => env_state
    else
      allocate(new_obj%env_states(1)%val)
      new_obj%owns_env_states = .true.
    end if

    ! Set up the environmental state array
    allocate(new_obj%env_var(CAMP_STATE_NUM_ENV_PARAM))

  end function constructor_one_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for camp_state_t
  function constructor_multi_cell(num_cells, env_states) result (new_obj)

    !> New model state
    type(camp_state_t), pointer :: new_obj
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: num_cells
    !> Environmental state
    type(env_state_ptr), target, intent(in), optional :: env_states(:)

    integer(kind=i_kind) :: i_cell

    ! Allocate space for the new object
    allocate(new_obj)

    ! Set the environmental state (if present)
    if (present(env_states)) then
      new_obj%env_states => env_states
    else
      allocate(new_obj%env_states(num_cells))
      do i_cell = 1, num_cells
        allocate(new_obj%env_states(i_cell)%val)
      end do
      new_obj%owns_env_states = .true.
    end if

    ! Set up the environmental state array
    allocate(new_obj%env_var(CAMP_STATE_NUM_ENV_PARAM*num_cells))

  end function constructor_multi_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the environmental state array
  subroutine update_env_state(this)

    !> Model state
    class(camp_state_t), intent(inout) :: this

    integer :: i_cell, grid_offset

    do i_cell = 1, size(this%env_states)
      grid_offset = (i_cell-1)*CAMP_STATE_NUM_ENV_PARAM
      this%env_var(grid_offset+1) = this%env_states(i_cell)%val%temp          ! Temperature (K)
      this%env_var(grid_offset+2) = this%env_states(i_cell)%val%pressure      ! Pressure (Pa)
    end do

  end subroutine update_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the state
  elemental subroutine finalize(this)

    !> CAMP model state
    type(camp_state_t), intent(inout) :: this

    integer(kind=i_kind) :: i_cell

    if (allocated(this%env_var)) deallocate(this%env_var)
    if (allocated(this%state_var)) deallocate(this%state_var)
    if (associated(this%env_states) .and. this%owns_env_states) then
      do i_cell = 1, size(this%env_states)
        deallocate(this%env_states(i_cell)%val)
      end do
      deallocate(this%env_states)
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Deference a pointer to a camp state
  elemental subroutine dereference(this)

    !> Pointer to the camp state
    class(camp_state_ptr), intent(inout) :: this

    this%val => null()

  end subroutine dereference

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a pointer to a camp state
  elemental subroutine ptr_finalize(this)

    !> Pointer to the camp state
    type(camp_state_ptr), intent(inout) :: this

    if (associated(this%val)) deallocate(this%val)

  end subroutine ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_camp_state
