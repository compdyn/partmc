! Copyright (C) 2017 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

! TODO Figure out a way to avoid type-specific value getter functions
! Seems like this should be possible with class(*)?

!> \file
!> The pmc_property module.

!> The property_t structure and associated subroutines.
module pmc_property

  use pmc_constants,                only : i_kind, dp
  use pmc_util,                     only : die_msg

  implicit none
  private

  public :: property_t, property_link_t

  !> Property data
  !!
  !! A set of physical properties, sub-model parameters and similar constants
  !! related to a chemical species, reaction, or other data object.
  type property_t
    private
    !> Number of elements
    integer(kind=i_kind) :: num_elem = 0
    !> First element in the set
    type(property_link_t), pointer :: first_link => null()
    !> Last element in the set
    type(property_link_t), pointer :: last_link => null()
    !> Iterator
    type(property_link_t), pointer :: curr_link => null()
  contains
    !> Put a value in the data set
    procedure :: put => pmc_property_put
    !> Find a key-value pair by key name
    procedure :: get => pmc_property_get
    !> Get the number of key-value pairs
    procedure :: size => pmc_property_size
    !> Reset the iterator
    procedure :: iter_reset => pmc_property_iter_reset
    !> Get the current key-value pair
    procedure :: iter_curr => pmc_property_iter_curr
    !> Increment the iterator
    procedure :: iter_next => pmc_property_iter_next
    !> Finalize
    final :: pmc_property_final
  end type property_t

  !> Constructor for property_t
  interface property_t
    procedure :: pmc_property_constructor
  end interface property_t

  !> Property link data
  !!
  !! An element of a property data set. Property values can be of any
  !! primitive type or be a pointer to a sub-set of property data
  type property_link_t
    private
    !> Key name
    character(:), allocatable :: key_name
    !> Value
    class(*), pointer :: val => null()
    !> Next link
    type(property_link_t), pointer :: next_link => null()
  contains
    !> Get the key name
    procedure :: key => pmc_property_link_key
    !> Set the value
    procedure, private :: set_value => pmc_property_link_set_value
    !> Get the int value
    procedure :: value_int => pmc_property_link_value_int
    !> Get the real value
    procedure :: value_real => pmc_property_link_value_real
    !> Get the logical value
    procedure :: value_logical => pmc_property_link_value_logical
    !> Get the string value
    procedure :: value_string => pmc_property_link_value_string
    !> Get the property sub-set
    procedure :: value_property_t => pmc_property_link_value_property_t
    !> Finalize
    final :: pmc_property_link_final
  end type property_link_t

  !> Constructor for link
  interface property_link_t
    procedure pmc_property_link_constructor
  end interface property_link_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for property_t
  function pmc_property_constructor() result(new_obj)

    !> A new property set
    type(property_t), pointer :: new_obj

    allocate(new_obj)

  end function pmc_property_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Put an element in the property data set
  subroutine pmc_property_put(this, key, val)

    !> Property data set
    class(property_t), intent(inout) :: this
    !> New key
    character(len=*), intent(in) :: key
    !> New value
    class(*), intent(in) :: val

    !> New link
    type(property_link_t), pointer :: new_link

    new_link => this%get(key)

    if (associated(new_link)) then
      if (associated(new_link%val)) deallocate(new_link%val)
      call new_link%set_value(val)
      return
    end if

    new_link => property_link_t(key, val)
    if (.not.associated(this%first_link)) then
      this%first_link => new_link
      this%last_link => this%first_link
    else
      this%last_link%next_link => new_link
      this%last_link => new_link
    end if

    this%num_elem = this%num_elem + 1

  end subroutine pmc_property_put

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find a key-value pair by key name. Returns a null pointer if the key name
  !! is not found.
  function pmc_property_get(this, key) result(found_pair)

    !> Pointer to property key-value pair
    type(property_link_t), pointer :: found_pair
    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name to search for
    character(len=*), intent(in) :: key

    type(property_link_t), pointer :: curr_link

    found_pair => null()
    curr_link => this%first_link
    do while (associated(curr_link))
      if (key .eq. curr_link%key()) then
        found_pair => curr_link
        return
      end if
      curr_link => curr_link%next_link
    end do

  end function pmc_property_get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of elements in the property set
  function pmc_property_size(this) result(size)

    !> Number of elements in the property set
    integer(kind=i_kind) :: size
    !> Property dataset
    class(property_t), intent(in) :: this
    
    type(property_link_t), pointer :: curr_link

    size = 0
    curr_link => this%first_link
    do while (associated(curr_link))
      size = size + 1
      curr_link => curr_link%next_link
    end do 

  end function pmc_property_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the iterator. It will now point to the first property in the
  !! dataset, or be NULL in the case of an empty property dataset
  subroutine pmc_property_iter_reset(this)

    !> Property dataset
    class(property_t), intent(inout) :: this

    this%curr_link => this%first_link

  end subroutine pmc_property_iter_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a pointer to the current property from the interator
  function pmc_property_iter_curr(this) result(curr_link)

    !> Pointer to current property key-value pair
    type(property_link_t), pointer :: curr_link
    !> Property dataset
    class(property_t), intent(in) :: this

    curr_link => this%curr_link

  end function pmc_property_iter_curr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Increment the interator
  subroutine pmc_property_iter_next(this)

    !> Property dataset
    class(property_t), intent(inout) :: this

    this%curr_link => this%curr_link%next_link

  end subroutine pmc_property_iter_next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a property_t variable
  elemental subroutine pmc_property_final(this)

    !> Property dataset
    type(property_t), intent(inout) :: this
    type(property_link_t), pointer :: curr, next

    next => this%first_link
    do while (associated(next))
      curr => next
      next => curr%next_link 
      deallocate(curr)
    end do

  end subroutine pmc_property_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for property_link_t
  function pmc_property_link_constructor(key, val) result(new_obj)

    !> Pointer to new property key-value pair
    type(property_link_t), pointer :: new_obj
    !> Key name
    character(len=*), intent(in) :: key
    !> New value
    class(*), intent(in) :: val

    allocate(new_obj)
    new_obj%key_name = trim(key)
    new_obj%next_link => null()
    call new_obj%set_value(val)

  end function pmc_property_link_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the key name of a property
  function pmc_property_link_key(this) result(key)

    !> Key name
    character(:), allocatable :: key
    !> Property key-value pair
    class(property_link_t), intent(in) :: this
    
    key = this%key_name

  end function pmc_property_link_key

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the value of a property key-value pair
  subroutine pmc_property_link_set_value(this, val)

    !> Property key-value pair
    class(property_link_t), intent(inout) :: this
    !> New value
    class(*), intent(in) :: val

    select type (val)
      type is (integer(kind=i_kind))
      type is (real(kind=dp))
      type is (logical)
      type is (character(len=*))
      type is (property_t)
      class default
        call die_msg(555557200,"Invalid property data type")
    end select

    allocate(this%val, source=val)

  end subroutine pmc_property_link_set_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the int value of a property
  function pmc_property_link_value_int(this) result(val)

    !> Value
    integer(kind=i_kind) :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type(this_val)
      type is (integer(kind=i_kind))
        val = this_val
    end select

  end function pmc_property_link_value_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the real value of a property
  function pmc_property_link_value_real(this) result(val)

    !> Value
    real(kind=dp) :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type(this_val)
      type is (real(kind=dp))
        val = this_val
      class default
        call die_msg(151463892, "Property type mismatch")
    end select

  end function pmc_property_link_value_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the logical value of a property
  function pmc_property_link_value_logical(this) result(val)

    !> Value
    logical :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type(this_val)
      type is (logical)
        val = this_val
      class default
        call die_msg(371288570, "Property type mismatch")
    end select

  end function pmc_property_link_value_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the string value of a property
  function pmc_property_link_value_string(this) result(val)

    !> Value
    character(:), pointer :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type(this_val)
      type is (character(len=*))
        val => this_val
      class default
        call die_msg(192508586, "Property type mismatch")
    end select

  end function pmc_property_link_value_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the property_t value of a property
  function pmc_property_link_value_property_t(this) result(val)

    !> Value
    type(property_t), pointer :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type(this_val)
      type is (property_t)
        val => this_val
      class default
        call die_msg(641781966, "Property type mismatch")
    end select

  end function pmc_property_link_value_property_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the property_link_t variable
  elemental subroutine pmc_property_link_final(this)

    !> Property key-value pair
    type(property_link_t), intent(inout) :: this

    class(*), pointer :: this_val

    if (associated(this%val)) then
      this_val => this%val
      select type (this_val)
        type is (property_t)
          deallocate(this_val)
      end select
    end if

    if (allocated(this%key_name)) deallocate(this%key_name)

  end subroutine pmc_property_link_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_property
