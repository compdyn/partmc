! Copyright (C) 2017 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_property module.

!> The property_t structure and associated subroutines.
module pmc_property

#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_constants,                only : i_kind, dp
  use pmc_util,                     only : die_msg, warn_msg, to_string, string_t


  implicit none
  private

  public :: property_t

  !> Property data
  !!
  !! A set of physical properties, sub-model parameters and similar constants
  !! related to a chemical species, reaction, or other data object. The \c
  !! property_t type can be used to build a set of data with a \c json -like
  !! structure.
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
    !> Load input data
    procedure :: load
    !> Put a value in the data set
    procedure :: put
    !> Get the current key name
    procedure :: get_key
    !> Get an integer value
    procedure :: get_int
    !> Get a real value
    procedure :: get_real
    !> Get a logical value
    procedure :: get_logical
    !> Get a string value
    procedure :: get_string
    !> Get a sub-set of properties
    procedure :: get_property_t
    !> Get the number of key-value pairs
    procedure :: size => get_size
    !> Reset the iterator
    procedure :: iter_reset
    !> Increment the iterator
    procedure :: iter_next
    !> Move property data from one property_t instance to another
    procedure :: move
    !> Update this property_t instance with data from another
    procedure :: update
    !> Print the contents of a property set
    procedure :: print => do_print
    !> Finalize
    final :: finalize
    
    !> Private functions
    !> Find a key-value pair by key name
    procedure, private :: get
  end type property_t

  ! Constructor for property_t
  interface property_t
    procedure :: constructor
  end interface property_t

  !> Property link data
  !!
  !! An element of a property data set. Property values can be of any
  !! primitive type or be a pointer to a sub-set of property data.
  !! The property_link_t object is for internal use in the pmc_property
  !! module. All interactions with property data should be made using
  !! property_t objects.
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
    procedure :: key
    !> Set the value
    procedure :: set_value
    !> Get the int value
    procedure :: value_int
    !> Get the real value
    procedure :: value_real
    !> Get the logical value
    procedure :: value_logical
    !> Get the string value
    procedure :: value_string
    !> Get the property sub-set
    procedure :: value_property_t
    !> Print the contents of a property key-value pair
    procedure :: print => link_do_print
    !> Finalize
    final :: link_finalize
  end type property_link_t

  ! Constructor for link
  interface property_link_t
    procedure link_constructor
  end interface property_link_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for property_t
  function constructor() result(new_obj)

    !> A new property set
    type(property_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a property set from input data
#ifdef PMC_USE_JSON
  recursive subroutine load(this, json, j_obj, as_object)

    !> Property dataset
    class(property_t), intent(inout) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj
    !> Set to true if j_obj is a json object to parse, adding all child 
    !! key-value pairs to the data set, or false if j_obj is a single
    !! key-value pair to add to the data set
    logical, intent(in) :: as_object

    type(json_value), pointer :: child, next
    type(property_t), pointer :: sub_prop
    character(kind=json_ck, len=:), allocatable :: unicode_prop_key
    character(len=:), allocatable :: prop_key

    character(kind=json_ck, len=:), allocatable :: unicode_val
    character(len=:), allocatable :: str_val
    logical(json_lk)     :: bool_val
    real(json_rk)        :: real_val
    integer(json_ik)     :: int_val

    integer(json_ik)     :: var_type

    ! initialize pointer to next object to parse
    next => null()

    ! determine whether to add parent or children key-value pairs
    if (as_object) then 
      call json%get_child(j_obj, child)
    else 
      child => j_obj
    end if
    
    ! loop through set of json objects to add to the property set
    do while (associated(child))

      ! get the key name and value
      call json%info(child, name=unicode_prop_key, var_type=var_type)
      prop_key = unicode_prop_key

      ! add key-value pair of appropriate type
      select case (var_type)

        ! skip null objects
        case (json_null)

        ! integer
        case (json_integer)
          call json%get(child, int_val)
          call this%put(prop_key, int(int_val, i_kind))
       
        ! double
        case (json_double)
          call json%get(child, real_val)
          call this%put(prop_key, real(real_val, dp))
        
        ! boolean
        case (json_logical)
          call json%get(child, bool_val)
          call this%put(prop_key, logical(bool_val))
        
        ! string
        case (json_string)
          call json%get(child, unicode_val)
          str_val = unicode_val
          call this%put(prop_key, str_val)
        
        ! sub-set of key-value pairs
        case (json_object)
          sub_prop => property_t()
          call sub_prop%load(json, child, .true.)
          call this%put(prop_key, sub_prop)
          deallocate(sub_prop)
        
        ! sub-set of values
        case (json_array)
          sub_prop => property_t()
          call sub_prop%load(json, child, .true.)
          call this%put(prop_key, sub_prop)
          deallocate(sub_prop)
        
        ! skip other types
        case default
      end select

      ! get the next object to add
      if (as_object) call json%get_next(child, next)
      child => next
    
    end do

#else
  subroutine load(this)

    !> Property dataset
    class(property_t), intent(inout) :: this

    call warn_msg(733896496, "No support for input files.")
#endif
  end subroutine load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Put an element in the property data set
  recursive subroutine put(this, key, val)

    !> Property data set
    class(property_t), intent(inout) :: this
    !> New key
    character(len=:), allocatable, intent(in) :: key
    !> New value
    class(*), intent(in) :: val

    type(property_link_t), pointer :: new_link, sub_link
    type(property_t), allocatable :: sub_prop_set
    class(*), pointer :: curr_val

    ! if this is an array element, the key will be empty
    if (allocated(key).and.len(key).ge.1) then
   
      ! look for the key in the existing properties
      new_link => this%get(key)

      ! do not allow overwrites of existing properties, but allow sub-sets
      ! of properties to be appended
      if (associated(new_link)) then
        curr_val => new_link%val
        select type (curr_val)
        class is (property_t)
          select type (val)
          class is (property_t)
            sub_link => val%first_link
            do while (associated(sub_link))
              call curr_val%put(sub_link%key_name, sub_link%val)
              sub_link => sub_link%next_link
            end do
          class default
            call die_msg(698012538, "Property type mismatch for "//key)
          end select
        class default
          call die_msg(359604264, "Trying to overwrite property "//key)
        end select
        return
      end if

    end if

    ! create a new link. for property_t sub-property sets,
    ! copy the passed value to a new object
    select type (val)
    class is (property_t)
      allocate(sub_prop_set)
      sub_link => val%first_link
      do while (associated(sub_link))
        call sub_prop_set%put(sub_link%key_name, sub_link%val)
        sub_link => sub_link%next_link
      end do
      new_link => property_link_t(key, sub_prop_set)
      sub_prop_set%first_link => null()
      sub_prop_set%last_link => null()
      deallocate(sub_prop_set)
    class default
      new_link => property_link_t(key, val)
    end select

    ! if the key does not exist in the property dataset,
    ! create a new link to add it.
    if (.not.associated(this%first_link)) then
      this%first_link => new_link
      this%last_link => this%first_link
    else
      this%last_link%next_link => new_link
      this%last_link => new_link
    end if

    this%num_elem = this%num_elem + 1

  end subroutine put

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the key name of the element currently pointed to by the iterator. 
  !! Returns true if the iterator points to a key-value pair; false indicates
  !! the list is empty, the iterator was never reset, or the end of the list
  !! has been reached. Array elements return true, but have an empty key name.
  logical function get_key(this, key) result (found)

    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name
    character(len=:), allocatable, intent(out) :: key

    found = .false.
    if (.not.associated(this%curr_link)) return
    key = this%curr_link%key()
    found = .true.

  end function get_key

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an integer value. The return value is true if the key-value pair
  !! was found, and false otherwise. If no key name is specified, the current
  !! value of the iterator is returned. In this case true indicates a current 
  !! key-value exists; false indicates the list is empty, the iterator was
  !! never reset, or the end of the list has been reached.
  logical function get_int(this, key, val) result(found)

    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name to search for
    character(len=:), allocatable, intent(in), optional :: key
    !> Property value
    integer(kind=i_kind), intent(out) :: val

    type(property_link_t), pointer :: link

    found = .false.
    if (present(key)) then
      link => get(this, key)
      if (.not. associated(link)) return
      val = link%value_int()
    else
      if (.not.associated(this%curr_link)) return
      val = this%curr_link%value_int()
    endif
    found = .true.

  end function get_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a real value. The return value is true if the key-value pair
  !! was found, and false otherwise. If no key name is specified, the current
  !! value of the iterator is returned. In this case true indicates a current 
  !! key-value exists; false indicates the list is empty, the iterator was
  !! never reset, or the end of the list has been reached.
  logical function get_real(this, key, val) result(found)

    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name to search for
    character(len=:), allocatable, intent(in), optional :: key
    !> Property value
    real(kind=dp), intent(out) :: val

    type(property_link_t), pointer :: link

    found = .false.
    if (present(key)) then
      link => get(this, key)
      if (.not. associated(link)) return
      val = link%value_real()
    else
      if (.not.associated(this%curr_link)) return
      val = this%curr_link%value_real()
    endif
    found = .true.

  end function get_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a logical value. The return value is true if the key-value pair
  !! was found, and false otherwise. If no key name is specified, the current
  !! value of the iterator is returned. In this case true indicates a current 
  !! key-value exists; false indicates the list is empty, the iterator was
  !! never reset, or the end of the list has been reached.
  logical function get_logical(this, key, val) result(found)

    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name to search for
    character(len=:), allocatable, intent(in), optional :: key
    !> Property value
    logical, intent(out) :: val

    type(property_link_t), pointer :: link

    found = .false.
    if (present(key)) then
      link => get(this, key)
      if (.not. associated(link)) return
      val = link%value_logical()
    else
      if (.not.associated(this%curr_link)) return
      val = this%curr_link%value_logical()
    endif
    found = .true.

  end function get_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a string value. The return value is true if the key-value pair
  !! was found, and false otherwise. If no key name is specified, the current
  !! value of the iterator is returned. In this case true indicates a current 
  !! key-value exists; false indicates the list is empty, the iterator was
  !! never reset, or the end of the list has been reached.
  logical function get_string(this, key, val) result(found)

    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name to search for
    character(len=:), allocatable, intent(in), optional :: key
    !> Property value
    character(len=:), allocatable, intent(out) :: val

    type(property_link_t), pointer :: link

    found = .false.
    if (present(key)) then
      link => get(this, key)
      if (.not. associated(link)) return
      val = link%value_string()
    else
      if (.not.associated(this%curr_link)) return
      val = this%curr_link%value_string()
    endif
    found = .true.

  end function get_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a property sub-set. The return value is true if the key-value pair
  !! was found, and false otherwise. If no key name is specified, the current
  !! value of the iterator is returned. In this case true indicates a current 
  !! key-value exists; false indicates the list is empty, the iterator was
  !! never reset, or the end of the list has been reached.
  logical function get_property_t(this, key, val) result(found)

    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name to search for
    character(len=:), allocatable, intent(in), optional :: key
    !> Property value
    type(property_t), pointer, intent(out) :: val

    type(property_link_t), pointer :: link

    found = .false.
    val => null()
    if (present(key)) then
      link => get(this, key)
      if (.not. associated(link)) return
      val => link%value_property_t()
    else
      if (.not. associated(this%curr_link)) return
      val => this%curr_link%value_property_t()
    end if
    found = .true.

  end function get_property_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of elements in the property set
  function get_size(this)

    !> Number of elements in the property set
    integer(kind=i_kind) :: get_size
    !> Property dataset
    class(property_t), intent(in) :: this
    
    type(property_link_t), pointer :: curr_link

    get_size = 0
    curr_link => this%first_link
    do while (associated(curr_link))
      get_size = get_size + 1
      curr_link => curr_link%next_link
    end do 

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the iterator. It will now point to the first property in the
  !! dataset, or be NULL in the case of an empty property dataset
  subroutine iter_reset(this)

    !> Property dataset
    class(property_t), intent(inout) :: this

    this%curr_link => this%first_link

  end subroutine iter_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Increment the interator
  subroutine iter_next(this)

    !> Property dataset
    class(property_t), intent(inout) :: this

    if (associated(this%curr_link)) then
      this%curr_link => this%curr_link%next_link
    else
      call warn_msg(365476096, "Trying to iterate NULL iterator.")
    end if

  end subroutine iter_next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Move data from one property_t instance to another
  elemental subroutine move(this, dest)

    !> Property dataset to move
    class(property_t), intent(inout) :: this
    !> Property dataset destination
    class(property_t), intent(inout) :: dest

    dest%first_link => this%first_link
    dest%curr_link => this%curr_link
    dest%last_link => this%last_link
    dest%num_elem = this%num_elem
    this%first_link => null()
    this%curr_link => null()
    this%last_link => null()
    this%num_elem = 0

  end subroutine move

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update this property_t instance with data from another instance
  subroutine update(this, source)

    !> Property dataset to update
    class(property_t), intent(inout) :: this
    !> Property dataset to update from
    class(property_t), intent(inout) :: source

    type(property_link_t), pointer :: curr_prop

    curr_prop => source%first_link
    do while (associated(curr_prop))
      call this%put(curr_prop%key_name, curr_prop%val)
      curr_prop => curr_prop%next_link
    end do

  end subroutine update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the contents of a property set
  recursive subroutine do_print(this, file_unit)

    !> Property dataset
    class(property_t), intent(in) :: this
    !> File unit for output
    integer(kind=i_kind), optional, intent(in) :: file_unit

    type(property_link_t), pointer :: curr_link
    integer(kind=i_kind) :: f_unit = 6

    if (present(file_unit)) f_unit = file_unit

    curr_link => this%first_link
    do while (associated(curr_link))
      if (associated(curr_link%next_link)) then
        call curr_link%print(",", f_unit)
      else
        call curr_link%print("", f_unit)
      endif
      curr_link => curr_link%next_link
    end do

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a property_t variable
  elemental subroutine finalize(this)

    !> Property dataset
    type(property_t), intent(inout) :: this

    type(property_link_t), pointer :: next

    next => null()
    do while (associated(this%first_link))
      next => this%first_link%next_link 
      deallocate(this%first_link)
      this%first_link => next
    end do

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Find a key-value pair by key name. Returns a null pointer if the key name
  !! is not found.
  function get(this, key) result(found_pair)

    !> Pointer to property key-value pair
    type(property_link_t), pointer :: found_pair
    !> Property dataset
    class(property_t), intent(in) :: this
    !> Key name to search for
    character(len=:), allocatable, intent(in) :: key

    type(property_link_t), pointer :: curr_link

    found_pair => null()
    if (.not. associated(this%first_link)) return
    curr_link => this%first_link
    do while (associated(curr_link))
      if (key .eq. curr_link%key()) then
        found_pair => curr_link
        return
      end if
      curr_link => curr_link%next_link
    end do

  end function get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! property_link_t functions
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for property_link_t
  function link_constructor(key, val) result(new_obj)

    !> Pointer to new property key-value pair
    type(property_link_t), pointer :: new_obj
    !> Key name
    character(len=:), allocatable, intent(in) :: key
    !> New value
    class(*), intent(in) :: val

    allocate(new_obj)
    new_obj%key_name = trim(key)
    new_obj%next_link => null()
    call new_obj%set_value(val)

  end function link_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the key name of a property
  function key(this)

    !> Key name
    character(:), allocatable :: key
    !> Property key-value pair
    class(property_link_t), intent(in) :: this
    
    key = this%key_name

  end function key

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the value of a property key-value pair
  subroutine set_value(this, val)

    !> Property key-value pair
    class(property_link_t), intent(inout) :: this
    !> New value
    class(*), intent(in) :: val

    type(string_t), pointer :: str_val

    ! determine the value type
    select type(val)
      
      ! add integers, reals, logicals, and string_t as-is
      type is (integer(kind=i_kind))
      type is (real(kind=dp))
      type is (logical)
      type is (string_t)

      ! handle empty sub-sets
      class is (property_t)
        if (.not.associated(val%first_link)) then
          this%val => property_t()
          return
        end if

      ! convert character arrays to string_t objects
      type is (character(len=*))
        allocate(str_val)
        str_val%string = val
        this%val => str_val
        return

      ! error on unsupported types
      class default
        call die_msg(728532218, "Unsupported property type")
    end select

    allocate(this%val, source=val)

  end subroutine set_value

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the int value of a property
  function value_int(this) result(val)

    !> Value
    integer(kind=i_kind) :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type(this_val)
      type is (integer(kind=i_kind))
        val = this_val
      class default
        call die_msg(509101133, "Property type mismatch for key "//&
                trim(this%key_name))
    end select

  end function value_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the real value of a property
  function value_real(this) result(val)

    !> Value
    real(kind=dp) :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type(this_val)
      type is (integer(kind=i_kind))
        val = real(this_val, kind=dp)
      type is (real(kind=dp))
        val = this_val
      class default
        call die_msg(151463892, "Property type mismatch for key "//&
                trim(this%key_name))
    end select

  end function value_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the logical value of a property
  function value_logical(this) result(val)

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
        call die_msg(371288570, "Property type mismatch for key "//&
                trim(this%key_name))
    end select

  end function value_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the string value of a property
  function value_string(this) result(val)

    !> Value
    character(:), pointer :: val
    !> Property key-value pair
    class(property_link_t), intent(in) :: this

    class(*), pointer :: this_val

    this_val => this%val
    select type (this_val)
      type is (string_t)
        val => this_val%string
      class default
        call die_msg(153505401, "Property type mismatch for key "//&
                trim(this%key_name))
    end select

  end function value_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the property_t value of a property
  function value_property_t(this) result(val)

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
        call die_msg(641781966, "Property type mismatch for key "//&
                trim(this%key_name))
    end select

  end function value_property_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the contents of a property key-value pair
  recursive subroutine link_do_print(this, suffix, file_unit)

    !> Property key-value pair
    class(property_link_t), intent(in) :: this
    !> Text to append to the end of the line
    character(len=*), intent(in) :: suffix
    !> File unit for output
    integer(kind=i_kind), optional, intent(in) :: file_unit

    class(*), pointer :: val
    integer(kind=i_kind) :: f_unit = 6

    if (present(file_unit)) f_unit = file_unit

    val => this%val
    select type(val)
      type is (integer(kind=i_kind))
        write(f_unit,*) '"'//this%key_name//'" : '//trim(to_string(val))// &
                suffix
      type is (real(kind=dp))
        write(f_unit,*) '"'//this%key_name//'" : '//trim(to_string(val))// &
                suffix
      type is (logical)
        write(f_unit,*) '"'//this%key_name//'" : '//trim(to_string(val))// &
                suffix
      type is (string_t)
        write(f_unit,*) '"'//this%key_name//'" : "'//val%string//'"'//suffix
      class is (property_t)
        write(f_unit,*) '"'//this%key_name//'" : {'
        call val%print(f_unit)
        write(f_unit,*) '}'//suffix
      class default
        call die_msg(711028956, "Unsupported property type")
    end select

  end subroutine link_do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the property_link_t variable
  elemental subroutine link_finalize(this)

    !> Property key-value pair
    type(property_link_t), intent(inout) :: this

    if (associated(this%val)) deallocate(this%val)

  end subroutine link_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_property
