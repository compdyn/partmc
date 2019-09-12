! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_factory module.

!> \page camp_aero_rep_add CAMP: Adding an Aerosol Representation
!!
!! NOTE: These instructions are out-of-date. TODO update
!!
!! Adding an \ref camp_aero_rep "aerosol representation" to the \ref
!! camp_chem "camp-chem" module can be done in the following steps:
!!
!! ## Step 1. Create a new aerosol representation module ##
!!   The module should be placed in the \c /src/aero_reps folder and
!!   extend the abstract \c pmc_aero_rep_data::aero_rep_data_t type,
!!   overriding all deferred functions, and providing a constructor that
!!   returns a pointer to a newly allocated instance of the new type:
!!
!! \code{.f90}
!! module pmc_aero_rep_foo
!!
!!   use ...
!!
!!   implicit none
!!   private
!!
!!   public :: aero_rep_foo_t
!!
!!   type, extends(aero_rep_data_t) :: aero_rep_foo_t
!!   contains
!!     ... (all deferred functions) ...
!!   end type :: aero_rep_foo_t
!!
!!   ! Constructor
!!   interface aero_rep_foo_t
!!     procedure :: constructor
!!   end interface aero_rep_foo_t
!!
!!   function constructor() result (new_obj)
!!     type(aero_rep_foo_t), pointer :: new_obj
!!     allocate(new_obj)
!!   end function constructor
!!
!!   ...
!!
!! end module pmc_aero_rep_foo
!! \endcode
!!
!! ## Step 2. Create a new aerosol representation state module ##
!!   The module should also be placed in the \c /src/aero_reps folder and
!!   should extend the abstract \c pmc_aero_rep_state::aero_rep_state_t
!!   type, overriding all deferred functions, and providing a constructor
!!   that returns a pointer to a newly allocated instance of the new type:
!!
!! \code{.f90}
!! module pmc_aero_rep_foo_state
!!
!!   use ...
!!
!!   implicit none
!!   private
!!
!!   public :: aero_rep_foo_state_t
!!
!!   type, extends(aero_rep_state_t) :: aero_rep_foo_state_t
!!   contains
!!     ... (all deferred functions) ...
!!   end type :: aero_rep_foo_state_t
!!
!!   ! Constructor
!!   interface aero_rep_foo_state_t
!!     procedure :: constructor
!!   end interface aero_rep_foo_state_t
!!
!!   function constructor() result (new_obj)
!!     type(aero_rep_foo_state_t), pointer :: new_obj
!!     allocate(new_obj)
!!   end function constructor
!!
!!   ...
!!
!! end module pmc_aero_rep_foo_state
!! \endcode
!!
!! ## Step 3. Add the aerosol representation to the \c  pmc_aero_rep_factory module ##
!!
!! \code{.f90}
!! module pmc_aero_rep_factory
!!
!! ...
!!
!!  ! Use all aerosol representation modules
!!  ...
!!  use pmc_aero_rep_foo
!!
!!  ...
!!
!!  !> Identifiers for aerosol representations - used by binary packing/unpacking
!!  !! functions
!!  ...
!!  integer(kind=i_kind), parameter :: AERO_REP_FOO = 32
!!
!!  ...
!!
!!  !> Create a new aerosol representation by type name
!!  function create(this, type_name) result (new_obj)
!!    ...
!!    select case (type_name)
!!      ...
!!      case ("AERO_REP_FOO")
!!        new_obj => aero_rep_foo_t()
!!    ...
!!  end function create
!!
!! ...
!!
!!  !> Pack the given value to the buffer, advancing position
!!  subroutine bin_pack(this, aero_rep, buffer, pos)
!!    ...
!!    select type (aero_rep)
!!      ...
!!      type is (aero_rep_foo_t)
!!        aero_rep_type = AERO_REP_FOO
!!    ...
!!  end subroutine bin_pack
!!
!! ...
!!
!!  !> Unpack the given value to the buffer, advancing position
!!  function bin_unpack(this, buffer, pos) result (aero_rep)
!!    ...
!!    select case (aero_rep_type)
!!      ...
!!      case (AERO_REP_FOO)
!!        aero_rep => aero_rep_foo_t()
!!    ...
!!  end function bin_unpack
!!
!! ...
!! end module pmc_aero_rep_factory
!! \endcode
!!
!! ## Step 4. Add the new module to the CMakeList file in the root directory. ##
!!
!! \code{.unparsed}
!! ...
!!
!! # partmc library
!! ...
!! set(AEROSOL_REPS_SRC
!!   ...
!!   src/aero_reps/aero_rep_foo.F90
!!   src/aero_reps/aero_rep_foo_state.F90
!! )
!!
!! ...
!! \endcode
!!
!! ## Step 5. Add unit tests for the new \c aero_rep_foo_t and \c aero_rep_foo_state_t types ##
!!
!! Unit testing should be based on the actual functions of the new module, but
!! in general 80% code coverage is recommended. Some examples can be found in
!! the \c /src/test folder
!!
!! ## Usage ##
!! The new \ref camp_aero_rep "aerosol representation" is now ready to use.
!! To include it in a \c pmc_camp_core::camp_core_t instance, add an \ref
!! input_format_aero_rep "aerosol representation object" to a new or existing
!! \ref input_format_camp_config "camp-chem configuration file" with a
!! \b type corresponding to the newly created type, along with any required
!! parameters:
!!
!! \code{.json}
!! { "pmc-data" : [
!!   {
!!     "type" : "AERO_REP_FOO",
!!     ...
!!   },
!!   ...
!! ]}
!! \endcode
!!

!> The aero_rep_factory_t type and associated subroutines.
module pmc_aero_rep_factory

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_aero_rep_data
  use pmc_constants,                  only : i_kind, dp
  use pmc_mpi
  use pmc_util,                       only : die_msg, string_t, assert_msg, &
                                             warn_msg

  ! Use all aerosol representation modules
  use pmc_aero_rep_modal_binned_mass
  use pmc_aero_rep_single_particle

  use iso_c_binding

  implicit none
  private

  public :: aero_rep_factory_t

  !> Identifiers for aerosol representations - used by binary
  !! packing/unpacking functions
  integer(kind=i_kind), parameter :: AERO_REP_SINGLE_PARTICLE   = 1
  integer(kind=i_kind), parameter :: AERO_REP_MODAL_BINNED_MASS = 2

  !> Factory type for aerosol representations
  !!
  !! Provides new instances of types extending aero_rep_data_t by name or
  !! from input file data
  type :: aero_rep_factory_t
  contains
    !> Create a new aerosol representation by type name
    procedure :: create
    !> Create a new aerosol representation from input data
    procedure :: load
    !> Get the aerosol representation type
    procedure :: get_type
    !> Initialize an update data object
    procedure :: initialize_update_data
    !> Determine the number of bytes required to pack an given aerosol
    !! representation
    procedure :: pack_size
    !> Pack a given aerosol representation to the buffer, advancing the
    !! position
    procedure :: bin_pack
    !> Unpack a aerosol representation from the buffer, advancing the position
    procedure :: bin_unpack
  end type aero_rep_factory_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Create a new aerosol representation by type name
  function create(this, type_name) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Name of the aerosol representation
    character(len=*), intent(in) :: type_name

    new_obj => null()

    select case (type_name)
      case ("AERO_REP_MODAL_BINNED_MASS")
        new_obj => aero_rep_modal_binned_mass_t()
      case ("AERO_REP_SINGLE_PARTICLE")
        new_obj => aero_rep_single_particle_t()
      case default
        call die_msg(792930166, "Unknown aerosol representation type: " &
                //type_name)
    end select

  end function create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load an aerosol represenation based on its type
#ifdef PMC_USE_JSON
  function load(this, json, j_obj) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> JSON core
    type(json_core), pointer, intent(in) :: json
    !> JSON object
    type(json_value), pointer, intent(in) :: j_obj

    character(kind=json_ck, len=:), allocatable :: unicode_name
    character(len=:), allocatable :: type_name
    logical(kind=json_lk) :: found

    new_obj => null()

    ! get the aerosol representation type
    call json%get(j_obj, "type", unicode_name, found)
    call assert_msg(904082395, found, 'Missing aerosol representation type.')
    type_name = unicode_name

    ! create a new aerosol representation
    new_obj => this%create(type_name)

    ! load aerosol representation parameters from the json object
    call new_obj%load(json, j_obj)

#else
  function load(this) result (new_obj)

    !> An aerosol representation
    class(aero_rep_data_t), pointer :: new_obj
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this

    new_obj => null()

    call warn_msg(723960750, "No support for input files.")
#endif
  end function load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the aerosol reaction type
  integer(kind=i_kind) function get_type(this, aero_rep) result(aero_rep_type)

    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Aerosol representation to get the type of
    class(aero_rep_data_t), intent(in) :: aero_rep

    select type (aero_rep)
      type is (aero_rep_modal_binned_mass_t)
        aero_rep_type = AERO_REP_MODAL_BINNED_MASS
      type is (aero_rep_single_particle_t)
        aero_rep_type = AERO_REP_SINGLE_PARTICLE
      class default
        call die_msg(865927801, "Unknown aerosol representation type")
    end select

  end function get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a new update data object
  subroutine initialize_update_data(this, aero_rep, update_data)

    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Aerosol representation to be updated
    class(aero_rep_data_t), intent(inout) :: aero_rep
    !> Update data object
    class(aero_rep_update_data_t), intent(out) :: update_data

    select type (update_data)
      type is (aero_rep_update_data_single_particle_radius_t)
        select type (aero_rep)
          type is (aero_rep_single_particle_t)
            call aero_rep%update_data_initialize_radius(update_data, &
                                                 AERO_REP_SINGLE_PARTICLE)
          class default
            call die_msg(761019956, "Update data <-> aero rep mismatch")
        end select
      type is (aero_rep_update_data_single_particle_number_t)
        select type (aero_rep)
          type is (aero_rep_single_particle_t)
            call aero_rep%update_data_initialize_number(update_data, &
                                                 AERO_REP_SINGLE_PARTICLE)
          class default
            call die_msg(584145506, "Update data <-> aero rep mismatch")
        end select
      type is (aero_rep_update_data_modal_binned_mass_GMD_t)
        select type (aero_rep)
          type is (aero_rep_modal_binned_mass_t)
            call aero_rep%update_data_initialize_GMD(update_data, &
                                                 AERO_REP_MODAL_BINNED_MASS)
          class default
            call die_msg(696463851, "Update data <-> aero rep mismatch")
        end select
      type is (aero_rep_update_data_modal_binned_mass_GSD_t)
        select type (aero_rep)
          type is (aero_rep_modal_binned_mass_t)
            call aero_rep%update_data_initialize_GSD(update_data, &
                                                 AERO_REP_MODAL_BINNED_MASS)
          class default
            call die_msg(526306947, "Update data <-> aero rep mismatch")
        end select
      class default
        call die_msg(916635086, "Internal error - update data type missing.")
    end select

  end subroutine initialize_update_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack an aerosol representation
  integer(kind=i_kind) function pack_size(this, aero_rep, comm)

    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Aerosol representation to pack
    class(aero_rep_data_t), intent(in) :: aero_rep
    !> MPI communicator
    integer, intent(in) :: comm

    pack_size =  pmc_mpi_pack_size_integer(int(1, kind=i_kind), comm) + &
                 aero_rep%pack_size(comm)

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, aero_rep, buffer, pos, comm)

    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Aerosol representation to pack
    class(aero_rep_data_t), intent(in) :: aero_rep
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: aero_rep_type, i_aero_rep, prev_position

    prev_position = pos
    select type (aero_rep)
      type is (aero_rep_modal_binned_mass_t)
        aero_rep_type = AERO_REP_MODAL_BINNED_MASS
      type is (aero_rep_single_particle_t)
        aero_rep_type = AERO_REP_SINGLE_PARTICLE
      class default
        call die_msg(278244560, &
                "Trying to pack aerosol representation of unknown type.")
    end select
    call pmc_mpi_pack_integer(buffer, pos, aero_rep_type, comm)
    call aero_rep%bin_pack(buffer, pos, comm)
    call assert(897674844, &
         pos - prev_position <= this%pack_size(aero_rep, comm))
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpack the given value to the buffer, advancing position
  function bin_unpack(this, buffer, pos, comm) result (aero_rep)

    !> Unpacked aerosol representation
    class(aero_rep_data_t), pointer :: aero_rep
    !> Aerosol representation factory
    class(aero_rep_factory_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos
    !> MPI communicator
    integer, intent(in) :: comm

#ifdef PMC_USE_MPI
    integer :: aero_rep_type, i_aero_rep, prev_position

    prev_position = pos
    call pmc_mpi_unpack_integer(buffer, pos, aero_rep_type, comm)
    select case (aero_rep_type)
      case (AERO_REP_MODAL_BINNED_MASS)
        aero_rep => aero_rep_modal_binned_mass_t()
      case (AERO_REP_SINGLE_PARTICLE)
        aero_rep => aero_rep_single_particle_t()
      case default
        call die_msg(106634417, &
                "Trying to unpack aerosol representation of unknown type:"// &
                to_string(aero_rep_type))
    end select
    call aero_rep%bin_unpack(buffer, pos, comm)
    call assert(948795857, &
         pos - prev_position <= this%pack_size(aero_rep, comm))
#endif

  end function bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_factory
