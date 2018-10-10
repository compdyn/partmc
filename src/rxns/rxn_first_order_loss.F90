! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_first_order_loss module.

!> \page phlex_rxn_first_order_loss Phlexible Module for Chemistry: First-Order Loss
!!
!! First-Order Loss reactions take the form:
!!
!! \f[\mbox{\ch{
!!   X ->
!! }}\f]
!!
!! where \f$\mbox{\ch{X}}\f$ is the species being lost.
!!
!! First-Order loss rate constants can be constant or set from an external
!! module using the
!! \c pmc_rxn_first_order_loss::rxn_update_data_first_order_loss_rate_t object.
!! External modules can use the
!! \c pmc_rxn_first_order_loss::rxn_first_order_loss_t::get_property_set()
!! function during initilialization to access any needed reaction parameters
!! to identify certain first-order loss reactions. The
!! \c pmc_rxn_first_order_loss::rxn_first_order_loss_t::set_rxn_id() function
!! can be used during initialization to set an integer id for a particular
!! reaction that can be used during solving to update the first-order loss
!! rate from an external module.
!!
!! Input data for first-order loss reactions have the following format :
!! \code{.json}
!!   {
!!     "type" : "FIRST_ORDER_LOSS",
!!     "species" : "species_name",
!!     "rate const" : 12.5,
!!     "scaling factor" : 1.2,
!!     ...
!!   }
!! \endcode
!! The key-value pairs \b species is required and its value must be the name
!! of the species being removed by the reaction. The \b rate \b const is
!! optional and can be used to set a rate constant that remains constant
!! throughout the model run. The \b scaling \b factor is also optional, and
!! can be used to set a constant scaling factor for the rate constant. When a
!! \b scaling \b factor is not provided, it is assumed to be 1.0. All other
!! data is optional and will be available to external modules during
!! initialization. Rate constants are in units of \f$s^{-1}\f$.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_first_order_loss_t type and associated functions.
module pmc_rxn_first_order_loss

  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_phlex_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, die_msg

  use iso_c_binding

  implicit none
  private

#define RXN_ID_ this%condensed_data_int(1)
#define REACT_ this%condensed_data_int(2)
#define DERIV_ID_ this%condensed_data_int(3)
#define JAC_ID_ this%condensed_data_int(4)
#define BASE_RATE_ this%condensed_data_real(1)
#define SCALING_ this%condensed_data_real(2)
#define RATE_CONSTANT_ this%condensed_data_real(3)
#define NUM_INT_PROP_ 4
#define NUM_REAL_PROP_ 3

public :: rxn_first_order_loss_t, rxn_update_data_first_order_loss_rate_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_first_order_loss_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Set the reaction id for this reaction
    procedure :: set_rxn_id
    !> Get the reaction property set
    procedure :: get_property_set
    !> Finalize the reaction
    final :: finalize
  end type rxn_first_order_loss_t

  !> Constructor for rxn_first_order_loss_t
  interface rxn_first_order_loss_t
    procedure :: constructor
  end interface rxn_first_order_loss_t

  !> First-Order Loss rate update object
  type, extends(rxn_update_data_t) :: rxn_update_data_first_order_loss_rate_t
  private
    logical :: is_malloced = .false.
  contains
    !> Initialize update data
    procedure :: initialize => update_data_rate_initialize
    !> Update the rate data
    procedure :: set_rate => update_data_rate_set
    !> Finalize the rate update data
    final :: update_data_rate_finalize
  end type rxn_update_data_first_order_loss_rate_t

  !> Interface to c reaction functions
  interface

    !> Allocate space for a rate update
    function rxn_first_order_loss_create_rate_update_data() &
              result (update_data) bind (c)
      use iso_c_binding
      !> Allocated update_data object
      type(c_ptr) :: update_data
    end function rxn_first_order_loss_create_rate_update_data

    !> Set a new first_order_loss rate
    subroutine rxn_first_order_loss_set_rate_update_data(update_data, &
              rxn_id, base_rate) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Reaction id from
      !! \c pmc_rxn_first_order_loss::rxn_first_order_loss_t::set_rxn_id
      integer(kind=c_int), value :: rxn_id
      !> New pre-scaling base first_order_loss rate
      real(kind=c_double), value :: base_rate
    end subroutine rxn_first_order_loss_set_rate_update_data

    !> Free an update rate data object
    pure subroutine rxn_free_update_data(update_data) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value, intent(in) :: update_data
    end subroutine rxn_free_update_data

  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for First-Order Loss reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_first_order_loss_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = GAS_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep)

    !> Reaction data
    class(rxn_first_order_loss_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props
    character(len=:), allocatable :: key_name, spec_name
    integer(kind=i_kind) :: i_spec, i_qty

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    call assert_msg(128411383, associated(this%property_set), &
            "Missing property set needed to initialize reaction")
    key_name = "species"
    call assert_msg(164644065, &
            this%property_set%get_string(key_name, spec_name), &
            "First-Order Loss reaction is missing species name")

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_))
    allocate(this%condensed_data_real(NUM_REAL_PROP_))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Get reaction parameters
    key_name = "rate const"
    if (.not. this%property_set%get_real(key_name, BASE_RATE_)) then
      BASE_RATE_ = real(0.0, kind=dp)
    end if
    key_name = "scaling factor"
    if (.not. this%property_set%get_real(key_name, SCALING_)) then
      SCALING_ = real(1.0, kind=dp)
    end if

    ! Save the index of this species in the state variable array
    REACT_ = chem_spec_data%gas_state_id(spec_name)

    ! Make sure the species exists
    call assert_msg(331442196, REACT_.gt.0, &
            "Missing first-order loss species: "//spec_name)


  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set an id for this reaction that can be used by an external
  !! module to update the base (unscaled) rate constant during the model run.
  subroutine set_rxn_id(this, rxn_id)

    !> Reaction data
    class(rxn_first_order_loss_t), intent(inout) :: this
    !> Reaction id
    integer(kind=i_kind), intent(in) :: rxn_id

    RXN_ID_ = rxn_id

  end subroutine set_rxn_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the reaction properties. (For use by external modules.)
  function get_property_set(this) result(prop_set)

    !> Reaction properties
    type(property_t), pointer :: prop_set
    !> Reaction data
    class(rxn_first_order_loss_t), intent(in) :: this

    prop_set => this%property_set

  end function get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_first_order_loss_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for first_order_loss rate constants
  subroutine update_data_rate_set(this, rxn_id, base_rate)

    !> Update data
    class(rxn_update_data_first_order_loss_rate_t), intent(inout) :: this
    !> Reaction id from
    !! \c pmc_rxn_first_order_loss::rxn_first_order_loss_t::set_rxn_id
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Updated pre-scaling first_order_loss rate
    real(kind=dp), intent(in) :: base_rate

    call rxn_first_order_loss_set_rate_update_data(this%get_data(), rxn_id, &
            base_rate)

  end subroutine update_data_rate_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize update data
  subroutine update_data_rate_initialize(this, rxn_type)

    !> Update data object
    class(rxn_update_data_first_order_loss_rate_t) :: this
    !> Reaction type id
    integer(kind=i_kind), intent(in) :: rxn_type

    this%rxn_type = int(rxn_type, kind=c_int)
    this%update_data = rxn_first_order_loss_create_rate_update_data()
    this%is_malloced = .true.

  end subroutine update_data_rate_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize an update data object
  elemental subroutine update_data_rate_finalize(this)

    !> Update data object to free
    type(rxn_update_data_first_order_loss_rate_t), intent(inout) :: this

    if (this%is_malloced) call rxn_free_update_data(this%update_data)

  end subroutine update_data_rate_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef RXN_ID_
#undef REACT_
#undef DERIV_ID_
#undef JAC_ID_
#undef BASE_RATE_
#undef SCALING_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
end module pmc_rxn_first_order_loss
