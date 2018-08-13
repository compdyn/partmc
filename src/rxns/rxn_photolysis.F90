! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_photolysis module.

!> \page phlex_rxn_photolysis Phlexible Module for Chemistry: Photolysis
!!
!! Photolysis reactions take the form:
!!
!! \f[\mbox{\ch{
!!   X + h $\nu$ -> Y_1 ( + Y_2 \dots )
!! }}\f]
!!
!! where \f$\mbox{\ch{X}}\f$ is the species being photolyzed, and
!! \f$\mbox{\ch{Y_n}}\f$ are the photolysis products.
!!
!! Photolysis rate constants (including the \f$\mbox{\ch{h $\nu$}}\f$ term)
!! can be constant or set from an external photolysis module using the
!! \c pmc_rxn_photolysis::rxn_update_data_photolysis_rate_t object.
!! External modules can use the
!! \c pmc_rxn_photolysis::rxn_photolysis_t::get_property_set() function during
!! initilialization to access any needed reaction parameters to identify
!! certain photolysis reactions. The 
!! \c pmc_rxn_photolysis::rxn_photolysis_t::set_photo_id() function can be
!! used during initialization to set an integer id for a particular reaction
!! that can be used during solving to update the photolysis rate from an
!! external module.
!!
!! Input data for photolysis reactions have the following format :
!! \code{.json}
!!   {
!!     "type" : "PHOTOLYSIS",
!!     "reactants" : {
!!       "spec1" : {}
!!     },
!!     "products" : {
!!       "spec2" : {},
!!       "spec3" : { "yield" : 0.65 },
!!       ...
!!     },
!!     "rate const" : 12.5,
!!     "scaling factor" : 1.2,
!!     ...
!!   }
!! \endcode
!! The key-value pairs \b reactants, and \b products are required. There must
!! be exactly one key-value pair in the \b reactants object whose name is the
!! species being photolyzed and whose value is an empty \c json object. Any
!! number of products may be present. Products without a specified \b yield
!! are assumed to have a \b yield of 1.0. The \b rate \b const is optional and
!! can be used to set a rate constant (including the \f$h\nu\f$ term) that
!! remains constant throughout the model run. The \b scaling \b factor is also
!! optional, and can be used to set a constant scaling factor for the rate
!! constant. When the \b scaling \b factor is not provided, it is assumed to
!! be 1.0. All other data is optional and will be available to external
!! photolysis modules during initialization. Rate constants are in units of 
!! \f$s^{-1}\f$.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_photolysis_t type and associated functions. 
module pmc_rxn_photolysis

  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_phlex_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: phlex_real, phlex_int, &
                                                  to_string, assert, &
                                                  assert_msg, die_msg, &
                                                  align_ratio

  use iso_c_binding

  implicit none
  private

#define NUM_REACT_ this%condensed_data_int(1)
#define NUM_PROD_ this%condensed_data_int(2)
#define PHOTO_ID_ this%condensed_data_int(3)
#define INT_DATA_SIZE_ this%condensed_data_int(4)
#define FLOAT_DATA_SIZE_ this%condensed_data_int(5)
#define BASE_RATE_ this%condensed_data_real(1)
#define SCALING_ this%condensed_data_real(2)
#define RATE_CONSTANT_ this%condensed_data_real(3)
#define NUM_INT_PROP_ 5
#define NUM_REAL_PROP_ 3
#define REACT_(x) this%condensed_data_int(NUM_INT_PROP_ + x)
#define PROD_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x)
#define YIELD_(x) this%condensed_data_real(NUM_REAL_PROP_ + x)

public :: rxn_photolysis_t, rxn_update_data_photolysis_rate_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_photolysis_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Set the photo id for this reaction
    procedure :: set_photo_id
    !> Get the reaction property set
    procedure :: get_property_set
    !> Finalize the reaction
    final :: finalize
  end type rxn_photolysis_t

  !> Constructor for rxn_photolysis_t
  interface rxn_photolysis_t
    procedure :: constructor
  end interface rxn_photolysis_t

  !> Photolysis rate update object
  type, extends(rxn_update_data_t) :: rxn_update_data_photolysis_rate_t
  private
    logical :: is_malloced = .false.
  contains
    !> Initialize update data
    procedure :: initialize => update_data_rate_initialize
    !> Update the rate data
    procedure :: set_rate => update_data_rate_set
    !> Finalize the rate update data
    final :: update_data_rate_finalize
  end type rxn_update_data_photolysis_rate_t

  !> Interface to c reaction functions
  interface

    !> Allocate space for a rate update
    function rxn_photolysis_create_rate_update_data() result (update_data) &
              bind (c)
      use iso_c_binding
      !> Allocated update_data object
      type(c_ptr) :: update_data
    end function rxn_photolysis_create_rate_update_data

    !> Set a new photolysis rate
    subroutine rxn_photolysis_set_rate_update_data(update_data, photo_id, &
              base_rate) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Photo id from pmc_rxn_photolysis::rxn_photolysis_t::set_photo_id
      integer(kind=c_int), value :: photo_id
      !> New pre-scaling base photolysis rate
      real(kind=PMC_F90_C_FLOAT), value :: base_rate
    end subroutine rxn_photolysis_set_rate_update_data
  
    !> Free an update rate data object
    pure subroutine rxn_free_update_data(update_data) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value, intent(in) :: update_data
    end subroutine rxn_free_update_data

  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Photolysis reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_photolysis_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = GAS_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep)
    
    !> Reaction data
    class(rxn_photolysis_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name
    integer(kind=phlex_int) :: i_spec, i_qty
    integer(kind=phlex_int) :: int_data_size, float_data_size

    integer(kind=phlex_int) :: temp_int
    real(kind=phlex_real) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(408416753, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(415586594, &
            this%property_set%get_property_t(key_name, reactants), &
            "Photolysis reaction is missing reactants")
    key_name = "products"
    call assert_msg(180421290, &
            this%property_set%get_property_t(key_name, products), &
            "Photolysis reaction is missing products")

    ! Count the number of reactants (including those with a qty specified)
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(240165383, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) i_spec = i_spec+temp_int-1
      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Calculate int and float array sizes with alignment spacing
    int_data_size = NUM_INT_PROP_ + &
            (i_spec + 2) * (i_spec + products%size())
    int_data_size = int_data_size + mod(int_data_size, align_ratio)
    float_data_size = NUM_REAL_PROP_ + products%size()

    ! Allocate space in the condensed data arrays
    ! Space in this example is allocated for two sets of inidices for the 
    ! reactants and products, one molecular property for each reactant, 
    ! yields for the products and three reaction parameters.
    allocate(this%condensed_data_int(int_data_size))
    allocate(this%condensed_data_real(float_data_size))
    this%condensed_data_int(:) = int(0, kind=phlex_int)
    this%condensed_data_real(:) = real(0.0, kind=phlex_real)
    INT_DATA_SIZE_ = int_data_size
    FLOAT_DATA_SIZE_ = float_data_size
    
    ! Save the size of the reactant and product arrays (for reactions where
    ! these can vary)
    NUM_REACT_ = i_spec
    NUM_PROD_ = products%size()

    ! Get reaction parameters (it might be easiest to keep these at the
    ! beginning of the condensed data array, so they can be accessed using 
    ! compliler flags)
    key_name = "rate const"
    if (.not. this%property_set%get_real(key_name, BASE_RATE_)) then
      BASE_RATE_ = real(0.0, kind=phlex_real)
    end if
    key_name = "scaling factor"
    if (.not. this%property_set%get_real(key_name, SCALING_)) then
      SCALING_ = real(1.0, kind=phlex_real)
    end if

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(747277322, REACT_(i_spec).gt.0, &
              "Missing photolysis reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(231542303, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) then
        do i_qty = 1, temp_int - 1
          REACT_(i_spec + i_qty) = REACT_(i_spec)
        end do
        i_spec = i_spec + temp_int - 1
      end if

      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Make sure exactly one reactant is present
    call assert_msg(908486656, i_spec.eq.2, "Incorrect number of reactants"//&
            " for Photolysis reaction: "//to_string(i_spec-1))

    ! Get the indices and chemical properties for the products
    call products%iter_reset()
    i_spec = 1
    do while (products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(360988742, PROD_(i_spec).gt.0, &
              "Missing photolysis product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(173703744, products%get_property_t(val=spec_props))
      key_name = "yield"
      if (spec_props%get_real(key_name, temp_real)) then
        YIELD_(i_spec) = temp_real
      else
        YIELD_(i_spec) = 1.0
      end if

      call products%iter_next()
      i_spec = i_spec + 1
    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set an id for this reaction that can be used by an external photolysis
  !! module to update the base (unscaled) rate constant during the model run.
  subroutine set_photo_id(this, photo_id)

    !> Reaction data 
    class(rxn_photolysis_t), intent(inout) :: this
    !> Photo id
    integer(kind=phlex_int), intent(in) :: photo_id

    PHOTO_ID_ = photo_id

  end subroutine set_photo_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the reaction properties. (For use by external photolysis modules.)
  function get_property_set(this) result(prop_set)

    !> Reaction properties
    type(property_t), pointer :: prop_set
    !> Reaction data
    class(rxn_photolysis_t), intent(in) :: this

    prop_set => this%property_set

  end function get_property_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_photolysis_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for photolysis rate constants
  subroutine update_data_rate_set(this, photo_id, base_rate)

    !> Update data
    class(rxn_update_data_photolysis_rate_t), intent(inout) :: this
    !> Photo id from pmc_rxn_photolysis::rxn_photolysis_t::set_photo_id
    integer(kind=phlex_int), intent(in) :: photo_id
    !> Updated pre-scaling photolysis rate
    real(kind=phlex_real), intent(in) :: base_rate

    call rxn_photolysis_set_rate_update_data(this%get_data(), photo_id, &
            base_rate)

  end subroutine update_data_rate_set

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize update data
  subroutine update_data_rate_initialize(this, rxn_type)

    !> Update data object
    class(rxn_update_data_photolysis_rate_t) :: this
    !> Reaction type id
    integer(kind=phlex_int), intent(in) :: rxn_type

    this%rxn_type = int(rxn_type, kind=c_int)
    this%update_data = rxn_photolysis_create_rate_update_data()
    this%is_malloced = .true.

  end subroutine update_data_rate_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize an update data object
  elemental subroutine update_data_rate_finalize(this)

    !> Update data object to free
    type(rxn_update_data_photolysis_rate_t), intent(inout) :: this

    if (this%is_malloced) call rxn_free_update_data(this%update_data)

  end subroutine update_data_rate_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef NUM_REACT_
#undef NUM_PROD_
#undef PHOTO_ID_
#undef BASE_RATE_
#undef SCALING_
#undef RATE_CONSTANT_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef REACT_
#undef PROD_
#undef DERIV_ID_
#undef JAC_ID_
#undef YIELD_
end module pmc_rxn_photolysis
