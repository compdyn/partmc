! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_photolysis module.

!> \page phlex_rxn_photolysis Phlexible Mechanism for Chemistry: Photolysis
!!
!! Photolysis reactions take the form:
!!
!! \f[
!!   \mbox{X} + h\nu \to \mbox{Y_1} ( + \mbox{Y_2} \dots )
!! \f]
!!
!! where \f$\mbox{X}\f$ is the species being photolyzed, and
!! \f$\mbox{Y_n}\f$ are the photolysis products.
!!
!! Photolysis rate constants (including the \f$h\nu\f$ term) can be constant
!! or set from an external photolysis module using the
!! \c pmc_rxn_photolysis::rxn_photolysis_t::set_rate_const() function.
!! External modules can use the
!! \c pmc_rxn_photolysis::rxn_photolysis_t::get_property_set() function during
!! initilialization to access any needed reaction parameters.
!!
!! Input data for Photolysis equations should take the form :
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
!!   }
!! \endcode
!! The key-value pairs \b reactants, and \b products are required. There must
!! be exactly one key-value pair in the \b reactants object whose name is the
!! species being photolyzed and whose value is an empty \c json object. Any
!! number of products may be present. Products without a specified \b yield
!! are assumed to have a \b yield of 1.0. The \b "rate const" is optional and
!! can be used to set a rate constant (including the \f$h\nu\f$ term) that
!! remains constant throughout the model run. All other data is optional and
!! will be available to external photolysis modules during initialization.
!! Rate constants should be in units of \f$s^{-1}\f$.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_photolysis_t type and associated functions. 
module pmc_rxn_photolysis

  use pmc_constants,                        only: const
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, die_msg
  use pmc_rxn_data
  use pmc_chem_spec_data
  use pmc_property
  use pmc_phlex_state

  implicit none
  private

#define _NUM_REACT_ this%condensed_data_int(1)
#define _NUM_PROD_ this%condensed_data_int(2)
#define _PHOTO_ID_ this%condensed_data_int(3)
#define _BASE_RATE_ this%condensed_data_real(1)
#define _SCALING_ this%condensed_data_real(2)
#define _RATE_CONSTANT_ this%condensed_data_real(3)
#define _NUM_INT_PROP_ 3
#define _NUM_REAL_PROP_ 3
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _DERIV_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + _NUM_PROD_ + x)
#define _JAC_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + 2*(_NUM_REACT_+_NUM_PROD_) + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)

public :: rxn_photolysis_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_photolysis_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
    !> Set the photo id for this reaction
    procedure :: set_photo_id
    !> Get the reaction property set
    procedure :: get_property_set
  end type rxn_photolysis_t

  !> Constructor for rxn_photolysis_t
  interface rxn_photolysis_t
    procedure :: constructor
  end interface rxn_photolysis_t

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
  !! any required information from reactant, product and reaction 
  !! property_t objects. This routine should be called once for each reaction
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensures all data required during the model run are included
  !! in the condensed data arrays.
  subroutine initialize(this, chem_spec_data)
    
    !> Reaction data
    class(rxn_photolysis_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name
    integer(kind=i_kind) :: i_spec, i_qty

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(408416753, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(415586594, this%property_set%get_property_t(key_name, reactants), &
            "Photolysis reaction is missing reactants")
    key_name = "products"
    call assert_msg(180421290, this%property_set%get_property_t(key_name, products), &
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

    ! Allocate space in the condensed data arrays
    ! Space in this example is allocated for two sets of inidices for the 
    ! reactants and products, one molecular property for each reactant, 
    ! yields for the products and three reaction parameters.
    allocate(this%condensed_data_int(_NUM_INT_PROP_ + &
            (i_spec * 3) * (i_spec + products%size())))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_ + products%size()))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)
    
    ! Save the size of the reactant and product arrays (for reactions where these
    ! can vary)
    _NUM_REACT_ = i_spec
    _NUM_PROD_ = products%size()

    ! Get reaction parameters (it might be easiest to keep these at the beginning
    ! of the condensed data array, so they can be accessed using compliler flags)
    key_name = "rate const"
    if (.not. this%property_set%get_real(key_name, _BASE_RATE_)) then
      _BASE_RATE_ = real(0.0, kind=dp)
    end if
    key_name = "scaling factor"
    if (.not. this%property_set%get_real(key_name, _SCALING_)) then
      _SCALING_ = real(1.0, kind=dp)
    end if

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(747277322, _REACT_(i_spec).gt.0, &
              "Missing photolysis reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(231542303, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) then
        do i_qty = 1, temp_int - 1
          _REACT_(i_spec + i_qty) = _REACT_(i_spec)
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
      _PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(360988742, _PROD_(i_spec).gt.0, &
              "Missing photolysis product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(173703744, products%get_property_t(val=spec_props))
      key_name = "yield"
      if (spec_props%get_real(key_name, temp_real)) then
        _yield_(i_spec) = temp_real
      else
        _yield_(i_spec) = 1.0
      end if

      call products%iter_next()
      i_spec = i_spec + 1
    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rate constant expression
  function build_rate_const_expr(this, rxn_id) result (expr)

    !> Rate constant expression
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_photolysis_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id

    expr = ""

  end function build_rate_const_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build time derivative expression
  function build_deriv_expr(this, rxn_id, spec_id, chem_spec_data) &
                  result (expr)

    !> Contribution to time derivative expression for species spec_id
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_photolysis_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    ! Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set an id for this reaction that can be used by an external photolysis
  !! module to update the base (unscaled) rate constant during the model run.
  subroutine set_photo_id(this, photo_id)

    !> Reaction data 
    class(rxn_photolysis_t), intent(inout) :: this
    !> Photo id
    integer(kind=i_kind), intent(in) :: photo_id

    _PHOTO_ID_ = photo_id

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

#undef _NUM_REACT_
#undef _NUM_PROD_
#undef _PHOTO_ID_
#undef _BASE_RATE_
#undef _SCALING_
#undef _RATE_CONSTANT_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _REACT_
#undef _PROD_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _yield_
end module pmc_rxn_photolysis
