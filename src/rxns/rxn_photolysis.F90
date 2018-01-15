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
!!     "rate const" : 12.5
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
#define _RATE_CONST_ this%condensed_data_real(1)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 1
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)

public :: rxn_photolysis_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_photolysis_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Time derivative contribution
    procedure :: func_contrib
    !> Jacobian matrix contribution
    procedure :: jac_contrib
    !> Calculate the rate constant
    procedure, private :: rate_const
    !> Set the photolysis rate constant
    procedure :: set_rate_const
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
    integer(kind=i_kind) :: i_spec

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(255324828, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(250060521, this%property_set%get_property_t(key_name, reactants), &
            "Photolysis reaction is missing reactants")
    key_name = "products"
    call assert_msg(304540307, this%property_set%get_property_t(key_name, products), &
            "Photolysis reaction is missing products")

    ! Count the number of reactants (including those with a qty specified)
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(243342975, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) i_spec = i_spec+temp_int-1
      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Allocate space in the condensed data arrays
    ! Space in this example is allocated for two sets of inidices for the 
    ! reactants and products, one molecular property for each reactant, 
    ! yields for the products and three reaction parameters.
    allocate(this%condensed_data_int(_NUM_INT_PROP_+i_spec+products%size()))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_+products%size()))
    
    ! Save the size of the reactant and product arrays (for reactions where these
    ! can vary)
    _NUM_REACT_ = i_spec
    _NUM_PROD_ = products%size()

    ! Get reaction parameters (it might be easiest to keep these at the beginning
    ! of the condensed data array, so they can be accessed using compliler flags)
    key_name = "rate const"
    if (.not. this%property_set%get_real(key_name, _RATE_CONST_)) then
      _RATE_CONST_ = real(0.0, kind=dp)
    end if

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(929298013, _REACT_(i_spec).gt.0, &
              "Missing photolysis reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(796763915, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) then
        do i_spec = i_spec + 1, i_spec + temp_int - 1
          _REACT_(i_spec) = _REACT_(i_spec-1)
        end do
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
      call assert(451185800, products%get_property_t(val=spec_props))
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

  !> Calculate the contribution to the time derivative vector. The current 
  !! model state is provided for species concentrations, aerosol state and 
  !! environmental variables. All other parameters must have been saved to the
  !! reaction data instance during initialization.
  !!
  !! Sample reaction rate equation:
  !!     [reactant1] * [reactant2] *A*exp(B+temp/theta)
  !!
  subroutine  func_contrib(this, phlex_state, func)

    !> Reaction data
    class(rxn_photolysis_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Time derivative vector. This vector may include contributions from
    !! other reactions, so the contributions from this reaction should
    !! append, not overwrite, the values already in the vector
    real(kind=dp), pointer, intent(inout) :: func(:)
   
    integer(kind=i_kind) :: i_spec
    real(kind=dp) :: rate

    rate = this%rate_const(phlex_state)
    do i_spec=1, _NUM_REACT_
      rate = rate * phlex_state%state_var(_REACT_(i_spec))
    end do
    
    if (rate.eq.real(0.0, kind=dp)) return

    do i_spec=1, _NUM_REACT_
      func(_REACT_(i_spec)) = &
        func(_REACT_(i_spec)) &
        - rate
    end do

    do i_spec=1, _NUM_PROD_
      func(_PROD_(i_spec)) = &
        func(_PROD_(i_spec)) &
        + _yield_(i_spec) * rate
    end do

  end subroutine func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution of this reaction to the Jacobian matrix.
  !! The current model state is provided for species concentrations and 
  !! aerosol state. All other parameters must have been saved to the reaction 
  !! data instance during initialization.
  subroutine jac_contrib(this, phlex_state, jac_matrix)

    !> Reaction data
    class(rxn_photolysis_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Jacobian matrix. This matrix may include contributions from other
    !! reactions, so the contributions from this reaction should append,
    !! not overwrite, the values already in the matrix.
    real(kind=dp), pointer, intent(inout) :: jac_matrix(:,:)

    character(len=16) :: out
    integer(kind=i_kind) :: i_spec_i, i_spec_d
    real(kind=dp) :: rate, rate_const

    rate_const = this%rate_const(phlex_state)
    do i_spec_i=1, _NUM_REACT_
      rate_const = rate_const * &
              phlex_state%state_var(_REACT_(i_spec_i))
    end do

    if (rate_const.eq.real(0.0, kind=dp)) return

    do i_spec_d=1, _NUM_REACT_
      do i_spec_i=1, _NUM_REACT_
        rate = rate_const
        rate = rate / phlex_state%state_var( &
                _REACT_(i_spec_i))
        jac_matrix(_REACT_(i_spec_d), &
                _REACT_(i_spec_i)) = &
              jac_matrix(_REACT_(i_spec_d) &
              , _REACT_(i_spec_i)) - &
              rate
      end do
    end do

    do i_spec_d=1, _NUM_PROD_
      do i_spec_i=1, _NUM_REACT_
        rate = rate_const
        rate = rate / phlex_state%state_var( &
                _REACT_(i_spec_i) )
        jac_matrix(_PROD_(i_spec_d), &
                _REACT_(i_spec_i) ) = &
              jac_matrix(_PROD_(i_spec_d) &
              , _REACT_(i_spec_i) ) + &
              _yield_(i_spec_d) * rate
      end do
    end do

  end subroutine jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the reaction rate constant
  real(kind=dp) function rate_const(this, phlex_state)

    !> Reaction data
    class(rxn_photolysis_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state

    rate_const = _RATE_CONST_

  end function rate_const

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the rate constant from an external photolysis module
  subroutine set_rate_const(this, rate_const)

    !> Reaction data
    class(rxn_photolysis_t), intent(inout) :: this
    !> Rate constant \f$k_{photo}*h\nu\f$ (\f$s^{-1}\f$)
    real(kind=dp), intent(in) ::rate_const

    _RATE_CONST_ = rate_const

  end subroutine set_rate_const

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

end module pmc_rxn_photolysis
