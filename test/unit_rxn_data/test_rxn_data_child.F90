! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_rxn_data module.

!> The rxn_test_t type and associated functions, used to test the
!! abstract rxn_data_t type. This module can be used as a template for
!! developing new reaction modules.
module pmc_test_rxn_data_child

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
#define _A_ this%condensed_data_real(1)
#define _B_ this%condensed_data_real(2)
#define _theta_ this%condensed_data_real(3)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 3
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _JACID_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + _NUM_PROD_ + x)
#define _MW_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + _NUM_REACT_ + x)

public :: rxn_test_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_test_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
  end type rxn_test_t

  !> Constructor for rxn_test_t
  interface rxn_test_t
    procedure :: constructor
  end interface rxn_test_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for test reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_test_t), pointer :: new_obj

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
    class(rxn_test_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name
    integer(kind=i_kind) :: i_spec, i_qty, spec_phase

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(751697804, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(278982963, this%property_set%get_property_t(key_name, reactants), &
            "Generic reaction is missing reactants")
    key_name = "products"
    call assert_msg(197632993, this%property_set%get_property_t(key_name, products), &
            "Generic reaction is missing products")

    ! Count the number of reactants (including those with a qty specified)
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(472956832, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) i_spec = i_spec+temp_int-1
      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Allocate space in the condensed data arrays
    ! Space in this example is allocated for two sets of inidices for the 
    ! reactants and products, one molecular property for each reactant, 
    ! yields for the products and three reaction parameters.
    allocate(this%condensed_data_int(_NUM_INT_PROP_+i_spec+products%size() &
            + i_spec * (i_spec + products%size())))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_+i_spec+products%size()))
    
    ! Save the size of the reactant and product arrays (for reactions where these
    ! can vary)
    _NUM_REACT_ = i_spec
    _NUM_PROD_ = products%size()

    ! Check the number of reactants
    call assert_msg(212395756, _NUM_REACT_.eq.2, "Expected 2 reactants, got "// &
            to_string(_NUM_REACT_))

    ! Get reaction parameters (it might be easiest to keep these at the beginning
    ! of the condensed data array, so they can be accessed using compliler flags)
    key_name = "A"
    call assert_msg(785292700, this%property_set%get_real(key_name, _A_), &
            "Missing reaction parameter A")
    key_name = "B"
    call assert_msg(132874546, this%property_set%get_real(key_name, _B_), &
            "Missing reaction parameter B")
    key_name = "theta"
    call assert_msg(127610239, this%property_set%get_real(key_name, _theta_), &
            "Missing reaction parameter theta")

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Check the species type
      call assert(545500490, chem_spec_data%get_phase(spec_name, spec_phase))
      call assert_msg(193108104, spec_phase.eq.CHEM_SPEC_GAS_PHASE, &
              "Expected gas-phase species for "//spec_name)

      ! Save the index of this species in the chem_spec_state_t variable
      _REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Get a chemical property for this species and save it in the condensed real array
      call assert(639119436, chem_spec_data%get_property_set(spec_name, spec_props))
      key_name = "MW"
      call assert_msg(542784410, spec_props%get_real(key_name, _MW_(i_spec)), &
              "Missing reaction parameter MW")

      ! Get properties included with this reactant in the reaction data
      call assert(193983497, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) then
        do i_qty = 1, temp_int - 1
          _REACT_(i_spec + i_qty) = _REACT_(i_spec)
          _MW_(i_spec + i_qty) = _MW_(i_spec)
        end do
        i_spec = i_spec + temp_int + 1
      end if

      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Get the indices and chemical properties for the products
    call products%iter_reset()
    i_spec = 1
    do while (products%get_key(spec_name))

      ! Check the species type
      call assert(293069092, chem_spec_data%get_phase(spec_name, spec_phase))
      call assert_msg(310787481, spec_phase.eq.CHEM_SPEC_GAS_PHASE, &
              "Expected gas-phase species for "//spec_name)

      ! Save the index of this species in the chem_spec_state_t variable
      _PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Get properties included with this product in the reaction data
      call assert(425790028, products%get_property_t(val=spec_props))
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
    class(rxn_test_t), intent(in) :: this
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
    class(rxn_test_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    ! Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Optional support function for calculating reactions rate constants
  !!
  !! Sample rate constant equation:
  !!     A*exp(B+temp/theta)
  !!
  real(kind=dp) function rate_const(this, phlex_state)

    !> Reaction data
    class(rxn_test_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state

    integer(kind=i_kind) :: i_spec

    rate_const = _A_*exp(_B_ + &
            phlex_state%env_state%temp/_theta_)

    do i_spec = 1, _NUM_REACT_
      rate_const = rate_const/_MW_(i_spec) 
    end do

  end function rate_const

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_test_rxn_data_child
