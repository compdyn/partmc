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
  use pmc_integration_data

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
    !> Indicate which Jacobian elements are populated by this reaction
    procedure :: get_used_jac_elem
    !> Update indexes used during integration
    procedure :: update_integration_ids
    !> Time derivative contribution
    procedure :: func_contrib
    !> Jacobian matrix contribution
    procedure :: jac_contrib
    !> Get test info
    procedure :: get_test_info
    !> Calculate the rate constant
    procedure, private :: rate_const
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
    integer(kind=i_kind) :: i_spec, i_qty

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
      call assert_msg(193108104, chem_spec_data%get_type(spec_name).eq.GAS_SPEC, &
              "Expected gas-phase species for "//spec_name)

      ! Save the index of this species in the chem_spec_state_t variable
      _REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Get a chemical property for this species and save it in the condensed real array
      spec_props => chem_spec_data%get_property_set(spec_name)
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
      call assert_msg(293069092, chem_spec_data%get_type(spec_name).eq.GAS_SPEC, &
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

  !> Indicate which Jacobian elements are populated by this reaction
  subroutine get_used_jac_elem(this, use_jac_elem)

    !> Reaction data
    class(rxn_test_t), intent(in) :: this
    !> Matrix of flags indicating whether to include Jacobian elements
    !! in the sparse matrix
    integer(kind=i_kind), allocatable, intent(inout) :: use_jac_elem(:,:)

    integer(kind=i_kind) :: i_spec_d, i_spec_i

    do i_spec_d = 1, _NUM_REACT_
      do i_spec_i = 1, _NUM_REACT_
        use_jac_elem(_REACT_(i_spec_d), &
            _REACT_(i_spec_i)) = 1
      end do
    end do

    do i_spec_d = 1, _NUM_PROD_
      do i_spec_i = 1, _NUM_REACT_
        use_jac_elem(_PROD_(i_spec_d), &
            _REACT_(i_spec_i)) = 1
      end do
    end do

  end subroutine get_used_jac_elem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update indexes used during integration
  subroutine update_integration_ids(this, integration_data)

    !> Reation data
    class(rxn_test_t), intent(inout) :: this
    !> Integration data
    class(integration_data_t), intent(in) :: integration_data

    integer(kind=i_kind) :: i_spec_d, i_spec_i, jac_id

    jac_id = 0
    do i_spec_d = 1, _NUM_REACT_
      do i_spec_i = 1, _NUM_REACT_
        jac_id = jac_id + 1  
        _JACID_(jac_id) = &
            integration_data%get_jac_elem_id( &
            _REACT_(i_spec_d), &
            _REACT_(i_spec_i))
      end do
    end do

    do i_spec_d = 1, _NUM_PROD_
      do i_spec_i = 1, _NUM_REACT_
        jac_id = jac_id + 1  
        _JACID_(jac_id) = &
            integration_data%get_jac_elem_id( &
            _PROD_(i_spec_d), &
            _REACT_(i_spec_i))
      end do
    end do

  end subroutine update_integration_ids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution to the time derivative vector. The current 
  !! model state is provided for species concentrations, aerosol state and 
  !! environmental variables. All other parameters must have been saved to the
  !! reaction data instance during initialization.
  !!
  !! Sample reaction rate equation:
  !!     [reactant1]/MW1 * [reactant2]/MW2 *A*exp(B+temp/theta)
  !!
  subroutine  func_contrib(this, integration_data)

    !> Reaction data
    class(rxn_test_t), intent(in) :: this
    !> Integration data
    class(integration_data_t), intent(inout) :: integration_data

    integer(kind=i_kind) :: i_spec
    real(kind=dp) :: rate

    rate = this%rate_const(integration_data%phlex_state)
    do i_spec=1, _NUM_REACT_
      rate = rate * integration_data%phlex_state%state_var(_REACT_(i_spec))
    end do
    
    if (rate.eq.real(0.0, kind=dp)) return

    do i_spec=1, _NUM_REACT_
      call integration_data%add_deriv_contrib( _REACT_(i_spec), &
        -rate * _MW_(i_spec))
    end do

    do i_spec=1, _NUM_PROD_
      call integration_data%add_deriv_contrib( _PROD_(i_spec), &
        _yield_(i_spec) * rate)
    end do

  end subroutine func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution of this reaction to the Jacobian matrix.
  !! The current model state is provided for species concentrations and 
  !! aerosol state. All other parameters must have been saved to the reaction 
  !! data instance during initialization.
  subroutine jac_contrib(this, integration_data)

    !> Reaction data
    class(rxn_test_t), intent(in) :: this
    !> Integration data
    class(integration_data_t), intent(inout) :: integration_data

    integer(kind=i_kind) :: i_spec_i, i_spec_d, jac_id
    real(kind=dp) :: rate, rate_const

    rate_const = this%rate_const(integration_data%phlex_state)
    do i_spec_i=1, _NUM_REACT_
      rate_const = rate_const * &
              integration_data%phlex_state%state_var(_REACT_(i_spec_i))
    end do

    if (rate_const.eq.real(0.0, kind=dp)) return

    jac_id = 0
    do i_spec_d=1, _NUM_REACT_
      do i_spec_i=1, _NUM_REACT_
        jac_id = jac_id + 1
        rate = rate_const
        rate = rate / integration_data%phlex_state%state_var( &
                _REACT_(i_spec_i))
        call integration_data%add_jac_contrib( &
                _JACID_(jac_id), &
                -rate*_MW_(i_spec_d)) 
      end do
    end do

    do i_spec_d=1, _NUM_PROD_
      do i_spec_i=1, _NUM_REACT_
        jac_id = jac_id + 1
        rate = rate_const
        rate = rate / integration_data%phlex_state%state_var( &
                _REACT_(i_spec_i))
        call integration_data%add_jac_contrib( &
              _JACID_(jac_id), &
              _yield_(i_spec_d) * rate)
      end do
    end do

  end subroutine jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the reaction rate and rate constant for a given model state. The
  !! definition of the rate and rate constant depends on the extending type.
  !! THIS FUNCTION IS ONLY FOR TESTING.
  subroutine get_test_info(this, phlex_state, rate, rate_const, property_set)

    !> Reaction data
    class(rxn_test_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Reaction rate (definition depends on extending type)
    real(kind=dp), intent(out) :: rate
    !> Rate constant (definition depends on extending type)
    real(kind=dp), intent(out) :: rate_const
    !> Reaction properties
    type(property_t), pointer, intent(out) :: property_set

    integer(kind=i_kind) :: i_spec

    rate_const = this%rate_const(phlex_state)
    rate = rate_const
    do i_spec=1, _NUM_REACT_
      rate = rate * phlex_state%state_var(_REACT_(i_spec))
    end do
    property_set => this%property_set

  end subroutine get_test_info

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
