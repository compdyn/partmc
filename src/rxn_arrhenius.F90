! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_arrhenius module.

!> The rxn_arrhenius_t type and associated functions. The Arrhenius reaction
!! rate constant equation is:
!!
!!   A*exp(-Ea/(kb*temp))
!!
!! where A is the pre-exponential factor ([concentration unit]^-(n-1) s^-1),
!! n is the number of reactants, Ea is the  activation energy (J), kb is the 
!! Boltzmann constant (J/K) and temp is the temperature (K).
!!
#ifdef PMC_USE_JSON
!! Input data for Arrhenius equations should take the form :
!!
!!   {
!!     "type" : "RXN_ARRHENIUS",
!!     "A" : 123.45,
!!     "Ea" : 123.45,
!!     "reactants" : {
!!       "spec1" : {},
!!       "spec2" : { "qty" : 2 },
!!       ...
!!     },
!!     "products" : {
!!       "spec3" : {},
!!       "spec4" : { "yield" : 0.65 },
!!       ...
!!     }
!!   }
!!
!! The key-value pairs A, Ea, reactants, and products are required. Reactants
!! without a qty value are assumed to appear once in the reaction equation.
!! Products without a specified yield are assumed to have a yield of 1.0.
#endif
module pmc_rxn_arrhenius

  use pmc_constants,                        only: const_t
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
#define _Ea_ this%condensed_data_real(2)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 2
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)

public :: rxn_arrhenius_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_arrhenius_t
  contains
    !> Reaction initialization
    procedure :: initialize => pmc_rxn_arrhenius_initialize
    !> Time derivative contribution
    procedure :: func_contrib => pmc_rxn_arrhenius_func_contrib
    !> Jacobian matrix contribution
    procedure :: jac_contrib => pmc_rxn_arrhenius_jac_contrib
    !> Calculate the rate constant
    procedure, private :: rate_const => pmc_rxn_arrhenius_rate_const
  end type rxn_arrhenius_t

  !> Constructor for rxn_arrhenius_t
  interface rxn_arrhenius_t
    procedure :: pmc_rxn_arrhenius_constructor
  end interface rxn_arrhenius_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Arrhenius reaction
  function pmc_rxn_arrhenius_constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_arrhenius_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = GAS_RXN

  end function pmc_rxn_arrhenius_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information from reactant, product and reaction 
  !! property_t objects. This routine should be called once for each reaction
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensures all data required during the model run are included
  !! in the condensed data arrays.
  subroutine pmc_rxn_arrhenius_initialize(this, chem_spec_data)
    
    !> Reaction data
    class(rxn_arrhenius_t), intent(inout) :: this
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
            "Arrhenius reaction is missing reactants")
    key_name = "products"
    call assert_msg(304540307, this%property_set%get_property_t(key_name, products), &
            "Arrhenius reaction is missing products")

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
    key_name = "A"
    call assert_msg(189315484, this%property_set%get_real(key_name, _A_), &
            "Missing Arrhenius reaction parameter A")
    key_name = "Ea"
    call assert_msg(861320020, this%property_set%get_real(key_name, _Ea_), &
            "Missing Arrhenius reaction parameter Ea")

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _REACT_(i_spec) = chem_spec_data%state_id(spec_name)

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

    ! Get the indices and chemical properties for the products
    call products%iter_reset()
    i_spec = 1
    do while (products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _PROD_(i_spec) = chem_spec_data%state_id(spec_name)

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

  end subroutine pmc_rxn_arrhenius_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution to the time derivative vector. The current 
  !! model state is provided for species concentrations, aerosol state and 
  !! environmental variables. All other parameters must have been saved to the
  !! reaction data instance during initialization.
  !!
  !! Sample reaction rate equation:
  !!     [reactant1] * [reactant2] *A*exp(B+temp/theta)
  !!
  subroutine  pmc_rxn_arrhenius_func_contrib(this, phlex_state, func)

    !> Reaction data
    class(rxn_arrhenius_t), intent(in) :: this
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

  end subroutine pmc_rxn_arrhenius_func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution of this reaction to the Jacobian matrix.
  !! The current model state is provided for species concentrations and 
  !! aerosol state. All other parameters must have been saved to the reaction 
  !! data instance during initialization.
  subroutine pmc_rxn_arrhenius_jac_contrib(this, phlex_state, jac_matrix)

    !> Reaction data
    class(rxn_arrhenius_t), intent(in) :: this
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

  end subroutine pmc_rxn_arrhenius_jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the reaction rate constant
  real(kind=dp) function pmc_rxn_arrhenius_rate_const(this, phlex_state) &
                  result(rate_const)

    !> Reaction data
    class(rxn_arrhenius_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state

    integer(kind=i_kind) :: i_spec
    type(const_t) :: const

    rate_const = _A_*exp(- _Ea_ / ( &
            phlex_state%env_state%temp * const%boltzmann))

  end function pmc_rxn_arrhenius_rate_const

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_arrhenius
