! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_arrhenius module.

!> \page phlex_rxn_arrhenius Phlexible Mechanism for Chemistry: Arrhenius Reaction
!!
!! Arrhenius-like reaction rate constant equations can take the form:
!!
!! \f[
!!   Ae^{(\frac{-E_a}{k_bT})}(\frac{T}{D})^B(1.0+E*P)
!! \f]
!!
!! where \f$A\f$ is the pre-exponential factor 
!! (\f$[\mbox{concentration unit}]^{-(n-1)} s^{-1}\f$),
!! \f$n\f$ is the number of reactants, \f$E_a\f$ is the activation energy
!! (J), \f$k_b\f$ is the Boltzmann constant (J/K), \f$D\f$ (K), \f$B\f$ 
!! (unitless) and \f$E\f$ (atm) are reaction parameters, \f$T\f$ is the
!! temperature (K), and \f$P\f$ is the pressure (atm). The first two terms
!! are described in Finlayson-Pitts and Pitts (2000). The final term is
!! included to accomodate CMAQ EBI solver type 7 rate constants.
!!
!! Input data for Arrhenius equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "ARRHENIUS",
!!     "A" : 123.45,
!!     "Ea" : 123.45,
!!     "B"  : 1.3,
!!     "D"  : 300.0,
!!     "E"  : 0.6,
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
!! \endcode
!! The key-value pairs \b reactants, and \b products are required. Reactants
!! without a \b qty value are assumed to appear once in the reaction equation.
!! Products without a specified \b yield are assumed to have a \b yield of
!! 1.0.
!!
!! Optionally, a parameter \b C may be included, and is taken to equal
!! \f$\frac{-E_a}{k_b}\f$. Note that either \b Ea or \b C may be included, but
!! not both. When neither \b Ea or \b C are included, they are assumed to be
!! 0.0. When \b A is not included, it is assumed to be 1.0, when \b D is not
!! included, it is assumed to be 300.0 K, when \b B is not included, it is
!! assumed to be 0.0, and when \b E is not included, it is assumed to be 0.0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_arrhenius_t type and associated functions. 
module pmc_rxn_arrhenius

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
#define _A_ this%condensed_data_real(1)
#define _C_ this%condensed_data_real(2)
#define _D_ this%condensed_data_real(3)
#define _B_ this%condensed_data_real(4)
#define _E_ this%condensed_data_real(5)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 5
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)

public :: rxn_arrhenius_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_arrhenius_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Time derivative contribution
    procedure :: func_contrib
    !> Jacobian matrix contribution
    procedure :: jac_contrib
    !> Calculate the rate constant
    procedure, private :: rate_const
  end type rxn_arrhenius_t

  !> Constructor for rxn_arrhenius_t
  interface rxn_arrhenius_t
    procedure :: constructor
  end interface rxn_arrhenius_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Arrhenius reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_arrhenius_t), pointer :: new_obj

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
    if (.not. this%property_set%get_real(key_name, _A_)) then
      _A_ = 1.0
    end if
    key_name = "Ea"
    if (this%property_set%get_real(key_name, temp_real)) then
      _C_ = -temp_real/const%boltzmann
      key_name = "C"
      call assert_msg(297370315, &
              .not.this%property_set%get_real(key_name, temp_real), &
              "Received both Ea and C parameter for Arrhenius equation")
    else
      key_name = "C"
      if (.not. this%property_set%get_real(key_name, _C_)) then
        _C_ = 0.0
      end if
    end if
    key_name = "D"
    if (.not. this%property_set%get_real(key_name, _D_)) then
      _D_ = 300.0
    end if
    key_name = "B"
    if (.not. this%property_set%get_real(key_name, _B_)) then
      _B_ = 0.0
    end if
    key_name = "E"
    if (.not. this%property_set%get_real(key_name, _E_)) then
      _E_ = 0.0
    end if
    
    call assert_msg(344705857, .not. ((_B_.ne.real(0.0, kind=dp)) &
            .and.(_D_.eq.real(0.0, kind=dp))), &
            "D cannot be zero if B is non-zero in Arrhenius reaction.")

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(751684145, _REACT_(i_spec).gt.0, &
              "Missing Arrhenius reactant: "//spec_name)

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
      _PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(234495887, _PROD_(i_spec).gt.0, &
              "Missing Arrhenius product: "//spec_name)

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

  end subroutine func_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the contribution of this reaction to the Jacobian matrix.
  !! The current model state is provided for species concentrations and 
  !! aerosol state. All other parameters must have been saved to the reaction 
  !! data instance during initialization.
  subroutine jac_contrib(this, phlex_state, jac_matrix)

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

  end subroutine jac_contrib

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the reaction rate constant
  real(kind=dp) function rate_const(this, phlex_state)

    !> Reaction data
    class(rxn_arrhenius_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state

    integer(kind=i_kind) :: i_spec

    if (_B_.eq.real(0.0, kind=dp)) then
      rate_const = _A_*exp(_C_/ &
          phlex_state%env_state%temp)* &
          (1.0 + _E_*const%air_std_press* &
            phlex_state%env_state%pressure)
    else
      rate_const = _A_*exp(_C_/ &
          phlex_state%env_state%temp)* &
          (phlex_state%env_state%temp/_D_)**_B_* &
          (1.0 + _E_*const%air_std_press* &
            phlex_state%env_state%pressure)
    end if

  end function rate_const

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_arrhenius
