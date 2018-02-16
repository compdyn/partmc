! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_CMAQ_OH_HNO3 module.

!> \page phlex_rxn_CMAQ_OH_HNO3 Phlexible Mechanism for Chemistry: Special CMAQ Reaction for OH+HNO3
!!
!! A special CMAQ rate constant for the reactions:
!!
!! \f[
!!   \mbox{OH} + \mbox{HNO3} -> \mbox{NO_3} + \mbox{H_2O}
!! \f]
!!
!! takes the form:
!!
!! \f[
!!   k=k_0+(\frac{k_3[\mbox{M}]}{1+k_3[\mbox{M}]/k_2})
!! \f]
!!
!! where \f$k_0\f$, \f$k_2\f$ and \f$k_3\f$ are \ref phlex_rxn_arrhenius
!! "Arrhenius" rate constants with \f$D=300\f$ and \f$E=0\f$, and
!! \f$[\mbox{M}]\f$ is the density of air (taken to be \f$10^6\f$ ppm;
!! Gipson and Young, 1999).
!!
!! Input data for CMAQ OH+HNO3 equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "CMAQ_OH_HNO3",
!!     "k0_A" : 5.6E-12,
!!     "k0_B" : -1.8,
!!     "k0_C" : 180.0,
!!     "k2_A" : 3.4E-12,
!!     "k2_B" : -1.6,
!!     "k2_C" : 104.1,
!!     "k3_A" : 3.2E-11,
!!     "k3_B" : -1.5,
!!     "k3_C" : 92.0,
!!     "time unit" : "MIN"
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
!! The three sets of parameters beginning with \b k0_, \b k2_, and \b k3_, are
!! the \ref phlex_rxn_arrhenius "Arrhenius" parameters for the \f$k_0\f$,
!! \f$k_2\f$ and \f$k_3\f$ rate constants, respectively. When not present,
!! \b _A parameters are assumed to be 1.0, \b _B to be 0.0, and \b _C to be
!! 0.0.
!!
!! The unit for time is assumed to be s, but inclusion of the optional
!! key-value pair \b "time unit" = "MIN" can be used to indicate a rate
!! with min as the time unit.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_CMAQ_OH_HNO3_t type and associated functions. 
module pmc_rxn_CMAQ_OH_HNO3

  use pmc_constants,                        only: const
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
#define _k0_A_ this%condensed_data_real(1)
#define _k0_B_ this%condensed_data_real(2)
#define _k0_C_ this%condensed_data_real(3)
#define _k2_A_ this%condensed_data_real(4)
#define _k2_B_ this%condensed_data_real(5)
#define _k2_C_ this%condensed_data_real(6)
#define _k3_A_ this%condensed_data_real(7)
#define _k3_B_ this%condensed_data_real(8)
#define _k3_C_ this%condensed_data_real(9)
#define _CONV_ this%condensed_data_real(10)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 10
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _JACID_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + _NUM_PROD_ + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)

public :: rxn_CMAQ_OH_HNO3_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_CMAQ_OH_HNO3_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Indicate which Jacobian elements are populated by this reactions
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
  end type rxn_CMAQ_OH_HNO3_t

  !> Constructor for rxn_CMAQ_OH_HNO3_t
  interface rxn_CMAQ_OH_HNO3_t
    procedure :: constructor
  end interface rxn_CMAQ_OH_HNO3_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for CMAQ OH+HNO3 reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_CMAQ_OH_HNO3_t), pointer :: new_obj

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
    class(rxn_CMAQ_OH_HNO3_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name, string_val
    integer(kind=i_kind) :: i_spec, i_qty

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(255324828, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(250060521, this%property_set%get_property_t(key_name, reactants), &
            "CMAQ OH+HNO3 reaction is missing reactants")
    key_name = "products"
    call assert_msg(304540307, this%property_set%get_property_t(key_name, products), &
            "CMAQ OH+HNO3 reaction is missing products")

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
    allocate(this%condensed_data_int(_NUM_INT_PROP_+i_spec+products%size() &
            + i_spec * (i_spec + products%size())))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_+products%size()))
    
    ! Save the size of the reactant and product arrays (for reactions where these
    ! can vary)
    _NUM_REACT_ = i_spec
    _NUM_PROD_ = products%size()

    ! Set the #/cc -> ppm conversion prefactor
    _CONV_ = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0)

    ! Get reaction parameters (it might be easiest to keep these at the beginning
    ! of the condensed data array, so they can be accessed using compliler flags)
    key_name = "k0_A"
    if (.not. this%property_set%get_real(key_name, _k0_A_)) then
      _k0_A_ = 1.0
    end if
    key_name = "k0_B"
    if (.not. this%property_set%get_real(key_name, _k0_B_)) then
      _k0_B_ = 0.0
    end if
    key_name = "k0_C"
    if (.not. this%property_set%get_real(key_name, _k0_C_)) then
      _k0_C_ = 0.0
    end if
    key_name = "k2_A"
    if (.not. this%property_set%get_real(key_name, _k2_A_)) then
      _k2_A_ = 1.0
    end if
    key_name = "k2_B"
    if (.not. this%property_set%get_real(key_name, _k2_B_)) then
      _k2_B_ = 0.0
    end if
    key_name = "k2_C"
    if (.not. this%property_set%get_real(key_name, _k2_C_)) then
      _k2_C_ = 0.0
    end if
    key_name = "k3_A"
    if (.not. this%property_set%get_real(key_name, _k3_A_)) then
      _k3_A_ = 1.0
    end if
    key_name = "k3_B"
    if (.not. this%property_set%get_real(key_name, _k3_B_)) then
      _k3_B_ = 0.0
    end if
    key_name = "k3_C"
    if (.not. this%property_set%get_real(key_name, _k3_C_)) then
      _k3_C_ = 0.0
    end if
    key_name = "time unit"
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        _k0_A_ = _k0_A_ / 60.0
        _k2_A_ = _k2_A_ / 60.0
        _k3_A_ = _k3_A_ / 60.0
      end if
    endif
    
    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(665522163, _REACT_(i_spec).gt.0, &
              "Missing CMAQ OH HNO3 reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(796763915, reactants%get_property_t(val=spec_props))
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

    ! Get the indices and chemical properties for the products
    call products%iter_reset()
    i_spec = 1
    do while (products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(881535773, _PROD_(i_spec).gt.0, &
              "Missing CMAQ OH HNO3 product: "//spec_name)

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

  !> Indicate which Jacobian elements are populated by this reaction
  subroutine get_used_jac_elem(this, use_jac_elem)

    !> Reaction data
    class(rxn_CMAQ_OH_HNO3_t), intent(in) :: this
    !> Matrix of flags indicating whether to include Jacobian elements
    !! in the sparse matrix
    integer(kind=i_kind), allocatable, intent(inout) :: use_jac_elem(:,:)

    integer(kind=i_kind) :: i_spec_d, i_spec_i

    do i_spec_d=1, _NUM_REACT_
      do i_spec_i=1, _NUM_REACT_
        use_jac_elem(_REACT_(i_spec_d), &
            _REACT_(i_spec_i)) = 1
      end do
    end do

    do i_spec_d=1, _NUM_PROD_
      do i_spec_i=1, _NUM_REACT_
        use_jac_elem(_PROD_(i_spec_d), &
            _REACT_(i_spec_i)) = 1
      end do
    end do

  end subroutine get_used_jac_elem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update indexes used during integration
  subroutine update_integration_ids(this, integration_data)

    !> Reaction data
    class(rxn_CMAQ_OH_HNO3_t), intent(inout) :: this
    !> Integration data
    class(integration_data_t), intent(in) :: integration_data

    integer(kind=i_kind) :: i_spec_d, i_spec_i, jac_id

    jac_id = 0
    do i_spec_d=1, _NUM_REACT_
      do i_spec_i=1, _NUM_REACT_
        jac_id = jac_id + 1
        _JACID_(jac_id) = &
            integration_data%get_jac_elem_id(&
            _REACT_(i_spec_d), &
            _REACT_(i_spec_i))
      end do
    end do

    do i_spec_d=1, _NUM_PROD_
      do i_spec_i=1, _NUM_REACT_
        jac_id = jac_id + 1
        _JACID_(jac_id) = &
            integration_data%get_jac_elem_id(&
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
  !!     [reactant1] * [reactant2] *A*exp(B+temp/theta)
  !!
  subroutine  func_contrib(this, integration_data)

    !> Reaction data
    class(rxn_CMAQ_OH_HNO3_t), intent(in) :: this
    !> Integration data
    class(integration_data_t), pointer, intent(inout) :: integration_data

    integer(kind=i_kind) :: i_spec
    real(kind=dp) :: rate

    rate = this%rate_const(integration_data%phlex_state)
    do i_spec=1, _NUM_REACT_
      rate = rate * integration_data%phlex_state%state_var(_REACT_(i_spec))
    end do
    
    if (rate.eq.real(0.0, kind=dp)) return

    do i_spec=1, _NUM_REACT_
      call integration_data%add_deriv_contrib(_REACT_(i_spec), -rate)
    end do

    do i_spec=1, _NUM_PROD_
      call integration_data%add_deriv_contrib(_PROD_(i_spec), &
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
    class(rxn_CMAQ_OH_HNO3_t), intent(in) :: this
    !> Integration data
    class(integration_data_t), pointer, intent(inout) :: integration_data

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
                -rate) 
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
    class(rxn_CMAQ_OH_HNO3_t), intent(in) :: this
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

  !> Calculate the reaction rate constant
  real(kind=dp) function rate_const(this, phlex_state)

    !> Reaction data
    class(rxn_CMAQ_OH_HNO3_t), intent(in) :: this
    !> Current model state
    type(phlex_state_t), intent(in) :: phlex_state

    integer(kind=i_kind) :: i_spec
    real(kind=dp) :: k0, k2, k3
    real(kind=dp), parameter :: M = 1.0d6
    real(kind=dp), parameter :: D = 300.0d0
    real(kind=dp) :: conv

    ! k0
    if (_k0_B_.eq.real(0.0, kind=dp)) then
      k0 = _k0_A_*exp(_k0_C_/ &
          phlex_state%env_state%temp)
    else
      k0 = _k0_A_*exp(_k0_C_/ &
          phlex_state%env_state%temp)* &
          (phlex_state%env_state%temp/D)**_k0_B_
    end if

    ! k2
    if (_k2_B_.eq.real(0.0, kind=dp)) then
      k2 = _k2_A_*exp(_k2_C_/ &
          phlex_state%env_state%temp)
    else
      k2 = _k2_A_*exp(_k2_C_/ &
          phlex_state%env_state%temp)* &
          (phlex_state%env_state%temp/D)**_k2_B_
    end if

    ! k3
    if (_k3_B_.eq.real(0.0, kind=dp)) then
      k3 = _k3_A_*exp(_k3_C_/ &
          phlex_state%env_state%temp)
    else
      k3 = _k3_A_*exp(_k3_C_/ &
          phlex_state%env_state%temp)* &
          (phlex_state%env_state%temp/D)**_k3_B_
    end if

    ! Convert from #/cc -> ppm
    conv = _CONV_ * &
            phlex_state%env_state%pressure / phlex_state%env_state%temp
    
    k3 = k3*M*conv
    rate_const = k0 + k3 / (1.0d0 + k3/k2)
    rate_const = rate_const * conv**(_NUM_REACT_-1)

  end function rate_const

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_CMAQ_OH_HNO3
