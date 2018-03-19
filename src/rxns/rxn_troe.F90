! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_troe module.

!> \page phlex_rxn_troe Phlexible Mechanism for Chemistry: Troe Reaction
!!
!! Troe (fall-off) reaction rate constant equations take the form:
!!
!! \f[
!!   \frac{k_0[\mbox{M}]}{1+k_0[\mbox{M}]/k_{\inf}}F_C^{1+(1/N[log_{10}(k_0[\mbox{M}]/k_{\inf})]^2)^{-1}}
!! \f]
!!
!! where \f$k_0\f$ is the low-pressure limiting rate constant, \f$k_{\inf}\f$
!! is the high-pressure limiting rate constant, \f$[\mbox{M}]\f$ is the
!! density of air (taken to be \f$10^6\f$ ppm), and \f$F_C\f$ and \f$N\f$
!! are parameters that determine the shape of the fall-off curve, and are
!! typically 0.6 and 1.0, respectively (Finalyson-Pitts and Pitts, 2000;
!! Gipson and Young, 1999). \f$k_0\f$ and \f$k_{\inf}\f$ are assumed to be
!! \ref phlex_rxn_arrhenius "Arrhenius" rate constants with \f$D=300\f$ and
!! \f$E=0\f$.
!!
!! Input data for Troe equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "TROE",
!!     "k0_A" : 5.6E-12,
!!     "k0_B" : -1.8,
!!     "k0_C" : 180.0,
!!     "kinf_A" : 3.4E-12,
!!     "kinf_B" : -1.6,
!!     "kinf_C" : 104.1,
!!     "Fc"  : 0.7,
!!     "N"  : 0.9,
!!     "time unit" : "MIN",
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
!! The two sets of parameters beginning with \b k0_ and \b kinf_ are the
!! \ref phlex_rxn_arrhenius "Arrhenius" parameters for the \f$k_0\f$ and
!! \f$k_{\inf}\f$ rate constants, respectively. When not present, \b _A
!! parameters are assumed to be 1.0, \b _B to be 0.0, \b _C to be 0.0, \b Fc
!! to be 0.6 and \b N to be 1.0.
!!
!! The unit for time is assumed to be s, but inclusion of the optional
!! key-value pair \b "time unit" = "MIN" can be used to indicate a rate
!! with min as the time unit.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_troe_t type and associated functions. 
module pmc_rxn_troe

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
#define _k0_A_ this%condensed_data_real(1)
#define _k0_B_ this%condensed_data_real(2)
#define _k0_C_ this%condensed_data_real(3)
#define _kinf_A_ this%condensed_data_real(4)
#define _kinf_B_ this%condensed_data_real(5)
#define _kinf_C_ this%condensed_data_real(6)
#define _Fc_ this%condensed_data_real(7)
#define _N_ this%condensed_data_real(8)
#define _SCALING_ this%condensed_data_real(9)
#define _CONV_ this%condensed_data_real(10)
#define _RATE_CONSTANT_ this%condensed_data_real(11)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 11
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _DERIV_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + _NUM_PROD_ + x)
#define _JAC_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + 2*(_NUM_REACT_+_NUM_PROD_) + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)

public :: rxn_troe_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_troe_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
  end type rxn_troe_t

  !> Constructor for rxn_troe_t
  interface rxn_troe_t
    procedure :: constructor
  end interface rxn_troe_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Troe reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_troe_t), pointer :: new_obj

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
    class(rxn_troe_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name, string_val
    integer(kind=i_kind) :: i_spec, i_qty

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(510658779, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(852878121, this%property_set%get_property_t(key_name, reactants), &
            "Troe reaction is missing reactants")
    key_name = "products"
    call assert_msg(965196466, this%property_set%get_property_t(key_name, products), &
            "Troe reaction is missing products")

    ! Count the number of reactants (including those with a qty specified)
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(844081716, reactants%get_property_t(val=spec_props))
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
            (i_spec + 2) * (i_spec + products%size())))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_ + products%size()))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)
    
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
    key_name = "kinf_A"
    if (.not. this%property_set%get_real(key_name, _kinf_A_)) then
      _kinf_A_ = 1.0
    end if
    key_name = "kinf_B"
    if (.not. this%property_set%get_real(key_name, _kinf_B_)) then
      _kinf_B_ = 0.0
    end if
    key_name = "kinf_C"
    if (.not. this%property_set%get_real(key_name, _kinf_C_)) then
      _kinf_C_ = 0.0
    end if
    key_name = "Fc"
    if (.not. this%property_set%get_real(key_name, _Fc_)) then
      _Fc_ = 0.6
    end if
    key_name = "N"
    if (.not. this%property_set%get_real(key_name, _N_)) then
      _N_ = 1.0
    end if
    key_name = "time unit"
    _SCALING_ = real(1.0, kind=dp)
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        _SCALING_ = real(1.0d0/60.0d0, kind=dp)
      end if
    endif
  
    ! Include [M] in _k0_A_
    _k0_A_ = _k0_A_ * real(1.0d6, kind=dp)

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      _REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(595701751, _REACT_(i_spec).gt.0, &
              "Missing Troe reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(221292659, reactants%get_property_t(val=spec_props))
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
      call assert_msg(480024633, _PROD_(i_spec).gt.0, &
              "Missing Troe product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(393355097, products%get_property_t(val=spec_props))
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
    class(rxn_troe_t), intent(in) :: this
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
    class(rxn_troe_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    ! Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef _NUM_REACT_
#undef _NUM_PROD_
#undef _k0_A_
#undef _k0_B_
#undef _k0_C_
#undef _kinf_A_
#undef _kinf_B_
#undef _kinf_C_
#undef _Fc_
#undef _N_
#undef _SCALING_
#undef _CONV_
#undef _RATE_CONSTANT_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _REACT_
#undef _PROD_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _yield_
end module pmc_rxn_troe
