! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_arrhenius module.

!> \page phlex_rxn_arrhenius Phlexible Module for Chemistry: Arrhenius Reaction
!!
!! Arrhenius-like reaction rate constant equations can take the form:
!!
!! \f[
!!   Ae^{(\frac{-E_a}{k_bT})}(\frac{T}{D})^B(1.0+E*P)
!! \f]
!!
!! where \f$A\f$ is the pre-exponential factor 
!! (\f$(\mbox{\SI[mode=text]{}{\#.cm^{-3}}})^{-(n-1)}\f$\SI{\per\second}),
!! \f$n\f$ is the number of
!! reactants, \f$E_a\f$ is the activation energy (J), \f$k_b\f$ is the
!! Boltzmann constant (J/K), \f$D\f$ (K), \f$B\f$ (unitless) and \f$E\f$
!! (\f$Pa^{-1}\f$) are reaction parameters, \f$T\f$ is the temperature (K),
!! and \f$P\f$ is the pressure (Pa). The first two terms are described in
!! Finlayson-Pitts and Pitts (2000). The final term is included to accomodate
!! CMAQ EBI solver type 7 rate constants.
!!
!! Input data for Arrhenius equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "ARRHENIUS",
!!     "A" : 123.45,
!!     "Ea" : 123.45,
!!     "B"  : 1.3,
!!     "D"  : 300.0,
!!     "E"  : 0.6E-5,
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
!! Optionally, a parameter \b C may be included, and is taken to equal
!! \f$\frac{-E_a}{k_b}\f$. Note that either \b Ea or \b C may be included, but
!! not both. When neither \b Ea or \b C are included, they are assumed to be
!! 0.0. When \b A is not included, it is assumed to be 1.0, when \b D is not
!! included, it is assumed to be 300.0 K, when \b B is not included, it is
!! assumed to be 0.0, and when \b E is not included, it is assumed to be 0.0.
!! The unit for time is assumed to be s, but inclusion of the optional
!! key-value pair \b "time unit" = "MIN" can be used to indicate a rate
!! with min as the time unit.

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
  use pmc_aero_rep_data

  implicit none
  private

#define _NUM_REACT_ this%condensed_data_int(1)
#define _NUM_PROD_ this%condensed_data_int(2)
#define _A_ this%condensed_data_real(1)
#define _B_ this%condensed_data_real(2)
#define _C_ this%condensed_data_real(3)
#define _D_ this%condensed_data_real(4)
#define _E_ this%condensed_data_real(5)
#define _CONV_ this%condensed_data_real(6)
#define _RATE_CONSTANT_ this%condensed_data_real(7)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 7
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + x)
#define _DERIV_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_REACT_ + _NUM_PROD_ + x)
#define _JAC_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + 2*(_NUM_REACT_+_NUM_PROD_) + x)
#define _yield_(x) this%condensed_data_real(_NUM_REAL_PROP_ + x)

  public :: rxn_arrhenius_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_arrhenius_t
  contains
    !> Reaction initialization
    procedure :: initialize
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
  subroutine initialize(this, chem_spec_data, aero_rep)
    
    !> Reaction data
    class(rxn_arrhenius_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

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
    key_name = "A"
    if (.not. this%property_set%get_real(key_name, _A_)) then
      _A_ = 1.0
    end if
    key_name = "time unit"
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        _A_ = _A_ / 60.0
      end if
    endif
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

#undef _NUM_REACT_
#undef _NUM_PROD_
#undef _A_
#undef _B_
#undef _C_
#undef _D_
#undef _E_
#undef _CONV_
#undef _RATE_CONSTANT_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _REACT_
#undef _PROD_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _yield_
end module pmc_rxn_arrhenius
