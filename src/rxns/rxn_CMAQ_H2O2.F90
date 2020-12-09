! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_CMAQ_H2O2 module.

!> \page camp_rxn_CMAQ_H2O2 CAMP: Special CMAQ Reaction for H2O2
!!
!! For the reactions:
!!
!! \f[\mbox{\ch{
!!   HO2 + HO2 ->[k] H2O2
!! }}\f]
!! \f[\mbox{\ch{
!!   HO2 + HO2 + H2O ->[k] H2O2
!! }}\f]
!!
!! CMAQ rate constants are calculated as:
!!
!! \f[
!!   k=k_1+k_2[\mbox{M}]
!! \f]
!!
!! where \f$k_1\f$ and \f$k_2\f$ are \ref camp_rxn_arrhenius "Arrhenius" rate
!! constants with \f$D=300\f$ and \f$E=0\f$, and \f$[\mbox{M}]\f$ is the
!! concentration of air (\f$10^6\f$ ppm) \cite Gipson.
!!
!! Input data for CMAQ \f$\mbox{\ch{H2O2}}\f$ reactions have the following
!! format:
!! \code{.json}
!!   {
!!     "type" : "CMAQ_H2O2",
!!     "k1_A" : 5.6E-12,
!!     "k1_B" : -1.8,
!!     "k1_C" : 180.0,
!!     "k2_A" : 3.4E-12,
!!     "k2_B" : -1.6,
!!     "k2_C" : 104.1,
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
!! The two sets of parameters beginning with \b k1_ and \b k2_ are the
!! \ref camp_rxn_arrhenius "Arrhenius" parameters for the \f$k_1\f$ and
!! \f$k_2\f$ rate constants, respectively. When not present, \b _A
!! parameters are assumed to be 1.0, \b _B to be 0.0, and \b _C to be 0.0.
!!
!! The unit for time is assumed to be s, but inclusion of the optional
!! key-value pair \b time \b unit = \b MIN can be used to indicate a rate
!! with min as the time unit.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_CMAQ_H2O2_t type and associated functions.
module pmc_rxn_CMAQ_H2O2

  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_camp_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, die_msg

  implicit none
  private

#define NUM_REACT_ this%condensed_data_int(1)
#define NUM_PROD_ this%condensed_data_int(2)
#define k1_A_ this%condensed_data_real(1)
#define k1_B_ this%condensed_data_real(2)
#define k1_C_ this%condensed_data_real(3)
#define k2_A_ this%condensed_data_real(4)
#define k2_B_ this%condensed_data_real(5)
#define k2_C_ this%condensed_data_real(6)
#define CONV_ this%condensed_data_real(7)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 7
#define NUM_ENV_PARAM_ 1
#define REACT_(x) this%condensed_data_int(NUM_INT_PROP_ + x)
#define PROD_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + 2*(NUM_REACT_ + NUM_PROD_) + x)
#define YIELD_(x) this%condensed_data_real(NUM_REAL_PROP_ + x)

public :: rxn_CMAQ_H2O2_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_CMAQ_H2O2_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
  end type rxn_CMAQ_H2O2_t

  !> Constructor for rxn_CMAQ_H2O2_t
  interface rxn_CMAQ_H2O2_t
    procedure :: constructor
  end interface rxn_CMAQ_H2O2_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for CMAQ H2O2 reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_CMAQ_H2O2_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = GAS_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)

    !> Reaction data
    class(rxn_CMAQ_H2O2_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name, string_val
    integer(kind=i_kind) :: i_spec, i_qty

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(255324828, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(940945637, &
            this%property_set%get_property_t(key_name, reactants), &
            "CMAQ H2O2 reaction is missing reactants")
    key_name = "products"
    call assert_msg(435739232, &
            this%property_set%get_property_t(key_name, products), &
            "CMAQ H2O2 reaction is missing products")

    ! Count the number of reactants (including those with a qty specified)
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(952475195, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) i_spec = i_spec+temp_int-1
      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Allocate space in the condensed data arrays
    ! Space in this example is allocated for two sets of inidices for the
    ! reactants and products, one molecular property for each reactant,
    ! yields for the products and three reaction parameters.
    allocate(this%condensed_data_int(NUM_INT_PROP_ + &
            (i_spec + 2) * (i_spec + products%size())))
    allocate(this%condensed_data_real(NUM_REAL_PROP_ + products%size()))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save space for the environment dependent parameters
    this%num_env_params = NUM_ENV_PARAM_

    ! Save the size of the reactant and product arrays (for reactions where
    ! these can vary)
    NUM_REACT_ = i_spec
    NUM_PROD_ = products%size()

    ! Set the #/cc -> ppm conversion prefactor
    CONV_ = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0)

    ! Get reaction parameters (it might be easiest to keep these at the
    ! beginning of the condensed data array, so they can be accessed using
    ! compliler flags)
    key_name = "k1_A"
    if (.not. this%property_set%get_real(key_name, k1_A_)) then
      k1_A_ = 1.0
    end if
    key_name = "k1_B"
    if (.not. this%property_set%get_real(key_name, k1_B_)) then
      k1_B_ = 0.0
    end if
    key_name = "k1_C"
    if (.not. this%property_set%get_real(key_name, k1_C_)) then
      k1_C_ = 0.0
    end if
    key_name = "k2_A"
    if (.not. this%property_set%get_real(key_name, k2_A_)) then
      k2_A_ = 1.0
    end if
    key_name = "k2_B"
    if (.not. this%property_set%get_real(key_name, k2_B_)) then
      k2_B_ = 0.0
    end if
    key_name = "k2_C"
    if (.not. this%property_set%get_real(key_name, k2_C_)) then
      k2_C_ = 0.0
    end if
    key_name = "time unit"
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        k1_A_ = k1_A_ / 60.0
        k2_A_ = k2_A_ / 60.0
      end if
    endif

    ! Include the multiplication of [M]*k2 into k2_A_
    k2_A_ = k2_A_ * real(1.0d6, kind=dp)

    ! Get the indices and chemical properties for the reactants
    !write(*,*) "CMAQ_H2O2 species,REACT_ID"
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      !write(*,*) chem_spec_data%gas_state_name(chem_spec_data%gas_state_id(spec_name)), &
      !        REACT_(i_spec)

      ! Make sure the species exists
      call assert_msg(345360993, REACT_(i_spec).gt.0, &
              "Missing CMAQ H2O2 reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(796763915, reactants%get_property_t(val=spec_props))
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

    ! Get the indices and chemical properties for the products
    call products%iter_reset()
    i_spec = 1
    do while (products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      !write(*,*) chem_spec_data%gas_state_name(chem_spec_data%gas_state_id(spec_name)), &
      !        PROD_(i_spec)

      ! Make sure the species exists
      call assert_msg(234948182, PROD_(i_spec).gt.0, &
              "Missing CMAQ H2O2 product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(267035567, products%get_property_t(val=spec_props))
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

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_CMAQ_H2O2_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_CMAQ_H2O2
