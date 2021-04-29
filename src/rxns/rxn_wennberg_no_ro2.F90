! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_wennberg_no_ro2 module.

!> \page camp_rxn_wennberg_no_ro2 CAMP: Wennberg NO + RO2 Reaction
!!
!! Wennberg NO + RO2 reactions are branched reactions with one branch forming
!! alkoxy radicals plus \f$\ce{NO2}\f$ and the other forming organic nitrates
!! \cite Wennberg2018 . The rate constants for each branch are based on an
!! Arrhenius rate constant and a temperature- and structure-dependent
!! branching ratio calculated as:
!!
!! \f{align}{
!!   k_{nitrate} & = \left(X e^{-Y/T}\right) \left(\frac{A(T, \mbox{[M]}, n)}{A(T, \mbox{[M]}, n) + Z}\right) \\
!!   k_{alkoxy} & = \left(X e^{-Y/T}\right)\left(\frac{Z}{Z + A(T, \mbox{[M]}, n)}\right) \\
!!   A(T, \mbox{[M]}, n) & = \frac{2 \times 10^{-22} e^n \mbox{[M]}}{1 + \frac{2 \times 10^{-22} e^n \mbox{[M]}}{0.43(T/298)^{-8}}} 0.41^{(1+[log( \frac{2 \times 10^{-22} e^n \mbox{[M]}}{0.43(T/298)^{-8}})]^2)^{-1}}
!! \f}
!!
!! where \f$T\f$ is temperature (K), [M] is the number density of air
!! (molecules \f$\mbox{cm}^{-3}\f$), \f$X\f$ and \f$Y\f$ are Arrhenius
!! parameters for the overall reaction, \f$n\f$ is the number of heavy atoms
!! in the \f$\ce{RO2}\f$ reacting species (excluding the peroxy moiety), and
!! \f$Z\f$ is defined as a function of two parameters (\f$\alpha_0, n\f$):
!!
!! \f[
!!   Z( \alpha_0, n) = A(T = 293 \mbox{K}, \mbox{[M]} = 2.45 \times 10^{19} \frac{\mbox{molec}}{\mbox{cm}^3}, n) \frac{(1-\alpha_0)}{\alpha_0}
!! \f]
!!
!! More details can be found in Wennberg et al. (2018) \cite Wennberg2018 .
!!
!! Input data for Wennberg NO + RO2 equations has the following format:
!! \code{.json}
!!   {
!!     "type" : "WENNBERG_NO_RO2",
!!     "X" : 123.45,
!!     "Y"  : 1200.0,
!!     "a0"  : 1.0e8,
!!     "n" : 6,
!!     "time unit" : "MIN",
!!     "reactants" : {
!!       "spec1" : {},
!!       "spec2" : { "qty" : 2 },
!!       ...
!!     },
!!     "alkoxy products" : {
!!       "spec3" : {},
!!       "spec4" : { "yield" : 0.65 },
!!       ...
!!     },
!!     "nitrate products" : {
!!       "spec5" : {},
!!       "spec6" : { "yield" : 0.32 },
!!       ...
!!     }
!!   }
!! \endcode
!! The key-value pairs \b reactants, and both sets of \b products are
!! required. Reactants without a \b qty value are assumed to appear once in
!! the reaction equation. Products without a specified \b yield are assumed
!! to have a \b yield of 1.0.
!!
!! When \b X is not included, it is assumed to be 1.0, when \b Y is not
!! included, it is assumed to be 0.0 K, when \b a0 is not included, it is
!! assumed to be 1.0, and when \b n is not included, it is assumed to be 0.
!! The unit for time is assumed to be s, but inclusion of the optional
!! key-value pair \b time \b unit = \b MIN can be used to indicate a rate
!! with min as the time unit.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_wennberg_no_ro2_t type and associated functions.
module pmc_rxn_wennberg_no_ro2

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
#define NUM_ALKOXY_PROD_ this%condensed_data_int(2)
#define NUM_NITRATE_PROD_ this%condensed_data_int(3)
#define X_ this%condensed_data_real(1)
#define Y_ this%condensed_data_real(2)
#define a0_ this%condensed_data_real(3)
#define n_ this%condensed_data_real(4)
#define CONV_ this%condensed_data_real(5)
#define NUM_INT_PROP_ 3
#define NUM_REAL_PROP_ 5
#define NUM_ENV_PARAM_ 2
#define REACT_(x) this%condensed_data_int(NUM_INT_PROP_ + x)
#define PROD_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_ + x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + 2*(NUM_REACT_ + NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_) + x)
#define YIELD_(x) this%condensed_data_real(NUM_REAL_PROP_ + x)

  public :: rxn_wennberg_no_ro2_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_wennberg_no_ro2_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
  end type rxn_wennberg_no_ro2_t

  !> Constructor for rxn_wennberg_no_ro2_t
  interface rxn_wennberg_no_ro2_t
    procedure :: constructor
  end interface rxn_wennberg_no_ro2_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Wennberg no_ro2 reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_wennberg_no_ro2_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = GAS_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)

    !> Reaction data
    class(rxn_wennberg_no_ro2_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells being solved simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

    type(property_t), pointer :: spec_props, reactants, alkoxy_products,      &
                                 nitrate_products
    character(len=:), allocatable :: key_name, spec_name, string_val
    integer(kind=i_kind) :: i_spec, i_qty

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real

    ! Get the species involved
    if (.not. associated(this%property_set)) call die_msg(573144252, &
            "Missing property set needed to initialize reaction")
    key_name = "reactants"
    call assert_msg(350413096, &
            this%property_set%get_property_t(key_name, reactants), &
            "Wennberg NO + RO2 reaction is missing reactants")
    key_name = "alkoxy products"
    call assert_msg(575049786, &
            this%property_set%get_property_t(key_name, alkoxy_products), &
            "Wennberg NO + RO2 reaction is missing alkoxy products")
    key_name = "nitrate products"
    call assert_msg(794422169, &
            this%property_set%get_property_t(key_name, nitrate_products), &
            "Wennberg NO + RO2 reaction is missing nitrate products")

    ! Count the number of reactants (including those with a qty specified)
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(626170799, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) i_spec = i_spec+temp_int-1
      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_ + &
            (i_spec + 2) * (i_spec + alkoxy_products%size() + &
                            nitrate_products%size())))
    allocate(this%condensed_data_real(NUM_REAL_PROP_ + &
                            alkoxy_products%size() + nitrate_products%size()))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save space for the environment dependent parameters
    this%num_env_params = NUM_ENV_PARAM_

    ! Save the size of the reactant and product arrays (for reactions where
    ! these can vary)
    NUM_REACT_ = i_spec
    NUM_ALKOXY_PROD_  = alkoxy_products%size()
    NUM_NITRATE_PROD_ = nitrate_products%size()

    ! Set the #/cc -> ppm conversion prefactor
    CONV_ = const%avagadro / const%univ_gas_const * 10.0d0**(-12.0d0)

    ! Get reaction parameters (it might be easiest to keep these at the
    ! beginning of the condensed data array, so they can be accessed using
    ! compliler flags)
    key_name = "X"
    if (.not. this%property_set%get_real(key_name, X_)) then
      X_ = 1.0
    end if
    key_name = "time unit"
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        X_ = X_ / 60.0
      end if
    endif
    key_name = "Y"
    if (.not. this%property_set%get_real(key_name, Y_)) then
      Y_ = 0.0
    end if
    key_name = "a0"
    if (.not. this%property_set%get_real(key_name, a0_)) then
      a0_ = 1.0
    end if
    key_name = "n"
    if (.not. this%property_set%get_real(key_name, n_)) then
      n_ = 0.0
    end if

    ! Get the indices and chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(716430972, REACT_(i_spec).gt.0, &
              "Missing Wennberg NO + RO2 reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(493699816, reactants%get_property_t(val=spec_props))
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
    call alkoxy_products%iter_reset()
    i_spec = 1
    do while (alkoxy_products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(825390544, PROD_(i_spec).gt.0, &
              "Missing Wennberg NO + RO2 alkoxy product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(879870330, alkoxy_products%get_property_t(val=spec_props))
      key_name = "yield"
      if (spec_props%get_real(key_name, temp_real)) then
        YIELD_(i_spec) = temp_real
      else
        YIELD_(i_spec) = 1.0
      end if

      call alkoxy_products%iter_next()
      i_spec = i_spec + 1
    end do
    call nitrate_products%iter_reset()
    do while (nitrate_products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      ! Make sure the species exists
      call assert_msg(590677535, PROD_(i_spec).gt.0, &
              "Missing Wennberg NO + RO2 nitrate product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(815314225, nitrate_products%get_property_t(val=spec_props))
      key_name = "yield"
      if (spec_props%get_real(key_name, temp_real)) then
        YIELD_(i_spec) = temp_real
      else
        YIELD_(i_spec) = 1.0
      end if

      call nitrate_products%iter_next()
      i_spec = i_spec + 1
    end do
    call assert(699637107, i_spec - 1 .eq. &
                           NUM_ALKOXY_PROD_ + NUM_NITRATE_PROD_ )

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_wennberg_no_ro2_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_wennberg_no_ro2
