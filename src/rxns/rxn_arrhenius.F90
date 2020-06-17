! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_arrhenius module.

!> \page camp_rxn_arrhenius CAMP: Arrhenius Reaction
!!
!! Arrhenius-like reaction rate constant equations are calculated as follows:
!!
!! \f[
!!   Ae^{(\frac{-E_a}{k_bT})}(\frac{T}{D})^B(1.0+E*P)
!! \f]
!!
!! where \f$A\f$ is the pre-exponential factor
!! (\f$(\mbox{\si{\#.cm^{-3}}})^{-(n-1)}\mbox{\si{\per\second}}\f$),
!! \f$n\f$ is the number of reactants, \f$E_a\f$ is the activation energy (J),
!! \f$k_b\f$ is the Boltzmann constant (J/K), \f$D\f$ (K), \f$B\f$ (unitless)
!! and \f$E\f$ (\f$Pa^{-1}\f$) are reaction parameters, \f$T\f$ is the
!! temperature (K), and \f$P\f$ is the pressure (Pa). The first two terms are
!! described in Finlayson-Pitts and Pitts (2000) \cite Finlayson-Pitts2000 .
!! The final term is included to accomodate CMAQ EBI solver type 7 rate
!! constants.
!!
!! Input data for Arrhenius equations has the following format:
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
!! key-value pair \b time \b unit = \b MIN can be used to indicate a rate
!! with min as the time unit.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_arrhenius_t type and associated functions.
module pmc_rxn_arrhenius

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
#define A_ this%condensed_data_real(1)
#define B_ this%condensed_data_real(2)
#define C_ this%condensed_data_real(3)
#define D_ this%condensed_data_real(4)
#define E_ this%condensed_data_real(5)
#define CONV_ this%condensed_data_real(6)
#define NUM_INT_PROP_ 2
#define NUM_REAL_PROP_ 6
#define NUM_ENV_PARAM_ 1
#define REACT_(x) this%condensed_data_int(NUM_INT_PROP_ + x)
#define PROD_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + NUM_REACT_ + NUM_PROD_ + x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_ + 2*(NUM_REACT_+NUM_PROD_) + x)
#define YIELD_(x) this%condensed_data_real(NUM_REAL_PROP_ + x)

  public :: rxn_arrhenius_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_arrhenius_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
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
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)

    !> Reaction data
    class(rxn_arrhenius_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells being solved simultaneously
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
    call assert_msg(250060521, &
            this%property_set%get_property_t(key_name, reactants), &
            "Arrhenius reaction is missing reactants")
    key_name = "products"
    call assert_msg(304540307, &
            this%property_set%get_property_t(key_name, products), &
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
    key_name = "rxn id"
    if (this%property_set%get_string(key_name, string_val)) then
      write (*,*) "Arrhenius RXN id:", string_val
    endif
    key_name = "A"
    if (.not. this%property_set%get_real(key_name, A_)) then
      A_ = 1.0
    end if
    key_name = "time unit"
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        A_ = A_ / 60.0
      end if
    endif
    key_name = "Ea"
    if (this%property_set%get_real(key_name, temp_real)) then
      C_ = -temp_real/const%boltzmann
      key_name = "C"
      call assert_msg(297370315, &
              .not.this%property_set%get_real(key_name, temp_real), &
              "Received both Ea and C parameter for Arrhenius equation")
    else
      key_name = "C"
      if (.not. this%property_set%get_real(key_name, C_)) then
        C_ = 0.0
      end if
    end if
    key_name = "D"
    if (.not. this%property_set%get_real(key_name, D_)) then
      D_ = 300.0
    end if
    key_name = "B"
    if (.not. this%property_set%get_real(key_name, B_)) then
      B_ = 0.0
    end if
    key_name = "E"
    if (.not. this%property_set%get_real(key_name, E_)) then
      E_ = 0.0
    end if

    write (*,*) "A_,B_,C_,D_,E_,CONV_"
    write (*,*) A_,B_,C_,D_,E_,CONV_

    call assert_msg(344705857, .not. ((B_.ne.real(0.0, kind=dp)) &
            .and.(D_.eq.real(0.0, kind=dp))), &
            "D cannot be zero if B is non-zero in Arrhenius reaction.")

    ! Get the indices and chemical properties for the reactants
    write(*,*) "arrhenius species,REACT_ID"
    call reactants%iter_reset()
    i_spec = 1
    do while (reactants%get_key(spec_name))

      ! Save the index of this species in the state variable array
      REACT_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      write(*,*) chem_spec_data%gas_state_name(chem_spec_data%gas_state_id(spec_name)), &
              REACT_(i_spec)

      ! Make sure the species exists
      call assert_msg(751684145, REACT_(i_spec).gt.0, &
              "Missing Arrhenius reactant: "//spec_name)

      ! Get properties included with this reactant in the reaction data
      call assert(796763915, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) then
        do i_qty = 1, temp_int - 1
          REACT_(i_spec + i_qty) = REACT_(i_spec)
        end do
        i_spec = i_spec + temp_int - 1
      end if

      !write(*,*) "after qty,", chem_spec_data%gas_state_name(chem_spec_data%gas_state_id(spec_name)), &
      !        REACT_(i_spec)

      call reactants%iter_next()
      i_spec = i_spec + 1
    end do

    write(*,*) "arrhenius species, PROD_ID"

    ! Get the indices and chemical properties for the products
    call products%iter_reset()
    i_spec = 1
    do while (products%get_key(spec_name))

      ! Save the index of this species in the state variable array
      PROD_(i_spec) = chem_spec_data%gas_state_id(spec_name)

      write(*,*) chem_spec_data%gas_state_name(chem_spec_data%gas_state_id(spec_name)), &
              PROD_(i_spec)

      ! Make sure the species exists
      call assert_msg(234495887, PROD_(i_spec).gt.0, &
              "Missing Arrhenius product: "//spec_name)

      ! Get properties included with this product in the reaction data
      call assert(451185800, products%get_property_t(val=spec_props))
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
    type(rxn_arrhenius_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_arrhenius
