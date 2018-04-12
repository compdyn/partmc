! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_aqueous_equilibrium module.

!> \page phlex_rxn_aqueous_equilibrium Phlexible Mechanism for Chemistry: Phase-Transfer Reaction
!!
!! Aqueous equilibrium reactions are calculated as forward and reverse
!! reactions, based on a provided reverse reaction rate constant and the
!! equilibrium constant that takes the form:
!!
!! \f[
!!   Ae^{C({1/T-1/298})}
!! \f]
!!
!! where \f$A\f$ is the pre-exponential factor (\f$s^{-1}\f$), \f$C\f$ is a
!! constant and \f$T\f$ is the temperature (K). Uptake kinetics are based on
!! the particle size, the gas-phase species diffusion coefficient and
!! molecular weight, and \f$N^{*}\f$, which is used to calculate the mass
!! accomodation coefficient. Details of the calculations can be found in:
!!
!! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended
!! and condensed tropospheric aqueous mechanism and its application."
!! J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202
!!
!! Input data for Aqueous equilibrium equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "AQUEOUS_EQUILIBRIUM",
!!     "A" : 123.45,
!!     "C" : 123.45,
!!     "k_reverse" : 123.45,
!!     "phase" : "my aqueous phase",
!!     "time unit" : "MIN",
!!     "aqueous-phase water" : "H2O_aq",
!!     "reactants" : {
!!       "spec1" : {},
!!       "spec2" : { "qty" : 2 },
!!       ...
!!     },
!!     "products" : {
!!       "spec3" : {},
!!       "spec4" : { "qty" : 0.65 },
!!       ...
!!     }
!!     ...
!!   }
!! \endcode
!! The key-value pairs \b reactants and \b products are required. Reactants
!! and products without a \b qty value are assumed to appear once in the
!! reaction equation. Reactant and product species must be present in the
!! specified phase and include a \b "molecular weight" parameter. The
!! parameter \b  "aqueous-phase water" is required and must be the name of the
!! aerosol-phase species that is used for water.
!!
!! When \b A is not included, it is assumed
!! to be 1.0, when \b C is not included, it is assumed to be 0.0. The reverse
!! reaction rate constant \b k_reverse is required.
!!
!! The unit for time is assumed to be s, but inclusion of the optional 
!! key-value pair \b "time unit" = "MIN" can be used to indicate a rate
!! with min as the time unit.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_aqueous_equilibrium_t type and associated functions. 
module pmc_rxn_aqueous_equilibrium

  use pmc_constants,                        only: const
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, die_msg, &
                                                  string_t
  use pmc_rxn_data
  use pmc_chem_spec_data
  use pmc_property
  use pmc_phlex_state
  use pmc_aero_rep_data
  use pmc_aero_phase_data

  implicit none
  private

#define _NUM_REACT_ this%condensed_data_int(1)
#define _NUM_PROD_ this%condensed_data_int(2)
#define _NUM_AERO_PHASE_ this%condensed_data_int(3)
#define _A_ this%condensed_data_real(1)
#define _C_ this%condensed_data_real(2)
#define _k_reverse_ this%condensed_data_real(3)
#define _k_forward_ this%condensed_data_real(4)
#define _NUM_INT_PROP_ 3
#define _NUM_REAL_PROP_ 4
#define _REACT_(x) this%condensed_data_int(_NUM_INT_PROP_+x)
#define _REACT_ACT_COEFF_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_REACT_*_NUM_AERO_PHASE_+x)
#define _PROD_(x) this%condensed_data_int(_NUM_INT_PROP_+2*_NUM_REACT_*_NUM_AERO_PHASE_+x)
#define _PROD_ACT_COEFF_(x) this%condensed_data_int(_NUM_INT_PROP_+(_NUM_PROD_+2*_NUM_REACT_)*_NUM_AERO_PHASE_+x)
#define _WATER_(x) this%condensed_data_int(_NUM_INT_PROP_+2*(_NUM_REACT_+_NUM_PROD_)*_NUM_AERO_PHASE_+x)
#define _WATER_ACT_COEFF_(x) this%condensed_data_int(_NUM_INT_PROP_+(2*(_NUM_REACT_+_NUM_PROD_)+1)*_NUM_AERO_PHASE_+x)
#define _DERIV_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+(2*(_NUM_REACT_+_NUM_PROD_)+2)*_NUM_AERO_PHASE_+x)
#define _JAC_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+(3*(_NUM_REACT_+_NUM_PROD_)+2)*_NUM_AERO_PHASE_+x)
#define _mass_frac_TO_M_(x) this%condensed_data_real(_NUM_REAL_PROP_+x)

  public :: rxn_aqueous_equilibrium_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_aqueous_equilibrium_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
  end type rxn_aqueous_equilibrium_t

  !> Constructor for rxn_aqueous_equilibrium_t
  interface rxn_aqueous_equilibrium_t
    procedure :: constructor
  end interface rxn_aqueous_equilibrium_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Aqueous equilibrium reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_aqueous_equilibrium_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

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
    class(rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name, water_name, phase_name, &
            string_val
    integer(kind=i_kind) :: i_spec, j_spec, i_qty, i_aero_rep, i_aero_phase, n_aero_ids
    integer(kind=i_kind) :: i_aero_id, num_spec_per_phase, num_phase, num_react, num_prod
    class(string_t), allocatable :: unique_names(:)
    class(string_t), allocatable :: react_names(:), prod_names(:)
    integer(kind=i_kind), allocatable :: aero_spec_ids(:)
    integer(kind=i_kind), allocatable :: water_spec_ids(:)

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real, N_star

    ! Get the property set
    if (.not. associated(this%property_set)) call die_msg(206493887, &
            "Missing property set needed to initialize reaction")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(138126714, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in aqueous-equilibrium reaction")

    ! Get the aerosol-phase water name
    key_name = "aerosol-phase water"
    call assert_msg(583327500, &
            this%property_set%get_string(key_name, water_name), &
            "Missing aerosol-phase water in aqueous-equilibrium reaction")

    ! Get reactant species
    key_name = "reactants"
    call assert_msg(237749385, &
            this%property_set%get_property_t(key_name, reactants), &
            "Missing reactant species in aqueous-equilibrium reaction")

    ! Get product species
    key_name = "products"
    call assert_msg(350520025, &
            this%property_set%get_property_t(key_name, products), &
            "Missing product species in aqueous-equilibrium reaction")

    ! Count the number of species per phase instance (including species
    ! with a "qty" specified)
    num_react = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(422080799, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) &
              num_react = num_react + temp_int - 1
      call reactants%iter_next()
      num_react = num_react + 1
    end do
    num_prod = 0
    do while (products%get_key(spec_name))
      ! Get properties included with this product in the reaction data
      call assert(971363961, products%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) &
              num_prod = num_prod + temp_int - 1
      call products%iter_next()
      num_prod = num_prod + 1
    end do
    num_spec_per_phase = num_prod + num_react

    ! Check for aerosol representations
    call assert_msg(191050890, associated(aero_rep), &
            "Missing aerosol representation for aqueous equilibrium reaction")
    call assert_msg(868319733, size(aero_rep).gt.0, &
            "Missing aerosol representation for aqueous equilibrium reaction")
    
    ! Count the instances of this phase/species pair
    num_phase = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Check for the specified phase in this aero rep
      if (aero_rep(i_aero_rep)%val%phase_id(phase_name).eq.0) cycle

      ! Get the unique names for aerosol-phase water
      unique_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = water_name)
      call assert_msg(416554966, allocated(unique_names), &
              "Missing aerosol-phase water species '"//water_name// &
              "' in phase '"//phase_name//"'")
      call assert_msg(344829020, size(unique_names).gt.0, &
              "Missing aerosol-phase water species '"//water_name// &
              "' in phase '"//phase_name//"'")

      ! Count the instances of the specified phase in this aero rep
      ! (One instance per unique water name)
      num_phase = num_phase + size(unique_names)

    end do

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(_NUM_INT_PROP_ + &
            num_phase * (num_spec_per_phase * (num_spec_per_phase + 4) + 2)))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_ + &
            num_spec_per_phase))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set the number of products, reactants and aerosol phase instances
    _NUM_REACT_ = num_react
    _NUM_PROD_ = num_prod
    _NUM_AERO_PHASE_ = num_phase

    ! Get the rate constant parameters
    key_name = "A"
    if (.not. this%property_set%get_real(key_name, _A_)) then
      _A_ = 1.0
    end if
    key_name = "C"
    if (.not. this%property_set%get_real(key_name, _C_)) then
      _C_ = 0.0
    end if
    key_name = "k_reverse"
    call assert_msg(519823533, &
            this%property_set%get_real(key_name, _k_reverse_), &
            "Missing 'k_reverse' for aqueous equilibrium reaction")
    key_name = "time_unit"
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        _k_reverse_ = _k_reverse_ / 60.0
      end if
    end if

    ! Set up an array to the reactant, product and water names
    allocate(react_names(_NUM_REACT_))
    allocate(prod_names(_NUM_PROD_))

    ! Get the chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))

      ! Get the reactant species properties
      call assert_msg(513131584, &
           chem_spec_data%get_property_set(spec_name, spec_props), &
           "Missing properties required for aqueous equilibrium "// &
           "reaction involving species '"//trim(spec_name)//"'")

      ! Get the molecular weight
      key_name = "molecular weight"
      call assert_msg(332898361, spec_props%get_real(key_name, temp_real), &
           "Missing 'molecular weight' for species '"//trim(spec_name)// &
           "' in aqueous equilibrium reaction.")
      
      ! Set properties for each occurance of a reactant in the rxn equation
      call assert(971363961, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (.not.spec_props%get_int(key_name, temp_int)) temp_int = 1
      do i_qty = 1, temp_int
        i_spec = i_spec + 1

        ! Add the reactant name to the list
        react_names(i_spec)%string = spec_name

        ! Use the MW to calculate the mass frac -> M conversion
        _mass_frac_TO_M_(i_spec) = 1.0d0/temp_real

      end do

      ! Go to the next reactant
      call reactants%iter_next()

    end do

    ! Get the chemical properties for the products
    call products%iter_reset()
    i_spec = 0
    do while (products%get_key(spec_name))

      ! Get the product species properties
      call assert_msg(513131584, &
           chem_spec_data%get_property_set(spec_name, spec_props), &
           "Missing properties required for aqueous equilibrium "// &
           "reaction involving species '"//trim(spec_name)//"'")

      ! Get the molecular weight
      key_name = "molecular weight"
      call assert_msg(332898361, spec_props%get_real(key_name, temp_real), &
           "Missing 'molecular weight' for species '"//trim(spec_name)// &
           "' in aqueous equilibrium reaction.")

      ! Set properties for each occurance of a reactant in the rxn equation
      call assert(294785742, products%get_property_t(val=spec_props))
      key_name = "qty"
      if (.not.spec_props%get_int(key_name, temp_int)) temp_int = 1
      do i_qty = 1, temp_int
        i_spec = i_spec + 1

        ! Add the product name to the list
        prod_names(i_spec)%string = spec_name

        ! Use the MW to calculate the mass frac -> M conversion
        _mass_frac_TO_M_(_NUM_REACT_ + i_spec) = 1.0d0/temp_real

      end do

      ! Go to the next product
      call products%iter_next()

    end do

    ! Set the state array indices for the reactants, products and water
    i_aero_phase = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Check for the specified phase in this aero rep
      if (aero_rep(i_aero_rep)%val%phase_id(phase_name).eq.0) cycle

      ! Get the unique names for aerosol-phase water
      unique_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = water_name)

      ! Save the number of instances to check for the presence of 
      ! each product and reactant
      num_phase = size(unique_names)

      ! Save the state ids for aerosol water concentration and activity coeff.
      do i_spec = 1, num_phase
        _WATER_(i_aero_phase + i_spec) = aero_rep(i_aero_rep)%val%spec_state_id( &
                unique_names(i_spec)%string)
        _WATER_ACT_COEFF_(i_aero_phase + i_spec) = &
                aero_rep(i_aero_rep)%val%activity_coeff_state_id( &
                unique_names(i_spec)%string)
      end do

      ! Loop through the reactants
      do i_spec = 1, _NUM_REACT_

        ! Get the unique names for the reactants
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = react_names(i_spec)%string)

        ! Make sure the right number of instances are present
        call assert_msg(280008989, size(unique_names).eq.num_phase, &
                "Incorrect instances of reactant '"//react_names(i_spec)%string// &
                "' in phase '"//phase_name//"' in an aqueous equilibrium reaction")

        ! Save the state ids for the reactant concentration and activity coeff.
        ! IDs are grouped by phase instance: R1(phase1), R2(phase1), ..., R1(phase2)...
        do j_spec = 1, num_phase
          _REACT_((i_aero_phase+j_spec-1)*_NUM_REACT_ + i_spec) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(j_spec)%string)
          _REACT_ACT_COEFF_((i_aero_phase+j_spec-1)*_NUM_REACT_ + i_spec) = &
                  aero_rep(i_aero_rep)%val%activity_coeff_state_id( &
                  unique_names(j_spec)%string)
        end do
      end do

      ! Loop through the products
      do i_spec = 1, _NUM_PROD_

        ! Get the unique names for the products
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = prod_names(i_spec)%string)

        ! Make sure the right number of instances are present
        call assert_msg(557959349, size(unique_names).eq.num_phase, &
                "Incorrect instances of product '"//prod_names(i_spec)%string// &
                "' in phase '"//phase_name//"' in an aqueous equilibrium reaction")

        ! Save the state ids for the product concentration and activity coeff.
        ! IDs are grouped by phase instance: P1(phase1), P2(phase1), ..., P1(phase2)...
        do j_spec = 1, num_phase
          _PROD_((i_aero_phase+j_spec-1)*_NUM_PROD_ + i_spec) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(j_spec)%string)
          _PROD_ACT_COEFF_((i_aero_phase+j_spec-1)*_NUM_PROD_ + i_spec) = &
                  aero_rep(i_aero_rep)%val%activity_coeff_state_id( &
                  unique_names(j_spec)%string)
        end do
      end do

      ! Increment the index offset for the next aerosol representation
      i_aero_phase = i_aero_phase + num_phase

    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rate constant expression
  function build_rate_const_expr(this, rxn_id) result (expr)

    !> Rate constant expression
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_aqueous_equilibrium_t), intent(in) :: this
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
    class(rxn_aqueous_equilibrium_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef _NUM_REACT_
#undef _NUM_PROD_
#undef _NUM_AERO_PHASE_
#undef _A_
#undef _C_
#undef _k_reverse_
#undef _k_forward_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _REACT_
#undef _REACT_ACT_COEFF_
#undef _PROD_
#undef _PROD_ACT_COEFF_
#undef _WATER_
#undef _WATER_ACT_COEFF_
#undef _DERIV_ID_
#undef _JAC_ID_
#undef _mass_frac_TO_M_
end module pmc_rxn_aqueous_equilibrium
