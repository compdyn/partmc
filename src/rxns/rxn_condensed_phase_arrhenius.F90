! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_condensed_phase_arrhenius module.

!> \page camp_rxn_condensed_phase_arrhenius CAMP: Condensed-Phase Arrhenius Reaction
!!
!! Condensed-phase Arrhenius reactions are calculated using an Arrhenius-like
!! rate constant that takes the form:
!!
!! \f[
!!   Ae^{(\frac{-E_a}{k_bT})}(\frac{T}{D})^B(1.0+E*P)
!! \f]
!!
!! where \f$A\f$ is the pre-exponential factor
!! (\f$[\mbox{U}]^{-(n-1)} s^{-1}\f$), \f$U\f$ is the unit of the reactants
!! and products, which can be \f$M\f$ for aqueous-phase reactions or
!! \f$\mbox{\si{\mole\per\cubic\metre}}\f$ for all other condensed-phase
!! reactions, \f$n\f$ is the number of reactants, \f$E_a\f$ is the activation
!! energy (J), \f$k_b\f$ is the Boltzmann constant (J/K), \f$D\f$ (K), \f$B\f$
!! (unitless) and \f$E\f$ (\f$Pa^{-1}\f$) are reaction parameters, \f$T\f$ is
!! the temperature (K), and \f$P\f$ is the pressure (Pa). The first two terms
!! are described in Finlayson-Pitts and Pitts (2000) \cite Finlayson-Pitts2000
!! . The final term is included to accomodate CMAQ EBI solver type 7 rate
!! constants \cite Gipson.
!!
!! Input data for condensed-phase Arrhenius reactions have the following
!! format:
!! \code{.json}
!!   {
!!     "type" : "CONDENSED_PHASE_ARRHENIUS",
!!     "A" : 123.45,
!!     "Ea" : 123.45,
!!     "B"  : 1.3,
!!     "D"  : 300.0,
!!     "E"  : 0.6E-5,
!!     "units" : "M",
!!     "time unit" : "MIN",
!!     "aerosol phase" : "my aqueous phase",
!!     "aerosol-phase water" : "H2O_aq",
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
!! Units for the reactants and products must be specified using the key
!! \b units and can be either \b M or \b mol \b m-3. If units of \b M are
!! specified, a key-value pair \b aerosol-phase \b water must also be included
!! whose value is a string specifying the name for water in the aerosol phase.
!!
!! The unit for time is assumed to be s, but inclusion of the optional
!! key-value pair \b time \b unit = \b MIN can be used to indicate a rate
!! with min as the time unit.
!!
!! The key-value pair \b aerosol \b phase is required and must specify the name
!! of the aerosol-phase in which the reaction occurs.
!!
!! Optionally, a parameter \b C may be included, and is taken to equal
!! \f$\frac{-E_a}{k_b}\f$. Note that either \b Ea or \b C may be included, but
!! not both. When neither \b Ea or \b C are included, they are assumed to be
!! 0.0. When \b A is not included, it is assumed to be 1.0, when \b D is not
!! included, it is assumed to be 300.0 K, when \b B is not included, it is
!! assumed to be 0.0, and when \b E is not included, it is assumed to be 0.0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_condensed_phase_arrhenius_t type and associated functions.
module pmc_rxn_condensed_phase_arrhenius

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_camp_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, &
                                                  die_msg, string_t

  implicit none
  private

#define NUM_REACT_ this%condensed_data_int(1)
#define NUM_PROD_ this%condensed_data_int(2)
#define NUM_AERO_PHASE_ this%condensed_data_int(3)
#define A_ this%condensed_data_real(1)
#define B_ this%condensed_data_real(2)
#define C_ this%condensed_data_real(3)
#define D_ this%condensed_data_real(4)
#define E_ this%condensed_data_real(5)
#define NUM_INT_PROP_ 3
#define NUM_REAL_PROP_ 5
#define NUM_ENV_PARAM_ 1
#define REACT_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define PROD_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_REACT_*NUM_AERO_PHASE_+x)
#define WATER_(x) this%condensed_data_int(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_)*NUM_AERO_PHASE_+x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+1)*NUM_AERO_PHASE_+x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+(2*(NUM_REACT_+NUM_PROD_)+1)*NUM_AERO_PHASE_+x)
#define YIELD_(x) this%condensed_data_real(NUM_REAL_PROP_+x)
#define UGM3_TO_MOLM3_(x) this%condensed_data_real(NUM_REAL_PROP_+NUM_PROD_+x)

  public :: rxn_condensed_phase_arrhenius_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_condensed_phase_arrhenius_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
  end type rxn_condensed_phase_arrhenius_t

  !> Constructor for rxn_condensed_phase_arrhenius_t
  interface rxn_condensed_phase_arrhenius_t
    procedure :: constructor
  end interface rxn_condensed_phase_arrhenius_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for Condensed-phase Arrhenius reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_condensed_phase_arrhenius_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep, n_cells)

    !> Reaction data
    class(rxn_condensed_phase_arrhenius_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)
    !> Number of grid cells to solve simultaneously
    integer(kind=i_kind), intent(in) :: n_cells

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name, water_name, &
            phase_name, temp_string
    integer(kind=i_kind) :: i_spec, i_phase_inst, i_qty, i_aero_rep, &
            i_aero_phase, num_spec_per_phase, num_phase, num_react, num_prod
    type(string_t), allocatable :: unique_names(:)
    type(string_t), allocatable :: react_names(:), prod_names(:)

    integer(kind=i_kind) :: temp_int
    real(kind=dp) :: temp_real
    logical :: is_aqueous

    is_aqueous = .false.

    ! Get the property set
    call assert_msg(852163263, associated(this%property_set), &
            "Missing property set needed to initialize reaction")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(845445717, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in condensed-phase Arrhenius reaction")

    ! Get reactant species
    key_name = "reactants"
    call assert_msg(501773136, &
            this%property_set%get_property_t(key_name, reactants), &
            "Missing reactant species in condensed-phase Arrhenius reaction")

    ! Get product species
    key_name = "products"
    call assert_msg(556252922, &
            this%property_set%get_property_t(key_name, products), &
            "Missing product species in condensed-phase Arrhenius reaction")

    ! Count the number of species per phase instance (including reactants
    ! with a "qty" specified)
    call reactants%iter_reset()
    num_react = 0
    do while (reactants%get_key(spec_name))
      ! Get properties included with this reactant in the reaction data
      call assert(428593951, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (spec_props%get_int(key_name, temp_int)) &
              num_react = num_react + temp_int - 1
      call reactants%iter_next()
      num_react = num_react + 1
    end do
    num_prod = products%size()
    num_spec_per_phase = num_prod + num_react

    ! Check for aerosol representations
    call assert_msg(370755392, associated(aero_rep), &
            "Missing aerosol representation for condensed-phase "// &
            "Arrhenius reaction")
    call assert_msg(483073737, size(aero_rep).gt.0, &
            "Missing aerosol representation for condensed-phase "// &
            "Arrhenius reaction")

    ! Count the instances of the specified aerosol phase
    num_phase = 0
    do i_aero_rep = 1, size(aero_rep)
      num_phase = num_phase + &
              aero_rep(i_aero_rep)%val%num_phase_instances(phase_name)
    end do

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(NUM_INT_PROP_ + &
            num_phase * (num_spec_per_phase * (num_react + 3) + 1)))
    allocate(this%condensed_data_real(NUM_REAL_PROP_ + &
            num_spec_per_phase + num_prod))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Save space for the environment-dependent parameters
    this%num_env_params = NUM_ENV_PARAM_

    ! Set the number of products, reactants and aerosol phase instances
    NUM_REACT_ = num_react
    NUM_PROD_ = num_prod
    NUM_AERO_PHASE_ = num_phase

    ! Get the rate constant parameters
    key_name = "A"
    if (.not. this%property_set%get_real(key_name, A_)) then
      A_ = 1.0
    end if
    key_name = "time unit"
    if (this%property_set%get_string(key_name, temp_string)) then
      if (trim(temp_string).eq."MIN") then
        A_ = A_ / 60.0
      else
        call assert_msg(390870843, trim(temp_string).eq."s", &
                "Received invalid time unit: '"//temp_string//"' in "// &
                "condnesed-phase Arrhenius reaction. Valid units are "// &
                "'MIN' and 's'.")
      end if
    end if
    key_name = "Ea"
    if (this%property_set%get_real(key_name, temp_real)) then
      C_ = -temp_real/const%boltzmann
      key_name = "C"
      call assert_msg(827485736, &
              .not.this%property_set%get_real(key_name, temp_real), &
              "Received both Ea and C parameter for condensed-phase "// &
              "Arrhenius equation.")
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

    ! Set up an array to the reactant and product names
    allocate(react_names(NUM_REACT_))
    allocate(prod_names(NUM_PROD_))

    ! Get the chemical properties for the reactants
    call reactants%iter_reset()
    i_spec = 0
    do while (reactants%get_key(spec_name))

      ! Get the reactant species properties
      call assert_msg(365229559, &
           chem_spec_data%get_property_set(spec_name, spec_props), &
           "Missing properties required for condensed-phase Arrhenius "// &
           "reaction involving species '"//trim(spec_name)//"'")

      ! Get the molecular weight
      key_name = "molecular weight [kg mol-1]"
      call assert_msg(409180731, spec_props%get_real(key_name, temp_real), &
           "Missing 'molecular weight' for species '"//trim(spec_name)// &
           "' in condensed-phase Arrhenius reaction.")

      ! Set properties for each occurance of a reactant in the rxn equation
      call assert(186449575, reactants%get_property_t(val=spec_props))
      key_name = "qty"
      if (.not.spec_props%get_int(key_name, temp_int)) temp_int = 1
      do i_qty = 1, temp_int
        i_spec = i_spec + 1

        ! Add the reactant name to the list
        react_names(i_spec)%string = spec_name

        ! Use the MW to calculate the ug/m3 -> mol/m3 conversion
        UGM3_TO_MOLM3_(i_spec) = 1.0d-9/temp_real

      end do

      ! Go to the next reactant
      call reactants%iter_next()

    end do

    ! Get the chemical properties for the products
    call products%iter_reset()
    i_spec = 0
    do while (products%get_key(spec_name))

      ! Get the product species properties
      call assert_msg(450225425, &
           chem_spec_data%get_property_set(spec_name, spec_props), &
           "Missing properties required for condensed-phase Arrhenius "// &
           "reaction involving species '"//trim(spec_name)//"'")

      ! Increment the product counter
      i_spec = i_spec + 1

      ! Get the molecular weight
      key_name = "molecular weight [kg mol-1]"
      call assert_msg(504705211, spec_props%get_real(key_name, temp_real), &
           "Missing 'molecular weight' for species '"//trim(spec_name)// &
           "' in condensed phase Arrhenius reaction.")

      ! Use the MW to calculate the ug/m3 -> mol/m3 conversion
      UGM3_TO_MOLM3_(NUM_REACT_+i_spec) = 1.0d-9/temp_real

      ! Set properties for each occurance of a reactant in the rxn equation
      call assert(846924553, products%get_property_t(val=spec_props))
      key_name = "yield"
      if (spec_props%get_real(key_name, temp_real)) then
        YIELD_(i_spec) = temp_real
      else
        YIELD_(i_spec) = 1.0d0
      end if

      ! Add the product name to the list
      prod_names(i_spec)%string = spec_name

      ! Go to the next product
      call products%iter_next()

    end do

    ! Get the units for the reactants and products and name for water
    ! if this is an aqueous reaction
    key_name = "units"
    call assert_msg(348722817, &
            this%property_set%get_string(key_name, temp_string), &
            "Missing units for condensed-phase Arrhenius reaction.")
    if (trim(temp_string).eq."mol m-3") then
      is_aqueous = .false.
      key_name = "aerosol-phase water"
      call assert_msg(767767240, &
              .not.this%property_set%get_string(key_name, temp_string), &
              "Aerosol-phase water specified for non-aqueous condensed-"// &
              "phase Arrhenius reaction. Change units to 'M' or remove "// &
              "aerosol-phase water")
    else if (trim(temp_string).eq."M") then
      is_aqueous = .true.
      key_name = "aerosol-phase water"
      call assert_msg(199910264, &
              this%property_set%get_string(key_name, water_name), &
              "Missing aerosol-phase water for aqeuous condensed-phase "// &
              "Arrhenius reaction.")
    else
      call die_msg(161772048, "Received invalid units for condensed-"// &
              "phase Arrhenius reaction: '"//temp_string//"'. Valid "// &
              "units are 'mol m-3' or 'M'.")
    end if

    ! Set the state array indices for the reactants, products and water
    i_aero_phase = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Check for the specified phase in this aero rep
      num_phase = aero_rep(i_aero_rep)%val%num_phase_instances(phase_name)
      if (num_phase.eq.0) cycle

      ! Save the state ids for aerosol-phase water for aqueous reactions.
      ! For non-aqueous reactions, set the water id to -1
      if (is_aqueous) then

        ! Get the unique names for aerosol-phase water
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = water_name)

        ! Make sure water is present
        call assert_msg(196838614, size(unique_names).eq.num_phase, &
                "Missing aerosol-phase water species '"//water_name// &
                "' in phase '"//phase_name//"' in aqueous condensed-"// &
                "phase Arrhenius reacion.")

        ! Save the ids for water in this phase
        do i_phase_inst = 1, num_phase
          WATER_(i_aero_phase + i_phase_inst) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(i_phase_inst)%string)
        end do

        deallocate(unique_names)

      else

        ! Set the water ids to -1 for non-aqueous condensed-phase reactions
        do i_phase_inst = 1, num_phase
          WATER_(i_aero_phase + i_phase_inst) = -1
        end do

      end if

      ! Loop through the reactants
      do i_spec = 1, NUM_REACT_

        ! Get the unique names for the reactants
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = react_names(i_spec)%string)

        ! Make sure the right number of instances are present
        call assert_msg(360730267, size(unique_names).eq.num_phase, &
                "Incorrect instances of reactant '"// &
                react_names(i_spec)%string//"' in phase '"//phase_name// &
                "' in a condensed-phase Arrhenius reaction")

        ! Save the state ids for the reactant concentration
        ! IDs are grouped by phase instance:
        !   R1(phase1), R2(phase1), ..., R1(phase2)...
        do i_phase_inst = 1, num_phase
          REACT_((i_aero_phase+i_phase_inst-1)*NUM_REACT_ + i_spec) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(i_phase_inst)%string)
        end do

        deallocate(unique_names)

      end do

      ! Loop through the products
      do i_spec = 1, NUM_PROD_

        ! Get the unique names for the products
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = prod_names(i_spec)%string)

        ! Make sure the right number of instances are present
        call assert_msg(399869427, size(unique_names).eq.num_phase, &
                "Incorrect instances of product '"// &
                prod_names(i_spec)%string//"' in phase '"//phase_name// &
                "' in a condensed-phase Arrhenius reaction")

        ! Save the state ids for the product concentration
        ! IDs are grouped by phase instance:
        !   P1(phase1), P2(phase1), ..., P1(phase2)...
        do i_phase_inst = 1, num_phase
          PROD_((i_aero_phase+i_phase_inst-1)*NUM_PROD_ + i_spec) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(i_phase_inst)%string)
        end do

        deallocate(unique_names)

      end do

      ! Increment the index offset for the next aerosol representation
      i_aero_phase = i_aero_phase + num_phase

    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_condensed_phase_arrhenius_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_rxn_condensed_phase_arrhenius
