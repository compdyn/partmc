! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_aqueous_equilibrium module.

!> \page phlex_rxn_aqueous_equilibrium Phlexible Module for Chemistry: Aqueous-Phase Equilibrium
!!
!! Aqueous equilibrium reactions are calculated as forward and reverse
!! reactions. The reverse rate constant must be provided as part of the input
!! data and the forward rate constant is calculated using the provided
!! reverse rate constant and an equilibrium constant of the following format:
!!
!! \f[
!!   Ae^{C({1/T-1/298})}
!! \f]
!!
!! where \f$A\f$ is the pre-exponential factor (\f$s^{-1}\f$), \f$C\f$ is a
!! constant and \f$T\f$ is the temperature (K).
!!
!! Input data for aqueous equilibrium equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "AQUEOUS_EQUILIBRIUM",
!!     "A" : 123.45,
!!     "C" : 123.45,
!!     "k_reverse" : 123.45,
!!     "phase" : "my aqueous phase",
!!     "time unit" : "MIN",
!!     "aqueous-phase water" : "H2O_aq",
!!     "ion pair" : "spec3-spec4",
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
!! specified phase and include a \b molecular \b weight parameter. The
!! parameter \b aqueous-phase \b water is required and must be the name of the
!! aerosol-phase species that is used for water. The parameter \b ion \b pair
!! is optional. When it is included its value must be the name of an ion pair
!! that is present in the specified aerosol phase. Its mean binary activity
!! coefficient will be applied to the reverse reaction.
!!
!! When \b A is not included, it is assumed to be 1.0, when \b C is not
!! included, it is assumed to be 0.0. The reverse reaction rate constant
!! \b k_reverse is required.
!!
!! The unit for time is assumed to be s, but inclusion of the optional 
!! key-value pair \b time \b unit = \b MIN can be used to indicate a rate with
!! min as the time unit.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_aqueous_equilibrium_t type and associated functions. 
module pmc_rxn_aqueous_equilibrium

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_phlex_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: phlex_real, phlex_int, &
                                                  to_string, assert, &
                                                  assert_msg, die_msg, &
                                                  string_t, align_ratio

  implicit none
  private

#define NUM_REACT_ this%condensed_data_int(1)
#define NUM_PROD_ this%condensed_data_int(2)
#define NUM_AERO_PHASE_ this%condensed_data_int(3)
#define INT_DATA_SIZE_ this%condensed_data_int(4)
#define FLOAT_DATA_SIZE_ this%condensed_data_int(5)
#define A_ this%condensed_data_real(1)
#define C_ this%condensed_data_real(2)
#define RATE_CONST_REVERSE_ this%condensed_data_real(3)
#define RATE_CONST_FORWARD_ this%condensed_data_real(4)
#define NUM_INT_PROP_ 5
#define NUM_REAL_PROP_ 4
#define REACT_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define PROD_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_REACT_*NUM_AERO_PHASE_+x)
#define WATER_(x) this%condensed_data_int(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_)*NUM_AERO_PHASE_+x)
#define ACTIVITY_COEFF_(x) this%condensed_data_int(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+1)*NUM_AERO_PHASE_+x)
#define DERIV_ID_(x) this%condensed_data_int(NUM_INT_PROP_+(NUM_REACT_+NUM_PROD_+2)*NUM_AERO_PHASE_+x)
#define JAC_ID_(x) this%condensed_data_int(NUM_INT_PROP_+(2*(NUM_REACT_+NUM_PROD_)+2)*NUM_AERO_PHASE_+x)
#define MASS_FRAC_TO_M_(x) this%condensed_data_real(NUM_REAL_PROP_+x)

  public :: rxn_aqueous_equilibrium_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_aqueous_equilibrium_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
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
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep)
    
    !> Reaction data
    class(rxn_aqueous_equilibrium_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props, reactants, products
    character(len=:), allocatable :: key_name, spec_name, water_name, &
            phase_name, string_val, ion_pair_name
    integer(kind=phlex_int) :: i_phase_inst, j_spec, i_qty, i_aero_rep, &
            i_aero_phase, num_spec_per_phase, num_phase, num_react, &
            num_prod, temp_int, tracer_type
    real(kind=phlex_real) :: temp_real
    type(string_t), allocatable :: unique_names(:), react_names(:), &
            prod_names(:)
    integer(kind=phlex_int) :: int_data_size, float_data_size

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
    call reactants%iter_reset()
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
    call products%iter_reset()
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
      if (aero_rep(i_aero_rep)%val%num_phase_instances(phase_name).eq.0) cycle

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

      deallocate(unique_names)

    end do

    ! Calculate int and float array sizes with alignment spacing
    int_data_size = NUM_INT_PROP_ + &
            num_phase * (num_spec_per_phase * (num_spec_per_phase + 3) + 2)
    int_data_size = int_data_size + mod(int_data_size, align_ratio)
    float_data_size = NUM_REAL_PROP_ + num_spec_per_phase

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(int_data_size))
    allocate(this%condensed_data_real(float_data_size))
    this%condensed_data_int(:) = int(0, kind=phlex_int)
    this%condensed_data_real(:) = real(0.0, kind=phlex_real)
    INT_DATA_SIZE_ = int_data_size
    FLOAT_DATA_SIZE_ = float_data_size

    ! Set the number of products, reactants and aerosol phase instances
    NUM_REACT_ = num_react
    NUM_PROD_ = num_prod
    NUM_AERO_PHASE_ = num_phase

    ! Get the rate constant parameters
    key_name = "A"
    if (.not. this%property_set%get_real(key_name, A_)) then
      A_ = 1.0
    end if
    key_name = "C"
    if (.not. this%property_set%get_real(key_name, C_)) then
      C_ = 0.0
    end if
    key_name = "k_reverse"
    call assert_msg(519823533, &
            this%property_set%get_real(key_name, RATE_CONST_REVERSE_), &
            "Missing 'k_reverse' for aqueous equilibrium reaction")
    key_name = "time_unit"
    if (this%property_set%get_string(key_name, string_val)) then
      if (trim(string_val).eq."MIN") then
        RATE_CONST_REVERSE_ = RATE_CONST_REVERSE_ / 60.0
      end if
    end if

    ! Set up an array to the reactant, product and water names
    allocate(react_names(NUM_REACT_))
    allocate(prod_names(NUM_PROD_))

    ! Get the chemical properties for the reactants
    call reactants%iter_reset()
    i_phase_inst = 0
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
        i_phase_inst = i_phase_inst + 1

        ! Add the reactant name to the list
        react_names(i_phase_inst)%string = spec_name

        ! Use the MW to calculate the mass frac -> M conversion
        MASS_FRAC_TO_M_(i_phase_inst) = 1.0d3/temp_real

      end do

      ! Go to the next reactant
      call reactants%iter_next()

    end do

    ! Get the chemical properties for the products
    call products%iter_reset()
    i_phase_inst = 0
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
        i_phase_inst = i_phase_inst + 1

        ! Add the product name to the list
        prod_names(i_phase_inst)%string = spec_name

        ! Use the MW to calculate the mass frac -> M conversion
        MASS_FRAC_TO_M_(NUM_REACT_ + i_phase_inst) = 1.0d3/temp_real

      end do

      ! Go to the next product
      call products%iter_next()

    end do

    ! Check for an ion pair to use for the activity coefficient of the reverse
    ! reaction
    key_name = "ion pair"
    if (.not. this%property_set%get_string(key_name, ion_pair_name)) then
        ion_pair_name = ""
    end if

    ! Set the state array indices for the reactants, products and water
    i_aero_phase = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Check for the specified phase in this aero rep
      if (aero_rep(i_aero_rep)%val%num_phase_instances(phase_name).eq.0) cycle

      ! Get the unique names for aerosol-phase water
      unique_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = water_name)

      ! Save the number of instances to check for the presence of 
      ! each product and reactant
      num_phase = size(unique_names)

      ! Save the state ids for aerosol water concentration
      do i_phase_inst = 1, num_phase
        WATER_(i_aero_phase + i_phase_inst) = &
                aero_rep(i_aero_rep)%val%spec_state_id( &
                unique_names(i_phase_inst)%string)
      end do

      deallocate(unique_names)

      ! If an ion pair was specified, save its id for each phase instance
      if (ion_pair_name.ne."") then
     
        ! Get the unique names for the ion pair
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = ion_pair_name)

        ! Make sure the right number of instances are present
        call assert_msg(298310186, size(unique_names).eq.num_phase, &
                "Incorrect instances of ion pair '"//ion_pair_name// &
                "' in phase '"//phase_name// &
                "' in an aqueous equilibrium reaction")

        ! Make sure the specified ion pair is of the right tracer type
        call assert(690720371, &
                chem_spec_data%get_type(ion_pair_name, tracer_type))
        call assert_msg(450743055, tracer_type.eq.CHEM_SPEC_ACTIVITY_COEFF, &
                "Ion pair '"//ion_pair_name//"' must have type "// &
                "'ION_PAIR' to be used as an ion pair in an aqueous "// &
                "equilibrium reaction.")

        ! Save the ion pair id
        do i_phase_inst = 1, num_phase
          ACTIVITY_COEFF_(i_aero_phase + i_phase_inst) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(i_phase_inst)%string)
        end do

        deallocate(unique_names)

      else
        do i_phase_inst = 1, num_phase
          ACTIVITY_COEFF_(i_aero_phase + i_phase_inst) = 0
        end do
      end if

      ! Loop through the reactants
      do i_phase_inst = 1, NUM_REACT_

        ! Get the unique names for the reactants
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, &
                spec_name = react_names(i_phase_inst)%string)

        ! Make sure the right number of instances are present
        call assert_msg(280008989, size(unique_names).eq.num_phase, &
                "Incorrect instances of reactant '"// &
                react_names(i_phase_inst)%string// &
                "' in phase '"//phase_name// &
                "' in an aqueous equilibrium reaction")

        ! Save the state ids for the reactant concentration
        ! IDs are grouped by phase instance: 
        !    R1(phase1), R2(phase1), ..., R1(phase2)...
        do j_spec = 1, num_phase
          REACT_((i_aero_phase+j_spec-1)*NUM_REACT_ + i_phase_inst) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(j_spec)%string)
        end do
        
        deallocate(unique_names)

      end do

      ! Loop through the products
      do i_phase_inst = 1, NUM_PROD_

        ! Get the unique names for the products
        unique_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, &
                spec_name = prod_names(i_phase_inst)%string)

        ! Make sure the right number of instances are present
        call assert_msg(557959349, size(unique_names).eq.num_phase, &
                "Incorrect instances of product '"// &
                prod_names(i_phase_inst)%string// &
                "' in phase '"//phase_name// &
                "' in an aqueous equilibrium reaction")

        ! Save the state ids for the product concentration
        ! IDs are grouped by phase instance:
        !    P1(phase1), P2(phase1), ..., P1(phase2)...
        do j_spec = 1, num_phase
          PROD_((i_aero_phase+j_spec-1)*NUM_PROD_ + i_phase_inst) = &
                  aero_rep(i_aero_rep)%val%spec_state_id( &
                  unique_names(j_spec)%string)
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
    type(rxn_aqueous_equilibrium_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef NUM_REACT_
#undef NUM_PROD_
#undef NUM_AERO_PHASE_
#undef A_
#undef C_
#undef RATE_CONST_REVERSE_
#undef RATE_CONST_FORWARD_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef REACT_
#undef PROD_
#undef WATER_
#undef ACTIVITY_COEFF_
#undef DERIV_ID_
#undef JAC_ID_
#undef MASS_FRAC_TO_M_
end module pmc_rxn_aqueous_equilibrium
