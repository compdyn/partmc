! Copyright (C) 2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_sub_model_PDFiTE module.

! FIXME Move to a sub model

!> \page phlex_sub_model_PDFiTE Phlexible Module for Chemistry: PD-FiTE Activity
!!
!! PD-FiTE activity calculates aerosol-phase species activities using
!! Taylor series to describe partial derivatives of mean activity coefficients
!! for ternary solutions, as described in Topping et al. (2009)
!! \cite Topping2009 . The mean binary activity coefficients for ion pairs are
!! calculated according to eq. 15 in \cite Topping2009 . The values are then
!! made available to aqueous-phase reactions during solving.
!!
!! Input data for PDFiTE activity equations have the following format :
!! \code{.json}
!!  { "pmc-data" : [
!!   {
!!     "name" : "my PDFiTE activity",
!!     "type" : "SUB_MODEL_PDFITE",
!!     "gas-phase water" : "H2O",
!!     "aerosol-phase water" : "H2O_aq",
!!     "aerosol phase" : "my aero phase",
!!     "calculate for" : {
!!       "H-NO3" : {
!!         "interactions" : [
!!           {
!!             "ion pair" : "H-NO3",
!!             "min RH" : 0.0,
!!             "max RH" : 0.1,
!!             "B" : [ 0.925113 ]
!!           },
!!           {
!!             "ion pair" : "H-NO3",
!!             "min RH" : 0.1,
!!             "max RH" : 0.4,
!!             "B" : [ 0.12091, 13.497, -67.771, 144.01, -117.97 ]
!!           },
!!           {
!!             "ion pair" : "H-NO3",
!!             "min RH" : 0.4,
!!             "max RH" : 0.9,
!!             "B" : [ 1.3424, -0.8197, -0.52983, -0.37335 ]
!!           },
!!           {
!!             "ion pair" : "H-NO3",
!!             "min RH" : 0.9,
!!             "max RH" : 1.0,
!!             "B" : [ -0.3506505 ]
!!           },
!!           {
!!             "ion pair" : "NH4-NO3",
!!             "min RH" : 0.0,
!!             "max RH" : 0.1,
!!             "B" : [ -11.93308 ]
!!           },
!!           {
!!             "ion pair" : "NH4-NO3",
!!             "min RH" : 0.1,
!!             "max RH" : 0.99,
!!             "B" : [ -17.0372, 59.232, -86.312, 44.04 ]
!!           },
!!           {
!!             "ion pair" : "NH4-NO3",
!!             "min RH" : 0.99,
!!             "max RH" : 1.0,
!!             "B" : [ -0.2599432 ]
!!           }
!!         ]
!!       }
!!       ...
!!     }
!!   }
!!  ]}
!! \endcode
!! The key-value pair \b aerosol \b phase is required to specify the aerosol
!! phase for which to calculate activity coefficients. The key-value pairs
!! \b gas-phase \b water and \b aerosol-phase \b water must also be present
!! and specify the names for the water species in each phase. The final
!! required key-value pair is \b calculated \b for, which should contain a set
!! of ion pairs that activity coefficients will be calculated for.
!!
!! The key names in this set must correspond to ion pairs that are present in
!! the specified aerosol phase. The values must contain a key-value pair
!! named \b interactions which includes an array of ion-pair interactions
!! used to calculate equation 15 in \cite Topping2009 .
!!
!! Each element in the \b interactions array must include an \b ion pair
!! that exists in the specified aerosol phase, a \b min \b RH and \b max \b RH
!! that specify the bounds for which the fitted curve is valid, and an array
!! of \b B values that specify the polynomial coefficients B0, B1, B2, ...
!! as in equation 19 in \cite Topping2009 . At least one polynomial
!! coefficient must be present.
!!
!! If at least one interaction with an ion pair is included, enough
!! interactions with that ion pair must be included to cover the entire RH
!! range (0.0--1.0). Interactions are assume to cover the range (minRH, maxRH],
!! except for the lowest RH interaction, which covers th range [0.0, maxRH].
!!
!! When the interacting ion pair is the same as the ion-pair for which the
!! mean binary activity coefficient is being calculated, the interaction
!! parameters are used to calculate \f$ln(\gamma_A^0(RH))\f$. Otherwise, the
!! parameters are used to calculate
!! \f$\frac{dln(gamma_A))}{d(N_{B,M}N_{B,x})}\f$.
!!
!! For the above example, the following input data should be present:
!! \code{.json}
!! {
!!   "name" : "H2O",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "GAS",
!! },
!! {
!!   "name" : "H2O_aq",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!! },
!! {
!!   "name" : "H_p",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : 1,
!!   "molecular weight [kg mol-1]" : 1.008
!! },
!! {
!!   "name" : "NH4_p",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : 1,
!!   "molecular weight" : 18.04
!! },
!! {
!!   "name" : "NO3_m",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : -1
!!   "molecular weight [kg mol-1]" : 62.0049
!! },
!! {
!!   "name" : "NH4-NO3",
!!   "type" : "CHEM_SPEC",
!!   "tracer type" : "ION_PAIR",
!!   "ions" : {
!!     "NH4_p" : {},
!!     "NO3_m" : {}
!!   }
!! },
!! {
!!   "name" : "H-NO3",
!!   "type" : "CHEM_SPEC",
!!   "tracer type" : "ION_PAIR",
!!   "ions" : {
!!     "H_p" : {},
!!     "NO3_m" : {}
!!   }
!! },
!! {
!!   "name" : "my aero phase",
!!   "type" : "AERO_PHASE",
!!   "species" : ["H_p", "NO3_m", "NH4_p", "NH4-NO3", "H-NO3", "H2O_aq"]
!! }
!! \endcode
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The sub_model_PDFiTE_t type and associated functions.
module pmc_sub_model_PDFiTE

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_phlex_state
  use pmc_property
  use pmc_sub_model_data
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, &
                                                  die_msg, string_t, &
                                                  warn_assert_msg

  implicit none
  private

#define NUM_PHASE_ this%condensed_data_int(1)
#define GAS_WATER_ID_ this%condensed_data_int(2)
#define NUM_ION_PAIRS_ this%condensed_data_int(3)
#define TOTAL_INT_PARAM_ this%condensed_data_int(4)
#define TOTAL_FLOAT_PARAM_ this%condensed_data_int(5)
#define PPM_TO_RH_ this%condensed_data_real(1)
#define NUM_INT_PROP_ 5
#define NUM_REAL_PROP_ 1
#define PHASE_ID_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define PAIR_INT_PARAM_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_PHASE_+x)
#define PAIR_FLOAT_PARAM_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_PHASE_+NUM_ION_PAIRS_+x)
#define ION_PAIR_ACT_ID_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x))
#define NUM_CATION_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+1)
#define NUM_ANION_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+2)
#define CATION_ID_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+3)
#define ANION_ID_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+4)
#define NUM_INTER_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+5)
#define NUM_B_(x,y) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+5+y)
#define INTER_SPEC_ID_(x,y) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+5+NUM_INTER_(x)+y)
#define INTER_SPEC_LOC_(x,y) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+5+2*(NUM_INTER_(x))+y)
#define CATION_MW_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x))
#define ANION_MW_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+1)
#define CATION_N_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+2)
#define ANION_N_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+3)
#define MIN_RH_(x,y) this%condensed_data_real(INTER_SPEC_LOC_(x,y))
#define MAX_RH_(x,y) this%condensed_data_real(INTER_SPEC_LOC_(x,y)+1)
#define B_Z_(x,y,z) this%condensed_data_real(INTER_SPEC_LOC_(x,y)+1+z)

  ! Update types (These must match values in sub_model_UNIFAC.c)
  ! (none for now)

  public :: sub_model_PDFiTE_t

  !> Generic test reaction data type
  type, extends(sub_model_data_t) :: sub_model_PDFiTE_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize the reaction
    final :: finalize
  end type sub_model_PDFiTE_t

  !> Constructor for sub_model_PDFiTE_t
  interface sub_model_PDFiTE_t
    procedure :: constructor
  end interface sub_model_PDFiTE_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for activity reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(sub_model_PDFiTE_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, aero_rep_set, aero_phase_set, chem_spec_data)

    !> Reaction data
    class(sub_model_PDFiTE_t), intent(inout) :: this
    !> The set of aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep_set(:)
    !> The set of aerosol phases
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    type(property_t), pointer :: spec_props, ion_pairs, ion_pair, &
            sub_props, ions, interactions, interaction, poly_coeffs
    character(len=:), allocatable :: key_name, spec_name, phase_name, &
            string_val, inter_spec_name
    integer(kind=i_kind) :: n_phase, n_ion_pair, n_int_param, n_float_param
    integer(kind=i_kind) :: i_aero_rep, i_phase, i_ion_pair, i_ion, i_spec, &
            i_poly_coeff, i_interaction, j_ion_pair, j_interaction
    integer(kind=i_kind) :: qty, int_val, charge, total_charge, tracer_type
    real(kind=dp) :: real_val, molecular_weight, min_RH, max_RH
    type(string_t), allocatable :: unique_spec_names(:)
    character(len=:), allocatable :: ion_pair_name, ion_name
    type(string_t), allocatable :: ion_pair_names(:), temp_ion_pair_names(:)
    integer(kind=i_kind), allocatable :: num_inter(:)
    real(kind=dp), allocatable :: rh_range(:)

    ! Get the reaction property set
    if (.not. associated(this%property_set)) call die_msg(101529793, &
            "Missing property set needed to initialize PDFiTE activity"// &
            " reaction.")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(429861748, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in PDFiTE activity reaction.")

    ! Count the instances of the aerosol phase
    n_phase = 0
    do i_aero_rep = 1, size(aero_rep_set)

      ! Get the number of instances of the phase in this representation
      n_phase = n_phase + &
              aero_rep_set(i_aero_rep)%val%num_phase_instances(phase_name)

    end do

    call assert_msg(656403972, n_phase.gt.0, &
            "Aerosol phase '"//phase_name//"' not present in any aerosol "// &
            "representation for PDFiTE activity reaction.")

    ! Get the ion pairs to calculate mean binary activity coefficients for
    key_name = "calculate for"
    call assert_msg(258251605, &
            this%property_set%get_property_t(key_name, ion_pairs), &
            "Missing ion pairs to calculate activity for in PDFiTE "// &
            " activity reaction.")

    ! Count the ion pairs to calculate activity coefficients for
    n_ion_pair = ion_pairs%size()
    call assert_msg(479529522, n_ion_pair .gt. 0, &
            "Empty ion pair set in PDFiTE activity reaction.")

    ! Allocate space for the ion pair names and array to check for their
    ! interaction parameters
    allocate(ion_pair_names(n_ion_pair))

    ! Set the ion pair names
    call ion_pairs%iter_reset()
    do i_ion_pair = 1, n_ion_pair

      ! Get the name of the ion_pair
      call assert(497564051, ion_pairs%get_key(ion_pair_name))
      ion_pair_names(i_ion_pair)%string = ion_pair_name

      call ion_pairs%iter_next()
    end do

    ! Get the number of parameters required for each ion pair
    ! Adding space for INT_PROPS and phase state ids
    n_int_param = NUM_INT_PROP_ + n_phase
    ! Adding space for REAL_PROPS
    n_float_param = NUM_REAL_PROP_
    call ion_pairs%iter_reset()
    do i_ion_pair = 1, n_ion_pair

      ! Get the name of the ion_pair
      call assert(680654801, ion_pairs%get_key(ion_pair_name))

      ! Get the ion_pair properties
      call assert_msg(287766741, ion_pairs%get_property_t(val=ion_pair), &
              "Missing ion pair properties for '"//ion_pair_name// &
              "' in PDFiTE activity reaction.")

      ! Get the interactions
      key_name = "interactions"
      call assert_msg(883233319, &
              ion_pair%get_property_t(key_name, interactions), &
              "Missing interaction parameters for '"//ion_pair_name// &
              "' in PDFiTE activity reaction.")

      ! Loop through the interaction set and count the parameters for each
      call interactions%iter_reset()
      do i_interaction = 1, interactions%size()

        ! Get the current interaction
        call assert(347348920, interactions%get_property_t(val=interaction))

        ! Get the name of the interacting ion pair for use in error messages
        key_name = "ion pair"
        call assert_msg(924945098, &
                interaction%get_string(key_name, inter_spec_name), &
                "Missing interaction species name for ion pair '"// &
                ion_pair_name//"' in PD-FiTE activity reaction")

        ! Get the set of polynomial coefficients
        key_name = "B"
        call assert_msg(157416082, &
                interaction%get_property_t(key_name, poly_coeffs), &
                "Missing 'B' array of polynomial coefficients for ion "// &
                "pair '"//inter_spec_name//"' in activity coefficient "// &
                "calculation for ion pair '"//ion_pair_name//"' in "// &
                "PD-FiTE activity reaction.")

        ! Check the number of polynomial coefficients
        call assert_msg(320927773, poly_coeffs%size().gt.0, &
                "Insufficient polynomial coefficients for PDFiTE activity "//&
                "calculation for ion_pair '"//ion_pair_name//"'.")

        ! Adding space for minRH, maxRH, B[]
        n_float_param = n_float_param + 2 + poly_coeffs%size()

        ! Adding space for size(B[]), interacting species id, location of
        ! interaction parameters
        n_int_param = n_int_param + 3

        ! Make sure that this ion pair is in the list of ion pairs. If not,
        ! add it to the list.
        do j_ion_pair = 1, size(ion_pair_names)
          if (inter_spec_name.eq.ion_pair_names(j_ion_pair)%string) exit
          if (j_ion_pair.eq.size(ion_pair_names)) then
            allocate(temp_ion_pair_names(j_ion_pair))
            forall (i_spec=1:j_ion_pair)
              temp_ion_pair_names(i_spec)%string = &
                    ion_pair_names(i_spec)%string
            end forall
            deallocate(ion_pair_names)
            allocate(ion_pair_names(j_ion_pair+1))
            forall (i_spec=1:j_ion_pair)
              ion_pair_names(i_spec)%string = &
                    temp_ion_pair_names(i_spec)%string
            end forall
            deallocate(temp_ion_pair_names)
            ion_pair_names(j_ion_pair+1)%string = inter_spec_name
          end if
        end do

        call interactions%iter_next()
      end do

      call ion_pairs%iter_next()
    end do

    ! Update the number of ion pairs to include those from interactions that
    ! activity coefficients are not calculated for
    n_ion_pair = size(ion_pair_names)

    ! Adding space for locations of int and float data for each ion pair,
    ! ion pair activity coefficient state id, number of cations and anions,
    ! state ids of cation and anion, and the number of interactions
    n_int_param = n_int_param + 8*n_ion_pair
    ! Adding space for MW of the cation and anion, and cation and anion
    ! concentrations (for use during solving)
    n_float_param = n_float_param + 4*n_ion_pair

    ! Allocate space for an array used to make sure all ion pair interactions
    ! are included and the RH range is covered for each
    allocate(num_inter(n_ion_pair))
    allocate(rh_range(n_ion_pair))

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(n_int_param))
    allocate(this%condensed_data_real(n_float_param))
    this%condensed_data_int(:) = int(9999, kind=i_kind)
    this%condensed_data_real(:) = real(9999.0, kind=dp)

    ! Set some data dimensions
    NUM_PHASE_  = n_phase
    NUM_ION_PAIRS_ = n_ion_pair
    TOTAL_INT_PARAM_ = n_int_param
    TOTAL_FLOAT_PARAM_ = n_float_param

    ! Set the gas-phase water species
    key_name = "gas-phase water"
    call assert_msg(901506020, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing gas-phase water species name in PDFiTE activity "// &
            "reaction.")

    GAS_WATER_ID_ = chem_spec_data%gas_state_id(spec_name)

    call assert_msg(442608616, GAS_WATER_ID_ .gt. 0, &
            "Cannot find gas-phase water species '"//spec_name//"' for "// &
            "PDFiTE activity reaction.")

    ! Set the aerosol-water species
    key_name = "aerosol-phase water"
    call assert_msg(214613153, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing aerosol-phase water species name in PDFiTE activity "// &
            "reaction.")

    ! Make the PHASE_ID_(x) hold the state id of aerosol water in each
    ! phase instance. Then the aerosol water id is 1, and the ion
    ! ids will be relative to the water id in each phase.
    i_phase = 1
    do i_aero_rep = 1, size(aero_rep_set)
      unique_spec_names = aero_rep_set(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)
      if (.not.allocated(unique_spec_names)) cycle
      do i_spec = 1, size(unique_spec_names)
        PHASE_ID_(i_phase) = aero_rep_set(i_aero_rep)%val%spec_state_id( &
                unique_spec_names(i_spec)%string)
        call assert(658622226, PHASE_ID_(i_phase).gt.0)
        i_phase = i_phase + 1
      end do
      deallocate(unique_spec_names)
    end do
    i_phase = i_phase - 1
    call assert_msg(653357919, i_phase.eq.NUM_PHASE_, &
            "Incorrect number of aerosol water instances in PDFiTE "// &
            "activity reaction. Expected "//trim(to_string(NUM_PHASE_))// &
            " but got "//trim(to_string(i_phase)))

    ! Save the ion_pair parameters
    n_int_param = NUM_INT_PROP_ + NUM_PHASE_ + 2*NUM_ION_PAIRS_
    n_float_param = NUM_REAL_PROP_
    call ion_pairs%iter_reset()
    do i_ion_pair = 1, n_ion_pair

      ! Get the name of the ion_pair
      ion_pair_name = ion_pair_names(i_ion_pair)%string

      ! Set the location of this ion_pair's parameters in the condensed data
      ! arrays.
      PAIR_INT_PARAM_LOC_(i_ion_pair) = n_int_param + 1
      PAIR_FLOAT_PARAM_LOC_(i_ion_pair) = n_float_param + 1

      ! Get the activity coefficient state ids for this ion_pair
      i_phase = 1
      do i_aero_rep = 1, size(aero_rep_set)
        unique_spec_names = aero_rep_set(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = ion_pair_name)
        if (.not.allocated(unique_spec_names)) cycle
        do i_spec = 1, size(unique_spec_names)
          if (i_phase.eq.1) then
            ION_PAIR_ACT_ID_(i_ion_pair) = &
                    aero_rep_set(i_aero_rep)%val%spec_state_id( &
                    unique_spec_names(i_spec)%string) - &
                    PHASE_ID_(i_phase)
          else
            call assert(142173386, ION_PAIR_ACT_ID_(i_ion_pair).eq. &
                    aero_rep_set(i_aero_rep)%val%spec_state_id( &
                    unique_spec_names(i_spec)%string) - &
                    PHASE_ID_(i_phase))
          end if
          i_phase = i_phase + 1
        end do
        deallocate(unique_spec_names)
      end do
      i_phase = i_phase - 1
      call assert_msg(700406338, i_phase.eq.NUM_PHASE_, &
              "Incorrect number of instances of ion pair '"// &
              ion_pair_name//"' in PDFiTE activity reaction. Expected "// &
              trim(to_string(NUM_PHASE_))//" but got "// &
              trim(to_string(i_phase)))

      ! Get the ion_pair species properties
      call assert_msg(737691158, &
              chem_spec_data%get_property_set(ion_pair_name, ion_pair), &
              "Missing species properties for ion pair '"// &
              ion_pair_name//"' in PDFiTE activity reaction.")

      ! Make sure the specified ion pair is of the right tracer type
      call assert(667033398, &
              chem_spec_data%get_type(ion_pair_name, tracer_type))
      call assert_msg(153203913, tracer_type.eq.CHEM_SPEC_ACTIVITY_COEFF, &
              "Ion pair '"//ion_pair_name//"' must have a tracer type '"// &
              "ION_PAIR' to participate in a PD-FiTE activity reaction.")

      ! Get the number and id of the ions and their molecular weights
      key_name = "ions"
      call assert_msg(974596824, ion_pair%get_property_t(key_name, ions), &
              "Mission ions for ion pair '"//ion_pair_name// &
              "' in PDFiTE activity reaction.")

      call assert_msg(852202160, ions%size().eq.2, &
              "Invalid number of unique ions specified for ion pair '"// &
              ion_pair_name//"' in for PDFiTE in activity reaction. "// &
              "Expected 2 got "//trim(to_string(ions%size())))

      ! TODO Consider moving some of the ion data validation to
      ! pmc_chem_spec_data
      call ions%iter_reset()
      total_charge = 0
      do i_ion = 1, 2

        ! Get the ion name
        call assert(622753458, ions%get_key(ion_name))

        ! Get the qty, if specified
        qty = 1
        if (ions%get_property_t(val=sub_props)) then
          key_name = "qty"
          if (sub_props%get_int(key_name, int_val)) qty = int_val
        end if

        ! Get the species properties
        call assert_msg(619394685, &
                chem_spec_data%get_property_set(ion_name, spec_props), &
                "Missing species properties for ion '"//ion_name// &
                "' in PDFiTE activity reaction.")

        ! Get the molecular weight
        key_name = "molecular weight [kg mol-1]"
        call assert_msg(951085413, &
                spec_props%get_real(key_name, molecular_weight), &
                "Missing molecular weight for ion '"//ion_name// &
                "' in PDFiTE activity reaction.")

        ! Add the charge from this species
        key_name = "charge"
        call assert_msg(663798152, spec_props%get_int(key_name, charge), &
                "Missing charge for ion '"//ion_name//"' in PDFiTE "// &
                "activity reaction.")

        if (charge.gt.0) then
          NUM_CATION_(i_ion_pair) = qty
          CATION_MW_(i_ion_pair) = molecular_weight
        else if (charge.lt.0) then
          NUM_ANION_(i_ion_pair) = qty
          ANION_MW_(i_ion_pair) = molecular_weight
        else
          call die_msg(939555855, "Neutral species '"//ion_name// &
                  "' not allowed in PDFiTE activity reaction ion pair")
        end if

        ! Add contribution to total charge
        total_charge = total_charge + qty * charge

        ! Get the state ids for this species
        i_phase = 1
        do i_aero_rep = 1, size(aero_rep_set)
          unique_spec_names = aero_rep_set(i_aero_rep)%val%unique_names( &
                  phase_name = phase_name, spec_name = ion_name)
          if (.not.allocated(unique_spec_names)) cycle
          do i_spec = 1, size(unique_spec_names)
            if (charge.gt.0) then
              if (i_phase.eq.1) then
                CATION_ID_(i_ion_pair) = &
                        aero_rep_set(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        PHASE_ID_(i_phase)
              else
                call assert(425726370, CATION_ID_(i_ion_pair).eq. &
                        aero_rep_set(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        PHASE_ID_(i_phase))
              end if
            else
              if (i_phase.eq.1) then
                ANION_ID_(i_ion_pair) = &
                        aero_rep_set(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        PHASE_ID_(i_phase)
              else
                call assert(192466600, ANION_ID_(i_ion_pair).eq. &
                        aero_rep_set(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        PHASE_ID_(i_phase))
              end if
            end if
            i_phase = i_phase + 1
          end do
          deallocate(unique_spec_names)
        end do
        i_phase = i_phase - 1
        call assert_msg(759322632, i_phase.eq.NUM_PHASE_, &
                "Incorrect number of instances of ion species '"// &
                ion_name//"' in PDFiTE activity reaction. Expected "// &
                trim(to_string(NUM_PHASE_))//" but got "// &
                trim(to_string(i_phase)))

        ! Get the next ion
        call ions%iter_next()

      end do

      call assert_msg(415650051, total_charge.eq.0, &
              "Charge imbalance for ion_pair '"//ion_pair_name// &
              " in PDFiTE activity reaction. Total charge: "// &
              trim(to_string(total_charge)))

      n_float_param = n_float_param + 4
      n_int_param = n_int_param + 6

      ! The first portion of ion_pair_names holds species for which
      ! activity coefficients are calculated
      if (i_ion_pair.le.ion_pairs%size()) then
        call assert(362061689, ion_pairs%get_key(string_val))
        call assert(444603372, string_val.eq.ion_pair_name)

        ! Get the ion_pair properties
        call assert(520886936, ion_pairs%get_property_t(val=ion_pair))

        ! Get the interactions
        key_name = "interactions"
        call assert(216229321, ion_pair%get_property_t(key_name,interactions))

        ! Set the number of interactions for this ion pair
        NUM_INTER_(i_ion_pair) = interactions%size()

        ! Get the interaction parameters
        num_inter(:) = 0
        call interactions%iter_reset()
        do i_interaction = 1, interactions%size()

          ! Set the location of the interaction float parameters
          INTER_SPEC_LOC_(i_ion_pair, i_interaction) = n_float_param + 1

          ! Get the current interaction
          call assert(105466070, interactions%get_property_t(val=interaction))

          ! Get the name of the interacting species
          key_name = "ion pair"
          call assert_msg(220815620, &
                  interaction%get_string(key_name, inter_spec_name), &
                  "Missing interacting species name for ion pair '"// &
                  ion_pair_name//"' in PDFiTE activity reaction.")

          ! Mark the species as having an interaction and get the id of the
          ! ion_pair in the 'calculate for' list
          do i_spec = 1, size(ion_pair_names)
            if (ion_pair_names(i_spec)%string.eq.inter_spec_name) then
              num_inter(i_spec) = num_inter(i_spec) + 1
              INTER_SPEC_ID_(i_ion_pair,i_interaction) = i_spec
              exit
            end if
          end do

          ! Set the minimum RH value
          key_name = "min RH"
          call assert_msg(519674084, &
                  interaction%get_real(key_name, &
                  MIN_RH_(i_ion_pair, i_interaction)), &
                  "Missing minimum RH value for ion pair '"// &
                  ion_pair_name//"' interaction with '"//inter_spec_name// &
                  "' in PD-FiTE activity reaction.")
          min_RH = MIN_RH_(i_ion_pair, i_interaction)
          call assert_msg(294172408, &
                  min_RH.ge.real(0.0, kind=dp).and. &
                  min_RH.lt.real(1.0, kind=dp), &
                  "Invalid value for minimum RH for ion pair '"// &
                  ion_pair_name//"' interaction with '"//inter_spec_name// &
                  "' in PD-FiTE activity reaction: "//to_string(min_RH))

          ! Set the maximum RH value
          key_name = "max RH"
          call assert_msg(649525712, &
                  interaction%get_real(key_name, &
                  MAX_RH_(i_ion_pair, i_interaction)), &
                  "Missing maximum RH value for ion pair '"// &
                  ion_pair_name//"' interaction with '"//inter_spec_name// &
                  "' in PD-FiTE activity reaction.")
          max_RH = MAX_RH_(i_ion_pair, i_interaction)
          call assert_msg(840423507, &
                  max_RH.gt.real(0.0, kind=dp).and. &
                  max_RH.le.real(1.0, kind=dp), &
                  "Invalid value for maximum RH for ion pair '"// &
                  ion_pair_name//"' interaction with '"//inter_spec_name// &
                  "' in PD-FiTE activity reaction: "//to_string(max_RH))

          ! Get the number of B_z parameters
          key_name = "B"
          call assert(981216440, &
                  interaction%get_property_t(key_name, poly_coeffs))
          NUM_B_(i_ion_pair, i_interaction) = poly_coeffs%size()

          ! Get the B_z parameters
          call poly_coeffs%iter_reset()
          do i_poly_coeff = 1, poly_coeffs%size()
            call assert_msg(988796931, poly_coeffs%get_real(val=real_val), &
                    "Invalid polynomial coefficient for ion_pair '"// &
                    ion_pair_name//"' interaction with '"//inter_spec_name// &
                    "' in PDFiTE activity reaction.")
            B_Z_(i_ion_pair, i_interaction, i_poly_coeff) = real_val
            call poly_coeffs%iter_next()
          end do

          n_float_param = n_float_param + 2 + NUM_B_(i_ion_pair, i_interaction)
          n_int_param = n_int_param + 3

          call interactions%iter_next()
        end do

        ! Check that all the interactions were included at least once
        do i_spec = 1, size(num_inter)
          call warn_assert_msg(793223082, num_inter(i_spec).ge.1, &
                  "Missing interaction parameters between ion_pair '"// &
                  ion_pair_name//"' and '"//ion_pair_names(i_spec)%string// &
                  "' in PDFiTE activity reaction: "// &
                  trim(to_string(num_inter(i_spec)))//".")
        end do

        ! Make sure no interactions with the same ion pair overlap in
        ! their RH ranges and that the entire RH range (0.0-1.0) is covered.
        ! Treat ranges as R = (minRH,maxRH] to avoid overlap at boundaries
        rh_range(:) = 0.0
        do i_interaction = 1, NUM_INTER_(i_ion_pair)

          ! Check for RH range overlaps with other interactions with this
          ! ion pair
          do j_interaction = i_interaction+1, NUM_INTER_(i_ion_pair)
            if (INTER_SPEC_ID_(i_ion_pair, i_interaction).ne. &
                    INTER_SPEC_ID_(i_ion_pair, j_interaction)) cycle
            call assert_msg(858420675, &
                    MIN_RH_(i_ion_pair, i_interaction).ge. &
                    MAX_RH_(i_ion_pair, j_interaction).or. &
                    MAX_RH_(i_ion_pair, i_interaction).le. &
                    MIN_RH_(i_ion_pair, j_interaction), &
                    "Overlapping RH range for interactions for ion pair '"// &
                    ion_pair_name//"' in PD-FiTE activity reaction.")
          end do

          ! Add the constribution from this interaction to RH coverage for
          ! the interacting ion pair
          rh_range(INTER_SPEC_ID_(i_ion_pair, i_interaction)) = &
                  rh_range(INTER_SPEC_ID_(i_ion_pair, i_interaction)) + &
                  (MAX_RH_(i_ion_pair, i_interaction) - &
                   MIN_RH_(i_ion_pair, i_interaction))
        end do

        ! Make sure the entire RH range is covered
        do i_spec = 1, size(rh_range)
          call assert_msg(370258071, &
                  rh_range(i_spec).eq.real(1.0, kind=dp).or. &
                  num_inter(i_spec).eq.0, &
                  "Incomplete RH coverage for interaction with ion pair '"// &
                  ion_pair_names(i_spec)%string//"' for '"//ion_pair_name// &
                  "' PD-FiTE activity coefficient calculation.")
        end do

        call ion_pairs%iter_next()

      ! The last portion of ion_pair_names includes ion pairs that are
      ! included in the interactions but for which acitivty coefficients
      ! are not calculated
      else

        ! Set the number of interactions to zero to indicate not to calculate
        ! activity coefficients for this ion pair
        NUM_INTER_(i_ion_pair) = 0

      end if

    end do

    call assert(938415336, n_int_param.eq.TOTAL_INT_PARAM_)
    call assert(433208931, n_float_param.eq.TOTAL_FLOAT_PARAM_)

    deallocate(ion_pair_names)
    deallocate(num_inter)
    deallocate(rh_range)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(sub_model_PDFiTE_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef NUM_PHASE_
#undef GAS_WATER_ID_
#undef NUM_ION_PAIRS_
#undef TOTAL_INT_PARAM_
#undef TOTAL_FLOAT_PARAM_
#undef PPM_TO_RH_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef PHASE_ID_
#undef PAIR_INT_PARAM_LOC_
#undef PAIR_FLOAT_PARAM_LOC_
#undef ION_PAIR_ACT_ID_
#undef NUM_CATION_
#undef NUM_ANION_
#undef CATION_ID_
#undef ANION_ID_
#undef NUM_INTER_
#undef NUM_B_
#undef INTER_SPEC_ID_
#undef INTER_SPEC_LOC_
#undef CATION_MW_
#undef ANION_MW_
#undef CATION_N_
#undef ANION_N_
#undef MIN_RH_
#undef MAX_RH_
#undef B_Z_

end module pmc_sub_model_PDFiTE
