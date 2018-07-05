! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_ZSR_aerosol_water module.

! TODO Incorporate deliquesence calculations
! FIXME Move to a sub-model

!> \page phlex_rxn_ZSR_aerosol_water Phlexible Module for Chemistry: ZSR Aerosol Water Reaction
!!
!! ZSR aerosol water reactions calculate equilibrium aerosol water content
!! based on the Zdanovski-Stokes-Robinson mixing rule \cite Jacobson1996 in
!! the following generalized format:
!!
!! \f[
!!   W = \sum\limits_{i=0}^{n}\frac{1000 M_i}{MW_i m_{i}(a_w)}
!! \f]
!!
!! where \f$M\f$ is the concentration of binary electrolyte \f$i\f$ 
!! (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$) with molecular weight
!! \f$MW_i\f$ (\f$\mbox{\si{\kilo\gram\per\mole}}\f$) and molality
!! \f$m_{i}\f$ at a given water activity \f$a_w\f$ (RH; 0--1)
!! contributing to the total aerosol water content \f$W\f$ 
!! (\f$\mbox{\si{\micro\gram\per\cubic\metre}}\f$).
!!
!! Input data for ZSR aerosol water calculations have the following format :
!! \code{.json}
!!   {
!!     "type" : "ZSR_AEROSOL_WATER",
!!     "aerosol phase" : "my aero phase",
!!     "gas-phase water" : "H2O",
!!     "aerosol-phase water" : "H2O_aq",
!!     "ion pairs" : {
!!       "Na2SO4" : {
!!         "type" : "JACOBSON",
!!         "ions" : {
!!           "Nap" : { "qty" : 2 },
!!           "SO4mm" : {}
!!         },
!!         "Y_j" : [-3.295311e3, 3.188349e4, -1.305168e5, 2.935608e5],
!!         "low RH" : 0.51
!!       },
!!       "H2SO4" : {
!!         "type" : "EQSAM",
!!         "ions" : {
!!           "SO4mm" : {}
!!         },
!!         "NW" : 4.5,
!!         "ZW" : 0.5,
!!         "MW" : 98.0
!!       }
!!       ...
!!     }
!!   }
!! \endcode
!! The key-value pair \b aerosol \b phase is required to specify the aerosol
!! phase for which to calculate water content. Key-value pairs 
!! \b gas-phase \b water and \b aerosol-phase \b water must also be present
!! and specify the names for the water species in each phase. The final
!! required key-value pair is \b ion \b pairs which should contain a set of
!! key-value pairs where the key of each member of the set is the name of a
!! binary  electrolyte and the contents contain parameters required to
!! estimate the contribution of the this electrolyte to total aerosol water.
!! The name of the electrolyte may or may not refer to an actual aerosol-phase
!! species.
!!
!! Each binary electrolyte must include a \b type that refers to a method
!! of calculating ion-pair contributions to aerosol water. Valid values for
!! \b type are \b JACOBSON and \b EQSAM. These are described next.
!!
!! Aerosol water from ion pairs with type \b JACOBSON use equations (28) and
!! (29) in Jacobson et al. (1996) \cite Jacobson1996 where experimentally
!! determined binary solution molalities are fit to a polynomial as:
!!
!! \f[
!!   \sqrt{m_{i}(a_w)} = Y_0 + Y_1 a_w + Y_2 a_w^2 + Y_3 a_w^3 + ...,
!! \f]
!!
!! where \f$Y_j\f$ are the fitting parameters. Thus, \f$m_i(a_w)\f$ is
!! calculated at each time step, assuming constant \f$a_w\f$. These values
!! must be included in a key-value pair \b Y_j whose value is an array
!! with the \f$Y_j\f$ parameters. The size of the array corresponds to the
!! order of the polynomial equation, which must be greater than 1. The
!! key-value pair \b low \b RH is required to specify the lowest RH (0--1) 
!! for which this fit is valid. This value for RH will be used for all lower
!! RH in calculations of \f$m_i(a_w)\f$ as per Jacobson et al. (1996) 
!! \cite Jacobson1996.
!!
!! The key-value pair \b ions must contain the set of ions this binary
!! electrolyte includes. Each species must correspond to a species present in
!! \b aerosol \b phase and  have a \b charge parameter that specifies their
!! charge (uncharged species are not permitted in this set) and a 
!! \b molecular \b weight (\f$\mbox{\si{\kilo\gram\per\mole}}\f$) property.
!! Ions without a \b qty specified are assumed to appear once in the binary
!! electrolyte. The total molecular weight for the binary electroly
!! \f$MW_i\f$ is calculated as a sum of its ionic components, and the ion
!! species concentrations are used to determine the \f$M_i\f$ during
!! integration.
!!
!! For the above example, the following input data should be present:
!! \code{.json}
!! {
!!   "name" : "H2O",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "GAS",
!! },  
!! {
!!   "name" : "Nap",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : 1,
!!   "molecular weight" : 0.0229898
!! },
!! {
!!   "name" : "SO4mm",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : -2,
!!   "molecular weight" : 0.09606
!! },
!! {
!!   "name" : "my aero phase",
!!   "type" : "AERO_PHASE",
!!   "species" : ["Nap", "SO4mm", H2O_aq"]
!! }
!! \endcode
!! Aerosol water from ion pairs with type \b EQSAM use the parameterization of
!! Metzger et al. (2002) \cite Metzger2002 for aerosol water content:
!!
!! \f[
!!   \sqrt{m_{i}(a_w)} = (NW_i MW_{H2O}/MW_i 1/(a_w-1))^{ZW_i}
!! \f]
!!
!! where \f$NW_i\f$ and \f$ZW_i\f$ are fitting parameters \cite Metzger2002,
!! and must be provided in key-value pairs \b NW and \b ZW, along with the
!! binary electrolyte molecular weight \b MW
!! (\f$\mbox{\si{\kilo\gram\per\mole}}\f$). The key-value pair \b ions must
!! contain a set of ions that can be summed to calculate \f$M_i\f$ at runtime.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! TODO Find a way to incorporate the "regimes" in EQSAM

!> The rxn_ZSR_aerosol_water_t type and associated functions. 
module pmc_rxn_ZSR_aerosol_water

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_constants,                        only: const
  use pmc_phlex_state
  use pmc_property
  use pmc_rxn_data
  use pmc_util,                             only: i_kind, dp, to_string, &
                                                  assert, assert_msg, &
                                                  die_msg, string_t, &
                                                  align_ratio

  implicit none
  private
#define ACT_CALCJACOBSON 1
#define ACT_CALCEQSAM 2

#define NUM_PHASE_ this%condensed_data_int(1)
#define GAS_WATER_ID_ this%condensed_data_int(2)
#define NUM_ION_PAIR_ this%condensed_data_int(3)
#define TOTAL_INT_PARAM_ this%condensed_data_int(4)
#define TOTAL_FLOAT_PARAM_ this%condensed_data_int(5)
#define INT_DATA_SIZE_ this%condensed_data_int(6)
#define FLOAT_DATA_SIZE_ this%condensed_data_int(7)
#define PPM_TO_RH_ this%condensed_data_real(1)
#define NUM_INT_PROP_ 7
#define NUM_REAL_PROP_ 1
#define PHASE_ID_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define PAIR_INT_PARAM_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_PHASE_+x)
#define PAIR_FLOAT_PARAM_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_PHASE_+NUM_ION_PAIR_+x)
#define TYPE_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x))
#define JACOB_NUM_CATION_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+1)
#define JACOB_NUM_ANION_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+2)
#define JACOB_CATION_ID_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+3)
#define JACOB_ANION_ID_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+4)
#define JACOB_NUM_Y_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+5)
#define EQSAM_NUM_ION_(x) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+1)
#define EQSAM_ION_ID_(x,y) this%condensed_data_int(PAIR_INT_PARAM_LOC_(x)+1+y)
#define JACOB_low_RH_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x))
#define JACOB_CATION_MW_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+1)
#define JACOB_ANION_MW_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+2)
#define JACOB_Y_(x,y) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+2+y)
#define EQSAM_NW_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x))
#define EQSAM_ZW_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+1)
#define EQSAM_ION_PAIR_MW_(x) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+2)
#define EQSAM_ION_MW_(x,y) this%condensed_data_real(PAIR_FLOAT_PARAM_LOC_(x)+2+y)

  public :: rxn_ZSR_aerosol_water_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_ZSR_aerosol_water_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Finalize
    final :: finalize
  end type rxn_ZSR_aerosol_water_t

  !> Constructor for rxn_ZSR_aerosol_water_t
  interface rxn_ZSR_aerosol_water_t
    procedure :: constructor
  end interface rxn_ZSR_aerosol_water_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for ZSR aerosol water reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_ZSR_aerosol_water_t), pointer :: new_obj

    allocate(new_obj)
    new_obj%rxn_phase = AERO_RXN

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the reaction data, validating component data and loading
  !! any required information into the condensed data arrays for use during
  !! solving
  subroutine initialize(this, chem_spec_data, aero_rep)
    
    !> Reaction data
    class(rxn_ZSR_aerosol_water_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props, ion_pairs, ion_pair, sub_props, &
            ions
    character(len=:), allocatable :: key_name, spec_name, phase_name
    integer(kind=i_kind) :: n_phase, n_ion_pair, n_int_param, n_float_param, &
            i_aero_rep, i_phase, i_ion_pair, i_ion, i_spec, i_sub_prop, &
            qty, int_val, charge, total_charge
    real(kind=dp) :: real_val, molecular_weight
    type(string_t), allocatable :: unique_spec_names(:)
    character(len=:), allocatable :: str_type, ion_pair_name, ion_name
    integer(kind=i_kind) :: int_data_size, float_data_size

    ! Get the reaction property set
    if (.not. associated(this%property_set)) call die_msg(344693903, &
            "Missing property set needed to initialize ZSR aerosol water"// &
            " reaction.")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(613734060, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in ZSR aerosol water reaction.")

    ! Count the instances of the aerosol phase
    n_phase = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Get the number of instances of the phase in this representation
      n_phase = n_phase + &
              aero_rep(i_aero_rep)%val%num_phase_instances(phase_name)

    end do

    call assert_msg(450206321, n_phase.gt.0, &
            "Aerosol phase '"//phase_name//"' not present in any aerosol "// &
            "representation for ZSR aerosol water reaction.")

    ! Get the ion pairs
    key_name = "ion pairs"
    call assert_msg(640916964, &
            this%property_set%get_property_t(key_name, ion_pairs), &
            "Missing ion pairs in ZSR aerosol water reaction.")

    ! Count the ion pairs
    n_ion_pair = ion_pairs%size()
    call assert_msg(743158990, n_ion_pair .gt. 0, &
            "Empty ion pair set in ZSR aerosol water reaction.")

    ! Get the number of parameters required for each ion pair
    n_int_param = NUM_INT_PROP_ + n_phase + 2*n_ion_pair
    n_float_param = NUM_REAL_PROP_
    call ion_pairs%iter_reset()
    do i_ion_pair = 1, n_ion_pair
  
      ! Get the name of the ion pair
      call assert(476976534, ion_pairs%get_key(ion_pair_name))

      ! Get the ion pair properties
      call assert_msg(280814432, ion_pairs%get_property_t(val=ion_pair), &
              "Missing ion pair properties for '"//ion_pair_name// &
              "' in ZSR aerosol water reaction.")

      ! Get the activity calculation type
      key_name = "type"
      call assert_msg(334930304, ion_pair%get_string(key_name, str_type), &
              "Missing activity calculation type for ion pair '"// &
              ion_pair_name//"' in ZSR aerosol water reaction.")

      ! Get the number of parameters according to activity calculation type
      if (str_type.eq."JACOBSON") then
        
        ! Get the number of Y_j parameters
        key_name = "Y_j"
        call assert_msg(286454243, &
                ion_pair%get_property_t(key_name, sub_props), &
                "Missing Y_j parameters for Jacobson activity calculation "//&
                "for ion pair '"//ion_pair_name//"'.")

        call assert_msg(495036486, sub_props%size().gt.0, &
                "Insufficient Y_j parameters for Jacobson activity "// &
                "calculation for ion pair '"//ion_pair_name//"' in "// &
                "ZSR aerosol water reaction.")

        n_float_param = n_float_param + 3 + sub_props%size()
        n_int_param = n_int_param + 6

      else if (str_type.eq."EQSAM") then

        ! Get the number of ions
        key_name = "ions"
        call assert_msg(244982851, &
                ion_pair%get_property_t(key_name, sub_props), &
                "Mission ions for EQSAM activity calculation for ion "// &
                "pair '"//ion_pair_name//"' in ZSR aerosol water "// &
                "reaction.")
        
        call assert_msg(849524804, sub_props%size().gt.0, &
                "Insufficient ions specified for EQSAM activity "// &
                "calculation for ion pair '"//ion_pair_name// &
                "' in ZSR aerosol water reaction.")

        n_float_param = n_float_param + 3 + sub_props%size()
        n_int_param = n_int_param + 2 + sub_props%size()

      else
        call die_msg(704759248, "Invalid activity type specified for ZSR "// &
                "aerosol water reaction: '"//str_type//"'")
      end if

      call ion_pairs%iter_next()
    end do

    ! Calculate int and float array sizes with alignment spacing
    int_data_size = n_int_param
    int_data_size = int_data_size + mod(int_data_size, align_ratio)
    float_data_size = n_float_param

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(int_data_size))
    allocate(this%condensed_data_real(float_data_size))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)
    INT_DATA_SIZE_ = int_data_size
    FLOAT_DATA_SIZE_ = float_data_size

    ! Set some data dimensions
    NUM_PHASE_  = n_phase
    NUM_ION_PAIR_ = n_ion_pair
    TOTAL_INT_PARAM_ = n_int_param
    TOTAL_FLOAT_PARAM_ = n_float_param

    ! Set the gas-phase water species
    key_name = "gas-phase water"
    call assert_msg(386389634, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing gas-phase water species name in ZSR aerosol water "// &
            "reaction.")

    GAS_WATER_ID_ = chem_spec_data%gas_state_id(spec_name)
    
    call assert_msg(709909577, GAS_WATER_ID_ .gt. 0, &
            "Cannot find gas-phase water species '"//spec_name//"' for "// &
            "ZSR aerosol water reaction.")
    
    ! Set the aerosol-water species
    key_name = "aerosol-phase water"
    call assert_msg(771445226, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing aerosol-phase water species name in ZSR aerosol "// &
            "water reaction.")

    ! Make the PHASE_ID_(x) hold the state id of aerosol water in each
    ! phase instance. Then the aerosol water id is 1, and the ion
    ! ids will be relative to the water id in each phase.
    i_phase = 1
    do i_aero_rep = 1, size(aero_rep)
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)
      if (.not.allocated(unique_spec_names)) cycle
      do i_spec = 1, size(unique_spec_names)
        PHASE_ID_(i_phase) = aero_rep(i_aero_rep)%val%spec_state_id( &
                unique_spec_names(i_spec)%string)
        call assert(204327668, PHASE_ID_(i_phase).gt.0)
        i_phase = i_phase + 1
      end do
      deallocate(unique_spec_names)
    end do
    i_phase = i_phase - 1
    call assert_msg(418435744, i_phase.eq.NUM_PHASE_, &
            "Incorrect number of aerosol water instances in ZSR aerosol "// &
            "water reaction. Expected "//trim(to_string(NUM_PHASE_))// &
            " but got "//trim(to_string(i_phase)))

    ! Save the ion-pair parameters
    n_int_param = NUM_INT_PROP_ + NUM_PHASE_ + 2*NUM_ION_PAIR_
    n_float_param = NUM_REAL_PROP_
    call ion_pairs%iter_reset()
    do i_ion_pair = 1, n_ion_pair
   
      ! Get the name of the ion pair
      call assert(476976534, ion_pairs%get_key(ion_pair_name))

      ! Get the ion pair properties
      call assert(660267400, ion_pairs%get_property_t(val=ion_pair))

      ! Set the location of this ion pair's parameters in the condensed data
      ! arrays.
      PAIR_INT_PARAM_LOC_(i_ion_pair) = n_int_param + 1
      PAIR_FLOAT_PARAM_LOC_(i_ion_pair) = n_float_param + 1

      ! Get the activity calculation type
      key_name = "type"
      call assert(288245799, ion_pair%get_string(key_name, str_type))

      ! Get the number of parameters according to activity calculation type
      if (str_type.eq."JACOBSON") then
       
        ! Set the type
        TYPE_(i_ion_pair) = ACT_CALCJACOBSON

        ! Get the Y_j parameters
        key_name = "Y_j"
        call assert(227500762, ion_pair%get_property_t(key_name, sub_props))
        JACOB_NUM_Y_(i_ion_pair) = sub_props%size()
        call sub_props%iter_reset()
        do i_sub_prop = 1, sub_props%size()
          call assert_msg(149509565, sub_props%get_real(val=real_val), &
                  "Invalid Y parameter for ion pair '"// &
                  ion_pair_name//"' in ZSR aerosol water reaction.")
          JACOB_Y_(i_ion_pair, i_sub_prop) = real_val
          call sub_props%iter_next()
        end do

        ! Get the low RH value
        key_name = "low RH"
        call assert_msg(462500894, ion_pair%get_real(key_name, real_val), &
                "Missing 'low RH' value for ion pair '"// &
                ion_pair_name//"' in ZSR aerosol water reaction.")
        JACOB_low_RH_(i_ion_pair) = real_val

        ! Get the number and id of the ions and the ion-pair molecular weight
        molecular_weight = 0.0
        key_name = "ions"
        call assert_msg(661006818, ion_pair%get_property_t(key_name, ions), &
                "Mission ions for Jacobson activity calculation for ion "// &
                "pair '"//ion_pair_name//"' in ZSR aerosol water "// &
                "reaction.")
        call assert_msg(880831496, ions%size().eq.2, &
                "Invalid number of unique ions specified for ion pair '"// &
                ion_pair_name//"' in for Jacobson activity "// &
                "calculation in ZSR aerosol water reaction. Expected 2 "// &
                "got "//trim(to_string(ions%size())))
        call ions%iter_reset()
        total_charge = 0
        do i_ion = 1, 2
        
          ! Get the ion name
          call assert(849711956, ions%get_key(ion_name))

          ! Get the qty, if specified
          qty = 1
          if (ions%get_property_t(val=sub_props)) then
            key_name = "qty"
            if (sub_props%get_int(key_name, int_val)) qty = int_val 
          end if 

          ! Get the species properties
          call assert_msg(315479897, &
                  chem_spec_data%get_property_set(ion_name, spec_props), &
                  "Missing species properties for ion '"//ion_name// &
                  "' in ZSR aerosol water reaction.")
          
          ! Add the molecular weight
          key_name = "molecular weight"
          call assert_msg(897812513, &
                  spec_props%get_real(key_name, molecular_weight), &
                  "Missing molecular weight for ion '"//ion_name// &
                  "' in ZSR aerosol water reaction.")

          ! Add the charge from this species
          key_name = "charge"
          call assert_msg(310667885, spec_props%get_int(key_name, charge), &
                  "Missing charge for ion '"//ion_name//"' in ZSR "// &
                  "aerosol water reaction.")

          if (charge.gt.0) then
            JACOB_NUM_CATION_(i_ion_pair) = qty
            JACOB_CATION_MW_(i_ion_pair) = molecular_weight
          else if (charge.lt.0) then
            JACOB_NUM_ANION_(i_ion_pair) = qty
            JACOB_ANION_MW_(i_ion_pair) = molecular_weight
          else
            call die_msg(899416917, "Neutral species '"//ion_name// &
                    "' not allowed in ZSR aerosol water reaction ion pair")
          end if

          ! Add contribution to total charge
          total_charge = total_charge + qty * charge

          ! Get the state ids for this species
          i_phase = 1
          do i_aero_rep = 1, size(aero_rep)
            unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
                    phase_name = phase_name, spec_name = ion_name)
            if (.not.allocated(unique_spec_names)) cycle
            do i_spec = 1, size(unique_spec_names)
              if (charge.gt.0) then
                if (i_phase.eq.1) then
                  JACOB_CATION_ID_(i_ion_pair) = &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          PHASE_ID_(i_phase)
                else
                  call assert(473680545, JACOB_CATION_ID_(i_ion_pair).eq. &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          PHASE_ID_(i_phase))
                end if
              else
                if (i_phase.eq.1) then
                  JACOB_ANION_ID_(i_ion_pair) = &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          PHASE_ID_(i_phase)
                else
                  call assert(234155524, JACOB_ANION_ID_(i_ion_pair).eq. &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          PHASE_ID_(i_phase))
                end if
              end if
              i_phase = i_phase + 1
            end do
            deallocate(unique_spec_names)
          end do
          i_phase = i_phase - 1
          call assert_msg(623684811, i_phase.eq.NUM_PHASE_, &
            "Incorrect number of instances of ion species '"// &
            ion_name//"' in ZSR aerosol water reaction. Expected "// &
            trim(to_string(NUM_PHASE_))//" but got "// &
            trim(to_string(i_phase)))

          ! Get the next ion
          call ions%iter_next()

        end do
    
        call assert_msg(319151390, total_charge.eq.0, &
                "Charge imbalance for ion pair '"//ion_pair_name// &
                " in ZSR aerosol water reaction. Total charge: "// &
                trim(to_string(total_charge)))

        n_float_param = n_float_param + 3 + JACOB_NUM_Y_(i_ion_pair)
        n_int_param = n_int_param + 6

      else if (str_type.eq."EQSAM") then

        ! Set the type
        TYPE_(i_ion_pair) = ACT_CALCEQSAM

        ! Get the required parameters for calculating activity
        key_name = "NW"
        call assert_msg(692339107, &
                ion_pair%get_real(key_name, EQSAM_NW_(i_ion_pair)), &
                "Missing parameter NW for ion pair '"//ion_pair_name// &
                "' in ZSR aerosol water reaction.")

        key_name = "ZW"
        call assert_msg(894917625, &
                ion_pair%get_real(key_name, EQSAM_ZW_(i_ion_pair)), &
                "Missing parameter ZW for ion pair '"//ion_pair_name// &
                "' in ZSR aerosol water reaction.")

        key_name = "MW"
        call assert_msg(272128568, &
                ion_pair%get_real(key_name, EQSAM_ION_PAIR_MW_(i_ion_pair)), &
                "Missing parameter MW for ion pair '"//ion_pair_name// &
                "' in ZSR aerosol water reaction.")

        ! Get the number and id of ions
        key_name = "ions"
        call assert(381088140, ion_pair%get_property_t(key_name, ions))
        EQSAM_NUM_ION_(i_ion_pair) = ions%size()
        call ions%iter_reset()
        do i_ion = 1, ions%size()
          
          ! Get the ion name
          call assert(849711956, ions%get_key(ion_name))

          ! Get the species properties
          call assert_msg(826137761, &
                  chem_spec_data%get_property_set(ion_name, spec_props), &
                  "Missing species properties for ion '"//ion_name// &
                  "' in ZSR aerosol water reaction.")
          
          ! Add the molecular weight
          key_name = "molecular weight"
          call assert_msg(598142298, &
                  spec_props%get_real(key_name, molecular_weight), &
                  "Missing molecular weight for ion '"//ion_name// &
                  "' in ZSR aerosol water reaction.")
          EQSAM_ION_MW_(i_ion_pair, i_ion) = molecular_weight

          ! Set the ion id (relative to water within the specified phase)
          i_phase = 1
          do i_aero_rep = 1, size(aero_rep)
            unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
                    phase_name = phase_name, spec_name = ion_name)
            if (.not.allocated(unique_spec_names)) cycle
            do i_spec = 1, size(unique_spec_names)
              if (i_phase.eq.1) then
                EQSAM_ION_ID_(i_ion_pair,i_ion) = &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        PHASE_ID_(i_phase)
              else
                call assert(973648240, EQSAM_ION_ID_(i_ion_pair,i_ion) .eq. &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        PHASE_ID_(i_phase))
              end if
              i_phase = i_phase + 1
            end do
            deallocate(unique_spec_names)
          end do
          i_phase = i_phase - 1
          call assert_msg(900921350, i_phase.eq.NUM_PHASE_, &
            "Incorrect number of instances of ion species '"// &
            ion_name//"' in ZSR aerosol water reaction. Expected "// &
            trim(to_string(NUM_PHASE_))//" but got "// &
            trim(to_string(i_phase)))

           ! Get the next ion
           call ions%iter_next()

        end do

        n_float_param = n_float_param + 3 + EQSAM_NUM_ION_(i_ion_pair)
        n_int_param = n_int_param + 2 + EQSAM_NUM_ION_(i_ion_pair)

      else
        call die_msg(186680407, "Internal error.")
      end if

      call ion_pairs%iter_next()
    end do

    call assert(859412771, n_int_param.eq.TOTAL_INT_PARAM_)
    call assert(568314442, n_float_param.eq.TOTAL_FLOAT_PARAM_)
      
  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the reaction
  elemental subroutine finalize(this)

    !> Reaction data
    type(rxn_ZSR_aerosol_water_t), intent(inout) :: this

    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef ACT_TYPE_JACOBSON
#undef ACT_TYPE_EQSAM

#undef NUM_PHASE_
#undef GAS_WATER_ID_
#undef NUM_ION_PAIR_
#undef TOTAL_INT_PARAM_
#undef TOTAL_FLOAT_PARAM_
#undef PPM_TO_RH_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef PHASE_ID_
#undef PAIR_INT_PARAM_LOC_
#undef PAIR_FLOAT_PARAM_LOC_
#undef TYPE_
#undef JACOB_NUM_CATION_
#undef JACOB_NUM_ANION_
#undef JACOB_CATION_ID_
#undef JACOB_ANION_ID_
#undef JACOB_NUM_Y_
#undef EQSAM_NUM_ION_
#undef EQSAM_ION_ID_
#undef JACOB_low_RH_
#undef JACOB_CATION_MW_
#undef JACOB_ANION_MW_
#undef JACOB_Y_
#undef EQSAM_NW_
#undef EQSAM_ZW_
#undef EQSAM_ION_PAIR_MW_
#undef EQSAM_ION_MW_

end module pmc_rxn_ZSR_aerosol_water
