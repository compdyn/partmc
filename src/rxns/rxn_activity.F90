! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_activity module.

! TODO Incorporate deliquesence calculations

!> \page phlex_rxn_activity Phlexible Mechanism for Chemistry: Activity Reaction
!!
!! activity reactions calculate equilibrium aerosol water content
!! based on the Zdanovski-Stokes-Robinson mixing rule \cite{Stokes1966,
!! Jacobson1996} in the following generalized format:
!!
!! \f[
!!   W = \sum\limits_{i=0}^{n}\frac{1000 M_i}{MW_i m_{i}(a_w)}
!! \f]
!!
!! where \f$M\f$ is the concentration of binary electrolyte \f$i\f$ 
!! (\f$\mu g m^{-3}\f$) with molecular weight \f$MW_i\f$ (g/mol) and 
!! molality \f$m_{i}\f$ at a given water activity \f$a_w\f$ (RH; 0-1)
!! contributing to the total aerosol water content \f$W\f$ 
!! (\f$\mu g m^{-3}\f$).
!!
!! Input data for activity equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "ACTIVITY",
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
!! The key-value pair \b "aerosol phase" is required to specify the aerosol
!! phase for which to calculate water content. Key-value pairs 
!! \b "gas-phase water" and \b "aerosol-phase water" must also be present and
!! specify the names for the water species in each phase. The final required
!! key-value pair is \b "ion pairs" which should contain a set of key-value
!! pairs where the key of each member of the set is the name of a binary
!! electrolyte and the contents contain parameters required to estimate
!! the contribution of the this electrolyte to total aerosol water. The
!! name of the electrolyte may or may not refer to an actual aerosol-phase
!! species.
!!
!! Each binary electrolyte must include a \b "type" that refers to a method
!! of calculating ion-pair contributions to aerosol water. Valid values for
!! \b "type" are "JACOBSON" and "EQSAM". These are described next.
!!
!! Aerosol water from ion pairs with type "JACOBSON" use equations (28) and
!! (29) in Jacobson et al. \cite{Jacobson1996} where experimentally determined
!! binary solution molalities are fit to a polynomial as:
!!
!! \f[
!!   \sqrt{m_{i}(a_w)} = Y_0 + Y_1 a_w + Y_2 a_w^2 + Y_3 a_w^3 + ...,
!! \f]
!!
!! where \f$Y_j\f$ are the fitting parameters. Thus, \f$m_i(a_w)\f$ is
!! calculated at each time step, assuming constant \f$a_w\f$. These values
!! must be included in a key-value pair \b "Y_j" whose value is an array
!! with the \f$Y_j\f$ parameters. The size of the array corresponds to the
!! order of the polynomial equation, which must be greater than 1. The
!! key-value pair \b "low RH" is required to specify the lowest RH (0-1) 
!! for which this fit is valid. This value for RH will be used for all lower
!! RH in calculations of \f$m_i(a_w)\f$ as per Jacobson et al. \cite{1996}.
!!
!! The key-value pair "ions" must contain the set of ions this binary
!! electrolyte includes. Each species must correspond to a species present in
!! \b "aerosol phase" and  have a \b "charge" parameter that specifies their
!! charge (uncharged species are not permitted in this set) and a 
!! \b "molecular weight" (g/mol) property. Ions without a \b "qty" specified
!! are assumed to appear once in the binary electrolyte. The total
!! molecular weight for the binary electrolye \f$MW_i\f$ is calculated as a 
!! sum of its ionic components, and the ion species concentrations are used
!! to determine the \f$M_i\f$ during integration.
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
!!   "molecular weight" : 22.9898
!! },
!! {
!!   "name" : "SO4mm",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : -2
!!   "molecular weight" : 96.06
!! },
!! {
!!   "name" : "my aero phase",
!!   "type" : "AERO_PHASE",
!!   "species" : ["Nap", "SO4mm", H2O_aq"]
!! }
!!
!! Aerosol water from ion pairs with type "EQSAM" use the parameterization of
!! Metzger et al. \cite{Metzget2002} for aerosol water content:
!!
!! \f[
!!   \sqrt{m_{i}(a_w)} = (NW_i MW_{H2O}/MW_i 1/(a_w-1))^{ZW_i}
!! \f]
!!
!! where \f$NW_i\f$ and \f$ZW_i\f$ are fitting parameters \cite{Metger2002},
!! and must be provided in key-value pairs \b "NW" and \b "ZW", along with the
!! binary electrolyte molecular weight \b "MW" (g/mol). The key-value pair
!! \b "ions" must contain a set of ions that can be summed to calculate
!! \f$M_i\f$ at runtime.
!!
!! TODO Find a way to incorporate the "regimes" in EQSAM
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_activity_t type and associated functions. 
module pmc_rxn_activity

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
#define ACT_CALC_JACOBSON 1
#define ACT_CALC_EQSAM 2

#define _NUM_PHASE_ this%condensed_data_int(1)
#define _GAS_WATER_ID_ this%condensed_data_int(2)
#define _NUM_ION_PAIR_ this%condensed_data_int(3)
#define _TOTAL_INT_PARAM_ this%condensed_data_int(4)
#define _TOTAL_FLOAT_PARAM_ this%condensed_data_int(5)
#define _ppm_TO_RH_ this%condensed_data_real(1)
#define _NUM_INT_PROP_ 5
#define _NUM_REAL_PROP_ 1
#define _PHASE_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+x)
#define _PAIR_INT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_PHASE_+x)
#define _PAIR_FLOAT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_PHASE_+_NUM_ION_PAIR_+x)
#define _TYPE_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x))
#define _JACOB_NUM_CATION_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+1)
#define _JACOB_NUM_ANION_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+2)
#define _JACOB_CATION_ID_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+3)
#define _JACOB_ANION_ID_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+4)
#define _JACOB_NUM_Y_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+5)
#define _EQSAM_NUM_ION_(x) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+1)
#define _EQSAM_ION_ID_(x,y) this%condensed_data_int(_PAIR_INT_PARAM_LOC_(x)+1+y)
#define _JACOB_low_RH_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x))
#define _JACOB_CATION_MW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+1)
#define _JACOB_ANION_MW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+2)
#define _JACOB_Y_(x,y) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+2+y)
#define _EQSAM_NW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x))
#define _EQSAM_ZW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+1)
#define _EQSAM_ION_PAIR_MW_(x) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+2)
#define _EQSAM_ION_MW_(x,y) this%condensed_data_real(_PAIR_FLOAT_PARAM_LOC_(x)+2+y)

  public :: rxn_activity_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_activity_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
  end type rxn_activity_t

  !> Constructor for rxn_activity_t
  interface rxn_activity_t
    procedure :: constructor
  end interface rxn_activity_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for activity reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_activity_t), pointer :: new_obj

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
    class(rxn_activity_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props, ion_pairs, ion_pair, sub_props, &
            ions
    character(len=:), allocatable :: key_name, water_name, spec_name, &
            phase_name, string_val
    integer(kind=i_kind) :: n_phase, n_ion_pair, n_int_param, n_float_param
    integer(kind=i_kind) :: i_aero_rep, i_phase, i_ion_pair, i_ion, i_spec, &
            i_sub_prop
    integer(kind=i_kind) :: qty, int_val, charge, total_charge
    real(kind=dp) :: real_val, molecular_weight
    class(string_t), allocatable :: unique_spec_names(:)
    character(len=:), allocatable :: str_type, ion_pair_name, ion_name

    ! Get the reaction property set
    if (.not. associated(this%property_set)) call die_msg(344693903, &
            "Missing property set needed to initialize activity"// &
            " reaction.")

    ! Get the aerosol phase name
    key_name = "aerosol phase"
    call assert_msg(613734060, &
            this%property_set%get_string(key_name, phase_name), &
            "Missing aerosol phase in activity reaction.")

    ! Count the instances of the aerosol phase
    n_phase = 0
    do i_aero_rep = 1, size(aero_rep)

      ! Get the number of instances of the phase in this representation
      n_phase = n_phase + &
              aero_rep(i_aero_rep)%val%num_phase_instances(phase_name)

    end do

    call assert_msg(450206321, n_phase.gt.0, &
            "Aerosol phase '"//phase_name//"' not present in any aerosol "// &
            "representation for activity reaction.")

    ! Get the ion pairs
    key_name = "ion pairs"
    call assert_msg(640916964, &
            this%property_set%get_property_t(key_name, ion_pairs), &
            "Missing ion pairs in activity reaction.")

    ! Count the ion pairs
    n_ion_pair = ion_pairs%size()
    call assert_msg(743158990, n_ion_pair .gt. 0, &
            "Empty ion pair set in activity reaction.")

    ! Get the number of parameters required for each ion pair
    n_int_param = _NUM_INT_PROP_ + n_phase + 2*n_ion_pair
    n_float_param = _NUM_REAL_PROP_
    call ion_pairs%iter_reset()
    do i_ion_pair = 1, n_ion_pair
  
      ! Get the name of the ion pair
      call assert(476976534, ion_pairs%get_key(ion_pair_name))

      ! Get the ion pair properties
      call assert_msg(280814432, ion_pairs%get_property_t(val=ion_pair), &
              "Missing ion pair properties for '"//ion_pair_name// &
              "' in activity reaction.")

      ! Get the activity calculation type
      key_name = "type"
      call assert_msg(334930304, ion_pair%get_string(key_name, str_type), &
              "Missing activity calculation type for ion pair '"// &
              ion_pair_name//"' in activity reaction.")

      ! Get the number of parameters according to activity calculation type
      if (str_type.eq."JACOBSON") then
        
        ! Get the number of Y_j parameters
        key_name = "Y_j"
        call assert_msg(286454243, &
                ion_pair%get_property_t(key_name, sub_props), &
                "Missing Y_j parameters for Jacobson activity calculation "// &
                "for ion pair '"//ion_pair_name//"'.")

        call assert_msg(495036486, sub_props%size().gt.0, &
                "Insufficient Y_j parameters for Jacobson activity "// &
                "calculation for ion pair '"//ion_pair_name//"' in "// &
                "activity reaction.")

        n_float_param = n_float_param + 3 + sub_props%size()
        n_int_param = n_int_param + 6

      else if (str_type.eq."EQSAM") then

        ! Get the number of ions
        key_name = "ions"
        call assert_msg(244982851, &
                ion_pair%get_property_t(key_name, sub_props), &
                "Mission ions for EQSAM activity calculation for ion "// &
                "pair '"//ion_pair_name//"' in activity "// &
                "reaction.")
        
        call assert_msg(849524804, sub_props%size().gt.0, &
                "Insufficient ions specified for EQSAM activity "// &
                "calculation for ion pair '"//ion_pair_name// &
                "' in activity reaction.")

        n_float_param = n_float_param + 3 + sub_props%size()
        n_int_param = n_int_param + 2 + sub_props%size()

      else
        call die_msg(704759248, "Invalid activity type specified for "// &
                "activity reaction: '"//str_type//"'")
      end if

      call ion_pairs%iter_next()
    end do


    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(n_int_param))
    allocate(this%condensed_data_real(n_float_param))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set some data dimensions
    _NUM_PHASE_  = n_phase
    _NUM_ION_PAIR_ = n_ion_pair
    _TOTAL_INT_PARAM_ = n_int_param
    _TOTAL_FLOAT_PARAM_ = n_float_param

    ! Set the gas-phase water species
    key_name = "gas-phase water"
    call assert_msg(386389634, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing gas-phase water species name in activity "// &
            "reaction.")

    _GAS_WATER_ID_ = chem_spec_data%gas_state_id(spec_name)
    
    call assert_msg(709909577, _GAS_WATER_ID_ .gt. 0, &
            "Cannot find gas-phase water species '"//spec_name//"' for "// &
            "activity reaction.")
    
    ! Set the aerosol-water species
    key_name = "aerosol-phase water"
    call assert_msg(771445226, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing aerosol-phase water species name in activity "// &
            "reaction.")

    ! Make the _PHASE_ID_(x) hold the state id of aerosol water in each
    ! phase instance. Then the aerosol water id is 1, and the ion
    ! ids will be relative to the water id in each phase.
    i_phase = 1
    do i_aero_rep = 1, size(aero_rep)
      unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
              phase_name = phase_name, spec_name = spec_name)
      if (.not.allocated(unique_spec_names)) cycle
      do i_spec = 1, size(unique_spec_names)
        _PHASE_ID_(i_phase) = aero_rep(i_aero_rep)%val%spec_state_id( &
                unique_spec_names(i_spec)%string)
        call assert(204327668, _PHASE_ID_(i_phase).gt.0)
        i_phase = i_phase + 1
      end do
    end do
    i_phase = i_phase - 1
    call assert_msg(418435744, i_phase.eq._NUM_PHASE_, &
            "Incorrect number of aerosol water instances in activity "// &
            "reaction. Expected "//trim(to_string(_NUM_PHASE_))// &
            " but got "//trim(to_string(i_phase)))

    ! Save the ion-pair parameters
    n_int_param = _NUM_INT_PROP_ + _NUM_PHASE_ + 2*_NUM_ION_PAIR_
    n_float_param = _NUM_REAL_PROP_
    call ion_pairs%iter_reset()
    do i_ion_pair = 1, n_ion_pair
   
      ! Get the name of the ion pair
      call assert(476976534, ion_pairs%get_key(ion_pair_name))

      ! Get the ion pair properties
      call assert(660267400, ion_pairs%get_property_t(val=ion_pair))

      ! Set the location of this ion pair's parameters in the condensed data
      ! arrays.
      _PAIR_INT_PARAM_LOC_(i_ion_pair) = n_int_param + 1
      _PAIR_FLOAT_PARAM_LOC_(i_ion_pair) = n_float_param + 1

      ! Get the activity calculation type
      key_name = "type"
      call assert(288245799, ion_pair%get_string(key_name, str_type))

      ! Get the number of parameters according to activity calculation type
      if (str_type.eq."JACOBSON") then
       
        ! Set the type
        _TYPE_(i_ion_pair) = ACT_CALC_JACOBSON

        ! Get the Y_j parameters
        key_name = "Y_j"
        call assert(227500762, ion_pair%get_property_t(key_name, sub_props))
        _JACOB_NUM_Y_(i_ion_pair) = sub_props%size()
        call sub_props%iter_reset()
        do i_sub_prop = 1, sub_props%size()
          call assert_msg(149509565, sub_props%get_real(val=real_val), &
                  "Invalid Y parameter for ion pair '"// &
                  ion_pair_name//"' in activity reaction.")
          _JACOB_Y_(i_ion_pair, i_sub_prop) = real_val
          call sub_props%iter_next()
        end do

        ! Get the low RH value
        key_name = "low RH"
        call assert_msg(462500894, ion_pair%get_real(key_name, real_val), &
                "Missing 'low RH' value for ion pair '"// &
                ion_pair_name//"' in activity reaction.")
        _JACOB_low_RH_(i_ion_pair) = real_val

        ! Get the number and id of the ions and the ion-pair molecular weight
        molecular_weight = 0.0
        key_name = "ions"
        call assert_msg(661006818, ion_pair%get_property_t(key_name, ions), &
                "Mission ions for Jacobson activity calculation for ion "// &
                "pair '"//ion_pair_name//"' in activity "// &
                "reaction.")
        call assert_msg(880831496, ions%size().eq.2, &
                "Invalid number of unique ions specified for ion pair '"// &
                ion_pair_name//"' in for Jacobson activity "// &
                "calculation in activity reaction. Expected 2 "// &
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
                  "' in activity reaction.")
          
          ! Add the molecular weight
          key_name = "molecular weight"
          call assert_msg(897812513, &
                  spec_props%get_real(key_name, molecular_weight), &
                  "Missing molecular weight for ion '"//ion_name// &
                  "' in activity reaction.")

          ! Add the charge from this species
          key_name = "charge"
          call assert_msg(310667885, spec_props%get_int(key_name, charge), &
                  "Missing charge for ion '"//ion_name//"' in activity "// &
                  "reaction.")

          if (charge.gt.0) then
            _JACOB_NUM_CATION_(i_ion_pair) = qty
            _JACOB_CATION_MW_(i_ion_pair) = molecular_weight
         else if (charge.lt.0) then
            _JACOB_NUM_ANION_(i_ion_pair) = qty
            _JACOB_ANION_MW_(i_ion_pair) = molecular_weight
          else
            call die_msg(899416917, "Neutral species '"//ion_name// &
                    "' not allowed in activity reaction ion pair")
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
              if (charge.gt.1) then
                if (i_phase.eq.1) then
                  _JACOB_CATION_ID_(i_ion_pair) = &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          _PHASE_ID_(i_phase)
                else
                  call assert(473680545, _JACOB_CATION_ID_(i_ion_pair).eq. &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          _PHASE_ID_(i_phase))
                end if
              else
                if (i_phase.eq.1) then
                  _JACOB_ANION_ID_(i_ion_pair) = &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          _PHASE_ID_(i_phase)
                else
                  call assert(234155524, _JACOB_ANION_ID_(i_ion_pair).eq. &
                          aero_rep(i_aero_rep)%val%spec_state_id( &
                          unique_spec_names(i_spec)%string) - &
                          _PHASE_ID_(i_phase))
                end if
              end if
              i_phase = i_phase + 1
            end do
          end do
          i_phase = i_phase - 1
          call assert_msg(623684811, i_phase.eq._NUM_PHASE_, &
            "Incorrect number of instances of ion species '"// &
            ion_name//"' in activity reaction. Expected "// &
            trim(to_string(_NUM_PHASE_))//" but got "// &
            trim(to_string(i_phase)))

          ! Get the next ion
          call ions%iter_next()

        end do
    
        call assert_msg(319151390, total_charge.eq.0, &
                "Charge imbalance for ion pair '"//ion_pair_name// &
                " in activity reaction. Total charge: "// &
                trim(to_string(total_charge)))

        n_float_param = n_float_param + 3 + _JACOB_NUM_Y_(i_ion_pair)
        n_int_param = n_int_param + 6

      else if (str_type.eq."EQSAM") then

        ! Set the type
        _TYPE_(i_ion_pair) = ACT_CALC_EQSAM

        ! Get the required parameters for calculating activity
        key_name = "NW"
        call assert_msg(692339107, &
                ion_pair%get_real(key_name, _EQSAM_NW_(i_ion_pair)), &
                "Missing parameter NW for ion pair '"//ion_pair_name// &
                "' in activity reaction.")

        key_name = "ZW"
        call assert_msg(894917625, &
                ion_pair%get_real(key_name, _EQSAM_ZW_(i_ion_pair)), &
                "Missing parameter ZW for ion pair '"//ion_pair_name// &
                "' in activity reaction.")

        key_name = "MW"
        call assert_msg(272128568, &
                ion_pair%get_real(key_name, _EQSAM_ION_PAIR_MW_(i_ion_pair)), &
                "Missing parameter MW for ion pair '"//ion_pair_name// &
                "' in activity reaction.")

        ! Get the number and id of ions
        key_name = "ions"
        call assert(381088140, ion_pair%get_property_t(key_name, ions))
        _EQSAM_NUM_ION_(i_ion_pair) = ions%size()
        call ions%iter_reset()
        do i_ion = 1, ions%size()
          
          ! Get the ion name
          call assert(849711956, ions%get_key(ion_name))

          ! Get the species properties
          call assert_msg(826137761, &
                  chem_spec_data%get_property_set(ion_name, spec_props), &
                  "Missing species properties for ion '"//ion_name// &
                  "' in activity reaction.")
          
          ! Add the molecular weight
          key_name = "molecular weight"
          call assert_msg(598142298, &
                  spec_props%get_real(key_name, molecular_weight), &
                  "Missing molecular weight for ion '"//ion_name// &
                  "' in activity reaction.")
          _EQSAM_ION_MW_(i_ion_pair, i_ion) = molecular_weight

          ! Set the ion id (relative to water within the specified phase)
          i_phase = 1
          do i_aero_rep = 1, size(aero_rep)
            unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
                    phase_name = phase_name, spec_name = ion_name)
            if (.not.allocated(unique_spec_names)) cycle
            do i_spec = 1, size(unique_spec_names)
              if (i_phase.eq.1) then
                _EQSAM_ION_ID_(i_ion_pair,i_ion) = &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        _PHASE_ID_(i_phase)
              else
                call assert(973648240, _EQSAM_ION_ID_(i_ion_pair,i_ion) .eq. &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        _PHASE_ID_(i_phase))
              end if
              i_phase = i_phase + 1
            end do
          end do
          i_phase = i_phase - 1
          call assert_msg(900921350, i_phase.eq._NUM_PHASE_, &
            "Incorrect number of instances of ion species '"// &
            ion_name//"' in activity reaction. Expected "// &
            trim(to_string(_NUM_PHASE_))//" but got "// &
            trim(to_string(i_phase)))

           ! Get the next ion
           call ions%iter_next()

        end do

        n_float_param = n_float_param + 3 + _EQSAM_NUM_ION_(i_ion_pair)
        n_int_param = n_int_param + 2 + _EQSAM_NUM_ION_(i_ion_pair)

      else
        call die_msg(186680407, "Internal error.")
      end if

      call ion_pairs%iter_next()
    end do

    call assert(859412771, n_int_param.eq._TOTAL_INT_PARAM_)
    call assert(568314442, n_float_param.eq._TOTAL_FLOAT_PARAM_)
      
  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rate constant expression
  function build_rate_const_expr(this, rxn_id) result (expr)

    !> Rate constant expression
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_activity_t), intent(in) :: this
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
    class(rxn_activity_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef ACT_TYPE_JACOBSON
#undef ACT_TYPE_EQSAM

#undef _NUM_PHASE_
#undef _GAS_WATER_ID_
#undef _NUM_ION_PAIR_
#undef _TOTAL_INT_PARAM_
#undef _TOTAL_FLOAT_PARAM_
#undef _ppm_TO_RH_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _PHASE_ID_
#undef _PAIR_INT_PARAM_LOC_
#undef _PAIR_FLOAT_PARAM_LOC_
#undef _TYPE_
#undef _JACOB_NUM_CATION_
#undef _JACOB_NUM_ANION_
#undef _JACOB_CATION_ID_
#undef _JACOB_ANION_ID_
#undef _JACOB_NUM_Y_
#undef _EQSAM_NUM_ION_
#undef _EQSAM_ION_ID_
#undef _JACOB_low_RH_
#undef _JACOB_CATION_MW_
#undef _JACOB_ANION_MW_
#undef _JACOB_Y_
#undef _EQSAM_NW_
#undef _EQSAM_ZW_
#undef _EQSAM_ION_PAIR_MW_
#undef _EQSAM_ION_MW_

end module pmc_rxn_activity
