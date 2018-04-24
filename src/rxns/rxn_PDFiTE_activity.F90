! Copyright (C) 2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_rxn_PDFiTE_activity module.

!> \page phlex_rxn_PDFiTE_activity Phlexible Mechanism for Chemistry: PDFiTE Activity Reaction
!!
!! PDFiTE activity reactions calculate aerosol-phase species activities using
!! Taylor series to describe partial derivatives of mean activity coefficients
!! for ternary solutions, as described in \cite{Topping2009}. Thus, the mean
!! binary activity coefficients for electrolytes are calculated according to
!! eq. 15 in \cite{Topping2009}. The values are then available to aqueous-
!! phase reactions during solving.
!!
!! Input data for PDFiTE activity equations should take the form :
!! \code{.json}
!!   {
!!     "type" : "PDFITE_ACTIVITY",
!!     "gas-phase water" : "H2O",
!!     "aerosol-phase water" : "H2O_aq",
!!     "aerosol phase" : "my aero phase",
!!     "electrolytes" : {
!!       "Na2SO4_aq" : {
!!         "ions" : {
!!           "Na_p" : { "qty" : 2 },
!!           "SO4_mm" : {}
!!         },
!!         "interactions" : {
!!           "Na2SO4_aq" : [3.188349e4, -1.305168e5, 2.935608e5],
!!           "H2SO4_aq" : [-3.295311e3, 3.188349e4, -1.305168e5, 2.935608e5],
!!           "HHSO4_aq" : [412.423566, 1.412e-3, 2.421e7]
!!         }
!!       },
!!       "H2SO4_aq" : {
!!         "ions" : {
!!           "H_pp" : { "qty" : 2},
!!           "SO4_mm" : {}
!!         },
!!         "interactions" : {
!!           "H2SO4_aq" : [4.3124e3, -1.35e2, -32.6254e7, 32.4120],
!!           "Na2SO4_aq" : [-1.35e2, -32.6254e7, 32.4120],
!!           "HHSO4_aq" : [-1.493857e4, 3.50193, -3.401929858e3, 2.92345023]
!!         }
!!       },
!!       "HHSO4_aq" : {
!!         "ions" : {
!!           "H_pp" : {},
!!           "HSO4m" : {}
!!         },
!!         "interactions" : {
!!           "HHSO4_aq" : [-1.35e2, -32.6254e7, 32.4120],
!!           "Na2SO4_aq" : [31.24205, -1.35e2, -32.6254e7, 32.4120],
!!           "H2SO4_aq" : [-2.312412,  3.50193, -3.401929858e3, 2.92345023]
!!         }
!!       }
!!       ...
!!     }
!!   }
!! \endcode
!! The key-value pair \b "aerosol phase" is required to specify the aerosol
!! phase for which to calculate activity coefficients. The key-value pairs
!! \b "gas-phase water" and \b "aerosol-phase water" must also be present
!! and specify the names for the water species in each phase. The final
!! required key-value pair is \b "electrolytes", which should contain a set of
!! key-value pairs where the key of each member of the set is the name of a
!! binary electrolyte and the contents contain parameters required to estimate
!! their activity coefficients. The name of the electrolyte must refer to an
!! actual species in the specified aerosol phase.
!!
!! Each binary electrolyte must include a set of exactly two \b "ions" that
!! must also be present in the specified aerosol phase, must have a non-zero
!! \b "charge" parameter and a \b "molecular weight". If a \b "qty" is not
!! specified, it is assumed to be 1. The net charge must be zero.
!!
!! Each binary electrolyte must also have an \b "interactions" parameter that
!! has exactly one entry for every other binary electrolyte in the the
!! PDFiTE_activity reaction, where the key is the an electrolyte name, and the
!! value is an array containing at least one polynomial coefficient. The list
!! must also include exactly one entry with a key name of the electrolyte
!! itself where the parameters can be used to calculate \f$\gamma_x_0\f$
!! (Table 5 in  \cite{Topping2009}).
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
!!   "name" : "Na_p",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : 1,
!!   "molecular weight" : 22.9898
!! },
!! {
!!   "name" : "SO4_mm",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : -2
!!   "molecular weight" : 96.06
!! },
!! {
!!   "name" : "HSO4_m",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "charge" : -1
!!   "molecular weight" : 97.06
!! },
!! {
!!   "name" : "H2SO4_aq",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "molecular weight" : 98.06
!! },
!! {
!!   "name" : "HHSO4_aq",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "note" : "this species is only present for activity calculations"
!!   "molecular weight" : 98.06
!! },
!! {
!!   "name" : "Na2SO4_aq",
!!   "type" : "CHEM_SPEC",
!!   "phase" : "AEROSOL",
!!   "molecular weight" : 142.04
!! },
!! {
!!   "name" : "my aero phase",
!!   "type" : "AERO_PHASE",
!!   "species" : ["Na_p", "SO4_mm", "HSO4_m", "Na2SO4_aq", "H2SO4_aq", "HHSO4_aq", "H2O_aq"]
!! }
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The rxn_PDFiTE_activity_t type and associated functions. 
module pmc_rxn_PDFiTE_activity

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

#define _NUM_PHASE_ this%condensed_data_int(1)
#define _GAS_WATER_ID_ this%condensed_data_int(2)
#define _NUM_ELECTROLYTES_ this%condensed_data_int(3)
#define _TOTAL_INT_PARAM_ this%condensed_data_int(4)
#define _TOTAL_FLOAT_PARAM_ this%condensed_data_int(5)
#define _ppm_TO_RH_ this%condensed_data_real(1)
#define _NUM_INT_PROP_ 5
#define _NUM_REAL_PROP_ 1
#define _PHASE_ID_(x) this%condensed_data_int(_NUM_INT_PROP_+x)
#define _ELECT_INT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_PHASE_+x)
#define _ELECT_FLOAT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_PHASE_+_NUM_ELECTROLYTES_+x)
#define _ELECTROLYTE_ID_(x) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x))
#define _ELECTROLYTE_ACT_ID_(x) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+1)
#define _NUM_CATION_(x) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+2)
#define _NUM_ANION_(x) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+3)
#define _CATION_ID_(x) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+4)
#define _ANION_ID_(x) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+5)
#define _NUM_B_(x,y) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+5+y)
#define _INTER_SPEC_ID_(x,y) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+5+_NUM_ELECTROLYTES_+y)
#define _B0_LOC_(x,y) this%condensed_data_int(_ELECT_INT_PARAM_LOC_(x)+5+2*(_NUM_ELECTROLYTES_)+y)
#define _ELECTROLYTE_MW_(x) this%condensed_data_real(_ELECT_FLOAT_PARAM_LOC_(x))
#define _CATION_MW_(x) this%condensed_data_real(_ELECT_FLOAT_PARAM_LOC_(x)+1)
#define _ANION_MW_(x) this%condensed_data_real(_ELECT_FLOAT_PARAM_LOC_(x)+2)
#define _CATION_N_(x) this%condensed_data_real(_ELECT_FLOAT_PARAM_LOC_(x)+3)
#define _ANION_N_(x) this%condensed_data_real(_ELECT_FLOAT_PARAM_LOC_(x)+4)
#define _B_z_(x,y,z) this%condensed_data_real(_B0_LOC_(x,y)-1+z)

  public :: rxn_PDFiTE_activity_t

  !> Generic test reaction data type
  type, extends(rxn_data_t) :: rxn_PDFiTE_activity_t
  contains
    !> Reaction initialization
    procedure :: initialize
    !> Build rate constant expression
    procedure :: build_rate_const_expr
    !> Build time derivative expression
    procedure :: build_deriv_expr
  end type rxn_PDFiTE_activity_t

  !> Constructor for rxn_PDFiTE_activity_t
  interface rxn_PDFiTE_activity_t
    procedure :: constructor
  end interface rxn_PDFiTE_activity_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for activity reaction
  function constructor() result(new_obj)

    !> A new reaction instance
    type(rxn_PDFiTE_activity_t), pointer :: new_obj

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
    class(rxn_PDFiTE_activity_t), intent(inout) :: this
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data
    !> Aerosol representations
    class(aero_rep_data_ptr), pointer, intent(in) :: aero_rep(:)

    type(property_t), pointer :: spec_props, electrolytes, electrolyte, &
            sub_props, ions, interactions
    character(len=:), allocatable :: key_name, water_name, spec_name, &
            phase_name, string_val, inter_spec_name
    integer(kind=i_kind) :: n_phase, n_electrolyte, n_int_param, n_float_param
    integer(kind=i_kind) :: i_aero_rep, i_phase, i_electrolyte, i_ion, i_spec, &
            i_sub_prop, i_interaction
    integer(kind=i_kind) :: qty, int_val, charge, total_charge
    real(kind=dp) :: real_val, molecular_weight
    class(string_t), allocatable :: unique_spec_names(:)
    character(len=:), allocatable :: electrolyte_name, ion_name
    type(string_t), allocatable :: electrolyte_names(:)
    integer(kind=i_kind), allocatable :: num_param(:)

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
    do i_aero_rep = 1, size(aero_rep)

      ! Get the number of instances of the phase in this representation
      n_phase = n_phase + &
              aero_rep(i_aero_rep)%val%num_phase_instances(phase_name)

    end do

    call assert_msg(656403972, n_phase.gt.0, &
            "Aerosol phase '"//phase_name//"' not present in any aerosol "// &
            "representation for PDFiTE activity reaction.")

    ! Get the electrolytes
    key_name = "electrolytes"
    call assert_msg(258251605, &
            this%property_set%get_property_t(key_name, electrolytes), &
            "Missing electrolytess in PDFiTE activity reaction.")

    ! Count the electrolytess
    n_electrolyte = electrolytes%size()
    call assert_msg(479529522, n_electrolyte .gt. 0, &
            "Empty ion pair set in PDFiTE activity reaction.")

    ! Allocate space for the electrolyte names and array to check for their
    ! interaction parameters
    allocate(electrolyte_names(n_electrolyte))
    allocate(num_param(n_electrolyte))

    ! Get the number of parameters required for each ion pair
    n_int_param = _NUM_INT_PROP_ + n_phase + 2*n_electrolyte
    n_float_param = _NUM_REAL_PROP_ + 5*n_electrolyte
    call electrolytes%iter_reset()
    do i_electrolyte = 1, n_electrolyte
  
      ! Get the name of the electrolyte
      call assert(680654801, electrolytes%get_key(electrolyte_name))
      electrolyte_names(i_electrolyte)%string = electrolyte_name

      ! Get the electrolyte properties
      call assert_msg(287766741, electrolytes%get_property_t(val=electrolyte), &
              "Missing electrolyte properties for '"//electrolyte_name// &
              "' in PDFiTE activity reaction.")

      ! Get the interactions
      key_name = "interactions"
      call assert_msg(883233319, &
              electrolyte%get_property_t(key_name, interactions), &
              "Missing interaction parameters for '"//electrolyte_name// &
              "' in PDFiTE activity reaction.")

      ! Check the number of interactions
      call assert_msg(977304560, &
              interactions%size().eq.n_electrolyte, &
              "Incorrect number of interaction parameters for species '"// &
              electrolyte_name//"' in PDFiTE activity reaction. Expected "// &
              to_string(n_electrolyte)//" but got "//to_string(interactions%size()))

      call interactions%iter_reset()
      do i_interaction = 1, interactions%size()

        ! Get the number of B_z parameters
        call assert_msg(552282009, &
                interactions%get_property_t(val=sub_props), &
                "Missing polynomial coefficients for PDFiTE activity calculation "// &
                "for electrolyte '"//electrolyte_name//"'.")

        call assert_msg(320927773, sub_props%size().gt.0, &
                "Insufficient polynomial coefficients for PDFiTE activity "// &
                "calculation for electrolyte '"//electrolyte_name//"'.")

        n_float_param = n_float_param + sub_props%size()
        n_int_param = n_int_param + 3

        call interactions%iter_next()
      end do

      n_int_param = n_int_param + 6

      call electrolytes%iter_next()
    end do

    ! Allocate space in the condensed data arrays
    allocate(this%condensed_data_int(n_int_param))
    allocate(this%condensed_data_real(n_float_param))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set some data dimensions
    _NUM_PHASE_  = n_phase
    _NUM_ELECTROLYTES_ = n_electrolyte
    _TOTAL_INT_PARAM_ = n_int_param
    _TOTAL_FLOAT_PARAM_ = n_float_param

    ! Set the gas-phase water species
    key_name = "gas-phase water"
    call assert_msg(901506020, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing gas-phase water species name in PDFiTE activity "// &
            "reaction.")

    _GAS_WATER_ID_ = chem_spec_data%gas_state_id(spec_name)
    
    call assert_msg(442608616, _GAS_WATER_ID_ .gt. 0, &
            "Cannot find gas-phase water species '"//spec_name//"' for "// &
            "PDFiTE activity reaction.")
    
    ! Set the aerosol-water species
    key_name = "aerosol-phase water"
    call assert_msg(214613153, &
            this%property_set%get_string(key_name, spec_name), &
            "Missing aerosol-phase water species name in PDFiTE activity "// &
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
        call assert(658622226, _PHASE_ID_(i_phase).gt.0)
        i_phase = i_phase + 1
      end do
    end do
    i_phase = i_phase - 1
    call assert_msg(653357919, i_phase.eq._NUM_PHASE_, &
            "Incorrect number of aerosol water instances in PDFiTE "// &
            "activity reaction. Expected "//trim(to_string(_NUM_PHASE_))// &
            " but got "//trim(to_string(i_phase)))

    ! Save the electrolyte parameters
    n_int_param = _NUM_INT_PROP_ + _NUM_PHASE_ + 2*_NUM_ELECTROLYTES_
    n_float_param = _NUM_REAL_PROP_
    call electrolytes%iter_reset()
    do i_electrolyte = 1, n_electrolyte
   
      ! Get the name of the electrolyte
      call assert(806720958, electrolytes%get_key(electrolyte_name))

      ! Get the electrolyte properties
      call assert(520886936, electrolytes%get_property_t(val=electrolyte))

      ! Set the location of this electrolyte's parameters in the condensed data
      ! arrays.
      _ELECT_INT_PARAM_LOC_(i_electrolyte) = n_int_param + 1
      _ELECT_FLOAT_PARAM_LOC_(i_electrolyte) = n_float_param + 1

      ! Get the state ids for this electrolyte
      i_phase = 1
      do i_aero_rep = 1, size(aero_rep)
        unique_spec_names = aero_rep(i_aero_rep)%val%unique_names( &
                phase_name = phase_name, spec_name = electrolyte_name)
        if (.not.allocated(unique_spec_names)) cycle
        do i_spec = 1, size(unique_spec_names)
          if (i_phase.eq.1) then
            _ELECTROLYTE_ID_(i_electrolyte) = &
                    aero_rep(i_aero_rep)%val%spec_state_id( &
                    unique_spec_names(i_spec)%string) - &
                    _PHASE_ID_(i_phase)
          else
            call assert(142173386, _ELECTROLYTE_ID_(i_electrolyte).eq. &
                    aero_rep(i_aero_rep)%val%spec_state_id( &
                    unique_spec_names(i_spec)%string) - &
                    _PHASE_ID_(i_phase))
          end if
          i_phase = i_phase + 1
        end do
      end do
      i_phase = i_phase - 1
      call assert_msg(700406338, i_phase.eq._NUM_PHASE_, &
              "Incorrect number of instances of ion species '"// &
              ion_name//"' in PDFiTE activity reaction. Expected "// &
              trim(to_string(_NUM_PHASE_))//" but got "// &
              trim(to_string(i_phase)))
        
      ! Get the electrolyte species properties
      call assert_msg(737691158, &
              chem_spec_data%get_property_set(electrolyte_name, spec_props), &
              "Missing species properties for electrolyte '"// &
              electrolyte_name//"' in PDFiTE activity reaction.")

      ! Set the electrolyte molecular weight
      key_name = "molecular weight"
      call assert_msg(380583485, &
              spec_props%get_real(key_name, _ELECTROLYTE_MW_(i_electrolyte)), &
              "Missing molecular weight for electrolyte '"// &
              electrolyte_name//"' in PDFiTE activity reaction.")

      ! Get the number and id of the ions and their molecular weights
      key_name = "ions"
      call assert_msg(974596824, electrolyte%get_property_t(key_name, ions), &
              "Mission ions for electrolyte '"//electrolyte_name// &
              "' in PDFiTE activity reaction.")

      call assert_msg(852202160, ions%size().eq.2, &
              "Invalid number of unique ions specified for electrolyte '"// &
              electrolyte_name//"' in for PDFiTE in activity reaction. "// &
              "Expected 2 got "//trim(to_string(ions%size())))

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
        key_name = "molecular weight"
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
          _NUM_CATION_(i_electrolyte) = qty
          _CATION_MW_(i_electrolyte) = molecular_weight
        else if (charge.lt.0) then
          _NUM_ANION_(i_electrolyte) = qty
          _ANION_MW_(i_electrolyte) = molecular_weight
        else
          call die_msg(939555855, "Neutral species '"//ion_name// &
                  "' not allowed in PDFiTE activity reaction ion pair")
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
                _CATION_ID_(i_electrolyte) = &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        _PHASE_ID_(i_phase)
              else
                call assert(425726370, _CATION_ID_(i_electrolyte).eq. &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        _PHASE_ID_(i_phase))
              end if
            else
              if (i_phase.eq.1) then
                _ANION_ID_(i_electrolyte) = &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        _PHASE_ID_(i_phase)
              else
                call assert(192466600, _ANION_ID_(i_electrolyte).eq. &
                        aero_rep(i_aero_rep)%val%spec_state_id( &
                        unique_spec_names(i_spec)%string) - &
                        _PHASE_ID_(i_phase))
              end if
            end if
            i_phase = i_phase + 1
          end do
        end do
        i_phase = i_phase - 1
        call assert_msg(759322632, i_phase.eq._NUM_PHASE_, &
                "Incorrect number of instances of ion species '"// &
                ion_name//"' in PDFiTE activity reaction. Expected "// &
                trim(to_string(_NUM_PHASE_))//" but got "// &
                trim(to_string(i_phase)))

        ! Get the next ion
        call ions%iter_next()

      end do
    
      call assert_msg(415650051, total_charge.eq.0, &
              "Charge imbalance for electrolyte '"//electrolyte_name// &
              " in PDFiTE activity reaction. Total charge: "// &
              trim(to_string(total_charge)))

      ! Get the interaction parameters
      n_float_param = n_float_param + 5
      num_param(:) = 0
      call interactions%iter_reset()
      do i_interaction = 1, interactions%size()

        ! Set the location of the B_z parameters
        _B0_LOC_(i_electrolyte, i_interaction) = n_float_param + 1

        ! Get the number of B_z parameters
        call assert(329488069, interactions%get_property_t(val=sub_props))
        _NUM_B_(i_electrolyte, i_interaction) = sub_props%size()
          
        ! Get the name of the interacting species
        call assert_msg(220815620, sub_props%get_key(inter_spec_name), &
                "Missing interacting species name for electrolyte '"// &
                electrolyte_name//"' in PDFiTE activity reaction.")

        ! Mark the species as having parameters and get the id of the
        ! electrolyte
        do i_spec = 1, size(electrolyte_names)
          if (electrolyte_names(i_spec)%string.eq.inter_spec_name) then
            num_param(i_spec) = num_param(i_spec) + 1
            _INTER_SPEC_ID_(i_electrolyte,i_interaction) = i_spec
            exit
          end if
        end do

        ! Get the B_z parameters
        call sub_props%iter_reset()
        do i_sub_prop = 1, sub_props%size()
          call assert_msg(988796931, sub_props%get_real(val=real_val), &
                  "Invalid polynomial coefficient for electrolyte '"// &
                  electrolyte_name//"' in PDFiTE activity reaction.")
          _B_z_(i_electrolyte, i_interaction, i_sub_prop) = real_val
          call sub_props%iter_next()
        end do
  
        n_float_param = n_float_param + _NUM_B_(i_electrolyte, i_interaction)
      
        call interactions%iter_next()
      end do
      
      ! Check that all the interactions were included exactly once
      do i_spec = 1, size(num_param)
        call assert_msg(793223082, num_param(i_spec).eq.1, &
                "Incorrect number of interaction parameters between "// &
                "electrolyte '"//electrolyte_name//"' and '"// &
                electrolyte_names(i_spec)%string//"' in PDFiTE activity "// &
                "reaction: "//to_string(num_param(i_spec))//".")
      end do

      n_int_param = n_int_param + 5 + 3*(_NUM_ELECTROLYTES_)

      call electrolytes%iter_next()
    end do

    call assert(938415336, n_int_param.eq._TOTAL_INT_PARAM_)
    call assert(433208931, n_float_param.eq._TOTAL_FLOAT_PARAM_)
      
  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rate constant expression
  function build_rate_const_expr(this, rxn_id) result (expr)

    !> Rate constant expression
    character(len=:), allocatable :: expr
    !> Reaction data
    class(rxn_PDFiTE_activity_t), intent(in) :: this
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
    class(rxn_PDFiTE_activity_t), intent(in) :: this
    !> Reaction id in mechanism
    integer(kind=i_kind), intent(in) :: rxn_id
    !> Species id to get contribution for
    integer(kind=i_kind), intent(in) :: spec_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    expr = ""

  end function build_deriv_expr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef _NUM_PHASE_
#undef _GAS_WATER_ID_
#undef _NUM_ELECTROLYTES_
#undef _TOTAL_INT_PARAM_
#undef _TOTAL_FLOAT_PARAM_
#undef _ppm_TO_RH_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _PHASE_ID_
#undef _ELECT_INT_PARAM_LOC_
#undef _ELECT_FLOAT_PARAM_LOC_
#undef _ELECTROLYTE_ID_
#undef _ELECTROLYTE_ACT_ID_
#undef _NUM_CATION_
#undef _NUM_ANION_
#undef _CATION_ID_
#undef _ANION_ID_
#undef _NUM_B_
#undef _INTER_SPEC_ID_
#undef _B0_LOC_
#undef _ELECTROLYTE_MW_
#undef _CATION_MW_
#undef _ANION_MW_
#undef _CATION_N_
#undef _ANION_N_
#undef _B_z_

end module pmc_rxn_PDFiTE_activity
