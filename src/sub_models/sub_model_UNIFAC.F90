! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_sub_model_UNIFAC module.

!> \page phlex_sub_model_UNIFAC Phlexible Module for Chemistry: UNIFAC Activity Coefficients
!!
!! The UNIFAC activity coefficient sub model calculates activity coefficients
!! for species in an aerosol phase based on the current aerosol phase
!! composition. The \c json object for this \ref phlex_sub_model "sub model"
!! is of the form:
!! \code{.json}
!!  { "pmc-data" : [
!!    {
!!      "type" : "SUB_MODEL_UNIFAC",
!!      "phases" : [
!!        "some phase",
!!        "some other phase"
!!      ],
!!      "functional groups" : {
!!        "CH3" : {
!!          "main group" : "CH2",
!!          "volume param" : 0.9011,
!!          "surface param" : 0.8480
!!        },
!!        "CH2" : {
!!          "main group" : "CH2",
!!          "volume param" : 0.6744,
!!          "surface param" : 0.5400
!!        },
!!        "CH=CH" : {
!!          "main group" : "C=C",
!!          "volume param" : 1.1167,
!!          "suface param" : 0.8670
!!        }
!!      },
!!      "main groups" : {
!!        "CH2" : {
!!          "interactions with" : {
!!            "C=C" : -35.36
!!          }
!!        },
!!        "C=C" : {
!!          "interactions with" : {
!!            "CH2" : 86.02
!!          }
!!        }
!!      }
!!    },
!!    ...
!!  ]}
!! \endcode
!! The key-value pair \b type is required and must be \b SUB_MODEL_UNIFAC.
!! The key-value pair \b phases is also required, and its value must be an
!! array of strings that correspond to valid 
!! \ref phlex_aero_phase "aerosol phases". The key-value pair \b "functional
!! groups" is also required, and must contain a set of key-value pairs whose
!! keys are the names of UNIFAC functions groups, and whose values are a set
!! of key value pairs that contain, at minimum:
!!   - \b "main group" : a string that corresponds to a key in the \b
!!                       "main groups" set.
!!   - \b "volume param" : the floating-point volume parameter for this 
!!                         functional group.
!!   - \b "surface param" : this floating-point surface parameter for this
!!                          functional group.
!! The last required key-value pair is \b "main groups" whose value must
!! be a set of key-value pairs whose keys are the names of the UNIFAC main
!! groups and whose values are a set key-pairs that contain, at minimum,
!! \b "interaction with" whose value is a set of key-value pairs whose keys
!! are the names of the other \b "main groups" and whose values are the
!! floating-point interation parameters for that interaction. Each main group
!! may contain up to one interaction with each other main group, and may
!! not contain an interaction with itself. Missing interactions are assumed
!! to be 0.0.
!!
!! Species in the specified phase for whom acitivity coefficients will be
!! calculated must contain a key-value pair \b "UNIFAC groups" whose value
!! is a set of key value pairs that correspond with members of the
!! \b "functional groups" set and whose values are the integer number of
!! instances of a particular functional group in this species. For the
!! above example UNIFAC model, the following species would be valid and
!! included in activity coefficient calculations:
!! \code{.json}
!! { "pmc-data" : [
!!   {
!!     "name" : "my species",
!!     "type" : "CHEM_SPEC",
!!     "phase" : "AEROSOL",
!!     "UNIFAC groups" : {
!!       "CH3" : 4,
!!       "C=C" : 1
!!     }
!!   },
!!   {
!!     "name" : "my other species",
!!     "type" : "CHEM_SPEC",
!!     "phase" : "AEROSOL",
!!     "UNIFAC groups" : {
!!       "CH3" : 2,
!!       "CH2" : 4
!!     },
!!   },
!!   {
!!     "name" : "some phase",
!!     "type" : "AERO_PHASE",
!!     "species" : { "my species", "my other species" }
!!   }
!! ]}
!! \endcode

!> The sub_model_UNIFAC_t type and assocatiated subroutines
!!
!! This module is based on the UNIFAC module of Alf Grini, CNRM, 2005, which
!! in turn, was based original code received from Pierre Tulet, who got it
!! from Betty Pun who (it seems) got it from Pradeep Saxena. The UNIFAC
!! module is part of the Model to Predict the Multi-phase Partitioning of
!! Organics (Griffif et al., JGR 110, D05304, 2005 doi: 10.1029/2004JD005219)
!!
!! Equations referenced are from Marcolli and Peter, ACP 5(2), 1501-1527, 
!! 2005, and variable names roughly follow their naming scheme.
!!
module pmc_sub_model_UNIFAC

  use pmc_util,                                 only : dp, i_kind, &
                                                       string_t, assert_msg, &
                                                       die_msg, to_string, &
                                                       assert
  use pmc_property
  use pmc_sub_model_data
  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_aero_phase_data
  use pmc_phlex_state

  implicit none
  private

#define _NUM_UNIQUE_PHASE_ this%condensed_data_int(1)
#define _NUM_GROUP_ this%condensed_data_int(2)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 0
#define _PHASE_INT_LOC_(p) this%condensed_data_int(_NUM_INT_PROP_+p)
#define _PHASE_REAL_LOC_(p) this%condensed_data_int(_NUM_INT_PROP_+_NUM_UNIQUE_PHASE_+p)
#define _NUM_PHASE_INSTANCE_(p) this%condensed_data_int(_PHASE_INT_LOC_(p))
#define _NUM_SPEC_(p) this%condensed_data_int(_PHASE_INT_LOC_(p)+1)
#define _PHASE_INST_REAL_LOC_(p,c) this%condensed_data_int(_PHASE_INT_LOC_(p)+1+c)
#define _PHASE_INST_ID_(p,c) this%condensed_data_int(_PHASE_INT_LOC_(p)+1+_NUM_PHASE_INSTANCE_(p)+c)
#define _SPEC_ID_(p,i) this%condensed_data_int(_PHASE_INT_LOC_(p)+1+2*_NUM_PHASE_INSTANCE_(p)+i)
#define _v_ik_(p,i,k) this%condensed_data_int(_PHASE_INT_LOC_(p)+1+2*_NUM_PHASE_INSTANCE_(p)+k*_NUM_SPEC_(p)+i)

#define _Q_k_(k) this%condensed_data_real(k)
#define _R_k_(k) this%condensed_data_real(_NUM_GROUP_+k)
#define _a_mn_(m,n) this%condensed_data_real((m+1)*_NUM_GROUP_+n)
#define _PSI_mn_(m,n) this%condensed_data_real((m+1+_NUM_GROUP_)*_NUM_GROUP_+n)
#define _r_i_(p,i) this%condensed_data_real(_PHASE_REAL_LOC_(p)+i-1)
#define _q_i_(p,i) this%condensed_data_real(_PHASE_REAL_LOC_(p)+_NUM_SPEC_(p)+i-1)
#define _l_i_(p,i) this%condensed_data_real(_PHASE_REAL_LOC_(p)+2*_NUM_SPEC_(p)+i-1)
#define _ln_GAMMA_ik_(p,i,k) this%condensed_data_real(_PHASE_REAL_LOC_(p)+(i-1)*_NUM_GROUP_+3*_NUM_SPEC_(p)+k-1)
#define _gamma_i_(p,c,i) this%condensed_data_real(_PHASE_INST_REAL_LOC_(p,c)+i-1)

  ! Update types (These must match values in sub_model_UNIFAC.c)
  ! (none for now)

  public :: sub_model_UNIFAC_t

  !> UNIFAC activity coefficient calculation
  !!
  !! Time-invariant data required by the UNIFAC activity coefficient sub model
  type, extends(sub_model_data_t) :: sub_model_UNIFAC_t
  contains
    !> Initialize the sub model data, validating input parameters and
    !! loading any required information form the \c
    !! sub_model_data_t::property_set. This routine should be called
    !! once for each sub model at the beginning of the model run after all
    !! the input files have been read in. It ensures all data required
    !! during the model run are included in the condensed data arrays.
    procedure :: initialize
  end type sub_model_UNIFAC_t

  ! Constructor for sub_model_UNIFAC_t
  interface sub_model_UNIFAC_t
    procedure :: constructor
  end interface sub_model_UNIFAC_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for sub_model_UNIFAC_t
  function constructor() result (new_obj)

    !> New sub model
    type(sub_model_UNIFAC_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  !> Initialize the sub model data, validating input parameters and
  !! loading any required information form the \c
  !! sub_model_data_t::property_set. This routine should be called
  !! once for each sub model at the beginning of the model run after all
  !! the input files have been read in. It ensures all data required
  !! during the model run are included in the condensed data arrays.
  subroutine initialize(this, aero_rep_set, aero_phase_set, chem_spec_data)

    !> Sub model data
    class(sub_model_UNIFAC_t), intent(inout) :: this
    !> The set of aerosol representations
    type(aero_rep_data_ptr), pointer, intent(in) :: aero_rep_set(:)
    !> The set of aerosol phases
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    type(property_t), pointer :: spec_props, phases
    type(property_t), pointer :: func_groups, func_group
    type(property_t), pointer :: main_groups, main_group
    type(property_t), pointer :: spec_groups, interactions
    character(len=:), allocatable :: key_name, spec_name, phase_name
    character(len=:), allocatable :: string_val, inter_group_name
    character(len=:), allocatable :: main_group_name, spec_group_name
    integer(kind=i_kind) :: i_spec, i_phase, i_rep, i_main_group, i_group
    integer(kind=i_kind) :: i_instance, i_inter, i_phase_inst, i_spec_group
    integer(kind=i_kind) :: i_inter_group
    integer(kind=i_kind) :: i_UNIFAC_phase
    integer(kind=i_kind) :: num_unique_phase, num_group, num_main_group
    integer(kind=i_kind) :: num_int_data, num_real_data, num_spec_group
    integer(kind=i_kind) :: curr_spec_id, curr_phase_inst_id
    integer(kind=i_kind) :: m, n
    integer(kind=i_kind), allocatable :: num_phase_inst(:)
    integer(kind=i_kind), allocatable :: num_phase_spec(:)
    integer(kind=i_kind), allocatable :: unique_phase_set_id(:)
    integer(kind=i_kind), allocatable :: main_group_id(:)
    type(string_t), allocatable :: phase_names(:)
    type(string_t), allocatable :: main_group_names(:)
    type(string_t), allocatable :: group_names(:)
    type(string_t), allocatable :: unique_names(:)
    real(kind=dp) :: q_i, r_i
    real(kind=dp) :: inter_param
    real(kind=dp), allocatable :: main_group_interactions(:,:)
    logical :: found, phase_ids_set

    ! Get the property set
    call assert_msg(403771584, associated(this%property_set), &
            "Missing property set needed to initialize UNIFAC model.")

    ! Get the aerosol phase names
    key_name = "phases"
    call assert_msg(114578789, &
            this%property_set%get_property_t(key_name, phases), &
            "Missing set of aerosol phase names for UNIFAC model.")
    num_unique_phase = phases%size()
    call assert_msg(780318074, num_unique_phase.gt.0, &
            "Received empty list of aerosol phase names for UNIFAC model.")

    ! Get the functional groups
    key_name = "functional groups"
    call assert_msg(372089388, &
            this%property_set%get_property_t(key_name, func_groups), &
            "Missing set of functional groups for UNIFAC model.")
    num_group = func_groups%size()
    call assert_msg(974837385, num_group.gt.0, &
            "Received empty set of functional groups for UNIFAC model.")

    ! Get the main groups
    key_name = "main groups"
    call assert_msg(956137986, &
            this%property_set%get_property_t(key_name, main_groups), &
            "Missing set of main groups for UNIFAC model.")
    num_main_group = main_groups%size()
    call assert_msg(556532380, num_main_group.gt.0, &
            "Received empty set of main groups for UNIFAC model.")

    ! Count the species in each phase, and the number of instances of each 
    ! phase
    allocate(num_phase_inst(num_unique_phase))
    allocate(num_phase_spec(num_unique_phase))
    allocate(unique_phase_set_id(num_unique_phase))
    allocate(phase_names(num_unique_phase))
    key_name = "UNIFAC groups"
    call phases%iter_reset()
    do i_UNIFAC_phase = 1, num_unique_phase
      call assert_msg(112836027, phases%get_string(val = phase_name), &
              "Received non-string phase name in UNIFAC model.")
      phase_names(i_UNIFAC_phase)%string = phase_name
      found = .false.
      do i_phase = 1, size(aero_phase_set)
        if (aero_phase_set(i_phase)%val%name().eq.phase_name) then
          found = .true.
          unique_phase_set_id(i_UNIFAC_phase) = i_phase
          num_phase_spec(i_UNIFAC_phase) = 0
          do i_spec = 1, aero_phase_set(i_phase)%val%size()
            spec_name = aero_phase_set(i_phase)%val%get_species_name(i_spec)
            call assert(698678581, &
                    chem_spec_data%get_property_set(spec_name, spec_props))
            if (spec_props%get_property_t(key_name, spec_groups)) then
              num_phase_spec(i_UNIFAC_phase) = &
                      num_phase_spec(i_UNIFAC_phase) + 1
            end if
          end do
        end if
      end do
      call assert_msg(835247755, found, "Cannot find aerosol phase '"// &
              phase_name//"' for UNIFAC model.")

      num_phase_inst(i_UNIFAC_phase) = 0
      do i_rep = 1, size(aero_rep_set)
        num_phase_inst(i_UNIFAC_phase) = num_phase_inst(i_UNIFAC_phase) + &
                aero_rep_set(i_rep)%val%num_phase_instances(phase_name)
      end do
      call assert_msg(187041753, num_phase_inst(i_UNIFAC_phase).gt.0, &
              "No instances of phase '"//phase_name//"' for UNIFAC model.")

      call phases%iter_next()
    end do

    ! Size of condensed data arrays
    num_int_data =   _NUM_INT_PROP_              & ! int props
                     + 2*num_unique_phase          ! PHASE_INT_LOC, PHASE_REAL_LOC
    num_real_data =  _NUM_REAL_PROP_             & ! real props
                     + 2*num_group               & ! Q_k, R_k
                     + 2*num_group*num_group       ! a_mn, PSI_mn
    do i_UNIFAC_phase = 1, num_unique_phase
      num_int_data = num_int_data + 2                    & ! NUM_PHASE_INSTANCE, NUM_SPEC
                     + 2*num_phase_inst(i_UNIFAC_phase)  & ! PHASE_INST_REAL_LOC, PHASE_INST_ID
                     + num_phase_spec(i_UNIFAC_phase)    & ! SPEC_ID
                     + num_phase_spec(i_UNIFAC_phase) * num_group          ! v_ik
      num_real_data = num_real_data &
                     + 3*num_phase_spec(i_UNIFAC_phase)  & ! r_i, q_i, l_i
                     + num_phase_spec(i_UNIFAC_phase) * num_group        & ! ln_GAMMA_ik
                     + num_phase_spec(i_UNIFAC_phase) * num_phase_inst(i_UNIFAC_phase)    ! gamma
    end do

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(num_int_data))
    allocate(this%condensed_data_real(num_real_data))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)

    ! Set sub model dimensions
    _NUM_UNIQUE_PHASE_ = num_unique_phase
    _NUM_GROUP_ = num_group
    
    ! Set data locations
    num_int_data =   _NUM_INT_PROP_              & ! int props
                     + 2*num_unique_phase          ! PHASE_INT_LOC, PHASE_REAL_LOC
    num_real_data =  _NUM_REAL_PROP_             & ! real props
                     + 2*num_group               & ! Q_k, R_k
                     + 2*num_group*num_group       ! a_mn, PSI_mn
    do i_UNIFAC_phase = 1, num_unique_phase
      _PHASE_INT_LOC_(i_UNIFAC_phase) = num_int_data + 1
      _PHASE_REAL_LOC_(i_UNIFAC_phase) = num_real_data + 1
      num_int_data = num_int_data + 2                    & ! NUM_PHASE_INSTANCE, NUM_SPEC
                     + 2*num_phase_inst(i_UNIFAC_phase)  & ! PHASE_INST_REAL_LOC, PHASE_INST_ID
                     + num_phase_spec(i_UNIFAC_phase)    & ! SPEC_ID
                     + num_phase_spec(i_UNIFAC_phase) * num_group          ! v_ik
      num_real_data = num_real_data &
                     + 3*num_phase_spec(i_UNIFAC_phase)  & ! r_i, q_i, l_i
                     + num_phase_spec(i_UNIFAC_phase) * num_group          ! ln_GAMMA_ik
      do i_phase_inst = 1, num_phase_inst(i_UNIFAC_phase)
        _PHASE_INST_REAL_LOC_(i_UNIFAC_phase, i_phase_inst) = num_real_data + 1
        num_real_data = num_real_data + num_phase_spec(i_UNIFAC_phase)     ! gamma
      end do
    end do

    ! Set phase dimensions
    do i_UNIFAC_phase = 1, _NUM_UNIQUE_PHASE_
      _NUM_SPEC_(i_UNIFAC_phase) = num_phase_spec(i_UNIFAC_phase)
      _NUM_PHASE_INSTANCE_(i_UNIFAC_phase) = num_phase_inst(i_UNIFAC_phase)
    end do

    ! Set the main group names
    allocate(main_group_names(main_groups%size()))
    call main_groups%iter_reset()
    do i_main_group = 1, main_groups%size()

      ! Save the main group name
      call assert(296339642, &
              main_groups%get_key(main_group_names(i_main_group)%string))

      call main_groups%iter_next()
    end do

    ! Set the main group interaction parameter matrix
    allocate(main_group_interactions(main_groups%size(), main_groups%size()))
    main_group_interactions(:,:) = real(0.0, kind=dp)
    call main_groups%iter_reset()
    do i_main_group = 1, main_groups%size()

      ! Get the main group properties
      call assert_msg(577361652, main_groups%get_property_t(val = main_group), &
              "Invalid main group '"//main_group_names(i_main_group)%string// &
              "' in UNIFAC model.")

      ! Get the interactions
      key_name = "interactions with"
      call assert_msg(208272126, &
              main_group%get_property_t(key_name, interactions), &
              "Missing interactions for main group '"// &
              main_group_names(i_main_group)%string//"' in UNIFAC model.")

      ! Set the interactions
      call interactions%iter_reset()
      do i_inter = 1, interactions%size()

        ! Get the interacting group name
        call assert(363540699, interactions%get_key(inter_group_name))

        ! Get the interaction parameter
        call assert_msg(976253437, interactions%get_real(val = inter_param), &
                "Invalid interaction parameter for interaction between "// &
                "main groups '"//main_group_names(i_main_group)%string// &
                "' and '"//trim(inter_group_name)//"' in UNIFAC model.")

        ! Set the interaction parameter
        do i_inter_group = 1, size(main_group_names)
          found = .false.
          if (main_group_names(i_inter_group)%string .eq. &
                  inter_group_name) then
            main_group_interactions(i_main_group, i_inter_group) = inter_param
            found = .true.
            exit
          end if
        end do
        call assert_msg(898262240, found, "Bad main group name '"// &
                inter_group_name//"' in interactions of '"// &
                main_group_names(i_main_group)%string//"' in UNIFAC model.")

        call interactions%iter_next()
      end do

      call main_groups%iter_next()
    end do

    ! Set the functional group parameters
    allocate(group_names(_NUM_GROUP_))
    allocate(main_group_id(_NUM_GROUP_))
    call func_groups%iter_reset()
    do i_group = 1, _NUM_GROUP_

      ! Save the group name
      call assert(803878279, func_groups%get_key(group_names(i_group)%string))

      ! Get the functional group
      call assert_msg(657972204, func_groups%get_property_t(val = func_group), &
              "Invalid functional group '"//group_names(i_group)%string// &
              "' in UNIFAC model.")

      ! Set the group volume parameter (R_k) Eq. 6
      key_name = "volume param"
      call assert_msg(549012632, &
              func_group%get_real(key_name, _R_k_(i_group)), &
              "Missing volume parameter in functional group '"// &
              group_names(i_group)%string//"' in UNIFAC model.")

      ! Set the group volume parameter (Q_k) Eq. 6
      key_name = "surface param"
      call assert_msg(127348854, &
              func_group%get_real(key_name, _Q_k_(i_group)), &
              "Missing surface parameter in functional group '"// &
              group_names(i_group)%string//"' in UNIFAC model.")

      ! Get the main group name
      key_name = "main group"
      call assert_msg(702688391, &
              func_group%get_string(key_name, main_group_name), &
              "Missing main group name in functional group '"// &
              group_names(i_group)%string//"' in UNIFAC model.")
            
      ! Set the main group id
      do i_main_group = 1, num_main_group
        found = .false.
        if (main_group_names(i_main_group)%string.eq.main_group_name) then
          main_group_id(i_group) = i_main_group
          found = .true.
          exit
        end if
      end do
      call assert_msg(752356165, found, "Missing main group '"// &
              main_group_name//"' needed by functional group '"// &
              group_names(i_group)%string//"' in UNIFAC model.")

      call func_groups%iter_next()
    end do

    ! Set the group interaction parameters
    do m = 1, _NUM_GROUP_
      do n = 1, _NUM_GROUP_
        _a_mn_(m,n) = main_group_interactions(main_group_id(m), &
                                              main_group_id(n))
      end do
    end do

    ! Set up parameters for each phase
    do i_UNIFAC_phase = 1, _NUM_UNIQUE_PHASE_
      phase_name = phase_names(i_UNIFAC_phase)%string
      i_phase = unique_phase_set_id(i_UNIFAC_phase)

      ! Set the properties for each species in the phase
      curr_spec_id = 0
      phase_ids_set = .false.
      do i_spec = 1, aero_phase_set(i_phase)%val%size()
        spec_name = aero_phase_set(i_phase)%val%get_species_name(i_spec)
        
        ! Get the species properties
        call assert(698678581, &
                chem_spec_data%get_property_set(spec_name, spec_props))
        
        ! Check if this is a UNIFAC species, and get its groups
        key_name = "UNIFAC groups"
        if (spec_props%get_property_t(key_name, spec_groups)) then
          curr_spec_id = curr_spec_id + 1

          ! Check the number of UNIFAC groups
          call assert_msg(511238330, spec_groups%size().gt.0, &
                  "Received empty set of UNIFAC groups for species '"// &
                  spec_name//"'")

          ! Set the functional groups (v_ik) Eq. 6
          ! and the r_i, q_i, and l_i parameters
          call spec_groups%iter_reset()
          do i_spec_group = 1, spec_groups%size()

            ! Get the group name
            call assert(649713038, spec_groups%get_key(spec_group_name))
            
            ! Get the number of this group for this species
            call assert_msg(429888360, spec_groups%get_int(val = num_spec_group), &
                    "Received non-integer number of UNIFAC groups for '"// &
                    spec_name//"'")

            ! Locate the group in the set of functional groups
            ! and set the v_ik parameter
            found = .false.
            do i_group = 1, size(group_names)
              if (group_names(i_group)%string.eq.spec_group_name) then
                found = .true.
                _v_ik_(i_UNIFAC_phase, curr_spec_id,i_group) = num_spec_group
                exit
              end if
            end do
            call assert_msg(175022713, found, &
                    "Invalid UNIFAC functional group specified for '"// &
                    spec_name//"'")

            ! Set the surface area (q_i) and volume (r_i) parameter for this
            ! species
            r_i = real(0.0, kind=dp)
            q_i = real(0.0, kind=dp)
            do i_group = 1, _NUM_GROUP_
              r_i = r_i + _R_k_(i_group) &
                      * real(_v_ik_(i_UNIFAC_phase ,curr_spec_id, i_group), &
                             kind=dp)
              q_i = q_i + _Q_k_(i_group) &
                      * real(_v_ik_(i_UNIFAC_phase, curr_spec_id, i_group), &
                             kind=dp)
            end do
            _r_i_(i_UNIFAC_phase, curr_spec_id) = r_i
            _q_i_(i_UNIFAC_phase, curr_spec_id) = q_i

            ! Set the l_i parameter for this species (Eq. 5)
            _l_i_(i_UNIFAC_phase, curr_spec_id) = 5.0d0 &
                    * ( r_i - q_i )  - ( r_i - 1.0d0 )

            call spec_groups%iter_next()
          end do

          ! Get the state ids of this species relative to the phase id
          ! If the phase ids have not been set, use this species to set them
          if (.not.phase_ids_set) then
            curr_phase_inst_id = 0
            do i_rep = 1, size(aero_rep_set)
              unique_names = aero_rep_set(i_rep)%val%unique_names( &
                      phase_name = phase_name, spec_name = spec_name )
              do i_instance = 1, size(unique_names)
                curr_phase_inst_id = curr_phase_inst_id + 1
                _PHASE_INST_ID_(i_UNIFAC_phase, curr_phase_inst_id) = &
                        aero_rep_set(i_rep)%val%spec_state_id( &
                        unique_names(i_instance)%string)
              end do
            end do
            _SPEC_ID_(i_UNIFAC_phase, curr_spec_id) = 0
            phase_ids_set = .true.
          else
            do i_rep = 1, size(aero_rep_set)
              unique_names = aero_rep_set(i_rep)%val%unique_names( &
                      phase_name = phase_name, spec_name = spec_name )
              if (size(unique_names).gt.0) then
                _SPEC_ID_(i_UNIFAC_phase, curr_spec_id) = &
                        aero_rep_set(i_rep)%val%spec_state_id( &
                        unique_names(1)%string) - &
                        _PHASE_INST_ID_(i_UNIFAC_phase, 1)
                exit
              end if
            end do 
          end if
        end if
      end do
    end do 
      
    call this%print()

    ! Clean up
    deallocate(num_phase_inst)
    deallocate(num_phase_spec)
    deallocate(unique_phase_set_id)
    deallocate(phase_names)
    deallocate(main_group_names)
    deallocate(main_group_interactions)
    deallocate(group_names)

  end subroutine initialize
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef _NUM_UNIQUE_PHASE_
#undef _NUM_GROUP_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _PHASE_INT_LOC_
#undef _PHASE_REAL_LOC_
#undef _NUM_PHASE_INSTANCE_
#undef _NUM_SPEC_
#undef _PHASE_INST_REAL_LOC_
#undef _PHASE_INST_ID_
#undef _SPEC_ID_
#undef _v_ik_

#undef _Q_k_
#undef _R_k_
#undef _a_mn_
#undef _PSI_mn_
#undef _r_i_
#undef _q_i_
#undef _l_i_
#undef _ln_GAMMA_ik_
#undef _gamma_i_
end module pmc_sub_model_UNIFAC
