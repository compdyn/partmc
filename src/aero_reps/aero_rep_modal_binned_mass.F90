! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_modal_binned_mass module.

!> \page phlex_aero_rep_modal_binned_mass Phlexible Module for Chemistry: Modal/Binned Mass Aerosol Representation
!!
!! The modal/binned mass aerosol representation includes a set of sections/bins
!! that are made up of one or more \ref phlex_aero_phase "aerosol phases." The
!! \c json object for this \ref phlex_aero_rep "aerosol representation" is of
!! the form:
!! \code{.json}
!!  { "pmc-data" : [
!!    {
!!      "type" : "AERO_REP_MODAL_BINNED_MASS",
!!      "sections/bins" : 
!!      {
!!        "dust" : 
!!        {
!!          "type" : "MODAL",
!!          "phases" : [ "insoluble", "organic", "aqueous" ],
!!          "shape" : "LOG_NORMAL",
!!          "geometric mean diameter" : 1.2e-6,
!!          "geometric standard deviation" : 1.2
!!        },
!!        "depeche" :
!!        {
!!          "type" : "BINNED",
!!          "phases" : [ "moody", "listless" ],
!!          "bins" : 8,
!!          "minimum diameter" : 0.8e-9,
!!          "maximum deviation" : 1.0e-6,
!!          "scale" : "LOG"
!!        }
!!      }
!!    },
!!    ...
!!  ]}
!! \endcode
!! The key-value pair \b type is required and must be
!! \b AERO_REP_MODAL_BINNED_MASS. The key-value pair \b "sections/bins" is also
!! required and must contain a set of at each one uniquely named mode or
!! bin-set key-value pairs whose values specify a \b type that must be either
!! 'MODAL' or 'BINNED' and an array of \b phases that correspond to existing
!! \ref phlex_aero_phase "aerosol phase" objects. Each phase will be present
!! once within a mode or once within each bin in a bin-set.
!!
!! Modes must also specify a distribution \b shape which must be 'LOG_NORMAL'
!! (this available shapes may be expanded in the future). Log-normal sections
!! must include a \b "geometric mean diameter" (m) and a 
!! \b "geometric standard deviation" (unitless) that will be used along with
!! the mass concetration of species in each phase and their densities to 
!! calculate a lognormal distribution for each mode at runtime.
!!
!! Bin sets must specify the number of \b bins, a \b "minimum diameter" and
!! \b "maximum diameter" and a \b scale, which must be 'LOG' or 'LINEAR'. The
!! number concentration will be calculated at run-time based on the total
!! mass (and thus volume) of each bin and the diameter of particles in that
!! bin.

!> The abstract aero_rep_modal_binned_mass_t structure and associated subroutines.
module pmc_aero_rep_modal_binned_mass

  use pmc_util,                               only: dp, i_kind, &
                                                    string_t, assert_msg, &
                                                    assert, die_msg, to_string
  use pmc_property
  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_aero_phase_data
  use pmc_phlex_state

  implicit none
  private

#define BINNED 1
#define MODAL 2

#define _NUM_SECTION_ this%condensed_data_int(1)
#define _INT_DATA_SIZE_ this%condensed_data_int(2)
#define _REAL_DATA_SIZE_ this%condensed_data_int(3)
#define _NUM_INT_PROP_ 3
#define _NUM_REAL_PROP_ 0
#define _MODE_INT_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+x)
#define _MODE_REAL_PARAM_LOC_(x) this%condensed_data_int(_NUM_INT_PROP_+_NUM_SECTION_+x)
#define _SECTION_TYPE_(x) this%condensed_data_int(_MODE_INT_PARAM_LOC_(x))

! For modes, _NUM_BINS_ = 1
#define _NUM_BINS_(x) this%condensed_data_int(_MODE_INT_PARAM_LOC_(x)+1)

#define _NUM_PHASE_(x) this%condensed_data_int(_MODE_INT_PARAM_LOC_(x)+2)
#define _PHASE_INT_PARAM_LOC_(x,y) this%condensed_data_int(_MODE_INT_PARAM_LOC_(x)+2+y)
#define _PHASE_REAL_PARAM_LOC_(x,y) this%condensed_data_int(_MODE_INT_PARAM_LOC_(x)+2+_NUM_PHASE_(x)+y)
#define _NUM_SPEC_(x,y) this%condensed_data_int(_PHASE_INT_PARAM_LOC_(x,y))

! Species state ids - for modes, b=1
#define _SPEC_STATE_ID_(x,y,b,z) this%condensed_data_int(_PHASE_INT_PARAM_LOC_(x,y)+(b-1)*_NUM_SPEC_(x,y)+z)

! GMD and bin diameter are stored in the same position - for modes, b=1
#define _GMD_(x,b) this%condensed_data_real(_MODE_REAL_PARAM_LOC_(x)+(b-1)*4)
#define _BIN_Dp_(x,b) this%condensed_data_real(_MODE_REAL_PARAM_LOC_(x)+(b-1)*4)

! GSD - only used for modes, b=1
#define _GSD_(x,b) this%condensed_data_real(_MODE_REAL_PARAM_LOC_(x)+(b-1)*4+1)

! Real-time number concetration - used for modes and bins - for modes, b=1
#define _NUMBER_CONC_(x,b) this%condensed_data_real(_MODE_REAL_PARAM_LOC_(x)+(b-1)*4+2)

! Real-time effective radius - only used for modesi, b=1
#define _EFFECTIVE_RADIUS_(x,b) this%condensed_data_real(_MODE_REAL_PARAM_LOC_(x)+(b-1)*4+3)

! Real-time aerosol phase mass - used for modes and bins - for modes, b=1
#define _AERO_PHASE_MASS_(x,y,b) this%condensed_data_real(_PHASE_REAL_PARAM_LOC_(x,y)-1+b)

! Species density
#define _DENSITY_(x,y,z) this%condensed_data_real(_PHASE_REAL_PARAM_LOC_(x,y)-1+_NUM_BINS_(x)+z)

  ! Update types (These must match values in aero_rep_modal_binned_mass.c)
  integer(kind=i_kind), parameter, public :: UPDATE_GMD = 0
  integer(kind=i_kind), parameter, public :: UPDATE_GSD = 1

  public :: aero_rep_modal_binned_mass_t

  !> Modal mass aerosol representation
  !!
  !! Time-invariant data related to a modal/binned mass aerosol representation. 
  type, extends(aero_rep_data_t) :: aero_rep_modal_binned_mass_t
    !> Mode names (only used during initialization)
    type(string_t), allocatable :: section_name(:)
    !> Phase state id (only used during initialization)
    integer(kind=i_kind), allocatable :: phase_state_id(:)
  contains
    !> Initialize the aerosol representation data, validating component data and
    !! loading any required information from the \c
    !! aero_rep_data_t::property_set. This routine should be called once for
    !! each aerosol representation at the beginning of a model run after all
    !! the input files have been read in. It ensures all data required during
    !! the model run are included in the condensed data arrays.
    procedure :: initialize
    !> Get the size of the section of the
    !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
    !! aerosol representation.
    !!
    !! For a modal/binned mass representation, the size will correspond to the
    !! the sum of the sizes of a single instance of each aerosol phase
    !! provided to \c aero_rep_modal_binned_mass::initialize()
    procedure :: size => get_size
    !> Get a list of unique names for each element on the
    !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
    !! representation. The list may be restricted to a particular phase and/or
    !! aerosol species by including the phase_name and spec_name arguments.
    !! 
    !! For a modal/binned mass representation, the unique names will be the
    !! phase name with the species name separated by a '.'
    procedure :: unique_names
    !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
    !! array by unique name. These are unique ids for each element on the
    !! state array for this \ref phlex_aero_rep  "aerosol representation" and
    !! are numbered:
    !!
    !!   \f$x_u = x_f ... (x_f+n-1)\f$
    !!
    !! where \f$x_u\f$ is the id of the element corresponding to concentration
    !! of the species with unique name \f$u\f$ on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
    !! of the first element for this aerosol representation on the state array
    !! and \f$n\f$ is the total number of variables on the state array from
    !! this aerosol representation.
    procedure :: spec_state_id
    !> Get the non-unique name of a species in this aerosol representation by
    !! id.
    procedure :: spec_name_by_id
    !> Get the number of instances of an aerosol phase in this representation
    procedure :: num_phase_instances

  end type aero_rep_modal_binned_mass_t

  ! Constructor for aero_rep_modal_binned_mass_t
  interface aero_rep_modal_binned_mass_t
    procedure :: constructor
  end interface aero_rep_modal_binned_mass_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_rep_modal_binned_mass_t
  function constructor() result (new_obj)

    !> New aerosol representation
    type(aero_rep_modal_binned_mass_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the \c
  !! aero_rep_data_t::property_set. This routine should be called once for
  !! each aerosol representation at the beginning of a model run after all
  !! the input files have been read in. It ensures all data required during
  !! the model run are included in the condensed data arrays.
  subroutine initialize(this, aero_phase_set, &
                  spec_state_id, chem_spec_data)

    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(inout) :: this
    !> The set of aerosol phases
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Beginning state id for this aerosol representationin the model species
    !! state array
    integer(kind=i_kind), intent(in) :: spec_state_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    type(property_t), pointer :: sections, section, phases, spec_props
    integer(kind=i_kind) :: i_section, i_phase, j_phase, k_phase, i_spec, &
                            i_bin
    integer(kind=i_kind) :: curr_spec_state_id
    integer(kind=i_kind) :: num_phase, num_bin
    integer(kind=i_kind) :: n_int_param, n_float_param
    integer(kind=i_kind) :: spec_type
    character(len=:), allocatable :: key_name, phase_name, sect_type, str_val
    real(kind=dp) :: min_Dp, max_Dp, d_log_Dp

    ! Determine the size of the condensed data arrays
    n_int_param = _NUM_INT_PROP_
    n_float_param = _NUM_REAL_PROP_

    ! Get the set of sections/bin-sets
    key_name = "modes/bins"
    call assert_msg(877855909, &
            this%property_set%get_property_t(key_name, sections), &
            "Missing sections/bins for modal/binned mass aerosol "// &
            "representation '"//this%rep_name//"'")
    call assert_msg(894962494, sections%size().gt.0, "No sections or bins "// &
            "specified for modal/binned mass aerosol representation '"// &
            this%rep_name//"'")

    ! Allocate space for the mode/bin names
    allocate(this%section_name(sections%size()))

    ! Loop through the sections, adding names and counting the spaces needed
    ! on the condensed data arrays, and counting the total phases instances
    num_phase = 0
    call sections%iter_reset()
    do i_section = 1, sections%size()
    
      ! Get the mode/bin name
      call assert(867378489, sections%get_key(key_name))
      call assert_msg(234513113, len(key_name).gt.0, "Missing mode/bin "// &
              "name in modal/binned mass aerosol representation '"// &
              this%rep_name//"'")
      this%section_name(i_section)%string = key_name

      ! Get the mode/bin properties
      call assert_msg(517138327, sections%get_property_t(val=section), &
              "Invalid structure for mode/bin '"// &
              this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"// &
              this%rep_name//"'")

      ! Get the section type
      key_name = "type"
      call assert_msg(742404898, section%get_string(key_name, sect_type), &
              "Missing mode/bin type in mode/bin '"// &
              this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"// &
              this%rep_name//"'")
      call assert_msg(618556995, &
              sect_type.eq."MODAL".or.sect_type.eq."BINNED", &
              "Invalid mode/bin type '"//sect_type//"' in mode/bin '"// &
              this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"// &
              this%rep_name//"'")
      
      ! Get the number of bins (or set to 1 for a mode)      
      num_bin = 1
      if (sect_type.eq."BINNED") then
      
        key_name = "bins"
        call assert_msg(824494286, section%get_int(key_name, num_bin), &
                "Missing number of bins in bin '"// &
              this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"// &
              this%rep_name//"'")
      end if 

      ! Add space for the mode/bin type, mode data locations, number of bins
      ! and phase count
      n_int_param = n_int_param + 5

      ! Add space for the GMD, GSD, number concentration, and effective radius
      n_float_param = n_float_param + 4*num_bin
 
      ! Get the set of phases
      key_name = "phases"
      call assert_msg(815518058, section%get_property_t(key_name, phases), &
              "Missing phases for mode '"//this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"//this%rep_name//"'")

      ! Add the phases to the counter
      call assert_msg(772593427, phases%size().gt.0, &
              "No phases specified for mode '"// &
              this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"//this%rep_name//"'")
      num_phase = num_phase + phases%size() * num_bin

      ! Loop through the phases and count the species
      call phases%iter_reset()
      do i_phase = 1, phases%size()

        ! Add space for phase data locations and number of species
        n_int_param = n_int_param + 3

        ! Get the phase name
        call assert_msg(393427582, phases%get_string(val=phase_name), &
                "Non-string phase name for mode '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")
        
        ! Find the aerosol phase and add space for the species variables
        do j_phase = 1, size(aero_phase_set)
          if (phase_name.eq.aero_phase_set(j_phase)%val%name()) then
            
            ! Add space for each species state id
            n_int_param = n_int_param + &
                    aero_phase_set(j_phase)%val%size() * num_bin

            ! Add space for total aerosol phase mass and each species density
            n_float_param = n_float_param + num_bin + &
                    aero_phase_set(j_phase)%val%size()

            exit
          else if (j_phase.eq.size(aero_phase_set)) then
            call die_msg(652391420, "Non-existant aerosol phase '"// &
                    phase_name//"' specified for mode '"// &
                    this%section_name(i_section)%string// &
                    "' in modal/binned mass aerosol representation '"// &
                    this%rep_name//"'")
          end if
        end do
        
        call phases%iter_next()
      end do

      call sections%iter_next()
    end do

    ! Allocate space for the aerosol phases and species state ids
    allocate(this%aero_phase(num_phase))
    allocate(this%phase_state_id(size(this%aero_phase)))

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(n_int_param))
    allocate(this%condensed_data_real(n_float_param))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)
    _INT_DATA_SIZE_ = n_int_param
    _REAL_DATA_SIZE_ = n_float_param

    ! Set the number of sections
    _NUM_SECTION_ = sections%size()

    ! Loop through the sections, adding names and distribution parameters and
    ! counting the phases in each section
    i_phase = 1
    curr_spec_state_id = spec_state_id
    n_int_param = _NUM_INT_PROP_+2*_NUM_SECTION_+1
    n_float_param = _NUM_REAL_PROP_+1
    call sections%iter_reset()
    do i_section = 1, _NUM_SECTION_
    
      ! Set the data locations for this mode
      _MODE_INT_PARAM_LOC_(i_section) = n_int_param
      _MODE_REAL_PARAM_LOC_(i_section) = n_float_param

      ! Get the mode/bin properties
      call assert(394743663, sections%get_property_t(val=section))

      ! Get the mode/bin type
      key_name = "type"
      call assert(667058653, section%get_string(key_name, sect_type))
      if (sect_type.eq."MODAL") then
        _SECTION_TYPE_(i_section) = MODAL
      else if (sect_type.eq."BINNED") then
        _SECTION_TYPE_(i_section) = BINNED
      else
        call die_msg(256924433, "Internal error")
      end if

      ! Get the number of bins (or set to 1 for a mode)      
      _NUM_BINS_(i_section) = 1
      if (_SECTION_TYPE_(i_section).eq.BINNED) then
        key_name = "bins"
        call assert(315215287, section%get_int(key_name, _NUM_BINS_(i_section)))
      end if

      ! Get mode parameters
      if (_SECTION_TYPE_(i_section).eq.MODAL) then
        
        ! Get the geometric mean diameter
        key_name = "geometric mean diameter"
        call assert_msg(414771933, &
                section%get_real(key_name, &
                _GMD_(i_section, _NUM_BINS_(i_section))), &
                "Missing geometric mean diameter in mode '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")

        ! Get the geometric standard deviation
        key_name = "geometric standard deviation"
        call assert_msg(163265059, &
                section%get_real(key_name, &
                _GSD_(i_section, _NUM_BINS_(i_section))), &
                "Missing geometric standard deviation in mode '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")

      ! Get bin parameters
      else if (_SECTION_TYPE_(i_section).eq.BINNED) then

        ! Get the minimum diameter (m)
        key_name = "minimum diameter"
        call assert_msg(548762180, section%get_real(key_name, min_Dp), &
                "Missing minimum diameter for bin '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")

        ! Get the maximum diameter (m)
        key_name = "maximum diameter"
        call assert_msg(288632226, section%get_real(key_name, max_Dp), &
                "Missing maximum diameter for bin '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")

        ! Get the scale
        key_name = "scale"
        call assert_msg(404761639, section%get_string(key_name, str_val), &
                "Missing bin scale for bin '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")

        ! Assign the bin diameters
        if (str_val.eq."LINEAR") then
          do i_bin = 1, _NUM_BINS_(i_section)
            _BIN_Dp_(i_section,i_bin) = min_Dp + &
                    (i_bin-1) * (max_Dp-min_Dp)/(_NUM_BINS_(i_section)-1)
          end do
        else if (str_val.eq."LOG") then
          d_log_Dp = (log10(max_Dp)-log10(min_Dp))/(_NUM_BINS_(i_section)-1)
          do i_bin = 1, _NUM_BINS_(i_section)
            _BIN_Dp_(i_section,i_bin) = 10.0d0**( log10(min_Dp) + &
                    (i_bin-1) * d_log_Dp )
          end do
        else
          call die_msg(236797392, "Invalid scale specified for bin '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")         
        end if

        ! Set the effective radius
        _EFFECTIVE_RADIUS_(i_section,i_bin) = _BIN_Dp_(i_section,i_bin)

      end if

      ! Get the set of phases
      key_name = "phases"
      call assert(712411046, section%get_property_t(key_name, phases))

      ! Save the number of phases
      _NUM_PHASE_(i_section) = phases%size()

      ! Add space for the mode/bin type, number of bins, number of phases and
      ! the phase data locations
      n_int_param = n_int_param + 3 + 2*_NUM_PHASE_(i_section)

      ! Add space for GMD, GSD, number concentration and effective radius
      n_float_param = n_float_param + 4 * _NUM_BINS_(i_section)

      ! Loop through the phase names, look them up, and add them to the list
      call phases%iter_reset()
      do j_phase = 1, phases%size()

        ! Set the data locations for this mode
        _PHASE_INT_PARAM_LOC_(i_section, j_phase) = n_int_param
        _PHASE_REAL_PARAM_LOC_(i_section, j_phase) = n_float_param

        ! Add space for the number of species
        n_int_param = n_int_param + 1

        ! Get the phase name
        call assert(775801035, phases%get_string(val=phase_name))
        
        ! Find the aerosol phase and add it to the list
        do k_phase = 1, size(aero_phase_set)
          if (phase_name.eq.aero_phase_set(k_phase)%val%name()) then
            
            ! Save the number of species
            _NUM_SPEC_(i_section, j_phase) = aero_phase_set(k_phase)%val%size()

            ! Loop through the bins
            do i_bin = 1, _NUM_BINS_(i_section)

              ! Add the aerosol phase to the list
              this%aero_phase(i_phase) = aero_phase_set(k_phase)

              ! Save the starting id for this phase on the state array
              this%phase_state_id(i_phase) = curr_spec_state_id

              ! Loop through the species in this phase
              do i_spec = 1, _NUM_SPEC_(i_section, j_phase)
            
                ! Save the species state ids for this phase
                _SPEC_STATE_ID_(i_section, j_phase, i_bin, i_spec) = curr_spec_state_id
                curr_spec_state_id = curr_spec_state_id + 1
            
                ! Get the property set for this species
                call assert_msg(741499484, chem_spec_data%get_property_set( &
                        this%aero_phase(i_phase)%val%get_species_name( i_spec), spec_props), &
                        "Missing property set for aerosol species '"// &
                        this%aero_phase(i_phase)%val%get_species_name(i_spec)// &
                        "' in phase '"//this%aero_phase(i_phase)%val%name()// &
                        "' of mode/bin '"// &
                        this%section_name(i_section)%string// &
                        "' in modal/binned mass aerosol representation '"//this%rep_name//"'")

                ! Get the species density (expect for ion-pair activity species)
                call assert(874736323, chem_spec_data%get_type( &
                        this%aero_phase(i_phase)%val%get_species_name( i_spec), spec_type))
                if (spec_type.eq.CHEM_SPEC_ACTIVITY_COEFF) then
                    _DENSITY_(i_section, j_phase, i_spec) = 0.0
                else
                  key_name = "density"
                  call assert_msg(821683338, spec_props%get_real(key_name, &
                          _DENSITY_(i_section, j_phase, i_spec)), &
                          "Missing density for aerosol species '"// &
                          this%aero_phase(i_phase)%val%get_species_name(i_spec)// &
                          "' in phase '"//this%aero_phase(i_phase)%val%name()// &
                          "' of mode/bin '"// &
                          this%section_name(i_section)%string// &
                          "' in modal/binned mass aerosol representation '"//this%rep_name//"'")
                end if

              end do

              ! Add space for species state ids
              n_int_param = n_int_param + this%aero_phase(i_phase)%val%size()

              i_phase = i_phase + 1
            end do

            ! Add space for aerosol phase mass and species densities
            n_float_param = n_float_param + _NUM_BINS_(i_section) + &
                    aero_phase_set(k_phase)%val%size()

            exit
          else if (k_phase.eq.size(aero_phase_set)) then
            call die_msg(652391420, "Internal error.")
          end if
        end do

        call phases%iter_next()
      end do

      call sections%iter_next()
    end do

    ! Check the data sizes
    call assert(951534966, i_phase-1.eq.num_phase)
    call assert(951534966, n_int_param.eq._INT_DATA_SIZE_+1)
    call assert(325387136, n_float_param.eq._REAL_DATA_SIZE_+1)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of the section of the
  !! \c pmc_phlex_state::phlex_state_t::state_var array required for this
  !! aerosol representation.
  !!
  !! For a modal/binned mass representation, the size will correspond to the
  !! the sum of the sizes of a single instance of each aerosol phase
  !! provided to \c aero_rep_modal_binned_mass::initialize()
  function get_size(this) result (state_size)

    !> Size on the state array
    integer(kind=i_kind) :: state_size
    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(in) :: this

    integer(kind=i_kind) :: i_phase

    ! Get the total number of species across all phases
    state_size = 0
    do i_phase = 1, size(this%aero_phase)
      state_size = state_size + this%aero_phase(i_phase)%val%size()
    end do

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of unique names for each element on the
  !! \c pmc_phlex_state::phlex_state_t::state_var array for this aerosol
  !! representation. The list may be restricted to a particular phase and/or
  !! aerosol species by including the phase_name and spec_name arguments.
  !! 
  !! For a modal/binned mass representation, the unique names will be the
  !! phase name with the species name separated by a '.'
  function unique_names(this, phase_name, tracer_type, spec_name)

    !> List of unique names
    type(string_t), allocatable :: unique_names(:)
    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(in) :: this
    !> Aerosol phase name
    character(len=:), allocatable, optional, intent(in) :: phase_name
    !> Aerosol-phase species tracer type
    integer(kind=i_kind), optional, intent(in) :: tracer_type
    !> Aerosol-phase species name
    character(len=:), allocatable, optional, intent(in) :: spec_name

    integer(kind=i_kind) :: num_spec, i_spec, j_spec, i_phase, j_phase, &
                            i_section, i_bin
    integer(kind=i_kind) :: curr_tracer_type
    character(len=:), allocatable :: curr_section_name, curr_phase_name, &
                                     curr_spec_name, curr_bin_str
    
    ! Count the number of unique names
    num_spec = 0
    do i_phase = 1, size(this%aero_phase)

      ! Filter by phase name
      if (present(phase_name)) then
        curr_phase_name = this%aero_phase(i_phase)%val%name()
        if (phase_name.ne.curr_phase_name) cycle
      end if

      ! Filter by spec name and/or tracer type
      if (present(spec_name).or.present(tracer_type)) then
        do j_spec = 1, this%aero_phase(i_phase)%val%size()
          curr_spec_name = &
                  this%aero_phase(i_phase)%val%get_species_name(j_spec)
          curr_tracer_type = &
                  this%aero_phase(i_phase)%val%get_species_type(j_spec)
          if (present(spec_name)) then
            if (spec_name.ne.curr_spec_name) cycle
          end if
          if (present(tracer_type)) then
            if (tracer_type.ne.curr_tracer_type) cycle
          end if
          num_spec = num_spec + 1
        end do
      else
        num_spec = num_spec + this%aero_phase(i_phase)%val%size()
      end if

    end do

    ! Allocate space for the unique names
    allocate(unique_names(num_spec))

    ! Loop through the modes/bin sets
    i_phase = 1
    i_spec = 1
    do i_section = 1, _NUM_SECTION_

      ! Get the current section name
      curr_section_name = this%section_name(i_section)%string

      ! Loop through the phases for this mode/bin set
      do j_phase = 1, _NUM_PHASE_(i_section)

        ! Set the current phase name
        curr_phase_name = this%aero_phase(i_phase)%val%name()
      
        ! Filter by phase name
        if (present(phase_name)) then
          if (phase_name.ne.curr_phase_name) then
            i_phase = i_phase + _NUM_BINS_(i_section)
            cycle
          end if
        end if

        ! Loop through the bins (one iteration for modes)
        do i_bin = 1, _NUM_BINS_(i_section)

          ! Set the current bin label (except for single bins or mode)
          if (_NUM_BINS_(i_section).gt.1) then
            curr_bin_str = trim(to_string(i_bin))//"."
          else
            curr_bin_str = ""
          end if

          ! Add species from this phase/bin
          num_spec = this%aero_phase(i_phase)%val%size()
          do j_spec = 1, num_spec

            ! Get the species name
            curr_spec_name = &
                  this%aero_phase(i_phase)%val%get_species_name(j_spec)
        
            ! Filter by species name
            if (present(spec_name)) then
              if (spec_name.ne.curr_spec_name) cycle
            end if

            ! Filter by species tracer type
            if (present(tracer_type)) then
              curr_tracer_type = &
                  this%aero_phase(i_phase)%val%get_species_type(j_spec)
              if (tracer_type.ne.curr_tracer_type) cycle
            end if

            ! Add the unique name for this species
            unique_names(i_spec)%string = curr_section_name//"."// &
                  curr_bin_str//curr_phase_name//'.'//curr_spec_name
          
            i_spec = i_spec + 1
          end do
          
          ! Move to the next phase instance
          i_phase = i_phase + 1
        
        end do
      end do
    end do

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
  !! array by unique name. These are unique ids for each element on the
  !! state array for this \ref phlex_aero_rep  "aerosol representation" and
  !! are numbered:
  !!
  !!   \f$x_u = x_f ... (x_f+n-1)\f$
  !!
  !! where \f$x_u\f$ is the id of the element corresponding to concentration
  !! of the species with unique name \f$u\f$ on the \c
  !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
  !! of the first element for this aerosol representation on the state array
  !! and \f$n\f$ is the total number of variables on the state array from
  !! this aerosol representation.
  function spec_state_id(this, unique_name) result (spec_id)

    !> Species state id
    integer(kind=i_kind) :: spec_id
    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(in) :: this
    !> Unique name
    character(len=:), allocatable, intent(in) :: unique_name

    type(string_t), allocatable, save :: unique_names(:)
    integer(kind=i_kind) :: i_spec

    spec_id = 0
    if (.not.allocated(unique_names)) unique_names = this%unique_names()
    do i_spec = 1, size(unique_names)
      if (unique_names(i_spec)%string .eq. unique_name) then
        spec_id = this%phase_state_id(1) + i_spec - 1
        return
      end if
    end do

  end function spec_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the non-unique name of a species in this aerosol representation by
  !! id.
  function spec_name_by_id(this, aero_rep_spec_id)

    !> Chemical species name
    character(len=:), allocatable :: spec_name_by_id
    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(in) :: this
    !> Indoex of species in this aerosol representation
    integer(kind=i_kind) :: aero_rep_spec_id

    ! Indices for iterators
    integer(kind=i_kind) :: i_spec, j_spec, i_phase
    
    i_spec = 1
    do i_phase = 1, size(this%aero_phase)
      do j_spec = 1, this%aero_phase(i_phase)%val%size()
        if (i_spec.eq.aero_rep_spec_id) then
          spec_name_by_id = this%aero_phase(i_phase)%val%get_species_name(j_spec)
        end if
        i_spec = i_spec + 1
      end do
    end do

  end function spec_name_by_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the number of instances of a specified aerosol phase.
  function num_phase_instances(this, phase_name)

    !> Number of instances of the aerosol phase
    integer(kind=i_kind) :: num_phase_instances
    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(in) :: this
    !> Aerosol phase name
    character(len=:), allocatable, intent(in) :: phase_name

    integer(kind=i_kind) :: i_phase

    num_phase_instances = 0
    do i_phase = 1, size(this%aero_phase)
      if (this%aero_phase(i_phase)%val%name().eq.phase_name) then
        num_phase_instances = 1
        return
      end if
    end do

  end function num_phase_instances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef BINNED
#undef MODAL
#undef _NUM_SECTION_
#undef _INT_DATA_SIZE_
#undef _REAL_DATA_SIZE_
#undef _NUM_INT_PROP_
#undef _NUM_REAL_PROP_
#undef _MODE_INT_PARAM_LOC_
#undef _MODE_REAL_PARAM_LOC_
#undef _SECTION_TYPE_
#undef _NUM_BINS_
#undef _NUM_PHASE_
#undef _PHASE_INT_PARAM_LOC_
#undef _PHASE_REAL_PARAM_LOC_
#undef _NUM_SPEC_
#undef _SPEC_STATE_ID_
#undef _GMD_
#undef _BIN_Dp_
#undef _GSD_
#undef _NUMBER_CONC_
#undef _EFFECTIVE_RADUIS_
#undef _AERO_PHASE_MASS_
#undef _DENSITY_

end module pmc_aero_rep_modal_binned_mass
