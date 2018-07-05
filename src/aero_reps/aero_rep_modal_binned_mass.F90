! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_modal_binned_mass module.

!> \page phlex_aero_rep_modal_binned_mass Phlexible Module for Chemistry: Modal/Binned Mass Aerosol Representation
!!
!! The modal/binned mass aerosol representation includes a set of sections/bins
!! that are made up of one or more \ref phlex_aero_phase "aerosol phases." The
!! \c json object for this \ref phlex_aero_rep "aerosol representation" has the
!! following format:
!! \code{.json}
!!  { "pmc-data" : [
!!    {
!!      "name" : "my modal/binned aero rep",
!!      "type" : "AERO_REP_MODAL_BINNED_MASS",
!!      "modes/bins" : 
!!      {
!!        "dust" : 
!!        {
!!          "type" : "BINNED",
!!          "phases" : [ "insoluble", "organic", "aqueous" ],
!!          "bins" : 8,
!!          "minimum diameter" : 0.8e-9,
!!          "maximum deviation" : 1.0e-6,
!!          "scale" : "LOG"
!!        },
!!        "depeche" :
!!        {
!!          "type" : "MODAL",
!!          "phases" : [ "moody", "listless" ],
!!          "shape" : "LOG_NORMAL",
!!          "geometric mean diameter" : 1.2e-6,
!!          "geometric standard deviation" : 1.2
!!        }
!!      }
!!    },
!!    ...
!!  ]}
!! \endcode
!! The key-value pair \b type is required and must be
!! \b AERO_REP_MODAL_BINNED_MASS. The key-value pair \b modes/bins is also
!! required and must contain a set of at least one uniquely named mode or
!! bin-set key-value pair whose value(s) specify a \b type that must be either
!! \b MODAL or \b BINNED and an array of \b phases that correspond to existing
!! \ref phlex_aero_phase "aerosol phase" objects. Each phase will be present
!! once within a mode or once within each bin in a bin-set.
!!
!! Modes must also specify a distribution \b shape which must be \b LOG_NORMAL
!! (the available shapes may be expanded in the future). Log-normal sections
!! must include a \b geometric \b mean \b diameter (m) and a \b geometric 
!! \b standard \b deviation (unitless) that will be used along with the mass
!! concentration of species in each phase and their densities to calculate a
!! lognormal distribution for each mode at runtime.
!!
!! Bin sets must specify the number of \b bins, a \b minimum \b diameter (m),
!! a \b maximum \b diameter (m) and a \b scale, which must be \b LOG or 
!! \b LINEAR. The number concentration will be calculated at run-time based on
!! the total mass of each bin, the species densities and the diameter of
!! particles in that bin.

!> The abstract aero_rep_modal_binned_mass_t structure and associated subroutines.
module pmc_aero_rep_modal_binned_mass

  use pmc_aero_phase_data
  use pmc_aero_rep_data
  use pmc_chem_spec_data
  use pmc_phlex_state
  use pmc_property
  use pmc_util,                               only: dp, i_kind, &
                                                    string_t, assert_msg, &
                                                    assert, die_msg, &
                                                    to_string, align_ratio

  use iso_c_binding

  implicit none
  private

#define BINNED 1
#define MODAL 2

#define NUM_SECTION_ this%condensed_data_int(1)
#define INT_DATA_SIZE_ this%condensed_data_int(2)
#define FLOAT_DATA_SIZE_ this%condensed_data_int(3)
#define AERO_REP_ID_ this%condensed_data_int(4)
#define NUM_INT_PROP_ 4
#define NUM_REAL_PROP_ 0
#define MODE_INT_PROP_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+x)
#define MODE_REAL_PROP_LOC_(x) this%condensed_data_int(NUM_INT_PROP_+NUM_SECTION_+x)
#define SECTION_TYPE_(x) this%condensed_data_int(MODE_INT_PROP_LOC_(x))

! For modes, NUM_BINS_ = 1
#define NUM_BINS_(x) this%condensed_data_int(MODE_INT_PROP_LOC_(x)+1)

! Number of aerosol phases in this mode/bin set
#define NUM_PHASE_(x) this%condensed_data_int(MODE_INT_PROP_LOC_(x)+2)

! Phase state and model data ids
#define PHASE_STATE_ID_(x,y,b) this%condensed_data_int(MODE_INT_PROP_LOC_(x)+2+(b-1)*NUM_PHASE_(x)+y)
#define PHASE_MODEL_DATA_ID_(x,y,b) this%condensed_data_int(MODE_INT_PROP_LOC_(x)+2+NUM_BINS_(x)*NUM_PHASE_(x)+(b-1)*NUM_PHASE_(x)+y)

! GMD and bin diameter are stored in the same position - for modes, b=1
#define GMD_(x,b) this%condensed_data_real(MODE_REAL_PROP_LOC_(x)+(b-1)*4)
#define BIN_DP_(x,b) this%condensed_data_real(MODE_REAL_PROP_LOC_(x)+(b-1)*4)

! GSD - only used for modes, b=1
#define GSD_(x,b) this%condensed_data_real(MODE_REAL_PROP_LOC_(x)+(b-1)*4+1)

! Real-time number concetration - used for modes and bins - for modes, b=1
#define NUMBER_CONC_(x,b) this%condensed_data_real(MODE_REAL_PROP_LOC_(x)+(b-1)*4+2)

! Real-time effective radius - only used for modesi, b=1
#define EFFECTIVE_RADIUS_(x,b) this%condensed_data_real(MODE_REAL_PROP_LOC_(x)+(b-1)*4+3)

! Real-time aerosol phase mass - used for modes and bins - for modes, b=1
#define PHASE_MASS_(x,y,b) this%condensed_data_real(MODE_REAL_PROP_LOC_(x)+4*NUM_BINS_(x)+(b-1)*NUM_PHASE_(x)+y-1)

! Real-time aerosol phase average MW - used for modes and bins - for modes, b=0
#define PHASE_AVG_MW_(x,y,b) this%condensed_data_real(MODE_REAL_PROP_LOC_(x)+(4+NUM_PHASE_(x))*NUM_BINS_(x)+(b-1)*NUM_PHASE_(x)+y-1)

  ! Update types (These must match values in aero_rep_modal_binned_mass.c)
  integer(kind=i_kind), parameter, public :: UPDATE_GMD = 0
  integer(kind=i_kind), parameter, public :: UPDATE_GSD = 1

  public :: aero_rep_modal_binned_mass_t, &
            aero_rep_update_data_modal_binned_mass_GMD_t, &
            aero_rep_update_data_modal_binned_mass_GSD_t

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
    !> Set an id for this aerosol representation for use with updates from
    !! external modules
    procedure :: set_id
    !> Get an id for a mode or bin in the aerosol representation by name for
    !! use with updates from external modules
    procedure :: get_section_id
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
    !! For a modal/binned mass representation, the unique names for bins are:
    !!   - "bin name.bin #.phase name.species name"
    !!
    !! ... and for modes are:
    !!   - "mode name.phase name.species name"
    procedure :: unique_names
    !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
    !! array by its unique name. These are unique ids for each element on the
    !! state array for this \ref phlex_aero_rep "aerosol representation" and
    !! are numbered:
    !!
    !!   \f[x_u \in x_f ... (x_f+n-1)\f]
    !!
    !! where \f$x_u\f$ is the id of the element corresponding to the species
    !! with unique name \f$u\f$ on the \c
    !! pmc_phlex_state::phlex_state_t::state_var array, \f$x_f\f$ is the index
    !! of the first element for this aerosol representation on the state array
    !! and \f$n\f$ is the total number of variables on the state array from
    !! this aerosol representation.
    procedure :: spec_state_id
    !> Get the non-unique name of a species by its unique name
    procedure :: spec_name
    !> Get the number of instances of an aerosol phase in this representation
    procedure :: num_phase_instances
    !> Finalize the aerosol representation
    final :: finalize

  end type aero_rep_modal_binned_mass_t

  ! Constructor for aero_rep_modal_binned_mass_t
  interface aero_rep_modal_binned_mass_t
    procedure :: constructor
  end interface aero_rep_modal_binned_mass_t

  !> Update GMD object
  type, extends(aero_rep_update_data_t) :: &
            aero_rep_update_data_modal_binned_mass_GMD_t
  private
    logical :: is_malloced = .false.
  contains
    !> Update the GMD
    procedure :: set_GMD => update_data_set_GMD
    !> Finalize the GMD update data
    final :: update_data_GMD_finalize
  end type aero_rep_update_data_modal_binned_mass_GMD_t

  !> Constructor for aero_rep_update_data_modal_binned_mass_GMD_t
  interface aero_rep_update_data_modal_binned_mass_GMD_t
    procedure :: update_data_GMD_constructor
  end interface aero_rep_update_data_modal_binned_mass_GMD_t

  !> Update GSD object
  type, extends(aero_rep_update_data_t) :: &
            aero_rep_update_data_modal_binned_mass_GSD_t
  private
    logical :: is_malloced = .false.
  contains
    !> Update the GSD
    procedure :: set_GSD => update_data_set_GSD
    !> Finalize the GSD update data
    final :: update_data_GSD_finalize
  end type aero_rep_update_data_modal_binned_mass_GSD_t

  !> Constructor for aero_rep_update_data_modal_binned_mass_GSD_t
  interface aero_rep_update_data_modal_binned_mass_GSD_t
    procedure :: update_data_GSD_constructor
  end interface aero_rep_update_data_modal_binned_mass_GSD_t

  !> Interface to c aerosol representation functions
  interface

    !> Allocate space for a GMD update object
    function aero_rep_modal_binned_mass_create_gmd_update_data() &
              result (update_data) bind (c)
      use iso_c_binding
      !> Allocated update data object
      type(c_ptr) :: update_data
    end function aero_rep_modal_binned_mass_create_gmd_update_data

    !> Set a new mode GMD
    subroutine aero_rep_modal_binned_mass_set_gmd_update_data(update_data, &
              aero_rep_id, section_id, gmd) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Aerosol representation id from
      !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::set_id
      integer(kind=c_int), value :: aero_rep_id
      !> Section id from
      !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::get_section_id
      integer(kind=c_int), value :: section_id
      !> New GMD (m)
      real(kind=c_double), value :: gmd
    end subroutine aero_rep_modal_binned_mass_set_gmd_update_data

    !> Allocate space for a GSD update object
    function aero_rep_modal_binned_mass_create_gsd_update_data() &
              result (update_data) bind (c)
      use iso_c_binding
      !> Allocated update data object
      type(c_ptr) :: update_data
    end function aero_rep_modal_binned_mass_create_gsd_update_data

    !> Set a new mode GSD
    subroutine aero_rep_modal_binned_mass_set_gsd_update_data(update_data, &
              aero_rep_id, section_id, gsd) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value :: update_data
      !> Aerosol representation id from
      !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::set_id
      integer(kind=c_int), value :: aero_rep_id
      !> Section id from
      !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::get_section_id
      integer(kind=c_int), value :: section_id
      !> New GSD (m)
      real(kind=c_double), value :: gsd
    end subroutine aero_rep_modal_binned_mass_set_gsd_update_data

    !> Free an update data object
    pure subroutine aero_rep_free_update_data(update_data) bind (c)
      use iso_c_binding
      !> Update data
      type(c_ptr), value, intent(in) :: update_data
    end subroutine aero_rep_free_update_data

  end interface

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
  subroutine initialize(this, aero_phase_set, spec_state_id)

    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(inout) :: this
    !> The set of aerosol phases
    type(aero_phase_data_ptr), pointer, intent(in) :: aero_phase_set(:)
    !> Beginning state id for this aerosol representationin the model species
    !! state array
    integer(kind=i_kind), intent(in) :: spec_state_id

    type(property_t), pointer :: sections, section, phases
    integer(kind=i_kind) :: i_section, i_phase, j_phase, k_phase, &
                            i_bin
    integer(kind=i_kind) :: curr_spec_state_id
    integer(kind=i_kind) :: num_phase, num_bin
    integer(kind=i_kind) :: n_int_param, n_float_param
    character(len=:), allocatable :: key_name, phase_name, sect_type, str_val
    real(kind=dp) :: min_Dp, max_Dp, d_log_Dp
    integer(kind=i_kind) :: int_data_size, float_data_size

    ! Determine the size of the condensed data arrays
    n_int_param = NUM_INT_PROP_
    n_float_param = NUM_REAL_PROP_

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

      ! Add space for the mode/bin type, number of bins, and phase count
      ! and parameter locations
      n_int_param = n_int_param + 5

      ! Add space for the GMD, GSD, number concentration, and effective radius
      n_float_param = n_float_param + 4*num_bin
 
      ! Get the set of phases
      key_name = "phases"
      call assert_msg(815518058, section%get_property_t(key_name, phases), &
              "Missing phases for mode '"// &
              this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"// &
              this%rep_name//"'")

      ! Add the phases to the counter
      call assert_msg(772593427, phases%size().gt.0, &
              "No phases specified for mode '"// &
              this%section_name(i_section)%string// &
              "' in modal/binned mass aerosol representation '"// &
              this%rep_name//"'")
      num_phase = num_phase + phases%size() * num_bin

      ! Loop through the phases and make sure they exist
      call phases%iter_reset()
      do i_phase = 1, phases%size()

        ! Get the phase name
        call assert_msg(393427582, phases%get_string(val=phase_name), &
                "Non-string phase name for mode '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")
        
        ! Find the aerosol phase and add space for its variables
        do j_phase = 1, size(aero_phase_set)
          if (phase_name.eq.aero_phase_set(j_phase)%val%name()) then
            
            ! Add space for the phase state and model data ids
            n_int_param = n_int_param + 2 * num_bin

            ! Add space for total aerosol phase mass and average MW,
            n_float_param = n_float_param + 2 * num_bin

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

    ! Calculate int and float array sizes with alignment spacing
    int_data_size = n_int_param
    int_data_size = int_data_size + mod(int_data_size, align_ratio)
    float_data_size = n_float_param

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(int_data_size))
    allocate(this%condensed_data_real(float_data_size))
    this%condensed_data_int(:) = int(0, kind=i_kind)
    this%condensed_data_real(:) = real(0.0, kind=dp)
    INT_DATA_SIZE_ = int_data_size
    FLOAT_DATA_SIZE_ = float_data_size

    ! Set the number of sections
    NUM_SECTION_ = sections%size()

    ! Loop through the sections, adding names and distribution parameters and
    ! counting the phases in each section
    i_phase = 1
    curr_spec_state_id = spec_state_id
    n_int_param = NUM_INT_PROP_+2*NUM_SECTION_+1
    n_float_param = NUM_REAL_PROP_+1
    call sections%iter_reset()
    do i_section = 1, NUM_SECTION_
    
      ! Set the data locations for this mode
      MODE_INT_PROP_LOC_(i_section) = n_int_param
      MODE_REAL_PROP_LOC_(i_section) = n_float_param

      ! Get the mode/bin properties
      call assert(394743663, sections%get_property_t(val=section))

      ! Get the mode/bin type
      key_name = "type"
      call assert(667058653, section%get_string(key_name, sect_type))
      if (sect_type.eq."MODAL") then
        SECTION_TYPE_(i_section) = MODAL
      else if (sect_type.eq."BINNED") then
        SECTION_TYPE_(i_section) = BINNED
      else
        call die_msg(256924433, "Internal error")
      end if

      ! Get the number of bins (or set to 1 for a mode)      
      NUM_BINS_(i_section) = 1
      if (SECTION_TYPE_(i_section).eq.BINNED) then
        key_name = "bins"
        call assert(315215287, section%get_int(key_name, NUM_BINS_(i_section)))
      end if

      ! Get mode parameters
      if (SECTION_TYPE_(i_section).eq.MODAL) then
        
        ! Get the geometric mean diameter
        key_name = "geometric mean diameter"
        call assert_msg(414771933, &
                section%get_real(key_name, &
                GMD_(i_section, NUM_BINS_(i_section))), &
                "Missing geometric mean diameter in mode '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")

        ! Get the geometric standard deviation
        key_name = "geometric standard deviation"
        call assert_msg(163265059, &
                section%get_real(key_name, &
                GSD_(i_section, NUM_BINS_(i_section))), &
                "Missing geometric standard deviation in mode '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")

        ! Calculate the effective radius (m)
        ! (See aero_rep_modal_binned_mass_get_effective_radius for details)
        EFFECTIVE_RADIUS_(i_section, NUM_BINS_(i_section)) = &
                GMD_(i_section, NUM_BINS_(i_section)) / 2.0d0 * &
                exp(5.0d0/2.0d0*(GSD_(i_section, NUM_BINS_(i_section)))**2) 

      ! Get bin parameters
      else if (SECTION_TYPE_(i_section).eq.BINNED) then

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
          do i_bin = 1, NUM_BINS_(i_section)
            BIN_DP_(i_section,i_bin) = min_Dp + &
                    (i_bin-1) * (max_Dp-min_Dp)/(NUM_BINS_(i_section)-1)
          end do
        else if (str_val.eq."LOG") then
          d_log_Dp = (log10(max_Dp)-log10(min_Dp))/(NUM_BINS_(i_section)-1)
          do i_bin = 1, NUM_BINS_(i_section)
            BIN_DP_(i_section,i_bin) = 10.0d0**( log10(min_Dp) + &
                    (i_bin-1) * d_log_Dp )
          end do
        else
          call die_msg(236797392, "Invalid scale specified for bin '"// &
                this%section_name(i_section)%string// &
                "' in modal/binned mass aerosol representation '"// &
                this%rep_name//"'")         
        end if

        ! Set the effective radius
        EFFECTIVE_RADIUS_(i_section,i_bin) = BIN_DP_(i_section,i_bin)

      end if

      ! Get the set of phases
      key_name = "phases"
      call assert(712411046, section%get_property_t(key_name, phases))

      ! Save the number of phases
      NUM_PHASE_(i_section) = phases%size()

      ! Add space for the mode/bin type, number of bins, and number of phases
      n_int_param = n_int_param + 3

      ! Add space for GMD, GSD, number concentration and effective radius
      n_float_param = n_float_param + 4 * NUM_BINS_(i_section)

      ! Loop through the phase names, look them up, and add them to the list
      call phases%iter_reset()
      do j_phase = 1, phases%size()

        ! Get the phase name
        call assert(775801035, phases%get_string(val=phase_name))
        
        ! Find the aerosol phase and add it to the list
        do k_phase = 1, size(aero_phase_set)
          if (phase_name.eq.aero_phase_set(k_phase)%val%name()) then
            
            ! Loop through the bins
            do i_bin = 1, NUM_BINS_(i_section)

              ! Add the aerosol phase to the list
              this%aero_phase(i_phase) = aero_phase_set(k_phase)

              ! Save the starting id for this phase on the state array
              this%phase_state_id(i_phase) = curr_spec_state_id
              PHASE_STATE_ID_(i_section, j_phase, i_bin) = curr_spec_state_id

              ! Increment the state id by the size of the phase
              curr_spec_state_id = curr_spec_state_id + &
                        aero_phase_set(k_phase)%val%size()

              ! Save the phase model data id
              PHASE_MODEL_DATA_ID_(i_section, j_phase, i_bin) = k_phase

              i_phase = i_phase + 1
            end do

            ! Add space for aerosol phase state and model data ids
            n_int_param = n_int_param + 2*NUM_BINS_(i_section)

            ! Add space for aerosol phase mass and average MW
            n_float_param = n_float_param + 2*NUM_BINS_(i_section)

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
    call assert(831761020, i_phase-1.eq.num_phase)
    int_data_size = n_int_param - 1 + mod(n_int_param - 1, align_ratio)
    float_data_size = n_float_param - 1
    call assert(951534966, int_data_size.eq.INT_DATA_SIZE_)
    call assert(325387136, float_data_size.eq.FLOAT_DATA_SIZE_)

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set an id for this aerosol representation that can be used by external
  !! modules to update the GMD or GSD of a mode
  subroutine set_id(this, new_id)

    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(inout) :: this
    !> New id
    integer(kind=i_kind), intent(in) :: new_id

    AERO_REP_ID_ = new_id

  end subroutine set_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an id for a mode or bin by name for use with updates from external
  !! modules
  function get_section_id(this, section_name, section_id) result (found)

    !> Flag indicating whether the mode/bin was found
    logical :: found
    !> Aerosol representation
    class(aero_rep_modal_binned_mass_t), intent(in) :: this
    !> Section name
    character(len=:), allocatable, intent(in) :: section_name
    !> Section id
    integer(kind=i_kind), intent(out) :: section_id

    integer(kind=i_kind) :: i_section

    call assert_msg(194186171, allocated(this%section_name), &
            "Trying to get section id of uninitialized aerosol "// &
            "representation.")

    found = .false.
    do i_section = 1, size(this%section_name)
      if (this%section_name(i_section)%string.eq.section_name) then
        found = .true.
        section_id = i_section
        return
      end if
    end do

  end function get_section_id

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
  !! For a modal/binned mass representation, the unique names for bins are:
  !!   - "bin name.bin #.phase name.species name"
  !!
  !! ... and for modes are:
  !!   - "mode name.phase name.species name"
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
                                     curr_bin_str
    type(string_t), allocatable :: spec_names(:)
    
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
        spec_names = this%aero_phase(i_phase)%val%get_species_names()
        do j_spec = 1, size(spec_names)
          curr_tracer_type = &
                  this%aero_phase(i_phase)%val%get_species_type( &
                  spec_names(j_spec)%string)
          if (present(spec_name)) then
            if (spec_name.ne.spec_names(j_spec)%string) cycle
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
    do i_section = 1, NUM_SECTION_

      ! Get the current section name
      curr_section_name = this%section_name(i_section)%string

      ! Loop through the phases for this mode/bin set
      do j_phase = 1, NUM_PHASE_(i_section)

        ! Set the current phase name
        curr_phase_name = this%aero_phase(i_phase)%val%name()
      
        ! Filter by phase name
        if (present(phase_name)) then
          if (phase_name.ne.curr_phase_name) then
            i_phase = i_phase + NUM_BINS_(i_section)
            cycle
          end if
        end if

        ! Get the species names in this phase
        spec_names = this%aero_phase(i_phase)%val%get_species_names()

        ! Loop through the bins (one iteration for modes)
        do i_bin = 1, NUM_BINS_(i_section)

          ! Set the current bin label (except for single bins or mode)
          if (NUM_BINS_(i_section).gt.1) then
            curr_bin_str = trim(to_string(i_bin))//"."
          else
            curr_bin_str = ""
          end if

          ! Add species from this phase/bin
          num_spec = this%aero_phase(i_phase)%val%size()
          do j_spec = 1, num_spec

            ! Filter by species name
            if (present(spec_name)) then
              if (spec_name.ne.spec_names(j_spec)%string) cycle
            end if

            ! Filter by species tracer type
            if (present(tracer_type)) then
              curr_tracer_type = &
                  this%aero_phase(i_phase)%val%get_species_type( &
                  spec_names(j_spec)%string)
              if (tracer_type.ne.curr_tracer_type) cycle
            end if

            ! Add the unique name for this species
            unique_names(i_spec)%string = curr_section_name//"."// &
                  curr_bin_str//curr_phase_name//'.'// &
                  spec_names(j_spec)%string
          
            i_spec = i_spec + 1
          end do
          
          ! Move to the next phase instance
          i_phase = i_phase + 1
        
        end do

        deallocate(spec_names)

      end do
    end do

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species id on the \c pmc_phlex_state::phlex_state_t::state_var
  !! array by its unique name. These are unique ids for each element on the
  !! state array for this \ref phlex_aero_rep "aerosol representation" and
  !! are numbered:
  !!
  !!   \f[x_u \in x_f ... (x_f+n-1)\f]
  !!
  !! where \f$x_u\f$ is the id of the element corresponding to the species
  !! with unique name \f$u\f$ on the \c
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

  !> Get the non-unique name of a species by its unique name
  function spec_name(this, unique_name)

    !> Chemical species name
    character(len=:), allocatable :: spec_name
    !> Aerosol representation data
    class(aero_rep_modal_binned_mass_t), intent(in) :: this
    !> Unique name of the species in this aerosol representation
    character(len=:), allocatable :: unique_name

    ! Indices for iterators
    integer(kind=i_kind) :: i_spec, j_spec, i_phase
   
    ! species names in the aerosol phase
    type(string_t), allocatable :: spec_names(:)

    ! unique name list
    type(string_t), allocatable :: unique_names(:)

    unique_names = this%unique_names()

    i_spec = 1
    do i_phase = 1, size(this%aero_phase)
      spec_names = this%aero_phase(i_phase)%val%get_species_names()
      do j_spec = 1, this%aero_phase(i_phase)%val%size()
        if (unique_name.eq.unique_names(i_spec)%string) then
          spec_name = spec_names(j_spec)%string
        end if
        i_spec = i_spec + 1
      end do
      deallocate(spec_names)
    end do

    deallocate(unique_names)

  end function spec_name

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

  !> Finalize the aerosol representation
  elemental subroutine finalize(this)

    !> Aerosol representation data
    type(aero_rep_modal_binned_mass_t), intent(inout) :: this

    if (allocated(this%rep_name)) deallocate(this%rep_name)
    if (allocated(this%aero_phase)) then
      ! The core will deallocate the aerosol phases
      call this%aero_phase(:)%dereference()
      deallocate(this%aero_phase)
    end if
    if (associated(this%property_set)) &
            deallocate(this%property_set)
    if (allocated(this%condensed_data_real)) &
            deallocate(this%condensed_data_real)
    if (allocated(this%condensed_data_int)) &
            deallocate(this%condensed_data_int)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_rep_modal_binned_mass_GMD_t
  function update_data_GMD_constructor(aero_rep_type) result(new_obj)

    !> New update data object
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: new_obj
    !> Aerosol representation id
    integer(kind=i_kind), intent(in) :: aero_rep_type

    new_obj%aero_rep_type = int(aero_rep_type, kind=c_int)
    new_obj%update_data = aero_rep_modal_binned_mass_create_gmd_update_data()
    new_obj%is_malloced = .true.

  end function update_data_GMD_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for mode GMD
  subroutine update_data_set_GMD(this, aero_rep_id, section_id, GMD)

    !> Update data
    class(aero_rep_update_data_modal_binned_mass_GMD_t), intent(inout) :: this
    !> Aerosol representation id from
    !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::set_id
    integer(kind=i_kind), intent(in) :: aero_rep_id
    !> Aerosol section id from
    !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::get_section_id
    integer(kind=i_kind), intent(in) :: section_id
    !> Updated GMD (m)
    real(kind=dp), intent(in) :: GMD

    call aero_rep_modal_binned_mass_set_gmd_update_data(this%get_data(), &
            aero_rep_id, section_id, GMD)

  end subroutine update_data_set_GMD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a GMD update data object
  elemental subroutine update_data_GMD_finalize(this)

    !> Update data object to free
    type(aero_rep_update_data_modal_binned_mass_GMD_t), intent(inout) :: this

    if (this%is_malloced) call aero_rep_free_update_data(this%update_data)

  end subroutine update_data_GMD_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_rep_modal_binned_mass_GSD_t
  function update_data_GSD_constructor(aero_rep_type) result(new_obj)

    !> New update data object
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: new_obj
    !> Aerosol representation id
    integer(kind=i_kind), intent(in) :: aero_rep_type

    new_obj%aero_rep_type = int(aero_rep_type, kind=c_int)
    new_obj%update_data = aero_rep_modal_binned_mass_create_gsd_update_data()
    new_obj%is_malloced = .true.

  end function update_data_GSD_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set packed update data for mode GSD
  subroutine update_data_set_GSD(this, aero_rep_id, section_id, GSD)

    !> Update data
    class(aero_rep_update_data_modal_binned_mass_GSD_t), intent(inout) :: this
    !> Aerosol representation id from
    !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::set_id
    integer(kind=i_kind), intent(in) :: aero_rep_id
    !> Aerosol section id from
    !! pmc_aero_rep_modal_binned_mass::aero_rep_modal_binned_mass_t::get_section_id
    integer(kind=i_kind), intent(in) :: section_id
    !> Updated GSD (m)
    real(kind=dp), intent(in) :: GSD

    call aero_rep_modal_binned_mass_set_gsd_update_data(this%get_data(), &
            aero_rep_id, section_id, GSD)

  end subroutine update_data_set_GSD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize a GSD update data object
  elemental subroutine update_data_GSD_finalize(this)

    !> Update data object to free
    type(aero_rep_update_data_modal_binned_mass_GSD_t), intent(inout) :: this

    if (this%is_malloced) call aero_rep_free_update_data(this%update_data)

  end subroutine update_data_GSD_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef BINNED
#undef MODAL
#undef NUM_SECTION_
#undef INT_DATA_SIZE_
#undef FLOAT_DATA_SIZE_
#undef AERO_REP_ID_
#undef NUM_INT_PROP_
#undef NUM_REAL_PROP_
#undef MODE_INT_PROP_LOC_
#undef MODE_REAL_PROP_LOC_
#undef SECTION_TYPE_
#undef NUM_BINS_
#undef NUM_PHASE_
#undef PHASE_STATE_ID_
#undef PHASE_MODEL_DATA_ID_
#undef SPEC_STATE_ID_
#undef GMD_
#undef BIN_DP_
#undef GSD_
#undef NUMBER_CONC_
#undef EFFECTIVE_RADIUS_
#undef PHASE_MASS_
#undef PHASE_AVG_MW_

end module pmc_aero_rep_modal_binned_mass
