! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_single_particle module.

!> The abstract aero_rep_single_particle_t structure and associated subroutines.
module pmc_aero_rep_single_particle

  use pmc_util,                                  only: dp, i_kind, &
                                                       string_t, assert_msg, &
                                                       die_msg
  use pmc_property
  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_aero_phase_data
  use pmc_aero_rep_state
  use pmc_aero_rep_single_particle_state
  use pmc_phlex_state

  implicit none
  private

#define _NUM_PHASE_ this%condensed_data_int(1)
#define _AERO_STATE_ID_ this%condensed_data_int(2)
#define _NUM_INT_PROP_ 2
#define _NUM_REAL_PROP_ 0
#define _PHASE_STATE_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + x)
#define _SPEC_STATE_ID_(y,x) _PHASE_STATE_ID_(y) + x - 1
#define _PHASE_SPEC_ID_(x) this%condensed_data_int(_NUM_INT_PROP_ + _NUM_PHASE_ + x)
#define _DENSITY_(y,x) this%condensed_data_real(_PHASE_SPEC_ID_(y) + x - 1)
#define _MASS_(y,x) phlex_state%state_var(_PHASE_STATE_ID_(y) + x - 1)
  
  public :: aero_rep_single_particle_t

  !> Single particle aerosol representation
  !!
  !! Time-invariant data related to a single particle aerosol representation. 
  type, extends(aero_rep_data_t) :: aero_rep_single_particle_t
  contains
    !> Aerosol representation initialization
    procedure :: initialize
    !> Get the size of this representation on the state variable array
    procedure :: size => get_size
    !> Get a list of unique names for each element on the state array
    procedure :: unique_names
    !> Get a species state id by its unique name
    procedure :: state_id_by_unique_name
    !> Get an instance of the state variable for this aerosol representation
    procedure :: new_state
    !> Get an aerosol species state id
    procedure :: species_state_id
    !> Get the surface area concentration (m^2/m^3)
    procedure :: surface_area_conc
    !> Get the surface area concentration for a specific aerosol species
    !! (m^2/m^3)
    procedure :: species_surface_area_conc
    !> Get the vapor pressure scaling for a particular species (unitless)
    procedure :: vapor_pressure_scaling

    !> Private functions
    !> Get the associated aero_rep_state_t variable
    procedure, private :: get_state => get_state
  end type aero_rep_single_particle_t

  !> Constructor for aero_rep_single_particle_t
  interface aero_rep_single_particle_t
    procedure :: constructor
  end interface aero_rep_single_particle_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_rep_single_particle_t
  function constructor() result (new_obj)

    !> New aerosol representation
    type(aero_rep_single_particle_t), pointer :: new_obj

    allocate(new_obj)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the  property_t object. This 
  !! routine should be called once for each aerosol representation
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensures all data required during the model run are included
  !! in the condensed data arrays.
  subroutine initialize(this, aero_phase_set, &
                  spec_state_id, aero_state_id, chem_spec_data)

    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(inout) :: this
    !> The set of aerosol phases
    type(aero_phase_data_t), pointer, intent(in) :: aero_phase_set(:)
    !> Beginning state id for this aerosol representationin the model species
    !! state array
    integer(kind=i_kind), intent(in) :: spec_state_id
    !> Index for this representation in the model state aero_rep_state_t
    !! array
    integer(kind=i_kind), intent(in) :: aero_state_id
    !> Chemical species data
    type(chem_spec_data_t), intent(in) :: chem_spec_data

    integer(kind=i_kind) :: i_phase, i_spec, curr_id, curr_spec_id, num_spec
    type(property_t), pointer :: spec_props
    type(string_t), allocatable :: species(:)
    real(kind=dp) :: density
    character(len=:), allocatable :: key

    ! Assume all phases will be applied to each particle
    allocate(this%aero_phase(size(aero_phase_set)))
    do i_phase = 1, size(aero_phase_set)
      this%aero_phase(i_phase) = aero_phase_set(i_phase)
    end do

    ! Get the total number of species across all phases
    num_spec = 0
    do i_phase = 1, size(this%aero_phase)
      num_spec = num_spec + this%aero_phase(i_phase)%size()
    end do

    ! Allocate condensed data arrays
    allocate(this%condensed_data_int(_NUM_INT_PROP_ + 2*size(this%aero_phase)))
    allocate(this%condensed_data_real(_NUM_REAL_PROP_ + num_spec))

    ! Set indexes
    _NUM_PHASE_ = size(this%aero_phase)
    _AERO_STATE_ID_ = aero_state_id
    curr_id = spec_state_id
    curr_spec_id = _NUM_REAL_PROP_ + 1
    do i_phase = 1, _NUM_PHASE_
      _PHASE_STATE_ID_(i_phase) = curr_id
      _PHASE_SPEC_ID_(i_phase) = curr_spec_id
      curr_spec_id = curr_spec_id + this%aero_phase(i_phase)%size()
      curr_id = curr_id + this%aero_phase(i_phase)%size()
    end do

    ! Set densities
    key = "density"
    do i_phase = 1, _NUM_PHASE_
      curr_spec_id = 1
      species = this%aero_phase(i_phase)%get_species()
      do i_spec = 1, size(species)
        spec_props = chem_spec_data%get_property_set(species(i_spec)%string)
        if (.not.associated(spec_props)) then
          call die_msg(204001989, "Missing properties for species "// &
                  species(i_spec)%string)
        end if
        if (.not.spec_props%get_real(key, density)) then
          call die_msg(532333944, "Missing density for species "// &
                  species(i_spec)%string)
        end if
        _DENSITY_(i_phase, curr_spec_id) = density
        curr_spec_id = curr_spec_id + 1
      end do
    end do

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the size of this aerosol representation on the state variable array
  function get_size(this) result (state_size)

    !> Size on the state array
    integer(kind=i_kind) :: state_size
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this

    integer(kind=i_kind) :: i_phase

    ! Get the total number of species across all phases
    state_size = 0
    do i_phase = 1, size(this%aero_phase)
      state_size = state_size + this%aero_phase(i_phase)%size()
    end do

  end function get_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a list of unique names for each species on the state array
  function unique_names(this)

    !> List of unique names
    type(string_t), allocatable :: unique_names(:)
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this

    integer(kind=i_kind) :: num_spec, i_spec, i_phase
    type(string_t), allocatable :: phase_names(:)
    
    num_spec = 0
    do i_phase = 1, size(this%aero_phase)
      num_spec = num_spec + this%aero_phase(i_phase)%size()
    end do
    allocate(unique_names(num_spec))
    i_spec = 1
    do i_phase = 1, size(this%aero_phase)
      num_spec = this%aero_phase(i_phase)%size()
      phase_names = this%aero_phase(i_phase)%get_species()
      unique_names(i_spec:i_spec-1+num_spec) = phase_names(1:num_spec)
      i_spec = i_spec + num_spec
    end do

  end function unique_names

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get a species state id by its unique name
  function state_id_by_unique_name(this, unique_name) result (spec_id)

    !> Species state id
    integer(kind=i_kind) :: spec_id
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Unique name
    character(len=:), allocatable, intent(in) :: unique_name

    type(string_t), allocatable :: unique_names(:)
    integer(kind=i_kind) :: i_spec

    spec_id = 0
    unique_names = this%unique_names()
    do i_spec = 1, size(unique_names)
      if (unique_names(i_spec)%string .eq. unique_name) then
        spec_id = i_spec
        return
      end if
    end do

  end function state_id_by_unique_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an instance of the state type for this aerosol representation
  !!
  !! For PartMC single particle runs, the aerosol state will be set by PartMC
  !! at the beginning of each chemistry integration
  function new_state(this) result (aero_rep_state)

    !> Aerosol representation state
    class(aero_rep_state_t), pointer :: aero_rep_state
    !> Aerosol representaiton data
    class(aero_rep_single_particle_t), intent(in) :: this

    ! Empty state variable
    aero_rep_state => aero_rep_single_particle_state_t()

  end function new_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !> Get an aerosol species state id(s). The returned array will have the
 !! index of the species in the specified phase in each aerosol group (e.g.,
 !! bin, mode). If the species is not present in a certain group-phase the
 !! index will be 0.
 !!
 !! This function should only be called during initialization.
 !!
 !! One of each aerosol phase exists for each single particle, so this array
 !! will always be of length 1.
 function species_state_id(this, phase_name, species_name) result(spec_index)
    use pmc_util,                                     only : i_kind

    !> Species index array
    integer(kind=i_kind), allocatable :: spec_index(:)
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase id
    character(len=:), allocatable, intent(in) :: phase_name
    !> Species name
    character(len=:), allocatable, intent(in) :: species_name

    integer(kind=i_kind) :: i_phase, i_spec

    allocate(spec_index(1))
    i_phase = this%phase_id(phase_name)
    call assert_msg(820720954, i_phase.gt.0, "Invalid phase requested: "// &
            phase_name)
    i_spec = this%aero_phase(i_phase)%state_id(species_name)
    call assert_msg(413945507, i_spec.gt.0, "Invalid species requested: "// &
            species_name//" for phase: "//phase_name)
    spec_index(1) = _PHASE_STATE_ID_(i_phase) + i_spec

  end function species_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get surface area concentration (m^2/m^3) between two phases. One phase
  !! may be set to 0 to indicate the gas-phase, or both phases may be aerosol
  !! phases, with a corresponding index
  function surface_area_conc(this, i_phase1, i_phase2, phlex_state, &
                  jac_contrib)
    use pmc_util,                                     only : dp
    use pmc_constants

    !> Surface area concentration
    real(kind=dp) :: surface_area_conc
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase1 id
    integer(kind=i_kind), intent(in) :: i_phase1
    !> Aerosol phase2 id
    integer(kind=i_kind), intent(in) :: i_phase2
    !> Model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Contribution to Jacobian matrix. An array of the same size as the
    !! state array that, when present, will be filled with the partial
    !! derivatives of the result of this calculation with each state
    !! variable.
    real(kind=dp), allocatable, intent(out), optional :: jac_contrib(:)

    real(kind=dp) :: volume, radius
    integer(kind=i_kind) :: i_phase, i_spec

    ! Currently, only gas-aerosol interfaces are supported
    if (i_phase1.ne.0 .and. i_phase2.ne.0) then
      call die_msg(768138678, "Internal aerosol interface surface area not "//&
              "currently available.")
    end if
    i_phase = i_phase1 + i_phase2

    ! Calculate the volume density (m^3_aerosol/m^3_air)
    volume = 0
    do i_spec = 1, this%aero_phase(i_phase)%size()
      volume = volume + _DENSITY_(i_phase, i_spec) &
              * _MASS_(i_phase, i_spec)
    end do
    radius = (3.0/4.0*volume/const%pi)**(1/3)

    ! Calculate the surface area density (m^2/m^3)
    surface_area_conc = 3.0 * volume / radius

    ! Calculate jac_contrib
    ! TODO check math
    ! dSA/dx = 2*(3/4*V/pi)^(-1/3)*dV/dx
    ! dV/dx = density(x)
    if (present(jac_contrib)) then
      jac_contrib(:) = real(0.0, kind=dp)
      do i_spec = 1, this%aero_phase(i_phase)%size()
        jac_contrib(_SPEC_STATE_ID_(i_phase, i_spec)) = &
                2.0/((3.0/4.0*volume/const%pi)**(1/3)) * &
                _DENSITY_(i_phase, i_spec)
      end do
    end if

  end function surface_area_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the surface area concentration for a specific aerosol species
  !! (m^2/m^3). It is assumed the surface is between the gas-phase and an
  !! aerosol phase.
  function species_surface_area_conc(this,  i_phase, i_spec, phlex_state, &
                jac_contrib) result(surface_area_conc)
    use pmc_util,                                     only : dp

    !> Surface area concentration
    real(kind=dp) :: surface_area_conc
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase index
    integer(kind=i_kind), intent(in) :: i_phase
    !> Species id
    integer(kind=i_kind), intent(in) :: i_spec
    !> Model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Contribution to Jacobian matrix. An array of the same size as the
    !! state array that, when present, will be filled with the partial
    !! derivatives of the result of this calculation with each state
    !! variable.
    real(kind=dp), allocatable, intent(out), optional :: jac_contrib(:)

    real(kind=dp) :: mass_frac, surface_area
    integer(kind=i_kind) :: j_spec
    
    ! Get the mass fraction of this species (kg/kg)
    mass_frac = real(0.0, kind=dp)
    do j_spec = 1, this%aero_phase(i_phase)%size()
      mass_frac = mass_frac + _MASS_(i_phase, j_spec)
    end do
    if (mass_frac.gt.0.0) mass_frac = _MASS_(i_phase, i_spec) / mass_frac
    
    ! Get the surface area density (m^2/m^3)
    if (present(jac_contrib)) then
      surface_area = this%surface_area_conc(0, i_phase, phlex_state, jac_contrib)
    else
      surface_area = this%surface_area_conc(0, i_phase, phlex_state)
    end if

    ! get the surface area for this species (m^2/m^3)
    surface_area_conc = mass_frac * surface_area

    ! Calculate jac_contrib
    ! TODO check math
    ! dSSA/dx = dSA/dx * MF + SA * dMF/dx
    if (present(jac_contrib)) then
      do j_spec = 1, this%aero_phase(i_phase)%size()
        jac_contrib(_SPEC_STATE_ID_(i_phase, j_spec)) = &
                jac_contrib(_SPEC_STATE_ID_(i_phase, j_spec)) * mass_frac + &
                surface_area_conc * merge(0, 1, i_spec.eq.j_spec)
      end do
    end if

  end function species_surface_area_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
  !> Get the vapor pressure scaling for a particular species (unitless)
  function vapor_pressure_scaling(this, i_spec, phlex_state, jac_contrib)
    use pmc_util,                                     only : dp

    !> Vapor pressure scaling
    real(kind=dp) :: vapor_pressure_scaling
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Species id
    integer(kind=i_kind), intent(in) :: i_spec
    !> Model state
    type(phlex_state_t), intent(in) :: phlex_state
    !> Contribution to Jacobian matrix. An array of the same size as the
    !! state array that, when present, will be filled with the partial
    !! derivatives of the result of this calculation with each state
    !! variable.
    real(kind=dp), allocatable, intent(out), optional :: jac_contrib(:)

    ! TODO Finish

    call die_msg(787876225, "Kelvin effect not available yet.")
    
  end function vapor_pressure_scaling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the associated aero_rep_state_t variable
  function get_state(this, phlex_state) &
                  result (aero_rep_state)

    !> Aerosol representation state
    type(aero_rep_single_particle_state_t), pointer :: aero_rep_state
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Model state
    type(phlex_state_t), intent(in) :: phlex_state

    select type (state_ptr => phlex_state%aero_rep_state(_AERO_STATE_ID_)%val)
      type is (aero_rep_single_particle_state_t)
        aero_rep_state => state_ptr
      class default
        call die_msg(814359298, "Received incorrect aero_rep_state_t variable "//&
                "for aero_rep_single_particle_t")
    end select

  end function get_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_single_particle
