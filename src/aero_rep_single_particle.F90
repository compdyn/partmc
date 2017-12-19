! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_rep_single_particle module.

!> The abstract aero_rep_single_particle_t structure and associated subroutines.
module pmc_aero_rep_single_particle

  use pmc_chem_spec_data
  use pmc_aero_rep_data
  use pmc_aero_phase_data

  implicit none
  private

  public :: aero_rep_single_particle_t

  !> Single particle aerosol representation
  !!
  !! Time-invariant data related to a single particle aerosol representation. 
  type, extends(aero_rep_data_t) :: aero_rep_single_particle_t
  contains
    !> Aerosol representation initialization
    procedure :: initialize => pmc_aero_rep_single_particle_initialize
    !> Get an aerosol species state id
    procedure :: species_state_id => pmc_aero_rep_single_particle_species_state_id
    !> Get the surface area concentration (m^2/m^3)
    procedure :: surface_area_conc => &
            pmc_aero_rep_single_particle_surface_area_conc
    !> Get the surface area concentration for a specific aerosol species
    !! (m^2/m^3)
    procedure :: species_surface_area_conc => &
            pmc_aero_rep_single_particle_species_surface_area_conc
    !> Get the vapor pressure scaling for a particular species (unitless)
    procedure :: vapor_pressure_scaling => &
            pmc_aero_rep_single_particle_vapor_pressure_scaling
  end type aero_rep_single_particle_t

  !> Constructor for aero_rep_single_particle_t
  interface aero_rep_single_particle_t
    procedure :: pmc_aero_rep_single_particle_constructor
  end interface aero_rep_single_particle_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for aero_rep_single_particle_t
  function pmc_aero_rep_single_particle_constructor() result (new_obj)

    !> New aerosol representation
    type(aero_rep_single_particle_t), pointer :: new_obj

    allocate(new_obj)

  end function pmc_aero_rep_single_particle_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the aerosol representation data, validating component data and
  !! loading any required information from the  property_t object. This 
  !! routine should be called once for each aerosol representation
  !! at the beginning of a model run after all the input files have been
  !! read in. It ensures all data required during the model run are included
  !! in the condensed data arrays.
  subroutine pmc_aero_rep_single_particle_initialize(this, aero_phase_set)

    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(inout) :: this
    !> Chemical species data
    type(aero_phase_data_t), pointer, intent(in) :: aero_phase_set(:)

  end subroutine pmc_aero_rep_single_particle_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get an aerosol species state id(s). The returned array will have the
  !! index of the species in the specified phase in each aerosol group (e.g.,
  !! bin, mode). If the species is not present in a certain group-phase the
  !! index will be 0.
 function pmc_aero_rep_single_particle_species_state_id(this, phase, &
                    species_name) result(spec_index)
    use pmc_util,                                     only : i_kind

    !> Species index array
    integer(kind=i_kind), allocatable :: spec_index(:)
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase
    character(len=:), allocatable, intent(in) :: phase
    !> Species name
    character(len=:), allocatable, intent(in) :: species_name

    ! TODO Finish
    allocate(spec_index(1))

  end function pmc_aero_rep_single_particle_species_state_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get surface area concentration (m^2/m^3) between two phases. One phase
  !! may be 'gas' to indicate the gas-phase, or both phases may be aerosol
  !! phases.
  function pmc_aero_rep_single_particle_surface_area_conc(this, phase1, phase2) &
                    result(surface_area_conc)
    use pmc_util,                                     only : dp

    !> Surface area concentration
    real(kind=dp) :: surface_area_conc
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase1
    character(len=:), allocatable, intent(in) :: phase1
    !> Aerosol phase2
    character(len=:), allocatable, intent(in) :: phase2

    ! TODO Finish
    surface_area_conc = real(0.0, kind=dp)

  end function pmc_aero_rep_single_particle_surface_area_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the surface area concentration for a specific aerosol species
  !! (m^2/m^3). It is assumed the surface is between the gas-phase and an
  !! aerosol phase.
  function pmc_aero_rep_single_particle_species_surface_area_conc(this, phase, species_name) &
                    result(surface_area_conc)
    use pmc_util,                                     only : dp

    !> Surface area concentration
    real(kind=dp) :: surface_area_conc
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Aerosol phase
    character(len=:), allocatable, intent(in) :: phase
    !> Species name
    character(len=:), allocatable, intent(in) :: species_name

    ! TODO Finish
    surface_area_conc = real(0.0, kind=dp)

  end function pmc_aero_rep_single_particle_species_surface_area_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
  !> Get the vapor pressure scaling for a particular species (unitless)
  function pmc_aero_rep_single_particle_vapor_pressure_scaling(this, species_name) &
                    result(vapor_pressure_scaling)
    use pmc_util,                                     only : dp

    !> Vapor pressure scaling
    real(kind=dp) :: vapor_pressure_scaling
    !> Aerosol representation data
    class(aero_rep_single_particle_t), intent(in) :: this
    !> Species name
    character(len=:), allocatable, intent(in) :: species_name

    ! TODO Finish
    vapor_pressure_scaling = real(0.0, kind=dp)

  end function pmc_aero_rep_single_particle_vapor_pressure_scaling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_rep_single_particle
