! Copyright (C) 2005-2012, 2016 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_particle module.

!> The aero_particle_t structure and associated subroutines.
module pmc_aero_particle

  use pmc_util
  use pmc_aero_data
  use pmc_spec_file
  use pmc_env_state
  use pmc_mpi
  use pmc_constants
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Single aerosol particle data structure.
  !!
  !! The \c vol array stores the total volumes of the different
  !! species that make up the particle. This array must have length
  !! equal to aero_data%%n_spec, so that \c vol(i) is the volume (in
  !! m^3) of the i'th aerosol species.
  type aero_particle_t
     !> Constituent species volumes [length aero_data_n_spec()] (m^3).
     real(kind=dp), allocatable :: vol(:)
     !> Number of original particles from each source that coagulated
     !> to form this one [length aero_data_n_source()].
     integer, allocatable :: n_orig_part(:)
     !> Weighting function group number (see \c aero_weight_array_t).
     integer :: weight_group
     !> Weighting function class number (see \c aero_weight_array_t).
     integer :: weight_class
     !> Absorption cross-section (m^2).
     real(kind=dp) :: absorb_cross_sect
     !> Scattering cross-section (m^2).
     real(kind=dp) :: scatter_cross_sect
     !> Asymmetry parameter (1).
     real(kind=dp) :: asymmetry
     !> Refractive index of the shell (1).
     complex(kind=dc) :: refract_shell
     !> Refractive index of the core (1).
     complex(kind=dc) :: refract_core
     !> Volume of the core (m^3).
     real(kind=dp) :: core_vol
     !> Water hysteresis curve section (0 = lower, 1 = upper)
     integer :: water_hyst_leg
     !> Unique ID number.
     integer :: id
     !> First time a constituent was created (s).
     real(kind=dp) :: least_create_time
     !> Last time a constituent was created (s).
     real(kind=dp) :: greatest_create_time
  end type aero_particle_t

  !> Next unique ID number to use for a particle.
  integer, save :: next_id = 1

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Shift data from one aero_particle_t to another and free the first
  !> one.
  subroutine aero_particle_shift(aero_particle_from, aero_particle_to)

    !> Reference particle (will be deallocated on return).
    type(aero_particle_t), intent(inout) :: aero_particle_from
    !> Destination particle (not allocated on entry).
    type(aero_particle_t), intent(inout) :: aero_particle_to

    call move_alloc(aero_particle_from%vol, aero_particle_to%vol)
    call move_alloc(aero_particle_from%n_orig_part, &
         aero_particle_to%n_orig_part)
    aero_particle_to%weight_group = aero_particle_from%weight_group
    aero_particle_to%weight_class = aero_particle_from%weight_class
    aero_particle_to%absorb_cross_sect = aero_particle_from%absorb_cross_sect
    aero_particle_to%scatter_cross_sect = &
         aero_particle_from%scatter_cross_sect
    aero_particle_to%asymmetry = aero_particle_from%asymmetry
    aero_particle_to%refract_shell = aero_particle_from%refract_shell
    aero_particle_to%refract_core = aero_particle_from%refract_core
    aero_particle_to%core_vol = aero_particle_from%core_vol
    aero_particle_to%water_hyst_leg = aero_particle_from%water_hyst_leg
    aero_particle_to%id = aero_particle_from%id
    aero_particle_to%least_create_time = aero_particle_from%least_create_time
    aero_particle_to%greatest_create_time = &
         aero_particle_from%greatest_create_time

  end subroutine aero_particle_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an aero_particle to be zero.
  subroutine aero_particle_zero(aero_particle, aero_data)

    !> Particle to zero.
    type(aero_particle_t), intent(inout) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    call ensure_real_array_size(aero_particle%vol, aero_data_n_spec(aero_data))
    aero_particle%vol = 0d0
    call ensure_integer_array_size(aero_particle%n_orig_part, &
         aero_data_n_source(aero_data))
    aero_particle%n_orig_part = 0
    aero_particle%weight_group = 0
    aero_particle%weight_class = 0
    aero_particle%absorb_cross_sect = 0d0
    aero_particle%scatter_cross_sect = 0d0
    aero_particle%asymmetry = 0d0
    aero_particle%refract_shell = (0d0, 0d0)
    aero_particle%refract_core = (0d0, 0d0)
    aero_particle%core_vol = 0d0
    aero_particle%water_hyst_leg = 0
    aero_particle%id = 0
    aero_particle%least_create_time = 0d0
    aero_particle%greatest_create_time = 0d0

  end subroutine aero_particle_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Assigns a globally-unique new ID number to the particle.
  subroutine aero_particle_new_id(aero_particle)

    !> Particle to set ID for.
    type(aero_particle_t), intent(inout) :: aero_particle

    aero_particle%id = (next_id - 1) * pmc_mpi_size() + pmc_mpi_rank() + 1
    next_id = next_id + 1

  end subroutine aero_particle_new_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the creation times for the particle.
  subroutine aero_particle_set_create_time(aero_particle, create_time)

    !> Particle to set time for.
    type(aero_particle_t), intent(inout) :: aero_particle
    !> Creation time.
    real(kind=dp), intent(in) :: create_time

    aero_particle%least_create_time = create_time
    aero_particle%greatest_create_time = create_time

  end subroutine aero_particle_set_create_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the aerosol particle volumes.
  subroutine aero_particle_set_vols(aero_particle, vols)

    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle
    !> New volumes.
    real(kind=dp), intent(in) :: vols(:)

    aero_particle%vol = vols

  end subroutine aero_particle_set_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the aerosol particle source.
  subroutine aero_particle_set_source(aero_particle, i_source)

    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle
    !> Source number for the particle.
    integer, intent(in) :: i_source

    aero_particle%n_orig_part = 0
    aero_particle%n_orig_part(i_source) = 1

  end subroutine aero_particle_set_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the aerosol particle weight group.
  subroutine aero_particle_set_weight(aero_particle, i_group, i_class)

    !> Particle.
    type(aero_particle_t), intent(inout) :: aero_particle
    !> Weight group number for the particle.
    integer, intent(in), optional :: i_group
    !> Weight class number for the particle.
    integer, intent(in), optional :: i_class

    if (present(i_group)) aero_particle%weight_group = i_group
    if (present(i_class)) aero_particle%weight_class = i_class

  end subroutine aero_particle_set_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total mass of the particle (kg).
  elemental real(kind=dp) function aero_particle_mass(aero_particle, aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_mass = sum(aero_particle%vol * aero_data%density)

  end function aero_particle_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mass of a single species in the particle (kg).
  elemental real(kind=dp) function aero_particle_species_mass(aero_particle, &
       i_spec, aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Species number to find mass of.
    integer, intent(in) :: i_spec
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_species_mass = aero_particle%vol(i_spec) &
         * aero_data%density(i_spec)

  end function aero_particle_species_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mass of all species in the particle (kg).
  function aero_particle_species_masses(aero_particle, aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Return value.
    real(kind=dp) :: aero_particle_species_masses(aero_data_n_spec(aero_data))

    aero_particle_species_masses = aero_particle%vol * aero_data%density

  end function aero_particle_species_masses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total moles in the particle (1).
  elemental real(kind=dp) function aero_particle_moles(aero_particle, &
       aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_moles = sum(aero_particle%vol * aero_data%density &
         / aero_data%molec_weight)

  end function aero_particle_moles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total volume of the particle (m^3).
  elemental real(kind=dp) function aero_particle_volume(aero_particle)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_particle_volume = sum(aero_particle%vol)

  end function aero_particle_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Volume of a single species in the particle (m^3).
  elemental real(kind=dp) function aero_particle_species_volume( &
       aero_particle, i_spec)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Species number to find volume of.
    integer, intent(in) :: i_spec

    aero_particle_species_volume = aero_particle%vol(i_spec)

  end function aero_particle_species_volume
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total dry volume of the particle (m^3).
  elemental real(kind=dp) function aero_particle_dry_volume(aero_particle, &
       aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle 
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: i_spec

    aero_particle_dry_volume = 0d0
    do i_spec = 1,aero_data_n_spec(aero_data)
       if (i_spec /= aero_data%i_water) then
          aero_particle_dry_volume = aero_particle_dry_volume &
               + aero_particle%vol(i_spec)
       end if
    end do

  end function aero_particle_dry_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total volume (dry or wet) of the particle (m^3).
  elemental real(kind=dp) function aero_particle_volume_maybe_dry( &
       aero_particle, aero_data, dry_volume)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether the desired volume is dry (otherwise wet).
    logical, intent(in) :: dry_volume

    if (dry_volume) then
       aero_particle_volume_maybe_dry &
            = aero_particle_dry_volume(aero_particle, aero_data)
    else
       aero_particle_volume_maybe_dry = aero_particle_volume(aero_particle)
    end if

  end function aero_particle_volume_maybe_dry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total radius of the particle (m).
  elemental real(kind=dp) function aero_particle_radius(aero_particle, &
       aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_radius = aero_data_vol2rad(aero_data, &
         aero_particle_volume(aero_particle))

  end function aero_particle_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total dry radius of the particle (m).
  elemental real(kind=dp) function aero_particle_dry_radius(aero_particle, &
       aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_dry_radius = aero_data_vol2rad(aero_data, &
         aero_particle_dry_volume(aero_particle, aero_data))

  end function aero_particle_dry_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total diameter of the particle (m).
  elemental real(kind=dp) function aero_particle_diameter(aero_particle, &
       aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_diameter = 2d0 * aero_particle_radius(aero_particle, &
         aero_data)

  end function aero_particle_diameter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Total dry diameter of the particle (m).
  elemental real(kind=dp) function aero_particle_dry_diameter(aero_particle, &
       aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_dry_diameter &
         = 2d0 * aero_particle_dry_radius(aero_particle, aero_data)

  end function aero_particle_dry_diameter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mobility diameter of the particle (m).
  real(kind=dp) function aero_particle_mobility_diameter(aero_particle, &
       aero_data, env_state)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: volume, mobility_radius

    volume = aero_particle_volume(aero_particle)
    mobility_radius = fractal_vol_to_mobility_rad(aero_data%fractal, &
         volume, env_state%temp, env_state%pressure)
    aero_particle_mobility_diameter = rad2diam(mobility_radius)

  end function aero_particle_mobility_diameter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Average density of the particle (kg/m^3).
  real(kind=dp) function aero_particle_density(aero_particle, aero_data)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_density = aero_particle_mass(aero_particle, aero_data) &
         / aero_particle_volume(aero_particle)

  end function aero_particle_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the volume-average of the non-water elements of quantity.
  real(kind=dp) function aero_particle_average_solute_quantity( &
       aero_particle, aero_data, quantity)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Quantity to average.
    real(kind=dp), intent(in) :: quantity(:)

    real(kind=dp) :: ones(aero_data_n_spec(aero_data))

    ones = 1d0
    aero_particle_average_solute_quantity = &
         aero_particle_total_solute_quantity(aero_particle, &
         aero_data, quantity) &
         / aero_particle_total_solute_quantity(aero_particle, aero_data, ones)

  end function aero_particle_average_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the volume-total of the non-water elements of quantity.
  real(kind=dp) function aero_particle_total_solute_quantity(aero_particle, &
       aero_data, quantity)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Quantity to total.
    real(kind=dp), intent(in) :: quantity(:)

    real(kind=dp) total
    integer i

    total = 0d0
    do i = 1,aero_data_n_spec(aero_data)
       if (i /= aero_data%i_water) then
          total = total + aero_particle%vol(i) * quantity(i)
       end if
    end do
    aero_particle_total_solute_quantity = total

  end function aero_particle_total_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the water element of quantity.
  real(kind=dp) function aero_particle_average_water_quantity(aero_particle, &
       aero_data, quantity)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Quantity to average.
    real(kind=dp), intent(in) :: quantity(:)

    call assert(420016623, aero_data%i_water > 0)
    aero_particle_average_water_quantity = quantity(aero_data%i_water)

  end function aero_particle_average_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the volume-total of the water element of quantity.
  real(kind=dp) function aero_particle_total_water_quantity(aero_particle, &
       aero_data, quantity)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Quantity to total.
    real(kind=dp), intent(in) :: quantity(:)

    call assert(223343210, aero_data%i_water > 0)
    aero_particle_total_water_quantity &
         = aero_particle%vol(aero_data%i_water) &
         * quantity(aero_data%i_water)

  end function aero_particle_total_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the water molecular weight.
  !> (kg/mole)
  real(kind=dp) function aero_particle_water_molec_weight(aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    call assert(772012490, aero_data%i_water > 0)
    aero_particle_water_molec_weight &
         = aero_data%molec_weight(aero_data%i_water)

  end function aero_particle_water_molec_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the average of the solute molecular weight (kg/mole).
  real(kind=dp) function aero_particle_solute_molec_weight(aero_particle, &
       aero_data)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_solute_molec_weight &
         = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, aero_data%molec_weight)

  end function aero_particle_solute_molec_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the average of the solute ion number (1).
  real(kind=dp) function aero_particle_solute_num_ions(aero_particle, &
       aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_particle_solute_num_ions &
         = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, real(aero_data%num_ions, kind=dp))

  end function aero_particle_solute_num_ions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the water density (kg/m^3).
  real(kind=dp) function aero_particle_water_density(aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    call assert(235482108, aero_data%i_water > 0)
    aero_particle_water_density = aero_data%density(aero_data%i_water)

  end function aero_particle_water_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the average of the solute densities (kg/m^3).
  real(kind=dp) function aero_particle_solute_density(aero_particle, &
       aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_particle_solute_density &
         = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, aero_data%density)

  end function aero_particle_solute_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the water mass (kg).
  real(kind=dp) function aero_particle_water_mass(aero_particle, aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle

    call assert(888636139, aero_data%i_water > 0)
    aero_particle_water_mass = aero_particle%vol(aero_data%i_water) &
         * aero_data%density(aero_data%i_water)

  end function aero_particle_water_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total solute mass (kg).
  real(kind=dp) function aero_particle_solute_mass(aero_particle, aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle

    aero_particle_solute_mass &
         = aero_particle_total_solute_quantity(aero_particle, &
         aero_data, aero_data%density)

  end function aero_particle_solute_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total solute volume (m^3).
  real(kind=dp) function aero_particle_solute_volume(aero_particle, &
       aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle

    real(kind=dp) :: ones(aero_data_n_spec(aero_data))

    ones = 1d0
    aero_particle_solute_volume &
         = aero_particle_total_solute_quantity(aero_particle, &
         aero_data, ones)

  end function aero_particle_solute_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total solute radius (m).
  real(kind=dp) function aero_particle_solute_radius(aero_particle, &
       aero_data)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    aero_particle_solute_radius &
         = aero_data_vol2rad(aero_data, &
         aero_particle_solute_volume(aero_particle, aero_data))

  end function aero_particle_solute_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the average of the solute kappas (1).
  real(kind=dp) function aero_particle_solute_kappa(aero_particle, aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle

    real(kind=dp) :: kappa(aero_data_n_spec(aero_data)), M_w, rho_w, M_a, rho_a
    integer :: i_spec

    M_w = aero_particle_water_molec_weight(aero_data)
    rho_w = aero_particle_water_density(aero_data)
    do i_spec = 1,aero_data_n_spec(aero_data)
       if (i_spec == aero_data%i_water) then
          kappa(i_spec) = 0d0
       elseif (aero_data%num_ions(i_spec) > 0) then
          call assert_msg(123681459, aero_data%kappa(i_spec) == 0d0, &
               "species has nonzero num_ions and kappa: " &
               // trim(aero_data%name(i_spec)))
          M_a = aero_data%molec_weight(i_spec)
          rho_a = aero_data%density(i_spec)
          kappa(i_spec) = M_w * rho_a / (M_a * rho_w) &
               * real(aero_data%num_ions(i_spec), kind=dp)
       else
          kappa(i_spec) = aero_data%kappa(i_spec)
       end if
    end do
    aero_particle_solute_kappa &
         = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, kappa)

  end function aero_particle_solute_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the approximate critical relative humidity (1).
  real(kind=dp) function aero_particle_approx_crit_rel_humid(aero_particle, &
       aero_data, env_state)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: kappa, diam, C, A

    A = env_state_A(env_state)
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    C = sqrt(4d0 * (A * const%water_surf_eng)**3 / 27d0)
    diam = aero_particle_diameter(aero_particle, aero_data)
    aero_particle_approx_crit_rel_humid = C / sqrt(kappa * diam**3) + 1d0

  end function aero_particle_approx_crit_rel_humid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the critical relative humidity (1).
  real(kind=dp) function aero_particle_crit_rel_humid(aero_particle, &
       aero_data, env_state)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: kappa, crit_diam, dry_diam, A

    A = env_state_A(env_state)
    dry_diam = aero_particle_dry_diameter(aero_particle, aero_data)
    crit_diam = aero_particle_crit_diameter(aero_particle, aero_data, &
         env_state)
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    if (kappa < 1d-30) then
       aero_particle_crit_rel_humid = exp(A * const%water_surf_eng / crit_diam)
    else
       aero_particle_crit_rel_humid = (crit_diam**3 - dry_diam**3) &
            / (crit_diam**3 - dry_diam**3 * (1 - kappa)) * exp(A * & 
            const%water_surf_eng / crit_diam)
    end if

  end function aero_particle_crit_rel_humid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! > Returns the critical relative humidity (1).
  real(kind=dp) function aero_particle_crit_rel_humid_varying_sigma(aero_particle, &
    aero_data, env_state)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp) :: A, dry_diam, kappa, d, sigma
    
    A = env_state_A(env_state)
    dry_diam = aero_particle_dry_diameter(aero_particle, aero_data)
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    d = aero_particle_crit_diameter_varying_sigma(aero_particle, aero_data, env_state, sigma)

    if (d == dry_diam) then
      aero_particle_crit_rel_humid_varying_sigma = exp(A * sigma / dry_diam)
    else 
      aero_particle_crit_rel_humid_varying_sigma = (d**3 - dry_diam**3) / & 
                  (d**3 - dry_diam**3 * (1d0 - kappa)) * exp(A * sigma / d) 
    end if     
 
  end function aero_particle_crit_rel_humid_varying_sigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the critical diameter (m).
  !!
  !! The method is as follows. We need to solve the polynomial
  !! equation \f$f(D)\f$ given in Eq. (7) of
  !!
  !! N. Riemer, M. West, R. Zaveri, D. Easter, "Estimating black
  !! carbon aging time-scales with a particle-resolved aerosol
  !! model", Aerosol Science 41, 143-158, 2010.
  !!
  !! This is the equation:
  !! \f[
  !!     f(D) = D^6 + c_4 D^4 + c_3 D^3 + c_0,
  !! \f]
  !! where \f$c_4 < 0\f$, \f$c_3 < 0\f$, and \f$c_0 > 0\f$. There is
  !! unique solution for \f$D > D_{\rm dry}\f$, as shown in the above
  !! paper.
  !!
  !! The polynomial has first derivative which factors like
  !! \f[
  !!     f'(D) = 6 D^5 + 4 c_4 D^3 + 3 c_3 D^2
  !!           = (3 D^2 + 4 c_4) D^3 + (D^3 + c_3) 3 D^2.
  !! \f]
  !! The first term is positive for \f$D > (-4 c_4 / 3)^{1/2}\f$ and
  !! the second is positive for \f$D > (-c_3)^{1/3}\f$. If we take
  !! \f[
  !!     D_0 = max((-4 c_4 / 3)^{1/2}, (-c_3)^{1/3})
  !! \f]
  !! then we have that \f$f'(D) > 0\f$ for \f$D > D_0\f$. Similarly,
  !! \f[
  !!     f''(D) = 30 D^4 + 12 c_4 D^2 + 6 c_3 D
  !!            = (5 D^2 + 4 c_4) 3 D^2 + (5 D^3 + 2 c_3) 3 D,
  !! \f]
  !! so \f$f''(D) > 0\f$ for \f$D > D_0\f$ (as this ensures that \f$D
  !! > (-4 c_4 / 5)^{1/2}\f$ and \f$D > (-2 c_3 / 5)^{1/3}\f$).
  !!
  !! Thus for $D > D_0$ we have that the first and second derivatives
  !! of $f(D)$ are positive, so Newton's method starting from
  !! \f$D_0\f$ will converge quadratically. This is the scheme used
  !! here.
  real(kind=dp) function aero_particle_crit_diameter(aero_particle, &
       aero_data, env_state)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    integer, parameter :: CRIT_DIAM_MAX_ITER = 100

    real(kind=dp) :: kappa, dry_diam, A, c4, c3, c0, d, f, df, dd
    integer :: i_newton

    A = env_state_A(env_state)
    dry_diam = aero_particle_dry_diameter(aero_particle, aero_data)
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    if (kappa < 1d-30) then
       ! bail out early for hydrophobic particles
       aero_particle_crit_diameter = dry_diam
       return
    end if

    c4 = - 3d0 * dry_diam**3 * kappa / (A * const%water_surf_eng)
    c3 = - dry_diam**3 * (2d0 - kappa)
    c0 = dry_diam**6 * (1d0 - kappa)

    ! Newton's method for f(d) = d^6 + c4 d^4 + c3 d^3 + c0
    d = max(sqrt(-4d0 / 3d0 * c4), (-c3)**(1d0/3d0))
    do i_newton = 1,CRIT_DIAM_MAX_ITER
       f = d**6 + c4 * d**4 + c3 * d**3 + c0
       df = 6 * d**5 + 4 * c4 * d**3 + 3 * c3 * d**2
       dd = f / df
       d = d - dd
       if (abs(dd / d) < 1d-14) then
          exit
       end if
    end do
    call warn_assert_msg(408545686, i_newton < CRIT_DIAM_MAX_ITER, &
         "critical diameter of old Newton loop failed to converge")
    call warn_assert_msg(353290871, d >= dry_diam, &
         "critical diameter of old Newton loop converged to invalid solution")
    aero_particle_crit_diameter = d

  end function aero_particle_crit_diameter
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! > Returns the critical diameter (m) for varying simga.
  real(kind=dp) function aero_particle_crit_diameter_varying_sigma(&
       aero_particle, aero_data, env_state, sigma)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp), intent(out) :: sigma

    integer, parameter :: CRIT_DIAM_MAX_ITER = 100
    integer :: i_newton
    real(kind=dp) :: delta_min = 1.6d-10
    real(kind=dp) :: kappa, dry_diam, A, T
    real(kind=dp) :: sigma_soluble, sigma_organic
    real(kind=dp) :: v_solid, v_sol, v_org, r_solid
    real(kind=dp) :: d, v_delta_solid, r_core, v_delta, d_v_delta
    real(kind=dp) :: d_sigma, d_2_sigma, R, d_R, f, df, dd
    real(kind=dp) :: c_1, c_2, c_3, c_4, c_5
    
    A = env_state_A(env_state)
    T = env_state%temp
    kappa = aero_particle_solute_kappa(aero_particle, aero_data)
    dry_diam = aero_particle_dry_diameter(aero_particle, aero_data)
    sigma_organic = aero_particle_sigma_organic(aero_particle, aero_data)
    v_solid = aero_particle_solid_volume(aero_particle, aero_data)
    v_org = aero_particle_organic_volume(aero_particle, aero_data)

    if (kappa < 1d-30) then
       ! bail out early for hydrophobic particles
      aero_particle_crit_diameter_varying_sigma = dry_diam
       return
    end if

    c_1 = 3d0 * dry_diam**3 * kappa / A 
    c_2= (2d0 - kappa) * dry_diam**3
    c_3 = (1d0 - kappa) * dry_diam**6
    c_4 = 2d0 * const%pi * delta_min ! d_2_v_delta_d 
  
    d = 50*sqrt(4d0 / 3d0 * c_1)

    if (v_solid == 0d0) then
      do i_newton = 1, CRIT_DIAM_MAX_ITER
        sigma_soluble = aero_particle_sigma_soluble(aero_particle, aero_data, env_state, d)
        c_5 = v_org * (sigma_organic - sigma_soluble)
        v_sol = const%pi * d**3 / 6d0 - v_org ! volume of soluble inorganics + water
        r_core = ((3d0 * v_sol)/(4d0 * const%pi))**(1d0 / 3d0)
        v_delta = (4d0 * const%pi / 3d0) * ((r_core + delta_min)**3 - (r_core)**3)
        if (v_org > v_delta) then
          sigma = sigma_organic
          f = d**6 - c_1 * d**4 / sigma - c_2 * d**3 + c_3  
          df = 6d0 * d**5 - 4d0 * c_1 * d**3 / sigma - 3d0 * c_2 * d**2
        else 
          if (v_org == 0) then
            sigma = sigma_soluble
            f = d**6 - c_1 * d**4 / sigma - c_2 * d**3 + c_3  
            df = 6d0 * d**5 - 4d0 * c_1 * d**3 / sigma - 3d0 * c_2 * d**2
          else
            sigma = sigma_soluble + c_5 / v_delta
            d_sigma = - c_5 * d_v_delta / v_delta**2
            d_v_delta = 2d0 * const%pi * delta_min * (d - delta_min)
            d_sigma = - c_5 * d_v_delta / v_delta**2
            R = sigma - d * d_sigma
            d_2_sigma = (c_5 / v_delta**3) * (2d0 * d_v_delta**2 - v_delta * c_4)
            d_R = - d * d_2_sigma
            f = R * (d**6 - c_2 * d**3 + c_3) - c_1 * d**4
            df = d_R * (d**6 - c_2 * d**3 + c_3) + R * (6d0 * d**5 - 3d0 * c_2 * d**2) &
                 - 4d0 * c_1 * d**3
          end if
        end if 
        dd = f / df
        d = d - dd
        if (abs(dd / d) < 1d-11) then
           exit
        end if
      end do
    else 
      r_solid = ((3d0 * v_solid) / (4d0 * const%pi))**(1d0 / 3d0)
      v_delta_solid = (4d0 * const%pi / 3d0) * ((r_solid + delta_min)**3 - (r_solid)**3)
      do i_newton = 1, CRIT_DIAM_MAX_ITER
        sigma_soluble = aero_particle_sigma_soluble(aero_particle, aero_data, env_state, d)
        c_5 = v_org * (sigma_organic - sigma_soluble)
        v_sol = const%pi * d**3 / 6 - v_solid - v_org
        if (v_sol + v_org > v_delta_solid) then 
          r_core = ((3d0 * (v_sol + v_solid))/(4 * const%pi))**(1d0 / 3d0) 
          v_delta = (4d0 * const%pi / 3d0) * ((r_core + delta_min)**3 - (r_core)**3)
          if (v_org > v_delta) then  
            sigma = sigma_organic
            f = d**6 - c_1 * d**4 / sigma - c_2 * d**3 + c_3  
            df = 6d0 * d**5 - 4d0 * c_1 * d**3 / sigma - 3d0 * c_2 * d**2
          else
            if (v_org == 0) then
              sigma = sigma_soluble
              f = d**6 - c_1 * d**4 / sigma - c_2 * d**3 + c_3  
              df = 6d0 * d**5 - 4d0 * c_1 * d**3 / sigma - 3d0 * c_2 * d**2
            else
              sigma = sigma_soluble + c_5 / v_delta
              d_sigma = - c_5 * d_v_delta / v_delta**2
              d_v_delta = 2d0 * const%pi * delta_min * (d - delta_min)
              d_sigma = - c_5 * d_v_delta / v_delta**2
              R = sigma - d * d_sigma
              d_2_sigma = (c_5 / v_delta**3) * (2d0 * d_v_delta**2 - v_delta * c_4)
              d_R = - d * d_2_sigma
              f = R * (d**6 - c_2 * d**3 + c_3) - c_1 * d**4
              df = d_R * (d**6 - c_2 * d**3 + c_3) + R * (6d0 * d**5 - 3d0 * c_2 * d**2) &
                   - 4d0 * c_1 * d**3
            end if
          end if
          dd = f / df
          d = d - dd
          if (abs(dd / d) < 1d-11) then
            exit
          end if
        else
          sigma = (v_org * sigma_organic + v_sol * sigma_soluble) / v_delta_solid
          d = dry_diam
          exit
        end if
      end do
    end if

    call warn_assert_msg(408545686, i_newton < CRIT_DIAM_MAX_ITER, &
         "critical diameter for new Newton loop failed to converge")
    call warn_assert_msg(353290871, d >= dry_diam, &
         "critical diameter for new Newton loop converged to invalid solution")
    aero_particle_crit_diameter_varying_sigma = d
    
  end function aero_particle_crit_diameter_varying_sigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Organic volume of a single species in the particle (m^3). 
  real(kind=dp) function aero_particle_organic_volume(aero_particle, aero_data) 

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    integer :: i_org_spec, n_org_spec
    character(len=AERO_NAME_LEN), parameter, dimension(10) :: &
    org_spec = ["MSA   ", "ARO1  ", "ARO2  ", "ALK1  ", "OLE1  ", &
                "API1  ", "API2  ", "LIM1  ", "LIM2  ", "OC    "]

    aero_particle_organic_volume = 0d0 
    
    do n_org_spec = 1, size(org_spec)
       i_org_spec = aero_data_spec_by_name(aero_data, org_spec(n_org_spec))
       aero_particle_organic_volume = aero_particle_organic_volume &
            + aero_particle%vol(i_org_spec)
    end do

  end function aero_particle_organic_volume
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! > soluble inrganic volume of a single species in the particle (m^3). 
  real(kind=dp) function aero_particle_solid_volume(aero_particle, aero_data) 

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.ls
    type(aero_data_t), intent(in) :: aero_data
    integer :: i_solid_spec, n_solid_spec
    character(len=AERO_NAME_LEN), parameter, dimension(2) :: &
    solid_spec = ["OIN   ", "BC    "]
     
    aero_particle_solid_volume = 0d0 

    do n_solid_spec = 1, size(solid_spec)
       i_solid_spec = aero_data_spec_by_name(aero_data, solid_spec(n_solid_spec))
       aero_particle_solid_volume = aero_particle_solid_volume &
            + aero_particle%vol(i_solid_spec)
    end do

  end function aero_particle_solid_volume
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the sigma of sol.
  real(kind=dp) function aero_particle_sigma_soluble(aero_particle, & 
              aero_data, env_state, d)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    
    integer :: i_sol_spec, n_sol_spec
    real(kind=dp) :: d, aero_particle_soluble_volume
    real(kind=dp) :: v_solid, v_org, v_sol, v_water
    character(len=AERO_NAME_LEN), parameter, dimension(7) :: &
    sol_spec = ["SO4   ", "NO3   ", "Cl    ", "NH4   ", "CO3   ", &
                "Na    ", "Ca    "]

    v_solid = aero_particle_solid_volume(aero_particle, aero_data)
    v_org = aero_particle_organic_volume(aero_particle, aero_data)

    aero_particle_soluble_volume = 0d0 

    do n_sol_spec = 1, size(sol_spec)
       i_sol_spec = aero_data_spec_by_name(aero_data, sol_spec(n_sol_spec))
       aero_particle_soluble_volume = aero_particle_soluble_volume &
            + aero_particle%vol(i_sol_spec)
    end do
  
    v_sol = d**3 * const%pi / 6 - v_solid - v_org
    v_water = v_sol - aero_particle_soluble_volume

    aero_particle_sigma_soluble = 0d0
    do n_sol_spec = 1, size(sol_spec)
      i_sol_spec = aero_data_spec_by_name(aero_data, sol_spec(n_sol_spec))
      aero_particle_sigma_soluble = aero_particle_sigma_soluble & 
              + aero_particle%vol(i_sol_spec) * &
              aero_data%sigma(i_sol_spec) / v_sol
    end do
    aero_particle_sigma_soluble = aero_particle_sigma_soluble + v_water &
              * const%water_surf_eng / v_sol

  end function aero_particle_sigma_soluble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !> Returns the sigma of shell
  real(kind=dp) function aero_particle_sigma_organic(aero_particle, aero_data)

    !> Aerosol particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    real(kind=dp) :: org_volume
    integer :: i_org_spec, n_org_spec
    character(len=AERO_NAME_LEN), parameter, dimension(10) :: &
    org_spec = ["MSA   ", "ARO1  ", "ARO2  ", "ALK1  ", "OLE1  ", &
                "API1  ", "API2  ", "LIM1  ", "LIM2  ", "OC    "]    

    org_volume = aero_particle_organic_volume(aero_particle, aero_data)
    aero_particle_sigma_organic = 0d0
    
    do n_org_spec = 1, size(org_spec)
      i_org_spec = aero_data_spec_by_name(aero_data, org_spec(n_org_spec))
      aero_particle_sigma_organic = aero_particle_sigma_organic + aero_particle%vol(i_org_spec) * &
              aero_data%sigma(i_org_spec) / org_volume
    end do 

  end function aero_particle_sigma_organic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Coagulate two particles together to make a new one. The new
  !> particle will not have its ID set.
  subroutine aero_particle_coagulate(aero_particle_1, &
       aero_particle_2, aero_particle_new)

    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Result particle.
    type(aero_particle_t), intent(inout) :: aero_particle_new

    call assert(203741686, size(aero_particle_1%vol) &
         == size(aero_particle_2%vol))
    call assert(483452167, size(aero_particle_1%n_orig_part) &
         == size(aero_particle_2%n_orig_part))
    aero_particle_new%vol = aero_particle_1%vol + aero_particle_2%vol
    aero_particle_new%n_orig_part = aero_particle_1%n_orig_part &
         + aero_particle_2%n_orig_part
    aero_particle_new%weight_group = 0
    aero_particle_new%weight_class = 0
    aero_particle_new%absorb_cross_sect = 0d0
    aero_particle_new%scatter_cross_sect = 0d0
    aero_particle_new%asymmetry = 0d0
    aero_particle_new%refract_shell = (0d0, 0d0)
    aero_particle_new%refract_core = (0d0, 0d0)
    aero_particle_new%core_vol = 0d0
    if ((aero_particle_1%water_hyst_leg == 1) &
         .and. (aero_particle_2%water_hyst_leg == 1)) then
       aero_particle_new%water_hyst_leg = 1
    else
       aero_particle_new%water_hyst_leg = 0
    end if
    aero_particle_new%id = 0
    aero_particle_new%least_create_time = &
         min(aero_particle_1%least_create_time, &
         aero_particle_2%least_create_time)
    aero_particle_new%greatest_create_time = &
         max(aero_particle_1%greatest_create_time, &
         aero_particle_2%greatest_create_time)

  end subroutine aero_particle_coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_particle(val)

    !> Value to pack.
    type(aero_particle_t), intent(in) :: val

    pmc_mpi_pack_size_aero_particle = &
         pmc_mpi_pack_size_real_array(val%vol) &
         + pmc_mpi_pack_size_integer_array(val%n_orig_part) &
         + pmc_mpi_pack_size_integer(val%weight_group) &
         + pmc_mpi_pack_size_integer(val%weight_class) &
         + pmc_mpi_pack_size_real(val%absorb_cross_sect) &
         + pmc_mpi_pack_size_real(val%scatter_cross_sect) &
         + pmc_mpi_pack_size_real(val%asymmetry) &
         + pmc_mpi_pack_size_complex(val%refract_shell) &
         + pmc_mpi_pack_size_complex(val%refract_core) &
         + pmc_mpi_pack_size_real(val%core_vol) &
         + pmc_mpi_pack_size_integer(val%water_hyst_leg) &
         + pmc_mpi_pack_size_integer(val%id) &
         + pmc_mpi_pack_size_real(val%least_create_time) &
         + pmc_mpi_pack_size_real(val%greatest_create_time)

  end function pmc_mpi_pack_size_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_particle(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_particle_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%vol)
    call pmc_mpi_pack_integer_array(buffer, position, val%n_orig_part)
    call pmc_mpi_pack_integer(buffer, position, val%weight_group)
    call pmc_mpi_pack_integer(buffer, position, val%weight_class)
    call pmc_mpi_pack_real(buffer, position, val%absorb_cross_sect)
    call pmc_mpi_pack_real(buffer, position, val%scatter_cross_sect)
    call pmc_mpi_pack_real(buffer, position, val%asymmetry)
    call pmc_mpi_pack_complex(buffer, position, val%refract_shell)
    call pmc_mpi_pack_complex(buffer, position, val%refract_core)
    call pmc_mpi_pack_real(buffer, position, val%core_vol)
    call pmc_mpi_pack_integer(buffer, position, val%water_hyst_leg)
    call pmc_mpi_pack_integer(buffer, position, val%id)
    call pmc_mpi_pack_real(buffer, position, val%least_create_time)
    call pmc_mpi_pack_real(buffer, position, val%greatest_create_time)
    call assert(810223998, position - prev_position &
         <= pmc_mpi_pack_size_aero_particle(val))
#endif

  end subroutine pmc_mpi_pack_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_particle(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_particle_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%vol)
    call pmc_mpi_unpack_integer_array(buffer, position, val%n_orig_part)
    call pmc_mpi_unpack_integer(buffer, position, val%weight_group)
    call pmc_mpi_unpack_integer(buffer, position, val%weight_class)
    call pmc_mpi_unpack_real(buffer, position, val%absorb_cross_sect)
    call pmc_mpi_unpack_real(buffer, position, val%scatter_cross_sect)
    call pmc_mpi_unpack_real(buffer, position, val%asymmetry)
    call pmc_mpi_unpack_complex(buffer, position, val%refract_shell)
    call pmc_mpi_unpack_complex(buffer, position, val%refract_core)
    call pmc_mpi_unpack_real(buffer, position, val%core_vol)
    call pmc_mpi_unpack_integer(buffer, position, val%water_hyst_leg)
    call pmc_mpi_unpack_integer(buffer, position, val%id)
    call pmc_mpi_unpack_real(buffer, position, val%least_create_time)
    call pmc_mpi_unpack_real(buffer, position, val%greatest_create_time)
    call assert(287447241, position - prev_position &
         <= pmc_mpi_pack_size_aero_particle(val))
#endif

  end subroutine pmc_mpi_unpack_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the particle data is consistent.
  subroutine aero_particle_check(aero_particle, aero_data, &
       continue_on_error)

    !> Aerosol particle to check.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Whether to continue despite error.
    logical, intent(in) :: continue_on_error

    if (allocated(aero_particle%vol)) then
       if (size(aero_particle%vol) /= aero_data_n_spec(aero_data)) then
          write(0, *) 'ERROR aero_particle A:'
          write(0, *) 'size(aero_particle%vol)', size(aero_particle%vol)
          write(0, *) 'aero_data_n_spec(aero_data)', &
               aero_data_n_spec(aero_data)
          call assert(185878626, continue_on_error)
       end if
    end if
    if (allocated(aero_particle%n_orig_part)) then
       if (size(aero_particle%n_orig_part) &
            /= aero_data_n_source(aero_data)) then
          write(0, *) 'ERROR aero_particle A:'
          write(0, *) 'size(aero_particle%n_orig_part)', &
               size(aero_particle%n_orig_part)
          write(0, *) 'aero_data_n_source(aero_data)', &
               aero_data_n_source(aero_data)
          call assert(625490639, continue_on_error)
       end if
    end if

  end subroutine aero_particle_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_particle
