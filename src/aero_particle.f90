! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Basic particle structure.

module pmc_aero_particle

  use pmc_util
  use pmc_aero_data
  use pmc_bin_grid
  use pmc_inout
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

  type aero_particle_t
     real*8, pointer :: vol(:)           ! constituent species volumes (m^3)
     integer :: n_orig_part              ! number of original particles
     real*8 :: absorb_cross_sect         ! absorption cross-section (m^2)
     real*8 :: scatter_cross_sect        ! scattering cross-section (m^2)
     real*8 :: asymmetry                 ! asymmetry parameter (1)
     complex*16 :: refract_shell         ! refractive index of the shell (1)
     complex*16 :: refract_core          ! refractive index of the core (1)
     real*8 :: core_vol                  ! volume of the core (m^3)
  end type aero_particle_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_alloc(aero_particle, n_spec)

    ! Allocates and initializes.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle to init
    integer, intent(in) :: n_spec       ! number of species

    allocate(aero_particle%vol(n_spec))
    call aero_particle_zero(aero_particle)

  end subroutine aero_particle_alloc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_free(aero_particle)

    ! Deallocates.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle to free
    
    deallocate(aero_particle%vol)

  end subroutine aero_particle_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_copy(aero_particle_from, aero_particle_to)

    ! Copies a particle.

    type(aero_particle_t), intent(in) :: aero_particle_from ! reference particle
    type(aero_particle_t), intent(inout) :: aero_particle_to ! already allocated
    
    integer :: n_spec

    n_spec = size(aero_particle_from%vol)
    if (n_spec /= size(aero_particle_to%vol)) then
       call aero_particle_free(aero_particle_to)
       call aero_particle_alloc(aero_particle_to, n_spec)
    end if
    call assert(651178226, size(aero_particle_from%vol) &
         == size(aero_particle_to%vol))
    aero_particle_to%vol = aero_particle_from%vol
    aero_particle_to%n_orig_part = aero_particle_from%n_orig_part
    aero_particle_to%absorb_cross_sect = aero_particle_from%absorb_cross_sect
    aero_particle_to%scatter_cross_sect = aero_particle_from%scatter_cross_sect
    aero_particle_to%asymmetry = aero_particle_from%asymmetry
    aero_particle_to%refract_shell = aero_particle_from%refract_shell
    aero_particle_to%refract_core = aero_particle_from%refract_core
    aero_particle_to%core_vol = aero_particle_from%core_vol

  end subroutine aero_particle_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_shift(aero_particle_from, aero_particle_to)

    ! Is the equivalent of aero_particle_copy(from, to) followed by
    ! aero_particle_free(from), but faster.

    type(aero_particle_t), intent(in) :: aero_particle_from ! reference particle
    type(aero_particle_t), intent(inout) :: aero_particle_to ! not allocated

    aero_particle_to%vol => aero_particle_from%vol
    nullify(aero_particle_from%vol)
    aero_particle_to%n_orig_part = aero_particle_from%n_orig_part
    aero_particle_to%absorb_cross_sect = aero_particle_from%absorb_cross_sect
    aero_particle_to%scatter_cross_sect = aero_particle_from%scatter_cross_sect
    aero_particle_to%asymmetry = aero_particle_from%asymmetry
    aero_particle_to%refract_shell = aero_particle_from%refract_shell
    aero_particle_to%refract_core = aero_particle_from%refract_core
    aero_particle_to%core_vol = aero_particle_from%core_vol
    
  end subroutine aero_particle_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_zero(aero_particle)

    ! Resets an aero_particle to be zero.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle to zero
    
    aero_particle%vol = 0d0
    aero_particle%n_orig_part = 1
    aero_particle%absorb_cross_sect = 0d0
    aero_particle%scatter_cross_sect = 0d0
    aero_particle%asymmetry = 0d0
    aero_particle%refract_shell = (0d0, 0d0)
    aero_particle%refract_core = (0d0, 0d0)
    aero_particle%core_vol = 0d0

  end subroutine aero_particle_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_set_vols(aero_particle, vols)

    ! Sets the aerosol particle volumes.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle
    real*8, intent(in) :: vols(size(aero_particle%vol)) ! new volumes

    aero_particle%vol = vols

  end subroutine aero_particle_set_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_mass(aero_particle, aero_data) ! (kg)

    ! Total mass of the particle.

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    
    aero_particle_mass = sum(aero_particle%vol * aero_data%density)

  end function aero_particle_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_moles(aero_particle, aero_data) ! (1)

    ! Total moles in the particle.

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    
    aero_particle_moles = sum(aero_particle%vol * aero_data%density &
         / aero_data%molec_weight)

  end function aero_particle_moles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_volume(aero_particle) ! m^3

    ! Total volume of the particle.

    type(aero_particle_t), intent(in) :: aero_particle ! particle

    aero_particle_volume = sum(aero_particle%vol)

  end function aero_particle_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer function aero_particle_in_bin(aero_particle, bin_grid)
    
    ! Find the bin number that contains a given particle.

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    
    aero_particle_in_bin = bin_grid_particle_in_bin(bin_grid, &
         aero_particle_volume(aero_particle))
    
  end function aero_particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_average_solute_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the volume-average of the non-water elements of quantity.

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to average

    real*8 :: ones(aero_data%n_spec)

    ones = 1d0
    aero_particle_average_solute_quantity = &
         aero_particle_total_solute_quantity(aero_particle, &
         aero_data, quantity) &
         / aero_particle_total_solute_quantity(aero_particle, aero_data, ones)

  end function aero_particle_average_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_total_solute_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the volume-total of the non-water elements of quantity.

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to total

    real*8 total
    integer i

    total = 0d0
    do i = 1,aero_data%n_spec
       if (i /= aero_data%i_water) then
          total = total + aero_particle%vol(i) * quantity(i)
       end if
    end do
    aero_particle_total_solute_quantity = total

  end function aero_particle_total_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_water_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the water element of quantity.

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to average

    call assert(420016623, aero_data%i_water > 0)
    average_water_quantity = quantity(aero_data%i_water)

  end function average_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_water_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the volume-total of the water element of quantity.

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to total

    call assert(223343210, aero_data%i_water > 0)
    total_water_quantity = aero_particle%vol(aero_data%i_water) &
         * quantity(aero_data%i_water)

  end function total_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_molec_weight(aero_data) ! (kg/mole)

    ! Returns the water molecular weight.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data

    call assert(772012490, aero_data%i_water > 0)
    aero_particle_water_molec_weight = aero_data%molec_weight(aero_data%i_water)

  end function aero_particle_water_molec_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_molec_weight(aero_particle, &
       aero_data) ! (kg/mole)

    ! Returns the average of the solute molecular weight.

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data ! aerosol data

    aero_particle_solute_molec_weight = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, aero_data%molec_weight)

  end function aero_particle_solute_molec_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_num_ions(aero_particle, aero_data) ! (1)

    ! Returns the average of the solute ion number.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_num_ions = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, dble(aero_data%num_ions))

  end function aero_particle_solute_num_ions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_solubility(aero_particle, &
       aero_data) ! (1)

    ! Returns the average of the solute solubilities.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_solubility = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, aero_data%solubility)

  end function aero_particle_solute_solubility

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_density(aero_data) ! (kg/m^3)

    ! Returns the water density.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data

    call assert(235482108, aero_data%i_water > 0)
    aero_particle_water_density = aero_data%density(aero_data%i_water)

  end function aero_particle_water_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_density(aero_particle, &
       aero_data) ! (kg/m^3)

    ! Returns the average of the solute densities.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_density = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, aero_data%density)

  end function aero_particle_solute_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_mass(aero_particle, aero_data) ! (kg)

    ! Returns the water mass.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    call assert(888636139, aero_data%i_water > 0)
    aero_particle_water_mass = aero_particle%vol(aero_data%i_water) &
         * aero_data%density(aero_data%i_water)

  end function aero_particle_water_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_mass(aero_particle, aero_data) ! (kg)

    ! Returns the total solute mass.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_mass = aero_particle_total_solute_quantity(aero_particle, &
         aero_data, aero_data%density)

  end function aero_particle_solute_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_kappa(aero_particle, aero_data) ! (1)

    ! Returns the average of the solute kappas.

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_kappa = aero_particle_average_solute_quantity(aero_particle, &
         aero_data, aero_data%kappa)

  end function aero_particle_solute_kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_coagulate(aero_particle_1, &
       aero_particle_2, aero_particle_new)

    ! Coagulate two particles together to make a new one.

    type(aero_particle_t), intent(in) :: aero_particle_1 ! first particle
    type(aero_particle_t), intent(in) :: aero_particle_2 ! second particle
    type(aero_particle_t), intent(inout) :: aero_particle_new ! result particle

    call assert(203741686, size(aero_particle_1%vol) &
         == size(aero_particle_new%vol))
    call assert(586181003, size(aero_particle_2%vol) &
         == size(aero_particle_new%vol))
    aero_particle_new%vol = aero_particle_1%vol + aero_particle_2%vol
    aero_particle_new%n_orig_part = aero_particle_1%n_orig_part &
         + aero_particle_2%n_orig_part
    aero_particle_new%absorb_cross_sect = 0d0
    aero_particle_new%scatter_cross_sect = 0d0
    aero_particle_new%asymmetry = 0d0
    aero_particle_new%refract_shell = (0d0, 0d0)
    aero_particle_new%refract_core = (0d0, 0d0)
    aero_particle_new%core_vol = 0d0

  end subroutine aero_particle_coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_write_aero_particle(file, aero_particle)
    
    ! Write full state.
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_t), intent(in) :: aero_particle ! aero_particle to write
    
    call inout_write_comment(file, "begin aero_particle")
    call inout_write_integer(file, "n_orig_part", aero_particle%n_orig_part)
    call inout_write_real(file, "absorb(m^2)", aero_particle%absorb_cross_sect)
    call inout_write_real(file, "scatter(m^2)", &
         aero_particle%scatter_cross_sect)
    call inout_write_real(file, "asymmetry(1)", aero_particle%asymmetry)
    call inout_write_complex(file, "refract_shell(1)", &
         aero_particle%refract_shell)
    call inout_write_complex(file, "refract_core(1)", &
         aero_particle%refract_core)
    call inout_write_real(file, "core_vol(m^3)", aero_particle%core_vol)
    call inout_write_real_array(file, "spec_vols(m^3)", aero_particle%vol)
    call inout_write_comment(file, "end aero_particle")
    
  end subroutine inout_write_aero_particle
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_read_aero_particle(file, aero_particle)
    
    ! Read full state.
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_t), intent(out) :: aero_particle ! aero_particle to read

    call inout_check_comment(file, "begin aero_particle")
    call inout_read_integer(file, "n_orig_part", aero_particle%n_orig_part)
    call inout_read_real(file, "absorb(m^2)", aero_particle%absorb_cross_sect)
    call inout_read_real(file, "scatter(m^2)", &
         aero_particle%scatter_cross_sect)
    call inout_read_real(file, "asymmetry(1)", aero_particle%asymmetry)
    call inout_read_complex(file, "refract_shell(1)", &
         aero_particle%refract_shell)
    call inout_read_complex(file, "refract_core(1)", aero_particle%refract_core)
    call inout_read_real(file, "core_vol(m^3)", aero_particle%core_vol)
    call inout_read_real_array(file, "spec_vols(m^3)", aero_particle%vol)
    call inout_check_comment(file, "end aero_particle")
    
  end subroutine inout_read_aero_particle
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_aero_particle(val)

    ! Determines the number of bytes required to pack the given value.

    type(aero_particle_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_size_aero_particle = &
         pmc_mpi_pack_size_real_array(val%vol) &
         + pmc_mpi_pack_size_integer(val%n_orig_part) &
         + pmc_mpi_pack_size_real(val%absorb_cross_sect) &
         + pmc_mpi_pack_size_real(val%scatter_cross_sect) &
         + pmc_mpi_pack_size_real(val%asymmetry) &
         + pmc_mpi_pack_size_complex(val%refract_shell) &
         + pmc_mpi_pack_size_complex(val%refract_core) &
         + pmc_mpi_pack_size_real(val%core_vol)
    
  end function pmc_mpi_pack_size_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_aero_particle(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_particle_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%vol)
    call pmc_mpi_pack_integer(buffer, position, val%n_orig_part)
    call pmc_mpi_pack_real(buffer, position, val%absorb_cross_sect)
    call pmc_mpi_pack_real(buffer, position, val%scatter_cross_sect)
    call pmc_mpi_pack_real(buffer, position, val%asymmetry)
    call pmc_mpi_pack_complex(buffer, position, val%refract_shell)
    call pmc_mpi_pack_complex(buffer, position, val%refract_core)
    call pmc_mpi_pack_real(buffer, position, val%core_vol)
    call assert(810223998, position - prev_position &
         == pmc_mpi_pack_size_aero_particle(val))
#endif

  end subroutine pmc_mpi_pack_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_aero_particle(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(aero_particle_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%vol)
    call pmc_mpi_unpack_integer(buffer, position, val%n_orig_part)
    call pmc_mpi_unpack_real(buffer, position, val%absorb_cross_sect)
    call pmc_mpi_unpack_real(buffer, position, val%scatter_cross_sect)
    call pmc_mpi_unpack_real(buffer, position, val%asymmetry)
    call pmc_mpi_unpack_complex(buffer, position, val%refract_shell)
    call pmc_mpi_unpack_complex(buffer, position, val%refract_core)
    call pmc_mpi_unpack_real(buffer, position, val%core_vol)
    call assert(287447241, position - prev_position &
         == pmc_mpi_pack_size_aero_particle(val))
#endif

  end subroutine pmc_mpi_unpack_aero_particle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_aero_particle
