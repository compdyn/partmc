! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Basic particle structure.

module mod_aero_particle

  type aero_particle_t
     real*8, pointer :: vol(:)           ! constituent species volumes (m^3)
     integer :: n_orig_part              ! number of original particles
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

    use mod_util

    type(aero_particle_t), intent(in) :: aero_particle_from ! reference particle
    type(aero_particle_t), intent(inout) :: aero_particle_to ! already allocated
    
    integer :: n_spec

    n_spec = size(aero_particle_from%vol)
    if (n_spec /= size(aero_particle_to%vol)) then
       call aero_particle_free(aero_particle_to)
       call aero_particle_alloc(aero_particle_to, n_spec)
    end if
    call assert(size(aero_particle_from%vol) == size(aero_particle_to%vol))
    aero_particle_to%vol = aero_particle_from%vol
    aero_particle_to%n_orig_part = aero_particle_from%n_orig_part

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
    
  end subroutine aero_particle_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_zero(aero_particle)

    ! Resets an aero_particle to be zero.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle to zero
    
    aero_particle%vol = 0d0
    aero_particle%n_orig_part = 1

  end subroutine aero_particle_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_set_vols(aero_particle, vols)

    ! Sets the aerosol particle volumes.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle
    real*8, intent(in) :: vols(size(aero_particle%vol)) ! new volumes

    aero_particle%vol = vols

  end subroutine aero_particle_set_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_mass(aero_particle, aero_data) ! kg

    ! Total mass of the particle.

    use mod_aero_data

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    
    aero_particle_mass = sum(aero_particle%vol * aero_data%density)

  end function aero_particle_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_volume(aero_particle) ! m^3

    ! Total volume of the particle.

    type(aero_particle_t), intent(in) :: aero_particle ! particle

    aero_particle_volume = sum(aero_particle%vol)

  end function aero_particle_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer function aero_particle_in_bin(aero_particle, bin_grid)
    
    ! Find the bin number that contains a given particle.

    use mod_bin_grid
    
    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    
    aero_particle_in_bin = &
         particle_in_bin(aero_particle_volume(aero_particle), bin_grid)
    
  end function aero_particle_in_bin
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_solute_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the volume-average of the non-water elements of quantity.

    use mod_aero_data

    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to average

    real*8 :: ones(aero_data%n_spec)

    ones = 1d0
    average_solute_quantity = &
         total_solute_quantity(aero_particle, aero_data, quantity) &
         / total_solute_quantity(aero_particle, aero_data, ones)

  end function average_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_solute_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the volume-total of the non-water elements of quantity.

    use mod_aero_data

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
    total_solute_quantity = total

  end function total_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_water_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the water element of quantity.

    use mod_util
    use mod_aero_data
    
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to average

    call assert(aero_data%i_water > 0)
    average_water_quantity = quantity(aero_data%i_water)

  end function average_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_water_quantity(aero_particle, &
       aero_data, quantity)

    ! Returns the volume-total of the water element of quantity.

    use mod_util
    use mod_aero_data
    
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to total

    call assert(aero_data%i_water > 0)
    total_water_quantity = aero_particle%vol(aero_data%i_water) &
         * quantity(aero_data%i_water)

  end function total_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_molec_weight(aero_data) ! (kg/mole)

    ! Returns the water molecular weight.

    use mod_util
    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data

    call assert(aero_data%i_water > 0)
    aero_particle_water_molec_weight = aero_data%molec_weight(aero_data%i_water)

  end function aero_particle_water_molec_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_molec_weight(aero_data, &
       aero_particle) ! (kg/mole)

    ! Returns the average of the solute molecular weight.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_molec_weight = average_solute_quantity(aero_particle, &
         aero_data, aero_data%molec_weight)

  end function aero_particle_solute_molec_weight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_num_ions(aero_data, aero_particle) ! (1)

    ! Returns the average of the solute ion number.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_num_ions = average_solute_quantity(aero_particle, &
         aero_data, dble(aero_data%num_ions))

  end function aero_particle_solute_num_ions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_solubility(aero_data, &
       aero_particle) ! (1)

    ! Returns the average of the solute solubilities.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_solubility = average_solute_quantity(aero_particle, &
         aero_data, aero_data%solubility)

  end function aero_particle_solute_solubility

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_density(aero_data) ! (kg/m^3)

    ! Returns the water density.

    use mod_util
    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data

    call assert(aero_data%i_water > 0)
    aero_particle_water_density = aero_data%density(aero_data%i_water)

  end function aero_particle_water_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_density(aero_data, &
       aero_particle) ! (kg/m^3)

    ! Returns the average of the solute densities.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_density = average_solute_quantity(aero_particle, &
         aero_data, aero_data%density)

  end function aero_particle_solute_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_mass(aero_data, aero_particle) ! (kg)

    ! Returns the water mass.

    use mod_util
    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    call assert(aero_data%i_water > 0)
    aero_particle_water_mass = aero_particle%vol(aero_data%i_water) &
         * aero_data%density(aero_data%i_water)

  end function aero_particle_water_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_mass(aero_data, aero_particle) ! (kg)

    ! Returns the total solute mass.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_mass = total_solute_quantity(aero_particle, &
         aero_data, aero_data%density)

  end function aero_particle_solute_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_coagulate(aero_particle_1, &
       aero_particle_2, aero_particle_new)

    ! Coagulate two particles together to make a new one.

    use mod_util

    type(aero_particle_t), intent(in) :: aero_particle_1 ! first particle
    type(aero_particle_t), intent(in) :: aero_particle_2 ! second particle
    type(aero_particle_t), intent(inout) :: aero_particle_new ! result particle

    call assert(size(aero_particle_1%vol) == size(aero_particle_new%vol))
    call assert(size(aero_particle_2%vol) == size(aero_particle_new%vol))
    aero_particle_new%vol = aero_particle_1%vol + aero_particle_2%vol
    aero_particle_new%n_orig_part = aero_particle_1%n_orig_part &
         + aero_particle_2%n_orig_part

  end subroutine aero_particle_coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_write_aero_particle(file, aero_particle)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_t), intent(in) :: aero_particle ! aero_particle to write
    
    call inout_write_integer(file, "n_orig_part", aero_particle%n_orig_part)
    call inout_write_real_array(file, "spec_vols(m^3)", aero_particle%vol)
    
  end subroutine inout_write_aero_particle
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_read_aero_particle(file, aero_particle)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_t), intent(out) :: aero_particle ! aero_particle to read

    call inout_read_integer(file, "n_orig_part", aero_particle%n_orig_part)
    call inout_read_real_array(file, "spec_vols(m^3)", aero_particle%vol)
    
  end subroutine inout_read_aero_particle
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_aero_particle
