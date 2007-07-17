! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Basic particle structure.

module mod_aero_particle

  type aero_particle_t
     real*8, pointer :: vols(:)          ! constituent species volumes
  end type aero_particle_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_alloc(n_spec, aero_particle)

    ! Allocates and initializes.

    integer, intent(in) :: n_spec       ! number of species
    type(aero_particle_t), intent(inout) :: aero_particle ! particle to init

    allocate(aero_particle%vols(n_spec))
    call aero_particle_zero(aero_particle)

  end subroutine aero_particle_alloc
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_free(aero_particle)

    ! Deallocates.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle to free
    
    deallocate(aero_particle%vols)

  end subroutine aero_particle_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_copy(aero_particle_from, aero_particle_to)

    ! Copies a particle.

    use mod_util

    type(aero_particle_t), intent(in) :: aero_particle_from ! reference particle
    type(aero_particle_t), intent(inout) :: aero_particle_to ! already allocated
    
    integer :: n_spec

    n_spec = size(aero_particle_from%vols)
    if (n_spec /= size(aero_particle_to%vols)) then
       call aero_particle_free(aero_particle_to)
       call aero_particle_alloc(n_spec, aero_particle_to)
    end if
    call assert(size(aero_particle_from%vols) == size(aero_particle_to%vols))
    aero_particle_to%vols = aero_particle_from%vols

  end subroutine aero_particle_copy
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_shift(aero_particle_from, aero_particle_to)

    ! Is the equivalent of aero_particle_copy(from, to) followed by
    ! aero_particle_free(from), but faster.

    type(aero_particle_t), intent(in) :: aero_particle_from ! reference particle
    type(aero_particle_t), intent(inout) :: aero_particle_to ! not allocated

    aero_particle_to%vols => aero_particle_from%vols
    nullify(aero_particle_from%vols)

  end subroutine aero_particle_shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_zero(aero_particle)

    ! Resets an aero_particle to be zero.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle to zero
    
    aero_particle%vols = 0d0

  end subroutine aero_particle_zero
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_set_vols(aero_particle, vols)

    ! Sets the aerosol particle volumes.

    type(aero_particle_t), intent(inout) :: aero_particle ! particle
    real*8, intent(in) :: vols(size(aero_particle%vols)) ! new volumes

    aero_particle%vols = vols

  end subroutine aero_particle_set_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_mass(aero_particle, aero_data) ! kg

    ! Total mass of the particle.

    use mod_aero_data

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    
    aero_particle_mass = sum(aero_particle%vols * aero_data%rho)

  end function aero_particle_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_volume(aero_particle) ! m^3

    ! Total volume of the particle.

    type(aero_particle_t), intent(in) :: aero_particle ! particle

    aero_particle_volume = sum(aero_particle%vols)

  end function aero_particle_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine aero_particle_in_bin(aero_particle, bin_grid, bin)
    
    ! Find the bin number that contains a given particle.

    use mod_bin_grid
    
    type(aero_particle_t), intent(in) :: aero_particle ! particle
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid
    integer, intent(out) :: bin         ! bin number containing particle
    
    call particle_in_bin(aero_particle_volume(aero_particle), bin_grid, bin)
    
  end subroutine aero_particle_in_bin
  
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
          total = total + aero_particle%vols(i) * quantity(i)
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
    total_water_quantity = aero_particle%vols(aero_data%i_water) &
         * quantity(aero_data%i_water)

  end function total_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_M_w(aero_data) ! (kg/mole)

    ! Returns the water molecular weight.

    use mod_util
    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data

    call assert(aero_data%i_water > 0)
    aero_particle_water_M_w = aero_data%M_w(aero_data%i_water)

  end function aero_particle_water_M_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_M_w(aero_data, aero_particle) ! (kg/mole)

    ! Returns the average of the solute molecular weight.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_M_w = average_solute_quantity(aero_particle, &
         aero_data, aero_data%M_w)

  end function aero_particle_solute_M_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_nu(aero_data, aero_particle) ! (1)

    ! Returns the average of the solute ion number.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_nu = average_solute_quantity(aero_particle, &
         aero_data, dble(aero_data%nu))

  end function aero_particle_solute_nu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_eps(aero_data, aero_particle) ! (1)

    ! Returns the average of the solute solubilities.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_eps = average_solute_quantity(aero_particle, &
         aero_data, aero_data%eps)

  end function aero_particle_solute_eps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_rho(aero_data) ! (kg/m^3)

    ! Returns the water density.

    use mod_util
    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data

    call assert(aero_data%i_water > 0)
    aero_particle_water_rho = aero_data%rho(aero_data%i_water)

  end function aero_particle_water_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_rho(aero_data, aero_particle) ! (kg/m^3)

    ! Returns the average of the solute densities.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_rho = average_solute_quantity(aero_particle, &
         aero_data, aero_data%rho)

  end function aero_particle_solute_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_water_mass(aero_data, aero_particle) ! (kg)

    ! Returns the water mass.

    use mod_util
    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    call assert(aero_data%i_water > 0)
    aero_particle_water_mass = aero_particle%vols(aero_data%i_water) &
         * aero_data%rho(aero_data%i_water)

  end function aero_particle_water_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_solute_mass(aero_data, aero_particle) ! (kg)

    ! Returns the total solute mass.

    use mod_aero_data

    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_particle_t), intent(in) :: aero_particle ! aerosol particle

    aero_particle_solute_mass = total_solute_quantity(aero_particle, &
         aero_data, aero_data%rho)

  end function aero_particle_solute_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_particle_coagulate(aero_particle_1, &
       aero_particle_2, aero_particle_new)

    ! Coagulate two particles together to make a new one.

    use mod_util

    type(aero_particle_t), intent(in) :: aero_particle_1 ! first particle
    type(aero_particle_t), intent(in) :: aero_particle_2 ! second particle
    type(aero_particle_t), intent(inout) :: aero_particle_new ! result particle

    call assert(size(aero_particle_1%vols) == size(aero_particle_new%vols))
    call assert(size(aero_particle_2%vols) == size(aero_particle_new%vols))
    aero_particle_new%vols = aero_particle_1%vols + aero_particle_2%vols

  end subroutine aero_particle_coagulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_write_aero_particle(file, aero_particle)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_t), intent(in) :: aero_particle ! aero_particle to write
    
    call inout_write_real_array(file, "spec_vols(m^3)", aero_particle%vols)
    
  end subroutine inout_write_aero_particle
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine inout_read_aero_particle(file, aero_particle)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_particle_t), intent(out) :: aero_particle ! aero_particle to read

    call inout_read_real_array(file, "spec_vols(m^3)", aero_particle%vols)
    
  end subroutine inout_read_aero_particle
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_aero_particle
