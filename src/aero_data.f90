! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Aerosol data 

module mod_aero_data

  type aero_data_t
     integer :: n_spec                  ! number of species
     integer :: i_water                 ! water species number
     character(len=10), pointer :: name(:) ! len n_spec, name of species
     integer, pointer :: mosaic_index(:) ! length n_spec, to_mosaic(i) is the
                                        ! mosaic index of species i, or 0 if
                                        ! there is no match
     real*8, pointer ::  rho(:)         ! len n_spec, densities (kg m^{-3})
     integer, pointer :: nu(:)          ! len n_spec, num ions in solute
     real*8, pointer :: eps(:)          ! len n_spec, solubilities (1)
     real*8, pointer :: M_w(:)          ! len n_spec, molec wghts (kg mole^{-1})
  end type aero_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_aero_data(aero_data, n_spec)

    ! Allocate storage for aero_data parameters given the number of
    ! species.

    type(aero_data_t), intent(inout) :: aero_data ! aerosol data
    integer, intent(in) :: n_spec       ! number of species

    aero_data%n_spec = n_spec
    allocate(aero_data%name(n_spec))
    allocate(aero_data%mosaic_index(n_spec))
    allocate(aero_data%rho(n_spec))
    allocate(aero_data%nu(n_spec))
    allocate(aero_data%eps(n_spec))
    allocate(aero_data%M_w(n_spec))
    aero_data%i_water = 0

  end subroutine allocate_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function aero_data_spec_by_name(aero_data, name)

    ! Returns the number of the species in aero_data with the given name, or
    ! returns 0 if there is no such species.

    type(aero_data_t), intent(in) :: aero_data     ! aero_data data
    character*10, intent(in) :: name      ! name of species to find

    integer i
    logical found

    found = .false.
    do i = 1,aero_data%n_spec
       if (index(name, aero_data%name(i)) == 1) then
          found = .true.
          exit
       end if
    end do
    if (found) then
       aero_data_spec_by_name = i
    else
       aero_data_spec_by_name = 0
    end if

  end function aero_data_spec_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_aero_data_water_index(aero_data)

    ! Fills in aero_data%i_water.

    type(aero_data_t), intent(inout) :: aero_data  ! aero_data data

    integer :: i

    do i = 1,aero_data%n_spec
       if (aero_data%name(i) == "H2O") then
          aero_data%i_water = i
       end if
    end do

  end subroutine set_aero_data_water_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_aero_data_mosaic_map(aero_data)

    ! Fills in aero_data%mosaic_index.

    type(aero_data_t), intent(inout) :: aero_data  ! aero_data data

    integer, parameter :: n_mosaic_species = 19
    character*10, parameter, dimension(n_mosaic_species) :: mosaic_species = [ &
         "SO4_a", "NO3_a", "Cl_a", "NH4_a", "CO3_a", "MSA_a", "Na_a", "Ca_a", &
         "OC_a", "BC_a", "OIN_a", "ARO1_a", "ARO2_a", "ALK1_a", "OLE1_a", "API1_a", &
         "API2_a", "LIM1_a", "LIM2_a"]

    integer spec, mosaic_spec, i

    aero_data%mosaic_index = 0
    do spec = 1,aero_data%n_spec
       mosaic_spec = 0
       do i = 1,n_mosaic_species
          if (aero_data%name(spec) == mosaic_species(i)) then
             mosaic_spec = i
          end if
       end do
       if (mosaic_spec > 0) then
          aero_data%mosaic_index(spec) = mosaic_spec
       end if
    end do

  end subroutine set_aero_data_mosaic_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function particle_mass(V, aero_data) ! kg

    ! Total mass of the particle.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    
    real*8 pm
    integer i

    pm = 0d0
    do i = 1,aero_data%n_spec
       pm = pm + V(i) * aero_data%rho(i)
    end do
    particle_mass = pm

  end function particle_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function particle_volume(V) ! m^3

    ! Total volume of the particle.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    
    real*8 pv
    integer i, n_spec

    n_spec = size(V)
    pv = 0d0
    do i = 1,n_spec
       pv = pv + V(i)
    end do
    particle_volume = pv

  end function particle_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_solute_quantity(V, aero_data, quantity)

    ! Returns the volume-average of the non-water elements of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to average

    real*8 :: ones(aero_data%n_spec)

    ones = 1d0
    average_solute_quantity = total_solute_quantity(V, aero_data, quantity) &
         / total_solute_quantity(V, aero_data, ones)

  end function average_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_solute_quantity(V, aero_data, quantity)

    ! Returns the volume-total of the non-water elements of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to total

    real*8 total
    integer i

    total = 0d0
    do i = 1,aero_data%n_spec
       if (i .ne. aero_data%i_water) then
          total = total + V(i) * quantity(i)
       end if
    end do
    total_solute_quantity = total

  end function total_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_water_quantity(V, aero_data, quantity)

    ! Returns the water element of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to average

    average_water_quantity = quantity(aero_data%i_water)

  end function average_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_water_quantity(V, aero_data, quantity)

    ! Returns the volume-total of the water element of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    real*8, intent(in) :: quantity(:)   ! quantity to total

    total_water_quantity = V(aero_data%i_water) * quantity(aero_data%i_water)

  end function total_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_aero_data
