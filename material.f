! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

module mod_material

  type material
     integer n_spec
     integer i_water ! water species number
     real*8, dimension(:), pointer ::  rho ! densities (kg m^{-3})
     integer, dimension(:), pointer :: nu ! number of ions in the solute
     real*8, dimension(:), pointer :: eps ! solubilities (1)
     real*8, dimension(:), pointer :: M_w ! molecular weights (kg mole^{-1})
  end type material

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_material(mat, n_spec)

    type(material), intent(inout) :: mat   ! material properties
    integer, intent(in) :: n_spec          ! number of species

    mat%n_spec = n_spec
    allocate(mat%rho(n_spec))
    allocate(mat%nu(n_spec))
    allocate(mat%eps(n_spec))
    allocate(mat%M_w(n_spec))

  end subroutine allocate_material

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function particle_mass(V, mat) ! kg

    ! total mass of the particle

    real*8, dimension(:), intent(in) :: V         ! species volumes (m^3)
    type(material), intent(in) :: mat             ! material properties
    
    real*8 pm
    integer i

    pm = 0d0
    do i = 1,mat%n_spec
       pm = pm + V(i) * mat%rho(i)
    end do
    particle_mass = pm

  end function particle_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function particle_volume(V, mat) ! m^3

    ! total volume of the particle

    real*8, dimension(:), intent(in) :: V         ! species volumes (m^3)
    type(material), intent(in) :: mat             ! material properties
    
    real*8 pv
    integer i

    pv = 0d0
    do i = 1,mat%n_spec
       pv = pv + V(i)
    end do
    particle_volume = pv

  end function particle_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_solute_quantity(V, mat, quantity)

    ! returns the volume-average of the non-water elements of quantity

    real*8, dimension(:), intent(in) :: V         ! species volumes (m^3)
    type(material), intent(in) :: mat             ! material properties
    real*8, dimension(:), intent(in) :: quantity  ! quantity to average

    real*8 :: ones(mat%n_spec)

    ones = 1d0
    average_solute_quantity = total_solute_quantity(V, mat, quantity) &
         / total_solute_quantity(V, mat, ones)

  end function average_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_solute_quantity(V, mat, quantity)

    ! returns the volume-total of the non-water elements of quantity

    real*8, dimension(:), intent(in) :: V         ! species volumes (m^3)
    type(material), intent(in) :: mat             ! material properties
    real*8, dimension(:), intent(in) :: quantity  ! quantity to total

    real*8 total
    integer i

    total = 0d0
    do i = 1,mat%n_spec
       if (i .ne. mat%i_water) then
          total = total + V(i) * quantity(i)
       end if
    end do
    total_solute_quantity = total

  end function total_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function average_water_quantity(V, mat, quantity)

    ! returns the water element of quantity

    real*8, dimension(:), intent(in) :: V         ! species volumes (m^3)
    type(material), intent(in) :: mat             ! material properties
    real*8, dimension(:), intent(in) :: quantity  ! quantity to average

    average_water_quantity = quantity(mat%i_water)

  end function average_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_water_quantity(V, mat, quantity)

    ! returns the volume-total of the water element of quantity

    real*8, dimension(:), intent(in) :: V         ! species volumes (m^3)
    type(material), intent(in) :: mat             ! material properties
    real*8, dimension(:), intent(in) :: quantity  ! quantity to total

    total_water_quantity = V(mat%i_water) * quantity(mat%i_water)

  end function total_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_material
