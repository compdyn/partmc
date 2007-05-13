! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Material parameters.

module mod_material

  type material
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
  end type material

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_material(mat, n_spec)

    ! Allocate storage for material parameters given the number of
    ! species.

    type(material), intent(inout) :: mat ! material properties
    integer, intent(in) :: n_spec       ! number of species

    mat%n_spec = n_spec
    allocate(mat%name(n_spec))
    allocate(mat%mosaic_index(n_spec))
    allocate(mat%rho(n_spec))
    allocate(mat%nu(n_spec))
    allocate(mat%eps(n_spec))
    allocate(mat%M_w(n_spec))
    mat%i_water = 0

  end subroutine allocate_material

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function material_spec_by_name(mat, name)

    ! Returns the number of the species in mat with the given name, or
    ! returns 0 if there is no such species.

    type(material), intent(in) :: mat     ! material data
    character*10, intent(in) :: name      ! name of species to find

    integer i
    logical found

    found = .false.
    do i = 1,mat%n_spec
       if (index(name, mat%name(i)) == 1) then
          found = .true.
          exit
       end if
    end do
    if (found) then
       material_spec_by_name = i
    else
       material_spec_by_name = 0
    end if

  end function material_spec_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_material_water_index(mat)

    ! Fills in mat%i_water.

    type(material), intent(inout) :: mat  ! material data

    integer :: i

    do i = 1,mat%n_spec
       if (mat%name(i) == "H2O") then
          mat%i_water = i
       end if
    end do

  end subroutine set_material_water_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_material_mosaic_map(mat)

    ! Fills in mat%mosaic_index.

    type(material), intent(inout) :: mat  ! material data

    ! FIXME: this list is wrong (it is the gas phase list, but should
    ! be the aerosol list)
    integer, parameter :: n_mosaic_species = 77
    character*10, parameter, dimension(n_mosaic_species) :: mosaic_species = [ &
         "H2SO4", "HNO3", "HCl", "NH3", "NO", "NO2", "NO3", "N2O5", &
         "HONO", "HNO4", "O3", "O1D", "O3P", "OH", "HO2", "H2O2", &
         "CO", "SO2", "CH4", "C2H6", "CH3O2", "ETHP", "HCHO", "CH3OH", &
         "ANOL", "CH3OOH", "ETHOOH", "ALD2", "HCOOH", "RCOOH", "C2O3", &
         "PAN", "ARO1", "ARO2", "ALK1", "OLE1", "API1", "API2", &
         "LIM1", "LIM2", "PAR", "AONE", "MGLY", "ETH", "OLET", "OLEI", &
         "TOL", "XYL", "CRES", "TO2", "CRO", "OPEN", "ONIT", "ROOH", &
         "RO2", "ANO2", "NAP", "XO2", "XPAR", "ISOP", "ISOPRD", &
         "ISOPP", "ISOPN", "ISOPO2", "API", "LIM", "DMS", "MSA", &
         "DMSO", "DMSO2", "CH3SO2H", "CH3SCH2OO", "CH3SO2", "CH3SO3", &
         "CH3SO2OO", "CH3SO2CH2OO", "SULFHOX"]

    integer spec, mosaic_spec, i

    mat%mosaic_index = 0
    do spec = 1,mat%n_spec
       mosaic_spec = 0
       do i = 1,n_mosaic_species
          if (mat%name(spec) == mosaic_species(i)) then
             mosaic_spec = i
          end if
       end do
       if (mosaic_spec > 0) then
          mat%mosaic_index(spec) = mosaic_spec
       end if
    end do

  end subroutine set_material_mosaic_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function particle_mass(V, mat) ! kg

    ! Total mass of the particle.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(material), intent(in) :: mat   ! material properties
    
    real*8 pm
    integer i

    pm = 0d0
    do i = 1,mat%n_spec
       pm = pm + V(i) * mat%rho(i)
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

  real*8 function average_solute_quantity(V, mat, quantity)

    ! Returns the volume-average of the non-water elements of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(material), intent(in) :: mat   ! material properties
    real*8, intent(in) :: quantity(:)   ! quantity to average

    real*8 :: ones(mat%n_spec)

    ones = 1d0
    average_solute_quantity = total_solute_quantity(V, mat, quantity) &
         / total_solute_quantity(V, mat, ones)

  end function average_solute_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_solute_quantity(V, mat, quantity)

    ! Returns the volume-total of the non-water elements of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(material), intent(in) :: mat   ! material properties
    real*8, intent(in) :: quantity(:)   ! quantity to total

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

    ! Returns the water element of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(material), intent(in) :: mat   ! material properties
    real*8, intent(in) :: quantity(:)   ! quantity to average

    average_water_quantity = quantity(mat%i_water)

  end function average_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function total_water_quantity(V, mat, quantity)

    ! Returns the volume-total of the water element of quantity.

    real*8, intent(in) :: V(:)          ! species volumes (m^3)
    type(material), intent(in) :: mat   ! material properties
    real*8, intent(in) :: quantity(:)   ! quantity to total

    total_water_quantity = V(mat%i_water) * quantity(mat%i_water)

  end function total_water_quantity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_material
