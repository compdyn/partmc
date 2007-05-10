! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Gas parameters.

module mod_gas

  type gas_chem
     integer :: n_spec                   ! number of species
     real*8, pointer :: conc(:)          ! length n_spec, concentration (ppb)
     character(len=10), pointer :: name(:) ! length n_spec, name of species
     integer, pointer :: mosaic_index(:) ! length n_spec, to_mosaic(i) is the
                                         ! mosaic index of species i, or 0 if
                                         ! there is no match
  end type gas_chem

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_gas(gas, n_spec)

    ! Allocate storage for gas species.

    type(gas_chem), intent(inout) :: gas  ! gas data
    integer, intent(in) :: n_spec         ! number of species

    gas%n_spec = n_spec
    allocate(gas%conc(n_spec))
    allocate(gas%name(n_spec))
    allocate(gas%mosaic_index(n_spec))

  end subroutine allocate_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function gas_spec_by_name(gas, name)

    ! Returns the number of the species in gas with the given name, or
    ! returns 0 if there is no such species.

    type(gas_chem), intent(in) :: gas     ! gas data
    character*10, intent(in) :: name      ! name of species to find

    integer i
    logical found

    found = .false.
    do i = 1,gas%n_spec
       if (index(name, gas%name(i)) == 1) then
          found = .true.
          exit
       end if
    end do
    if (found) then
       gas_spec_by_name = i
    else
       gas_spec_by_name = 0
    end if

  end function gas_spec_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_gas_mosaic_map(gas)

    ! Fills in gas%mosaic_index.

    type(gas_chem), intent(inout) :: gas  ! gas data

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

    gas%mosaic_index = 0
    do spec = 1,gas%n_spec
       mosaic_spec = 0
       do i = 1,n_mosaic_species
          if (gas%name(spec) == mosaic_species(i)) then
             mosaic_spec = i
          end if
       end do
       if (mosaic_spec > 0) then
          gas%mosaic_index(spec) = mosaic_spec
       end if
    end do

  end subroutine set_gas_mosaic_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_gas
