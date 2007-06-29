! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Gas parameters.

module mod_gas_data

  type gas_data_t
     integer :: n_spec                   ! number of species
     real*8, pointer :: M_w(:)           ! molecular weight (kg mole^{-1})
     character(len=10), pointer :: name(:) ! length n_spec, name of species
     integer, pointer :: mosaic_index(:) ! length n_spec, to_mosaic(i) is the
                                         ! mosaic index of species i, or 0 if
                                         ! there is no match
  end type gas_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_gas_data(n_spec, gas_data)

    ! Allocate storage for gas species.

    integer, intent(in) :: n_spec         ! number of species
    type(gas_data_t), intent(out) :: gas_data ! gas data

    gas_data%n_spec = n_spec
    allocate(gas_data%M_w(n_spec))
    allocate(gas_data%name(n_spec))
    allocate(gas_data%mosaic_index(n_spec))

  end subroutine allocate_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function gas_spec_by_name(gas_data, name)

    ! Returns the number of the species in gas with the given name, or
    ! returns 0 if there is no such species.

    type(gas_data_t), intent(in) :: gas_data ! gas data
    character*10, intent(in) :: name      ! name of species to find

    integer i
    logical found

    found = .false.
    do i = 1,gas_data%n_spec
       if (index(name, gas_data%name(i)) == 1) then
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

  subroutine set_gas_mosaic_map(gas_data)

    ! Fills in gas_data%mosaic_index.

    type(gas_data_t), intent(inout) :: gas_data ! gas data

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

    gas_data%mosaic_index = 0
    do spec = 1,gas_data%n_spec
       mosaic_spec = 0
       do i = 1,n_mosaic_species
          if (gas_data%name(spec) == mosaic_species(i)) then
             mosaic_spec = i
          end if
       end do
       if (mosaic_spec > 0) then
          gas_data%mosaic_index(spec) = mosaic_spec
       end if
    end do

  end subroutine set_gas_mosaic_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_gas_data
