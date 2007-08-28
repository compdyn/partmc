! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Gas parameters.

module pmc_gas_data

  integer, parameter :: GAS_NAME_LEN = 15

  type gas_data_t
     integer :: n_spec                   ! number of species
     real*8, pointer :: molec_weight(:)  ! molecular weight (kg mole^{-1})
     character(len=GAS_NAME_LEN), pointer :: name(:) ! len n_spec, species name
     integer, pointer :: mosaic_index(:) ! length n_spec, to_mosaic(i) is the
                                         ! mosaic index of species i, or 0 if
                                         ! there is no match
  end type gas_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_data_alloc(gas_data, n_spec)

    ! Allocate storage for gas species.

    type(gas_data_t), intent(out) :: gas_data ! gas data
    integer, intent(in) :: n_spec         ! number of species

    gas_data%n_spec = n_spec
    allocate(gas_data%molec_weight(n_spec))
    allocate(gas_data%name(n_spec))
    allocate(gas_data%mosaic_index(n_spec))

  end subroutine gas_data_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_data_free(gas_data)

    ! Free all storage.

    type(gas_data_t), intent(out) :: gas_data ! gas data

    deallocate(gas_data%molec_weight)
    deallocate(gas_data%name)
    deallocate(gas_data%mosaic_index)

  end subroutine gas_data_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function gas_data_spec_by_name(gas_data, name)

    ! Returns the number of the species in gas with the given name, or
    ! returns 0 if there is no such species.

    type(gas_data_t), intent(in) :: gas_data ! gas data
    character(len=GAS_NAME_LEN), intent(in) :: name ! name of species to find

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
       gas_data_spec_by_name = i
    else
       gas_data_spec_by_name = 0
    end if

  end function gas_data_spec_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gas_data_set_mosaic_map(gas_data)

    ! Fills in gas_data%mosaic_index.

    type(gas_data_t), intent(inout) :: gas_data ! gas data

    integer, parameter :: n_mosaic_species = 77
    character(len=GAS_NAME_LEN), parameter, dimension(n_mosaic_species) &
         :: mosaic_species = [ &
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

  end subroutine gas_data_set_mosaic_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_write_gas_data(file, gas_data)
    
    ! Write full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(gas_data_t), intent(in) :: gas_data ! gas_data to write

    call inout_write_comment(file, "begin gas_data")
    call inout_write_integer(file, "n_spec", gas_data%n_spec)
    call inout_write_real_array(file, "molec_wght(kg/mole)", &
         gas_data%molec_weight)
    call inout_write_string_array(file, "species_names", gas_data%name)
    call inout_write_integer_array(file, "mosaic_indices", &
         gas_data%mosaic_index)
    call inout_write_comment(file, "end gas_data")
    
  end subroutine inout_write_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_gas_data(file, gas_data)
    
    ! Read full state.
    
    use pmc_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(gas_data_t), intent(out) :: gas_data ! gas_data to read

    call inout_check_comment(file, "begin gas_data")
    call inout_read_integer(file, "n_spec", gas_data%n_spec)
    call inout_read_real_array(file, "molec_wght(kg/mole)", &
         gas_data%molec_weight)
    call inout_read_string_array(file, "species_names", gas_data%name)
    call inout_read_integer_array(file, "mosaic_indices", &
         gas_data%mosaic_index)
    call inout_check_comment(file, "end gas_data")
    
  end subroutine inout_read_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_gas_data(file, gas_data)

    ! Read gas data from a .spec file.

    use pmc_inout

    type(inout_file_t), intent(inout) :: file ! spec file
    type(gas_data_t), intent(out) :: gas_data ! gas data

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)

    ! read the gas data from the specified file
    call inout_read_string(file, 'gas_data', read_name)
    call inout_open_read(read_name, read_file)
    call inout_read_real_named_array(read_file, 0, species_name, species_data)
    call inout_close(read_file)

    ! check the data size
    if (size(species_data, 2) /= 1) then
       write(0,*) 'ERROR: each line in ', trim(read_name), &
            ' should only contain one value'
       call exit(1)
    end if

    ! allocate and copy over the data
    n_species = size(species_data, 1)
    call gas_data_alloc(gas_data, n_species)
    do i = 1,n_species
       gas_data%name(i) = species_name(i)(1:GAS_NAME_LEN)
       gas_data%molec_weight(i) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)
    
    call gas_data_set_mosaic_map(gas_data)

  end subroutine spec_read_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_gas_data(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_mpi

    type(gas_data_t), intent(in) :: val ! value to pack

    pmc_mpi_pack_size_gas_data = &
         pmc_mpi_pack_size_integer(val%n_spec) &
         + pmc_mpi_pack_size_real_array(val%molec_weight) &
         + pmc_mpi_pack_size_string_array(val%name) &
         + pmc_mpi_pack_size_integer_array(val%mosaic_index)

  end function pmc_mpi_pack_size_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_gas_data(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(gas_data_t), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_spec)
    call pmc_mpi_pack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_pack_string_array(buffer, position, val%name)
    call pmc_mpi_pack_integer_array(buffer, position, val%mosaic_index)
    call assert(position - prev_position == pmc_mpi_pack_size_gas_data(val))
#endif

  end subroutine pmc_mpi_pack_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_gas_data(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

#ifdef PMC_USE_MPI
    use mpi
    use pmc_mpi
    use pmc_util
#endif

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    type(gas_data_t), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_spec)
    call pmc_mpi_unpack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_unpack_string_array(buffer, position, val%name)
    call pmc_mpi_unpack_integer_array(buffer, position, val%mosaic_index)
    call assert(position - prev_position == pmc_mpi_pack_size_gas_data(val))
#endif

  end subroutine pmc_mpi_unpack_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_gas_data
