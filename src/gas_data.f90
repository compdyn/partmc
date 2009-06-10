! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_gas_data module.

!> The gas_data_t structure and associated subroutines.
module pmc_gas_data

  use pmc_spec_read
  use pmc_mpi
  use pmc_util
  use pmc_netcdf
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Maximum length of the name of a gas.
  integer, parameter :: GAS_NAME_LEN = 100

  !> Constant gas data.
  !!
  !! Each gas species is identified by an integer \c i between 1 and
  !! \c n_spec. Species \c i has name \c gas_data%%name(i), molecular
  !! weight gas_data%%molec_weight(i), etc. The variable gas data
  !! describing the current concentrations is stored in the
  !! gas_state_t structure, so the concentration of species \c i is
  !! gas_state%%conc(i).
  type gas_data_t
     !> Number of species.
     integer :: n_spec
     !> Molecular weight (kg mole^{-1}).
     real*8, pointer :: molec_weight(:)
     !> Species name [length \c n_spec].
     character(len=GAS_NAME_LEN), pointer :: name(:)
     !> Index of the corresponding MOSAIC species [length \c
     !> n_spec]. \c to_mosaic(i) is the mosaic index of species \c i,
     !> or 0 if there is no match.
     integer, pointer :: mosaic_index(:)
  end type gas_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for gas species.
  subroutine gas_data_alloc(gas_data, n_spec)

    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data
    !> Number of species.
    integer, intent(in) :: n_spec

    gas_data%n_spec = n_spec
    allocate(gas_data%molec_weight(n_spec))
    allocate(gas_data%name(n_spec))
    allocate(gas_data%mosaic_index(n_spec))

  end subroutine gas_data_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine gas_data_free(gas_data)

    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data

    deallocate(gas_data%molec_weight)
    deallocate(gas_data%name)
    deallocate(gas_data%mosaic_index)

  end subroutine gas_data_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of the species in gas with the given name, or
  !> returns 0 if there is no such species.
  integer function gas_data_spec_by_name(gas_data, name)

    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Name of species to find.
    character(len=GAS_NAME_LEN), intent(in) :: name

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

  !> Fills in gas_data%mosaic_index.
  subroutine gas_data_set_mosaic_map(gas_data)

    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data

    integer, parameter :: n_mosaic_species = 77
    character(len=GAS_NAME_LEN), parameter, dimension(n_mosaic_species) &
         :: mosaic_species = [ &
         "H2SO4        ", "HNO3         ", "HCl          ", &
         "NH3          ", "NO           ", "NO2          ", &
         "NO3          ", "N2O5         ", "HONO         ", &
         "HNO4         ", "O3           ", "O1D          ", &
         "O3P          ", "OH           ", "HO2          ", &
         "H2O2         ", "CO           ", "SO2          ", &
         "CH4          ", "C2H6         ", "CH3O2        ", &
         "ETHP         ", "HCHO         ", "CH3OH        ", &
         "ANOL         ", "CH3OOH       ", "ETHOOH       ", &
         "ALD2         ", "HCOOH        ", "RCOOH        ", &
         "C2O3         ", "PAN          ", "ARO1         ", &
         "ARO2         ", "ALK1         ", "OLE1         ", &
         "API1         ", "API2         ", "LIM1         ", &
         "LIM2         ", "PAR          ", "AONE         ", &
         "MGLY         ", "ETH          ", "OLET         ", &
         "OLEI         ", "TOL          ", "XYL          ", &
         "CRES         ", "TO2          ", "CRO          ", &
         "OPEN         ", "ONIT         ", "ROOH         ", &
         "RO2          ", "ANO2         ", "NAP          ", &
         "XO2          ", "XPAR         ", "ISOP         ", &
         "ISOPRD       ", "ISOPP        ", "ISOPN        ", &
         "ISOPO2       ", "API          ", "LIM          ", &
         "DMS          ", "MSA          ", "DMSO         ", &
         "DMSO2        ", "CH3SO2H      ", "CH3SCH2OO    ", &
         "CH3SO2       ", "CH3SO3       ", "CH3SO2OO     ", &
         "CH3SO2CH2OO  ", "SULFHOX      "]

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

  !> Read gas data from a .spec file.
  subroutine spec_read_gas_data(file, gas_data)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data

    integer :: n_species, species, i
    character(len=MAX_VAR_LEN) :: read_name
    type(spec_file_t) :: read_file
    character(len=MAX_VAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)

    ! read the gas data from the specified file
    call spec_read_string(file, 'gas_data', read_name)
    call spec_read_open(read_name, read_file)
    call spec_read_real_named_array(read_file, 0, species_name, species_data)
    call spec_read_close(read_file)

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

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_gas_data(val)

    !> Value to pack.
    type(gas_data_t), intent(in) :: val

    pmc_mpi_pack_size_gas_data = &
         pmc_mpi_pack_size_integer(val%n_spec) &
         + pmc_mpi_pack_size_real_array(val%molec_weight) &
         + pmc_mpi_pack_size_string_array(val%name) &
         + pmc_mpi_pack_size_integer_array(val%mosaic_index)

  end function pmc_mpi_pack_size_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_gas_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(gas_data_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_spec)
    call pmc_mpi_pack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_pack_string_array(buffer, position, val%name)
    call pmc_mpi_pack_integer_array(buffer, position, val%mosaic_index)
    call assert(449872094, &
         position - prev_position == pmc_mpi_pack_size_gas_data(val))
#endif

  end subroutine pmc_mpi_pack_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_gas_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(gas_data_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_spec)
    call pmc_mpi_unpack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_unpack_string_array(buffer, position, val%name)
    call pmc_mpi_unpack_integer_array(buffer, position, val%mosaic_index)
    call assert(137879163, &
         position - prev_position == pmc_mpi_pack_size_gas_data(val))
#endif

  end subroutine pmc_mpi_unpack_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the gas species dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine gas_data_netcdf_dim_gas_species(gas_data, ncid, &
       dimid_gas_species)

    !> Gas_data structure.
    type(gas_data_t), intent(in) :: gas_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_gas_species

    integer :: status, i_spec
    integer :: varid_gas_species
    integer :: gas_species_centers(gas_data%n_spec)
    character(len=((GAS_NAME_LEN + 2) * gas_data%n_spec)) :: gas_species_names

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "gas_species", dimid_gas_species)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "gas_species", &
         gas_data%n_spec, dimid_gas_species))
    gas_species_names = ""
    do i_spec = 1,gas_data%n_spec
       gas_species_names((len_trim(gas_species_names) + 1):) &
            = trim(gas_data%name(i_spec))
       if (i_spec < gas_data%n_spec) then
          gas_species_names((len_trim(gas_species_names) + 1):) = ","
       end if
    end do
    call pmc_nc_check(nf90_def_var(ncid, "gas_species", NF90_INT, &
         dimid_gas_species, varid_gas_species))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species, "unit", "1"))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species, "names", &
         gas_species_names))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_spec = 1,gas_data%n_spec
       gas_species_centers(i_spec) = i_spec
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_gas_species, &
         gas_species_centers))

  end subroutine gas_data_netcdf_dim_gas_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine gas_data_output_netcdf(gas_data, ncid)
    
    !> Gas_data to write.
    type(gas_data_t), intent(in) :: gas_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_gas_species

    call gas_data_netcdf_dim_gas_species(gas_data, ncid, &
         dimid_gas_species)

    call pmc_nc_write_integer_1d(ncid, gas_data%mosaic_index, &
         "gas_mosaic_index", "1", (/ dimid_gas_species /))
    call pmc_nc_write_real_1d(ncid, gas_data%molec_weight, &
         "gas_molec_weight", "kg/mole", (/ dimid_gas_species /))

  end subroutine gas_data_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine gas_data_input_netcdf(gas_data, ncid)
    
    !> Gas_data to read.
    type(gas_data_t), intent(inout) :: gas_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    character(len=1000) :: unit, name
    integer :: dimid_gas_species, n_spec, varid_gas_species, i_spec, i
    character(len=((GAS_NAME_LEN + 2) * 1000)) :: gas_species_names

    call pmc_nc_check(nf90_inq_dimid(ncid, "gas_species", dimid_gas_species))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_gas_species, name, n_spec))
    call gas_data_free(gas_data)
    call gas_data_alloc(gas_data, n_spec)
    call assert(719237193, n_spec < 1000)

    call pmc_nc_read_integer_1d(ncid, gas_data%mosaic_index, &
         "gas_mosaic_index", unit)
    call pmc_nc_read_real_1d(ncid, gas_data%molec_weight, &
         "gas_molec_weight", unit)

    call pmc_nc_check(nf90_inq_varid(ncid, "gas_species", varid_gas_species))
    call pmc_nc_check(nf90_get_att(ncid, varid_gas_species, "names", gas_species_names))
    ! gas_species_names are comma-separated, so unpack them
    do i_spec = 1,gas_data%n_spec
       i = 1
       do while ((gas_species_names(i:i) /= " ") &
            .and. (gas_species_names(i:i) /= ","))
          i = i + 1
       end do
       call assert(173021381, i > 1)
       gas_data%name(i_spec) = gas_species_names(1:(i-1))
       gas_species_names = gas_species_names((i+1):)
    end do
    call assert(729138192, gas_species_names == "")

  end subroutine gas_data_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_gas_data
