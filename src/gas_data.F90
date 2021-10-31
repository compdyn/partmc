! Copyright (C) 2005-2012, 2021 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_gas_data module.

!> The gas_data_t structure and associated subroutines.
module pmc_gas_data

  use pmc_spec_file
  use pmc_mpi
  use pmc_util
  use pmc_netcdf
#ifdef PMC_USE_CAMP
  use camp_camp_core
  use camp_chem_spec_data
  use camp_property
  use camp_util, only: string_t, split_string, split_char
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Maximum length of the name of a gas.
  integer, parameter :: GAS_NAME_LEN = 100

  !> Constant gas data.
  !!
  !! Each gas species is identified by an integer \c i between 1 and
  !! \c gas_data_n_spec(gas_data). Species \c i has name \c gas_data%%name(i).
  !! The variable gas data describing the current mixing ratios is stored
  !! in the gas_state_t structure, so the mixing ratio of species \c i
  !! is gas_state%%mix_rat(i).
  type gas_data_t
     !> Species name [length \c gas_data_n_spec(gas_data)].
     character(len=GAS_NAME_LEN), allocatable :: name(:)
     !> Index of the corresponding MOSAIC species [length \c
     !> gas_data_n_spec(gas_data)]. \c to_mosaic(i) is the mosaic index of
     !> species \c i, or 0 if there is no match.
     integer, allocatable :: mosaic_index(:)
     !> Index for gas-phase water in the CAMP state array
     integer :: i_camp_water = 0
  end type gas_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PMC_USE_CAMP
  !> Initialize the gas_data_t instance from camp_core data
  subroutine gas_data_initialize(gas_data, camp_core)

    !> Gas-phase species data
    class(gas_data_t), intent(inout) :: gas_data 
    !> CAMP core
    class(camp_core_t), intent(in) :: camp_core

    type(chem_spec_data_t), pointer :: chem_spec_data
    integer :: i_spec
    type(string_t), allocatable :: gas_spec_names(:)
    type(property_t), pointer :: property_set
    character(len=:), allocatable :: prop_name
    logical :: bool_val

    ! Get the chemical species data
    call assert_msg(139566827, &
            camp_core%get_chem_spec_data(chem_spec_data), &
            "No chemical species data in camp core.")

    ! Get the gas-phase species names
    gas_spec_names = chem_spec_data%get_spec_names( &
            spec_phase = CHEM_SPEC_GAS_PHASE)

    ! Allocate space for the gas-phase species
    allocate(gas_data%name(size(gas_spec_names)))

    ! Set the species names and locate gas-phase water
    prop_name = "is gas-phase water"
    do i_spec = 1, size(gas_spec_names)
      gas_data%name(i_spec) = gas_spec_names(i_spec)%string
      call assert_msg(990037352, &
                      chem_spec_data%get_property_set( &
                        gas_spec_names(i_spec)%string, &
                        property_set), &
                      "Missing property set for gas species "// &
                      gas_spec_names(i_spec)%string)
      if (property_set%get_logical(prop_name, bool_val)) then
        call assert_msg(423633615, gas_data%i_camp_water == 0, &
                        "More than one gas-phase water species specified")
        gas_data%i_camp_water = i_spec
      end if
    end do

    call assert_msg(134440820, gas_data%i_camp_water /= 0, &
                    "No gas-phase water species specified.")

    ! Allocate the mosaic index array and set to zero
    allocate(gas_data%mosaic_index(size(gas_spec_names)))
    gas_data%mosaic_index(:) = 0

  end subroutine gas_data_initialize
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the number of gas species.
  elemental integer function gas_data_n_spec(gas_data)

    !> Aero data structure.
    type(gas_data_t), intent(in) :: gas_data

    if (allocated(gas_data%name)) then
       gas_data_n_spec = size(gas_data%name)
    else
       gas_data_n_spec = 0
    end if

  end function gas_data_n_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of the species in gas with the given name, or
  !> returns 0 if there is no such species.
  integer function gas_data_spec_by_name(gas_data, name)

    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Name of species to find.
    character(len=*), intent(in) :: name

    integer i
    logical found

    found = .false.
    do i = 1,gas_data_n_spec(gas_data)
       if (name == gas_data%name(i)) then
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fills in gas_data%%mosaic_index.
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
    do spec = 1,gas_data_n_spec(gas_data)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read gas data from a .spec file.
  subroutine spec_file_read_gas_data(file, gas_data)

    !> Spec file to read data from.
    type(spec_file_t), intent(inout) :: file
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data

    integer :: n_species, species, i
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: species_name(:)
    real(kind=dp), allocatable :: species_data(:,:)

    !> \page input_format_gas_data Input File Format: Gas Material Data
    !!
    !! A gas material data file must consist of one line per gas
    !! species, with each line having the species name. This specifies
    !! which species are to be recognized as gases. For example, a \c
    !! gas_data file could contain:
    !! <pre>
    !! H2SO4
    !! HNO3
    !! HCl
    !! NH3
    !! </pre>
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_gas_data --- the corresponding output format

    ! read the gas data from the specified file
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    if (size(species_data, 2) /= 0) then
       call die_msg(614290516, 'each line in ' // trim(file%name) &
            // ' must only contain the species name')
    end if

    ! allocate and copy over the data
    n_species = size(species_data, 1)
    call ensure_string_array_size(gas_data%name, n_species)
    do i = 1,n_species
       gas_data%name(i) = species_name(i)(1:GAS_NAME_LEN)
    end do

    call ensure_integer_array_size(gas_data%mosaic_index, n_species)
    call gas_data_set_mosaic_map(gas_data)

  end subroutine spec_file_read_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_gas_data(val)

    !> Value to pack.
    type(gas_data_t), intent(in) :: val

    pmc_mpi_pack_size_gas_data = &
         pmc_mpi_pack_size_string_array(val%name) &
         + pmc_mpi_pack_size_integer_array(val%mosaic_index)

  end function pmc_mpi_pack_size_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call pmc_mpi_pack_string_array(buffer, position, val%name)
    call pmc_mpi_pack_integer_array(buffer, position, val%mosaic_index)
    call assert(449872094, &
         position - prev_position <= pmc_mpi_pack_size_gas_data(val))
#endif

  end subroutine pmc_mpi_pack_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_gas_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(gas_data_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_string_array(buffer, position, val%name)
    call pmc_mpi_unpack_integer_array(buffer, position, val%mosaic_index)
    call assert(137879163, &
         position - prev_position <= pmc_mpi_pack_size_gas_data(val))
#endif

  end subroutine pmc_mpi_unpack_gas_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    integer :: gas_species_centers(gas_data_n_spec(gas_data))
    character(len=((GAS_NAME_LEN + 2) * gas_data_n_spec(gas_data))) &
         :: gas_species_names

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "gas_species", dimid_gas_species)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "gas_species", &
         gas_data_n_spec(gas_data), dimid_gas_species))
    gas_species_names = ""
    do i_spec = 1,gas_data_n_spec(gas_data)
       gas_species_names((len_trim(gas_species_names) + 1):) &
            = trim(gas_data%name(i_spec))
       if (i_spec < gas_data_n_spec(gas_data)) then
          gas_species_names((len_trim(gas_species_names) + 1):) = ","
       end if
    end do
    call pmc_nc_check(nf90_def_var(ncid, "gas_species", NF90_INT, &
         dimid_gas_species, varid_gas_species))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species, "names", &
         gas_species_names))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species, "description", &
         "dummy dimension variable (no useful value) - read species names " &
         // "as comma-separated values from the 'names' attribute"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_spec = 1,gas_data_n_spec(gas_data)
       gas_species_centers(i_spec) = i_spec
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_gas_species, &
         gas_species_centers))

  end subroutine gas_data_netcdf_dim_gas_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine gas_data_output_netcdf(gas_data, ncid)

    !> Gas_data to write.
    type(gas_data_t), intent(in) :: gas_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_gas_species

    !> \page output_format_gas_data Output File Format: Gas Material Data
    !!
    !! The gas material data NetCDF dimensions are:
    !!   - \b gas_species: number of gas species
    !!
    !! The gas material data NetCDF variables are:
    !!   - \b gas_species (dim \c gas_species): dummy dimension variable
    !!     (no useful value) - read species names as comma-separated values
    !!     from the 'names' attribute
    !!   - \b gas_mosaic_index (dim \c gas_species): MOSAIC indices of
    !!     gas species
    !!
    !! See also:
    !!   - \ref input_format_gas_data --- the corresponding input format

    call gas_data_netcdf_dim_gas_species(gas_data, ncid, &
         dimid_gas_species)

    call pmc_nc_write_integer_1d(ncid, gas_data%mosaic_index, &
         "gas_mosaic_index", (/ dimid_gas_species /), &
         long_name="MOSAIC indices of gas species")

  end subroutine gas_data_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine gas_data_input_netcdf(gas_data, ncid)

    !> Gas_data to read.
    type(gas_data_t), intent(inout) :: gas_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    character(len=1000) :: name
    integer :: dimid_gas_species, n_spec, varid_gas_species, i_spec, i
    character(len=:), allocatable :: gas_species_names

    call pmc_nc_check(nf90_inq_dimid(ncid, "gas_species", dimid_gas_species))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_gas_species, name, &
         n_spec))
    call assert(719237193, n_spec < 1000)

    if (allocated(gas_data%name)) deallocate(gas_data%name)
    if (allocated(gas_data%mosaic_index)) deallocate(gas_data%mosaic_index)
    allocate(gas_data%name(n_spec))
    allocate(gas_data%mosaic_index(n_spec))

    call pmc_nc_read_integer_1d(ncid, gas_data%mosaic_index, &
         "gas_mosaic_index")

    call pmc_nc_check(nf90_inq_varid(ncid, "gas_species", varid_gas_species))
    allocate(character(len=((GAS_NAME_LEN + 2) * 1000)) :: gas_species_names)
    call pmc_nc_check(nf90_get_att(ncid, varid_gas_species, "names", &
         gas_species_names))
    ! gas_species_names are comma-separated, so unpack them
    do i_spec = 1,gas_data_n_spec(gas_data)
       i = 1
       do while ((gas_species_names(i:i) /= " ") &
            .and. (gas_species_names(i:i) /= ","))
          i = i + 1
       end do
       call assert(173021381, i > 1)
       gas_data%name(i_spec) = gas_species_names(1:(i-1))
       gas_species_names = gas_species_names((i+1):)
    end do
    call assert(469721220, gas_species_names == "")

  end subroutine gas_data_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_gas_data
