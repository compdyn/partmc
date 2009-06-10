! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_data module.

!> The aero_data_t structure and associated subroutines.
module pmc_aero_data

  use pmc_spec_read
  use pmc_mpi
  use pmc_util
  use pmc_netcdf
#ifdef PMC_USE_MPI
  use mpi
#endif

  integer, parameter :: AERO_NAME_LEN = 15

  !> Aerosol material properties and associated data.
  !!
  !! The data in this structure is constant, as it represents physical
  !! quantities that cannot change over time.
  !!
  !! Each aerosol species is identified by an index <tt>i =
  !! 1,...,n_spec</tt>. Then \c name(i) is the name of that species,
  !! \c density(i) is its density, etc. The ordering of the species is
  !! arbitrary and should not be relied upon (currently it is the
  !! order in the species data file). The only exception is that it is
  !! possible to find out which species is water from the \c i_water
  !! variable.
  !!
  !! The names of the aerosol species are not important to PartMC, as
  !! only the material properties are used. The names are used for
  !! input and output, and also for communication with MOSAIC. For the
  !! MOSAIC interface to work correctly the species must be named the
  !! same, but without the \c _a suffix.
  type aero_data_t
     !> Number of species.
     integer :: n_spec
     !> Water species number.
     integer :: i_water
     !> Len n_spec, species.
     character(len=AERO_NAME_LEN), pointer :: name(:)
     !> Length n_spec, mosaic_index(i) a positive integer giving the
     !> mosaic index of species i, or 0 if there is no match.
     integer, pointer :: mosaic_index(:)
     !> Len n_spec, densities (kg m^{-3}).
     real*8, pointer ::  density(:)
     !> Len n_spec, num ions in solute.
     integer, pointer :: num_ions(:)
     !> Len n_spec, solubilities (1).
     real*8, pointer :: solubility(:)
     !> Len n_spec, molc wghts (kg mole^{-1}).
     real*8, pointer :: molec_weight(:)
     !> Len n_spec, kappas (1).
     real*8, pointer :: kappa(:)
  end type aero_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for aero_data parameters given the number of
  !> species.
  subroutine aero_data_alloc(aero_data, n_spec)

    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Number of species.
    integer, intent(in) :: n_spec

    aero_data%n_spec = n_spec
    allocate(aero_data%name(n_spec))
    allocate(aero_data%mosaic_index(n_spec))
    allocate(aero_data%density(n_spec))
    allocate(aero_data%num_ions(n_spec))
    allocate(aero_data%solubility(n_spec))
    allocate(aero_data%molec_weight(n_spec))
    allocate(aero_data%kappa(n_spec))
    aero_data%i_water = 0

  end subroutine aero_data_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Frees all storage.
  subroutine aero_data_free(aero_data)

    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data

    deallocate(aero_data%name)
    deallocate(aero_data%mosaic_index)
    deallocate(aero_data%density)
    deallocate(aero_data%num_ions)
    deallocate(aero_data%solubility)
    deallocate(aero_data%molec_weight)
    deallocate(aero_data%kappa)

  end subroutine aero_data_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of the species in aero_data with the given name, or
  !> returns 0 if there is no such species.
  integer function aero_data_spec_by_name(aero_data, name)

    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Name of species to find.
    character(len=AERO_NAME_LEN), intent(in) :: name

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

  !> Fills in aero_data%i_water.
  subroutine aero_data_set_water_index(aero_data)

    !> Aero_data data.
    type(aero_data_t), intent(inout) :: aero_data

    integer :: i

    do i = 1,aero_data%n_spec
       if (aero_data%name(i) == "H2O") then
          aero_data%i_water = i
       end if
    end do

  end subroutine aero_data_set_water_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fills in aero_data%mosaic_index.
  subroutine aero_data_set_mosaic_map(aero_data)

    !> Aero_data data.
    type(aero_data_t), intent(inout) :: aero_data

    integer, parameter :: n_mosaic_spec = 19
    character(AERO_NAME_LEN), parameter, dimension(n_mosaic_spec) :: &
         mosaic_spec_name = [ &
         "SO4   ", "NO3   ", "Cl    ", "NH4   ", "MSA   ", "ARO1  ", &
         "ARO2  ", "ALK1  ", "OLE1  ", "API1  ", "API2  ", "LIM1  ", &
         "LIM2  ", "CO3   ", "Na    ", "Ca    ", "OIN   ", "OC    ", &
         "BC    "]

    integer :: i_spec, i_mosaic_spec, i

    aero_data%mosaic_index = 0
    do i_spec = 1,aero_data%n_spec
       i_mosaic_spec = 0
       do i = 1,n_mosaic_spec
          if (aero_data%name(i_spec) == mosaic_spec_name(i)) then
             i_mosaic_spec = i
          end if
       end do
       if (i_mosaic_spec > 0) then
          aero_data%mosaic_index(i_spec) = i_mosaic_spec
       end if
    end do

  end subroutine aero_data_set_mosaic_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aero_data specification from a spec file.
  subroutine spec_read_aero_data(file, aero_data)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(out) :: aero_data

    integer :: n_species, species, i
    character(len=MAX_VAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)

    call spec_read_real_named_array(file, 0, species_name, species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (.not. ((size(species_data, 2) == 5) .or. (n_species == 0))) then
       write(0,*) 'ERROR: each line in ', trim(file%name), &
            ' should contain exactly 5 values'
       call exit(1)
    end if

    ! allocate and copy over the data
    call aero_data_alloc(aero_data, n_species)
    do i = 1,n_species
       aero_data%name(i) = species_name(i)(1:AERO_NAME_LEN)
       if (species_name(i) == "H2O") then
          aero_data%i_water = i
       end if
       aero_data%density(i) = species_data(i,1)
       aero_data%num_ions(i) = nint(species_data(i,2))
       aero_data%solubility(i) = species_data(i,3)
       aero_data%molec_weight(i) = species_data(i,4)
       aero_data%kappa(i) = species_data(i,5)
    end do
    deallocate(species_name)
    deallocate(species_data)
    call aero_data_set_water_index(aero_data)
    call aero_data_set_mosaic_map(aero_data)

  end subroutine spec_read_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aero_data specification from a spec file.
  subroutine spec_read_aero_data_filename(file, aero_data)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(out) :: aero_data

    character(len=MAX_VAR_LEN) :: read_name
    type(spec_file_t) :: read_file

    ! read the aerosol data from the specified file
    call spec_read_string(file, 'aerosol_data', read_name)
    call spec_read_open(read_name, read_file)
    call spec_read_aero_data(read_file, aero_data)
    call spec_read_close(read_file)

  end subroutine spec_read_aero_data_filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a list of species from the given file with the given name.
  subroutine spec_read_species_list(file, name, aero_data, species_list)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name of line.
    character(len=*), intent(in) :: name
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> List of species numbers.
    integer, pointer :: species_list(:)

    type(spec_line_t) :: line
    integer :: i, spec

    call spec_read_line_no_eof(file, line)
    call spec_read_check_line_name(file, line, name)
    allocate(species_list(size(line%data)))
    do i = 1,size(line%data)
       spec = aero_data_spec_by_name(aero_data, line%data(i))
       if (spec == 0) then
          write(0,*) 'ERROR: unknown species ', trim(line%data(i)), &
               ' on line ', file%line_num, ' of file ', trim(file%name)
          call exit(1)
       end if
       species_list(i) = spec
    end do
    call spec_line_free(line)

  end subroutine spec_read_species_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_data(val)

    !> Value to pack.
    type(aero_data_t), intent(in) :: val

    pmc_mpi_pack_size_aero_data = &
         pmc_mpi_pack_size_integer(val%n_spec) &
         + pmc_mpi_pack_size_integer(val%i_water) &
         + pmc_mpi_pack_size_string_array(val%name) &
         + pmc_mpi_pack_size_integer_array(val%mosaic_index) &
         + pmc_mpi_pack_size_real_array(val%density) &
         + pmc_mpi_pack_size_integer_array(val%num_ions) &
         + pmc_mpi_pack_size_real_array(val%solubility) &
         + pmc_mpi_pack_size_real_array(val%molec_weight) &
         + pmc_mpi_pack_size_real_array(val%kappa)

  end function pmc_mpi_pack_size_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_data_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_spec)
    call pmc_mpi_pack_integer(buffer, position, val%i_water)
    call pmc_mpi_pack_string_array(buffer, position, val%name)
    call pmc_mpi_pack_integer_array(buffer, position, val%mosaic_index)
    call pmc_mpi_pack_real_array(buffer, position, val%density)
    call pmc_mpi_pack_integer_array(buffer, position, val%num_ions)
    call pmc_mpi_pack_real_array(buffer, position, val%solubility)
    call pmc_mpi_pack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_pack_real_array(buffer, position, val%kappa)
    call assert(183834856, &
         position - prev_position == pmc_mpi_pack_size_aero_data(val))
#endif

  end subroutine pmc_mpi_pack_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_data_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_spec)
    call pmc_mpi_unpack_integer(buffer, position, val%i_water)
    call pmc_mpi_unpack_string_array(buffer, position, val%name)
    call pmc_mpi_unpack_integer_array(buffer, position, val%mosaic_index)
    call pmc_mpi_unpack_real_array(buffer, position, val%density)
    call pmc_mpi_unpack_integer_array(buffer, position, val%num_ions)
    call pmc_mpi_unpack_real_array(buffer, position, val%solubility)
    call pmc_mpi_unpack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_unpack_real_array(buffer, position, val%kappa)
    call assert(188522823, &
         position - prev_position == pmc_mpi_pack_size_aero_data(val))
#endif

  end subroutine pmc_mpi_unpack_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the aero species dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine aero_data_netcdf_dim_aero_species(aero_data, ncid, &
       dimid_aero_species)

    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_aero_species

    integer :: status, i_spec
    integer :: varid_aero_species
    integer :: aero_species_centers(aero_data%n_spec)
    character(len=(AERO_NAME_LEN * aero_data%n_spec)) :: aero_species_names

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_species", dimid_aero_species)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_species", &
         aero_data%n_spec, dimid_aero_species))
    aero_species_names = ""
    do i_spec = 1,aero_data%n_spec
       aero_species_names((len_trim(aero_species_names) + 1):) &
            = trim(aero_data%name(i_spec))
       if (i_spec < aero_data%n_spec) then
          aero_species_names((len_trim(aero_species_names) + 1):) = ","
       end if
    end do
    call pmc_nc_check(nf90_def_var(ncid, "aero_species", NF90_INT, &
         dimid_aero_species, varid_aero_species))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_species, "unit", "1"))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_species, "names", &
         aero_species_names))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_spec = 1,aero_data%n_spec
       aero_species_centers(i_spec) = i_spec
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_species, &
         aero_species_centers))

  end subroutine aero_data_netcdf_dim_aero_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_data_output_netcdf(aero_data, ncid)
    
    !> Aero_data to write.
    type(aero_data_t), intent(in) :: aero_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_species

    call aero_data_netcdf_dim_aero_species(aero_data, ncid, &
         dimid_aero_species)

    call pmc_nc_write_integer_1d(ncid, aero_data%mosaic_index, &
         "aero_mosaic_index", "1", (/ dimid_aero_species /))
    call pmc_nc_write_real_1d(ncid, aero_data%density, &
         "aero_density", "kg/m^3", (/ dimid_aero_species /))
    call pmc_nc_write_integer_1d(ncid, aero_data%num_ions, &
         "aero_num_ions", "1", (/ dimid_aero_species /))
    call pmc_nc_write_real_1d(ncid, aero_data%solubility, &
         "aero_solubility", "1", (/ dimid_aero_species /))
    call pmc_nc_write_real_1d(ncid, aero_data%molec_weight, &
         "aero_molec_weight", "kg/mole", (/ dimid_aero_species /))
    call pmc_nc_write_real_1d(ncid, aero_data%kappa, &
         "aero_kappa", "1", (/ dimid_aero_species /))

  end subroutine aero_data_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_data_input_netcdf(aero_data, ncid)
    
    !> Aero_data to read.
    type(aero_data_t), intent(inout) :: aero_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    character(len=1000) :: unit, name
    integer :: dimid_aero_species, n_spec, varid_aero_species, i_spec, i
    character(len=((AERO_NAME_LEN + 2) * 1000)) :: aero_species_names

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_species", dimid_aero_species))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_species, name, n_spec))
    call aero_data_free(aero_data)
    call aero_data_alloc(aero_data, n_spec)
    call assert(739238793, n_spec < 1000)

    call pmc_nc_read_integer_1d(ncid, aero_data%mosaic_index, &
         "aero_mosaic_index", unit)
    call pmc_nc_read_real_1d(ncid, aero_data%density, &
         "aero_density", unit)
    call pmc_nc_read_integer_1d(ncid, aero_data%num_ions, &
         "aero_num_ions", unit)
    call pmc_nc_read_real_1d(ncid, aero_data%solubility, &
         "aero_solubility", unit)
    call pmc_nc_read_real_1d(ncid, aero_data%molec_weight, &
         "aero_molec_weight", unit)
    call pmc_nc_read_real_1d(ncid, aero_data%kappa, &
         "aero_kappa", unit)

    call pmc_nc_check(nf90_inq_varid(ncid, "aero_species", varid_aero_species))
    call pmc_nc_check(nf90_get_att(ncid, varid_aero_species, "names", aero_species_names))
    ! aero_species_names are comma-separated, so unpack them
    do i_spec = 1,aero_data%n_spec
       i = 1
       do while ((aero_species_names(i:i) /= " ") &
            .and. (aero_species_names(i:i) /= ","))
          i = i + 1
       end do
       call assert(852937292, i > 1)
       aero_data%name(i_spec) = aero_species_names(1:(i-1))
       aero_species_names = aero_species_names((i+1):)
    end do
    call assert(729138192, aero_species_names == "")

    call aero_data_set_water_index(aero_data)

  end subroutine aero_data_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_data
