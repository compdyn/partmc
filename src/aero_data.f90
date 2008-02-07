! Copyright (C) 2005-2008 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_data module.

!> The aero_data_t structure and associated subroutines.
module pmc_aero_data

  use pmc_inout
  use pmc_mpi
  use pmc_util
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
         "SO4", "NO3", "Cl", "NH4", "MSA", "ARO1", &
         "ARO2", "ALK1", "OLE1", "API1", "API2", "LIM1", &
         "LIM2", "CO3", "Na", "Ca", "OIN", "OC", "BC" ]

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

  !> Write full state.
  subroutine inout_write_aero_data(file, aero_data)
    
    !> File to write to.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data to write.
    type(aero_data_t), intent(in) :: aero_data

    call inout_write_comment(file, "begin aero_data")
    call inout_write_integer(file, "n_spec", aero_data%n_spec)
    call inout_write_integer(file, "i_water", aero_data%i_water)
    call inout_write_string_array(file, "species_names", aero_data%name)
    call inout_write_integer_array(file, "mosaic_indices", &
         aero_data%mosaic_index)
    call inout_write_real_array(file, "rho(kg/m^3)", aero_data%density)
    call inout_write_integer_array(file, "nu", aero_data%num_ions)
    call inout_write_real_array(file, "eps(1)", aero_data%solubility)
    call inout_write_real_array(file, "molec_wght(kg/mole)", &
         aero_data%molec_weight)
    call inout_write_real_array(file, "kappa(1)", aero_data%kappa)
    call inout_write_comment(file, "end aero_data")
    
  end subroutine inout_write_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine inout_read_aero_data(file, aero_data)
    
    !> File to read from.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data to read.
    type(aero_data_t), intent(out) :: aero_data

    call inout_check_comment(file, "begin aero_data")
    call inout_read_integer(file, "n_spec", aero_data%n_spec)
    call inout_read_integer(file, "i_water", aero_data%i_water)
    call inout_read_string_array(file, "species_names", aero_data%name)
    call inout_read_integer_array(file, "mosaic_indices", &
         aero_data%mosaic_index)
    call inout_read_real_array(file, "rho(kg/m^3)", aero_data%density)
    call inout_read_integer_array(file, "nu", aero_data%num_ions)
    call inout_read_real_array(file, "eps(1)", aero_data%solubility)
    call inout_read_real_array(file, "molec_wght(kg/mole)", &
         aero_data%molec_weight)
    call inout_read_real_array(file, "kappa(1)", aero_data%kappa)
    call inout_check_comment(file, "end aero_data")
    
  end subroutine inout_read_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aero_data specification from a inout file.
  subroutine spec_read_aero_data(file, aero_data)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(out) :: aero_data

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)

    call inout_read_real_named_array(file, 0, species_name, species_data)

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

  !> Read aero_data specification from a inout file.
  subroutine spec_read_aero_data_filename(file, aero_data)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(out) :: aero_data

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file

    ! read the aerosol data from the specified file
    call inout_read_string(file, 'aerosol_data', read_name)
    call inout_open_read(read_name, read_file)
    call spec_read_aero_data(read_file, aero_data)
    call inout_close(read_file)

  end subroutine spec_read_aero_data_filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a list of species from the given file with the given name.
  subroutine inout_read_species_list(file, name, aero_data, species_list)

    !> Inout file.
    type(inout_file_t), intent(inout) :: file
    !> Name of line.
    character(len=*), intent(in) :: name
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> List of species numbers.
    integer, pointer :: species_list(:)

    type(inout_line_t) :: line
    integer :: i, spec

    call inout_read_line_no_eof(file, line)
    call inout_check_line_name(file, line, name)
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
    call inout_line_free(line)

  end subroutine inout_read_species_list

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

end module pmc_aero_data
