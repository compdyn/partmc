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

  subroutine aero_data_alloc(aero_data, n_spec)

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

  end subroutine aero_data_alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine aero_data_free(aero_data)

    ! Frees all storage.

    type(aero_data_t), intent(inout) :: aero_data ! aerosol data

    deallocate(aero_data%name)
    deallocate(aero_data%mosaic_index)
    deallocate(aero_data%rho)
    deallocate(aero_data%nu)
    deallocate(aero_data%eps)
    deallocate(aero_data%M_w)

  end subroutine aero_data_free

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

  subroutine set_aero_particle_water_index(aero_data)

    ! Fills in aero_data%i_water.

    type(aero_data_t), intent(inout) :: aero_data  ! aero_data data

    integer :: i

    do i = 1,aero_data%n_spec
       if (aero_data%name(i) == "H2O") then
          aero_data%i_water = i
       end if
    end do

  end subroutine set_aero_particle_water_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_aero_data_mosaic_map(aero_data)

    ! Fills in aero_data%mosaic_index.

    type(aero_data_t), intent(inout) :: aero_data  ! aero_data data

    integer, parameter :: n_mosaic_species = 19
    character*10, parameter, dimension(n_mosaic_species) :: mosaic_species = [ &
         "SO4_a", "NO3_a", "Cl_a", "NH4_a", "CO3_a", "MSA_a", "Na_a", "Ca_a", &
         "OC_a", "BC_a", "OIN_a", "ARO1_a", "ARO2_a", "ALK1_a", "OLE1_a", &
         "API1_a", "API2_a", "LIM1_a", "LIM2_a"]

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

  subroutine inout_write_aero_data(file, aero_data)
    
    ! Write full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to write to
    type(aero_data_t), intent(in) :: aero_data ! aero_data to write

    call inout_write_integer(file, "n_spec", aero_data%n_spec)
    call inout_write_integer(file, "i_water", aero_data%i_water)
    call inout_write_string_array(file, "species_names", aero_data%name)
    call inout_write_integer_array(file, "mosaic_indices", &
         aero_data%mosaic_index)
    call inout_write_real_array(file, "rho(kg/m^3)", aero_data%rho)
    call inout_write_integer_array(file, "nu", aero_data%nu)
    call inout_write_real_array(file, "eps(1)", aero_data%eps)
    call inout_write_real_array(file, "M_w(kg/mole)", aero_data%M_w)
    
  end subroutine inout_write_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inout_read_aero_data(file, aero_data)
    
    ! Read full state.
    
    use mod_inout
    
    type(inout_file_t), intent(inout) :: file ! file to read from
    type(aero_data_t), intent(out) :: aero_data ! aero_data to read

    call inout_read_integer(file, "n_spec", aero_data%n_spec)
    call inout_read_integer(file, "i_water", aero_data%i_water)
    call inout_read_string_array(file, "species_names", aero_data%name)
    call inout_read_integer_array(file, "mosaic_indices", &
         aero_data%mosaic_index)
    call inout_read_real_array(file, "rho(kg/m^3)", aero_data%rho)
    call inout_read_integer_array(file, "nu", aero_data%nu)
    call inout_read_real_array(file, "eps(1)", aero_data%eps)
    call inout_read_real_array(file, "M_w(kg/mole)", aero_data%M_w)
    
  end subroutine inout_read_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_data(file, aero_data)

    ! Read aero_data specification from a inout file.

    use mod_inout

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(out) :: aero_data  ! aero_data data

    integer :: n_species, species, i
    character(len=MAX_CHAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)

    call inout_read_real_named_array(file, 0, species_name, species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (.not. ((size(species_data, 2) == 4) .or. (n_species == 0))) then
       write(0,*) 'ERROR: each line in ', trim(file%name), &
            ' should contain exactly 4 values'
       call exit(1)
    end if

    ! allocate and copy over the data
    call aero_data_alloc(aero_data, n_species)
    do i = 1,n_species
       aero_data%name(i) = species_name(i)
       if (species_name(i) == "H2O") then
          aero_data%i_water = i
       end if
       aero_data%rho(i) = species_data(i,1)
       aero_data%nu(i) = nint(species_data(i,2))
       aero_data%eps(i) = species_data(i,3)
       aero_data%M_w(i) = species_data(i,4)
    end do
    deallocate(species_name)
    deallocate(species_data)
    call set_aero_particle_water_index(aero_data)
    call set_aero_data_mosaic_map(aero_data)

  end subroutine spec_read_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_read_aero_data_filename(file, aero_data)

    ! Read aero_data specification from a inout file.

    use mod_inout

    type(inout_file_t), intent(inout) :: file ! inout file
    type(aero_data_t), intent(out) :: aero_data  ! aero_data data

    character(len=MAX_CHAR_LEN) :: read_name
    type(inout_file_t) :: read_file

    ! read the aerosol data from the specified file
    call inout_read_string(file, 'aerosol_data', read_name)
    call inout_open_read(read_name, read_file)
    call spec_read_aero_data(read_file, aero_data)
    call inout_close(read_file)

  end subroutine spec_read_aero_data_filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_aero_data
