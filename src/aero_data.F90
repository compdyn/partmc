! Copyright (C) 2005-2012, 2016, 2021 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_data module.

!> The aero_data_t structure and associated subroutines.
module pmc_aero_data

  use pmc_spec_file
  use pmc_mpi
  use pmc_util
  use pmc_fractal
  use pmc_netcdf
#ifdef PMC_USE_CAMP
  use camp_util, only: string_t
  use camp_camp_core
  use camp_chem_spec_data
  use camp_aero_rep_data
  use camp_aero_rep_single_particle
  use camp_property
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif

  integer, parameter :: AERO_NAME_LEN = 50
  integer, parameter :: AERO_SOURCE_NAME_LEN = 100

  !> Aerosol material properties and associated data.
  !!
  !! The data in this structure is constant, as it represents physical
  !! quantities that cannot change over time.
  !!
  !! Each aerosol species is identified by an index <tt>i =
  !! 1,...,aero_data_n_spec(aero_data)</tt>. Then \c name(i) is the name of
  !! that species, \c density(i) is its density, etc. The ordering of the
  !! species is arbitrary and should not be relied upon (currently it is the
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
     !> Water species number (0 if water is not a species).
     integer :: i_water
     !> Length \c aero_data_n_spec(aero_data), species names.
     character(len=AERO_NAME_LEN), allocatable :: name(:)
     !> Length \c aero_data_n_spec(aero_data), mosaic_index(i) a positive
     !> integer giving the mosaic index of species i, or 0 if there is no match.
     integer, allocatable :: mosaic_index(:)
     !> Length \c aero_data_n_spec(aero_data), densities (kg m^{-3}).
     real(kind=dp), allocatable ::  density(:)
     !> Length \c aero_data_n_spec(aero_data), number of ions in solute.
     integer, allocatable :: num_ions(:)
     !> Length \c aero_data_n_spec(aero_data), molecular weights (kg mol^{-1}).
     real(kind=dp), allocatable :: molec_weight(:)
     !> Length \c aero_data_n_spec(aero_data), kappas (1).
     real(kind=dp), allocatable :: kappa(:)
     !> Length \c aero_data_n_source(aero_data), source names.
     character(len=AERO_SOURCE_NAME_LEN), allocatable :: source_name(:)
     !> Fractal particle parameters.
     type(fractal_t) :: fractal
#ifdef PMC_USE_CAMP
     !> CAMP aerosol representation pointer
     class(aero_rep_data_t), pointer :: aero_rep_ptr
     !> Aerosol species ids on the camp chem state array for the first
     !! computational particle
     integer, allocatable :: camp_particle_spec_id(:)
     !> Number of elements on the camp chem state array per computational
     !! particle
     integer :: camp_particle_state_size = -1
  contains
     !> Get the index on the CAMP state array for a specified species and
     !! computation particle
     procedure :: camp_spec_id
#endif
  end type aero_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to geometric radius
  !> \f$R_{\rm geo}\f$ (m).
  real(kind=dp) elemental function aero_data_vol2rad(aero_data, v)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    aero_data_vol2rad = fractal_vol2rad(aero_data%fractal, v)

  end function aero_data_vol2rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to geometric diameter
  !> \f$D_{\rm geo}\f$ (m).
  real(kind=dp) elemental function aero_data_vol2diam(aero_data, v)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    aero_data_vol2diam = fractal_vol2diam(aero_data%fractal, v)

  end function aero_data_vol2diam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric radius \f$R_{\rm geo}\f$ (m) to mass-equivalent volume
  !> \f$V\f$ (m^3).
  real(kind=dp) elemental function aero_data_rad2vol(aero_data, r)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Geometric radius (m).
    real(kind=dp), intent(in) :: r

    aero_data_rad2vol = fractal_rad2vol(aero_data%fractal, r)

  end function aero_data_rad2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert geometric diameter \f$D_{\rm geo}\f$ (m) to
  !> mass-equivalent volume \f$V\f$ (m^3).
  real(kind=dp) elemental function aero_data_diam2vol(aero_data, d)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Geometric diameter (m).
    real(kind=dp), intent(in) :: d

    aero_data_diam2vol = fractal_rad2vol(aero_data%fractal, diam2rad(d))

  end function aero_data_diam2vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to number of
  !> monomers \f$N\f$ in a fractal particle cluster.
  !!
  !! Based on Eq. 5 in Naumann [2003].
  real(kind=dp) elemental function aero_data_vol_to_num_of_monomers( &
       aero_data, v)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v

    aero_data_vol_to_num_of_monomers = fractal_vol_to_num_of_monomers( &
         aero_data%fractal, v)

  end function aero_data_vol_to_num_of_monomers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mass-equivalent volume \f$V\f$ (m^3) to mobility equivalent
  !> radius \f$R_{\rm me}\f$ (m).
  !!
  !! Based on Eq. 5, 21 and 30 in Naumann [2003].
  real(kind=dp) function aero_data_vol_to_mobility_rad(aero_data, v, temp, &
       pressure)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Particle mass-equivalent volume (m^3).
    real(kind=dp), intent(in) :: v
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    aero_data_vol_to_mobility_rad = fractal_vol_to_mobility_rad( &
         aero_data%fractal, v, temp, pressure)

  end function aero_data_vol_to_mobility_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius \f$R_{\rm me}\f$ (m) to
  !> mass-equivalent volume \f$V\f$ (m^3).
  real(kind=dp) function aero_data_mobility_rad_to_vol(aero_data, &
       mobility_rad, temp, pressure)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    aero_data_mobility_rad_to_vol = fractal_mobility_rad_to_vol( &
         aero_data%fractal, mobility_rad, temp, pressure)

  end function aero_data_mobility_rad_to_vol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert mobility equivalent radius \f$R_{\rm me}\f$ (m) to
  !> geometric radius \f$R_{\rm geo}\f$ (m^3).
  real(kind=dp) function aero_data_mobility_rad_to_geometric_rad(aero_data, &
       mobility_rad, temp, pressure)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Mobility equivalent radius (m).
    real(kind=dp), intent(in) :: mobility_rad
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Pressure (Pa).
    real(kind=dp), intent(in) :: pressure

    aero_data_mobility_rad_to_geometric_rad &
         = fractal_mobility_rad_to_geometric_rad(aero_data%fractal, &
         mobility_rad, temp, pressure)

  end function aero_data_mobility_rad_to_geometric_rad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the number of aerosol species, or -1 if uninitialized.
  elemental integer function aero_data_n_spec(aero_data)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data

    if (allocated(aero_data%name)) then
       aero_data_n_spec = size(aero_data%name)
    else
       aero_data_n_spec = -1
    end if

  end function aero_data_n_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the number of aerosol sources, or -1 if uninitialized.
  elemental integer function aero_data_n_source(aero_data)

    !> Aero data structure.
    type(aero_data_t), intent(in) :: aero_data

    if (allocated(aero_data%source_name)) then
       aero_data_n_source = size(aero_data%source_name)
    else
       aero_data_n_source = -1
    end if

  end function aero_data_n_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of the species in aero_data with the given name, or
  !> returns 0 if there is no such species.
  integer function aero_data_spec_by_name(aero_data, name)

    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Name of species to find.
    character(len=*), intent(in) :: name

    integer i
    logical found

    found = .false.
    do i = 1,aero_data_n_spec(aero_data)
       if (name == aero_data%name(i)) then
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of the source in aero_data with the given name, or
  !> adds the source if it doesn't exist yet.
  integer function aero_data_source_by_name(aero_data, name)

    !> Aero_data data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Name of source to find.
    character(len=*), intent(in) :: name

    if (.not. allocated(aero_data%source_name)) then
       aero_data%source_name = [name(1:AERO_SOURCE_NAME_LEN)]
       aero_data_source_by_name = 1
       return
    end if
    aero_data_source_by_name = string_array_find(aero_data%source_name, name)
    if (aero_data_source_by_name > 0) return
    aero_data%source_name = [aero_data%source_name, &
         name(1:AERO_SOURCE_NAME_LEN)]
    aero_data_source_by_name = size(aero_data%source_name)

  end function aero_data_source_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fills in aero_data%%i_water.
  subroutine aero_data_set_water_index(aero_data)

    !> Aero_data data.
    type(aero_data_t), intent(inout) :: aero_data

    integer :: i

    aero_data%i_water = string_array_find(aero_data%name, "H2O")

  end subroutine aero_data_set_water_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fills in aero_data%%mosaic_index.
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
    do i_spec = 1,aero_data_n_spec(aero_data)
       aero_data%mosaic_index(i_spec) = string_array_find(mosaic_spec_name, &
            aero_data%name(i_spec))
    end do

  end subroutine aero_data_set_mosaic_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the index on the CAMP state array for a specified species and
  !! computational particle
  integer function camp_spec_id(aero_data, i_part, i_spec)

    !> Aerosol data.
    class(aero_data_t), intent(in) :: aero_data
    !> Computational particle index (1...aero_state_t%n_part).
    integer, intent(in) :: i_part
    !> Aerosol species index in aero_particle_t%vol(:) array.
    integer, intent(in) :: i_spec

#ifdef PMC_USE_CAMP
    call assert(106669451, allocated(aero_data%camp_particle_spec_id))
    call assert(278731889, aero_data%camp_particle_state_size .ge. 0)
    camp_spec_id = (i_part - 1) * aero_data%camp_particle_state_size + &
                   aero_data%camp_particle_spec_id(i_spec)
#else
    camp_spec_id = 0
#endif

  end function camp_spec_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aero_data specification from a spec file.
  subroutine spec_file_read_aero_data(file, aero_data)

    !> Spec file to read data from.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(inout) :: aero_data

    integer :: n_species, species, i
    character(len=SPEC_LINE_MAX_VAR_LEN), allocatable :: species_name(:)
    real(kind=dp), allocatable :: species_data(:,:)

    !> \page input_format_aero_data Input File Format: Aerosol Material Data
    !!
    !! A aerosol material data file must consist of one line per
    !! aerosol species, with each line having:
    !!   - species name (string)
    !!   - density (real, unit kg/m^3)
    !!   - ions per fully dissociated molecule (integer) - used to
    !!     compute kappa value if the corresponding kappa value is
    !!     zero
    !!   - molecular weight (real, unit kg/mol)
    !!   - kappa hygroscopicity parameter (real, dimensionless) - if
    !!     zero, then inferred from the ions value
    !!
    !! This specifies both which species are to be recognized as
    !! aerosol consituents, as well as their physical properties. For
    !! example, an aerosol material data file could contain:
    !! <pre>
    !! # species  dens (kg/m^3)   ions (1)    molec wght (kg/mole)   kappa (1)
    !! SO4        1800            0           96e-3                  0.65
    !! NO3        1800            0           62e-3                  0.65
    !! Cl         2200            1           35.5e-3                0
    !! NH4        1800            0           18e-3                  0.65
    !! </pre>
    !!
    !! Note that it is an error to specify a non-zero number of ions
    !! and a non-zero kappa value for a species. If both values are
    !! zero then that species has zero hygroscopicity parameter. If
    !! exactly one of kappa or ions is non-zero then the non-zero
    !! value is used and the zero value is ignored.
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref output_format_aero_data --- the corresponding output format

    call spec_file_read_real_named_array(file, 0, species_name, species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (.not. ((size(species_data, 2) == 4) .or. (n_species == 0))) then
       call die_msg(428926381, 'each line in ' // trim(file%name) &
            // ' should contain exactly 5 values')
    end if

    ! allocate and copy over the data
    call ensure_string_array_size(aero_data%name, n_species)
    call ensure_integer_array_size(aero_data%mosaic_index, n_species)
    call ensure_real_array_size(aero_data%density, n_species)
    call ensure_integer_array_size(aero_data%num_ions, n_species)
    call ensure_real_array_size(aero_data%molec_weight, n_species)
    call ensure_real_array_size(aero_data%kappa, n_species)
    do i = 1,n_species
       aero_data%name(i) = species_name(i)(1:AERO_NAME_LEN)
       aero_data%density(i) = species_data(i,1)
       aero_data%num_ions(i) = nint(species_data(i,2))
       aero_data%molec_weight(i) = species_data(i,3)
       aero_data%kappa(i) = species_data(i,4)
       call assert_msg(232362742, &
            (aero_data%num_ions(i) == 0) .or. (aero_data%kappa(i) == 0d0), &
            "ions and kappa both non-zero for species " &
            // trim(aero_data%name(i)) // " in " // trim(file%name))
       if (species_name(i) == "H2O") then
          aero_data%i_water = i
          call warn_assert_msg(945800387, &
               aero_data%density(i) == const%water_density, &
               "input H2O density not equal to const%water_density (" &
               // trim(real_to_string(aero_data%density(i))) // " /= " &
               // trim(real_to_string(const%water_density)) // ")")
          call warn_assert_msg(233517437, &
               aero_data%molec_weight(i) == const%water_molec_weight, &
               "input H2O molec_weight not equal " &
               // "to const%water_molec_weight (" &
               // trim(real_to_string(aero_data%molec_weight(i))) // " /= " &
               // trim(real_to_string(const%water_molec_weight)) // ")")
       end if
    end do
    call aero_data_set_water_index(aero_data)
    call aero_data_set_mosaic_map(aero_data)

  end subroutine spec_file_read_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a list of species from the given file with the given name.
  subroutine spec_file_read_species_list(file, name, aero_data, species_list)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Name of line.
    character(len=*), intent(in) :: name
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> List of species numbers.
    integer, allocatable :: species_list(:)

    type(spec_line_t) :: line
    integer :: i, spec

    call spec_file_read_line_no_eof(file, line)
    call spec_file_check_line_name(file, line, name)
    call ensure_integer_array_size(species_list, size(line%data))
    do i = 1,size(line%data)
       spec = aero_data_spec_by_name(aero_data, line%data(i))
       if (spec == 0) then
          call spec_file_die_msg(964771833, file, &
               'unknown species: ' // trim(line%data(i)))
       end if
       species_list(i) = spec
    end do

  end subroutine spec_file_read_species_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_data(val)

    !> Value to pack.
    type(aero_data_t), intent(in) :: val

    pmc_mpi_pack_size_aero_data = &
         pmc_mpi_pack_size_integer(val%i_water) &
         + pmc_mpi_pack_size_string_array(val%name) &
         + pmc_mpi_pack_size_integer_array(val%mosaic_index) &
         + pmc_mpi_pack_size_real_array(val%density) &
         + pmc_mpi_pack_size_integer_array(val%num_ions) &
         + pmc_mpi_pack_size_real_array(val%molec_weight) &
         + pmc_mpi_pack_size_real_array(val%kappa) &
         + pmc_mpi_pack_size_string_array(val%source_name) &
         + pmc_mpi_pack_size_fractal(val%fractal)

  end function pmc_mpi_pack_size_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    call pmc_mpi_pack_integer(buffer, position, val%i_water)
    call pmc_mpi_pack_string_array(buffer, position, val%name)
    call pmc_mpi_pack_integer_array(buffer, position, val%mosaic_index)
    call pmc_mpi_pack_real_array(buffer, position, val%density)
    call pmc_mpi_pack_integer_array(buffer, position, val%num_ions)
    call pmc_mpi_pack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_pack_real_array(buffer, position, val%kappa)
    call pmc_mpi_pack_string_array(buffer, position, val%source_name)
    call pmc_mpi_pack_fractal(buffer, position, val%fractal)
    call assert(183834856, &
         position - prev_position <= pmc_mpi_pack_size_aero_data(val))
#endif

  end subroutine pmc_mpi_pack_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_data_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%i_water)
    call pmc_mpi_unpack_string_array(buffer, position, val%name)
    call pmc_mpi_unpack_integer_array(buffer, position, val%mosaic_index)
    call pmc_mpi_unpack_real_array(buffer, position, val%density)
    call pmc_mpi_unpack_integer_array(buffer, position, val%num_ions)
    call pmc_mpi_unpack_real_array(buffer, position, val%molec_weight)
    call pmc_mpi_unpack_real_array(buffer, position, val%kappa)
    call pmc_mpi_unpack_string_array(buffer, position, val%source_name)
    call pmc_mpi_unpack_fractal(buffer, position, val%fractal)
    call assert(188522823, &
         position - prev_position <= pmc_mpi_pack_size_aero_data(val))
#endif

  end subroutine pmc_mpi_unpack_aero_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    integer :: aero_species_centers(aero_data_n_spec(aero_data))
    character(len=(AERO_NAME_LEN * aero_data_n_spec(aero_data))) :: &
         aero_species_names

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_species", dimid_aero_species)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_species", &
         aero_data_n_spec(aero_data), dimid_aero_species))
    aero_species_names = ""
    do i_spec = 1,aero_data_n_spec(aero_data)
       aero_species_names((len_trim(aero_species_names) + 1):) &
            = trim(aero_data%name(i_spec))
       if (i_spec < aero_data_n_spec(aero_data)) then
          aero_species_names((len_trim(aero_species_names) + 1):) = ","
       end if
    end do
    call pmc_nc_check(nf90_def_var(ncid, "aero_species", NF90_INT, &
         dimid_aero_species, varid_aero_species))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_species, "names", &
         aero_species_names))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_species, "description", &
         "dummy dimension variable (no useful value) - read species names " &
         // "as comma-separated values from the 'names' attribute"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_spec = 1,aero_data_n_spec(aero_data)
       aero_species_centers(i_spec) = i_spec
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_species, &
         aero_species_centers))

  end subroutine aero_data_netcdf_dim_aero_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the aero source dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine aero_data_netcdf_dim_aero_source(aero_data, ncid, &
       dimid_aero_source)

    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the source dimension.
    integer, intent(out) :: dimid_aero_source

    integer :: status, i_source
    integer :: varid_aero_source
    integer :: aero_source_centers(aero_data_n_source(aero_data))
    character(len=(AERO_SOURCE_NAME_LEN * aero_data_n_source(aero_data))) &
         :: aero_source_names

    ! try to get the dimension ID
    status = nf90_inq_dimid(ncid, "aero_source", dimid_aero_source)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_source", &
         aero_data_n_source(aero_data), dimid_aero_source))
    aero_source_names = ""
    do i_source = 1,aero_data_n_source(aero_data)
       aero_source_names((len_trim(aero_source_names) + 1):) &
            = trim(aero_data%source_name(i_source))
       if (i_source < aero_data_n_source(aero_data)) then
          aero_source_names((len_trim(aero_source_names) + 1):) = ","
       end if
    end do
    call pmc_nc_check(nf90_def_var(ncid, "aero_source", NF90_INT, &
         dimid_aero_source, varid_aero_source))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_source, "names", &
         aero_source_names))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_source, "description", &
         "dummy dimension variable (no useful value) - read source names " &
         // "as comma-separated values from the 'names' attribute"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_source = 1,aero_data_n_source(aero_data)
       aero_source_centers(i_source) = i_source
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_source, &
         aero_source_centers))

  end subroutine aero_data_netcdf_dim_aero_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write full state.
  subroutine aero_data_output_netcdf(aero_data, ncid)

    !> Aero_data to write.
    type(aero_data_t), intent(in) :: aero_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer :: dimid_aero_species, dimid_aero_source

    !> \page output_format_aero_data Output File Format: Aerosol Material Data
    !!
    !! The aerosol material data NetCDF dimensions are:
    !!   - \b aero_species: number of aerosol species
    !!   - \b aero_source: number of aerosol sources
    !!
    !! The aerosol material data NetCDF variables are:
    !!   - \b aero_species (dim \c aero_species): dummy dimension variable
    !!     (no useful value) - read species names as comma-separated values
    !!     from the 'names' attribute
    !!   - \b aero_source (dim \c aero_source): dummy dimension variable
    !!     (no useful value) - read source names as comma-separated values
    !!     from the 'names' attribute
    !!   - \b aero_mosaic_index (dim \c aero_species): indices of species
    !!     in MOSAIC
    !!   - \b aero_density (unit kg/m^3, dim \c aero_species): densities
    !!     of aerosol species
    !!   - \b aero_num_ions (dim \c aero_species): number of ions produced
    !!     when one molecule of each species fully dissociates in water
    !!   - \b aero_molec_weight (unit kg/mol, dim \c aero_species): molecular
    !!     weights of aerosol species
    !!   - \b aero_kappa (unit kg/mol, dim \c aero_species): hygroscopicity
    !!     parameters of aerosol species
    !!   - \b fractal parameters, see \ref output_format_fractal
    !!
    !! See also:
    !!   - \ref input_format_aero_data --- the corresponding input format

    call aero_data_netcdf_dim_aero_species(aero_data, ncid, &
         dimid_aero_species)
    call aero_data_netcdf_dim_aero_source(aero_data, ncid, &
         dimid_aero_source)

    call pmc_nc_write_integer_1d(ncid, aero_data%mosaic_index, &
         "aero_mosaic_index", (/ dimid_aero_species /), &
         long_name="MOSAIC indices of aerosol species")
    call pmc_nc_write_real_1d(ncid, aero_data%density, &
         "aero_density", (/ dimid_aero_species /), unit="kg/m^3", &
         long_name="densities of aerosol species")
    call pmc_nc_write_integer_1d(ncid, aero_data%num_ions, &
         "aero_num_ions", (/ dimid_aero_species /), &
         long_name="number of ions after dissociation of aerosol species")
    call pmc_nc_write_real_1d(ncid, aero_data%molec_weight, &
         "aero_molec_weight", (/ dimid_aero_species /), unit="kg/mol", &
         long_name="molecular weights of aerosol species")
    call pmc_nc_write_real_1d(ncid, aero_data%kappa, &
         "aero_kappa", (/ dimid_aero_species /), unit="1", &
         long_name="hygroscopicity parameters (kappas) of aerosol species")
    call fractal_output_netcdf(aero_data%fractal, ncid)

  end subroutine aero_data_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read full state.
  subroutine aero_data_input_netcdf(aero_data, ncid)

    !> Aero_data to read.
    type(aero_data_t), intent(inout) :: aero_data
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    integer, parameter :: MAX_SPECIES = 1000
    integer, parameter :: MAX_SOURCES = 1000

    character(len=1000) :: name
    integer :: dimid_aero_species, n_spec, varid_aero_species, i_spec, i
    integer :: dimid_aero_source, n_source, varid_aero_source, i_source
    character(len=((AERO_NAME_LEN + 2) * MAX_SPECIES)) :: aero_species_names
    character(len=:), allocatable :: aero_source_names

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_species", &
         dimid_aero_species))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, &
         dimid_aero_species, name, n_spec))
    call assert(141013948, n_spec < MAX_SPECIES)

    call pmc_nc_check(nf90_inq_dimid(ncid, "aero_source", &
         dimid_aero_source))
    call pmc_nc_check(nf90_Inquire_Dimension(ncid, &
         dimid_aero_source, name, n_source))
    call assert(739238793, n_source < MAX_SOURCES)

    call pmc_nc_read_integer_1d(ncid, aero_data%mosaic_index, &
         "aero_mosaic_index")
    call pmc_nc_read_real_1d(ncid, aero_data%density, "aero_density")
    call pmc_nc_read_integer_1d(ncid, aero_data%num_ions, "aero_num_ions")
    call pmc_nc_read_real_1d(ncid, aero_data%molec_weight, "aero_molec_weight")
    call pmc_nc_read_real_1d(ncid, aero_data%kappa, "aero_kappa")

    call pmc_nc_check(nf90_inq_varid(ncid, "aero_species", &
         varid_aero_species))
    call pmc_nc_check(nf90_get_att(ncid, varid_aero_species, "names", &
         aero_species_names))
    ! aero_species_names are comma-separated, so unpack them
    call ensure_string_array_size(aero_data%name, n_spec)
    do i_spec = 1,aero_data_n_spec(aero_data)
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

    call pmc_nc_check(nf90_inq_varid(ncid, "aero_source", &
         varid_aero_source))
    allocate(character(len=((AERO_SOURCE_NAME_LEN + 2) * MAX_SPECIES)) &
         :: aero_source_names)
    call pmc_nc_check(nf90_get_att(ncid, varid_aero_source, "names", &
         aero_source_names))
    ! aero_source_names are comma-separated, so unpack them
    call ensure_string_array_size(aero_data%source_name, n_source)
    do i_source = 1,aero_data_n_source(aero_data)
       i = 1
       do while ((aero_source_names(i:i) /= " ") &
            .and. (aero_source_names(i:i) /= ","))
          i = i + 1
       end do
       call assert(840982478, i > 1)
       aero_data%source_name(i_source) = aero_source_names(1:(i-1))
       aero_source_names = aero_source_names((i+1):)
    end do
    call assert(377166446, aero_source_names == "")

    call aero_data_set_water_index(aero_data)

    call fractal_input_netcdf(aero_data%fractal, ncid)

  end subroutine aero_data_input_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef PMC_USE_CAMP
  !> Initialize the aero_data_t variable with camp chem data
  subroutine aero_data_initialize(aero_data, camp_core)

    !> Aerosol data.
    class(aero_data_t), intent(inout) :: aero_data
    !> CAMP core.
    type(camp_core_t), intent(in) :: camp_core

    character(len=:), allocatable :: rep_name, prop_name, str_val
    type(string_t), allocatable :: spec_names(:), tmp_spec_names(:)
    integer :: num_spec, i_spec, spec_type
    type(chem_spec_data_t), pointer :: chem_spec_data
    type(property_t), pointer :: property_set

    rep_name = "PartMC single particle"
    if (.not. camp_core%get_aero_rep(rep_name, aero_data%aero_rep_ptr)) then
      call die_msg(418509983, "Missing 'PartMC single particle' aerosol "// &
           "representation.")
    end if

    call assert_msg(935419266, camp_core%get_chem_spec_data(chem_spec_data), &
         "No chemical species data in camp_core.")

    ! Only include real aerosol species (no activity coefficients)
    spec_names = aero_data%aero_rep_ptr%unique_names()
    allocate(tmp_spec_names(size(spec_names)))
    num_spec = 0
    do i_spec = 1, size(spec_names)
       call assert(496388827, chem_spec_data%get_type( &
            aero_data%aero_rep_ptr%spec_name(spec_names(i_spec)%string), &
            spec_type))
       if (spec_type /= CHEM_SPEC_VARIABLE .and. &
           spec_type /= CHEM_SPEC_CONSTANT .and. &
           spec_type /= CHEM_SPEC_PSSA) cycle
       if (spec_names(i_spec)%string(1:3) /= "P1.") exit
       num_spec = num_spec + 1
       tmp_spec_names(num_spec)%string = spec_names(i_spec)%string(4:)
    end do
    deallocate(spec_names)
    allocate(spec_names(num_spec))
    spec_names(:) = tmp_spec_names(1:num_spec)
    deallocate(tmp_spec_names)

    allocate(aero_data%name(num_spec))
    allocate(aero_data%mosaic_index(num_spec))
    allocate(aero_data%density(num_spec))
    allocate(aero_data%num_ions(num_spec))
    allocate(aero_data%molec_weight(num_spec))
    allocate(aero_data%kappa(num_spec))
    allocate(aero_data%camp_particle_spec_id(num_spec))

    ! Assume no aerosol water
    aero_data%i_water = 0

    do i_spec = 1, num_spec
       aero_data%name(i_spec) = spec_names(i_spec)%string
       if (.not. chem_spec_data%get_property_set( &
            aero_data%aero_rep_ptr%spec_name("P1." &
            // spec_names(i_spec)%string), property_set)) then
          call die_msg(934844845, "Missing property set for aerosol species " &
               // spec_names(i_spec)%string)
       end if
       prop_name = "density [kg m-3]"
       if (.not. property_set%get_real(prop_name, &
            aero_data%density(i_spec))) then
          call die_msg(547508215, "Missing density for aerosol species " &
               // spec_names(i_spec)%string)
       end if
       prop_name = "num_ions"
       if (.not. property_set%get_int(prop_name, &
            aero_data%num_ions(i_spec))) then
          call die_msg(324777059, "Missing num_ions for aerosol species " &
               // spec_names(i_spec)%string)
       end if
       prop_name = "molecular weight [kg mol-1]"
       if (.not. property_set%get_real(prop_name, &
            aero_data%molec_weight(i_spec))) then
          call die_msg(549413749, "Missing molec_weight for aerosol species " &
               // spec_names(i_spec)%string)
       end if
       prop_name = "kappa"
       if (.not. property_set%get_real(prop_name, &
            aero_data%kappa(i_spec))) then
         call die_msg(944207343, "Missing kappa for aerosol species " &
              // spec_names(i_spec)%string)
       end if
       prop_name = "PartMC name"
       if (property_set%get_string(prop_name, str_val)) then
          if (str_val == "H2O") then
             call assert_msg(227489086, aero_data%i_water == 0, &
                  "Multiple aerosol water species")
             aero_data%i_water = i_spec
          end if
       end if
       aero_data%camp_particle_spec_id(i_spec) = &
            aero_data%aero_rep_ptr%spec_state_id("P1." &
            // spec_names(i_spec)%string)
    end do

    select type( aero_rep => aero_data%aero_rep_ptr)
       type is(aero_rep_single_particle_t)

          ! Get the number of elements per-particle on the CAMP state array
          aero_data%camp_particle_state_size = aero_rep%per_particle_size()

       class default
          call die_msg(281737350, "Wrong aerosol representation type")
    end select

  end subroutine aero_data_initialize
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_data
