! Copyright (C) 2024 Jeff Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_tchem_interface module.

!> An interface between PartMC and TChem
module pmc_tchem_interface

  use pmc_aero_data
  use pmc_aero_particle
  use pmc_aero_state
  use pmc_constants
  use pmc_gas_data
  use pmc_gas_state
#ifdef PMC_USE_TCHEM
  use iso_c_binding
#endif
  use pmc_util, only : die_msg, warn_assert_msg, assert_msg

interface
  subroutine initialize(arg_chemfile, arg_aerofile, arg_numericsfile, &
       n_batch) &
       bind(c, name="initialize")
    use iso_c_binding
    character(kind=c_char), intent(in) :: arg_chemfile(*)
    character(kind=c_char), intent(in) :: arg_aerofile(*)
    character(kind=c_char), intent(in) :: arg_numericsfile(*)
    integer(c_int), intent(in), value :: n_batch
  end subroutine initialize
  subroutine finalize() bind(c, name="finalize")
  end subroutine finalize
  function TChem_getNumberOfSpecies() bind(c, name="TChem_getNumberOfSpecies")
    use iso_c_binding
    integer(kind=c_int) :: TChem_getNumberOfSpecies
  end function
  function TChem_getNumberOfAeroSpecies() bind(c, &
       name="TChem_getNumberOfAeroSpecies")
    use iso_c_binding
    integer(kind=c_int) :: TChem_getNumberOfAeroSpecies
  end function
  function TChem_getAerosolSpeciesDensity(i_spec) bind(c, &
       name="TChem_getAerosolSpeciesDensity")
    use iso_c_binding
    integer(kind=c_int) :: i_spec
    real(kind=c_double) :: TChem_getAerosolSpeciesDensity
  end function
  function TChem_getAerosolSpeciesMW(i_spec) bind(c, &
       name="TChem_getAerosolSpeciesMW")
    use iso_c_binding
    integer(kind=c_int) :: i_spec
    real(kind=c_double) :: TChem_getAerosolSpeciesMW
  end function
  function TChem_getAerosolSpeciesKappa(i_spec) bind(c, &
       name="TChem_getAerosolSpeciesKappa")
    use iso_c_binding
    integer(kind=c_int) :: i_spec
    real(kind=c_double) :: TChem_getAerosolSpeciesKappa
  end function
  function TChem_getLengthOfStateVector() bind(c, &
      name="TChem_getLengthOfStateVector")
    use iso_c_binding
    integer(kind=c_int) :: TChem_getLengthOfStateVector
  end function
  subroutine TChem_getStateVector(array, i_batch) bind(c, &
       name="TChem_getStateVector")
    use iso_c_binding
    real(kind=c_double) :: array(*)
    integer(c_int), value :: i_batch
  end subroutine
  subroutine TChem_setStateVector(array, i_batch) bind(c, &
       name="TChem_setStateVector")
    use iso_c_binding
    real(kind=c_double) :: array(*)
    integer(c_int), value :: i_batch
  end subroutine
  subroutine TChem_setNumberConcentrationVector(array, i_batch) bind(c, &
       name="TChem_setNumberConcentrationVector")
    use iso_c_binding
    real(kind=c_double) :: array(*)
    integer(c_int), value :: i_batch
  end subroutine
  integer(kind=c_size_t) function TChem_getSpeciesName(index, result, &
       buffer_size) bind(C, name="TChem_getSpeciesName")
     use iso_c_binding
     integer(kind=c_int), intent(in) :: index
     character(kind=c_char), intent(out) :: result(*)
     integer(kind=c_size_t), intent(in), value :: buffer_size
  end function
  integer(kind=c_size_t) function TChem_getAerosolSpeciesName(index, result, &
       buffer_size) bind(C, name="TChem_getAerosolSpeciesName")
     use iso_c_binding
     integer(kind=c_int), intent(in) :: index
     character(kind=c_char), intent(out) :: result(*)
     integer(kind=c_size_t), intent(in), value :: buffer_size
  end function
  subroutine TChem_doTimestep(del_t) bind(C, name="TChem_doTimestep")
     use iso_c_binding
     real(kind=c_double), intent(in) :: del_t
  end subroutine
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_TCHEM
  !> Run the CAMP module for the current PartMC state
  subroutine pmc_tchem_interface_solve(env_state, aero_data, aero_state, &
       gas_data, gas_state, del_t)

    !> Environment data.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time step (s).
    real(kind=dp), intent(in) :: del_t

    call tchem_from_partmc(aero_data, aero_state, gas_data, gas_state, &
         env_state)

    call tchem_timestep(del_t)

    call tchem_to_partmc(aero_data, aero_state, gas_data, gas_state, env_state)

  end subroutine pmc_tchem_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize TChem and PartMC gas and aerosol data.
  subroutine pmc_tchem_initialize(gas_config_filename, aero_config_filename, &
       solver_config_filename, gas_data, aero_data, n_grid_cells)
    use iso_c_binding

    !> Gas configuration filename.
    character(len=*), intent(in) :: gas_config_filename
    !> Aerosol configuration filename.
    character(len=*), intent(in) :: aero_config_filename
    !> Numerical configuration filename.
    character(len=*), intent(in) :: solver_config_filename
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Number of cells to solve.
    integer, intent(in) :: n_grid_cells

    integer(kind=c_int) :: n_aero_spec, n_gas_spec
    integer :: i
    real(kind=c_double), dimension(:), allocatable :: array 
    character(:), allocatable ::  val

    ! initialize the model
    call TChem_initialize(trim(gas_config_filename), &
         trim(aero_config_filename), trim(solver_config_filename), &
         n_grid_cells)

    ! Get size that gas_data should be
    n_gas_spec = TChem_getNumberOfSpecies()
    call ensure_string_array_size(gas_data%name, n_gas_spec)

    ! Populate gas_data with gas species from TChem
    do i = 1,n_gas_spec
       val = TChem_species_name(i-1, .true.)
       gas_data%name(i) = trim(val)
    end do   

    ! For output and MPI, this needs to be allocated (for now)
    allocate(gas_data%mosaic_index(n_gas_spec))
    gas_data%mosaic_index(:) = 0

    ! TODO: Create aero_data based on TChem input.
    ! From TChem we need:
    !   Species names
    !   Species properties - density, kappa, molecular weight
    n_aero_spec = TChem_getNumberOfAeroSpecies()
    call ensure_string_array_size(aero_data%name, n_aero_spec)
    call ensure_integer_array_size(aero_data%mosaic_index, n_aero_spec)
    call ensure_real_array_size(aero_data%wavelengths, n_swbands)
    call ensure_real_array_size(aero_data%density, n_aero_spec)
    call ensure_integer_array_size(aero_data%num_ions, n_aero_spec)
    call ensure_real_array_size(aero_data%molec_weight, n_aero_spec)
    call ensure_real_array_size(aero_data%kappa, n_aero_spec)
    do i = 1,n_aero_spec
       val = TChem_species_name(i-1, .false.)
       aero_data%name(i) =  trim(val)
       aero_data%density(i) = TChem_getAerosolSpeciesDensity(i-1)
       aero_data%molec_weight(i) = TChem_getAerosolSpeciesMW(i-1)
       aero_data%kappa(i) = TChem_getAerosolSpeciesKappa(i-1)
    end do

  end subroutine pmc_tchem_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clean up TChem.
  subroutine pmc_tchem_cleanup()

    call finalize()

  end subroutine pmc_tchem_cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map all data TChem -> PartMC.
  subroutine tchem_to_partmc(aero_data, aero_state, gas_data, gas_state, &
       env_state)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    integer(c_int) :: nSpec, stateVecDim
    integer :: i_part, i_spec
    real(kind=c_double), dimension(:), allocatable :: stateVector 
    integer :: n_gas_spec, n_aero_spec

    n_gas_spec = gas_data_n_spec(gas_data)
    n_aero_spec = aero_data_n_spec(aero_data)

    ! Get gas array
    stateVecDim = TChem_getLengthOfStateVector()
    allocate(stateVector(stateVecDim))
    call TChem_getStateVector(stateVector, 0)

    gas_state%mix_rat = 0.0
    ! Convert from ppm to ppb.
    gas_state%mix_rat = stateVector(4:n_gas_spec+3) * 1000.d0

    do i_part = 1,aero_state_n_part(aero_state)
       do i_spec = 1,n_aero_spec
          aero_state%apa%particle(i_part)%vol(i_spec) = stateVector( &
               n_gas_spec + 3 + i_spec + (i_part -1) * n_aero_spec) &
               / aero_data%density(i_spec)
       end do
    end do

  end subroutine tchem_to_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map all data PartMC -> TChem.
  subroutine tchem_from_partmc(aero_data, aero_state, gas_data, gas_state, &
       env_state)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas State.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp), allocatable :: stateVector(:), number_concentration(:)
    integer :: stateVecDim, tchem_n_part, i_spec
    integer :: i_part
    integer :: i_water

    integer :: n_gas_spec, n_aero_spec

    n_gas_spec = gas_data_n_spec(gas_data)
    n_aero_spec = aero_data_n_spec(aero_data)

    ! Get size of stateVector
    stateVecDim = TChem_getLengthOfStateVector()
    allocate(stateVector(stateVecDim))

    ! Get size of number concentration
    tchem_n_part = 20 
    allocate(number_concentration(tchem_n_part))

    ! First three elements are density, pressure and temperature
    stateVector(1) = env_state_air_den(env_state)
    stateVector(2) = env_state%pressure
    stateVector(3) = env_state%temp

    ! PartMC uses relative humidity and not H2O mixing ratio.
    i_water = gas_data_spec_by_name(gas_data, "H2O")
    gas_state%mix_rat(i_water) = env_state_rel_humid_to_mix_rat(env_state)
    ! Add gas species to state vector. Convert from ppb to ppm.
    stateVector(4:n_gas_spec + 3) = gas_state%mix_rat / 1000.d0

    do i_part = 1,aero_state_n_part(aero_state)
      do i_spec = 1,n_aero_spec
         stateVector(n_gas_spec + 3 + i_spec + (i_part - 1) * n_aero_spec) = &
              aero_particle_species_mass(aero_state%apa%particle(i_part), &
              i_spec, aero_data)
      end do
      number_concentration(i_part) = aero_state_particle_num_conc( &
            aero_state, aero_state%apa%particle(i_part), aero_data)
    end do

    ! FIXME: What do we have to do here for this to work well?
    do i_part = aero_state_n_part(aero_state)+1,tchem_n_part
       do i_spec = 1,n_aero_spec
          stateVector(n_gas_spec + 3 + i_spec + (i_part-1) * n_aero_spec) = &
1d-10
       end do
       number_concentration(i_part) = 0.0d0
    end do

    call TChem_setStateVector(stateVector, 0)

    call TChem_setNumberConcentrationVector(number_concentration, 0)

  end subroutine tchem_from_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a single timestep of TChem chemistry.
  subroutine TChem_timestep(del_t)
    use iso_c_binding

    !> Time step (s).
    real(kind=c_double) :: del_t

    call TChem_doTimestep(del_t)

  end subroutine TChem_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize TChem.
  subroutine TChem_initialize(chemFile, aeroFile, NumericsFile, n_batch)
    use iso_c_binding

    !> Chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: chemFile
    !> Chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: aeroFile
    !> Chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: numericsFile
    !> Number of systems to solve.
    integer(kind=c_int), intent(in) :: n_batch

    call initialize(chemFile//c_null_char, aeroFile//c_null_char, &
         numericsFile//c_null_char, n_batch)

  end subroutine TChem_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get species name from TChem for a given index.
  function TChem_species_name(i_spec, is_gas) result(species_name)
    use iso_c_binding

    !> Index of species.
    integer(kind=c_int), intent(in) :: i_spec
    !> Logical for if the species is a gas
    logical :: is_gas
    !> Species name.
    character(:), allocatable :: species_name

    character(kind=c_char, len=:), allocatable :: cbuf
    integer(kind=c_size_t) :: N

    allocate(character(256) :: cbuf)
    N = len(cbuf)
    if (is_gas) then
       N = TChem_getSpeciesName(i_spec, cbuf, N)
    else
       N = TChem_getAerosolSpeciesName(i_spec, cbuf, N)
    end if
    allocate(character(N) :: species_name)
    species_name = cbuf(:N)

  end function TChem_species_name
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_tchem_interface
