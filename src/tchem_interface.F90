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
  subroutine initialize(arg_chemfile, arg_aerofile, arg_numericsfile) &
       bind(c, name="initialize")
    use iso_c_binding
    character(kind=c_char), intent(in) :: arg_chemfile(*)
    character(kind=c_char), intent(in) :: arg_aerofile(*)
    character(kind=c_char), intent(in) :: arg_numericsfile(*)
  end subroutine initialize
  subroutine finalize() bind(c, name="finalize")
  end subroutine finalize
  function TChem_getNumberOfSpecies() bind(c, name="TChem_getNumberOfSpecies")
    use iso_c_binding
    integer(kind=c_int) :: TChem_getNumberOfSpecies
  end function
  function TChem_getLengthOfStateVector() bind(c, &
      name="TChem_getLengthOfStateVector")
    use iso_c_binding
    integer(kind=c_int) :: TChem_getLengthOfStateVector
  end function
  subroutine TChem_getStateVector(array) bind(c, name="TChem_getStateVector")
    use iso_c_binding
    real(kind=c_double) :: array(*)
  end subroutine
  subroutine TChem_setStateVector(array) bind(c, name="TChem_setStateVector")
    use iso_c_binding
    real(kind=c_double) :: array(*)
  end subroutine
  integer(kind=c_size_t) function TChem_getSpeciesName(index, result, buffer_size) &
       bind(C, name="TChem_getSpeciesName")
     use iso_c_binding
     integer(kind=c_int), intent(in) :: index
     character(kind=c_char), intent(out) :: result(*)
     integer(kind=c_size_t), intent(in), value :: buffer_size
  end function
  subroutine TChem_doTimestep(del_t) bind(C, name="TChem_doTimestep")
     use iso_c_binding
     real(kind=c_double), intent(in) :: del_t
  end subroutine
  subroutine TChem_getAllStateVectorHost(arg_state_vector) bind(C, &
      name="TChem_getAllStateVectorHost")
    use iso_c_binding, only: c_ptr
    type(c_ptr), value :: arg_state_vector
  end subroutine TChem_getAllStateVectorHost
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
       solver_config_filename, gas_data, aero_data)
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

    integer(kind=c_int) :: nSpec, nAeroSpec
    integer :: n_species
    integer :: i
    real(kind=c_double), dimension(:), allocatable :: array 
    character(:), allocatable ::  val

    ! initialize the model
    call TChem_initialize(trim(gas_config_filename), &
         trim(aero_config_filename), trim(solver_config_filename))

    ! what we need:
    !  - number of gas species for gas_data
    !  - name of gas species to set gas_data
    nSpec = TChem_getNumberOfSpecies() 
    call ensure_string_array_size(gas_data%name, nSpec)

    ! Populate gas_data with gas species from TChem
    do i = 1,nSpec
       val = TChem_species_name(i-1) 
       gas_data%name(i) = trim(val)
    end do   

    ! FIXME: We need this still
    allocate(gas_data%mosaic_index(gas_data_n_spec(gas_data)))
    gas_data%mosaic_index(:) = 0
 
    print*, 'in partmc', nSpec, trim(gas_config_filename), &
         gas_data_n_spec(gas_data)

    ! Get aerosol data
    ! Species names
    ! Species properties - density, kappa, molecular weight
    n_species = 10
    call ensure_string_array_size(aero_data%name, n_species)
    call ensure_integer_array_size(aero_data%mosaic_index, n_species)
    call ensure_real_array_size(aero_data%wavelengths, n_swbands)
    call ensure_real_array_size(aero_data%density, n_species)
    call ensure_integer_array_size(aero_data%num_ions, n_species)
    call ensure_real_array_size(aero_data%molec_weight, n_species)
    call ensure_real_array_size(aero_data%kappa, n_species)
    do i = 1,n_species

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
    integer :: i_part
    real(kind=c_double), dimension(:), allocatable :: stateVector 

    real(c_double), dimension(:,:), allocatable, target :: state_vector

    ! Get gas array
    stateVecDim = TChem_getLengthOfStateVector()
    nSpec = TChem_getNumberOfSpecies()
    allocate(stateVector(stateVecDim))
    call TChem_getStateVector(stateVector)

    allocate(state_vector(stateVecDim, 1))
    call TChem_getAllStateVectorHost(c_loc(state_vector))
    gas_state%mix_rat = 0.0
    ! FIXME: adjust this range later
    gas_state%mix_rat = stateVector(4:nSpec+3)
    ! Convert gas_state from ppm to ppb.
    call gas_state_scale(gas_state, 1000.d0)
    ! Map aerosols
    do i_part = 1,aero_state_n_part(aero_state)

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

    real(kind=dp), allocatable :: stateVector(:)
    integer :: stateVecDim
    real(kind=dp), parameter :: t_steam = 373.15 ! steam temperature (K)
    real(kind=dp) :: a, water_vp
    integer :: i_part
    integer :: i_water

    ! Get size of stateVector
    stateVecDim = TChem_getLengthOfStateVector()
    allocate(stateVector(stateVecDim))
    ! First three elements are density, pressure and temperature
    stateVector(1) = env_state_air_den(env_state)
    stateVector(2) = env_state%pressure
    stateVector(3) = env_state%temp

    ! PartMC uses relative humidity and not H2O mixing ratio.
    i_water = gas_data_spec_by_name(gas_data, "H2O")
    a = 1.0 - t_steam / env_state%temp
    a = (((-0.1299 * a - 0.6445) * a - 1.976) * a + 13.3185) * a
    water_vp = 101325.0 * exp(a)  ! (Pa)
    gas_state%mix_rat(i_water) = env_state%rel_humid * water_vp * 1.0e9 &
         / env_state%pressure ! (ppb)
    ! Convert from ppb to ppm.
    call gas_state_scale(gas_state, 1.0d0 / 1000.d0)
    ! Add gas species to state vector.
    stateVector(4:gas_data_n_spec(gas_data)+3) = gas_state%mix_rat

    ! TODO: Map aerosols
    do i_part = 1,aero_state_n_part(aero_state)

    end do

    call TChem_setStateVector(stateVector)

  end subroutine tchem_from_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a single timestep of TChem chemistry.
  subroutine TChem_timestep(del_t)
    use iso_c_binding
    real(kind=c_double) :: del_t

    call TChem_doTimestep(del_t)

  end subroutine TChem_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize TChem.
  subroutine TChem_initialize(chemFile, aeroFile, NumericsFile)
    use iso_c_binding

    !> Chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: chemFile
    !> Chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: aeroFile
    !> Chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: numericsFile

    call initialize(chemFile//c_null_char, aeroFile//c_null_char, &
         numericsFile//c_null_char)

  end subroutine TChem_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get species name from TChem for a given index.
  function TChem_species_name(i_spec) result(species_name)
    use iso_c_binding
    ! Species name.
    character(:), allocatable :: species_name

    integer(kind=c_int), intent(in) :: i_spec
    character(kind=c_char, len=:), allocatable :: cbuf
    integer(kind=c_size_t) :: N

    allocate(character(256) :: cbuf)
    N = len(cbuf)
    N = TChem_getSpeciesName(i_spec, cbuf, N)
    allocate(character(N) :: species_name)
    species_name = cbuf(:N)

  end function TChem_species_name
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_tchem_interface
