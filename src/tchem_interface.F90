! Copyright (C) 2017-2018 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_camp_interface module.

!> An interface between PartMC and the CAMP
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
  subroutine initialize(arg_chemfile) bind(c, name="initialize")
    use iso_c_binding
    character(kind=c_char) :: arg_chemfile(*)
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
  subroutine TChem_doTimestep() bind(C, name="TChem_doTimestep")
  end subroutine
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_TCHEM
  !> Run the CAMP module for the current PartMC state
  subroutine pmc_tchem_interface_solve(env_state, aero_data, aero_state, &
       gas_data, gas_state)

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

    call tchem_from_partmc(gas_data, gas_state, env_state)
    call tchem_timestep()
    call tchem_to_partmc(gas_data, gas_state, env_state)

  end subroutine pmc_tchem_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize TChem and PartMC gas and aerosol data.
  subroutine pmc_tchem_initialize(config_filename, gas_data, aero_data)
    use iso_c_binding

    !>
    character(len=*), intent(in) :: config_filename
    !> Gas data.
    type(gas_data_t), intent(inout) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data

    integer(kind=c_int) :: nSpec
    integer :: i
    real(kind=c_double), dimension(:), allocatable :: array 
    character(:), allocatable ::  val

    ! initialize the model
    call TChem_initialize(trim(config_filename))

    ! what we need:
    !  - number of gas species for gas_data
    !  - name of gas species to set gas_data
    nSpec = TChem_getNumberOfSpecies() 
    call ensure_string_array_size(gas_data%name, nSpec)

    ! name of gas species
    do i = 1,nSpec
       val = TChem_species_name(i-1) 
       gas_data%name(i) = trim(val)
    end do   

    ! We need this
    allocate(gas_data%mosaic_index(gas_data_n_spec(gas_data)))
    gas_data%mosaic_index(:) = 0
 
    print*, 'in partmc', nSpec, trim(config_filename), gas_data_n_spec(gas_data)
  
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clean up TChem.
  subroutine pmc_tchem_cleanup()

    call finalize()

  end subroutine pmc_tchem_cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map all data TChem -> PartMC.
  subroutine tchem_to_partmc(gas_data, gas_state, env_state)

    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    integer(c_int) :: nSpec, stateVecDim
    real(kind=c_double), dimension(:), allocatable :: stateVector 

    ! Get gas array
    stateVecDim = TChem_getLengthOfStateVector()
    nSpec = TChem_getNumberOfSpecies()
    allocate(stateVector(stateVecDim))
    array = 0.0d0
    call TChem_getStateVector(stateVector)

    ! FIXME: adjust this range later
    gas_state%mix_rat = stateVector(4:nSpec+3)
    ! Convert gas_state from mol m^-3 to ppb.
    call gas_state_mole_dens_to_ppb(gas_state, env_state)

  end subroutine tchem_to_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map all data PartMC -> TChem.
  subroutine tchem_from_partmc(gas_data, gas_state, env_state)

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

    ! Make a state vector
    ! Density
    ! Pressure
    ! Temperature
    ! Concentrations
    ! get size of stateVector
    stateVecDim = TChem_getLengthOfStateVector()
    allocate(stateVector(stateVecDim))
    stateVector(1) = env_state_air_den(env_state)
    stateVector(2) = env_state%pressure
    stateVector(3) = env_state%temp

    ! FIXME: we have to do something about water here
    i_water = gas_data_spec_by_name(gas_data, "H2O")
    a = 1.0 - t_steam / env_state%temp
    a = (((-0.1299 * a - 0.6445) * a - 1.976) * a + 13.3185) * a
    water_vp = 101325.0 * exp(a)  ! (Pa)
    gas_state%mix_rat(i_water) = env_state%rel_humid * water_vp * 1.0e9 &
         / env_state%pressure ! (ppb)

    call gas_state_ppb_to_mole_dens(gas_state, env_state)
    stateVector(4:gas_data_n_spec(gas_data)+3) = gas_state%mix_rat

    call TChem_setStateVector(stateVector)

  end subroutine tchem_from_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a single timestep of TChem chemistry.
  subroutine TChem_timestep()
    use iso_c_binding

    call TChem_doTimestep()

  end subroutine TChem_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize TChem.
  subroutine TChem_initialize(chemFile)
    use iso_c_binding

    !> Chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: chemFile

    print*, 'in TChemDriver: ', chemFile

    call initialize(chemFile//c_null_char)

  end subroutine TChem_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get i_spec species name from TChem.
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
