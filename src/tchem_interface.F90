! Copyright (C) 2024-2025 Jeff Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_tchem_interface module.

!> An interface between PartMC and TChem
module pmc_tchem_interface

  use pmc_aero_data
  use pmc_aero_particle
  use pmc_aero_state
  use pmc_aero_binned
  use pmc_bin_grid
  use pmc_constants
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env_state
#ifdef PMC_USE_TCHEM
  use iso_c_binding
#endif
  use pmc_util, only : die_msg, warn_assert_msg, assert_msg

  integer, parameter :: STATE_VEC_ENV_OFFSET = 3  ! density, pressure, temp
  real(kind=dp), parameter :: PPM_TO_PPB = 1000.0d0
  integer, parameter :: DEFAULT_BATCH_INDEX = 0

#ifdef PMC_USE_TCHEM
  interface
    subroutine tchem_c_initialize(arg_chemfile, arg_aerofile, &
         arg_numericsfile, n_batch) &
         bind(c, name="initialize")
      use iso_c_binding
      character(kind=c_char), intent(in) :: arg_chemfile(*)
      character(kind=c_char), intent(in) :: arg_aerofile(*)
      character(kind=c_char), intent(in) :: arg_numericsfile(*)
      integer(kind=c_int), intent(in), value :: n_batch
    end subroutine tchem_c_initialize
    subroutine tchem_c_finalize() bind(c, name="finalize")
    end subroutine tchem_c_finalize
    function TChem_getNumberOfSpecies() bind(c, &
         name="TChem_getNumberOfSpecies")
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
      integer(kind=c_int), value :: i_batch
    end subroutine
    subroutine TChem_setStateVector(array, i_batch) bind(c, &
         name="TChem_setStateVector")
      use iso_c_binding
      real(kind=c_double) :: array(*)
      integer(kind=c_int), value :: i_batch
    end subroutine
    function TChem_getNumberConcentrationVectorSize() bind(c, &
         name="TChem_getNumberConcentrationVectorSize")
      use iso_c_binding
      integer(kind=c_int) :: TChem_getNumberConcentrationVectorSize
    end function
    subroutine TChem_setNumberConcentrationVector(array, i_batch) bind(c, &
         name="TChem_setNumberConcentrationVector")
      use iso_c_binding
      real(kind=c_double) :: array(*)
      integer(kind=c_int), value :: i_batch
    end subroutine
    subroutine TChem_setNParticlesTrack(n) bind(c, &
         name="TChem_setNParticlesTrack")
      use iso_c_binding
      integer(kind=c_int), value :: n
    end subroutine
    integer(kind=c_size_t) function TChem_getSpeciesName(index, result, &
         buffer_size) bind(c, name="TChem_getSpeciesName")
      use iso_c_binding
      integer(kind=c_int), intent(in) :: index
      character(kind=c_char), intent(out) :: result(*)
      integer(kind=c_size_t), intent(in), value :: buffer_size
    end function
    integer(kind=c_size_t) function TChem_getAerosolSpeciesName(index, &
         result, buffer_size) bind(c, name="TChem_getAerosolSpeciesName")
      use iso_c_binding
      integer(kind=c_int), intent(in) :: index
      character(kind=c_char), intent(out) :: result(*)
      integer(kind=c_size_t), intent(in), value :: buffer_size
    end function
    subroutine TChem_doTimestep(del_t) bind(c, name="TChem_doTimestep")
      use iso_c_binding
      real(kind=c_double), intent(in) :: del_t
    end subroutine
  end interface
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_TCHEM
  !> Run chemistry using TChem for the current PartMC state.
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

  !> Run chemistry using TChem for the current PartMC sectional (binned) state.
  !!
  !! Each bin is mapped to a TChem "particle" (the bin's mean particle),
  !! mirroring pmc_tchem_interface_solve() with bins in place of particles.
  !! TChem holds the per-bin number concentration fixed and only changes the
  !! per-species masses, so on return only aero_binned%%vol_conc is updated and
  !! aero_binned%%num_conc is unchanged.
  subroutine pmc_tchem_interface_solve_sect(env_state, aero_data, bin_grid, &
       aero_binned, gas_data, gas_state, del_t)

    !> Environment data.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Binned aerosol state.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time step (s).
    real(kind=dp), intent(in) :: del_t

    type(aero_binned_t) :: aero_binned_pregrow
    integer, save :: redist_method = -1
    character(len=64) :: env_val
    integer :: env_status

    ! Select the sectional bin-remap method once (default moving-center; set
    ! PARTMC_SECTIONAL_REDISTRIBUTE=two_moment for the lower-diffusion
    ! linear-discrete scheme).
    if (redist_method < 0) then
       call get_environment_variable("PARTMC_SECTIONAL_REDISTRIBUTE", &
            env_val, status=env_status)
       if ((env_status == 0) .and. ((trim(env_val) == "two_moment") &
            .or. (trim(env_val) == "2"))) then
          redist_method = AERO_BINNED_REDIST_TWO_MOMENT
          write(*,'(a)') 'tchem sectional bin remap: TWO-MOMENT'
       else
          redist_method = AERO_BINNED_REDIST_MOVING_CENTER
          write(*,'(a)') 'tchem sectional bin remap: MOVING-CENTER'
       end if
    end if

    ! The two-moment remap reconstructs the sub-bin distribution from the
    ! pre-growth state, so capture it before TChem changes the masses.
    if (redist_method == AERO_BINNED_REDIST_TWO_MOMENT) then
       aero_binned_pregrow = aero_binned
    end if

    call tchem_from_partmc_sect(aero_data, bin_grid, aero_binned, gas_data, &
         gas_state, env_state)

    call tchem_timestep(del_t)

    call tchem_to_partmc_sect(aero_data, aero_binned, gas_data, gas_state, &
         env_state)

    ! Remap bins whose mean particle volume has grown (or shrunk) past their
    ! grid edges back onto the fixed bin grid.
    if (redist_method == AERO_BINNED_REDIST_TWO_MOMENT) then
       call aero_binned_redistribute_two_moment(aero_binned, &
            aero_binned_pregrow, bin_grid, aero_data)
    else
       call aero_binned_redistribute_moving_center(aero_binned, bin_grid, &
            aero_data)
    end if

  end subroutine pmc_tchem_interface_solve_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map PartMC binned aerosol and gas state into TChem, treating each bin as
  !> a single TChem particle (the bin's number-mean particle).
  subroutine tchem_from_partmc_sect(aero_data, bin_grid, aero_binned, &
       gas_data, gas_state, env_state)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Binned aerosol state.
    type(aero_binned_t), intent(in) :: aero_binned
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp), allocatable :: state_vector(:), number_concentration(:)
    integer :: state_vec_dim, tchem_n_part
    integer :: i_bin, n_bin, i_spec, i_water
    integer :: n_gas_spec, n_aero_spec, aero_offset
    real(kind=dp), parameter :: SMALL_MASS_PLACEHOLDER = 1.0d-30

    n_gas_spec = gas_data_n_spec(gas_data)
    n_aero_spec = aero_data_n_spec(aero_data)
    n_bin = size(aero_binned%num_conc)

    ! Get size of state vector in TChem
    state_vec_dim = TChem_getLengthOfStateVector()
    allocate(state_vector(state_vec_dim))
    state_vector = 0.0d0

    ! Get size of number concentration vector in TChem
    tchem_n_part = TChem_getNumberConcentrationVectorSize()
    allocate(number_concentration(tchem_n_part))

    call assert_msg(502348107, n_bin <= tchem_n_part, &
         "number of sectional bins exceeds the TChem particle capacity")

    ! First three elements are density, pressure and temperature
    state_vector(1) = env_state_air_den(env_state)
    state_vector(2) = env_state%pressure
    state_vector(3) = env_state%temp

    ! PartMC uses relative humidity and not H2O mixing ratio.
    i_water = gas_data_spec_by_name(gas_data, "H2O")
    gas_state%mix_rat(i_water) = env_state_rel_humid_to_mix_rat(env_state)
    ! Add gas species to state vector. Convert from ppb to ppm.
    state_vector(STATE_VEC_ENV_OFFSET+1:n_gas_spec + STATE_VEC_ENV_OFFSET) = &
         gas_state%mix_rat / PPM_TO_PPB

    ! Each bin becomes one TChem particle. aero_binned%num_conc and %vol_conc
    ! are per-log-width densities (dN/dlnD, dV/dlnD), so the per-particle mass
    ! is vol_conc/num_conc (the dlnD cancels), while the actual number
    ! concentration (#/m^3) passed to TChem is num_conc * the bin log-width.
    aero_offset = n_gas_spec + STATE_VEC_ENV_OFFSET
    do i_bin = 1,n_bin
       if (aero_binned%num_conc(i_bin) > 0.0d0) then
          do i_spec = 1,n_aero_spec
             state_vector(aero_offset + i_spec + (i_bin - 1) * n_aero_spec) = &
                  aero_binned%vol_conc(i_bin, i_spec) &
                  * aero_data%density(i_spec) / aero_binned%num_conc(i_bin)
          end do
          number_concentration(i_bin) = aero_binned%num_conc(i_bin) &
               * bin_grid%widths(i_bin)
       else
          do i_spec = 1,n_aero_spec
             state_vector(aero_offset + i_spec + (i_bin - 1) * n_aero_spec) = &
                  SMALL_MASS_PLACEHOLDER
          end do
          number_concentration(i_bin) = 0.0d0
       end if
    end do

    ! Pad any unused TChem particle slots.
    do i_bin = n_bin+1,tchem_n_part
       do i_spec = 1,n_aero_spec
          state_vector(aero_offset + i_spec + (i_bin - 1) * n_aero_spec) = &
               SMALL_MASS_PLACEHOLDER
       end do
       number_concentration(i_bin) = 0.0d0
    end do

    call TChem_setStateVector(state_vector, DEFAULT_BATCH_INDEX)

    call TChem_setNumberConcentrationVector(number_concentration, &
         DEFAULT_BATCH_INDEX)

    call TChem_setNParticlesTrack(n_bin)

  end subroutine tchem_from_partmc_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Map TChem results back into the PartMC binned state.
  !!
  !! Only the per-bin per-species masses are updated (into vol_conc); the bin
  !! number concentrations are held fixed by TChem and are left unchanged.
  subroutine tchem_to_partmc_sect(aero_data, aero_binned, gas_data, &
       gas_state, env_state)

    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Binned aerosol state.
    type(aero_binned_t), intent(inout) :: aero_binned
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp), allocatable :: state_vector(:)
    integer :: state_vec_dim, i_bin, n_bin, i_spec
    integer :: n_gas_spec, n_aero_spec, aero_offset

    n_gas_spec = gas_data_n_spec(gas_data)
    n_aero_spec = aero_data_n_spec(aero_data)
    n_bin = size(aero_binned%num_conc)

    state_vec_dim = TChem_getLengthOfStateVector()
    allocate(state_vector(state_vec_dim))
    call TChem_getStateVector(state_vector, DEFAULT_BATCH_INDEX)

    ! Convert gas from ppm to ppb.
    gas_state%mix_rat = &
         state_vector(STATE_VEC_ENV_OFFSET+1:n_gas_spec+STATE_VEC_ENV_OFFSET) &
         * PPM_TO_PPB

    ! Per-bin mean-particle masses back to per-bin volume concentration,
    ! holding the bin number concentration fixed (so empty bins stay empty).
    aero_offset = n_gas_spec + STATE_VEC_ENV_OFFSET
    do i_bin = 1,n_bin
       do i_spec = 1,n_aero_spec
          aero_binned%vol_conc(i_bin, i_spec) = &
               state_vector(aero_offset + i_spec + (i_bin - 1) * n_aero_spec) &
               * aero_binned%num_conc(i_bin) / aero_data%density(i_spec)
       end do
    end do

  end subroutine tchem_to_partmc_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize TChem and PartMC gas and aerosol data.
  subroutine pmc_tchem_initialize(gas_config_filename, aero_config_filename, &
       solver_config_filename, gas_data, aero_data, n_grid_cells)

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
    character(:), allocatable :: val
    integer :: i_spec
    logical :: is_gas

    ! initialize the model
    call tchem_c_initialize(trim(gas_config_filename), &
         trim(aero_config_filename), trim(solver_config_filename), &
         n_grid_cells)

    ! Get size of gas_data from TChem
    n_gas_spec = TChem_getNumberOfSpecies()
    call ensure_string_array_size(gas_data%name, n_gas_spec)

    ! Set gas_data with gas species from TChem
    is_gas = .true.
    do i_spec = 1,n_gas_spec
       val = tchem_species_name(i_spec - 1, is_gas)
       gas_data%name(i_spec) = trim(val)
    end do

    ! For output and MPI, this needs to be allocated (for now)
    allocate(gas_data%mosaic_index(n_gas_spec))
    gas_data%mosaic_index(:) = 0

    ! Get size of aero_data from TChem
    n_aero_spec = TChem_getNumberOfAeroSpecies()
    call ensure_string_array_size(aero_data%name, n_aero_spec)
    call ensure_integer_array_size(aero_data%mosaic_index, n_aero_spec)
    call ensure_real_array_size(aero_data%wavelengths, n_swbands)
    call ensure_real_array_size(aero_data%density, n_aero_spec)
    call ensure_integer_array_size(aero_data%num_ions, n_aero_spec)
    call ensure_real_array_size(aero_data%molec_weight, n_aero_spec)
    call ensure_real_array_size(aero_data%kappa, n_aero_spec)
    call ensure_real_array_size(aero_data%abifm_m, n_aero_spec)
    call ensure_real_array_size(aero_data%abifm_c, n_aero_spec)

    is_gas = .false.
    do i_spec = 1,n_aero_spec
       val = tchem_species_name(i_spec - 1, is_gas)
       aero_data%name(i_spec) = trim(val)
       aero_data%density(i_spec) = TChem_getAerosolSpeciesDensity(i_spec - 1)
       aero_data%molec_weight(i_spec) = TChem_getAerosolSpeciesMW(i_spec - 1)
       aero_data%kappa(i_spec) = TChem_getAerosolSpeciesKappa(i_spec - 1)
       aero_data%abifm_m(i_spec) = 0.0d0
       aero_data%abifm_c(i_spec) = 0.0d0
    end do

  end subroutine pmc_tchem_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clean up TChem.
  subroutine pmc_tchem_cleanup()

    call tchem_c_finalize()

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

    integer(kind=c_int) :: state_vec_dim
    integer :: i_part, i_spec
    real(kind=c_double), allocatable :: state_vector(:)
    integer :: n_gas_spec, n_aero_spec
    real(kind=dp) :: reweight_num_conc(aero_state_n_part(aero_state))
    integer :: aero_offset

    n_gas_spec = gas_data_n_spec(gas_data)
    n_aero_spec = aero_data_n_spec(aero_data)

    ! Get gas array
    state_vec_dim = TChem_getLengthOfStateVector()
    allocate(state_vector(state_vec_dim))
    call TChem_getStateVector(state_vector, DEFAULT_BATCH_INDEX)

    ! Convert from ppm to ppb.
    gas_state%mix_rat = &
       state_vector(STATE_VEC_ENV_OFFSET+1:n_gas_spec+STATE_VEC_ENV_OFFSET) &
       * PPM_TO_PPB

    call aero_state_num_conc_for_reweight(aero_state, aero_data, &
         reweight_num_conc)

    aero_offset = n_gas_spec + STATE_VEC_ENV_OFFSET
    do i_part = 1,aero_state_n_part(aero_state)
       do i_spec = 1,n_aero_spec
          aero_state%apa%particle(i_part)%vol(i_spec) = state_vector( &
               aero_offset + i_spec + (i_part - 1) * n_aero_spec) &
               / aero_data%density(i_spec)
       end do
    end do

    call aero_state_reweight(aero_state, aero_data, reweight_num_conc)

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
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state

    real(kind=dp), allocatable :: state_vector(:), number_concentration(:)
    integer :: state_vec_dim, tchem_n_part, i_spec
    integer :: i_part, n_part
    integer :: i_water
    integer :: n_gas_spec, n_aero_spec
    integer :: aero_offset
    real(kind=dp), parameter :: SMALL_MASS_PLACEHOLDER = 1.0d-30

    n_gas_spec = gas_data_n_spec(gas_data)
    n_aero_spec = aero_data_n_spec(aero_data)
    n_part = aero_state_n_part(aero_state)

    ! Get size of state vector in TChem
    state_vec_dim = TChem_getLengthOfStateVector()
    allocate(state_vector(state_vec_dim))
    state_vector = 0.0d0

    ! Get size of number concentration vector in TChem
    tchem_n_part = TChem_getNumberConcentrationVectorSize()
    allocate(number_concentration(tchem_n_part))

    ! First three elements are density, pressure and temperature
    state_vector(1) = env_state_air_den(env_state)
    state_vector(2) = env_state%pressure
    state_vector(3) = env_state%temp

    ! PartMC uses relative humidity and not H2O mixing ratio.
    i_water = gas_data_spec_by_name(gas_data, "H2O")
    gas_state%mix_rat(i_water) = env_state_rel_humid_to_mix_rat(env_state)
    ! Add gas species to state vector. Convert from ppb to ppm.
    state_vector(STATE_VEC_ENV_OFFSET+1:n_gas_spec + STATE_VEC_ENV_OFFSET) = &
         gas_state%mix_rat / PPM_TO_PPB

    aero_offset = n_gas_spec + STATE_VEC_ENV_OFFSET
    do i_part = 1,n_part
       do i_spec = 1,n_aero_spec
          state_vector(aero_offset + i_spec + (i_part - 1) * n_aero_spec) = &
               aero_particle_species_mass(aero_state%apa%particle(i_part), &
               i_spec, aero_data)
       end do
       number_concentration(i_part) = aero_state_particle_num_conc( &
            aero_state, aero_state%apa%particle(i_part), aero_data)
    end do

    do i_part = n_part+1,tchem_n_part
       do i_spec = 1,n_aero_spec
          state_vector(aero_offset + i_spec + (i_part-1) * n_aero_spec) = &
               SMALL_MASS_PLACEHOLDER
       end do
       number_concentration(i_part) = 0.0d0
    end do

    call TChem_setStateVector(state_vector, DEFAULT_BATCH_INDEX)

    call TChem_setNumberConcentrationVector(number_concentration, &
         DEFAULT_BATCH_INDEX)

    call TChem_setNParticlesTrack(n_part)

  end subroutine tchem_from_partmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a single timestep of TChem chemistry.
  subroutine tchem_timestep(del_t)

    !> Time step (s).
    real(kind=c_double), intent(in) :: del_t

    call TChem_doTimestep(del_t)

  end subroutine tchem_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize TChem.
  subroutine tchem_initialize(chem_file, aero_file, numerics_file, n_batch)

    !> Gas chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: chem_file
    !> Aerosol chemistry configuration file.
    character(kind=c_char,len=*), intent(in) :: aero_file
    !> Numerics configuration file.
    character(kind=c_char,len=*), intent(in) :: numerics_file
    !> Number of systems to solve.
    integer(kind=c_int), intent(in) :: n_batch

    call tchem_c_initialize(chem_file//c_null_char, aero_file//c_null_char, &
         numerics_file//c_null_char, n_batch)

  end subroutine tchem_initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get species name from TChem for a given index.
  function tchem_species_name(i_spec, is_gas) result(species_name)

    !> Index of species.
    integer(kind=c_int), intent(in) :: i_spec
    !> Logical for if the species is a gas.
    logical, intent(in) :: is_gas
    !> Species name.
    character(:), allocatable :: species_name

    character(kind=c_char, len=:), allocatable :: cbuf
    integer(kind=c_size_t) :: N

    if (is_gas) then
       N = GAS_NAME_LEN
    else
       N = AERO_NAME_LEN
    end if
    allocate(character(N) :: cbuf)
    if (is_gas) then
       N = TChem_getSpeciesName(i_spec, cbuf, N)
    else
       N = TChem_getAerosolSpeciesName(i_spec, cbuf, N)
    end if
    allocate(character(N) :: species_name)
    species_name = cbuf(:N)

  end function tchem_species_name
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_tchem_interface
