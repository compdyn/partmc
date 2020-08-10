! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_cb05cl_ae5 program

!> Test for the cb05cl_ae5 mechanism from MONARCH. This program runs the
!! MONARCH CB5 code and the CAMP-chem version and compares the output.


!todo: regroup cells by stiffness in mpi cell division
!todo improve test to handle all cases

program pmc_test_cb05cl_ae5

  use pmc_constants,                    only: const
  use pmc_util,                         only: i_kind, dp, assert, assert_msg, &
                                              almost_equal, string_t, &
                                              to_string, warn_assert_msg
  use pmc_mpi
  use pmc_camp_core
  use pmc_camp_state
  use pmc_camp_solver_data
  use pmc_solver_stats
  use pmc_chem_spec_data
  use pmc_mechanism_data
  use pmc_rxn_data
  use pmc_rxn_photolysis
  use pmc_rxn_factory
  use pmc_property
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_netcdf

#define PMC_MONARCH_INPUT

  ! EBI Solver
  use module_bsc_chem_data

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)
  ! EBI solver output file unit
  integer(kind=i_kind), parameter :: EBI_FILE_UNIT = 10
  ! KPP solver output file unit
  integer(kind=i_kind), parameter :: KPP_FILE_UNIT = 11
  ! CAMP-chem output file unit
  integer(kind=i_kind), parameter :: CAMP_FILE_UNIT = 12
  ! CAMP-chem output profiling stats file unit
  integer(kind=i_kind), parameter :: CAMP_FILE_UNIT_PROFILE = 13
#ifdef PMC_MONARCH_INPUT
  ! EBI solver output file unit
  integer(kind=i_kind), parameter :: CAMP_EBI_FILE_UNIT = 14
  ! file unit
  integer(kind=i_kind), parameter :: CAMP_KPP_FILE_UNIT = 15
  ! file unit
  integer(kind=i_kind), parameter :: EBI_KPP_FILE_UNIT = 16
#endif
  ! Number of timesteps to integrate over
  integer(kind=i_kind), parameter :: NUM_TIME_STEPS = 1
  ! Number of cells
  integer(kind=i_kind), parameter :: NUM_CELLS= 10800
  !100
  !1125
  !3375
  !5625
  !7875
  !10800
  !100000
  ! Number of EBI-solver species
  integer(kind=i_kind), parameter :: NUM_EBI_SPEC = 72
  ! Number of EBI-solever photolysis reactions
  integer(kind=i_kind), parameter :: NUM_EBI_PHOTO_RXN = 23
  ! Small number for minimum concentrations
  real(kind=dp), parameter :: SMALL_NUM = 1.0d-30
  ! Used to check availability of a solver
  type(camp_solver_data_t), pointer :: camp_solver_data

  !Command arguments mapping
  integer(kind=i_kind), parameter :: MAP_I_START = 1
  integer(kind=i_kind), parameter :: MAP_J_START = 2
  integer(kind=i_kind), parameter :: MAP_K_START = 3
  integer(kind=i_kind), parameter :: MAP_T_START = 4
  integer(kind=i_kind), parameter :: MAP_I_N = 5
  integer(kind=i_kind), parameter :: MAP_J_N = 6
  integer(kind=i_kind), parameter :: MAP_K_N = 7
  integer(kind=i_kind), parameter :: MAP_T_N = 8

  !integer(kind=i_kind), parameter :: N_COMMAND_ARGUMENTS = 9

#ifdef DEBUG
  integer(kind=i_kind), parameter :: DEBUG_UNIT = 13

  open(unit=DEBUG_UNIT, file="out/debug_cb05cl_ae.txt", status="replace", action="write")
#endif

  call pmc_mpi_init()

  camp_solver_data => camp_solver_data_t()

  if (.not.camp_solver_data%is_solver_available()) then
    write(*,*) "CB5 mechanism test - no solver available - PASS"
  else if (run_cb05cl_ae5_tests()) then
    write(*,*) "Finish test_cb05cl_ae5_big"
  else
    write(*,*) "CB5 mechanism tests - FAIL"
  end if

  deallocate(camp_solver_data)

#ifdef DEBUG
  close(DEBUG_UNIT)
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all CB5 tests
  logical function run_cb05cl_ae5_tests() result(passed)

    passed = run_standard_cb05cl_ae5_test()

  end function run_cb05cl_ae5_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run the cb05cl_ae5 mechanism under standard conditions using the original
  !! MONARCH ebi-solver code, the KPP CB5 module and the CAMP-chem version
  logical function run_standard_cb05cl_ae5_test() result(passed)

    ! EBI Solver
    use EXT_HRDATA
    use EXT_RXCM,                               only : NRXNS, RXLABEL

    ! KPP Solver
    use cb05cl_ae5_Initialize,                  only : KPP_Initialize => Initialize
    use cb05cl_ae5_Model,                       only : KPP_NSPEC => NSPEC, &
                                                       KPP_STEPMIN => STEPMIN, &
                                                       KPP_STEPMAX => STEPMAX, &
                                                       KPP_RTOL => RTOL, &
                                                       KPP_ATOL => ATOL, &
                                                       KPP_TIME => TIME, &
                                                       KPP_C => C, &
                                                       KPP_RCONST => RCONST, &
                                                       KPP_Update_RCONST => Update_RCONST, &
                                                       KPP_INTEGRATE => INTEGRATE, &
                                                       KPP_SPC_NAMES => SPC_NAMES, &
                                                       KPP_PHOTO_RATES => PHOTO_RATES, &
                                                       KPP_TEMP => TEMP, &
                                                       KPP_PRESS => PRESS, &
                                                       KPP_SUN => SUN, &
                                                       KPP_M => M, &
                                                       KPP_N2 => N2, &
                                                       KPP_O2 => O2, &
                                                       KPP_H2 => H2, &
                                                       KPP_H2O => H2O, &
                                                       KPP_N2O => N2O, &
                                                       KPP_CH4 => CH4, &
                                                       KPP_NVAR => NVAR, &
                                                       KPP_NREACT => NREACT, &
                                                       KPP_DT => DT
    use cb05cl_ae5_Parameters,                  only : KPP_IND_O2 => IND_O2
    use cb05cl_ae5_Initialize, ONLY: Initialize

    ! EBI-solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: ebi_spec_names
#ifdef PMC_MONARCH_INPUT
    ! EBI-solver species names in MONARCH order
    type(string_t), dimension(NUM_EBI_SPEC) :: ebi_monarch_spec_names
#endif

    ! KPP reaction labels
    type(string_t), allocatable :: kpp_rxn_labels(:)
    ! KPP rstate
    real(kind=dp) :: KPP_RSTATE(20)
    ! KPP control variables
    integer :: KPP_ICNTRL(20)
    ! #/cc -> ppm conversion factor
    real(kind=dp) :: conv

    ! Flag for sunlight
    logical :: is_sunny
    ! Photolysis rates (\min)
    real, allocatable :: photo_rates(:)
#ifdef PMC_MONARCH_INPUT
    real, allocatable :: new_rates(:)
#endif
    ! Temperature (K)
    real :: temperature
    ! Pressure (atm)
    real :: pressure
    ! Water vapor concentration (ppmV)
    real :: water_conc

    ! CAMP-chem core
    type(camp_core_t), pointer :: camp_core
    ! CAMP-chem state
    type(camp_state_t), pointer :: camp_state, camp_state_comp
    ! CAMP-chem species names
    type(string_t), allocatable :: camp_spec_names(:)
    ! EBI -> CAMP-chem species map
    integer(kind=i_kind), dimension(NUM_EBI_SPEC) :: spec_map

    ! Computation timer variables
    real(kind=dp) :: comp_start, comp_end, comp_ebi, comp_kpp, comp_camp

    type(chem_spec_data_t), pointer :: chem_spec_data
    class(rxn_data_t), pointer :: rxn
    type(property_t), pointer :: prop_set
    character(len=:), allocatable :: key, spec_name, string_val, camp_input_file
#ifdef PMC_MONARCH_INPUT
    character(500) :: name_cell
#endif
    real(kind=dp) :: real_val, camp_rate, camp_rate_const
    integer(kind=i_kind) :: i_spec, j_spec, i_rxn, i_ebi_rxn, i_kpp_rxn, &
            i_time, i_repeat, n_gas_spec

    integer(kind=i_kind) :: i_M, i_O2, i_N2, i_H2O, i_CH4, i_H2, i_DUMMY
    integer(kind=i_kind), allocatable :: ebi_rxn_map(:), kpp_rxn_map(:)
    integer(kind=i_kind), allocatable :: ebi_spec_map(:), kpp_spec_map(:)
    type(string_t) :: str_temp
    type(string_t), allocatable :: spec_names(:)
    type(solver_stats_t), target :: solver_stats

    ! Pointer to the mechanism
    type(mechanism_data_t), pointer :: mechanism

    ! Variables to set photolysis rates
    type(rxn_factory_t) :: rxn_factory
    integer(kind=i_kind) :: n_photo_rxn, i_photo_rxn
    type(rxn_update_data_photolysis_t), allocatable :: rate_update(:)
    type(rxn_update_data_photolysis_t) :: jo2_rate_update

    ! Arrays to hold starting concentrations
    real(kind=dp), allocatable :: ebi_init(:), kpp_init(:), camp_init(:)
    real(kind=dp), allocatable :: temperatures(:), pressures(:)
#ifdef PMC_MONARCH_INPUT
    real(kind=dp), dimension(NUM_EBI_SPEC) :: ebi_monarch_init
    integer, dimension(NUM_EBI_SPEC) :: map_ebi_monarch
#endif
    real(kind=dp), allocatable :: model_conc(:)
    integer(kind=i_kind) :: n_cells, n_blocks, n_cells_block, compare_results, i, j, k, s, state_size_cell, i_cell, i_block
    real(kind=dp) :: offset_conc
    real(kind=dp) :: offset_temp
    integer(kind=i_kind) :: n_repeats
    !netcdf
    integer(kind=i_kind) :: input_id, varid, offset
    character*8 :: t
    character(len=50) :: aux_arg
    integer :: status_code
    integer :: pmc_multicells

#ifdef PMC_USE_MPI
    character, allocatable :: buffer(:), buffer_copy(:)
    integer(kind=i_kind) :: pack_size, pos, i_elem, results, mpi_threads
#endif

    ! initialize MPI
     !call pmc_mpi_init()

    ! Read from the .sh script the command arguments
    call get_command_argument(1, aux_arg, status=status_code)
    read(aux_arg,*)n_cells_block
    call get_command_argument(2, aux_arg, status=status_code)
    read(aux_arg,*)n_blocks
    call get_command_argument(3, aux_arg, status=status_code)
    read(aux_arg,*)offset_conc
    call get_command_argument(4, aux_arg, status=status_code)
    read(aux_arg,*)offset_temp
    call get_command_argument(5, aux_arg, status=status_code)
    read(aux_arg,*)pmc_multicells

    n_cells=n_cells_block*n_blocks
    !n_blocks=n_cells/n_cells_block
    write(*,*) "n_cells", n_cells, "n_cells_block", n_cells_block, "n_blocks", n_blocks

    call assert_msg(921735481, (pmc_multicells.ne.0 .or. pmc_multicells.ne.1), "Wrong pmc_multicells config value (use 0 or 1)")

    !Default init
    !n_cells = NUM_CELLS !i_n*j_n*k_n
    KPP_ICNTRL( : ) = 0
    temperature = 297.93 !v9:202.9565 !v48:297.93 !orig:272.5
    pressure = 0.998 !v9:0.1456779 !v48:0.998 !orig:0.8
    water_conc = 0.0 ! (Set by CAMP-chem initial concentration)

    compare_results = 1
    n_repeats = 1!2000

    passed = .false.

    ! Set the #/cc -> ppm conversion factor
    conv = 1.0d0/ (const%avagadro /const%univ_gas_const * 10.0d0**(-12.0d0) * &
            (pressure*101325.d0) /temperature)

    ! Load the EBI solver species names
    call set_ebi_species(ebi_spec_names)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize the EBI solver !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call cpu_time(comp_start)
    ! Set the BSC chem parameters
    call init_bsc_chem_data()
    ! Set the output unit
    LOGDEV = 6
    ! Set the aerosol flag
    L_AE_VRSN = .false.
    ! Set the aq. chem flag
    L_AQ_VRSN = .false.
    ! Initialize the solver
    call EXT_HRINIT
    RKI(:) = 0.0
    RXRAT(:) = 0.0
    YC(:) = 0.0
    YC0(:) = 0.0
    YCP(:) = 0.0
    PROD(:) = 0.0
    LOSS(:) = 0.0
    PNEG(:) = 0.0
    ! Set the timestep (min)
#ifdef PMC_MONARCH_INPUT
    EBI_TMSTEP = 3. !monarch:3 !orig:0.1
#else
    EBI_TMSTEP = 0.1
#endif
    ! Set the number of timesteps
    N_EBI_STEPS = 1
    ! Set the number of internal timesteps
    N_INR_STEPS = 1
    call cpu_time(comp_end)
    write(*,*) "EBI initialization time: ", comp_end-comp_start," s"

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize the KPP CB5 module !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call cpu_time(comp_start)
    ! Set the step limits
    KPP_STEPMIN = 0.0d0
    KPP_STEPMAX = 0.0d0
    KPP_SUN = 1.0
    ! Set the tolerances
    do i_spec = 1, KPP_NVAR
      KPP_RTOL(i_spec) = 1.0d-4
      KPP_ATOL(i_spec) = 1.0d-3
    end do
    CALL KPP_Initialize()
    call cpu_time(comp_end)
    write(*,*) "KPP initialization time: ", comp_end-comp_start," s"

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize camp-chem !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_MPI
    if( pmc_mpi_rank( ) .eq. 0 ) then
#endif
    call cpu_time(comp_start)
    camp_input_file = "config_cb05cl_ae5_big.json"
    if(pmc_multicells.eq.1) then
      camp_core => camp_core_t(camp_input_file, n_cells=n_cells_block)
    else
      camp_core => camp_core_t(camp_input_file)
    endif
    ! Initialize the model
    call camp_core%initialize()

    ! Find the CB5 mechanism
    key = "cb05cl_ae5"
    call assert(418262750, camp_core%get_mechanism(key, mechanism))

    ! Set the photolysis rate ids
    key = "rxn id"
    n_photo_rxn = 0
    do i_rxn = 1, mechanism%size()
      rxn => mechanism%get_rxn(i_rxn)
      select type(rxn)
        type is (rxn_photolysis_t)
          call assert(265614917, rxn%property_set%get_string(key, string_val))
          if (trim(string_val).ne."jo2") then
            n_photo_rxn = n_photo_rxn + 1
          end if
      end select
    end do
    allocate(rate_update(n_photo_rxn))
    i_photo_rxn = 0
    do i_rxn = 1, mechanism%size()
      rxn => mechanism%get_rxn(i_rxn)
      select type(rxn)
        type is (rxn_photolysis_t)
          call assert(265614917, rxn%property_set%get_string(key, string_val))
          if (trim(string_val).eq."jo2") then
            ! Set O2 + hv rate constant to 0 (not present in ebi version)
            call camp_core%initialize_update_object(rxn, jo2_rate_update)
          else
            i_photo_rxn = i_photo_rxn + 1
            call camp_core%initialize_update_object(rxn, &
                                                    rate_update(i_photo_rxn))
          end if
      end select
    end do
    call assert(322300770, n_photo_rxn.eq.i_photo_rxn)

#ifdef PMC_USE_MPI

    ! Pack the cores and the unit test
    pack_size = camp_core%pack_size( )
    allocate( buffer( pack_size ) )
    pos = 0
    call camp_core%bin_pack(  buffer, pos )
    call assert( 881897913, pos .eq. pack_size )

    end if

    ! Broadcast the buffer size
    call pmc_mpi_bcast_integer( pack_size )

    if (pmc_mpi_rank().eq.1) then
      ! allocate the buffer to receive data
      allocate(buffer(pack_size))
    end if

    ! broadcast the data
    call pmc_mpi_bcast_packed(buffer)

    if( pmc_mpi_rank( ) .eq. 1 ) then
      write(*,*) "Rank 1 exists"

      ! unpack the data
      if(pmc_multicells.eq.1) then
        camp_core => camp_core_t(n_cells=n_cells_block)
      else
        camp_core => camp_core_t()
      endif
      !camp_core  => camp_core_t( )
      pos = 0
      call camp_core%bin_unpack(  buffer, pos )

      ! Try repacking the data and making sure it stays the same
      allocate( buffer_copy( pack_size ) )
      pos = 0
      call camp_core%bin_pack(  buffer_copy, pos )
      call assert( 276642139, pos .eq. pack_size )
      do i_elem = 1, pack_size
        call assert_msg( 443440270, buffer( i_elem ) .eq. &
                buffer_copy( i_elem ), &
                "Mismatch in element "//trim( to_string( i_elem ) ) )
      end do
      deallocate( buffer_copy )

    end if

#endif

    ! Initialize the solver
    call camp_core%solver_initialize()

    ! Get an new state variable
    camp_state => camp_core%new_state()

    ! Set the environmental conditions
    !call camp_state%env_states(1)%set_temperature_K( real( temperature, kind=dp ) )
    !call camp_state%env_states(1)%set_pressure_Pa( pressure * const%air_std_press )

    call cpu_time(comp_end)
    write(*,*) "CAMP-chem initialization time: ", comp_end-comp_start," s"

    ! Get a camp-state for rate comparisons
    !camp_state_comp => camp_core%new_state()
    !call camp_state_comp%env_states(1)%set_temperature_K( &
    !  camp_state%env_states(1)%val%temp )
    !call camp_state_comp%env_states(1)%set_pressure_Pa(   &
    !  camp_state%env_states(1)%val%pressure )

    ! Get the chemical species data
    call assert(298481296, camp_core%get_chem_spec_data(chem_spec_data))

    ! Find the constant species in the CB5 mechanism
    !todo this species are not plot, maybe there are on state but not in deriv I guess?
    spec_name = "M"
    i_M   = chem_spec_data%gas_state_id(spec_name)
    spec_name = "O2"
    i_O2  = chem_spec_data%gas_state_id(spec_name)
    spec_name = "N2"
    i_N2  = chem_spec_data%gas_state_id(spec_name)
    spec_name = "H2O"
    i_H2O = chem_spec_data%gas_state_id(spec_name)
    spec_name = "CH4"
    i_CH4 = chem_spec_data%gas_state_id(spec_name)
    spec_name = "H2"
    i_H2  = chem_spec_data%gas_state_id(spec_name)

    ! Set the photolysis rates (dummy values for solver comparison)
    is_sunny = .true.
    allocate(photo_rates(NUM_EBI_PHOTO_RXN))
#ifdef PMC_MONARCH_INPUT
    allocate(new_rates(NUM_EBI_PHOTO_RXN))
#endif

    !v48
    !new_rates = (/0.550871,3.0273698E-2,2.5476560E-3,12.23329,1.575956,0.120056,&
    !        4.4206658E-4,2.7387924E-4,4.3699056E-5,3.0446171E-3,1.7465037E-4,&
    !        3.3178381E-4,1.8682621E-3,2.8112135E-3,3.3792580E-4,5.1333202E-5,&
    !        1.6838829E-5,29445188E-3,5.0825302E-2,10457919E-4,0.0,0.0,0.0/)

    !v9 wrong
    !new_rates = (/0.8685482,3.6514506E-02,4.0941574E-03,13.91660,1.792809,&
    !0.1908291,6.7898177E-04,5.7340757E-04,5.6609770E-05,2.2438637E-03,&
    !1.6916823E-04,6.2022154E-04,3.8634937E-03,5.5694710E-03,7.1730785E-04,&
    !4.1930823E-05,3.3923301E-05,5.8016707E-03,6.8607539E-02,1.7641661E-04,&
    !0.0000000E+00,0.0000000E+00,0.0000000E+00/)

    !v9 good
    !new_rates = (/0.8685482,3.6514506E-02,4.0941574E-03,13.91660,1.792809,&
    !        0.1908291,6.7898177E-04,5.7340757E-04,5.6609770E-05,2.2438637E-03,&
    !        1.6916823E-04,6.2022154E-04,6.2022154E-04,3.8634937E-03,5.5694710E-03,&
    !        7.1730785E-04,4.1930823E-05,3.3923301E-05,5.8016707E-03,4.1930823E-05,&
    !        3.8634937E-03,6.8607539E-02,1.7641661E-04/)

#ifdef PMC_MONARCH_INPUT
    photo_rates(:) = 0.0
    !photo_rates(:) = new_rates(:)
    KPP_PHOTO_RATES(:) = photo_rates(:)/60
#else
    !todo: set photo_rates to 0 for night and X for day
    photo_rates(:) = 0.0001 * 60.0 ! EBI solver wants rates in min^-1
    KPP_PHOTO_RATES(:) = 0.0001
    !KPP_PHOTO_RATES(:) = 0.0
#endif
    !write(*,*) "photo_rates/60", new_rates(:)/60
    ! Set O2 + hv rate constant to 0 in KPP (not present in ebi version)
    KPP_PHOTO_RATES(1) = 0.0
    write(*,*) "n_photo_rxn", n_photo_rxn
    ! Set the remaining rates
    if(pmc_multicells.eq.1) then
      do i_cell = 1, n_cells_block
        ! Set the O2 + hv rate constant to 0 (not present in ebi version)
        call jo2_rate_update%set_rate(real(0.0, kind=dp))
        call camp_core%update_data(jo2_rate_update,i_cell+1)
        do i_photo_rxn = 1, n_photo_rxn
#ifdef PMC_MONARCH_INPUT
          call rate_update(i_photo_rxn)%set_rate(real(photo_rates(i_photo_rxn)/60, kind=dp))
#else
          call rate_update(i_photo_rxn)%set_rate(real(0.0001+0.0000001*(i_cell-1), kind=dp))
#endif
          call camp_core%update_data(rate_update(i_photo_rxn),i_cell)
        end do
      end do
    else
      ! Set the O2 + hv rate constant to 0 (not present in ebi version)
      call jo2_rate_update%set_rate(real(0.0, kind=dp))
      call camp_core%update_data(jo2_rate_update)
      do i_photo_rxn = 1, n_photo_rxn
#ifdef PMC_MONARCH_INPUT
        call rate_update(i_photo_rxn)%set_rate(real(photo_rates(i_photo_rxn)/60, kind=dp))
#else
        call rate_update(i_photo_rxn)%set_rate(real(0.1, kind=dp))
#endif
        call camp_core%update_data(rate_update(i_photo_rxn))
      end do
    end if

    ! Make sure the right number of reactions is present
    ! (KPP includes two Cl rxns with rate constants set to zero that are not
    !  present in camp-chem)
    call assert_msg(396732632, mechanism%size().eq.186, &
            "Wrong number of camp-chem reactions: "// &
                    trim(to_string(mechanism%size())))

    ! Set the initial concentrations
    key = "init conc"
    YC(:) = 0.0
    KPP_C(:) = 0.0
    camp_state%state_var(:) = 0.0

    !todo use all time n_cells *n_blocks is bored, group in one variable

    state_size_cell = size(chem_spec_data%get_spec_names()) !size(camp_state%state_var) / n_cells

    !if(pmc_multicells.eq.0) then

    allocate(model_conc(state_size_cell*n_cells))
    model_conc(:)=0.0
    !endif

    allocate(temperatures(n_cells))
    temperatures(:) = temperature
    allocate(pressures(n_cells))
    pressures(:) = pressure*const%air_std_press

#ifdef PMC_MONARCH_INPUT
    write(*,*) "size(camp_state%state_var)",size(camp_state%state_var), "state_size_cell", state_size_cell

    !open(30, file="../../../../test/chemistry/cb05cl_ae5/files/ebi_output_1_1_48_kss_kae.txt", status="old")
    !open(30, file="../../../../test/chemistry/cb05cl_ae5/files/ebi_output_1_1_9_ksskse_lemisF_ldrydepF_lcldchemF.txt", status="old")
    !open(30, file="../../../../test/chemistry/cb05cl_ae5/files/ebi_output_1_1_9_ksskse_photo0_lemisF_ldrydepF_lcldchemF.txt", status="old")
    !open(30, file="../../../../test/chemistry/cb05cl_ae5/files/ebi_input&
    !        _1_1_48_ksskse_photo0_lemisF_ldrydepF_lcldchemF.txt", status="old")

    !open(30, file="../../../../test/chemistry/cb05cl_ae5/files/ebi_input&
    !        _all_all_48_ksskse_photo0_lemisF_ldrydepF_lcldchemF_reverse.txt", status="old")
    open(30, file="../../../../test/chemistry/cb05cl_ae5/files/ebi_input123&
            _all_all_all_ksskse_photo0_lemisF_ldrydepF_lcldchemF.txt", status="old")

    open(31, file="../../../../test/chemistry/cb05cl_ae5/files/ebi123_temp_press&
            _all_all_all.txt", status="old")

    !Offset
    offset=0 !690 !1
    do i_cell = 1, offset

      !write (*,*) "offset concs hola"
      read(30,*) name_cell !First line indicating layer
      do i_spec = 1, NUM_EBI_SPEC
        read(30,*) real_val!ebi_monarch_init(i_spec)
      end do
      do i_spec = 1, 6 !Monarch extra species (NH3 and others, cut on 72)
        read(30,*) real_val !avoid them atm
      end do
      read(31,*) name_cell
      read(31,*) real_val!temperatures(i_cell+1)
      read(31,*) real_val!pressures(i_cell+1)

    end do

    ! Set the initial concentrations in each module
    !todo suppose not multi_cells
    do i_cell = 0, n_cells-1

      call set_ebi_monarch_species(ebi_monarch_spec_names)
      read(30,*) name_cell !First line indicating layer

      !write(*,*) name_cell

      ! Get ebi in MONARCH order
      do i_spec = 1, NUM_EBI_SPEC
        read(30,*) ebi_monarch_init(i_spec)
      end do

      !do i_spec = 1, 6 !Monarch extra species (NH3 and others, cut on 72)
      !  read(30,*) real_val !avoid them atm
      !end do
      read(30,*) model_conc(i_cell*state_size_cell+i_H2O)

      read(31,*) name_cell
      read(31,*) temperatures(i_cell+1)
      read(31,*) pressures(i_cell+1)

      do i_spec = 1, NUM_EBI_SPEC !72
        ! Get initial concentrations from camp-chem input data
        call assert(787326679, chem_spec_data%get_property_set( &
                ebi_spec_names(i_spec)%string, prop_set))
        if (prop_set%get_real(key, real_val)) then

          !DEBUG Set real_val from ebi_monarch file
          !Loop to search the same name and id from ebi to monarch names (todo there's a better way to do this)
          do j_spec = 1, NUM_EBI_SPEC
            if (trim(ebi_spec_names(i_spec)%string).eq.trim(ebi_monarch_spec_names(j_spec)%string)) then
              real_val = ebi_monarch_init(j_spec)
              !save id mapping ebi-monarch (difference between i_spec and j_spec
              map_ebi_monarch(j_spec)=i_spec
            end if
          end do

          ! Set the EBI solver concetration (ppm)
          YC(i_spec) = real_val

          ! Set the camp-chem concetration (ppm)
          model_conc(chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)+i_cell*state_size_cell) = real_val

        end if

        ! Set KPP species concentrations (#/cc)
        do j_spec = 1, KPP_NSPEC
          if (trim(ebi_spec_names(i_spec)%string).eq.trim(KPP_SPC_NAMES(j_spec))) then
            KPP_C(j_spec) = YC(i_spec) / conv
          end if
        end do
      end do
    end do
    close(30)
    close(31)
#else

    ! Set the initial concentrations in each module
    do i_spec = 1, NUM_EBI_SPEC

      ! Get initial concentrations from camp-chem input data
      call assert(787326679, chem_spec_data%get_property_set( &
              ebi_spec_names(i_spec)%string, prop_set))
      if (prop_set%get_real(key, real_val)) then

        ! Set the EBI solver concetration (ppm)
        YC(i_spec) = real_val

        ! Set the camp-chem concetration (ppm)
        camp_state%state_var( &
                chem_spec_data%gas_state_id( &
                        ebi_spec_names(i_spec)%string)) = real_val

      end if

      ! Set KPP species concentrations (#/cc)
      do j_spec = 1, KPP_NSPEC
        if (trim(ebi_spec_names(i_spec)%string).eq.trim(KPP_SPC_NAMES(j_spec))) then
          KPP_C(j_spec) = YC(i_spec) / conv
        end if
      end do
    end do
#endif

    write (*,*) "Setup model_concs not ebi species"

    ! Set EBI solver constant species concentrations in CAMP-chem
    spec_name = "M"
    call assert(273497194, chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(740666066, associated(prop_set))
    call assert(907464197, prop_set%get_real(key, real_val))
    camp_state%state_var(i_M) = real_val
    KPP_M = real_val / conv
    spec_name = "O2"
    call assert(557877977, chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(729136508, associated(prop_set))
    call assert(223930103, prop_set%get_real(key, real_val))
    camp_state%state_var(i_O2) = real_val
    KPP_O2 = real_val / conv
    KPP_C(KPP_IND_O2) = real_val / conv
    ! KPP has variable O2 concentration
    do j_spec = 1, KPP_NSPEC
      if (trim(KPP_SPC_NAMES(j_spec)).eq.'O2') then
        KPP_C(j_spec) = real_val / conv
      end if
    end do
    spec_name = "N2"
    call assert(329882514, chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(553715297, associated(prop_set))
    call assert(666033642, prop_set%get_real(key, real_val))
    camp_state%state_var(i_N2) = real_val
    KPP_N2 = real_val / conv
    spec_name = "H2O"
    call assert(101887051, chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(160827237, associated(prop_set))
    call assert(273145582, prop_set%get_real(key, real_val))
    camp_state%state_var(i_H2O) = real_val
    KPP_H2O = real_val / conv
    spec_name = "CH4"
    call assert(208941089, chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(667939176, associated(prop_set))
    call assert(780257521, prop_set%get_real(key, real_val))
    camp_state%state_var(i_CH4) = real_val
    KPP_CH4 = real_val / conv
    ! KPP has variable CH4 concentration
    do j_spec = 1, KPP_NSPEC
      if (trim(KPP_SPC_NAMES(j_spec)).eq.'CH4') then
        KPP_C(j_spec) = real_val / conv
      end if
    end do
    spec_name = "H2"
    call assert(663478776, chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(892575866, associated(prop_set))
    call assert(722418962, prop_set%get_real(key, real_val))
    camp_state%state_var(i_H2) = real_val
    KPP_H2 = real_val / conv

    ! Set the water concentration for EBI solver (ppmV)
    water_conc = camp_state%state_var(i_H2O)

    !if(pmc_multicells.eq.0) then
    spec_name = "DUMMY"
    i_DUMMY = chem_spec_data%gas_state_id(spec_name)
    call assert(663478276, chem_spec_data%get_property_set(spec_name, prop_set))
    call assert(892572866, associated(prop_set))
    call assert(722411962, prop_set%get_real(key, real_val))

#ifdef PMC_MONARCH_INPUT
    do i_cell = 0, n_cells-1
      model_conc(i_cell*state_size_cell+i_DUMMY) = real_val
      model_conc(i_cell*state_size_cell+i_M)=camp_state%state_var(i_M)
      model_conc(i_cell*state_size_cell+i_O2)=camp_state%state_var(i_O2)
      model_conc(i_cell*state_size_cell+i_N2)=camp_state%state_var(i_N2)
      !model_conc(i_cell*state_size_cell+i_H2O)=camp_state%state_var(i_H2O)
      model_conc(i_cell*state_size_cell+i_CH4)=camp_state%state_var(i_CH4) !interesting, when setting wrong conc like i_h20 it takes soo long
      model_conc(i_cell*state_size_cell+i_H2)=camp_state%state_var(i_H2)
    end do
#endif

    ! Set up the output files
    open(EBI_FILE_UNIT, file="out/cb05cl_ae5_ebi_results.txt", status="replace", &
            action="write")
    open(KPP_FILE_UNIT, file="out/cb05cl_ae5_kpp_results.txt", status="replace", &
            action="write")
    open(CAMP_FILE_UNIT, file="out/cb05cl_ae5_camp_results.txt", status="replace", &
            action="write")
#ifdef PMC_MONARCH_INPUT
    open(CAMP_EBI_FILE_UNIT, file="out/cb05cl_ae5_camp_ebi_diff.txt", status="replace", &
            action="write")
    open(CAMP_KPP_FILE_UNIT, file="out/cb05cl_ae5_camp_kpp_diff.txt", status="replace", &
            action="write")
    open(EBI_KPP_FILE_UNIT, file="out/cb05cl_ae5_ebi_kpp_diff.txt", status="replace", &
            action="write")
#endif
    n_gas_spec = chem_spec_data%size(spec_phase=CHEM_SPEC_GAS_PHASE)
    allocate(camp_spec_names(n_gas_spec))
    do i_spec = 1, n_gas_spec
      camp_spec_names(i_spec)%string = chem_spec_data%gas_state_name(i_spec)
    end do
    write(EBI_FILE_UNIT,*) "time ", (ebi_spec_names(i_spec)%string//" ", i_spec=1, &
    NUM_EBI_SPEC)
    write(KPP_FILE_UNIT,*) "time ", (trim(KPP_SPC_NAMES(i_spec))//" ", i_spec=1, &
    KPP_NSPEC)
    write(CAMP_FILE_UNIT,*) "time ", (camp_spec_names(i_spec)%string//" ", i_spec=1, &
            size(camp_spec_names))
#ifdef PMC_MONARCH_INPUT
    write(CAMP_EBI_FILE_UNIT,*) "CAMP order: ", (camp_spec_names(i_spec)%string//" ", i_spec=1, &
            NUM_EBI_SPEC)
#endif

    ! Set up the reaction map between camp-chem, kpp and ebi solvers
    key = "rxn id"
    allocate(ebi_rxn_map(mechanism%size()))
    ebi_rxn_map(:) = 0
    allocate(kpp_rxn_map(mechanism%size()))
    kpp_rxn_map(:) = 0
    call get_kpp_rxn_labels(kpp_rxn_labels)
    do i_rxn = 1, mechanism%size()
      rxn => mechanism%get_rxn(i_rxn)
      call assert_msg(917216189, associated(rxn), "Missing rxn "//to_string(i_rxn))
      call assert(656034097, rxn%property_set%get_string(key, string_val))
      do i_ebi_rxn = 1, NRXNS
        if (trim(RXLABEL(i_ebi_rxn)).eq.trim(string_val)) then
          ebi_rxn_map(i_rxn) = i_ebi_rxn
          exit
        end if
      end do
      if (trim(string_val).ne.'jo2') then ! jo2 rxn O2 + hv is not in EBI solver
        call assert_msg(921715481, ebi_rxn_map(i_rxn).ne.0, "EBI missing rxn "//string_val)
      end if
      do i_kpp_rxn = 1, KPP_NREACT
        if (trim(kpp_rxn_labels(i_kpp_rxn)%string).eq.trim(string_val)) then
          kpp_rxn_map(i_rxn) = i_kpp_rxn
          exit
        end if
      end do
      call assert_msg(360769001, kpp_rxn_map(i_rxn).ne.0, "KPP missing rxn "//string_val)
    end do

    ! Set up the species map between camp-chem, kpp and ebi solvers
    allocate(ebi_spec_map(chem_spec_data%size()))
    allocate(kpp_spec_map(chem_spec_data%size()))
    ebi_spec_map(:) = 0
    kpp_spec_map(:) = 0
    do i_spec = 1, NUM_EBI_SPEC
      j_spec = chem_spec_data%gas_state_id(ebi_spec_names(i_spec)%string)
      call assert_msg(194404050, j_spec.gt.0, "Missing EBI species: "//trim(ebi_spec_names(i_spec)%string))
      ebi_spec_map(j_spec) = i_spec
    end do
    do i_spec = 1, KPP_NSPEC
      spec_name = trim(KPP_SPC_NAMES(i_spec))
      j_spec = chem_spec_data%gas_state_id(spec_name)
      call assert_msg(194404050, j_spec.gt.0, "Missing KPP species: "//trim(KPP_SPC_NAMES(i_spec)))
      kpp_spec_map(j_spec) = i_spec
    end do

    ! Reset the computation timers
    comp_ebi = 0.0
    comp_kpp = 0.0
    comp_camp = 0.0

    ! Compare the rates for the initial conditions
    call EXT_HRCALCKS( NUM_EBI_PHOTO_RXN,       & ! Number of EBI solver photolysis reactions
            is_sunny,                & ! Flag for sunlight
            photo_rates,             & ! Photolysis rates
            temperature,             & ! Temperature (K)
            pressure,                & ! Air pressure (atm)
            water_conc,              & ! Water vapor concentration (ppmV)
            RKI)                       ! Rate constants
    KPP_DT = EBI_TMSTEP * 60.0d0 ! KPP Time step in seconds
    KPP_TIME = 0.0
    KPP_TEMP = temperature
    KPP_PRESS = pressure * 1013.25 ! KPP pressure in hPa
    CALL KPP_Update_RCONST()

#ifdef PMC_MONARCH_INPUT
    state_size_cell = size(chem_spec_data%get_spec_names())
#else
    state_size_cell = size(chem_spec_data%get_spec_names()) !size(camp_state%state_var) / n_cells
#endif
    spec_names = chem_spec_data%get_spec_names()

    !todo netcdf with monarch input instead txt file
    !Netcdf n cells exp values
#ifdef NETCDF_INPUT
    call set_input_from_netcdf(ncfile, state_size_cell, spec_names&
            ,model_conc, temperatures, pressures, n_cells)
#endif

    !write(*,*) "First setup camp_state_var"

    ! Set same conc per n_cells
    !do i_block = 0, n_blocks-1
      do i_cell = 0, n_cells_block-1
        do j = 1, state_size_cell
          if(pmc_multicells.eq.1) then
            !todo is this necessary when everything is stored in model_conc and updated later?
#ifdef PMC_MONARCH_INPUT
            camp_state%state_var(i_cell*state_size_cell+j) = camp_state%state_var(i_cell*state_size_cell+j) + offset_conc*i_cell !0.1*j
            !camp_state%state_var(i_cell*state_size_cell+j) = model_conc(i_block*(state_size_cell*n_cells_block)+i_cell*state_size_cell+j) + offset_conc*i_cell
#else
            camp_state%state_var(i_cell*state_size_cell+j) = camp_state%state_var(j) + offset_conc*i_cell !0.1*j !todo this should be i_cell...repeat tests
#endif
          else
#ifdef PMC_MONARCH_INPUT
            camp_state%state_var(j) = model_conc(j) + offset_conc*i_cell !+ 0.1*i_cell
#else
            model_conc(i_cell*state_size_cell+j) = camp_state%state_var(j) + offset_conc*i_cell !+ 0.1*i_cell !todo improve this part (make it more clear)
#endif
          endif
        end do
      end do
    !end do

    if(pmc_multicells.eq.1) then
      !Set default temperature & pressure to all cells
      do i_cell = 1, n_cells_block
        call camp_state%env_states(i_cell)%set_temperature_K( real( temperatures(i_cell)+offset_temp*(i_cell-1), kind=dp ) )
        call camp_state%env_states(i_cell)%set_pressure_Pa( pressures(i_cell) )
      end do
    else
      call camp_state%env_states(1)%set_temperature_K( real( temperatures(1), kind=dp ) )
      call camp_state%env_states(1)%set_pressure_Pa( pressures(1))
    endif
    ! Save the initial states for repeat calls
    allocate(ebi_init(size(YC)))
    allocate(kpp_init(size(KPP_C)))

!#ifdef PMC_USE_MPI
    !allocate(camp_init(size(camp_state%state_var/mpi_threads)))
    allocate(camp_init(size(camp_state%state_var)))
!#else
!    allocate(camp_init(size(camp_state%state_var)))
!#endif

    ebi_init(:) = YC(:)
    kpp_init(:) = KPP_C(:)
    camp_init(:) = camp_state%state_var(:)

#ifdef PRINT_STATE_INPUT
    do i_spec = 1, NUM_EBI_SPEC
      write(EBI_KPP_FILE_UNIT,*) "(id), spec_name, camp input"
      associate (camp_var=>camp_state%state_var( &
              chem_spec_data%gas_state_id( &
                      ebi_monarch_spec_names(i_spec)%string)))
        write(*,*) i_spec, ebi_monarch_spec_names(i_spec)%string, &
                camp_var
      end associate
    end do
#endif

    !write(*,*) model_conc(:)
    !write(*,*) temperatures(:)

    ! Repeatedly solve the mechanism
    do i_repeat = 1, n_repeats

      write(*,*) "Repeat", i_repeat
      !print*, "running"

      YC(:) = ebi_init(:)
      KPP_C(:) = kpp_init(:)
      camp_state%state_var(:) = camp_init(:)

      ! Solve the mechanism
      do i_time = 1,NUM_TIME_STEPS

        ! Set minimum concentrations in all solvers
        YC(:) = MAX(YC(:), SMALL_NUM)
        KPP_C(:) = MAX(KPP_C(:), SMALL_NUM/conv)
        camp_state%state_var(:) = max(camp_state%state_var(:), SMALL_NUM)

        ! EBI solver
        call cpu_time(comp_start)
        call EXT_HRCALCKS( NUM_EBI_PHOTO_RXN,       & ! Number of EBI solver photolysis reactions
                is_sunny,                & ! Flag for sunlight
                photo_rates,             & ! Photolysis rates
                temperature,             & ! Temperature (K)
                pressure,                & ! Air pressure (atm)
                water_conc,              & ! Water vapor concentration (ppmV)
                RKI)                       ! Rate constants
        call EXT_HRSOLVER( 2018012, 070000, 1, 1, 1)  ! These dummy variables are just for output
        call cpu_time(comp_end)
        comp_ebi = comp_ebi + (comp_end-comp_start)

        ! Set KPP and camp-chem concentrations to those of EBI at first time step to match steady-state
        ! EBI species
        if (i_time.eq.1) then
          ! Set KPP species in #/cc
          do i_spec = 1, NUM_EBI_SPEC
            do j_spec = 1, KPP_NSPEC
              if (trim(ebi_spec_names(i_spec)%string).eq.trim(KPP_SPC_NAMES(j_spec))) then
                KPP_C(j_spec) = YC(i_spec) / conv
              end if
            end do
        if(pmc_multicells.eq.1) then
            do i = 0, n_cells_block-1
#ifdef PMC_MONARCH_INPUT
#else
              !camp_state%state_var( &
              !  chem_spec_data%gas_state_id( &
              !  ebi_spec_names(i_spec)%string) &
              !  + i*state_size_cell)= YC(i_spec)
#endif
            end do
        else
            !camp_state%state_var( &
            !        chem_spec_data%gas_state_id( &
            !                ebi_spec_names(i_spec)%string)) = YC(i_spec)
        endif
          end do
        end if

        ! KPP module
        call cpu_time(comp_start)
        KPP_TIME = (i_time-1)*KPP_DT
        KPP_TEMP = temperature
        KPP_PRESS = pressure * 1013.25 ! KPP pressure in hPa
        CALL KPP_Update_RCONST()
        CALL KPP_INTEGRATE( TIN = KPP_TIME, TOUT = (KPP_TIME+KPP_DT), &
                RSTATUS_U = KPP_RSTATE, ICNTRL_U = KPP_ICNTRL )
        call cpu_time(comp_end)
        comp_kpp = comp_kpp + (comp_end-comp_start)

        if(pmc_multicells.eq.1) then
          do i_block = 0, n_blocks-1
            do i_cell = 0, n_cells_block-1
              call camp_state%env_states(i_cell+1)%set_temperature_K( real(temperatures( i_block*n_cells_block+i_cell+1), kind=dp ))
              call camp_state%env_states(i_cell+1)%set_pressure_Pa( pressures(i_block*n_cells_block+i_cell+1))

              do j = 1, state_size_cell
                camp_state%state_var(state_size_cell*i_cell+j) = model_conc(i_block*(state_size_cell*n_cells_block)+state_size_cell*i_cell+j)
              end do
            end do

            call cpu_time(comp_start)
            call camp_core%solve(camp_state, real(EBI_TMSTEP*60.0, kind=dp), &
                    solver_stats = solver_stats)
            !write (*,*) "a", mod(n_cells_block,i_cell)
            call cpu_time(comp_end)
            comp_camp = comp_camp + (comp_end-comp_start)
          end do
        else
          do i_cell = 0, n_cells-1

            call camp_state%env_states(1)%set_temperature_K( real( temperatures(i_cell+1)+offset_temp, kind=dp ) )
            call camp_state%env_states(1)%set_pressure_Pa( pressures(i_cell+1))

            do j = 1, state_size_cell
              camp_state%state_var(j) = model_conc(i_cell*state_size_cell+j)
            end do

            ! CAMP-chem
            call cpu_time(comp_start)
            !call date_and_time(time=t)
            !write(*, "(' The current time is ', A8 )")  t
            call camp_core%solve(camp_state, real(EBI_TMSTEP*60.0, kind=dp), &
                    solver_stats = solver_stats)
            !call date_and_time(time=t)
            !write(*, "(' The current time is ', A8 )")  t
            call cpu_time(comp_end)
            comp_camp = comp_camp + (comp_end-comp_start)
          end do
        endif
      end do
    end do

    ! Output the computational time
    write(*,*) "CAMP-chem calculation time: ", comp_camp," s"

    ! Output final timestep
    write(EBI_FILE_UNIT,*) i_time*EBI_TMSTEP, YC(:)
    write(KPP_FILE_UNIT,*) i_time*EBI_TMSTEP, KPP_C(:)*conv
    write(CAMP_FILE_UNIT,*) i_time*EBI_TMSTEP, camp_state%state_var(:)

    ! Output the computational time
    write(*,*) "EBI calculation time: ", comp_ebi," s"
    write(*,*) "KPP calculation time: ", comp_kpp," s"
    write(*,*) "CAMP-chem calculation time: ", comp_camp," s"

#ifdef PMC_MONARCH_INPUT
    write(CAMP_EBI_FILE_UNIT,*) "Repeat", i_repeat, "timestep ", NUM_TIME_STEPS
    write(CAMP_EBI_FILE_UNIT,*) "spec_name, concentrations rel. error [(camp_state-ebi)/(camp_state+ebi)], camp_state, ebi"
    write(CAMP_KPP_FILE_UNIT,*) "Repeat", i_repeat, "timestep ", NUM_TIME_STEPS
    write(CAMP_KPP_FILE_UNIT,*) "spec_name, concentrations rel. error [(camp_state-kpp)/(camp_state+kpp)], camp_state, kpp"
    write(EBI_KPP_FILE_UNIT,*) "Repeat", i_repeat, "timestep ", NUM_TIME_STEPS
    write(EBI_KPP_FILE_UNIT,*) "spec_name, concentrations rel. error [(èbi-kpp)/(ebi+kpp)], ebi, kpp"
#endif

    ! Compare the results
    ! EBI <-> CAMP-chem
    do i_spec = 1, NUM_EBI_SPEC
#ifdef PMC_COMPARE
      call assert_msg(749090387, almost_equal(real(YC(i_spec), kind=dp), &
          camp_state%state_var( &
                  chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)), 1.0d-4, 1.0d-3),& !.or. & !, 1.0d-3
          !(YC(i_spec).lt.ebi_init(i_spec)*1.0d-2 .and. &
           !camp_state%state_var(chem_spec_data%gas_state_id( &
          !        ebi_spec_names(i_spec)%string)) .lt. &
          ! camp_init(chem_spec_data%gas_state_id( &
          !        ebi_spec_names(i_spec)%string))*1.0d-2), &
          "Species "//ebi_spec_names(i_spec)%string//" has different result. "// &
          "EBI solver: "//trim(to_string(real(YC(i_spec), kind=dp)))// &
          "; CAMP-chem: "// &
          trim(to_string( camp_state%state_var( &
                  chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)))) // "; ebi init: "// &
          trim(to_string(real(ebi_init(i_spec), kind=dp)))//"; camp init: "// &
          trim(to_string(camp_init(chem_spec_data%gas_state_id( &
                  ebi_spec_names(i_spec)%string)))))
#endif

#ifdef PMC_MONARCH_INPUT
      associate (camp_var=>camp_state%state_var( &
              chem_spec_data%gas_state_id( &
                      ebi_monarch_spec_names(i_spec)%string)))
        write(CAMP_EBI_FILE_UNIT,*) ebi_monarch_spec_names(i_spec)%string, &
                (abs(camp_var) - abs(YC(map_ebi_monarch(i_spec)))) / &
                        (abs(camp_var) + abs(YC(map_ebi_monarch(i_spec)))+ &
                                1.0d-30), & !avoid division by zero, &
                camp_var, &
                YC(map_ebi_monarch(i_spec))
      end associate
#endif
    end do

#ifdef PMC_COMPARE_KPP
    ! KPP <-> CAMP-chem
    do i_spec = 1, KPP_NSPEC
      str_temp%string = trim(KPP_SPC_NAMES(i_spec))
      call assert_msg(749090436, almost_equal(real(KPP_C(i_spec)*conv, kind=dp), &
          camp_state%state_var( &
                  chem_spec_data%gas_state_id( &
                  str_temp%string)), 1.0d-4, 1.0d-3),& !.or. &
          !(KPP_C(i_spec) .lt. KPP_init(i_spec)*1.0d-2 .and. &
          ! camp_state%state_var(chem_spec_data%gas_state_id( &
          !        str_temp%string)) .lt. &
          ! camp_init(chem_spec_data%gas_state_id( &
          !        str_temp%string))*1.0d-2), &
          "Species "//str_temp%string//" has different result. "// &
          "KPP solver: "//trim(to_string(real(KPP_C(i_spec)*conv, kind=dp)))// &
          "; CAMP-chem: "// &
          trim(to_string( camp_state%state_var( &
                  chem_spec_data%gas_state_id( &
                  str_temp%string))))//"; KPP init: "// &
          trim(to_string(real(kpp_init(i_spec)*conv, kind=dp)))// &
          "; camp init: "//trim(to_string(camp_init( &
                  chem_spec_data%gas_state_id( &
                  str_temp%string)))))

      associate (camp_var=>camp_state%state_var( &
              chem_spec_data%gas_state_id( &
                      str_temp%string)))
        write(CAMP_KPP_FILE_UNIT,*) str_temp%string, &
                (abs(camp_var) - abs(real(KPP_C(i_spec)*conv, kind=dp))) / &
                        (abs(camp_var) + abs(real(KPP_C(i_spec)*conv, kind=dp))+ &
                                1.0d-30), & !avoid division by zero, &
                camp_var, &
                real(KPP_C(i_spec)*conv, kind=dp)
      end associate
    end do
#endif

    ! Close the output files
    close(EBI_FILE_UNIT)
    close(KPP_FILE_UNIT)
    close(CAMP_FILE_UNIT)

#ifdef PMC_MONARCH_INPUT
    close(CAMP_EBI_FILE_UNIT)
    close(CAMP_KPP_FILE_UNIT)
    close(EBI_KPP_FILE_UNIT)
#endif

    ! Create the gnuplot script
    open(unit=CAMP_FILE_UNIT, file="out/plot_cb05cl_ae.conf", status="replace", action="write")
    write(CAMP_FILE_UNIT,*) "# plot_cb05cl_ae5.conf"
    write(CAMP_FILE_UNIT,*) "# Run as: gnuplot plot_cb05cl_ae5.conf"
    write(CAMP_FILE_UNIT,*) "set terminal png truecolor"
    write(CAMP_FILE_UNIT,*) "set autoscale"
    spec_names = chem_spec_data%get_spec_names()
    do i_spec = 1, size(spec_names)
      write(CAMP_FILE_UNIT,*) "set output 'cb05cl_ae5_"//trim(spec_names(i_spec)%string)//".png'"
      write(CAMP_FILE_UNIT,*) "plot\"
      if (ebi_spec_map(i_spec).gt.0) then
        write(CAMP_FILE_UNIT,*) " 'cb05cl_ae5_ebi_results.txt'\"
        write(CAMP_FILE_UNIT,*) " using 1:"//trim(to_string(ebi_spec_map(i_spec)+1))//" title '"// &
                trim(spec_names(i_spec)%string)//" (ebi)',\"
      end if
      if (kpp_spec_map(i_spec).gt.0) then
        write(CAMP_FILE_UNIT,*) " 'cb05cl_ae5_kpp_results.txt'\"
        write(CAMP_FILE_UNIT,*) " using 1:"//trim(to_string(kpp_spec_map(i_spec)+1))//" title '"// &
                trim(spec_names(i_spec)%string)//" (kpp)',\"
      end if
      write(CAMP_FILE_UNIT,*) " 'cb05cl_ae5_camp_results.txt'\"
      write(CAMP_FILE_UNIT,*) " using 1:"//trim(to_string(i_spec+1))//" title '"// &
              trim(spec_names(i_spec)%string)//" (camp)'"
    end do
    close(CAMP_FILE_UNIT)

    deallocate(photo_rates)
    deallocate(rate_update)
    deallocate(camp_spec_names)
    deallocate(ebi_rxn_map)
    deallocate(kpp_rxn_map)
    deallocate(kpp_rxn_labels)
    deallocate(ebi_spec_map)
    deallocate(kpp_spec_map)
    deallocate(ebi_init)
    deallocate(kpp_init)
    deallocate(camp_init)
    !deallocate(camp_state_comp)
    deallocate(camp_state)
    deallocate(camp_core)

    passed = .true.

    call pmc_mpi_finalize( )

  end function run_standard_cb05cl_ae5_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_input_from_netcdf(camp_state, state_size_cell, spec_names&
          ,model_conc, temperatures, pressures, n_cells)

    !todo fix arguments
    type(camp_state_t), intent(inout) :: camp_state
    type(string_t), intent(in) :: spec_names(:)
    integer(kind=i_kind), intent(in) :: state_size_cell, n_cells
    real(kind=dp), intent(inout) :: model_conc(:)
    real(kind=dp), intent(inout) :: temperatures(:)
    real(kind=dp), intent(inout) :: pressures(:)

    real(kind=dp), allocatable :: aux_array(:,:)
    real(kind=dp), allocatable :: water_concs(:)
    integer(kind=i_kind) :: stat
    real(kind=dp) :: comp_start, comp_end, comp_camp

    integer(kind=i_kind) :: i, j, k, s, i_cell, i_spec
    integer(kind=i_kind) :: n_repeats
    !netcdf
    integer(kind=i_kind) :: input_id, varid
    character(len=:), allocatable :: ncfile

    integer :: max_k = 46
    integer :: max_j = 15
    integer :: max_i = 15
    integer :: i_start = 1 !136
    integer :: j_start= 1 !134
    integer :: k_start = 3 !1 !cell 1 and 2 are teorically non-chemistry computant
    integer :: t_start = 1 !14
    integer :: i_n = 5 !20
    integer :: j_n = 5 !20
    integer :: k_n = 5 !48
    integer :: t_n = 1 !1

#ifdef PMC_USE_MPI
    if( pmc_mpi_rank( ) .eq. 0 ) then
#endif

#ifdef ARGS_NETCDF
    ! Get from the .sh script the command arguments
    call read_args(i_start, MAP_I_START)
    call read_args(j_start, MAP_J_START)
    call read_args(k_start, MAP_K_START)
    call read_args(t_start, MAP_T_START)
    call read_args(i_n, MAP_I_N)
    call read_args(j_n, MAP_J_N)
    call read_args(k_n, MAP_K_N)
    call read_args(t_n, MAP_T_N)
#else
    !todo: we don't know how many cells has each rank (15 or 16). We will suppose 15 atm
    !todo check everything is multiple of max i, j and k and compute the non-multiple cells
    !set to 46 if >=46
    if(n_cells.ge.max_i) then
      k_n=max_k
      if(n_cells.ge.max_k*max_j) then
        j_n=max_j
        if(n_cells.ge.max_k*max_j*max_i) then
          i_n=n_cells/(max_k*max_j)
          write (*,*) "WARNING More cells than maxs ", max_k, max_j, max_i, ", not computing the rest"
        else
          i_n=n_cells/(max_k*max_j)
        end if
      else
        j_n=n_cells/max_k
      end if
    else
      k_n = n_cells
    end if

#endif

    !ncfile = '/esarchive/exp/monarch/a2bk/original_files/000/2016083012/MONARCH_d01_2016083012.nc'
    ncfile = '/gpfs/scratch/bsc32/bsc32815/a2s8/nmmb-monarch/OUTPUT/regional/000/20160801/CURRENT_RUN/&
            nmmb_hst24hAllBlockMulticell_01_nc4_0000h_00m_00.00s_copy.nc'

    stat =  nf90_open &
            (ncfile, &
            NF90_NOWRITE,input_id)

    allocate(aux_array(n_cells,state_size_cell))
    allocate(water_concs(n_cells))

    !CONCS
    call cpu_time(comp_start)
    do i_spec = 1 ,state_size_cell

      !Save concs in temporal array
      stat =  nf90_inq_varid(input_id,spec_names(i_spec)%string,varid)
      stat =  nf90_get_var(input_id,varid,aux_array(:,i_spec), &
              start=(/i_start,j_start,k_start,t_start/),count=(/i_n,j_n,k_n,t_n/))
      !start=(/i_start,j_start,k_start,t_start/),count=(/1,1,1,1/))


      !print *, "spec_names(i_spec)%string ",spec_names(i_spec)%string
    end do
    !print *, "spec_names(i_spec)%string, aux_array(1,1) ",spec_names(i_spec)%string, aux_array(1,1)
    !todo: change prints for tests that check the correct reading of netcdf variables

    !Set in state array
    do i_cell=0, n_cells-1
      do i_spec = 1, state_size_cell
        !camp_state%state_var(i_cell*state_size_cell+i_spec) = aux_array(i_cell,i_spec)
        model_conc(i_cell*state_size_cell+i_spec) = aux_array(i_cell,i_spec)
        !print*, i_cell, spec_names(i_spec)%string, camp_state%state_var(i_cell*state_size_cell+i_spec)
      end do
    end do

    !temps
    stat =  nf90_inq_varid(input_id,'T',varid)
    stat =  nf90_get_var(input_id,varid,temperatures, &
            start=(/i_start,j_start,k_start,t_start/),count=(/i_n,j_n,k_n,t_n/))
    !start=(/i_start,j_start,k_start,t_start/),count=(/1,1,1,1/))
    !return temperatures and later will set them to correct env_state

    !todo: pressure
    !do i_cell = 1, n_cells
    !  call camp_state%env_states(i_cell)%set_temperature_K( real( temperatures(i_cell), kind=dp ) )
    !  call camp_state%env_states(i_cell)%set_pressure_Pa( pressure * const%air_std_press )
    !end do

    call cpu_time(comp_end)
    comp_camp = comp_camp + (comp_end-comp_start)
    print*, "netcdf reading time", comp_camp
    comp_camp = 0.0

    stat = nf90_close(input_id)
    deallocate (aux_array)
    deallocate(water_concs)

#ifdef PMC_USE_MPI
    else
      write (*,*) "Not rank 0!"
    end if
#endif

  end subroutine set_input_from_netcdf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the EBI-solver species names
  subroutine set_ebi_species(spec_names)

    !> EBI solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: spec_names

    spec_names(1)%string = "NO2"
    spec_names(2)%string = "NO"
    spec_names(3)%string = "O"
    spec_names(4)%string = "O3"
    spec_names(5)%string = "NO3"
    spec_names(6)%string = "O1D"
    spec_names(7)%string = "OH"
    spec_names(8)%string = "HO2"
    spec_names(9)%string = "N2O5"
    spec_names(10)%string = "HNO3"
    spec_names(11)%string = "HONO"
    spec_names(12)%string = "PNA"
    spec_names(13)%string = "H2O2"
    spec_names(14)%string = "XO2"
    spec_names(15)%string = "XO2N"
    spec_names(16)%string = "NTR"
    spec_names(17)%string = "ROOH"
    spec_names(18)%string = "FORM"
    spec_names(19)%string = "ALD2"
    spec_names(20)%string = "ALDX"
    spec_names(21)%string = "PAR"
    spec_names(22)%string = "CO"
    spec_names(23)%string = "MEO2"
    spec_names(24)%string = "MEPX"
    spec_names(25)%string = "MEOH"
    spec_names(26)%string = "HCO3"
    spec_names(27)%string = "FACD"
    spec_names(28)%string = "C2O3"
    spec_names(29)%string = "PAN"
    spec_names(30)%string = "PACD"
    spec_names(31)%string = "AACD"
    spec_names(32)%string = "CXO3"
    spec_names(33)%string = "PANX"
    spec_names(34)%string = "ROR"
    spec_names(35)%string = "OLE"
    spec_names(36)%string = "ETH"
    spec_names(37)%string = "IOLE"
    spec_names(38)%string = "TOL"
    spec_names(39)%string = "CRES"
    spec_names(40)%string = "TO2"
    spec_names(41)%string = "TOLRO2"
    spec_names(42)%string = "OPEN"
    spec_names(43)%string = "CRO"
    spec_names(44)%string = "MGLY"
    spec_names(45)%string = "XYL"
    spec_names(46)%string = "XYLRO2"
    spec_names(47)%string = "ISOP"
    spec_names(48)%string = "ISPD"
    spec_names(49)%string = "ISOPRXN"
    spec_names(50)%string = "TERP"
    spec_names(51)%string = "TRPRXN"
    spec_names(52)%string = "SO2"
    spec_names(53)%string = "SULF"
    spec_names(54)%string = "SULRXN"
    spec_names(55)%string = "ETOH"
    spec_names(56)%string = "ETHA"
    spec_names(57)%string = "CL2"
    spec_names(58)%string = "CL"
    spec_names(59)%string = "HOCL"
    spec_names(60)%string = "CLO"
    spec_names(61)%string = "FMCL"
    spec_names(62)%string = "HCL"
    spec_names(63)%string = "TOLNRXN"
    spec_names(64)%string = "TOLHRXN"
    spec_names(65)%string = "XYLNRXN"
    spec_names(66)%string = "XYLHRXN"
    spec_names(67)%string = "BENZENE"
    spec_names(68)%string = "BENZRO2"
    spec_names(69)%string = "BNZNRXN"
    spec_names(70)%string = "BNZHRXN"
    spec_names(71)%string = "SESQ"
    spec_names(72)%string = "SESQRXN"

  end subroutine set_ebi_species

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the EBI-solver species names in MONARCH order
  subroutine set_ebi_monarch_species(spec_names)

    !> EBI solver species names
    type(string_t), dimension(NUM_EBI_SPEC) :: spec_names

    !Monarch order
    spec_names(1)%string = "NO2"
    spec_names(2)%string = "NO"
    spec_names(3)%string = "O3"
    spec_names(4)%string = "NO3"
    spec_names(5)%string = "N2O5"
    spec_names(6)%string = "HNO3"
    spec_names(7)%string = "HONO"
    spec_names(8)%string = "PNA"
    spec_names(9)%string = "H2O2"
    spec_names(10)%string = "NTR"
    spec_names(11)%string = "ROOH"
    spec_names(12)%string = "FORM"
    spec_names(13)%string = "ALD2"
    spec_names(14)%string = "ALDX"
    spec_names(15)%string = "PAR"
    spec_names(16)%string = "CO"
    spec_names(17)%string = "MEPX"
    spec_names(18)%string = "MEOH"
    spec_names(19)%string = "FACD"
    spec_names(20)%string = "PAN"
    spec_names(21)%string = "PACD"
    spec_names(22)%string = "AACD"
    spec_names(23)%string = "PANX"
    spec_names(24)%string = "OLE"
    spec_names(25)%string = "ETH"
    spec_names(26)%string = "IOLE"
    spec_names(27)%string = "TOL"
    spec_names(28)%string = "CRES"
    spec_names(29)%string = "OPEN"
    spec_names(30)%string = "MGLY"
    spec_names(31)%string = "XYL"
    spec_names(32)%string = "ISOP"
    spec_names(33)%string = "ISPD"
    spec_names(34)%string = "TERP"
    spec_names(35)%string = "SO2"
    spec_names(36)%string = "SULF"
    spec_names(37)%string = "ETOH"
    spec_names(38)%string = "ETHA"
    spec_names(39)%string = "CL2"
    spec_names(40)%string = "HOCL"
    spec_names(41)%string = "FMCL"
    spec_names(42)%string = "HCL"
    spec_names(43)%string = "BENZENE"
    spec_names(44)%string = "SESQ"
    spec_names(45)%string = "O"
    spec_names(46)%string = "O1D"
    spec_names(47)%string = "OH"
    spec_names(48)%string = "HO2"
    spec_names(49)%string = "XO2"
    spec_names(50)%string = "XO2N"
    spec_names(51)%string = "MEO2"
    spec_names(52)%string = "HCO3"
    spec_names(53)%string = "C2O3"
    spec_names(54)%string = "CXO3"
    spec_names(55)%string = "ROR"
    spec_names(56)%string = "TO2"
    spec_names(57)%string = "TOLRO2"
    spec_names(58)%string = "CRO"
    spec_names(59)%string = "XYLRO2"
    spec_names(60)%string = "ISOPRXN"
    spec_names(61)%string = "TRPRXN"
    spec_names(62)%string = "SULRXN"
    spec_names(63)%string = "CL"
    spec_names(64)%string = "CLO"
    spec_names(65)%string = "TOLNRXN"
    spec_names(66)%string = "TOLHRXN"
    spec_names(67)%string = "XYLNRXN"
    spec_names(68)%string = "XYLHRXN"
    spec_names(69)%string = "BENZRO2"
    spec_names(70)%string = "BNZNRXN"
    spec_names(71)%string = "BNZHRXN"
    spec_names(72)%string = "SESQRXN"

  end subroutine set_ebi_monarch_species

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Load a string array with KPP reaction labels
  subroutine get_kpp_rxn_labels(kpp_rxn_labels)

    use cb05cl_ae5_Model,                       only : NREACT

    type(string_t), allocatable :: kpp_rxn_labels(:)
    integer(kind=i_kind) :: i_rxn

    i_rxn = 1

    allocate(kpp_rxn_labels(NREACT))

    kpp_rxn_labels(i_rxn)%string = 'R1'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R2'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R3'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R4'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R5'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R6'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R7'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R8'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R9'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R10'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R11'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R12'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R13'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R14'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R15'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R16'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R17'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R18'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R19'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R20'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R21'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R22'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R23'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R24'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R25'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R26'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R27'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R28'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R29'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R30'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R31'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R32'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R33'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R34'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R35'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R36'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R37'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R38'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R39'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R40'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R41'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R42'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R43'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R44'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R45'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R46'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R47'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R48'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R49'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R50'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R51'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R52'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R53'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R54'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R55'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R56'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R57'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R58'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R59'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R60'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R61'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R62'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R63'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R64'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R65'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R66'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R67'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R68'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R69'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R70'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R71'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R72'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R73'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R74'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R75'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R76'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R77'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R78'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R79'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R80'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R81'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R82'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R83'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R84'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R85'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R86'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R87'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R88'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R89'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R90'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R91'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R92'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R93'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R94'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R95'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R96'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R97'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R98'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R99'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R100'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R101'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R102'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R103'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R104'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R105'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R106'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R107'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R108'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R109'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R110'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R111'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R112'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R113'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R114'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R115'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R116'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R117'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R118'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R119'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R120'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R121'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R122'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R123'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R124'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R125'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R126'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R127'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R128'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R129'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R130'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R131'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R132'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R133'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R134'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R135'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R136'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R137'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R138'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R139'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R140'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R141'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R142'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R143'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R144'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R145'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R146'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R147'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R148'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R149'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R150'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R151'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R152'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R153'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R154'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R155'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'R156'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL1'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL2'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL3'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL4'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL5'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL6'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL7'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL8'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL9'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL10'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL11'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL12'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL13'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL14'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL15'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL16'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL17'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL18'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL19'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL20'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'CL21'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA01'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA02'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA03'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA04'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA05'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA06'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA07'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA08'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA09'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'SA10'
    i_rxn = i_rxn + 1
    kpp_rxn_labels(i_rxn)%string = 'jo2'

    call assert_msg(313903826, i_rxn.eq.NREACT, "Labeled "//trim(to_string(i_rxn))// &
            " out of "//trim(to_string(NREACT))//" KPP reactions")

  end subroutine get_kpp_rxn_labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compare calculated rates between the modules
  subroutine compare_rates(camp_core, camp_state, ebi_spec_names, conv, &
          ebi_rxn_map, kpp_rxn_map)

    use EXT_RXCM,                               only : NRXNS, RXLABEL
    use EXT_HRDATA,                             only : EBI_PROD => PROD, &
            EBI_LOSS => LOSS, &
            EBI_YC => YC, &
            EBI_RXRAT => RXRAT
    use cb05cl_ae5_Function,                    only : KPP_Fun => Fun, &
            KPP_A => A
    use cb05cl_ae5_Model,                       only : KPP_NVAR => NVAR, &
            KPP_SPC_NAMES => SPC_NAMES, &
            KPP_VAR => VAR, &
            KPP_FIX => FIX, &
            KPP_RCONST => RCONST

    !> CAMP-chem core
    type(camp_core_t), pointer :: camp_core
    !> CAMP-chem state
    type(camp_state_t), pointer :: camp_state
    !> Species map CAMP-chem <-> EBI
    type(string_t), intent(in) :: ebi_spec_names(:)
    !> Coversion factor #/cc -> ppm
    real(kind=dp), intent(in) :: conv
    !> Reaction map EBI <-> camp
    integer(kind=i_kind), allocatable :: ebi_rxn_map(:)
    !> Reaction map KPP <-> camp
    integer(kind=i_kind), allocatable :: kpp_rxn_map(:)

    real(kind=dp) :: KPP_VDOT(KPP_NVAR)
    integer(kind=i_kind) :: i_ebi_spec, i_kpp_spec, i_camp_spec, i_rxn
    type(string_t) :: spec_name

    character(len=:), allocatable :: key

    write(*,*) "Comparing rates"

    call KPP_Fun ( KPP_VAR, KPP_FIX, KPP_RCONST, KPP_VDOT )

    call EXT_HRRATES()
    call EXT_HRPRODLOSS()

    ! Compare the calculated rates
    ! EBI <-> KPP
    do i_ebi_spec = 1, NUM_EBI_SPEC
      do i_kpp_spec = 1, KPP_NVAR
        ! Skip EBI PSSA species
        if (trim(ebi_spec_names(i_ebi_spec)%string).eq.'NO2' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'NO' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'O' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'O3' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'NO3' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'O1D' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'OH' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'HO2' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'N2O5' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'HONO' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'C2O3' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'PAN' .or. &
                trim(ebi_spec_names(i_ebi_spec)%string).eq.'PNA') cycle
        if (trim(ebi_spec_names(i_ebi_spec)%string).eq.trim(KPP_SPC_NAMES(i_kpp_spec))) then
          call warn_assert_msg(616862348, almost_equal(real(KPP_VDOT(i_kpp_spec)*60.0d0*conv, kind=dp), &
                  real(EBI_PROD(i_ebi_spec)-EBI_LOSS(i_ebi_spec), kind=dp), 1.0d-2), &
                  "Species "//ebi_spec_names(i_ebi_spec)%string//" has different time deriv. "// &
                          "EBI solver: "//trim(to_string(real((EBI_PROD(i_ebi_spec)-EBI_LOSS(i_ebi_spec))/60.0d0, kind=dp)))// &
                          "; KPP solver: "//trim(to_string(real(KPP_VDOT(i_kpp_spec)*conv, kind=dp))))
        end if
      end do
    end do

  end subroutine compare_rates

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_args(input_int, id)

    integer, intent(inout) :: input_int
    integer, intent(in), optional :: id
    integer :: status_code
    character(len=50) :: arg

    call get_command_argument(id, arg, status=status_code)

    if (LEN_TRIM(arg) == 0) then
#ifdef PMC_DEBUG
      print*, "Argument ", id, " not present, setting default values..."
#endif
    else
      arg = trim(arg)
      read(arg,*)input_int
    end if

  end subroutine read_args

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_cb05cl_ae5
