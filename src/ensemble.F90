! SYNCED 9/16/26 2:13 PM

module pmc_ensemble

    use pmc_tchem_interface
    use pmc_aero_state
    use pmc_gas_state
    use pmc_env_state
    use pmc_aero_data
    use pmc_gas_data
    use pmc_scenario
    use pmc_spec_file 
    use pmc_aero_dist
    use pmc_fractal
    use pmc_output 
    use pmc_util
    use pmc_rand

    type ensemble_opt_t
        !> Number of ensemble members
        integer :: ensemble_size
        !> Gas configuration filename
        character(len=300) :: gas_config_filename
        !> Aerosol configuration filename.
        character(len=300) :: aero_config_filename
        !> Numerical configuration filename.
        character(len=300) :: solver_config_filename
        !> Seed for first ensemble member
        integer :: reference_seed
        !> Per-member target particle count
        real(kind=dp) :: num_particles 
        
        real(kind=dp) :: t_max, del_t, t_output 
        !> Whether the weight classes for each source are specified in inputs.
        logical :: read_aero_weight_classes
        !> Whether to allow doubling of the population.
        logical :: allow_doubling
        !> Whether to allow halving of the population.
        logical :: allow_halving
        !> Name of directory where member spec files are located
        character(len=300) :: member_dir
        character(len=300) :: member_prefix
        !> Prefix of output files
        character(len=300) :: output_prefix 

        ! The following are included to match signature for calling output subroutine, 
        ! for now won't be actually using...
        !> Repeat number of run. 
        integer :: i_repeat
        !> Whether to compute optical properties. 
        logical :: do_optical
        !> Fractal parameters
        type(fractal_t) :: fractal
        !> Whether to record particle removal information.
        logical :: record_removals
        !> Parallel output type.
        integer :: output_type                        
    end type

    type ensemble_state_t 
        !> Ensemble-wide aerosol data
        type(aero_data_t) :: aero_data
        !> Ensemble-wide gas data 
        type(gas_data_t) :: gas_data
        !> Array of member aerosol states.
        type(aero_state_t), allocatable, dimension(:) :: aero_states
        !> Array of member scenario configs.
        type(scenario_t), allocatable, dimension(:) :: scenarios
        !> Array of member environmental states.
        type(env_state_t), allocatable, dimension(:) :: env_states
        !> Array of member gas_states.
        type(gas_state_t), allocatable, dimension(:) :: gas_states
        !> Array of member random seeds
        integer, allocatable, dimension(:) :: seeds
        !> UUID
        character(len=PMC_UUID_LEN) :: uuid
        !> Allocation upper bound for per-member number of particles (required 
        !> by TChem's fixed allocation size)
        integer :: tchem_max_num_particles
        !> Time (in sec) of previous write to output
        real(kind=dp) :: last_output_time
        !> Output index, increments with each call to the output subroutine
        integer :: i_output
    end type 

    contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef PMC_USE_TCHEM
    !> Allocate PartMC ensemble data structures
    subroutine ensemble_allocate(ensmb_opt, ensmb_state)

        type(ensemble_opt_t), intent(in) :: ensmb_opt
        type(ensemble_state_t), intent(inout) :: ensmb_state

        allocate(ensmb_state%aero_states(ensmb_opt%ensemble_size))
        allocate(ensmb_state%gas_states(ensmb_opt%ensemble_size))
        allocate(ensmb_state%env_states(ensmb_opt%ensemble_size))
        allocate(ensmb_state%scenarios(ensmb_opt%ensemble_size))
        allocate(ensmb_state%seeds(ensmb_opt%ensemble_size))

    end subroutine ensemble_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Read in ensemble-wide spec 
    subroutine spec_file_read_ensemble_opt(file, ensmb_opt)
        type(spec_file_t), intent(inout) :: file 
        type(ensemble_opt_t), intent(inout) :: ensmb_opt

        ! functions for loading various types:
        ! - Boolean variables: call spec_file_read_logical 
        ! - Real type: call spec_file_read_real
        ! - Integer type: call spec_file_read_integer
        ! - String type: call spec_file_read_string'
        ! - custom type: call_spec_file_read_[custom_type] (e.g., fractal)

        ! IMPORTANT NOTE: The reading of these attributes is POSITONAL, so 
        ! the order matters in how they are read in (i.e., read order should match 
        ! the order in the .ensemble file)

        call spec_file_read_integer(file, 'ensemble_size', &
            ensmb_opt%ensemble_size)
        call spec_file_read_integer(file, 'reference_seed', &
            ensmb_opt%reference_seed)  ! TODO check if reference seed is != 0 (clock based, dont want)
        
        call spec_file_read_string(file, 'member_dir', &
            ensmb_opt%member_dir)
        call spec_file_read_string(file, 'member_prefix', & ! assume something simple like "m"; meant to signify that the first block of digits in the output filename corresponds to the member index, e.g., out/ens_m0001_....nc
            ensmb_opt%member_prefix)
        call spec_file_read_string(file, 'output_prefix', &
            ensmb_opt%output_prefix)

        call spec_file_read_real(file, 't_max', ensmb_opt%t_max)
        call spec_file_read_real(file, 'del_t', ensmb_opt%del_t)
        call spec_file_read_real(file, 't_output', ensmb_opt%t_output)
        
        call spec_file_read_string(file, 'tchem_gas_config', &
            ensmb_opt%gas_config_filename)
        call spec_file_read_string(file, 'tchem_aero_config', &
            ensmb_opt%aero_config_filename)
        call spec_file_read_string(file, 'tchem_numerics_config', &
            ensmb_opt%solver_config_filename)
        
        call spec_file_read_real(file, 'n_part', ensmb_opt%num_particles)
        ! These should all be false by default
        call spec_file_read_logical(file, 'allow_doubling', &
            ensmb_opt%allow_doubling)
        call spec_file_read_logical(file, 'allow_halving', &
            ensmb_opt%allow_halving)
        call spec_file_read_logical(file, 'read_aero_weight_classes', &
            ensmb_opt%read_aero_weight_classes)
        call spec_file_read_logical(file, 'record_removals', &
            ensmb_opt%record_removals)
        call spec_file_read_logical(file, 'do_optical', &
            ensmb_opt%do_optical)
        call spec_file_read_fractal(file, ensmb_opt%fractal)
        ! TODO coagulation (not called yet)

        ! only need to read in output type from spec file if do_parallel is true
        ensmb_opt%output_type = OUTPUT_TYPE_SINGLE
        ensmb_opt%i_repeat = 1 ! not allowing repeats for ensemble members (yet), technically maybe this should be in ensmb_state since it (could) be allowed to change


    end subroutine spec_file_read_ensemble_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Read in ensemble member spec
    subroutine spec_file_read_ensemble_member(file, scenario, env_state, &
        gas_state, aero_state, gas_data, aero_data, num_particles, &
        read_aero_weight_classes, allow_doubling, allow_halving)

        type(spec_file_t), intent(inout) :: file 
        type(scenario_t), intent(inout) :: scenario 
        type(env_state_t), intent(inout) :: env_state 
        type(gas_state_t), intent(inout) :: gas_state
        type(aero_state_t), intent(inout) :: aero_state
        type(gas_data_t), intent(in) :: gas_data 
        type(aero_data_t), intent(inout) :: aero_data ! inout b/c reading in new mode sources are appended to aero_data
        real(kind=dp), intent(in) :: num_particles
        logical, intent(in) :: read_aero_weight_classes
        logical, intent(in) :: allow_doubling, allow_halving
        
        ! local variables
        !> Per-member aerosol distribution
        type(aero_dist_t) :: aero_dist
        real(kind=dp), parameter :: sample_prop = 1d0
        !> Factor to scale current sample to achieve characteristic sample.
        real(kind=dp), parameter :: characteristic_factor = 1d0
        !> Creation time for new particles (s).
        real(kind=dp), parameter :: create_time = 0d0
        character(len=300) :: sub_filename
        type(spec_file_t) :: sub_file

        ! debugging
        integer :: n_part_add

        ! IMPORTANT NOTE: The reading of these attributes is POSITONAL, so 
        ! the order matters in how they are read in (i.e., read order should match 
        ! the order in the .spec file)
        
        ! Read in gas IC -> gas_state
        call spec_file_read_string(file, 'gas_init', sub_filename)
        call spec_file_open(sub_filename, sub_file)
        call spec_file_read_gas_state(sub_file, gas_data, gas_state)
        call spec_file_close(sub_file)

        ! Read in aerosol IC -> aero_dist
        call spec_file_read_string(file, 'aerosol_init', sub_filename)
        call spec_file_open(sub_filename, sub_file)
        call spec_file_read_aero_dist(sub_file, aero_data, &
            read_aero_weight_classes, aero_dist)
        call spec_file_close(sub_file)

        ! Set the initial aero state and sample the aerosol distribution
        call aero_state_zero(aero_state)
        call aero_state_set_weight(aero_state, aero_data, &
            AERO_STATE_WEIGHT_NUMMASS_SOURCE) ! Using same scheme as partmc-tchem test case 
        call aero_state_set_n_part_ideal(aero_state, num_particles)
        call aero_state_add_aero_dist_sample(aero_state, &
            aero_data, aero_dist, sample_prop, & 
            characteristic_factor, create_time, & 
            allow_doubling, allow_halving)

        ! NOTE: the number of sources is lazily read in per ensemble member, and the number of 
        ! ideal particles for sampling each member's distribution relies on 1 / (n_group*n_class), 
        ! (which are determined by the # of sources). This means that if ensemble members which 
        ! are read in later have more sources than members read in earlier, n_part_ideal will differ
        ! despite the ensemble-wide n_part being fixed/the same for each member. 

        ! Read in profiles for temp/press/height/gas emiss/gas bkgnd/
        ! aero emiss/aero bkgnd/loss -> scenario
        call spec_file_read_scenario(file, gas_data, aero_data, &
         read_aero_weight_classes, scenario)
        ! Read RH/lat/lon/alt/start time/start day -> env_state
        call spec_file_read_env_state(file, env_state)
       
        ! set elapsed time and initialize the env state (initial temp/pressure/height)
        env_state%elapsed_time = 0d0
        call scenario_init_env_state(scenario, env_state, &
            env_state%elapsed_time)

        
    end subroutine spec_file_read_ensemble_member

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine ensemble_member_output_prefix(ensmb_opt, i_member, member_output_prefix)
        type(ensemble_opt_t), intent(in) :: ensmb_opt 
        integer, intent(in) :: i_member 
        character(len=:), allocatable, intent(out) :: member_output_prefix ! use allocatable length to avoid truncation when forming output prefix
        character(len=4) :: idx_char ! NOTE assumes max of 9999 members

        write(idx_char, '(I4.4)') i_member 
        member_output_prefix = trim(ensmb_opt%output_prefix) // "_" // &
            trim(ensmb_opt%member_prefix) // idx_char

    end subroutine ensemble_member_output_prefix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Initialize ensemble of PartMC-TChem runs 
    subroutine ensemble_init(ensmb_opt, ensmb_state)

        type(ensemble_opt_t), intent(in) :: ensmb_opt ! read only!
        type(ensemble_state_t), intent(inout) :: ensmb_state

        ! Local variables
        integer :: i
        character(len=300) :: member_spec_name
        character(len=4) :: member_idx_char ! NOTE assumes max of 9999 members
        type(spec_file_t):: file

        ! Initialize RNG
        call pmc_srand(0, 0)

        ! initialize PartMC-TChem (set up Kokkos instance, populate aerosol_data, 
        ! gas_data)
        call pmc_tchem_initialize(ensmb_opt%gas_config_filename, &
            ensmb_opt%aero_config_filename, ensmb_opt%solver_config_filename, & 
            ensmb_state%gas_data, ensmb_state%aero_data, ensmb_opt%ensemble_size)

        ! set fractal parameters
        ensmb_state%aero_data%fractal = ensmb_opt%fractal

        ! Max number of particles for fixed allocation size dictated by TChem
        ensmb_state%tchem_max_num_particles = TChem_getNumberConcentrationVectorSize()

        ! UUID -- use as a unique identifer for ensemble run
        call uuid4_str(ensmb_state%uuid)
        
        ! Set system state for each ensemble member
        do i = 1, ensmb_opt%ensemble_size 

            ! initialize member seed
            call pmc_rand_finalize()
            ensmb_state%seeds(i) = ensmb_opt%reference_seed + (i - 1)
            call pmc_srand(ensmb_state%seeds(i), 0)

            ! member spec name in member directory, e.g., members/m0001.spec
            write(member_idx_char, '(I0.4)') i ! NOTE assumes max of 9999 members
            member_spec_name = trim(ensmb_opt%member_dir) // "/" // &
                trim(ensmb_opt%member_prefix) // member_idx_char // ".spec"

            call spec_file_open(member_spec_name, file) 
            call spec_file_read_ensemble_member(file, ensmb_state%scenarios(i), &
                ensmb_state%env_states(i), ensmb_state%gas_states(i), &
                ensmb_state%aero_states(i), ensmb_state%gas_data, &
                ensmb_state%aero_data, ensmb_opt%num_particles, &
                ensmb_opt%read_aero_weight_classes, ensmb_opt%allow_doubling, &
                ensmb_opt%allow_halving)
            call spec_file_close(file)
   
        end do 

    end subroutine ensemble_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Run an ensemble of PartMC-TChem box model simulations
    subroutine ensemble_run(ensmb_opt, ensmb_state)

        type(ensemble_opt_t), intent(in) :: ensmb_opt
        type(ensemble_state_t), intent(inout) :: ensmb_state

        ! local variables
        integer :: i
        integer :: i_start, i_end, i_time, n_time
        real(kind=dp) :: elapsed_time
        character(len=:), allocatable :: member_output_prefix

        call check_time_multiple("t_max", ensmb_opt%t_max, &
         "del_t", ensmb_opt%del_t)
        call check_time_multiple("t_output", ensmb_opt%t_output, &
            "del_t", ensmb_opt%del_t)
        
        ensmb_state%i_output = 1
        ensmb_state%last_output_time = 0d0
        ! Save initial condition to output
        elapsed_time = 0d0
        if (ensmb_opt%t_output > 0d0) then
            do i = 1, ensmb_opt%ensemble_size
          
                ! member output file out/ensemble_[gas-mechanism-aero-mechanism]_m0001_....nc
                call ensemble_member_output_prefix(ensmb_opt, i, member_output_prefix)
                ! NOTE: Currently, this writes a new NetCDF file per-member and per output timestep.
                ! This will result in a lot of files at large ensemble sizes, so this should be 
                ! replaced with an approach akin to WRF-PartMC where the ensemble members are all 
                ! written to the same file at a given output time stamp, so we have n_output_times 
                ! total files intead of ensemble_size * n_output_times files.
                call output_state(member_output_prefix, &
                        ensmb_opt%output_type, ensmb_state%aero_data, &
                        ensmb_state%aero_states(i), ensmb_state%gas_data, &
                        ensmb_state%gas_states(i), ensmb_state%env_states(i), &
                        ensmb_state%i_output, elapsed_time, ensmb_opt%del_t, &
                        ensmb_opt%i_repeat, ensmb_opt%record_removals, &
                        ensmb_opt%do_optical, ensmb_state%uuid, ensmb_state%seeds(i))
                call aero_info_array_zero(ensmb_state%aero_states(i)%aero_info_array)
            end do 
        end if

        ! could consolidate do loops if ensmb_opt%t_output > 0d0 check always true
        do i = 1, ensmb_opt%ensemble_size   
            call aero_state_rebalance(ensmb_state%aero_states(i), &
                ensmb_state%aero_data, ensmb_opt%allow_doubling, &
                ensmb_opt%allow_halving, initial_state_warning=.true.)
        end do
        
        ! integer number of timesteps to take
        n_time = nint(ensmb_opt%t_max / ensmb_opt%del_t)
        
        i_start = 1
        i_end = n_time
        ! Primary time integration loop
        do i_time = i_start, i_end 
            call ensemble_timestep(ensmb_opt, ensmb_state, i_time)
        end do

    end subroutine ensemble_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Run chemistry from t to t + del_t for all ensemble members
    subroutine ensemble_timestep(ensmb_opt, ensmb_state, i_time)

        type(ensemble_opt_t), intent(in) :: ensmb_opt
        type(ensemble_state_t), intent(inout) :: ensmb_state
        integer, intent(in) :: i_time

        ! local variables
        integer :: i
        logical :: do_output
        type(env_state_t) :: old_env_state
        real(kind=dp) :: time
        character(len=:), allocatable :: member_output_prefix
        integer :: n_emit, n_dil_in, n_dil_out ! number of emitted/dilute-in/dilute-out particles (not used for anything currently)

        time = i_time * ensmb_opt%del_t

        ! update the system state for each member prior to chemistry integration
        do i = 1, ensmb_opt%ensemble_size
            old_env_state = ensmb_state%env_states(i)
            call scenario_update_env_state(ensmb_state%scenarios(i), ensmb_state%env_states(i), &
                time) ! Note that if adding restart, time arg becomes time + t_start where t_start is the time of restart (either per member array or scalar across all members)
            call scenario_update_gas_state(ensmb_state%scenarios(i), ensmb_opt%del_t, &
                ensmb_state%env_states(i), old_env_state, ensmb_state%gas_data, &
                ensmb_state%gas_states(i))
            call scenario_update_aero_state(ensmb_state%scenarios(i), ensmb_opt%del_t, &
                ensmb_state%env_states(i), old_env_state, ensmb_state%aero_data, &
                ensmb_state%aero_states(i), n_emit, n_dil_in, n_dil_out, & ! these local particle vars are intent(out)--assigned by subroutine
                ensmb_opt%allow_doubling, ensmb_opt%allow_halving)

            ! Map PartMC -> TChem 
            call tchem_from_partmc(ensmb_state%aero_data, ensmb_state%aero_states(i), &
                ensmb_state%gas_data, ensmb_state%gas_states(i), ensmb_state%env_states(i), &
                i - 1)
        end do

        ! Single timestep of TChem chemistry
        call tchem_timestep(ensmb_opt%del_t)

        ! update PartMC data structures and rebalance each member aero state
        do i = 1, ensmb_opt%ensemble_size
            ! Map TChem -> PartMC 
            call tchem_to_partmc(ensmb_state%aero_data, ensmb_state%aero_states(i), &
                ensmb_state%gas_data, ensmb_state%gas_states(i), ensmb_state%env_states(i), &
                i - 1)

            call aero_state_rebalance(ensmb_state%aero_states(i), ensmb_state%aero_data, &
                ensmb_opt%allow_doubling, ensmb_opt%allow_halving, &
                initial_state_warning=.false.) ! seems like this is always false since this is called after the initial state is set?
            
            if (.not. ensmb_opt%record_removals) then 
                ! flush removals array to avoid unnecessary memory accumulation 
                call aero_info_array_zero(ensmb_state%aero_states(i)%aero_info_array)
            end if
            ! Assert that particle count is <= max number of particles
        end do

        if (ensmb_opt%t_output > 0d0) then
            call check_event(time, ensmb_opt%del_t, ensmb_opt%t_output, &
                ensmb_state%last_output_time, do_output)
            if (do_output) then 
                ensmb_state%i_output = ensmb_state%i_output + 1
                do i = 1, ensmb_opt%ensemble_size
                    ! member output file out/ensemble_[gas-mechanism-aero-mechanism]_m0001_....nc
                    call ensemble_member_output_prefix(ensmb_opt, i, member_output_prefix)
                    call output_state(member_output_prefix, &
                            ensmb_opt%output_type, ensmb_state%aero_data, &
                            ensmb_state%aero_states(i), ensmb_state%gas_data, &
                            ensmb_state%gas_states(i), ensmb_state%env_states(i), &
                            ensmb_state%i_output, time, ensmb_opt%del_t, &
                            ensmb_opt%i_repeat, ensmb_opt%record_removals, &
                            ensmb_opt%do_optical, ensmb_state%uuid, ensmb_state%seeds(i))
                    ! zero out the record of removals (if record_removals then contains removals
                    ! between time - t_output and time)
                    call aero_info_array_zero(ensmb_state%aero_states(i)%aero_info_array)
                end do
            end if
        end if

    end subroutine ensemble_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Close out the ensemble run
    subroutine ensemble_finalize(ensmb_state)

        type(ensemble_state_t), intent(inout) :: ensmb_state 
        ! TODO deallocate manually?

        call pmc_rand_finalize()
        call pmc_tchem_cleanup()

    end subroutine ensemble_finalize
#endif

end module pmc_ensemble
