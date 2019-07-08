! Copyright (C) 2015 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> Test program for the PartMC aqueous chemistry mechanism
program test_aq

    use pmc_aq_rxn_file
    use pmc_aq_mech_data
    use pmc_aq_spec_data
    use pmc_aq_integrate_f


    !> Name of file to open.
    character(len=1000) :: filename
    !> Aqueous reaction file.
    type(aq_rxn_file_t) :: file_capram
    !> Aqueous species init conc input file
    type(spec_file_t) :: file_init
    !> Aqueous species data
    type(spec_file_t) :: file_spec_data
    !> Number of reactions
    integer :: num_rxns

    !> Aqueous chemistry mechanism
    !> This will be passed via pointers to the ODE solver, so
    !> it should be set as a target
    type(aq_mech_data_t), target :: aq_mech_data
    !> Aqueous chemistry related species
    !> This will be passed via pointers to the ODE solver, so
    !> it should be set as a target
    type(aq_spec_data_t), target :: aq_spec_data
    !> Aqueous chemistry state
    type(aq_state_t) :: aq_state
    !> Initial environmental state
    type(env_state_t) :: env_state_initial
    !> Final environmental state
    type(env_state_t) :: env_state_final
    !> Integration time (s)
    real(kind=dp) :: del_t

    character(len=1000) :: temp_str
    integer :: ios

    if (command_argument_count() /= 4) then
        write(*,*) "Invalid number of arguments"
        call die_msg(611054666, "invalid commandline arguments")
    end if

    !> get CAPRAM input file name and open it
    call get_command_argument(1, filename)
    call aq_rxn_file_open(filename, file_capram)

    !> Output number of reactions in input file
    call aq_rxn_file_count_rxns(file_capram, num_rxns)
    write(*,*) "Number of reactions: ", num_rxns

    !> Allocate space for mechanism and species data
    call aq_mech_data_allocate(aq_mech_data)
    call aq_spec_data_allocate(aq_spec_data)

    !> Allocate space for env_state variables
    call env_state_allocate(env_state_initial)
    call env_state_allocate(env_state_final)

    !> Set important environment parameters
    !> (Testing)
    env_state_initial%temp = 298.0
    env_state_final%temp = 298.0
    env_state_initial%solar_zenith_angle = 0.7
    env_state_final%solar_zenith_angle = 0.8


    !> Populate mechanism and species variables with input data
    call aq_mech_file_read_data(file_capram, aq_mech_data, aq_spec_data)

    !> Allocate space for the aqueous chemistry state variable
    call aq_state_allocate_size(aq_state, aq_spec_data%n_spec)

    !> Get aq_spec_data file and open it
    call get_command_argument(3, filename)
    call spec_file_open(filename, file_spec_data)

    !> Load aq. chemistry species data
    call aq_spec_file_read_data(file_spec_data, aq_spec_data)

    !> Make sure all the necessary parameters have been read in for the mechanism
    call aq_mech_data_check_const(aq_mech_data, aq_spec_data)

    !> Set integration time (s)
    !> (For SUNDIALS cvRoberts_dns.c comparison)
    ! del_t = 100.0
    call get_command_argument(4, temp_str)
    read(temp_str, '(g30.0)', iostat=ios) del_t

    !> Print species list
    call aq_spec_data_print(aq_spec_data)

    !> Print mechanism
    call aq_mech_data_print(aq_mech_data, aq_spec_data)

    !> Get init conc input file name and open it
    call get_command_argument(2, filename)
    call spec_file_open(filename, file_init)

    !> Load init conc data
    call aq_spec_file_read_aq_state(file_init, aq_spec_data, aq_state)

    !> Particle radius (m)
    aq_state%radius = 1.0e-5

    !> Particle number concentration (#/cc)
    aq_state%n_particle = 70.0

    !> Print initial species concentrations
    write(*,*) "Before Integration"
    call pmc_aq_state_print(aq_state, aq_spec_data)

    call aq_integrate(aq_state, aq_mech_data, aq_spec_data, env_state_initial, &
            env_state_final, del_t)

    !> Print final species concentrations
    write(*,*) " "
    write(*,*) "After Integration"
    call pmc_aq_state_print(aq_state, aq_spec_data)
    write(*,*) " "


end program test_aq
