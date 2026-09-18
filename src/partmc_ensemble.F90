program partmc_ensemble
    
    use pmc_ensemble
    use pmc_util ! including but redundant since used by pmc_ensemble
#ifdef PMC_USE_TCHEM
    type(ensemble_opt_t) :: ensmb_opt
    type(ensemble_state_t) :: ensmb_state
    
    character(len=300) :: ensmb_spec_name
    type(spec_file_t) :: file 
    integer :: i
#endif

#ifdef PMC_USE_TCHEM
    ! Read in ensemble file as the ONLY command line arg
    if (command_argument_count() /= 1) then
        call print_usage()
        call die_msg(739173192, "invalid commandline arguments") !TODO should be a unique integer, copied from from partmc.F90
    end if

    call get_command_argument(1, ensmb_spec_name)

    i = len_trim(ensmb_spec_name)
    if (ensmb_spec_name((i-8):i) /= '.ensemble') then
       call die_msg(710381938, "input filename must end in .ensemble") !TODO should be a unique integer, copied from from partmc.F90
    end if

    call spec_file_open(ensmb_spec_name, file)
    call spec_file_read_ensemble_opt(file, ensmb_opt)
    call spec_file_close(file)

    ! Allocate partmc arrays
    call ensemble_allocate(ensmb_opt, ensmb_state)

    ! Initialize the ensemble, create TChem instance
    call ensemble_init(ensmb_opt, ensmb_state)

    ! Run an ensemble of PartMC-TChem simulations
    call ensemble_run(ensmb_opt, ensmb_state)

    ! Close out the ensemble run
    call ensemble_finalize(ensmb_state)
#else
    ! TODO use a unique ID for the die message, this is taken from run_part.F90:580
    call die_msg(648994111, "cannot run partmc_ensemble, TChem support " & 
        // "not compiled in")
#endif

    contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the usage text to stdout.
  subroutine print_usage()

    write(*,*) 'Usage: partmc_ensemble <ensemble-file>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program partmc_ensemble