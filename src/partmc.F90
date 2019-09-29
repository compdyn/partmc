! Copyright (C) 2007-2012, 2016, 2017 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The partmc program.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \mainpage PartMC Code Documentation
!!
!! \subpage input_format - Input file format description.
!!
!! \subpage output_format - Output file format description.
!!
!! \subpage module_diagram - Diagram of modules and dependencies.
!!
!! \subpage coding_style - Description of code conventions and style.
!!
!! \subpage publications - Publications about PartMC.
!!
!! \subpage camp_chem - Description of the module for incorporation of chemical
!!   mechanisms in PartMC.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \page input_format Input File Format
!! \section ss_spec_file Spec File Format
!!
!! The input file format is plain text:
!!
!! \subpage spec_file_format "Spec File Format"
!!
!! When running PartMC with the command <tt>partmc input.spec</tt> the
!! first line of the <tt>input.spec</tt> file must define the \c
!! run_type with:
!! <pre>
!! run_type &lt;type&gt;
!! </pre>
!! where <tt>&lt;type&gt;</tt> is one of \c particle, \c exact, or \c
!! sectional. This determines the type of run as well as the format of the
!! remainder of the spec file. The rest of the spec file is described by:
!!
!! \subpage input_format_particle "Particle-resolved simulation"
!!
!! \subpage input_format_exact "Exact (analytical) solution"
!!
!! \subpage input_format_sectional "Sectional model simulation"
!!
!! \section ss_json JSON File Format
!!
!! Beginning with the \ref camp_chem "Chemistry Across Multiple Phases (CAMP)",
!! \c json format is used as an alternative input file format:
!!
!! <a href="https://www.json.org">www.json.org</a>
!!
!! Two types of \c json input files are used by \ref camp_chem
!! "camp-chem":
!!
!! \ref input_format_camp_file_list "CAMP file list"
!!
!! \ref input_format_camp_config "CAMP configuration data"
!!
!! Typically, one \ref input_format_camp_file_list "file list" file is
!! used for a PartMC run, which includes paths to multiple \ref
!! input_format_camp_config "configuration" files containing the \ref
!! camp_chem "camp-chem" configuration data. When running stand-alone
!! PartMC, the path to the \ref input_format_camp_file_list "file list"
!! file is included in the spec file. When  using the PartMC library, the path
!! to the \ref input_format_camp_file_list "file list" file can be passed
!! as an argument to the \c pmc_camp_core::camp_core_t constructor.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \page module_diagram Module Diagram
!!
!! \dotfile partmc_modules.gv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \page coding_style Coding Style
!!
!! The code is mainly modern Fortran, with a few parts still clearly
!! showing their Fortran 77 heritage. Fortran 2003 features are used
!! heavily (especially allocatable array features). The code needs to
!! be processed with \c cpp or a compatible pre-processor.
!!
!! \section oo_fortran Object Oriented Fortran
!!
!! Extensive use is made of Fortran 90 derived types. Derived types
!! are named \c my_type_t and are generally defined in modules named
!! \c pmc_my_type within files named \c my_type.F90. Almost all
!! subroutines and function in each \c my_type.F90 file have names of
!! the form \c my_type_*() and take an object of type \c my_type_t
!! (called \c my_type) as the first argument on which to operate.
!!
!! Module names are always the same as the name of the containing
!! file, but prefixed with \c pmc_. Thus the module \c
!! pmc_condense is contained in the file \c condense.F90.
!!
!! \subsection update_2017 2017 Updates
!!
!! Some changes were made beginning with the \ref camp_chem
!! "Chemistry Across Multiple Phases" module to take advantage of Fortran
!! 2003 and 2008 features. In general, modules are \c private and expose
!! functionality through a single `public` derived type and its
!! type-bound procedures. \c public module functions,
!! variables and parameters are discouraged.
!!
!! Derived-type constructors return a pointer to a newly allocated instance
!! of a module's primary derived type, and are typically included in each
!! module. In
!! addition to the primary derived type, a pointer type may be included for
!! building arrays of pointers to the primary derived type. For example, the
!! \c my_type_t type may have an associated \c my_type_ptr type in the same
!! module whose only member is a pointer to a \c my_type_t object named \c
!! val. (This makes it easier to build arrays of mixed extending
!! types of an abstract type.) When abstract types are used, for example \c
!! my_general_data_t, a factory type in a separate module, \c
!! my_general_factory_t, is provided to create instances of extending types
!! by name or from input files. Extending types are located in a \c
!! \\my_general folder in the \c \\src directory.
!!
!! \section mem_manage Memory Management
!!
!! The memory allocation policy to is to always use \c allocatable
!! arrays and to do the allocation in the lowest-level routine
!! possible. Explicit \c allocate() and \c deallocate() statements are
!! discouraged in favor of automatic memory management, where
!! possible.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \page publications PartMC Publications
!!
!!   - J.&nbsp;Ching, N.&nbsp;Riemer, and M.&nbsp;West (2016) Black
!!     carbon mixing state impacts on cloud microphysical properties:
!!     Effects of aerosol plume and environmental conditions,
!!     <i>Journal of Geophysical Research</i> <b>121</b>(10),
!!     5990-6013, DOI: <a
!!     href="http://dx.doi.org/10.1002/2016JD024851">10.1002/2016JD024851</a>.
!!   - J.&nbsp;H.&nbsp;Curtis, M.&nbsp;D.&nbsp;Michelotti,
!!     N.&nbsp;Riemer, M.&nbsp;Heath, and M.&nbsp;West (2016)
!!     Accelerated simulation of stochastic particle removal processes
!!     in particle-resolved aerosol models, <i>Journal of
!!     Computational Physics</i> <b>322</b>, 21-32, DOI: <a
!!     href="http://dx.doi.org/10.1016/j.jcp.2016.06.029">10.1016/j.jcp.2016.06.029</a>.
!!   - R.&nbsp;M.&nbsp;Healy, N.&nbsp;Riemer, J.&nbsp;C.&nbsp;Wenger,
!!     M.&nbsp;Murphy, M.&nbsp;West, L.&nbsp;Poulain,
!!     A.&nbsp;Wiedensohler, I.&nbsp;P.&nbsp;O'Connor,
!!     E.&nbsp;McGillicuddy, J.&nbsp;R.&nbsp;Sodeau, and
!!     G.&nbsp;J.&nbsp;Evans, Single particle diversity and mixing
!!     state measurements, <i>Atmospheric Chemistry and Physics</i>
!!     <b>14</b>, 6289-6299, DOI: <a
!!     href="http://dx.doi.org/10.5194/acp-14-6289-2014">10.5194/acp-14-6289-2014</a>.
!!   - J.&nbsp;Tian, N.&nbsp;Riemer, M.&nbsp;West,
!!     L.&nbsp;Pfaffenberger, H.&nbsp;Schlager, and A.&nbsp;Petzold
!!     (2014) Modeling the evolution of aerosol particles in a ship
!!     plume using PartMC-MOSAIC, <i>Atmospheric Chemistry and
!!     Physics</i> 14, 5327-5347, DOI: <a
!!     href="http://dx.doi.org/10.5194/acp-14-5327-2014">10.5194/acp-14-5327-2014</a>.
!!   - N.&nbsp;Riemer and M.&nbsp;West (2013) Quantifying aerosol mixing
!!     state with entropy and diversity measures, <i>Atmospheric
!!     Chemistry and Physics</i> 13, 11423-11439, DOI: <a
!!     href="http://dx.doi.org/10.5194/acp-13-11423-2013">10.5194/acp-13-11423-2013</a>.
!!   - M.&nbsp;D.&nbsp;Michelotti, M.&nbsp;T.&nbsp;Heath, and
!!     M.&nbsp;West (2013) Binning for efficient stochastic multiscale
!!     particle simulations, <i>Atmospheric Chemistry and Physics</i>
!!     11(4), 1071-1096, DOI: <a
!!     href="http://dx.doi.org/10.1137/130908038">10.1137/130908038</a>.
!!   - J.&nbsp;Ching, N.&nbsp;Riemer, and M.&nbsp;West (2012) Impacts of
!!     black carbon mixing state on black carbon nucleation scavenging:
!!     Insights from a particle-resolved model, <i>Journal of
!!     Geophysical Research</i> 117(D23209), DOI: <a
!!     href="http://dx.doi.org/10.1029/2012JD018269">10.1029/2012JD018269</a>.
!!   - R.&nbsp;E.&nbsp;L.&nbsp;DeVille, N.&nbsp;Riemer, and
!!     M.&nbsp;West (2011) Weighted Flow Algorithms (WFA) for
!!     stochastic particle coagulation, <i>Journal of Computational
!!     Physics</i> 230(23), 8427-8451, DOI: <a
!!     href="http://dx.doi.org/10.1016/j.jcp.2011.07.027">10.1016/j.jcp.2011.07.027</a>.
!!   - R.&nbsp;A.&nbsp;Zaveri, J.&nbsp;C.&nbsp;Barnard,
!!     R.&nbsp;C.&nbsp;Easter, N.&nbsp;Riemer, and M.&nbsp;West (2010)
!!     Particle-resolved simulation of aerosol size, composition,
!!     mixing state, and the associated optical and cloud condensation
!!     nuclei activation properties in an evolving urban plume,
!!     <i>Journal of Geophysical Research</i> 115(D17210), DOI: <a
!!     href="http://dx.doi.org/10.1029/2009JD013616">10.1029/2009JD013616</a>.
!!   - N.&nbsp;Riemer, M.&nbsp;West, R.&nbsp;A.&nbsp;Zaveri, and
!!     R.&nbsp;C.&nbsp;Easter (2010) Estimating black carbon aging
!!     time-scales with a particle-resolved aerosol model, <i>Journal
!!     of Aerosol Science</i> 41(1), 143-158, DOI: <a
!!     href="http://dx.doi.org/10.1016/j.jaerosci.2009.08.009">10.1016/j.jaerosci.2009.08.009</a>
!!   - N.&nbsp;Riemer, M.&nbsp;West, R.&nbsp;A.&nbsp;Zaveri, and
!!     R.&nbsp;C.&nbsp;Easter (2009) Simulating the evolution of soot
!!     mixing state with a particle-resolved aerosol model, <i>Journal
!!     of Geophysical Research</i> 114(D09202), DOI: <a
!!     href="http://dx.doi.org/10.1029/2008JD011073">10.1029/2008JD011073</a>
!!   - R.&nbsp;McGraw, L.&nbsp;Leng, W.&nbsp;Zhu, N.&nbsp;Riemer, and
!!     M.&nbsp;West (2008) Aerosol dynamics using the quadrature
!!     method of moments: Comparing several quadrature schemes with
!!     particle-resolved simulation, <i>Journal of Physics: Conference
!!     Series</i> 125(012020), DOI: <a
!!     href="http://dx.doi.org/10.1088/1742-6596/125/1/012020">10.1088/1742-6596/125/1/012020</a>

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Top level driver.
program partmc

  use pmc_mpi
  use pmc_bin_grid
  use pmc_aero_state
  use pmc_aero_dist
  use pmc_aero_binned
  use pmc_coag_kernel
  use pmc_aero_data
  use pmc_scenario
  use pmc_env_state
  use pmc_run_part
  use pmc_run_exact
  use pmc_run_sect
  use pmc_spec_file
  use pmc_gas_data
  use pmc_gas_state
  use pmc_util
  use pmc_camp_core
#ifdef PMC_USE_SUNDIALS
  use pmc_condense
#endif

  character(len=300) :: spec_name

  call pmc_mpi_init()

  if (pmc_mpi_rank() == 0) then
     ! only the root process accesses the commandline

     if (command_argument_count() /= 1) then
        call print_usage()
        call die_msg(739173192, "invalid commandline arguments")
     end if

     call get_command_argument(1, spec_name)
  end if

  call pmc_mpi_bcast_string(spec_name)
  call partmc_run(spec_name)

  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the usage text to stderr.
  subroutine print_usage()

    write(*,*) 'Usage: partmc <spec-file>'

  end subroutine print_usage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do a PartMC run.
  subroutine partmc_run(spec_name)

    !> Spec filename.
    character(len=*), intent(in) :: spec_name

    type(spec_file_t) :: file
    character(len=100) :: run_type
    integer :: i

    ! check filename (must be "filename.spec")
    i = len_trim(spec_name)
    if (spec_name((i-4):i) /= '.spec') then
       call die_msg(710381938, "input filename must end in .spec")
    end if

    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O
       call spec_file_open(spec_name, file)
       call spec_file_read_string(file, 'run_type', run_type)
    end if

    call pmc_mpi_bcast_string(run_type)
    if (trim(run_type) == 'particle') then
       call partmc_part(file)
    elseif (trim(run_type) == 'exact') then
       call partmc_exact(file)
    elseif (trim(run_type) == 'sectional') then
       call partmc_sect(file)
    else
       call die_msg(719261940, "unknown run_type: " // trim(run_type))
    end if

  end subroutine partmc_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run a Monte Carlo simulation.
  subroutine partmc_part(file)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file

    type(gas_data_t) :: gas_data
    type(gas_state_t) :: gas_state
    type(gas_state_t) :: gas_state_init
    type(aero_data_t) :: aero_data
    type(aero_dist_t) :: aero_dist_init
    type(aero_state_t) :: aero_state
    type(aero_state_t) :: aero_state_init
    type(scenario_t) :: scenario
    type(env_state_t), target  :: env_state
    type(env_state_t) :: env_state_init
    type(run_part_opt_t) :: run_part_opt
    type(camp_core_t), pointer :: camp_core
    integer :: i_repeat, i_group
    integer :: rand_init
    character, allocatable :: buffer(:)
    integer :: buffer_size, max_buffer_size
    integer :: position
    logical :: do_restart, do_init_equilibriate, aero_mode_type_exp_present
    character(len=PMC_MAX_FILENAME_LEN) :: restart_filename
    integer :: dummy_index, dummy_i_repeat
    real(kind=dp) :: dummy_time, dummy_del_t, n_part
    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file
    character(len=PMC_MAX_FILENAME_LEN) :: camp_config_filename

    !> \page input_format_particle Input File Format: Particle-Resolved Simulation
    !!
    !! See \ref spec_file_format for the input file text format.
    !!
    !! A particle-resolved simulation spec file has the parameters:
    !! - \b run_type (string): must be \c particle
    !! - \b output_prefix (string): prefix of the output filenames
    !!   --- see \ref output_format for the full name format
    !! - \b n_repeat (integer): number of repeats
    !! - \b n_part (integer): number of computational particles to
    !!   simulate (actual number used will vary between <tt>n_part /
    !!   2</tt> and <tt>n_part * 2</tt> if \c allow_doubling and \c
    !!   allow_halving are \c yes)
    !! - \b restart (logical): whether to restart the simulation from
    !!   a saved output data file. If \c restart is \c yes, then the
    !!   following parameters must also be provided:
    !!   - \b restart_file (string): name of file from which to load
    !!     restart data, which must be a PartMC output NetCDF file
    !! - \b t_max (real, unit s): total simulation time
    !! - \b del_t (real, unit s): timestep size
    !! - \b t_output (real, unit s): the interval on which to
    !!   output data to disk (see \ref output_format)
    !! - \b t_progress (real, unit s): the interval on which to
    !!   write summary information to the screen while running
    !! - \b do_camp_chem (logical): whether to run <b>CAMP
    !!   </b> (requires JSON and SUNDIALS support to be compiled
    !!   in). If \c do_camp_chem is \c yes, then the following parameters
    !!   must also be provided:
    !!   - \b camp_config (string): name of file containing a list of \b
    !!     camp-chem configuration files. File format should be \ref
    !!     input_format_camp_config
    !! - \b gas_data (string): name of file from which to read the gas
    !!   material data (only provide if \c restart and \c do_camp_chem
    !!   are \c no) --- the file format should be \subpage
    !!   input_format_gas_data
    !! - \b gas_init (string): name of file from which to read the
    !!   initial gas state at the start of the simulation (only
    !!   provide option if \c restart is \c no) --- the file format
    !!   should be \subpage input_format_gas_state
    !! - \b aerosol_data (string): name of file from which to read the
    !!   aerosol material data (only provide if \c restart and do_camp_chem
    !!   are \c no) --- the file format should be \subpage
    !!   input_format_aero_data
    !! - \b do_fractal (logical): whether to consider particles
    !!   as fractal agglomerates. If \c do_fractal is \c no, then all the
    !!   particles are treated as spherical. If \c do_fractal is \c yes,
    !!   then the following parameters must also be provided:
    !!   - \subpage input_format_fractal
    !! - \b aerosol_init (string): filename containing the initial
    !!   aerosol state at the start of the simulation (only provide
    !!   option if \c restart is \c no) --- the file format should
    !!   be \subpage input_format_aero_dist
    !! - \subpage input_format_scenario
    !! - \subpage input_format_env_state
    !! - \b do_coagulation (logical): whether to perform particle
    !!   coagulation. If \c do_coagulation is \c yes, then the
    !!   following parameters must also be provided:
    !!   - \subpage input_format_coag_kernel
    !! - \b do_condensation (logical): whether to perform explicit
    !!   water condensation (requires SUNDIALS support to be compiled
    !!   in; cannot be used simultaneously with MOSAIC). If \c
    !!   do_condensation is \c yes, then the following parameters must
    !!   also be provided:
    !!   - \b do_init_equilibriate (logical): whether to equilibriate
    !!     the water content of each particle before starting the
    !!     simulation
    !! - \b do_mosaic (logical): whether to use the MOSAIC chemistry
    !!   code (requires support to be compiled in; cannot be used
    !!   simultaneously with condensation). If \c do_mosaic is \c
    !!   yes, then the following parameters must also be provided:
    !!   - \b do_optical (logical): whether to compute optical
    !!     properties of the aersol particles for the output files ---
    !!     see output_format_aero_state
    !! - \b do_nucleation (logical): whether to perform particle
    !!   nucleation. If \c do_nucleation is \c yes, then the following
    !!   parameters must also be provided:
    !!   - \subpage input_format_nucleate
    !! - \b rand_init (integer): if greater than zero then use as
    !!   the seed for the random number generator, or if zero then
    !!   generate a random seed for the random number generator ---
    !!   two simulations on the same machine with the same seed
    !!   (greater than 0) will produce identical output
    !! - \b allow_doubling (logical): if \c yes, then whenever the
    !!   number of simulated particles falls below <tt>n_part /
    !!   2</tt>, every particle is duplicated to give better
    !!   statistics
    !! - \b allow_halving (logical): if \c yes, then whenever the
    !!   number of simulated particles rises above <tt>n_part *
    !!   2</tt>, half of the particles are removed (chosen randomly)
    !!   to reduce the computational expense
    !! - \b do_select_weighting (logical): whether to explicitly select
    !!   the weighting scheme. If \c do_select_weighting is \c yes, then the
    !!   following parameters must also be provided:
    !!   - \subpage input_format_weight_type
    !! - \b record_removals (logical): whether to record information
    !!   about aerosol particles removed from the simulation --- see
    !!   \ref output_format_aero_removed
    !! - \b do_parallel (logical): whether to run in parallel mode
    !!   (requires MPI support to be compiled in). If \c do_parallel
    !!   is \c yes, then the following parameters must also be
    !!   provided:
    !!   - \subpage input_format_output
    !!   - \b mix_timescale (real, unit s): timescale on which to mix
    !!     aerosol particle information amongst processes in an
    !!     attempt to keep the aerosol state consistent (the mixing
    !!     rate is inverse to \c mix_timescale)
    !!   - \b gas_average (logical): whether to average the gas state
    !!     amongst processes each timestep, to ensure uniform gas
    !!     concentrations
    !!   - \b env_average (logical): whether to average the
    !!     environment state amongst processes each timestep, to
    !!     ensure a uniform environment
    !!   - \subpage input_format_parallel_coag

    if (pmc_mpi_rank() == 0) then
       ! only the root process does I/O

       call spec_file_read_string(file, 'output_prefix', &
            run_part_opt%output_prefix)
       call spec_file_read_integer(file, 'n_repeat', run_part_opt%n_repeat)
       call spec_file_read_real(file, 'n_part', n_part)
       call spec_file_read_logical(file, 'restart', do_restart)
       if (do_restart) then
          call spec_file_read_string(file, 'restart_file', restart_filename)
       end if

       call spec_file_read_real(file, 't_max', run_part_opt%t_max)
       call spec_file_read_real(file, 'del_t', run_part_opt%del_t)
       call spec_file_read_real(file, 't_output', run_part_opt%t_output)
       call spec_file_read_real(file, 't_progress', run_part_opt%t_progress)

       call spec_file_read_logical(file, 'do_camp_chem', &
               run_part_opt%do_camp_chem)
       if (run_part_opt%do_camp_chem) then
#ifndef PMC_USE_JSON
         call spec_file_die_msg(581685398, file, &
                 'cannot do camp chem, JSON support not compiled in')
#endif
#ifndef PMC_USE_SUNDIALS
         call spec_file_die_msg(905205341, file, &
                 'cannot do camp chem, SUNDIALS support not compiled in')
#endif
         call spec_file_read_string(file, 'camp_config', &
                 camp_config_filename)
         camp_core => camp_core_t(camp_config_filename)
         call camp_core%initialize()
         ! FIXME: Temporary print state of the data
         call camp_core%print()
       end if

       if (do_restart) then
          call input_state(restart_filename, dummy_index, dummy_time, &
               dummy_del_t, dummy_i_repeat, run_part_opt%uuid, aero_data, &
               aero_state_init, gas_data, gas_state_init, env_state_init)
       end if

       if (.not. do_restart) then
          env_state_init%elapsed_time = 0d0

          if (.not. run_part_opt%do_camp_chem) then
            call spec_file_read_string(file, 'gas_data', sub_filename)
            call spec_file_open(sub_filename, sub_file)
            call spec_file_read_gas_data(sub_file, gas_data)
            call spec_file_close(sub_file)
          else
            call gas_data%initialize(camp_core)
          end if

          call spec_file_read_string(file, 'gas_init', sub_filename)
          call spec_file_open(sub_filename, sub_file)
          call spec_file_read_gas_state(sub_file, gas_data, &
               gas_state_init)
          call spec_file_close(sub_file)

          if (.not. run_part_opt%do_camp_chem) then
            call spec_file_read_string(file, 'aerosol_data', sub_filename)
            call spec_file_open(sub_filename, sub_file)
            call spec_file_read_aero_data(sub_file, aero_data)
            call spec_file_close(sub_file)
          else
            call aero_data%initialize(camp_core)
          end if

          call spec_file_read_fractal(file, aero_data%fractal)

          call spec_file_read_string(file, 'aerosol_init', sub_filename)
          call spec_file_open(sub_filename, sub_file)
          call spec_file_read_aero_dist(sub_file, aero_data, aero_dist_init)
          call spec_file_close(sub_file)
       end if

       call spec_file_read_scenario(file, gas_data, aero_data, scenario)
       call spec_file_read_env_state(file, env_state_init)

       call spec_file_read_logical(file, 'do_coagulation', &
            run_part_opt%do_coagulation)
       if (run_part_opt%do_coagulation) then
          call spec_file_read_coag_kernel_type(file, &
               run_part_opt%coag_kernel_type)
       else
          run_part_opt%coag_kernel_type = COAG_KERNEL_TYPE_INVALID
       end if

       call spec_file_read_logical(file, 'do_condensation', &
            run_part_opt%do_condensation)
#ifndef PMC_USE_SUNDIALS
       call assert_msg(121370218, &
            run_part_opt%do_condensation .eqv. .false., &
            "cannot use condensation, SUNDIALS support is not compiled in")
#endif
       if (run_part_opt%do_condensation .and. run_part_opt%do_camp_chem) then
          call spec_file_die_msg(556720748, file, &
              'cannot do condensation with camp chem')
       end if
       if (run_part_opt%do_condensation) then
          call spec_file_read_logical(file, 'do_init_equilibriate', &
               do_init_equilibriate)
       else
          do_init_equilibriate = .false.
       end if

       call spec_file_read_logical(file, 'do_mosaic', run_part_opt%do_mosaic)
       if (run_part_opt%do_mosaic .and. (.not. mosaic_support())) then
          call spec_file_die_msg(230495365, file, &
               'cannot use MOSAIC, support is not compiled in')
       end if
       if (run_part_opt%do_mosaic .and. run_part_opt%do_condensation) then
          call spec_file_die_msg(599877804, file, &
               'cannot use MOSAIC and condensation simultaneously')
       end if
       if (run_part_opt%do_mosaic .and. run_part_opt%do_camp_chem) then
          call spec_file_die_msg(952967581, file, &
               'cannot use MOSAIC with camp chem')
       end if
       if (run_part_opt%do_mosaic) then
          call spec_file_read_logical(file, 'do_optical', &
               run_part_opt%do_optical)
       else
          run_part_opt%do_optical = .false.
       end if

       call spec_file_read_logical(file, 'do_nucleation', &
            run_part_opt%do_nucleation)
       if (run_part_opt%do_nucleation) then
          call spec_file_read_nucleate_type(file, aero_data, &
               run_part_opt%nucleate_type, run_part_opt%nucleate_source)
       else
          run_part_opt%nucleate_type = NUCLEATE_TYPE_INVALID
       end if

       call spec_file_read_integer(file, 'rand_init', rand_init)
       call spec_file_read_logical(file, 'allow_doubling', &
            run_part_opt%allow_doubling)
       call spec_file_read_logical(file, 'allow_halving', &
            run_part_opt%allow_halving)
       if (.not. do_restart) then
          call spec_file_read_logical(file, 'do_select_weighting', &
               run_part_opt%do_select_weighting)
          if (run_part_opt%do_select_weighting) then
             call spec_file_read_aero_state_weighting_type(file, &
                  run_part_opt%weighting_type, run_part_opt%weighting_exponent)
          else
             run_part_opt%weighting_type = AERO_STATE_WEIGHT_NUMMASS_SOURCE
             run_part_opt%weighting_exponent = 0.0d0
          end if
       end if
       call spec_file_read_logical(file, 'record_removals', &
            run_part_opt%record_removals)

       call spec_file_read_logical(file, 'do_parallel', &
            run_part_opt%do_parallel)
       if (run_part_opt%do_parallel) then
#ifndef PMC_USE_MPI
          call spec_file_die_msg(929006383, file, &
               'cannot use parallel mode, support is not compiled in')
#endif
          call spec_file_read_output_type(file, run_part_opt%output_type)
          call spec_file_read_real(file, 'mix_timescale', &
               run_part_opt%mix_timescale)
          call spec_file_read_logical(file, 'gas_average', &
               run_part_opt%gas_average)
          call spec_file_read_logical(file, 'env_average', &
               run_part_opt%env_average)
          call spec_file_read_parallel_coag_type(file, &
               run_part_opt%parallel_coag_type)
       else
          run_part_opt%output_type = OUTPUT_TYPE_SINGLE
          run_part_opt%mix_timescale = 0d0
          run_part_opt%gas_average = .false.
          run_part_opt%env_average = .false.
          run_part_opt%parallel_coag_type = PARALLEL_COAG_TYPE_LOCAL
       end if

       call spec_file_close(file)
    end if

    ! finished reading .spec data, now broadcast data

    ! initialize RNG with random seed for UUID generation
    call pmc_srand(0, pmc_mpi_rank())

    if (.not. do_restart) then
       call uuid4_str(run_part_opt%uuid)
    end if

#ifdef PMC_USE_MPI
    if (pmc_mpi_rank() == 0) then
       ! root process determines size
       max_buffer_size = 0
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_run_part_opt(run_part_opt)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_real(n_part)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_gas_data(gas_data)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_gas_state(gas_state_init)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_aero_data(aero_data)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_aero_dist(aero_dist_init)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_scenario(scenario)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_env_state(env_state_init)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_integer(rand_init)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_logical(do_restart)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_logical(do_init_equilibriate)
       max_buffer_size = max_buffer_size &
            + pmc_mpi_pack_size_aero_state(aero_state_init)
       if (run_part_opt%do_camp_chem) then
         max_buffer_size = max_buffer_size &
              + camp_core%pack_size()
       end if

       allocate(buffer(max_buffer_size))

       position = 0
       call pmc_mpi_pack_run_part_opt(buffer, position, run_part_opt)
       call pmc_mpi_pack_real(buffer, position, n_part)
       call pmc_mpi_pack_gas_data(buffer, position, gas_data)
       call pmc_mpi_pack_gas_state(buffer, position, gas_state_init)
       call pmc_mpi_pack_aero_data(buffer, position, aero_data)
       call pmc_mpi_pack_aero_dist(buffer, position, aero_dist_init)
       call pmc_mpi_pack_scenario(buffer, position, scenario)
       call pmc_mpi_pack_env_state(buffer, position, env_state_init)
       call pmc_mpi_pack_integer(buffer, position, rand_init)
       call pmc_mpi_pack_logical(buffer, position, do_restart)
       call pmc_mpi_pack_logical(buffer, position, do_init_equilibriate)
       call pmc_mpi_pack_aero_state(buffer, position, aero_state_init)
       if (run_part_opt%do_camp_chem) then
         call camp_core%bin_pack(buffer, position)
       end if
       call assert(181905491, position <= max_buffer_size)
       buffer_size = position ! might be less than we allocated
    end if

    ! tell everyone the size
    call pmc_mpi_bcast_integer(buffer_size)

    if (pmc_mpi_rank() /= 0) then
       ! non-root processes allocate space
       allocate(buffer(buffer_size))
    end if

    ! broadcast data to everyone
    call pmc_mpi_bcast_packed(buffer)

    if (pmc_mpi_rank() /= 0) then
       ! non-root processes unpack data
       position = 0
       call pmc_mpi_unpack_run_part_opt(buffer, position, run_part_opt)
       call pmc_mpi_unpack_real(buffer, position, n_part)
       call pmc_mpi_unpack_gas_data(buffer, position, gas_data)
       call pmc_mpi_unpack_gas_state(buffer, position, gas_state_init)
       call pmc_mpi_unpack_aero_data(buffer, position, aero_data)
       call pmc_mpi_unpack_aero_dist(buffer, position, aero_dist_init)
       call pmc_mpi_unpack_scenario(buffer, position, scenario)
       call pmc_mpi_unpack_env_state(buffer, position, env_state_init)
       call pmc_mpi_unpack_integer(buffer, position, rand_init)
       call pmc_mpi_unpack_logical(buffer, position, do_restart)
       call pmc_mpi_unpack_logical(buffer, position, do_init_equilibriate)
       call pmc_mpi_unpack_aero_state(buffer, position, aero_state_init)
       if (run_part_opt%do_camp_chem) then
         ! set up the camp chem core
         camp_core => camp_core_t()
         ! upack the camp chem core
         call camp_core%bin_unpack(buffer, position)
       end if
       call assert(143770146, position == buffer_size)
    end if

    ! free the buffer
    deallocate(buffer)
#endif

    ! initialize the chemistry solver
    if (run_part_opt%do_camp_chem) then
      call camp_core%solver_initialize()
    end if

    ! re-initialize RNG with the given seed
    call pmc_rand_finalize()
    call pmc_srand(rand_init, pmc_mpi_rank())

    call cpu_time(run_part_opt%t_wall_start)

    do i_repeat = 1,run_part_opt%n_repeat
       run_part_opt%i_repeat = i_repeat

       gas_state = gas_state_init
       if (do_restart) then
          aero_state = aero_state_init
          call aero_state_set_n_part_ideal(aero_state, n_part)
       else
          call aero_state_zero(aero_state)
          aero_mode_type_exp_present &
               = aero_dist_contains_aero_mode_type(aero_dist_init, &
               AERO_MODE_TYPE_EXP) &
               .or. scenario_contains_aero_mode_type(scenario, &
               AERO_MODE_TYPE_EXP)
          if (aero_mode_type_exp_present) then
             call warn_msg(245301880, "using flat weighting only due to " &
                  // "presence of exp aerosol mode")
             call aero_state_set_weight(aero_state, aero_data, &
                  AERO_STATE_WEIGHT_FLAT)
          else
             call aero_state_set_weight(aero_state, aero_data, &
                  run_part_opt%weighting_type, run_part_opt%weighting_exponent)
          end if
          call aero_state_set_n_part_ideal(aero_state, n_part)
          call aero_state_add_aero_dist_sample(aero_state, aero_data, &
               aero_dist_init, 1d0, 0d0, run_part_opt%allow_doubling, &
               run_part_opt%allow_halving)
       end if
       env_state = env_state_init
       call scenario_init_env_state(scenario, env_state, &
            env_state_init%elapsed_time)

#ifdef PMC_USE_SUNDIALS
       if (do_init_equilibriate) then
          call condense_equilib_particles(env_state, aero_data, aero_state)
       end if
#endif
       if (run_part_opt%do_camp_chem) then
          call run_part(scenario, env_state, aero_data, aero_state, gas_data, &
               gas_state, run_part_opt, camp_core=camp_core)
       else
          call run_part(scenario, env_state, aero_data, aero_state, gas_data, &
               gas_state, run_part_opt)
       end if
    end do

    call pmc_rand_finalize()

  end subroutine partmc_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run an exact solution simulation.
  subroutine partmc_exact(file)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file

    character(len=100) :: soln_name
    type(aero_data_t) :: aero_data
    type(scenario_t) :: scenario
    type(env_state_t) :: env_state
    type(aero_dist_t) :: aero_dist_init
    type(run_exact_opt_t) :: run_exact_opt
    type(bin_grid_t) :: bin_grid
    type(gas_data_t) :: gas_data
    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file

    !> \page input_format_exact Exact (Analytical) Solution
    !!
    !! The coagulation kernel and initial distribution must be matched
    !! for an exact solution to exist. The valid choices are:
    !!
    !! <table>
    !! <tr><th>Coagulation kernel</th>
    !!     <th>Initial aerosol distribution</th></tr>
    !! <tr><td>Additive</td>
    !!     <td>Single exponential mode</td></tr>
    !! <tr><td>Constant</td>
    !!     <td>Single exponential mode</td></tr>
    !! <tr><td>Zero</td>
    !!     <td>Anything</td></tr>
    !! </table>
    !!
    !! See \ref spec_file_format for the input file text format.
    !!
    !! An exact (analytical) simulation spec file has the parameters:
    !! - \b run_type (string): must be \c exact
    !! - \b output_prefix (string): prefix of the output filenames ---
    !!   the filenames will be of the form \c PREFIX_SSSSSSSS.nc where
    !!   \c SSSSSSSS is is the eight-digit output index (starting at 1
    !!   and incremented each time the state is output)
    !! - \b t_max (real, unit s): total simulation time
    !! - \b t_output (real, unit s): the interval on which to output
    !!   data to disk and to print progress information to the screen
    !!   (see \ref output_format)
    !! - \subpage input_format_diam_bin_grid
    !! - \b gas_data (string): name of file from which to read the
    !!   gas material data --- the file format should be
    !!   \subpage input_format_gas_data
    !! - \b aerosol_data (string): name of file from which to read the
    !!   aerosol material data --- the file format should be
    !!   \subpage input_format_aero_data
    !! - \b do_fractal (logical): whether to consider particles
    !!   as fractal agglomerates. If \c do_fractal is \c no, then all the
    !!   particles are treated as spherical. If \c do_fractal is \c yes,
    !!   then the following parameters must also be provided:
    !!   - \subpage input_format_fractal
    !! - \b aerosol_init (string): filename containing the initial
    !!   aerosol state at the start of the simulation --- the file
    !!   format should be \subpage input_format_aero_dist
    !! - \subpage input_format_scenario
    !! - \subpage input_format_env_state
    !! - \b do_coagulation (logical): whether to perform particle
    !!   coagulation.  If \c do_coagulation is \c yes, then the
    !!   following parameters must also be provided:
    !!   - \subpage input_format_coag_kernel
    !!
    !! Example:
    !! <pre>
    !! run_type exact                  # exact solution
    !! output_prefix additive_exact    # prefix of output files
    !!
    !! t_max 600                       # total simulation time (s)
    !! t_output 60                     # output interval (0 disables) (s)
    !!
    !! n_bin 160                       # number of bins
    !! d_min 1e-8                      # minimum diameter (m)
    !! d_max 1e-3                      # maximum diameter (m)
    !!
    !! gas_data gas_data.dat           # file containing gas data
    !!
    !! aerosol_data aero_data.dat      # file containing aerosol data
    !! do_fractal no                   # whether to do fractal treatment
    !! aerosol_init aero_init_dist.dat # aerosol initial condition file
    !!
    !! temp_profile temp.dat           # temperature profile file
    !! height_profile height.dat       # height profile file
    !! gas_emissions gas_emit.dat      # gas emissions file
    !! gas_background gas_back.dat     # background gas mixing ratios file
    !! aero_emissions aero_emit.dat    # aerosol emissions file
    !! aero_background aero_back.dat   # aerosol background file
    !!
    !! rel_humidity 0.999              # initial relative humidity (1)
    !! pressure 1e5                    # initial pressure (Pa)
    !! latitude 0                      # latitude (degrees, -90 to 90)
    !! longitude 0                     # longitude (degrees, -180 to 180)
    !! altitude 0                      # altitude (m)
    !! start_time 0                    # start time (s since 00:00 UTC)
    !! start_day 1                     # start day of year (UTC)
    !!
    !! do_coagulation yes              # whether to do coagulation (yes/no)
    !! kernel additive                 # Additive coagulation kernel
    !! </pre>

    ! only serial code here
    if (pmc_mpi_rank() /= 0) then
       return
    end if

    call spec_file_read_string(file, 'output_prefix', run_exact_opt%prefix)

    call spec_file_read_real(file, 't_max', run_exact_opt%t_max)
    call spec_file_read_real(file, 't_output', run_exact_opt%t_output)

    call spec_file_read_radius_bin_grid(file, bin_grid)

    call spec_file_read_string(file, 'gas_data', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_data(sub_file, gas_data)
    call spec_file_close(sub_file)

    call spec_file_read_string(file, 'aerosol_data', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_data(sub_file, aero_data)
    call spec_file_close(sub_file)

    call spec_file_read_fractal(file, aero_data%fractal)

    call spec_file_read_string(file, 'aerosol_init', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dist(sub_file, aero_data, aero_dist_init)
    call spec_file_close(sub_file)

    call spec_file_read_scenario(file, gas_data, aero_data, scenario)
    call spec_file_read_env_state(file, env_state)

    call spec_file_read_logical(file, 'do_coagulation', &
         run_exact_opt%do_coagulation)
    if (run_exact_opt%do_coagulation) then
       call spec_file_read_coag_kernel_type(file, &
            run_exact_opt%coag_kernel_type)
    else
       run_exact_opt%coag_kernel_type = COAG_KERNEL_TYPE_INVALID
    end if

    call spec_file_close(file)

    ! finished reading .spec data, now do the run

    call pmc_srand(0, 0)

    call uuid4_str(run_exact_opt%uuid)

    call scenario_init_env_state(scenario, env_state, 0d0)

    call run_exact(bin_grid, scenario, env_state, aero_data, &
         aero_dist_init, gas_data, run_exact_opt)

    call pmc_rand_finalize()

  end subroutine partmc_exact

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run a sectional code simulation.
  subroutine partmc_sect(file)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file

    type(run_sect_opt_t) :: run_sect_opt
    type(aero_data_t) :: aero_data
    type(aero_dist_t) :: aero_dist_init
    type(aero_state_t) :: aero_init
    type(scenario_t) :: scenario
    type(env_state_t) :: env_state
    type(bin_grid_t) :: bin_grid
    type(gas_data_t) :: gas_data
    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    type(spec_file_t) :: sub_file

    !> \page input_format_sectional Sectional Model Simulation
    !!
    !! See \ref spec_file_format for the input file text format.
    !!
    !! A sectional simulation spec file has the parameters:
    !! - \b run_type (string): must be \c sectional
    !! - \b output_prefix (string): prefix of the output filenames ---
    !!   the filenames will be of the form \c PREFIX_SSSSSSSS.nc where
    !!   \c SSSSSSSS is is the eight-digit output index (starting at 1
    !!   and incremented each time the state is output)
    !! - \b del_t (real, unit s): timestep size
    !! - \b t_output (real, unit s): the interval on which to
    !!   output data to disk (see \ref output_format)
    !! - \b t_progress (real, unit s): the interval on which to
    !!   write summary information to the screen while running
    !! - \subpage input_format_diam_bin_grid
    !! - \b gas_data (string): name of file from which to read the
    !!   gas material data --- the file format should be
    !!   \subpage input_format_gas_data
    !! - \b aerosol_data (string): name of file from which to read the
    !!   aerosol material data --- the file format should be
    !!   \subpage input_format_aero_data
    !! - \b do_fractal (logical): whether to consider particles
    !!   as fractal agglomerates. If \c do_fractal is \c no, then all the
    !!   particles are treated as spherical. If \c do_fractal is \c yes,
    !!   then the following parameters must also be provided:
    !!   - \subpage input_format_fractal
    !! - \b aerosol_init (string): filename containing the initial
    !!   aerosol state at the start of the simulation --- the file
    !!   format should be \subpage input_format_aero_dist
    !! - \subpage input_format_scenario
    !! - \subpage input_format_env_state
    !! - \b do_coagulation (logical): whether to perform particle
    !!   coagulation.  If \c do_coagulation is \c yes, then the
    !!   following parameters must also be provided:
    !!   - \subpage input_format_coag_kernel
    !!
    !! Example:
    !! <pre>
    !! run_type sectional              # sectional code run
    !! output_prefix brown_sect        # prefix of output files
    !!
    !! t_max 86400                     # total simulation time (s)
    !! del_t 60                        # timestep (s)
    !! t_output 3600                   # output interval (0 disables) (s)
    !! t_progress 600                  # progress printing interval (0 disables) (s)
    !!
    !! n_bin 220                       # number of bins
    !! d_min 1e-10                     # minimum diameter (m)
    !! d_max 1e-4                      # maximum diameter (m)
    !!
    !! gas_data gas_data.dat           # file containing gas data
    !! aerosol_data aero_data.dat      # file containing aerosol data
    !! do_fractal no                   # whether to do fractal treatment
    !! aerosol_init aero_init_dist.dat # initial aerosol distribution
    !!
    !! temp_profile temp.dat           # temperature profile file
    !! height_profile height.dat       # height profile file
    !! gas_emissions gas_emit.dat      # gas emissions file
    !! gas_background gas_back.dat     # background gas mixing ratios file
    !! aero_emissions aero_emit.dat    # aerosol emissions file
    !! aero_background aero_back.dat   # aerosol background file
    !!
    !! rel_humidity 0.999              # initial relative humidity (1)
    !! pressure 1e5                    # initial pressure (Pa)
    !! latitude 0                      # latitude (degrees_north, -90 to 90)
    !! longitude 0                     # longitude (degrees_east, -180 to 180)
    !! altitude 0                      # altitude (m)
    !! start_time 0                    # start time (s since 00:00 UTC)
    !! start_day 1                     # start day of year (UTC)
    !!
    !! do_coagulation yes              # whether to do coagulation (yes/no)
    !! kernel brown                    # coagulation kernel
    !! </pre>

    ! only serial code here
    if (pmc_mpi_rank() /= 0) then
       return
    end if

    call spec_file_read_string(file, 'output_prefix', run_sect_opt%prefix)

    call spec_file_read_real(file, 't_max', run_sect_opt%t_max)
    call spec_file_read_real(file, 'del_t', run_sect_opt%del_t)
    call spec_file_read_real(file, 't_output', run_sect_opt%t_output)
    call spec_file_read_real(file, 't_progress', run_sect_opt%t_progress)

    call spec_file_read_radius_bin_grid(file, bin_grid)

    call spec_file_read_string(file, 'gas_data', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_gas_data(sub_file, gas_data)
    call spec_file_close(sub_file)

    call spec_file_read_string(file, 'aerosol_data', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_data(sub_file, aero_data)
    call spec_file_close(sub_file)

    call spec_file_read_fractal(file, aero_data%fractal)

    call spec_file_read_string(file, 'aerosol_init', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dist(sub_file, aero_data, aero_dist_init)
    call spec_file_close(sub_file)

    call spec_file_read_scenario(file, gas_data, aero_data, scenario)
    call spec_file_read_env_state(file, env_state)

    call spec_file_read_logical(file, 'do_coagulation', &
         run_sect_opt%do_coagulation)
    if (run_sect_opt%do_coagulation) then
       call spec_file_read_coag_kernel_type(file, &
            run_sect_opt%coag_kernel_type)
    else
       run_sect_opt%coag_kernel_type = COAG_KERNEL_TYPE_INVALID
    end if

    call spec_file_close(file)

    ! finished reading .spec data, now do the run

    call pmc_srand(0, 0)

    call uuid4_str(run_sect_opt%uuid)

    call scenario_init_env_state(scenario, env_state, 0d0)

    call run_sect(bin_grid, gas_data, aero_data, aero_dist_init, scenario, &
         env_state, run_sect_opt)

    call pmc_rand_finalize()

  end subroutine partmc_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program partmc
