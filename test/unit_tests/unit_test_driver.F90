! Copyright (C) 2019 Christian Guzman and Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The unit test driver program.
program unit_test_driver

  use pmc_camp_core
  use pmc_camp_state
  use pmc_mpi
  use pmc_rand
  use pmc_solver_stats
  use pmc_unit_test_data
  use pmc_util
  use UNIT_TEST_MODULE_

  implicit none

  ! Number of grid cells to solve
  integer(kind=i_kind), parameter :: N_CELLS = 10

  ! Cell to output results for
  integer(kind=i_kind), parameter :: OUTPUT_CELL = 3

  ! File unit for the results output
  integer(kind=i_kind), parameter :: OUTPUT_FILE_UNIT = 7

  ! Size of the state array for one grid cell
  integer(kind=i_kind) :: n_state_var_one_cell

  ! Size of the environmental state array for one grid cell
  integer(kind=i_kind), parameter :: N_ENV_STATE_VAR_ONE_CELL = 2

  ! State pointer for building array of states
  type :: state_ptr
    class(camp_state_t), pointer :: val
  end type state_ptr

  ! Unit test object
  class(unit_test_data_t), pointer :: unit_test
  ! CAMP core for individual and combined solving
  class(camp_core_t), pointer :: one_cell_core, multicell_core

  ! State and related objects
  type(state_ptr) :: grid_cell_state(N_CELLS)
  type(camp_state_t), pointer :: multicell_state, multicell_cell_state
  integer(kind=i_kind) :: grid_cell_state_id(N_CELLS)
  integer(kind=i_kind) :: state_start_idx, state_end_idx
  integer(kind=i_kind) :: env_start_idx, env_end_idx
  real(kind=dp) :: multicell_val, one_cell_val

  ! Solver evaluation object
  type(solver_stats_t) :: solver_stats
  real(kind=dp) :: rel_tol, abs_tol

  integer(kind=i_kind) :: i_cell, i_time, i_spec
  real(kind=dp) :: time
  logical :: passed

#ifdef PMC_USE_MPI
  character, allocatable :: buffer(:), buffer_copy(:)
  integer(kind=i_kind) :: pack_size, pos, i_elem, results
#endif

  ! initialize MPI
  call pmc_mpi_init()

  passed = .true.

  ! Seed the random number generator
  call pmc_srand(0,0)

#ifdef PMC_USE_MPI
  ! Load the model data and initialize the cores on the primary process, then
  ! pass the cores to process 1 for solving
  if( pmc_mpi_rank( ) .eq. 0 ) then
#endif

    ! Create a unit test
    unit_test => UNIT_TEST_TYPE_

    ! Initialize the model for solving individual cells
    one_cell_core => camp_core_t( unit_test%input_file_name( ) )
    call one_cell_core%initialize( )

    ! Initialize the model for solving multiple cells
    multicell_core => camp_core_t( unit_test%input_file_name( ), N_CELLS )
    call multicell_core%initialize( )

    ! Initialize the unit test
    call unit_test%initialize( one_cell_core )

#ifdef PMC_USE_MPI
    ! Pack the cores and the unit test
    pack_size = one_cell_core%pack_size( )  + &
                multicell_core%pack_size( ) + &
                unit_test%pack_size( )
    allocate( buffer( pack_size ) )
    pos = 0
    call one_cell_core%bin_pack(  buffer, pos )
    call multicell_core%bin_pack( buffer, pos )
    call unit_test%bin_pack(      buffer, pos )
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

  ! Upack the objects on process 1
  if( pmc_mpi_rank( ) .eq. 1 ) then

    ! unpack the data
    one_cell_core  => camp_core_t( )
    multicell_core => camp_core_t( )
    unit_test      => UNIT_TEST_TYPE_
    pos = 0
    call one_cell_core%bin_unpack(  buffer, pos )
    call multicell_core%bin_unpack( buffer, pos )
    call unit_test%bin_unpack(      buffer, pos )

    ! Try repacking the data and making sure it stays the same
    allocate( buffer_copy( pack_size ) )
    pos = 0
    call one_cell_core%bin_pack(  buffer_copy, pos )
    call multicell_core%bin_pack( buffer_copy, pos )
    call unit_test%bin_pack(      buffer_copy, pos )
    call assert( 276642139, pos .eq. pack_size )
    do i_elem = 1, pack_size
      call assert_msg( 443440270, buffer( i_elem ) .eq. &
                                  buffer_copy( i_elem ), &
                       "Mismatch in element "//trim( to_string( i_elem ) ) )
    end do
    deallocate( buffer_copy )

    ! Solve the system on process 1
#endif

    ! Open a file to output results
    open( unit=OUTPUT_FILE_UNIT, file="out/"//unit_test%output_file_name( ), &
          status="replace", action="write" )

    ! Initialize the solvers
    call one_cell_core%solver_initialize( )
    call multicell_core%solver_initialize( )

    ! **********************************
    ! *** Set the initial conditions ***
    ! **********************************

    ! Set the initial states (unit_test_data_t classes should provide a certain
    ! number of unique initial states that can be analyzed during solving)
    do i_cell=1, N_CELLS
      grid_cell_state( i_cell )%val => one_cell_core%new_state( )
      grid_cell_state_id( i_cell ) = &
        pmc_rand_int( unit_test%num_unique_states( ) )

      ! Set the initial state for this grid cell to grid_cell_state_id
      call unit_test%initialize_state( i_cell, one_cell_core, &
                                       grid_cell_state( i_cell )%val, &
                                       grid_cell_state_id( i_cell ) )
    end do

    ! Get a multi-cell state
    ! (and one single-cell state to use during analysis)
    multicell_state      => multicell_core%new_state( )
    multicell_cell_state => one_cell_core%new_state( )

    ! Set the size of the state array for one grid cell
    n_state_var_one_cell = size( grid_cell_state( 1 )%val%state_var )

    ! Make sure the multi-cell state is the right size
    call assert( 429906732, size( multicell_state%state_var ) .eq. &
                            n_state_var_one_cell * N_CELLS )
    call assert( 869103793, size( multicell_state%env_var ) .eq. &
                            N_ENV_STATE_VAR_ONE_CELL * N_CELLS )

    ! Set the multicell state using the individual single-cell states
    do i_cell=1, N_CELLS

      ! Chemical state
      state_start_idx = ( i_cell - 1 ) * n_state_var_one_cell + 1
      state_end_idx   = i_cell * n_state_var_one_cell
      multicell_state%state_var( state_start_idx:state_end_idx ) = &
        grid_cell_state( i_cell )%val%state_var( 1:n_state_var_one_cell )

      ! Environmental state
      call multicell_state%env_states( i_cell )%set_temperature_K( &
        grid_cell_state( i_cell )%val%env_states( 1 )%val%temp )
      call multicell_state%env_states( i_cell )%set_pressure_Pa( &
        grid_cell_state( i_cell )%val%env_states( 1 )%val%pressure )

    end do

#ifdef PMC_DEBUG
    ! Evaluate the Jacobian during solving
    solver_stats%eval_Jac = .true.
#endif

    ! *************************************************
    ! *** Solve & analyze results at each time step ***
    ! *************************************************

    ! Output the intial conditions for the output grid cell
    call unit_test%output_results( one_cell_core, &
                                   grid_cell_state( OUTPUT_CELL )%val, &
                                   grid_cell_state_id( OUTPUT_CELL ), &
                                   0, OUTPUT_FILE_UNIT )

    do i_time = 1, unit_test%num_time_steps( )

      ! Solve the multicell system first
      call multicell_core%solve( multicell_state, &
                                 unit_test%time_step_size( ), &
                                 solver_stats = solver_stats )

#ifdef PMC_DEBUG
      ! Check the Jacobian evaluations
      call assert_msg( 673820325, solver_stats%Jac_eval_fails.eq.0, &
                       trim( to_string( solver_stats%Jac_eval_fails ) )// &
                       " Jacobian evaluation failures in multi-cell "// &
                       "solver at time step "//trim( to_string( i_time ) ) )
#endif

      ! Loop over the grid cells to do the single-cell solving and to compare
      ! and analyze the results from multicell and single-cell solving for
      ! each cell
      do i_cell = 1, N_CELLS

        ! Put the multi-cell results in the temporary single-cell state object
        state_start_idx = ( i_cell - 1 ) * n_state_var_one_cell + 1
        state_end_idx   = i_cell * n_state_var_one_cell
        multicell_cell_state%state_var( 1:n_state_var_one_cell ) = &
          multicell_state%state_var( state_start_idx:state_end_idx )
        env_start_idx   = ( i_cell - 1 ) * N_ENV_STATE_VAR_ONE_CELL + 1
        env_end_idx     = i_cell * N_ENV_STATE_VAR_ONE_CELL
        multicell_cell_state%env_var( 1:N_ENV_STATE_VAR_ONE_CELL ) = &
          multicell_state%env_var( env_start_idx:env_end_idx )

        ! Solve the single-cell system
        call one_cell_core%solve( grid_cell_state( i_cell )%val, &
                                  unit_test%time_step_size( ), &
                                  solver_stats = solver_stats )

#ifdef PMC_DEBUG
        ! Check the Jacobian evaluations
        call assert_msg( 781326658, solver_stats%Jac_eval_fails.eq.0, &
                         trim( to_string( solver_stats%Jac_eval_fails ) )// &
                         " Jacobian evaluation failures in single-cell "// &
                         "solver at time step "//trim( to_string( i_time ) ) )
#endif

        ! Make sure the results are similar between the two solvers
        do i_spec = 1, n_state_var_one_cell
          multicell_val = multicell_cell_state%state_var( i_spec )
          one_cell_val  = grid_cell_state( i_cell )%val%state_var( i_spec )
          rel_tol       = one_cell_core%get_rel_tol( ) * i_time * 100.0
          abs_tol       = one_cell_core%get_abs_tol( i_spec ) * i_time * 100.0
          call warn_assert_msg( 294102573, &
                           almost_equal( multicell_val, one_cell_val, &
                                         rel_tol, abs_tol ), &
                           "Different results for single- and multi-cell "// &
                           "solving in cell "//trim( to_string( i_cell ) )// &
                           " for species "//trim( to_string( i_spec ) )// &
                           " at time step "//trim( to_string( i_time ) )// &
                           ". Multicell value: "// &
                           trim( to_string( multicell_val ) )// &
                           ", single-cell value: "// &
                           trim( to_string( one_cell_val ) )//". Relative "//&
                           "tolerance: "//trim( to_string( rel_tol ) )// &
                           ". Absolute tolerance: "// &
                           trim( to_string( abs_tol ) ) )
          passed = passed .and. &
                   almost_equal( multicell_val, one_cell_val, &
                                 rel_tol, abs_tol )
        end do

        ! Output the results
        if( i_cell .eq. OUTPUT_CELL ) &
          call unit_test%output_results( one_cell_core, &
                                         grid_cell_state( i_cell )%val, &
                                         grid_cell_state_id( i_cell ), &
                                         i_time, OUTPUT_FILE_UNIT )

        ! Do the system-specific analysis
        passed = passed .and. &
          unit_test%analyze_state( one_cell_core, &
                                   grid_cell_state( i_cell )%val, &
                                   grid_cell_state_id( i_cell ), &
                                   i_time )
        if( .not.passed ) exit
      end do

      if( .not.passed ) exit
    end do

    ! Close the output file
    close( OUTPUT_FILE_UNIT )

    ! Deallocate states
    do i_cell = 1, N_CELLS
      deallocate( grid_cell_state( i_cell )%val )
    end do
    deallocate( multicell_state )
    deallocate( multicell_cell_state )

#ifdef PMC_USE_MPI
    ! convert the results to an integer
    if( passed ) then
      results = 0
    else
      results = 1
    end if
  end if

  ! Send the results back to the primary process
  call pmc_mpi_transfer_integer(results, results, 1, 0)

  ! Convert the results back to a logical
  if( pmc_mpi_rank( ) .eq. 0 ) then
    if( results .eq. 0 ) then
      passed = .true.
    else
      passed = .false.
    end if
  end if

  deallocate( buffer )
#endif

  ! deallocate the cores and the unit test
  deallocate( unit_test )
  deallocate( one_cell_core )
  deallocate( multicell_core )

  ! Finalize the test
  if( passed ) then
    if( pmc_mpi_rank( ).eq.0 ) write(*,*) "Unit test - PASS"
    call pmc_mpi_finalize( )
  else
    if( pmc_mpi_rank( ).eq.0 ) write(*,*) "Unit test - FAIL"
    call pmc_mpi_finalize( )
    stop 3
  end if

end program unit_test_driver
