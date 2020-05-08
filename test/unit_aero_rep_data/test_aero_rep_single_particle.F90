! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_aero_rep_data program

!> Test class for the aero_rep_data_t extending types
program pmc_test_aero_rep_data

#ifdef PMC_USE_JSON
  use json_module
#endif
#ifdef PMC_USE_MPI
  use mpi
#endif
  use pmc_aero_rep_data
  use pmc_aero_rep_factory
  use pmc_aero_rep_single_particle
  use pmc_mpi
  use pmc_camp_core
  use pmc_camp_state
  use pmc_property
  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal

  use iso_c_binding
  implicit none

  ! Externally set properties
  real(kind=dp), parameter :: PART_NUM_CONC = 1.23e3
  real(kind=dp), parameter :: PART_RADIUS   = 2.43e-7

  !> Interface to c ODE solver and test functions
  interface
    !> Run the c function tests
    integer(kind=c_int) function run_aero_rep_single_particle_c_tests(solver_data, &
        state, env)  bind (c)
      use iso_c_binding
      !> Pointer to the initialized solver data
      type(c_ptr), value :: solver_data
      !> Pointer to the state array
      type(c_ptr), value :: state
      !> Pointer to the environmental state array
      type(c_ptr), value :: env
    end function run_aero_rep_single_particle_c_tests
  end interface

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  !> initialize mpi
  call pmc_mpi_init()

  if (run_pmc_aero_rep_data_tests()) then
    if (pmc_mpi_rank().eq.0) write(*,*) "Aerosol representation tests - PASS"
  else
    if (pmc_mpi_rank().eq.0) write(*,*) "Aerosol representation tests - FAIL"
    stop 3
  end if

  !> finalize mpi
  call pmc_mpi_finalize()

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_aero_rep_data tests
  logical function run_pmc_aero_rep_data_tests() result(passed)

    use pmc_camp_solver_data

    type(camp_solver_data_t), pointer :: camp_solver_data

    camp_solver_data => camp_solver_data_t()

    if (camp_solver_data%is_solver_available()) then
      ! The MPI tests only involve packing and unpacking the aero rep
      ! from a buffer on the primary task
      if (pmc_mpi_rank().eq.0) then
        passed = build_aero_rep_data_set_test()
      else
        passed = .true.
      end if
    else
      call warn_msg(594028423, "No solver available")
      passed = .true.
    end if

    deallocate(camp_solver_data)

  end function run_pmc_aero_rep_data_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build aero_rep_data set
  logical function build_aero_rep_data_set_test()

    type(camp_core_t), pointer :: camp_core
    type(camp_state_t), pointer :: camp_state
    class(aero_rep_data_t), pointer :: aero_rep

#ifdef PMC_USE_JSON

    integer(kind=i_kind) :: i_spec, j_spec, i_rep, i_phase
    type(string_t), allocatable :: rep_names(:)
    character(len=:), allocatable :: rep_name, spec_name, phase_name
    type(string_t), allocatable :: file_list(:), unique_names(:)
#ifdef PMC_USE_MPI
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_data_ptr), allocatable :: aero_rep_passed_data_set(:)
    character, allocatable :: buffer(:)
    integer(kind=i_kind) :: pos, pack_size, i_prop
#endif
    build_aero_rep_data_set_test = .false.

    camp_core => camp_core_t()

    allocate(file_list(1))
    file_list(1)%string = &
            'test_run/unit_aero_rep_data/test_aero_rep_single_particle.json'

    call camp_core%load(file_list)
    call camp_core%initialize()
    camp_state => camp_core%new_state()

    ! Set up the list of aerosol representation names
    ! !!! Add new aero_rep_data_t extending types here !!!
    allocate(rep_names(1))
    rep_names(1)%string = "AERO_REP_SINGLE_PARTICLE"

    ! Loop through all the aerosol representations
    do i_rep = 1, size(rep_names)

      ! Check the aerosol representation getter functions
      rep_name = rep_names(i_rep)%string
      call assert_msg(253854173, &
              camp_core%get_aero_rep(rep_name, aero_rep), rep_name)
      call assert_msg(362813745, associated(aero_rep), rep_name)
      select type (aero_rep)
        type is (aero_rep_single_particle_t)
        class default
          call die_msg(519535557, rep_name)
      end select

      ! Check the unique name functions
      unique_names = aero_rep%unique_names()
      call assert_msg(885541843, allocated(unique_names), rep_name)
      call assert_msg(206819761, size(unique_names).eq.8, rep_name)
      do i_spec = 1, size(unique_names)
        call assert_msg(142263656, aero_rep%spec_state_id(&
                unique_names(i_spec)%string).gt.0, rep_name)
        do j_spec = 1, size(unique_names)
          if (i_spec.eq.j_spec) cycle
          call assert_msg(414662586, aero_rep%spec_state_id(&
                  unique_names(i_spec)%string) .ne. aero_rep%spec_state_id(&
                  unique_names(j_spec)%string), rep_name)
        end do
      end do

      ! Set the species concentrations
      phase_name = "my test phase one"
      spec_name = "species a"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(258227897, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 1.5
      spec_name = "species b"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(418308482, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 2.5
      spec_name = "species c"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(420214016, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 3.5
      phase_name = "my test phase two"
      spec_name = "species c"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(416855243, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 4.5
      spec_name = "species d"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(578389067, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 5.5
      spec_name = "species e"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(147314014, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 6.5
      phase_name = "my last test phase"
      spec_name = "species b"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(401514617, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 7.5
      spec_name = "species e"
      unique_names = aero_rep%unique_names(phase_name = phase_name, &
              spec_name = spec_name)
      i_spec = aero_rep%spec_state_id(unique_names(1)%string)
      call assert_msg(291101806, i_spec.gt.0, rep_name)
      camp_state%state_var(i_spec) = 8.5

    end do

    rep_name = "AERO_REP_BAD_NAME"
    call assert(676257369, .not.camp_core%get_aero_rep(rep_name, aero_rep))
    call assert(453526213, .not.associated(aero_rep))

#ifdef PMC_USE_MPI
    pack_size = 0
    do i_rep = 1, size(rep_names)
      call assert(778520709, &
              camp_core%get_aero_rep(rep_names(i_rep)%string, aero_rep))
      pack_size = pack_size + aero_rep_factory%pack_size(aero_rep, MPI_COMM_WORLD)
    end do
    allocate(buffer(pack_size))
    pos = 0
    do i_rep = 1, size(rep_names)
      call assert(543807700, &
              camp_core%get_aero_rep(rep_names(i_rep)%string, aero_rep))
      call aero_rep_factory%bin_pack(aero_rep, buffer, pos, MPI_COMM_WORLD)
    end do
    allocate(aero_rep_passed_data_set(size(rep_names)))
    pos = 0
    do i_rep = 1, size(rep_names)
      aero_rep_passed_data_set(i_rep)%val => &
              aero_rep_factory%bin_unpack(buffer, pos, MPI_COMM_WORLD)
    end do
    do i_rep = 1, size(rep_names)
      associate (passed_aero_rep => aero_rep_passed_data_set(i_rep)%val)
        call assert(744932979, &
                camp_core%get_aero_rep(rep_names(i_rep)%string, aero_rep))
        call assert(860610097, size(aero_rep%condensed_data_real) .eq. &
                size(passed_aero_rep%condensed_data_real))
        do i_prop = 1, size(aero_rep%condensed_data_real)
          call assert(205686549, aero_rep%condensed_data_real(i_prop).eq. &
                  passed_aero_rep%condensed_data_real(i_prop))
        end do
        call assert(823663594, size(aero_rep%condensed_data_int) .eq. &
                size(passed_aero_rep%condensed_data_int))
        do i_prop = 1, size(aero_rep%condensed_data_int)
          call assert(318457189, aero_rep%condensed_data_int(i_prop).eq. &
                  passed_aero_rep%condensed_data_int(i_prop))
        end do
      end associate
    end do

    aero_rep => null()

    deallocate(buffer)
    deallocate(aero_rep_passed_data_set)
#endif

    ! Evaluate the aerosol representation c functions
    build_aero_rep_data_set_test = eval_c_func(camp_core)

    deallocate(camp_state)
    deallocate(camp_core)

#endif

  end function build_aero_rep_data_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Evaluate the aerosol representation c functions
  logical function eval_c_func(camp_core) result(passed)

    !> CAMP-core
    type(camp_core_t), intent(inout) :: camp_core

    class(aero_rep_data_t), pointer :: aero_rep
    type(camp_state_t), pointer :: camp_state
    integer(kind=i_kind), allocatable :: phase_ids(:)
    character(len=:), allocatable :: rep_name, phase_name
    type(aero_rep_factory_t) :: aero_rep_factory
    type(aero_rep_update_data_single_particle_radius_t) :: update_radius
    type(aero_rep_update_data_single_particle_number_t) :: update_number

    rep_name = "AERO_REP_SINGLE_PARTICLE"

    call assert_msg(264314298, camp_core%get_aero_rep(rep_name, aero_rep), &
                    rep_name)

    select type( aero_rep )
      type is(aero_rep_single_particle_t)
        call camp_core%initialize_update_object( aero_rep, update_radius )
        call camp_core%initialize_update_object( aero_rep, update_number )
      class default
        call die_msg(766425873, "Wrong aero rep type")
    end select

    call camp_core%solver_initialize()

    camp_state => camp_core%new_state()

    camp_state%state_var(:) = 0.0
    call camp_state%env_states(1)%set_temperature_K(  298.0d0 )
    call camp_state%env_states(1)%set_pressure_Pa( 101325.0d0 )

    ! Update external properties
    call update_radius%set_radius( PART_RADIUS )
    call update_number%set_number( 12.3d0 )
    call camp_core%update_data( update_radius )
    call camp_core%update_data( update_number )

    ! Test re-setting number concentration
    call update_number%set_number( PART_NUM_CONC )
    call camp_core%update_data( update_number )

    passed = run_aero_rep_single_particle_c_tests(                           &
                 camp_core%solver_data_gas_aero%solver_c_ptr,                &
                 c_loc(camp_state%state_var),                                &
                 c_loc(camp_state%env_var)                                   &
                 ) .eq. 0

    deallocate(camp_state)

  end function eval_c_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_aero_rep_data
