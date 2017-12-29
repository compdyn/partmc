! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_aero_rep_data program

!> Test class for the aero_rep_data_t extending types
program pmc_test_aero_rep_data

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal
  use pmc_property
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_aero_rep_data
  use pmc_aero_rep_single_particle
#ifdef PMC_USE_JSON
  use json_module
#endif
  use pmc_mpi

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  !> initialize mpi
  call pmc_mpi_init()

  if (run_pmc_aero_rep_data_tests()) then
    write(*,*) "Aerosol rep data tests - PASS"
  else
    write(*,*) "Aerosol rep data tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_aero_rep_data tests
  logical function run_pmc_aero_rep_data_tests() result(passed)

    use pmc_integration_data

    type(integration_data_t), pointer :: integration_data

    integration_data => integration_data_t()

    if (integration_data%is_solver_available()) then
      passed = build_aero_rep_data_set_test()
    else
      call warn_msg(594028423, "No solver available")
      passed = .true.
    end if

  end function run_pmc_aero_rep_data_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build aero_rep_data set
  logical function build_aero_rep_data_set_test()

    type(phlex_core_t), pointer :: phlex_core
    type(phlex_state_t), pointer :: phlex_state
    class(aero_rep_data_t), pointer :: aero_rep

#ifdef PMC_USE_JSON

    integer(kind=i_kind) :: i_rep, i_spec, j_spec, rep_id, i_phase
    integer(kind=i_kind), allocatable :: i_spec_list(:)
    type(string_t), allocatable :: rep_names(:)
    character(len=:), allocatable :: rep_name, spec_name, phase_name
    type(string_t), allocatable :: file_list(:), unique_names(:)
    real(kind=dp) :: surf_conc, v_t, v_p, r_t, jac_true
    real(kind=dp), allocatable :: jac_contrib(:)

    character, allocatable :: buffer(:)
    integer(kind=i_kind) :: pos, pack_size

    build_aero_rep_data_set_test = .false.

    phlex_core => phlex_core_t()

    allocate(file_list(1))
    file_list(1)%string = 'test_run/unit_aero_rep_data/test_aero_rep_data.json'

    call phlex_core%load(file_list)
    call phlex_core%initialize()
    phlex_state => phlex_core%new_state()
    allocate(jac_contrib(size(phlex_state%state_var)))

    ! Set up the list of aerosol representation names
    ! !!! Add new aero_rep_data_t extending types here !!!
    allocate(rep_names(1))
    rep_names(1)%string = "AERO_REP_SINGLE_PARTICLE"

    ! Check the number of aerosol representations in the core
    call assert(154970920, size(phlex_core%aero_rep) .eq. size(rep_names))

    ! Loop through all the aerosol representations
    do i_rep = 1, size(rep_names)
    
      ! Check the aerosol representation getter functions
      rep_name = rep_names(i_rep)%string
      call assert_msg(253854173, phlex_core%find_aero_rep(rep_name, rep_id), rep_name)
      call assert_msg(362813745, rep_id .gt. 0, rep_name)
      call assert_msg(589355969, phlex_core%find_aero_rep(rep_name, aero_rep), rep_name)
      call assert_msg(191203602, associated(aero_rep), rep_name)
      select type (aero_rep)
        type is (aero_rep_single_particle_t)
        class default
          call die_msg(519535557, rep_name)
      end select
      aero_rep => phlex_core%aero_rep(rep_id)%val
      call assert_msg(240871376, associated(aero_rep), rep_name)
      select type (aero_rep)
        type is (aero_rep_single_particle_t)
        class default
          call die_msg(625136356, rep_name)
      end select

      ! Check the unique name functions
      unique_names = aero_rep%unique_names()
      call assert_msg(885541843, allocated(unique_names), rep_name)
      call assert_msg(206819761, size(unique_names).eq.8, rep_name)
      do i_spec = 1, size(unique_names)
        call assert_msg(142263656, aero_rep%state_id_by_unique_name(&
                unique_names(i_spec)%string).gt.0, rep_name)
        do j_spec = 1, size(unique_names)
          if (i_spec.eq.j_spec) cycle
          call assert_msg(414662586, aero_rep%state_id_by_unique_name(&
                  unique_names(i_spec)%string) .ne. aero_rep%state_id_by_unique_name(&
                  unique_names(j_spec)%string), rep_name)
        end do
      end do

      ! Set the species concentrations
      phase_name = "my test phase one"
      spec_name = "species a"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(258227897, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 1.5
      spec_name = "species b"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(418308482, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 2.5
      spec_name = "species c"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(420214016, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 3.5
      phase_name = "my test phase two"
      spec_name = "species c"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(416855243, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 4.5
      spec_name = "species d"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(578389067, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 5.5
      spec_name = "species e"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(147314014, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 6.5
      phase_name = "my last test phase"
      spec_name = "species b"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(401514617, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 7.5
      spec_name = "species e"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(291101806, i_spec.gt.0, rep_name)
      phlex_state%state_var(i_spec) = 8.5

      ! Test the surface area 
      jac_contrib(:) = 0.0
      phase_name = "my test phase two"
      i_phase = aero_rep%phase_id(phase_name)
      call assert_msg(626765157, i_phase.gt.0, rep_name)
      surf_conc = aero_rep%surface_area_conc(0, i_phase, phlex_state)
      call assert_msg(790204515, surf_conc .eq. &
            aero_rep%surface_area_conc(i_phase, 0, phlex_state, jac_contrib) &
            , rep_name)
      v_p = 3.0 * 4.5 + 4.0 * 5.5 + 5.0 * 6.5 
      v_t = 1.0 * 1.5 + 2.0 * 2.5 + 3.0 * 3.5 + v_p + &
            2.0 * 7.5 + 5.0 * 8.5
      r_t = (3.0/4.0 * v_t/const%pi)**(1.0/3.0)
      call assert_msg(528621322, almost_equal(surf_conc, 3.0*v_p/r_t), rep_name)
      phase_name = "my test phase two"
      spec_name = "species c"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      jac_true = (3.0 * r_t * 3.0 - 3.0 * v_p * &
              1.0/3.0 * (3.0/4.0 * v_t/const%pi)**(1.0/3.0) * 3.0)/r_t**2
      call assert_msg(593916845, jac_true .eq. jac_contrib(i_spec), rep_name)
      phase_name = "my last test phase"
      spec_name = "species e"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      jac_true =  (-3.0 * v_p * &
              1.0/3.0 * (3.0/4.0 * v_t/const%pi)**(1.0/3.0) * 5.0)/r_t**2
      call assert_msg(638320312, jac_true .eq. jac_contrib(i_spec), rep_name)

      ! Test the species surface area 
      jac_contrib(:) = 0.0
      phase_name = "my test phase two"
      spec_name = "species c"
      i_phase = aero_rep%phase_id(phase_name)
      i_spec_list = aero_rep%aero_phase(i_phase)%val%spec_id(spec_name)
      i_spec = i_spec_list(1)
      call assert_msg(136045982, i_phase.gt.0, rep_name)
      surf_conc = aero_rep%species_surface_area_conc(i_phase, i_spec, phlex_state)
      call assert_msg(530839576, surf_conc .eq. &
            aero_rep%species_surface_area_conc(i_phase, i_spec, phlex_state, &
            jac_contrib), rep_name)
      v_p = 3.0 * 4.5
      v_t = 1.0 * 1.5 + 2.0 * 2.5 + 3.0 * 3.5 + &
            v_p + 4.0 * 5.5 + 5.0 * 6.5 + &  
            2.0 * 7.5 + 5.0 * 8.5
      r_t = (3.0/4.0 * v_t/const%pi)**(1.0/3.0)
      call assert_msg(299485340, surf_conc .eq. 3.0*v_p/r_t, rep_name)
      phase_name = "my test phase two"
      spec_name = "species c"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      jac_true = (3.0 * r_t * 3.0 - 3.0 * v_p * &
              1.0/3.0 * (3.0/4.0 * v_t/const%pi)**(1.0/3.0) * 3.0)/r_t**2
      call assert_msg(413709219, jac_true .eq. jac_contrib(i_spec), rep_name)
      phase_name = "my last test phase"
      spec_name = "species e"
      i_spec_list = aero_rep%species_state_id(phase_name, spec_name)
      i_spec = i_spec_list(1)
      jac_true =  (-3.0 * v_p * &
              1.0/3.0 * (3.0/4.0 * v_t/const%pi)**(1.0/3.0) * 5.0)/r_t**2
      call assert_msg(190978063, jac_true .eq. jac_contrib(i_spec), rep_name)


    end do

    rep_name = "AERO_REP_BAD_NAME"
    call assert(676257369, .not.phlex_core%find_aero_rep(rep_name, rep_id))
    call assert(453526213, rep_id .eq. 0)
    call assert(848319807, .not.phlex_core%find_aero_rep(rep_name, aero_rep))
    call assert(343113402, .not.associated(aero_rep))



#ifdef PMC_USE_MPI
#endif

    ! If condensed data arrays are used for aerosol phases in the future, put
    ! tests for passed info here

#endif  
    build_aero_rep_data_set_test = .true.

  end function build_aero_rep_data_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_aero_rep_data
