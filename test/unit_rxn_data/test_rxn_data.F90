! Copyright (C) 2017 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_test_rxn_data program

!> Test class for the rxn_test_t type, which extends the abstract
!! rxn_data_t type.
program pmc_test_rxn_data

  use pmc_util,                         only: i_kind, dp, assert, &
                                              almost_equal
  use pmc_test_rxn_data_child
  use pmc_phlex_core
  use pmc_phlex_state
  use pmc_chem_spec_data
  use pmc_rxn_data
  use pmc_integration_data
  use pmc_mpi
#ifdef PMC_USE_JSON
  use json_module
#endif
  use iso_c_binding

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

  call pmc_mpi_init()

  if (run_pmc_rxn_data_tests()) then
    write(*,*) "Rxn data tests - PASS"
  else
    write(*,*) "Rxn data tests - FAIL"
  end if

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run all pmc_rxn_data tests
  logical function run_pmc_rxn_data_tests() result(passed)

    passed = build_rxn_data_set_test()

  end function run_pmc_rxn_data_tests

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Build rxn_data set
  logical function build_rxn_data_set_test()

    type(rxn_data_ptr) :: rxn_test(3)
    type(phlex_core_t), pointer :: phlex_core
    type(c_ptr) :: phlex_core_c_ptr
    procedure(integration_data_deriv_func), pointer :: deriv_func_ptr => null()
    procedure(integration_data_jac_func), pointer :: jac_func_ptr => null()
    type(phlex_state_t), pointer :: phlex_state
    type(chem_spec_data_t) :: spec_data
    type(integration_data_t), pointer :: integration_data => null()

#ifdef PMC_USE_JSON
    type(json_core), pointer :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next, j_rxn

    character(kind=json_ck, len=:), allocatable :: str_val, key_name
    character(len=:), allocatable :: spec_name, other_spec_name
    character(len=:), allocatable :: json_string

    integer(kind=i_kind) :: i, j, i_spec, j_spec, i_elem
    real(kind=dp), pointer :: abs_tol(:)
    real(kind=dp) :: test_real

    integer(kind=i_kind) :: n_jac_elem
    integer(kind=i_kind), pointer :: n_jac_col_elem(:), jac_row_ids(:)
    integer(kind=i_kind), allocatable ::  use_jac_elem(:,:)

    real(kind=c_double) :: curr_time
    real(kind=dp), pointer :: deriv(:)
    real(kind=dp), pointer :: jac(:)
    type(c_ptr) :: state_array_c_p, deriv_c_p, jac_c_p, integration_data_c_p

    build_rxn_data_set_test = .false.

    json_string = '{ "pmc-data" : ['//new_line//&
            '{  "name" : "peanut butter", "type" : "GAS_SPEC", "MW" : 12.75 },'//new_line//&
            '{  "name" : "jelly", "type" : "GAS_SPEC", "MW" : 42.5 },'//new_line//&
            '{  "name" : "sandwich", "type" : "GAS_SPEC", "MW" : 13.58 },'//new_line//&
            '{  "name" : "snack", "type" : "GAS_SPEC", "MW" : 185.39 },'//new_line//&
            '{  "name" : "oreo", "type" : "GAS_SPEC", "MW" : 12.45 },'//new_line//&
            '{  "name" : "extra lid", "type" : "GAS_SPEC" },'//new_line//&
            '{  "name" : "lunch", "type" : "GAS_SPEC" },'//new_line//&
            '{  "name" : "my mechanism",'//new_line//&
            '  "type" : "MECHANISM",'//new_line//&
            '  "reactions" : [{'//new_line//&
            '    "rxn type" : "TEST",'//new_line//&
            '    "reactants" : {'//new_line//&
            '      "peanut butter" : {},'//new_line//&
            '      "jelly" : {}'//new_line//&
            '    },'//new_line//&
            '    "products" : {'//new_line//&
            '      "sandwich" : {}'//new_line//&
            '    },'//new_line//&
            '    "A" : 1.2e3,'//new_line//&
            '    "B" : 1.2,'//new_line//&
            '    "theta" : 298.15'//new_line//&
            '  },'//new_line//&
            '  {'//new_line//&
            '    "rxn type" : "TEST",'//new_line//&
            '    "reactants" : {'//new_line//&
            '      "oreo" : { "qty" : 2 }'//new_line//&
            '    },'//new_line//&
            '    "products" : {'//new_line//&
            '      "snack" : {},'//new_line//&
            '      "extra lid" : { "yield" : 2.0 }'//new_line//&
            '    },'//new_line//&
            '    "A" : 3056,'//new_line//&
            '    "B" : 0.9,'//new_line//&
            '    "theta" : 275.0'//new_line//&
            '  },'//new_line//&
            '  {'//new_line//&
            '    "rxn type" : "TEST",'//new_line//&
            '    "reactants" : {'//new_line//&
            '      "snack" : {},'//new_line//&
            '      "sandwich" : {}'//new_line//&
            '    },'//new_line//&
            '    "products" : {'//new_line//&
            '      "lunch" : { "yield" : 0.73 }'//new_line//&
            '    },'//new_line//&
            '    "A" : 5.2e4,'//new_line//&
            '    "B" : 1.0,'//new_line//&
            '    "theta" : 298'//new_line//&
            '  }]'//new_line//&
            '}'//new_line//&
            ']}'//new_line
 
    ! Set up the json core
    allocate(json)
    call j_file%load_from_string(json_string)

    ! Get the pmc-data object
    call j_file%get('pmc-data(1)', j_obj)

    ! Set up the phlex chem state
    phlex_state => phlex_state_t()

    ! Set up the chemical species data 
    spec_data = chem_spec_data_t()
    call assert(782321038, associated(j_obj))

    ! Build rxn list
    rxn_test(1)%val => rxn_test_t()
    rxn_test(2)%val => rxn_test_t()
    rxn_test(3)%val => rxn_test_t()

    do while (associated(j_obj))
      call json%get(j_obj, 'type', str_val)
      if (str_val.eq."GAS_SPEC" .or. str_val.eq."AERO_SPEC") then
        call spec_data%load(json, j_obj)
      else if (str_val.eq."MECHANISM") then
        call json%get(j_obj, 'reactions(1)', j_rxn)
        call assert(767719830, associated(j_rxn))
        i=1
        do while(associated(j_rxn))
          call rxn_test(i)%val%load(json, j_rxn)
          i=i+1
          j_next => j_rxn
          call json%get_next(j_next, j_rxn)
        end do
        call assert(366660985, i.eq.4)
      else
        call assert(476334378, .false.)
      end if
      j_next => j_obj
      call json%get_next(j_next, j_obj)
    end do

    ! Get a species state variable
    phlex_state%env_state%temp = real(301.15, kind=dp)
    allocate(phlex_state%state_var(spec_data%size()))
    call assert(760288463, size(phlex_state%state_var).eq.7)

    ! Initialize the reactions
    do i=1,3
      call rxn_test(i)%val%initialize(spec_data)
    end do

    ! Get the Jacobian elements used
    allocate(use_jac_elem(spec_data%size(), spec_data%size()))
    use_jac_elem(:,:) = 0
    do i=1,3
      call rxn_test(i)%val%get_used_jac_elem(use_jac_elem)
    end do
    n_jac_elem = 0
    do i=1,7
      do j=1,7
        if (use_jac_elem(i,j).gt.0) n_jac_elem = n_jac_elem + 1
      end do
    end do
    call assert(122462114, n_jac_elem.eq.15)
    allocate(n_jac_col_elem(spec_data%size()))
    allocate(jac_row_ids(n_jac_elem))
    n_jac_col_elem(:) = 0
    i_elem = 0
    do j=1,7
      do i=1,7
        if (use_jac_elem(i,j).gt.0) then
          i_elem = i_elem + 1
          n_jac_col_elem(j) = n_jac_col_elem(j) + 1
          jac_row_ids(i_elem) = i
        end if
      end do
    end do
    deallocate(use_jac_elem)
    call assert(269559902, n_jac_col_elem(4).eq.3)
    call assert(318775381, n_jac_col_elem(5).eq.3)
    call assert(990779917, jac_row_ids(6).eq.3)
    call assert(205003797, jac_row_ids(7).eq.3)
    call assert(317322142, jac_row_ids(8).eq.4)

    ! Get an integration data object
    allocate(abs_tol(spec_data%size()))
    abs_tol(:) = 1d-30
    allocate(phlex_core)
    deriv_func_ptr => test_deriv_func
    jac_func_ptr => test_jac_func
    phlex_core_c_ptr = c_loc(phlex_core)
    integration_data => integration_data_t(phlex_core_c_ptr, deriv_func_ptr, &
            jac_func_ptr, abs_tol, n_jac_elem, n_jac_col_elem, jac_row_ids)
    integration_data%phlex_state => phlex_state
    
    ! Update the reaction integration indexes
    do i=1,3
      call rxn_test(i)%val%update_integration_ids(integration_data)
    end do

    ! Initialize the deriv and jac arrays and associate them with the integration
    ! data object
    curr_time = real(0.0, kind=c_double)
    allocate(deriv(spec_data%size()))
    allocate(jac(n_jac_elem))
    state_array_c_p = c_loc(phlex_state%state_var)
    deriv_c_p = c_loc(deriv)
    jac_c_p = c_loc(jac)
    integration_data_c_p = c_loc(integration_data)
    call deriv_func(int(spec_data%size(), kind=c_int), curr_time, &
            state_array_c_p, deriv_c_p, integration_data_c_p)
    call jac_func(int(spec_data%size(), kind=c_int), curr_time, &
            state_array_c_p, jac_c_p, integration_data_c_p)

    ! Set the species concentrations
    phlex_state%state_var(:) = real(0.0, kind=dp)
    spec_name = "peanut butter"
    phlex_state%state_var(spec_data%gas_state_id(spec_name)) = real(200.0, kind=dp)
    spec_name = "jelly"
    phlex_state%state_var(spec_data%gas_state_id(spec_name)) = real(150.0, kind=dp)
    spec_name = "oreo"
    phlex_state%state_var(spec_data%gas_state_id(spec_name)) = real(100.0, kind=dp)
    spec_name = "sandwich"
    phlex_state%state_var(spec_data%gas_state_id(spec_name)) = real(75.0, kind=dp)
    spec_name = "snack"
    phlex_state%state_var(spec_data%gas_state_id(spec_name)) = real(50.0, kind=dp)

    ! Calculate contributions from the test reaction to the time derivative 
    ! and the Jacobian matric
    do i=1, 3
      call rxn_test(i)%val%func_contrib(integration_data)
      call rxn_test(i)%val%jac_contrib(integration_data)
    end do

    ! ******************************************
    ! * Validate time derivative contributions *
    ! ******************************************

    ! Peanut butter is consumed in its reaction with jelly
    test_real = -200.0 * 150.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    spec_name = "peanut butter"
    i_spec = spec_data%gas_state_id(spec_name)
    call assert(445736741, abs(integration_data%get_deriv_elem(i_spec)-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Jelly is consumed in its reaction with peanut butter
    test_real = -200.0/12.75 * 150.0 * 1.2e3 * exp(1.2 + 301.15/298.15)
    spec_name = "jelly"
    i_spec = spec_data%gas_state_id(spec_name)
    call assert(509063946, abs(integration_data%get_deriv_elem(i_spec)-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Sandwiches are mode from peanut butter and jelly and destroyed in the lunch-
    ! forming reaction
    test_real = 200.0/12.75 * 150.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    test_real = test_real - 75.0 * 50.0/185.39 * 5.2e4 * exp(1.0 + 301.15/298)
    spec_name = "sandwich"
    i_spec = spec_data%gas_state_id(spec_name)
    call assert(286332790, abs(integration_data%get_deriv_elem(i_spec)-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Oreos are consumed in a self-reaction (double stuffing)
    test_real = -2.0 * 100.0 * 100.0/12.45 * 3056 * exp(0.9 + 301.15/275.0)
    spec_name = "oreo"
    i_spec = spec_data%gas_state_id(spec_name)
    call assert(398651135, abs(integration_data%get_deriv_elem(i_spec)-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Oreo lids accumulate from the double stuffing process
    test_real = 2.0 * 100.0/12.45 * 100.0/12.45 * 3056 * exp(0.9 + 301.15/275.0)
    spec_name = "extra lid"
    i_spec = spec_data%gas_state_id(spec_name)
    call assert(423701766, abs(integration_data%get_deriv_elem(i_spec)-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Snacks are formed from oreo double stuffing and react with sandwiches to
    ! form a lunch
    test_real = 100.0/12.45 * 100.0/12.45 * 3056 * exp(0.9 + 301.15/275.0)
    test_real = test_real - 50.0 * 75.0/13.58 * 5.2e4 * exp(1.0 + 301.15/298) 
    spec_name = "snack"
    i_spec = spec_data%gas_state_id(spec_name)
    call assert(121073879, abs(integration_data%get_deriv_elem(i_spec)-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! A partial lunch is made from a sandwich and a snack
    test_real = 0.73 * 50.0/185.39 * 75.0/13.58 * 5.2e4 * exp(1.0 + 301.15/298) 
    spec_name = "lunch"
    i_spec = spec_data%gas_state_id(spec_name)
    call assert(277795691, abs(integration_data%get_deriv_elem(i_spec)-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! ******************************************
    ! * Validate Jacobian matrix contributions *
    ! ******************************************

    ! Validate row for sandwich
    spec_name = "sandwich"
    i_spec = spec_data%gas_state_id(spec_name)

    ! Peanut butter contributes to sandwich production
    test_real = 1.0/12.75 * 150.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    other_spec_name = "peanut butter"
    j_spec = spec_data%gas_state_id(other_spec_name)
    call assert(244207961, abs(integration_data%get_jac_elem(i_spec, &
            j_spec)-test_real)/abs(test_real) .lt. 1e-6)

    ! Jelly contributes to sandwich production
    test_real = 200.0/12.75 * 1.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    other_spec_name = "jelly"
    j_spec = spec_data%gas_state_id(other_spec_name)
    call assert(993489873, abs(integration_data%get_jac_elem(i_spec, &
            j_spec)-test_real)/abs(test_real) .lt. 1e-6)

    ! Oreos do not directly affect sandwich loading
    test_real = 0.0 
    other_spec_name = "oreo"
    j_spec = spec_data%gas_state_id(other_spec_name)
    call assert(821879730, integration_data%get_jac_elem(i_spec, &
            j_spec).eq.test_real)

    ! Oreo lids also do not directly affect sandwich loading
    test_real = 0.0 
    other_spec_name = "extra lid"
    j_spec = spec_data%gas_state_id(other_spec_name)
    call assert(529328162, integration_data%get_jac_elem(i_spec, &
            j_spec).eq.test_real)

    ! Snacks react with sandwiches to form a lunch
    test_real = - 1.0/185.39 * 75.0 * 5.2e4 * exp(1.0 + 301.15/298) 
    other_spec_name = "snack"
    j_spec = spec_data%gas_state_id(other_spec_name)
    call assert(121099476, abs(integration_data%get_jac_elem(i_spec, &
            j_spec)-test_real)/abs(test_real) .lt. 1e-6)

    ! Sandwich loading affect the lunch forming reaction
    test_real = - 50.0/185.39 * 1.0 * 5.2e4 * exp(1.0 + 301.15/298) 
    other_spec_name = "sandwich"
    j_spec = spec_data%gas_state_id(other_spec_name)
    call assert(849489332, abs(integration_data%get_jac_elem(i_spec, &
            j_spec)-test_real)/abs(test_real) .lt. 1e-6)

    ! Because we have omitted the reverse (i.e. lunch decomposition) reaction,
    ! sandwich loading is not affected by total lunches
    test_real = 0.0 
    other_spec_name = "lunch"
    j_spec = spec_data%gas_state_id(other_spec_name)
    call assert(272556981, integration_data%get_jac_elem(i_spec, &
            j_spec).eq.test_real)

#endif
    build_rxn_data_set_test = .true.

  end function build_rxn_data_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Time derivative function
  subroutine test_deriv_func(integration_data)

    !> Pointer to integration data
    type(integration_data_t), pointer, intent(inout) :: integration_data

  end subroutine test_deriv_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Time derivative function
  subroutine test_jac_func(integration_data)

    !> Pointer to integration data
    type(integration_data_t), pointer, intent(inout) :: integration_data

  end subroutine test_jac_func

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_rxn_data
