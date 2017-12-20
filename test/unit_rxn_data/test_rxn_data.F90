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
  use pmc_model_state
  use pmc_chem_spec_data
  use pmc_chem_spec_state
  use pmc_rxn_data
#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)

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
    type(model_state_t) :: model_state
    type(chem_spec_data_t) :: spec_data

    real(kind=dp), pointer :: func(:)
    real(kind=dp), pointer :: jac_matrix(:,:)

#ifdef PMC_USE_JSON
    type(json_core), pointer :: json
    type(json_file) :: j_file
    type(json_value), pointer :: j_obj, j_next, j_rxn

    character(kind=json_ck, len=:), allocatable :: str_val, key_name
    character(len=:), allocatable :: spec_name, other_spec_name
    character(len=:), allocatable :: json_string

    integer(kind=i_kind) :: i
    real(kind=dp) :: test_real

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
    model_state%env_state%temp = real(301.15, kind=dp)
    allocate(model_state%state_var(spec_data%size()))
    call assert(760288463, size(model_state%state_var).eq.7)

    ! Set up the time derivative and Jacobian matrix arrays
    allocate(func(size(model_state%state_var)))
    allocate(jac_matrix(size(model_state%state_var), &
            size(model_state%state_var)))

    func (:) = real(0.0, kind=dp)
    jac_matrix(:,:) = real(0.0, kind=dp)

    ! Initialize the reactions
    do i=1,3
      call rxn_test(i)%val%initialize(spec_data)
    end do

    ! Set the species concentrations
    model_state%state_var(:) = real(0.0, kind=dp)
    spec_name = "peanut butter"
    model_state%state_var(spec_data%state_id(spec_name)) = real(200.0, kind=dp)
    spec_name = "jelly"
    model_state%state_var(spec_data%state_id(spec_name)) = real(150.0, kind=dp)
    spec_name = "oreo"
    model_state%state_var(spec_data%state_id(spec_name)) = real(100.0, kind=dp)
    spec_name = "sandwich"
    model_state%state_var(spec_data%state_id(spec_name)) = real(75.0, kind=dp)
    spec_name = "snack"
    model_state%state_var(spec_data%state_id(spec_name)) = real(50.0, kind=dp)

    ! Calculate contributions from the test reaction to the time derivative 
    ! and the Jacobian matric
    do i=1, 3
      call rxn_test(i)%val%func_contrib(model_state, func)
      call rxn_test(i)%val%jac_contrib(model_state, jac_matrix)
    end do

    ! ******************************************
    ! * Validate time derivative contributions *
    ! ******************************************

    ! Peanut butter is consumed in its reaction with jelly
    test_real = -200.0 * 150.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    spec_name = "peanut butter"
    call assert(445736741, abs(func(spec_data%state_id(spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Jelly is consumed in its reaction with peanut butter
    test_real = -200.0/12.75 * 150.0 * 1.2e3 * exp(1.2 + 301.15/298.15)
    spec_name = "jelly"
    call assert(509063946, abs(func(spec_data%state_id(spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Sandwiches are mode from peanut butter and jelly and destroyed in the lunch-
    ! forming reaction
    test_real = 200.0/12.75 * 150.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    test_real = test_real - 75.0 * 50.0/185.39 * 5.2e4 * exp(1.0 + 301.15/298)
    spec_name = "sandwich"
    call assert(286332790, abs(func(spec_data%state_id(spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Oreos are consumed in a self-reaction (double stuffing)
    test_real = -2.0 * 100.0 * 100.0/12.45 * 3056 * exp(0.9 + 301.15/275.0)
    spec_name = "oreo"
    call assert(398651135, abs(func(spec_data%state_id(spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Oreo lids accumulate from the double stuffing process
    test_real = 2.0 * 100.0/12.45 * 100.0/12.45 * 3056 * exp(0.9 + 301.15/275.0)
    spec_name = "extra lid"
    call assert(423701766, abs(func(spec_data%state_id(spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Snacks are formed from oreo double stuffing and react with sandwiches to
    ! form a lunch
    test_real = 100.0/12.45 * 100.0/12.45 * 3056 * exp(0.9 + 301.15/275.0)
    test_real = test_real - 50.0 * 75.0/13.58 * 5.2e4 * exp(1.0 + 301.15/298) 
    spec_name = "snack"
    call assert(121073879, abs(func(spec_data%state_id(spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! A partial lunch is made from a sandwich and a snack
    test_real = 0.73 * 50.0/185.39 * 75.0/13.58 * 5.2e4 * exp(1.0 + 301.15/298) 
    spec_name = "lunch"
    call assert(277795691, abs(func(spec_data%state_id(spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! ******************************************
    ! * Validate Jacobian matrix contributions *
    ! ******************************************

    ! Validate row for sandwich
    spec_name = "sandwich"

    ! Peanut butter contributes to sandwich production
    test_real = 1.0/12.75 * 150.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    other_spec_name = "peanut butter"
    call assert(244207961, abs(jac_matrix(spec_data%state_id(spec_name), &
            spec_data%state_id(other_spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Jelly contributes to sandwich production
    test_real = 200.0/12.75 * 1.0/42.5 * 1.2e3 * exp(1.2 + 301.15/298.15)
    other_spec_name = "jelly"
    call assert(993489873, abs(jac_matrix(spec_data%state_id(spec_name), &
            spec_data%state_id(other_spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Oreos do not directly affect sandwich loading
    test_real = 0.0 
    other_spec_name = "oreo"
    call assert(821879730, jac_matrix(spec_data%state_id(spec_name), &
            spec_data%state_id(other_spec_name)).eq.test_real)

    ! Oreo lids also do not directly affect sandwich loading
    test_real = 0.0 
    other_spec_name = "extra lid"
    call assert(529328162, jac_matrix(spec_data%state_id(spec_name), &
            spec_data%state_id(other_spec_name)).eq.test_real)

    ! Snacks react with sandwiches to form a lunch
    test_real = - 1.0/185.39 * 75.0 * 5.2e4 * exp(1.0 + 301.15/298) 
    other_spec_name = "snack"
    call assert(121099476, abs(jac_matrix(spec_data%state_id(spec_name), &
            spec_data%state_id(other_spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Sandwich loading affect the lunch forming reaction
    test_real = - 50.0/185.39 * 1.0 * 5.2e4 * exp(1.0 + 301.15/298) 
    other_spec_name = "sandwich"
    call assert(849489332, abs(jac_matrix(spec_data%state_id(spec_name), &
            spec_data%state_id(other_spec_name))-test_real)/ &
            abs(test_real) .lt. 1e-6)

    ! Because we have omitted the reverse (i.e. lunch decomposition) reaction,
    ! sandwich loading is not affected by total lunches
    test_real = 0.0 
    other_spec_name = "lunch"
    call assert(272556981, jac_matrix(spec_data%state_id(spec_name), &
            spec_data%state_id(other_spec_name)).eq.test_real)

#endif
    build_rxn_data_set_test = .true.

  end function build_rxn_data_set_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program pmc_test_rxn_data
