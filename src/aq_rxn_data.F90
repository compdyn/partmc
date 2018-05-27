! Copyright (C) 2015 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_rxn_data module.

!> The aq_rxn_data_t structure and associated subroutines.

module pmc_aq_rxn_data

  use pmc_aq_spec_data
  use pmc_constants
#ifdef PMC_USE_MPI
  use mpi
#endif

  implicit none

  integer, parameter :: AQ_RXN_DATA_ID_MAX_LEN = 100
  integer, parameter :: AQ_RXN_DATA_NUM_CLASS = 3
  integer, parameter :: AQ_RXN_DATA_NUM_RATE = 4
  integer, parameter :: AQ_RXN_STRING_MAX_LEN = 1000

  !> Class names
  character(len=AQ_RXN_DATA_ID_MAX_LEN), parameter, dimension(AQ_RXN_DATA_NUM_CLASS) :: &
                                                AQ_RXN_DATA_CLASS_NAME = &
                                                            (/"HENRY       ", &
                                                              "AQUA        ",  &
                                                              "DISS        "/)
  !> Rate names
  character(len=AQ_RXN_DATA_ID_MAX_LEN), parameter, dimension(AQ_RXN_DATA_NUM_RATE) :: &
                                                AQ_RXN_DATA_RATE_NAME = &
                                                            (/"TEMP3:      ",   &
                                                              "PHOTABC:    ", &
                                                              "DTEMP:      ",   &
                                                              "DCONST:     "/)
  !> Number of rate parameters for each type
  integer, parameter , dimension(AQ_RXN_DATA_NUM_RATE):: AQ_RXN_DATA_RATE_N_PARAM = &
                                                                (/ 2, & ! TEMP3
                                                                   3, & ! PHOTABC
                                                                   3, & ! DTEMP
                                                                   2/)  ! DCONST

  !> Allow product yield (only when a backward reaction is not included in the rate constant)
  !! 0 = allow yields; 1 = do not allow yields
  integer, parameter , dimension(AQ_RXN_DATA_NUM_RATE):: AQ_RXN_DATA_RATE_YIELD_OK = &
                                                                (/ 0, & ! TEMP3
                                                                   0, & ! PHOTABC
                                                                   1, & ! DTEMP
                                                                   1/)  ! DCONST


  !> Reaction data
  !! 
  !! Constant data for a single reaction
  type aq_rxn_data_t
     !> Class Index
     integer :: class_index
     !> Reactant species index in aq_spec_t structure
     integer, pointer :: reactant(:)
     !> Product species index in aq_spec_t structure
     integer, pointer :: product(:)
     !> Product yields (unitless)
     real(kind=dp), pointer :: prod_yield(:)
     !> Rate constant type index
     integer :: rate_index
     !> Rate constant parameters (must correspond to number of
     !! parameters required by rate type)
     real(kind=dp), pointer :: rate_param(:)
  end type aq_rxn_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for reaction.
  subroutine aq_rxn_data_allocate(aq_rxn_data)

    !> Reaction data.
    type(aq_rxn_data_t), intent(out) :: aq_rxn_data

    allocate(aq_rxn_data%reactant(0))
    allocate(aq_rxn_data%product(0))
    allocate(aq_rxn_data%prod_yield(0))
    allocate(aq_rxn_data%rate_param(0))

  end subroutine aq_rxn_data_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for reaction with the given dimensions.
  subroutine aq_rxn_data_allocate_size(aq_rxn_data, n_reactant, n_product, n_rate_param)

    !> Reaction data.
    type(aq_rxn_data_t), intent(out) :: aq_rxn_data
    !> Number of reactants.
    integer, intent(in) :: n_reactant
    !> Number of products.
    integer, intent(in) :: n_product
    !> Number of rate parameters.
    integer, intent(in) :: n_rate_param

    allocate(aq_rxn_data%reactant(n_reactant))
    allocate(aq_rxn_data%product(n_product))
    allocate(aq_rxn_data%prod_yield(n_product))
    allocate(aq_rxn_data%rate_param(n_rate_param))

  end subroutine aq_rxn_data_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aq_rxn_data_deallocate(aq_rxn_data)

    !> Reaction data.
    type(aq_rxn_data_t), intent(inout) :: aq_rxn_data

    deallocate(aq_rxn_data%reactant)
    deallocate(aq_rxn_data%product)
    deallocate(aq_rxn_data%prod_yield)
    deallocate(aq_rxn_data%rate_param)

  end subroutine aq_rxn_data_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get class index
  integer function aq_rxn_data_get_class_index(class_name)

    !> Class name
    character(len=*) :: class_name

    integer :: i

    aq_rxn_data_get_class_index = 0

    do i=1,AQ_RXN_DATA_NUM_CLASS
        if (trim(class_name) .eq. trim(AQ_RXN_DATA_CLASS_NAME(i))) then
            aq_rxn_data_get_class_index = i
        endif
    enddo

  end function aq_rxn_data_get_class_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get class name
  character(len=AQ_RXN_DATA_ID_MAX_LEN) function aq_rxn_data_get_class_name(class_index)

    !> Class index
    integer :: class_index

    if (class_index .gt. 0 .and. class_index .le. AQ_RXN_DATA_NUM_CLASS) then
        aq_rxn_data_get_class_name = AQ_RXN_DATA_CLASS_NAME(class_index)
    else
        aq_rxn_data_get_class_name = ""
    endif

  end function aq_rxn_data_get_class_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get rate index
  integer function aq_rxn_data_get_rate_index(rate_name)

    !> Rate name
    character(len=*) :: rate_name

    integer :: i

    aq_rxn_data_get_rate_index = 0

    do i=1,AQ_RXN_DATA_NUM_RATE
        if (trim(rate_name) .eq. trim(AQ_RXN_DATA_RATE_NAME(i))) then
            aq_rxn_data_get_rate_index = i
        endif
    enddo

  end function aq_rxn_data_get_rate_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get rate name
  character(len=AQ_RXN_DATA_ID_MAX_LEN) function aq_rxn_data_get_rate_name(rate_index)

    !> Rate index
    integer :: rate_index

    if (rate_index .gt. 0 .and. rate_index .le. AQ_RXN_DATA_NUM_RATE) then
        aq_rxn_data_get_rate_name = AQ_RXN_DATA_RATE_NAME(rate_index)
    else
        aq_rxn_data_get_rate_name = ""
    endif

  end function aq_rxn_data_get_rate_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get number of rate parameters for a particular type
  integer function aq_rxn_data_get_num_param(rate_index)

    !> Rate index
    integer, intent(in) :: rate_index

    if (rate_index .gt. 0 .and. rate_index .le. AQ_RXN_DATA_NUM_RATE) then
        aq_rxn_data_get_num_param = AQ_RXN_DATA_RATE_N_PARAM(rate_index)
    else
        aq_rxn_data_get_num_param = 0
    endif

  end function aq_rxn_data_get_num_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine if a rxn includes a backwards rate
  logical function aq_rxn_data_is_backward_rxn(aq_rxn_data)

    !> Reaction data
    type(aq_rxn_data_t), intent(in) :: aq_rxn_data

    aq_rxn_data_is_backward_rxn = .false.

    select case(trim(aq_rxn_data_get_rate_name(aq_rxn_data%rate_index)))
        case("DTEMP:")
            aq_rxn_data_is_backward_rxn = .true.
        case("DCONST:")
            aq_rxn_data_is_backward_rxn = .true.
    end select

    select case(trim(aq_rxn_data_get_class_name(aq_rxn_data%class_index)))
        case("HENRY")
            aq_rxn_data_is_backward_rxn = .true.
        case("DISS")
            aq_rxn_data_is_backward_rxn = .true.
    end select

  end function aq_rxn_data_is_backward_rxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get reaction rate constant for a particular reaction
  integer function aq_rxn_data_get_rate_constant(rc_forward, rc_backward, &
                    aq_rxn_data, aq_spec_data, temp, solar_zenith_angle, radius)
    !> Forward reaction rate constant
    real(kind=dp), intent(out) :: rc_forward
    !> Backward reaction rate constant
    real(kind=dp), intent(out) :: rc_backward
    !> Reaction data
    type(aq_rxn_data_t), intent(in) :: aq_rxn_data
    !> Aq. chemistry species data
    type(aq_spec_data_t), intent(in) :: aq_spec_data
    !> Temperature (K)
    real(kind=dp), intent(in) :: temp
    !> Solar Zenith Angle (radians) 
    real(kind=dp), intent(in) :: solar_zenith_angle
    !> Particle radius (m)
    real(kind=dp), intent(in) :: radius

    ! Rate parameters
    real(kind=dp) :: A, B, C

    ! parameters for phase transfer reactions
    real(kind=dp) :: k_mt, alpha
    integer :: species_index

    rc_forward = 0.0
    rc_backward = 0.0
    aq_rxn_data_get_rate_constant = 0

    select case(trim(aq_rxn_data_get_rate_name(aq_rxn_data%rate_index)))

        case("TEMP3:")
            A = aq_rxn_data%rate_param(1)
            B = aq_rxn_data%rate_param(2)
            rc_forward = A*exp(B*(1.0/temp - 1.0/298.0))

        case("PHOTABC:")
            A = aq_rxn_data%rate_param(1)
            B = aq_rxn_data%rate_param(2)
            C = aq_rxn_data%rate_param(3)
            rc_forward = A*exp(B*(1.0-1.0/(cos(C*solar_zenith_angle))))
            ! Use the same cutoff solar zenith angle as in MOSAIC
            if (solar_zenith_angle.gt.1.55334) rc_forward = 0.0

        case("DTEMP:")
            A = aq_rxn_data%rate_param(1)
            B = aq_rxn_data%rate_param(2)
            C = aq_rxn_data%rate_param(3)
            rc_forward = A*exp(B*(1.0/temp - 1.0/298.0))
            rc_backward = C

        case("DCONST:")
            A = aq_rxn_data%rate_param(1)
            B = aq_rxn_data%rate_param(2)
            rc_forward = A
            rc_backward = B

        case default

            aq_rxn_data_get_rate_constant = 1

    end select

    select case(trim(aq_rxn_data_get_class_name(aq_rxn_data%class_index)))

        ! HENRY
        case("HENRY")
            ! If this is a Henry's Law phase transfer,
            ! rc_forward currently holds the equilibrium constant

            ! calculate the uptake rate constant (including mass accomodation coefficient)
            species_index = aq_rxn_data%reactant(1)
            alpha = aq_rxn_data_get_mass_accom(aq_spec_data%N_star(species_index), temp)
            k_mt = aq_rxn_data_get_mass_transfer_rc(aq_spec_data%Dg(species_index), &
                    radius, alpha, aq_spec_data%MW(species_index), temp)

            ! At equilibrium the evaporation rate equals the condensation rate
            ! If rc_forward is infinite, then there is no evaportation
            if (rc_forward .eq. rc_forward-1) then
                rc_backward = 0.0
                rc_forward = k_mt
            else
                rc_backward = k_mt / rc_forward
                rc_forward = k_mt
            endif

        ! DISS
        case("DISS")
            ! If this is a dissociation reaction,
            ! rc_forward currently holds the equilibrium constant, and
            ! rc_backward holds the backwards rate constant

            rc_forward = rc_backward * rc_forward

    end select

  end function aq_rxn_data_get_rate_constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a reaction to a string.
  character(len=AQ_RXN_STRING_MAX_LEN) function aq_rxn_to_string(val, &
       aq_spec_data)

    !> Value to convert.
    type(aq_rxn_data_t), intent(in) :: val
    !> Aqueous chemistry related species data.
    type(aq_spec_data_t), intent(in) :: aq_spec_data

    character(len=AQ_RXN_STRING_MAX_LEN) :: ret_val
    character(len=PMC_UTIL_CONVERT_STRING_LEN) :: yield
    integer :: i

    ret_val = ""
    ret_val = trim(ret_val) &
         // trim(aq_rxn_data_get_class_name(val%class_index)) // ":"

    do i = 1,size(val%reactant)
       if (i > 1) then
          ret_val = trim(ret_val) // " +"
       end if
       ret_val = trim(ret_val) // " " &
            // trim(aq_spec_data%name(val%reactant(i))) &
            // "(" // trim(integer_to_string(val%reactant(i))) // ")"
    end do

    ret_val = trim(ret_val) // " ==>"

    do i = 1,size(val%product)
       if (i > 1) then
          ret_val = trim(ret_val) // " +"
       end if
       if (val%prod_yield(i) == 1d0) then
          yield = ""
       else
          yield = trim(real_or_integer_to_string(val%prod_yield(i))) // "*"
       end if
       ret_val = trim(ret_val) // " " &
            // trim(yield) // trim(aq_spec_data%name(val%product(i))) &
            // "(" // trim(integer_to_string(val%product(i))) // ")"
    end do

    ret_val = trim(ret_val) // " ; rate " &
         // trim(aq_rxn_data_get_rate_name(val%rate_index)) // ":"

    do i = 1,size(val%rate_param)
       if (i > 1) then
          ret_val = trim(ret_val) // ","
       end if
       ret_val = trim(ret_val) // " " &
            // trim(real_or_integer_to_string(val%rate_param(i)))
    end do
    
    aq_rxn_to_string = adjustl(ret_val)

  end function aq_rxn_to_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aq_rxn_data(val)

    !> Value to pack.
    type(aq_rxn_data_t), intent(in) :: val

    pmc_mpi_pack_size_aq_rxn_data = &
         pmc_mpi_pack_size_integer(val%class_index) &
         + pmc_mpi_pack_size_integer_array(val%reactant) &
         + pmc_mpi_pack_size_integer_array(val%product) &
         + pmc_mpi_pack_size_real_array(val%prod_yield) &
         + pmc_mpi_pack_size_integer(val%rate_index) &
         + pmc_mpi_pack_size_real_array(val%rate_param)

  end function pmc_mpi_pack_size_aq_rxn_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aq_rxn_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_rxn_data_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%class_index)
    call pmc_mpi_pack_integer_array(buffer, position, val%reactant)
    call pmc_mpi_pack_integer_array(buffer, position, val%product)
    call pmc_mpi_pack_real_array(buffer, position, val%prod_yield)
    call pmc_mpi_pack_integer(buffer, position, val%rate_index)
    call pmc_mpi_pack_real_array(buffer, position, val%rate_param)
    call assert(602214184, &
         position - prev_position <= pmc_mpi_pack_size_aq_rxn_data(val))
#endif

  end subroutine pmc_mpi_pack_aq_rxn_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aq_rxn_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_rxn_data_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%class_index)
    call pmc_mpi_unpack_integer_array(buffer, position, val%reactant)
    call pmc_mpi_unpack_integer_array(buffer, position, val%product)
    call pmc_mpi_unpack_real_array(buffer, position, val%prod_yield)
    call pmc_mpi_unpack_integer(buffer, position, val%rate_index)
    call pmc_mpi_unpack_real_array(buffer, position, val%rate_param)
    call assert(602348640, &
         position - prev_position <= pmc_mpi_pack_size_aq_rxn_data(val))
#endif

  end subroutine pmc_mpi_unpack_aq_rxn_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aq_rxn_data





