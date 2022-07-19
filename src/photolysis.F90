! Copyright (C) 2020, 2021 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> PartMC interface to a photolysis module

!> The photolysis_t type and related functions
module pmc_photolysis

  use pmc_constants
#ifdef PMC_USE_CAMP
  use camp_camp_core
  use camp_rxn_data
  use camp_rxn_photolysis
#endif
  use pmc_util

#ifdef PMC_USE_CAMP

  implicit none
  private

  public :: photolysis_t

  !> Provider of photolysis rates
  type :: photolysis_t
  private
    !> Pointer to the CAMP core to update rates for
    type(camp_core_t), pointer :: camp_core => null()
    !> Base photolysis rates (s-1)
    real(kind=dp), allocatable :: base_rates(:)
    !> Rate update objects
    type(rxn_update_data_photolysis_t), allocatable :: photo_rxns(:)
  contains
    !> Update the CAMP photolysis rate constants
    procedure :: update_rate_constants
    !> Determine the number of bytes required to pack the given value
    procedure :: pack_size
    !> Packs the given value into the buffer, advancing position
    procedure :: bin_pack
    !> Unpacks the given value from the buffer, advancing position
    procedure :: bin_unpack
    !> Print the photolysis module data
    procedure :: print => do_print
  end type photolysis_t

  !> Constructor for photolysis_t
  interface photolysis_t
    procedure :: constructor
  end interface photolysis_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor for photolysis_t
  function constructor(camp_core) result(new_obj)

    !> A new photolysis rate calculator
    type(photolysis_t), pointer :: new_obj
    !> The CAMP core containing photolysis reactions
    type(camp_core_t), pointer, intent(in), optional :: camp_core

    character(len=:), allocatable :: rxn_key, rxn_val, rate_key, str_val
    real(kind=dp) :: rate_val
    integer :: i_mech, i_rxn, i_photo_rxn, n_photo_rxns
    class(rxn_data_t), pointer :: rxn

    ! If no core is provided, return an empty object that can be
    ! unpacked from an MPI buffer
    new_obj => NULL()
    if (.not. present(camp_core)) return

    ! Strings needed to find all the photolysis reactions
    rxn_key = "type"
    rxn_val = "PHOTOLYSIS"
    rate_key = "base rate"

    call assert(254347663, associated(camp_core))
    call assert(689045432, camp_core%is_initialized())
    call assert(256253197, associated(camp_core%mechanism))

    allocate(new_obj)
    new_obj%camp_core => camp_core

    ! Count the photolysis reactions
    n_photo_rxns = 0
    do i_mech = 1, size(camp_core%mechanism)
      do i_rxn = 1, camp_core%mechanism(i_mech)%val%size()
        rxn => camp_core%mechanism(i_mech)%val%get_rxn(i_rxn)
        call assert(106297725, rxn%property_set%get_string(rxn_key, str_val))
        if (trim(str_val) == rxn_val) n_photo_rxns = n_photo_rxns + 1
      end do
    end do

    ! Create update rate objects for photolysis reactions in all CAMP
    ! mechanisms
    i_photo_rxn  = 0
    allocate(new_obj%photo_rxns(n_photo_rxns))
    allocate(new_obj%base_rates(n_photo_rxns))
    do i_mech = 1, size(camp_core%mechanism)

      ! Loop through all reactions in mechanism to find photolysis reactions
      do i_rxn = 1, camp_core%mechanism(i_mech)%val%size()
        rxn => camp_core%mechanism(i_mech)%val%get_rxn(i_rxn)
        call assert(799145523, rxn%property_set%get_string(rxn_key, str_val))

        ! Is this a photolysis reaction?
        if (trim(str_val) /= rxn_val) cycle
        i_photo_rxn = i_photo_rxn + 1

        ! Get the base photolysis rate
        call assert_msg(501329648, &
             rxn%property_set%get_real(rate_key, rate_val), &
             "Missing 'base rate' for photolysis reaction "// &
             trim(integer_to_string(i_photo_rxn)))
        new_obj%base_rates(i_photo_rxn) = rate_val

        ! Create an update rate object for this photolysis reaction
        select type (rxn_photo => rxn)
          class is (rxn_photolysis_t)
          call camp_core%initialize_update_object(rxn_photo, &
               new_obj%photo_rxns(i_photo_rxn))
          class default
            call die(722633162)
        end select
      end do
    end do

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the photolysis rates
  subroutine update_rate_constants(this)

    !> Photolysis calculator
    class(photolysis_t), intent(inout) :: this

    integer :: i_rxn

    ! NOTE: Update this after connection to a real photolysis module
    do i_rxn = 1, size(this%photo_rxns)
      call this%photo_rxns(i_rxn)%set_rate(this%base_rates(i_rxn))
      call this%camp_core%update_data(this%photo_rxns(i_rxn))
    end do

  end subroutine update_rate_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine the size of a binary required to pack the photolysis data
  integer function pack_size(this)

    use pmc_mpi

    !> Photolysis calculator
    class(photolysis_t), intent(in) :: this

    integer :: i_rxn

#ifdef PMC_USE_MPI

    call assert(127027009, allocated(this%base_rates))
    call assert(634138948, allocated(this%photo_rxns))

    pack_size = &
      pmc_mpi_pack_size_real_array(this%base_rates) + &
      pmc_mpi_pack_size_integer(size(this%photo_rxns))

    do i_rxn = 1, size(this%photo_rxns)
      pack_size = pack_size + this%photo_rxns(i_rxn)%pack_size()
    end do
#else
  pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_pack(this, buffer, pos)

    use pmc_mpi

    !> Photolysis calculator
    class(photolysis_t), intent(in) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: i_rxn, prev_position

    call assert(971093983, allocated(this%base_rates))
    call assert(913255424, allocated(this%photo_rxns))

    prev_position = pos
    call pmc_mpi_pack_real_array(buffer, pos, this%base_rates)
    call pmc_mpi_pack_integer(buffer, pos, size(this%photo_rxns))
    do i_rxn = 1, size(this%photo_rxns)
      call this%photo_rxns(i_rxn)%bin_pack(buffer, pos)
    end do
    call assert(234533342, pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Pack the given value to the buffer, advancing position
  subroutine bin_unpack(this, buffer, pos)

    use pmc_mpi

    !> Photolysis calculator
    class(photolysis_t), intent(out) :: this
    !> Memory buffer
    character, intent(inout) :: buffer(:)
    !> Current buffer position
    integer, intent(inout) :: pos

#ifdef PMC_USE_MPI
    integer :: i_rxn, n_rxns, prev_position

    prev_position = pos
    call pmc_mpi_unpack_real_array(buffer, pos, this%base_rates)
    call pmc_mpi_unpack_integer(buffer, pos, n_rxns)
    allocate(this%photo_rxns(n_rxns))
    do i_rxn = 1, size(this%photo_rxns)
      call this%photo_rxns(i_rxn)%bin_unpack(buffer, pos)
    end do
    call assert(391255154, pos - prev_position <= this%pack_size())
#endif

  end subroutine bin_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the photolysis data
  subroutine do_print(this, file_unit)

    !> Photolysis calculator
    class(photolysis_t), intent(in) :: this
    !> File unit for output
    integer, optional :: file_unit

    integer :: f_unit, i_rxn

    f_unit = 6
    if (present(file_unit)) f_unit = file_unit

    write(f_unit,*) "***************************"
    write(f_unit,*) "***   Photolysis Data   ***"
    write(f_unit,*) "***************************"
    write(f_unit,*) ""
    if (allocated(this%base_rates)) then
      do i_rxn = 1, size(this%base_rates)
        write(f_unit,*) " photo rxn(",i_rxn,") = ", this%base_rates(i_rxn)
      end do
    else
      write(f_unit,*) "  No photolysis data"
    end if
    write(f_unit,*) ""
    write(f_unit,*) "***************************"
    write(f_unit,*) "*** End Photolysis Data ***"
    write(f_unit,*) "***************************"

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
end module pmc_photolysis
