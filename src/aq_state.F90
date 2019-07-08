! Copyright (C) 2015 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_state module.

!> The aq_state_t structure and assocated subroutines.
module pmc_aq_state

  use pmc_aq_spec_data
  use pmc_spec_file
#ifdef PMC_USE_MPI
  use mpi
#endif

  implicit none

  character(len=AQ_SPEC_NAME_LEN), parameter :: AQ_STATE_OH_SPEC_NAME = "OHm"

  !> The current state of aerosol properties and aqueous-phase
  !! chemistry related species
  !!
  !! The aq. chemistry related species are defined by the aq_spec_data_t
  !! structure, so that \c aq_state%%mix_rat(i) is the current mixing 
  !! ratio of the species with name \c aq_spec_data%%name(i), etc.
  type aq_state_t
     !> Length n_spec, mixing ratio (ppb).
     real(kind=dp), pointer :: mix_rat(:)
    !> Particle radius (m)
    real(kind=dp) :: radius
    !> Particle number concentration (#/cc)
    real(kind=dp) :: n_particle
  end type aq_state_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for aq. chemistry related species.
  subroutine aq_state_allocate(aq_state)

    !> Aq. chemistry state to be allocated.
    type(aq_state_t), intent(out) :: aq_state

    allocate(aq_state%mix_rat(0))

  end subroutine aq_state_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for aq. chemistry related species of the given size.
  subroutine aq_state_allocate_size(aq_state, n_spec)

    !> Aq. chemistry state to be allocated.
    type(aq_state_t), intent(out) :: aq_state
    !> Number of species.
    integer, intent(in) :: n_spec

    allocate(aq_state%mix_rat(n_spec))
    call aq_state_zero(aq_state)

  end subroutine aq_state_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aq_state_deallocate(aq_state)

    !> Aq. chemistry state to be freed.
    type(aq_state_t), intent(inout) :: aq_state

    deallocate(aq_state%mix_rat)

  end subroutine aq_state_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Zeros the state.
  subroutine aq_state_zero(aq_state)

    !> Aq. chemistry state.
    type(aq_state_t), intent(inout) :: aq_state

    aq_state%mix_rat(:) = 0d0
    aq_state%radius = 0.0
    aq_state%n_particle = 0.0

  end subroutine aq_state_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy to an already allocated to_state.
  subroutine aq_state_copy(from_state, to_state)

    !> Existing aq. chemistry state.
    type(aq_state_t), intent(in) :: from_state
    !> Must be allocated already.
    type(aq_state_t), intent(inout) :: to_state

    integer :: n_spec

    n_spec = size(from_state%mix_rat)
    deallocate(to_state%mix_rat)
    allocate(to_state%mix_rat(n_spec))
    to_state%mix_rat(:) = from_state%mix_rat(:)
    to_state%radius = from_state%radius
    to_state%n_particle = from_state%n_particle

  end subroutine aq_state_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aq. chemistry state from the file named on the line read from file.
  subroutine aq_spec_file_read_aq_state(file, aq_spec_data, aq_state)

    !> File to read aq. chemistry state from.
    type(spec_file_t), intent(inout) :: file
    !> Aq. chemistry related species data.
    type(aq_spec_data_t), intent(in) :: aq_spec_data
    !> Aq. chemistry state data to read.
    type(aq_state_t), intent(inout) :: aq_state

    integer :: n_species, species, i
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_name(:)
    real(kind=dp), pointer :: species_data(:,:)

    !> \page input_format_aq_state Input File Format: Aqueous Chemistry Initial State
    !!
    !! This file specifies the intial concentrations of certain species
    !! in the aqueous chemistry mechanism. IMPORTANT: These species
    !! concentrations will be reset to these values at EACH call to
    !! the \ref aq_chem_timestep subroutine. This is intended for species
    !! that are included in the aqueous chemistry mechanism, but not
    !! in the PartMC gas or aerosol set of species, and whose 
    !! concentrations can be considered constant. Any species included
    !! here should not be included in the aqueous chemistry species map file.
    !!
    !! An aqueous chemistry state input file must consist of one line per
    !! species, with each line having the species name followed by the
    !! species mixing ratio in atm (for gas-phase species) or M (for
    !! aqueous-phase species). The valid species names are those 
    !! specfied in the CAPRAM mechanism input file, but not all 
    !! species have to be listed. Any missing species will have mixing 
    !! ratios of zero. For example, an aqueous chemistry state file could
    !! contain:
    !! <pre>
    !! # species    initial concentration (atm or M)
    !! CO2          5.0E-04
    !! [O2]         2.1E-01
    !! [aH2O]       5.5E+01    
    !! </pre>
    !!
    !! The intitial state file is specified by the parameter:
    !!   - \b aq_init (string): name of file from which to read the
    !!   initial aqueous chemistry state
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aq_spec_map --- the Aqueous chemistry
    !!       to PartMC map file

    allocate(species_name(0))
    allocate(species_data(0,0))
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    n_species = size(species_data, 1)
    if (.not. ((size(species_data, 2) == 1) .or. (n_species == 0))) then
       call die_msg(608970598, 'each line in ' // trim(file%name) &
            // ' must contain exactly one data value')
    end if

    ! copy over the data
    call aq_state_deallocate(aq_state)
    call aq_state_allocate_size(aq_state, aq_spec_data%n_spec)
    aq_state%mix_rat = 0d0
    do i = 1,n_species
       species = aq_spec_data_spec_by_name(aq_spec_data, species_name(i))
       if (species == 0) then
          call die_msg(609088247, 'unknown species ' // &
               trim(species_name(i)) // ' in file ' // trim(file%name))
       end if
       aq_state%mix_rat(species) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)

  end subroutine aq_spec_file_read_aq_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print aq. chemistry state
  subroutine pmc_aq_state_print(aq_state, aq_spec_data)

    !> Aqueous chemistry state
    type(aq_state_t) :: aq_state
    !> Aqueous chemistry species data
    type(aq_spec_data_t) :: aq_spec_data

    integer :: i

    write(*,*) " "
    write(*,*) "Aqeuous Chemistry Related Species Concentrations"
    write(*,*) " "

    do i=1, aq_spec_data%n_spec
        write(*,'(A20,G10.2)') trim(aq_spec_data%name(i)), aq_state%mix_rat(i)
    enddo

  end subroutine pmc_aq_state_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aq_state(val)

    !> Value to pack.
    type(aq_state_t), intent(in) :: val

    pmc_mpi_pack_size_aq_state = &
         pmc_mpi_pack_size_real_array(val%mix_rat) &
         + pmc_mpi_pack_size_real(val%radius) &
         + pmc_mpi_pack_size_real(val%n_particle)

  end function pmc_mpi_pack_size_aq_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aq_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_state_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_real_array(buffer, position, val%mix_rat)
    call pmc_mpi_pack_real(buffer, position, val%radius)
    call pmc_mpi_pack_real(buffer, position, val%n_particle)
    call assert(609441194, &
         position - prev_position <= pmc_mpi_pack_size_aq_state(val))
#endif

  end subroutine pmc_mpi_pack_aq_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aq_state(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_state_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_real_array(buffer, position, val%mix_rat)
    call pmc_mpi_unpack_real(buffer, position, val%radius)
    call pmc_mpi_unpack_real(buffer, position, val%n_particle)
    call assert(609592457, &
         position - prev_position <= pmc_mpi_pack_size_aq_state(val))
#endif

  end subroutine pmc_mpi_unpack_aq_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Charge balance the system using [OH-]aq
  subroutine aq_state_charge_balance(aq_spec_data, aq_state)

    !> Aq. chemistry related species data.
    type(aq_spec_data_t), intent(in) :: aq_spec_data
    !> Aq. chemistry state data to read.
    type(aq_state_t), intent(inout) :: aq_state

    real(kind=dp) :: n_charge ! (mol/L)
    integer :: i_OH, i

    i_OH = aq_spec_data_spec_by_name(aq_spec_data, trim(AQ_STATE_OH_SPEC_NAME))

    if (i_OH.eq.0) then
       call die_msg(610785754, 'Cannot find species OH- in aqueous chemistry mechanism')
    endif

    n_charge = 0.0
    do i=1, aq_spec_data%n_spec
        n_charge = aq_state%mix_rat(i) * aq_spec_data%charge(i)
    enddo

    if (n_charge.ge.0.0) then
        aq_state%mix_rat(i_OH) = n_charge
    else
        aq_state%mix_rat(i_OH) = 0.0
    endif

  end subroutine aq_state_charge_balance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aq_state




