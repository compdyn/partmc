! Copyright (C) 2005-2015 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_spec_data module.

!> The aq_spec_data_t structure and associated subroutines.
module pmc_aq_spec_data

  use pmc_spec_file
  use pmc_constants
  use pmc_gas_data
  use pmc_aero_data
#ifdef PMC_USE_MPI
  use mpi
#endif

  implicit none

  !> Maximum length of the name of a species.
  integer, parameter :: AQ_SPEC_NAME_LEN = 100
  !> Default absolute tolerance for integration
  !> Units match those of aq_state%mix_rat
  real(kind=dp), parameter :: AQ_SPEC_DEFAULT_ABSTOL = 1.0E-30

  !> Constant aqueous chemistry species data.
  !!
  !! This includes gas and particle phase species involved in aqueous 
  !! chemistry.
  !!
  !! Each gas species is identified by an integer \c i between 1 and
  !! \c n_spec. Species \c i has name \c aq_spec_data%%name(i). The
  !! variable aqueous species data describing the current mixing ratios
  !! are stored in the aq_spec_state_t structure, so the mixing ratio
  !! of species \c i is aq_spec_state%%mix_rat(i).
  type aq_spec_data_t
     !> Number of species.
     integer :: n_spec
     !> Species name [length \c n_spec].
     character(len=AQ_SPEC_NAME_LEN), pointer :: name(:)
     !> Index of the corresponding gas-phase PartMC species [length
     !> \c n_spec]. \c pmc_gas_index (i) is the index of species \c i in
     !> the gas_data_t variable, or 0 if there is no match.
     integer, pointer :: pmc_gas_index(:)
     !> Index of the corresponding aerosol-phase PartMC species [length
     !> \c n_spec]. \c pmc_aero_index (i) is the index of species \c i 
     !> in the aero_data_t variable, or 0 if there is no match.
     integer, pointer :: pmc_aero_index(:)
     !> Absolute tolerance for integration (units matching those of 
     !> aq_state%mix_rat
     real(kind=dp), pointer :: abstol(:)
     !> Constant concentration species?
     logical, pointer :: const_conc(:)
     !> Gas-phase diffusion coefficient (m^2/s)
     real(kind=dp), pointer :: Dg(:)
     !> N_star (see Ervens et al. JGR 2003)
     real(kind=dp), pointer :: N_star(:)
     !> Molecular weight (kg/mol)
     real(kind=dp), pointer :: MW(:)
     !> Net charge (+e)
     integer, pointer :: charge(:)
  end type aq_spec_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for aqueous chemistry related species.
  subroutine aq_spec_data_allocate(aq_spec_data)

    !> Aqueous chemistry related species data.
    type(aq_spec_data_t), intent(out) :: aq_spec_data

    aq_spec_data%n_spec = 0
    allocate(aq_spec_data%name(0))
    allocate(aq_spec_data%pmc_gas_index(0))
    allocate(aq_spec_data%pmc_aero_index(0))
    allocate(aq_spec_data%abstol(0))
    allocate(aq_spec_data%const_conc(0))
    allocate(aq_spec_data%Dg(0))
    allocate(aq_spec_data%N_star(0))
    allocate(aq_spec_data%MW(0))
    allocate(aq_spec_data%charge(0))

  end subroutine aq_spec_data_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for aqueous chemistry related species with the given size.
  subroutine aq_spec_data_allocate_size(aq_spec_data, n_spec)

    !> Aqueous chemistry related data.
    type(aq_spec_data_t), intent(out) :: aq_spec_data
    !> Number of species.
    integer, intent(in) :: n_spec

    aq_spec_data%n_spec = n_spec
    allocate(aq_spec_data%name(n_spec))
    allocate(aq_spec_data%pmc_gas_index(n_spec))
    allocate(aq_spec_data%pmc_aero_index(n_spec))
    allocate(aq_spec_data%abstol(n_spec))
    allocate(aq_spec_data%const_conc(n_spec))
    allocate(aq_spec_data%Dg(n_spec))
    allocate(aq_spec_data%N_star(n_spec))
    allocate(aq_spec_data%MW(n_spec))
    allocate(aq_spec_data%charge(n_spec))

    aq_spec_data%pmc_gas_index(n_spec) = 0
    aq_spec_data%pmc_aero_index(n_spec) = 0
    aq_spec_data%abstol(:) = AQ_SPEC_DEFAULT_ABSTOL
    aq_spec_data%const_conc(:) = .false.
    aq_spec_data%Dg(:) = 0.0
    aq_spec_data%N_star(:) = 0.0
    aq_spec_data%MW(:) = 0.0
    aq_spec_data%charge(:) = 0

  end subroutine aq_spec_data_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aq_spec_data_deallocate(aq_spec_data)

    !> Aqueous chemistry related data.
    type(aq_spec_data_t), intent(inout) :: aq_spec_data

    deallocate(aq_spec_data%name)
    deallocate(aq_spec_data%pmc_gas_index)
    deallocate(aq_spec_data%pmc_aero_index)
    deallocate(aq_spec_data%abstol)
    deallocate(aq_spec_data%const_conc)
    deallocate(aq_spec_data%Dg)
    deallocate(aq_spec_data%N_star)
    deallocate(aq_spec_data%MW)
    deallocate(aq_spec_data%charge)

  end subroutine aq_spec_data_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy data from one aq_spec_data_t variable to another
  subroutine aq_spec_data_copy(orig_data, copy_of_orig, n_spec)

    !> Original species data
    type(aq_spec_data_t), intent(in) :: orig_data
    !> Copy of original species data
    type(aq_spec_data_t), intent(inout) :: copy_of_orig
    !> Number of records to copy
    integer, intent(in) :: n_spec

    copy_of_orig%name(1:n_spec)           = orig_data%name(1:n_spec)
    copy_of_orig%pmc_gas_index(1:n_spec)  = orig_data%pmc_gas_index(1:n_spec)
    copy_of_orig%pmc_aero_index(1:n_spec) = orig_data%pmc_aero_index(1:n_spec)
    copy_of_orig%abstol(1:n_spec)         = orig_data%abstol(1:n_spec)
    copy_of_orig%const_conc(1:n_spec)     = orig_data%const_conc(1:n_spec)
    copy_of_orig%Dg(1:n_spec)             = orig_data%Dg(1:n_spec)
    copy_of_orig%N_star(1:n_spec)         = orig_data%N_star(1:n_spec)
    copy_of_orig%MW(1:n_spec)             = orig_data%MW(1:n_spec)
    copy_of_orig%charge(1:n_spec)         = orig_data%charge(1:n_spec)

  end subroutine aq_spec_data_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy data from one aq_spec_data_t variable to another in a different order
  subroutine aq_spec_data_copy_reorder(orig_data, copy_of_orig, spec_map, n_spec)

    !> Original species data
    type(aq_spec_data_t), intent(in) :: orig_data
    !> Copy of original species data
    type(aq_spec_data_t), intent(inout) :: copy_of_orig
    !> Map of species indices from orig_data to copy_of_orig
    !! [copy_of_orig%prop(spec_map(i)) = orig_data%prop(i)]
    integer, intent(in) :: spec_map(:)
    !> Number of records to copy
    integer, intent(in) :: n_spec

    integer :: i

    do i=1,n_spec
        copy_of_orig%name(spec_map(i))           = orig_data%name(i)
        copy_of_orig%pmc_gas_index(spec_map(i))  = orig_data%pmc_gas_index(i)
        copy_of_orig%pmc_aero_index(spec_map(i)) = orig_data%pmc_aero_index(i)
        copy_of_orig%abstol(spec_map(i))         = orig_data%abstol(i)
        copy_of_orig%const_conc(spec_map(i))     = orig_data%const_conc(i)
        copy_of_orig%Dg(spec_map(i))             = orig_data%Dg(i)
        copy_of_orig%N_star(spec_map(i))         = orig_data%N_star(i)
        copy_of_orig%MW(spec_map(i))             = orig_data%MW(i)
        copy_of_orig%charge(spec_map(i))         = orig_data%charge(i)
    enddo
  end subroutine aq_spec_data_copy_reorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Increase the size of an aq_spec_data_t variable by 1 and add a new species
  integer function aq_spec_data_add_spec(aq_spec_data, name)

    ! FIXME: It seems like there should be a way to allocate a
    ! new aq_spec_data_t variable, set the new values and then set
    ! aq_spec_data to point to the new variable. 

    !> Aqueous chemistry related data.
    type(aq_spec_data_t), intent(inout) :: aq_spec_data
    !> Name of new species
    character(len=*), intent(in) :: name

    ! Temporary species data
    type(aq_spec_data_t) :: aq_spec_data_temp
    ! Number of species
    integer :: n_spec

    n_spec = aq_spec_data%n_spec

    call aq_spec_data_allocate_size(aq_spec_data_temp, n_spec+1)

    call aq_spec_data_copy(aq_spec_data, aq_spec_data_temp, n_spec)

    call aq_spec_data_deallocate(aq_spec_data)

    n_spec = n_spec + 1

    call aq_spec_data_allocate_size(aq_spec_data, n_spec)

    call aq_spec_data_copy(aq_spec_data_temp, aq_spec_data, n_spec)

    call aq_spec_data_deallocate(aq_spec_data_temp)

    aq_spec_data%name(n_spec) = name

    aq_spec_data_add_spec = aq_spec_data%n_spec

  end function aq_spec_data_add_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of the species with the given name, or
  !> returns 0 if there is no such species.
  integer function aq_spec_data_spec_by_name(aq_spec_data, name)

    !> Aqueous chemistry related species data.
    type(aq_spec_data_t), intent(in) :: aq_spec_data
    !> Name of species to find.
    character(len=*), intent(in) :: name

    integer i
    logical found

    found = .false.
    do i = 1,aq_spec_data%n_spec
       if (name == aq_spec_data%name(i)) then
          found = .true.
          exit
       end if
    end do
    if (found) then
       aq_spec_data_spec_by_name = i
    else
       aq_spec_data_spec_by_name = 0
    end if

  end function aq_spec_data_spec_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the mass accomodation factor based on N* and Temperature (K)
  real(kind=dp) function aq_rxn_data_get_mass_accom(N_star, temp)

    !> N* (see Ervens et al. for description, unitless)
    real(kind=dp), intent(in) :: N_star
    !> temperature (K)
    real(kind=dp), intent(in) :: temp

    real(kind=dp), parameter :: J_per_kcal = 4184.0 ! (J/kcal)
    real(kind=dp) :: del_H, del_S, del_G, beta

    !> Based on equations in:
    !! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended
    !! and condensed tropospheric aqueous phase mechanism and its
    !! application."" J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202

    ! enthalpy change (kcal mol-1)
    del_H = - 10.0*(N_star-1.0) + 7.53*(N_star**(2.0/3.0)-1.0) - 1.0

    ! entropy change (cal mol-1)
    del_S = - 32.0*(N_star-1.0) + 9.21*(N_star**(2.0/3.0)-1.0) - 1.3

    ! Gibb's free energy change (J mol-1)
    del_G = (del_H - temp*del_S/1000.0) * J_per_kcal

    ! Calculate mass accomodation coefficient
    beta = exp(-del_G/(const%univ_gas_const*temp))
    aq_rxn_data_get_mass_accom = beta / (1.0 + beta)

  end function aq_rxn_data_get_mass_accom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the mass transfer rate constant for a gas-phase species
  real(kind=dp) function aq_rxn_data_get_mass_transfer_rc(Dg, radius, alpha, &
                            MW, temp)

    !> gas-phase diffusion coefficient (m^2/s)
    real(kind=dp), intent(in) :: Dg
    !> particle radius (m)
    real(kind=dp), intent(in) :: radius
    !> mass accomodation coefficient (unitless)
    real(kind=dp), intent(in) :: alpha
    !> molecular weight (kg/mol)
    real(kind=dp), intent(in) :: MW
    !> temperature (K)
    real(kind=dp), intent(in) :: temp

    ! Gas constant (L*atm/K/mol)
    real(kind=dp), parameter :: R_gas = 0.08205736

    ! mean velocity of a gas-phase species (m/s)
    real(kind=dp) :: c_rms

    !> Based on equations in:
    !! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended
    !! and condensed tropospheric aqueous phase mechanism and its
    !! application."" J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202

    ! calculate average velocity (m/s)
    c_rms = sqrt(8.0*const%univ_gas_const*temp/(const%pi*MW))

    ! calculate rate constant for mass transfer (diffusion limited) (1/s)
    aq_rxn_data_get_mass_transfer_rc = (radius**2/(3.0*Dg) + &
        4.0*radius/(3.0*c_rms*alpha))**(-1)

    ! convert to units of (M/atm/s)
    aq_rxn_data_get_mass_transfer_rc = aq_rxn_data_get_mass_transfer_rc &
        / (R_gas * temp)

  end function aq_rxn_data_get_mass_transfer_rc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read species data from a .spec file.
  subroutine aq_spec_file_read_data(file, aq_spec_data)

    !> Spec file to read data from.
    type(spec_file_t), intent(inout) :: file
    !> Aqueous species data.
    type(aq_spec_data_t), intent(inout) :: aq_spec_data

    integer :: species_index, i
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_name(:)
    real(kind=dp), pointer :: species_data(:,:)

    !> \page input_format_aq_spec_data Input File Format: Aqueous Chemistry Species Data
    !!
    !! An aqueous chemistry species data file must consist of one line per
    !! species, with each line having:
    !!   - species name (string)
    !!   - gas-phase diffusion coefficient (real, unit m^2/s)
    !!   - molecular weight (real, unit kg/mol)
    !!   - N* parameter (unitless). 
    !!   - absolute integration tolerance (atm or M)
    !!   - number of positive charges
    !!   - number of negative charges
    !!
    !! It is not necessary to include every species in the aqueous
    !! chemistry mechanism in this file, only the ones that partition
    !! between phases or use a non-standard integration tolerance.
    !! For partitioning species that use the standard integration 
    !! tolerance, the absolute tolerance can be set to zero.
    !! For non-partitioning species that need a non-standard absolute
    !! integration tolerance, the other parameters can be set to zero.
    !! For example, an aq_spec_data file could contain:
    !! <pre>
    !! # species  diff coeff (m^2/s) MW (kg/mol) N* (unitless) Abs Tol (atm or M)  +Chg   -Chg
    !! CO2        1.55e5             4.401e-2    4.493         0.0                  0       0
    !! HCl        1.89e5             3.646e-2    1.841         1e-7                 0       0
    !! NH3        1.3e5              1.7031e-2   1.924         1e-7                 0       0
    !! FEpp       0.0                0.0         0.0           1e-10                2       0
    !! </pre>
    !!
    !! The aqueous chemistry species data file is specified by the parameter:
    !!   - \b aq_spec (string): name of file from which to read the
    !!   aqueous chemistry related species
    !!
    !! The reference describing the N* parameter is:
    !! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended 
    !! and condensed tropospheric aqueous phase mechanism and its 
    !! application."" J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    ! read the gas data from the specified file
    allocate(species_name(0))
    allocate(species_data(0,0))
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    if (size(species_data, 2) /= 6) then
       call die_msg(607508389, 'each line in ' // trim(file%name) &
            // ' must contain the species name, gas-phase diffusion coefficient,' &
            // ' molecular weight, N* parameter and absolute integration tolerance.')
    end if

    ! copy over the data
    do i = 1,size(species_data,1)
        species_index = aq_spec_data_spec_by_name(aq_spec_data, species_name(i))
        if (species_index.eq.0) then
            call die_msg(607642845, 'species ' // trim(species_name(i)) // ' in ' &
              // trim(file%name) // ' is not a species in the aq. chemistry mechanism.')
        endif
        aq_spec_data%Dg(species_index) = species_data(i,1)
        aq_spec_data%MW(species_index) = species_data(i,2)
        aq_spec_data%N_star(species_index) = species_data(i,3)
        if (species_data(i,4).gt.0.0) then
            aq_spec_data%abstol(species_index) = species_data(i,4)
        endif
        aq_spec_data%charge(species_index) = int(species_data(i,5) - species_data(i,6))

    end do
    deallocate(species_name)
    deallocate(species_data)

  end subroutine aq_spec_file_read_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read species map from a .spec file.
  subroutine aq_spec_file_read_map(file, aq_spec_data, gas_data, aero_data)

    !> Spec file to read data from.
    type(spec_file_t), intent(inout) :: file
    !> Aqueous species data.
    type(aq_spec_data_t), intent(inout) :: aq_spec_data
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data

    integer :: species_index, i
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_name(:)
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_data(:,:)

    !> \page input_format_aq_spec_map Input File Format: Aqueous Chemistry Species Map
    !!
    !! An aqueous chemistry species map file must consist of one line per
    !! species, with each line having
    !!   - species name in aqueous chemistry mechanism (string)
    !!   - phase (string, "gas" or "aero")
    !!   - species name in PartMC (string)
    !!
    !! It is not necessary to include every species in the aqueous
    !! chemistry mechanism in this file, only the ones whose concentration
    !! should be preserved between calls to the aqueous chemistry code.
    !! Species included in this file should not be included in the
    !! the aqueous chemistry initial state file. For example, an
    !! aq_spec_data file could contain:
    !! <pre>
    !! # CAPRAM name   Phase   PartMC name
    !! SO2             gas     SO2
    !! OP1             gas     CH3OOH
    !! SO4             aero    SO4
    !! </pre>
    !!
    !! The aqueous chemistry species map file is specified by the parameter:
    !!   - \b aq_spec (string): name of file from which to read the
    !!   aqueous chemistry related species map
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format
    !!   - \ref input_format_aq_state --- the aqueous chemistry initial
    !!           state file

    ! read the mapping data from the specified file
    allocate(species_name(0))
    allocate(species_data(0,0))
    call spec_file_read_string_named_array(file, 0, species_name, &
         species_data)

    ! check the data size
    if (size(species_data, 2) /= 2) then
       call die_msg(607861336, 'each line in ' // trim(file%name) &
            // ' must contain the species name in the aqueous chemistry mechanism' &
            // ' followed by ''gas'' or ''aero'' and then the species name in PartMC.')
    end if

    ! copy over the data
    do i = 1,size(species_data,1)
        species_index = aq_spec_data_spec_by_name(aq_spec_data, species_name(i))
        if (species_index.eq.0) then
            call die_msg(607995792, 'species ' // trim(species_name(i)) // ' in ' &
              // trim(file%name) // ' is not a species in the aq. chemistry mechanism.')
        endif
        select case (species_data(i,1))
            case ("gas")
                aq_spec_data%pmc_gas_index(species_index) = &
                    gas_data_spec_by_name(gas_data, trim(species_data(i,2)))
            case ("aero")
                aq_spec_data%pmc_aero_index(species_index) = &
                    aero_data_spec_by_name(aero_data, trim(species_data(i,2)))
            case default
                call die_msg(608130248, 'species ' // trim(species_name(i)) // ' in ' &
                    // trim(file%name) // ' needs to be either ''gas'' or ''aero''')
        end select
        if (aq_spec_data%pmc_gas_index(species_index).eq.0 &
            .and. aq_spec_data%pmc_aero_index(species_index).eq.0) then
            call die_msg(608247897, 'species ' // trim(species_name(i)) // ' in ' &
              // trim(file%name) // ' is not a species in PartMC.')
        endif
    end do
    deallocate(species_name)
    deallocate(species_data)

  end subroutine aq_spec_file_read_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out a species list
  subroutine aq_spec_data_print(aq_spec_data)

    !> Aqueous chemistry related species data.
    type(aq_spec_data_t), intent(in) :: aq_spec_data

    integer :: i

    real(kind=dp) :: alpha_288

    write(*,*) "Aqueous-Phase Chemistry Related Species"
    write(*,*) " "
    write(*,*) "Species Name   Constant?  Abs Tol (M or atm)  Dg (m2 s-1)   MW (kg mol-1)" &
                // " N_star (unitless)  Charge  alpha(288K)"

    do i=1,aq_spec_data%n_spec
        if (aq_spec_data%N_star(i).ne.0.0) then
            alpha_288 = aq_rxn_data_get_mass_accom(aq_spec_data%N_star(i), real(288.0,dp))
            write(*,'(A25,g10.5,A,g10.5,A,g10.5,A,g10.5,A,I3,A,g10.5)') trim(aq_spec_data%name(i)) // &
                "   " // trim(logical_to_string(aq_spec_data%const_conc(i))) // "     ", &
                aq_spec_data%abstol(i),"        ", &
                aq_spec_data%Dg(i), "      ", aq_spec_data%MW(i), "        ", &
                aq_spec_data%N_star(i), "      ",aq_spec_data%charge(i), "    ", alpha_288
        else
            write(*,'(A25,g10.5,A,g10.5,A,g10.5,A,g10.5,A,I3)') trim(aq_spec_data%name(i)) // &
                "   " // trim(logical_to_string(aq_spec_data%const_conc(i))) // "     ", &
                aq_spec_data%abstol(i),"        ", &
                aq_spec_data%Dg(i), "      ", aq_spec_data%MW(i), "        ", &
                aq_spec_data%N_star(i), "      ",aq_spec_data%charge(i)
        endif
    enddo

    write(*,*) " "

  end subroutine aq_spec_data_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aq_spec_data(val)

    !> Value to pack.
    type(aq_spec_data_t), intent(in) :: val

    pmc_mpi_pack_size_aq_spec_data = &
         pmc_mpi_pack_size_integer(val%n_spec) &
         + pmc_mpi_pack_size_string_array(val%name) &
         + pmc_mpi_pack_size_integer_array(val%pmc_gas_index) &
         + pmc_mpi_pack_size_integer_array(val%pmc_aero_index) &
         + pmc_mpi_pack_size_real_array(val%abstol) &
         + pmc_mpi_pack_size_logical_array(val%const_conc) &
         + pmc_mpi_pack_size_real_array(val%Dg) &
         + pmc_mpi_pack_size_real_array(val%N_star) &
         + pmc_mpi_pack_size_real_array(val%MW) &
         + pmc_mpi_pack_size_integer_array(val%charge)

  end function pmc_mpi_pack_size_aq_spec_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aq_spec_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_spec_data_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_spec)
    call pmc_mpi_pack_string_array(buffer, position, val%name)
    call pmc_mpi_pack_integer_array(buffer, position, val%pmc_gas_index)
    call pmc_mpi_pack_integer_array(buffer, position, val%pmc_aero_index)
    call pmc_mpi_pack_real_array(buffer, position, val%abstol)
    call pmc_mpi_pack_logical_array(buffer, position, val%const_conc)
    call pmc_mpi_pack_real_array(buffer, position, val%Dg)
    call pmc_mpi_pack_real_array(buffer, position, val%N_star)
    call pmc_mpi_pack_real_array(buffer, position, val%MW)
    call pmc_mpi_pack_integer_array(buffer, position, val%charge)
    call assert(608533616, &
         position - prev_position <= pmc_mpi_pack_size_aq_spec_data(val))
#endif

  end subroutine pmc_mpi_pack_aq_spec_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aq_spec_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_spec_data_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_spec)
    call pmc_mpi_unpack_string_array(buffer, position, val%name)
    call pmc_mpi_unpack_integer_array(buffer, position, val%pmc_gas_index)
    call pmc_mpi_unpack_integer_array(buffer, position, val%pmc_aero_index)
    call pmc_mpi_unpack_real_array(buffer, position, val%abstol)
    call pmc_mpi_unpack_logical_array(buffer, position, val%const_conc)
    call pmc_mpi_unpack_real_array(buffer, position, val%Dg)
    call pmc_mpi_unpack_real_array(buffer, position, val%N_star)
    call pmc_mpi_unpack_real_array(buffer, position, val%MW)
    call pmc_mpi_unpack_integer_array(buffer, position, val%charge)
    call assert(608651265, &
         position - prev_position <= pmc_mpi_pack_size_aq_spec_data(val))
#endif

  end subroutine pmc_mpi_unpack_aq_spec_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aq_spec_data





