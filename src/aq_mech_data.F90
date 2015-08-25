! Copyright (C) 2015 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aq_mech_data module.

!> The aq_mech_data_t structure and associated subroutines.

module pmc_aq_mech_data

  use pmc_aq_rxn_data
  use pmc_aq_rxn_file
#ifdef PMC_USE_MPI
  use mpi
#endif

  implicit none

  !> Constant aqueous reaction mechanism data
  !!
  !! Each aqueous reaction is identified by an integer \c i between 1 
  !! and \c n_rxn.
  type aq_mech_data_t
     !> Number of reactions.
     integer :: n_rxn
     !> Reaction data
     type(aq_rxn_data_t), pointer :: rxn(:)
  end type aq_mech_data_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for aqueous chemistry mechanism
  subroutine aq_mech_data_allocate(aq_mech_data)

    !> Aqueous chemistry mechanism data.
    type(aq_mech_data_t), intent(out) :: aq_mech_data

    aq_mech_data%n_rxn = 0
    allocate(aq_mech_data%rxn(0))

  end subroutine aq_mech_data_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocate storage for aqueous chemistry mechanism with the given size.
  subroutine aq_mech_data_allocate_size(aq_mech_data, n_rxn)

    !> Aqueous chemistry mechanism data.
    type(aq_mech_data_t), intent(out) :: aq_mech_data
    !> Number of reactions.
    integer, intent(in) :: n_rxn

    aq_mech_data%n_rxn = n_rxn
    allocate(aq_mech_data%rxn(n_rxn))

  end subroutine aq_mech_data_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aq_mech_data_deallocate(aq_mech_data)

    !> Aqueous chemistry mechanism data.
    type(aq_mech_data_t), intent(inout) :: aq_mech_data

    integer :: i

    ! Deallocate each reaction in mechanism
    do i=1,size(aq_mech_data%rxn)
        call aq_rxn_data_deallocate(aq_mech_data%rxn(i))
    enddo

    deallocate(aq_mech_data%rxn)

  end subroutine aq_mech_data_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update an aqueous mechanism with new species indices
  subroutine aq_mech_data_spec_update(aq_mech_data, spec_map)

    !> Aqueous chemistry mechanism data.
    type(aq_mech_data_t), intent(inout) :: aq_mech_data
    !> Species map
    !! [(spec_map(i))_new = (i)_old]
    integer, intent(in) :: spec_map(:)

    integer :: i, j

    do i=1,size(aq_mech_data%rxn)
        do j=1,size(aq_mech_data%rxn(i)%reactant)
            aq_mech_data%rxn(i)%reactant(j) = spec_map(aq_mech_data%rxn(i)%reactant(j))
        enddo
        do j=1,size(aq_mech_data%rxn(i)%product)
            aq_mech_data%rxn(i)%product(j) = spec_map(aq_mech_data%rxn(i)%product(j))
        enddo
    enddo

  end subroutine aq_mech_data_spec_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read aqueous chemistry mechanism from a CAPRAM file
  subroutine aq_mech_file_read_data(file, aq_mech_data, aq_spec_data)

    !> CAPRAM file to read data from
    type(aq_rxn_file_t), intent(inout) :: file
    !> Aqueous chemistry mechanism
    type(aq_mech_data_t), intent(inout) :: aq_mech_data
    !> Aqueous chemistry related species data.
    type(aq_spec_data_t), intent(inout) :: aq_spec_data

    integer :: n_rxn, i

    !> \page input_format_aq_mech_data Input File Format: Aqueous Chemistry Mechanism Data
    !!
    !! Aqueous-Phase reaction files do not follow the standard 
    !! \ref spec_file_format but are instead in the CAPRAM text file 
    !! format, described here. Individual reactions require exactly
    !! three lines each, and follow the format:
    !! <pre>
    !!   CLASS: class_name
    !!   rct1 + rct2 + rct3 = [yield1] prod1 + [yield] prod2 + [yield3] prod3
    !!   rc_type: A: 123.4e-15 B: 123.0
    !!</pre>
    !! Here, lowercase terms and numeric values should be replaced by
    !! reaction parameters. Uppercase terms are keywords and must appear
    !! as shown. Terms enclosed in square brackets are optional.
    !! Text parameters are case sensitive and cannot
    !! contain spaces. By CAPRAM convention, text parameters are 
    !! uppercase. There should be no extra lines in between the three lines
    !! required for each reaction. The parameters are:
    !!
    !!   - \c class_name:       Type of reaction. Determines how mass
    !!                       transfer is calculated from this reaction.
    !!   - \c rct1, \c rct2, etc.: Reactant species. There can be any number
    !!                       of reactants.
    !!   - \c prod1, \c prod2, etc.: Product species. There can be any number
    !!                       of products.
    !!   - \c yield1, \c yield2, etc. : Yields may be included when the
    !!                       rate constant type does not include a backward 
    !!                       reaction
    !!   - \c rc_type:       Rate constant type. Determines form of rate
    !!                       constant expression. Must be accompanied by
    !!                       correct number of terms following "A:" "B:"
    !!                       "C:" etc.
    !!
    !! The keywords \c CLASS:, \c +, \c =, \c A:, and \c B: must be present as
    !! shown, although the number of rate constant parameters must
    !! match that required by \c rc_type and follow the format \c A:, \c B:,
    !! \c C:, etc. Keyword \c + can be replaced by \c - in product terms to
    !! specify consumption of a reactant without its concentration
    !! affecting the reaction rate.
    !!
    !! Species names that are surrounded by square brackets will be
    !! treated as constant in the aqueous chemistry mechanism.
    !!
    !! Comments may be included after the \c COMMENT keyword.
    !!
    !! For example, an aqueous chemistry mechanism could include the 
    !! following lines for a reaction:
    !! <pre>
    !! CLASS: AQUA
    !! O2m  +  CUp   = CUpp  + aH2O2   + 2 OHm - 2 [aH2O]
    !! TEMP3:   A:  1e10  B: 0.0
    !!</pre>
    !! (Note that it is not really necessary to include the consumption term
    !! for \c [aH2O] since it is a constant concentration species.)
    !!
    !! The aqueous chemistry mechanism file is specified by the parameter:
    !!   - \b aq_mech (string): name of file from which to read the
    !!   aqueous chemistry mechanism
    !!
    !! Reference:
    !! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended 
    !! and condensed tropospheric aqueous phase mechanism and its 
    !! application."" J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202


    ! Get number of reactions in CAPRAM file
    call aq_rxn_file_count_rxns(file, n_rxn)

    ! Allocate space for the whole mechanism
    call aq_mech_data_allocate_size(aq_mech_data, n_rxn)

    ! Read in reactions
    do i=1,n_rxn
        call aq_rxn_file_read_rxn(file, aq_mech_data%rxn(i), aq_spec_data)
    enddo

    ! Reorder species for efficient solving
    call aq_mech_data_reorder(aq_mech_data, aq_spec_data)

  end subroutine aq_mech_file_read_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print out a reaction mechanism
  subroutine aq_mech_data_print(aq_mech_data, aq_spec_data)

    !> Aqueous chemistry mechanism
    type(aq_mech_data_t), intent(in) :: aq_mech_data
    !> Aqueous chemistry related species data.
    type(aq_spec_data_t), intent(in) :: aq_spec_data

    character(len=100) :: temp_string, temp_elem_string

    integer :: i, j

    real(kind=dp) :: rc_forward, rc_backward
    real(kind=dp) :: ret_val

    write(*,*) "Aqueous-Phase Chemistry Mechanism"
    write(*,*) " "

    do i=1,aq_mech_data%n_rxn

        write(*,'(A,I4,A,A)') "Reaction ",i, "   Class: ", &
            trim(aq_rxn_data_get_class_name(aq_mech_data%rxn(i)%class_index))

        write(*,'(A)') "Reactants: "
        temp_string = ""
        do j=1,size(aq_mech_data%rxn(i)%reactant)
            temp_elem_string = trim(aq_spec_data%name(aq_mech_data%rxn(i)%reactant(j)))
            if ((len(trim(temp_elem_string))+len(trim(temp_string))).gt.99) then
                write (*,'   (A)') trim(temp_string)
                temp_string = trim(temp_elem_string)
            else
                temp_string = trim(temp_string) // " " // trim(temp_elem_string)
            endif
        enddo
        write (*,'   (A)') trim(temp_string)

        write(*,'(A)') "Products: "
        temp_string = ""
        do j=1,size(aq_mech_data%rxn(i)%product)
            write(temp_elem_string, '(g12.5,A,A)') aq_mech_data%rxn(i)%prod_yield(j), " ", &
                        trim(aq_spec_data%name(aq_mech_data%rxn(i)%product(j)))
            if (len(trim(temp_elem_string))+len(trim(temp_string)).gt.99) then
                write (*,'   (A)') temp_string
                temp_string = trim(temp_elem_string)
            else
                temp_string = trim(temp_string) // " " // trim(temp_elem_string)
            endif
        enddo
        write (*,'   (A)') trim(temp_string)

        write(*,'(A,A,A)') "Rate Parameters for type ", &
                        trim(aq_rxn_data_get_rate_name(aq_mech_data%rxn(i)%rate_index))
        temp_string = ""
        do j=1,size(aq_mech_data%rxn(i)%rate_param)
            write(temp_elem_string, '(g12.5)') aq_mech_data%rxn(i)%rate_param(j)
            if (len(trim(temp_elem_string))+len(trim(temp_string)).gt.99) then
                write (*,'   (A)') temp_string
                temp_string = trim(temp_elem_string)
            else
                temp_string = trim(temp_string) // " " // trim(temp_elem_string)
            endif
        enddo
        write (*,'   (A)') trim(temp_string)

        write(*,*) "k_forward      k_backward     (298K, SZA=0, r=100nm)"
        ret_val = aq_rxn_data_get_rate_constant(rc_forward, rc_backward, &
            aq_mech_data%rxn(i), aq_spec_data, real(298.0,dp), real(0.0,dp), real(1.0e-7,dp))
        write(*,"(g10.5,A6,g10.5)") rc_forward, " ", rc_backward
        if (trim(aq_rxn_data_get_class_name(aq_mech_data%rxn(i)%class_index)) &
            .eq. "HENRY") then
            write(*,"(A,g10.5)") "Henry's Law Eq. Constant @ 298K (M/atm) = ", rc_forward/rc_backward
        endif
        write(*,*) " "

    enddo

  end subroutine aq_mech_data_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that all reactions that require species physical constants have them
  subroutine aq_mech_data_check_const(aq_mech_data, aq_spec_data)

    !> Aqueous chemistry mechanism
    type(aq_mech_data_t), intent(in) :: aq_mech_data
    !> Aqueous chemistry related species data.
    type(aq_spec_data_t), intent(in) :: aq_spec_data

    integer :: i, species_index

    do i=1,aq_mech_data%n_rxn

        select case(trim(aq_rxn_data_get_class_name(aq_mech_data%rxn(i)%class_index)))

            case("HENRY")
                species_index = aq_mech_data%rxn(i)%reactant(1)
                if (aq_spec_data%Dg(species_index).eq.0.0) then
                    call die_msg(614298401, 'species ' // trim(aq_spec_data%name(species_index)) &
                        // ' is a partitioning species, but has no gas-phase diffusion constant.')
                endif
                if (aq_spec_data%MW(species_index).eq.0.0) then
                    call die_msg(611038401, 'species ' // trim(aq_spec_data%name(species_index)) &
                        // ' is a partitioning species, but has no molecular weight.')
                endif
                if (aq_spec_data%N_star(species_index).eq.0.0) then
                    call die_msg(610378401, 'species ' // trim(aq_spec_data%name(species_index)) &
                        // ' is a partitioning species, but has no N* parameter.')
                endif

            case default

        end select

    enddo

  end subroutine aq_mech_data_check_const

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aq_mech_data(val)

    !> Value to pack.
    type(aq_mech_data_t), intent(in) :: val

    integer :: i

    pmc_mpi_pack_size_aq_mech_data = &
         pmc_mpi_pack_size_integer(val%n_rxn)

    do i=1, val%n_rxn
        pmc_mpi_pack_size_aq_mech_data = pmc_mpi_pack_size_aq_mech_data &
            + pmc_mpi_pack_size_aq_rxn_data(val%rxn(i))
    enddo

  end function pmc_mpi_pack_size_aq_mech_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aq_mech_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_mech_data_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_pack_integer(buffer, position, val%n_rxn)

    do i=1, val%n_rxn
        call pmc_mpi_pack_aq_rxn_data(buffer, position, val%rxn(i))
    enddo

    call assert(599625906, &
         position - prev_position <= pmc_mpi_pack_size_aq_mech_data(val))
#endif

  end subroutine pmc_mpi_pack_aq_mech_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aq_mech_data(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aq_mech_data_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, val%n_rxn)
    allocate (val%rxn(val%n_rxn))

    do i=1, val%n_rxn
        call aq_rxn_data_allocate(val%rxn(i))
        call pmc_mpi_unpack_aq_rxn_data(buffer, position, val%rxn(i))
    enddo
    call assert(599777169, &
         position - prev_position <= pmc_mpi_pack_size_aq_mech_data(val))
#endif

  end subroutine pmc_mpi_unpack_aq_mech_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reorder species to minimize fill-in during LU decomposition
  subroutine aq_mech_data_reorder(aq_mech_data, aq_spec_data)

    !> Aqueous chemistry mechanism
    type(aq_mech_data_t), intent(inout) :: aq_mech_data
    !> Aqueous chemistry related species data
    type(aq_spec_data_t), intent(inout) :: aq_spec_data

    ! Old species indices
    integer, allocatable :: old_spec_index(:)
    ! New species indices
    integer, allocatable :: new_spec_index(:)
    ! Partitioning species
    logical, allocatable :: partitioning_spec(:)
    ! Sparsity pattern
    integer, allocatable :: sparsity(:,:)
    ! Temporary aq_spec_data_t variable for reordering
    type(aq_spec_data_t) :: temp_aq_spec_data

    integer :: i, j, k, n_spec, n_rxn, beta, beta_low, beta_low_index
    integer :: c_k, r_k

    !> Reordering is based on method described in:
    !! Wolke et al, (2005) SPACCIM: A parcel model with detailed 
    !! microphysics and complex multiphase chemistry, <i>Atmos. 
    !! Environ.</i> 39, 4375-4388.
    !! and in:
    !! Sandu, et al. (1996) Efficient implementation of fully 
    !! implicit methods for atmospheric chemical kinetics,
    !! <i>J. Comput. Phys.</i> 129, 101-110.

    n_spec = aq_spec_data%n_spec
    n_rxn = aq_mech_data%n_rxn

    allocate(old_spec_index(n_spec))
    allocate(new_spec_index(n_spec))
    allocate(partitioning_spec(n_spec))
    allocate(sparsity(n_spec, n_spec))

    do i=1, n_spec
        old_spec_index(i) = i
    enddo
    partitioning_spec(:) = .false.
    sparsity(:,:) = 0

    ! Calculate sparsity pattern
    do i=1, n_rxn
        do j=1, size(aq_mech_data%rxn(i)%reactant)
            sparsity(aq_mech_data%rxn(i)%reactant(j), aq_mech_data%rxn(i)%reactant(:)) = 1
        enddo
        do j=1, size(aq_mech_data%rxn(i)%product)
            sparsity(aq_mech_data%rxn(i)%product(j), aq_mech_data%rxn(i)%reactant(:)) = 1
        enddo
        if (aq_rxn_data_is_backward_rxn(aq_mech_data%rxn(i))) then
            do j=1, size(aq_mech_data%rxn(i)%reactant)
                sparsity(aq_mech_data%rxn(i)%reactant(j), aq_mech_data%rxn(i)%product(:)) = 1
            enddo
            do j=1, size(aq_mech_data%rxn(i)%product)
                sparsity(aq_mech_data%rxn(i)%product(j), aq_mech_data%rxn(i)%product(:)) = 1
            enddo
        endif
        ! identify gas-phase species that participate in HL partitioning
        if (trim(aq_rxn_data_get_class_name(aq_mech_data%rxn(i)%class_index)).eq."HENRY") then
            partitioning_spec(aq_mech_data%rxn(i)%reactant(1)) = .true.
        endif
    enddo

    ! Diagonal Markowitz
    do k=1, n_spec

        beta_low = 0
        beta_low_index = 0
        do i=k, n_spec

            ! calculate beta_k (Sandu et al. pg 102)
            r_k = 0
            do j=k, n_spec
                r_k = r_k + sparsity(j,i)
            enddo
            c_k = 0
            do j=k, n_spec
                c_k = c_k + sparsity(i,j)
            enddo
            beta = (r_k-1) * (c_k-1)

            ! ensure gas-phase species that participate in HL partitioning
            ! are at the end of the pivot order
            if (partitioning_spec(old_spec_index(i))) then
                beta = 1e6
            endif

            ! save index corresponding to lowest beta
            if (beta_low_index.eq.0 .or. beta.lt.beta_low) then
                beta_low = beta
                beta_low_index = i
            endif

        enddo

        ! Swap species i and k if necessary
        if (beta_low_index.ne.k) then

            ! update species map
            j = old_spec_index(beta_low_index)
            old_spec_index(beta_low_index) = old_spec_index(k)
            old_spec_index(k) = j

            ! swap matrix columns and rows
            do i=1,n_spec
                j = sparsity(beta_low_index,i)
                sparsity(beta_low_index,i) = sparsity(k,i)
                sparsity(k,i) = j
                j = sparsity(i,beta_low_index)
                sparsity(i,beta_low_index) = sparsity(i,k)
                sparsity(i,k) = j
            enddo

        endif

    enddo

    ! create new_species_index
    new_spec_index(:) = 0
    do i=1,n_spec
        do j=1,n_spec
            if (old_spec_index(j).eq.i) then
                if (new_spec_index(i).ne.0) then
                    call die_msg(600701554, 'internal error reordering species.')
                endif
                new_spec_index(i) = j
            endif
        enddo
        if (new_spec_index(i).eq.0) then
            call die_msg(600886431, 'internal error reordering species.')
        endif
    enddo

    ! reorder species in aq_spec_data
    call aq_spec_data_allocate_size(temp_aq_spec_data, n_spec)
    call aq_spec_data_copy_reorder(aq_spec_data, temp_aq_spec_data, new_spec_index, n_spec)
    call aq_spec_data_copy(temp_aq_spec_data, aq_spec_data, n_spec)
    call aq_spec_data_deallocate(temp_aq_spec_data)

    ! update mechanism
    call aq_mech_data_spec_update(aq_mech_data, new_spec_index)

    deallocate(old_spec_index)
    deallocate(new_spec_index)
    deallocate(partitioning_spec)
    deallocate(sparsity)

  end subroutine aq_mech_data_reorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aq_mech_data









