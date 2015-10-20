! Copyright (C) 2011-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_integer_rmap2 module.

!> The integer_rmap2_t structure and assocated subroutines.
module pmc_integer_rmap2

  use pmc_integer_varray
  use pmc_util
  use pmc_mpi

  !> A map \f$\mathbb{Z} \to \mathbb{Z} \times \mathbb{Z}\f$, together
  !> with its multi-valued inverse.
  !!
  !! The forward map takes an integer \f$i\f$ in the domain \f$\{1,
  !! \ldots, N_{\rm d}\}\f$ to integers \f$(j_1,j_2)\f$ in the range
  !! \f$\{1, \ldots, N_{\rm r,1} \times \{1, \ldots, N_{\rm
  !! r,2}\}\f$. This is stored with <tt>j_1 =
  !! integer_rmap2%%forward1%%entry(i)</tt> and <tt>j_2 =
  !! integer_rmap2%%forward2%%entry(i)</tt>. This map is not assumed to
  !! be one-to-one or onto.
  !!
  !! The inverse map is multi-valued, with
  !! <tt>integer_rmap2%%inverse(j_1, j_2)</tt> containing all the inverses of
  !! \f$(j_1, j_2)\f$. The entry numbers in the inverse map are given by
  !! <tt>inverse_rmap%%index</tt>. The relationships between the
  !! forward and reverse maps are as follows.
  !!
  !! Given \c i, let:
  !! <pre>
  !! j_1 = integer_rmap2%%forward1%%entry(i)
  !! j_2 = integer_rmap2%%forward2%%entry(i)
  !! k = integer_rmap2%%index%%entry(i)
  !! </pre>
  !! Then:
  !! <pre>
  !! integer_rmap2%%inverse(j_1, j_2)%%entry(k) == i
  !! </pre>
  !!
  !! Alternatively, given \c j_1, \c j_2 and \c k, let:
  !! <pre>
  !! i = integer_rmap2%%inverse(j_1, j_2)%%entry(k)
  !! </pre>
  !! Then:
  !! <pre>
  !! integer_rmap2%%forward1%%entry(i) == j_1
  !! integer_rmap2%%forward2%%entry(i) == j_2
  !! integer_rmap2%%index%%entry(i) == k
  !! </pre>
  type integer_rmap2_t
     !> Forward map to first range (single valued).
     type(integer_varray_t) :: forward1
     !> Forward map to second range (single valued).
     type(integer_varray_t) :: forward2
     !> Inverse map (multi-valued).
     type(integer_varray_t), allocatable :: inverse(:, :)
     !> Forward map to inverse map entries (single valued).
     type(integer_varray_t) :: index
  end type integer_rmap2_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the maximum ranges of the forward map.
  elemental subroutine integer_rmap2_set_ranges(integer_rmap2, n_range_1, &
       n_range_2)

    !> Structure to initialize.
    type(integer_rmap2_t), intent(out) :: integer_rmap2
    !> Size of first range space.
    integer, intent(in) :: n_range_1
    !> Size of second range space.
    integer, intent(in) :: n_range_2

    if (allocated(integer_rmap2%inverse)) then
       deallocate(integer_rmap2%inverse)
    end if
    allocate(integer_rmap2%inverse(n_range_1, n_range_2))

  end subroutine integer_rmap2_set_ranges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an integer_rmap2 to have no mappings.
  elemental subroutine integer_rmap2_zero(integer_rmap2)

    !> Structure to zero.
    type(integer_rmap2_t), intent(inout) :: integer_rmap2

    call integer_varray_zero(integer_rmap2%forward1)
    call integer_varray_zero(integer_rmap2%forward2)
    if (allocated(integer_rmap2%inverse)) then
       call integer_varray_zero(integer_rmap2%inverse)
    end if
    call integer_varray_zero(integer_rmap2%index)

  end subroutine integer_rmap2_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the map value of the next free domain value to
  !> <tt>(i_range_1, i_range_2</tt>.
  subroutine integer_rmap2_append(integer_rmap2, i_range_1, i_range_2)

    !> Map to append to.
    type(integer_rmap2_t), intent(inout) :: integer_rmap2
    !> First range value.
    integer, intent(in) :: i_range_1
    !> Second range value.
    integer, intent(in) :: i_range_2

    call assert(438521606, allocated(integer_rmap2%inverse))
    call assert(708651144, i_range_1 >= 1)
    call assert(779828769, i_range_1 <= size(integer_rmap2%inverse, 1))
    call assert(978259336, i_range_2 >= 1)
    call assert(238981205, i_range_2 <= size(integer_rmap2%inverse, 2))

    ! grow map by one element
    call integer_varray_append(integer_rmap2%forward1, i_range_1)
    call integer_varray_append(integer_rmap2%forward2, i_range_2)
    call integer_varray_append(integer_rmap2%inverse(i_range_1, i_range_2), &
         integer_varray_n_entry(integer_rmap2%forward1))
    call integer_varray_append(integer_rmap2%index, &
         integer_varray_n_entry(integer_rmap2%inverse(i_range_1, i_range_2)))

  end subroutine integer_rmap2_append

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Change the map value of \c i_domain to <tt>(i_range_1, i_range_2)</tt>.
  subroutine integer_rmap2_change(integer_rmap2, i_domain, i_range_1, &
       i_range_2)

    !> Map to change.
    type(integer_rmap2_t), intent(inout) :: integer_rmap2
    !> Domain value.
    integer, intent(in) :: i_domain
    !> First range value.
    integer, intent(in) :: i_range_1
    !> Second range value.
    integer, intent(in) :: i_range_2

    integer :: i_range_1_old, i_range_2_old, i_index_old, i_domain_shifted

    call assert(897948211, allocated(integer_rmap2%inverse))
    call assert(191141591, i_domain >= 1)
    call assert(240079303, &
         i_domain <= integer_varray_n_entry(integer_rmap2%forward1))
    call assert(671426897, i_range_1 >= 1)
    call assert(311976942, i_range_1 <= size(integer_rmap2%inverse, 1))
    call assert(383129645, i_range_2 >= 1)
    call assert(771283685, i_range_2 <= size(integer_rmap2%inverse, 2))

    i_range_1_old = integer_rmap2%forward1%entry(i_domain)
    i_range_2_old = integer_rmap2%forward2%entry(i_domain)
    if ((i_range_1_old == i_range_1) .and. (i_range_1_old == i_range_1)) return
    i_index_old = integer_rmap2%index%entry(i_domain)

    ! remove the old inverse map
    call integer_varray_remove_entry( &
         integer_rmap2%inverse(i_range_1_old, i_range_2_old), i_index_old)
    if (i_index_old &
         <= integer_varray_n_entry(integer_rmap2%inverse(i_range_1_old, &
         i_range_2_old))) then
       ! the removed entry wasn't the last one, so the last entry
       ! was moved and needs fixing
       i_domain_shifted = integer_rmap2%inverse(i_range_1_old, &
            i_range_2_old)%entry(i_index_old)
       integer_rmap2%index%entry(i_domain_shifted) = i_index_old
    end if

    ! set the new map and inverse
    integer_rmap2%forward1%entry(i_domain) = i_range_1
    integer_rmap2%forward2%entry(i_domain) = i_range_2
    call integer_varray_append(integer_rmap2%inverse(i_range_1, i_range_2), &
         i_domain)
    integer_rmap2%index%entry(i_domain) &
         = integer_varray_n_entry(integer_rmap2%inverse(i_range_1, i_range_2))

  end subroutine integer_rmap2_change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Replace the map at the given \c i_domain with the map value of
  !> the last entry, and delete the last entry.
  subroutine integer_rmap2_remove(integer_rmap2, i_domain)

    !> Map to remove from.
    type(integer_rmap2_t), intent(inout) :: integer_rmap2
    !> Domain value to replace.
    integer, intent(in) :: i_domain

    integer :: i_range_1_old, i_range_2_old, i_index_old, i_domain_shifted
    integer :: i_range_1_fix, i_range_2_fix, i_index_fix, i_domain_fix

    call assert(902132756, allocated(integer_rmap2%inverse))
    call assert(242566612, i_domain >= 1)
    call assert(110569289, &
         i_domain <= integer_varray_n_entry(integer_rmap2%forward1))

    ! Deleting particles shifts the end particles into the empty slots
    ! in the aero_particle_array and the aero_sorted forward and
    ! reverse indexes. All must be fixed in the right order to
    ! maintain consistency.

    i_range_1_old = integer_rmap2%forward1%entry(i_domain)
    i_range_2_old = integer_rmap2%forward2%entry(i_domain)
    i_index_old = integer_rmap2%index%entry(i_domain)

    ! old shifted value loc
    i_domain_shifted = integer_varray_n_entry(integer_rmap2%forward1)
    if (i_domain_shifted /= i_domain) then
       i_range_1_fix = integer_rmap2%forward1%entry(i_domain_shifted)
       i_range_2_fix = integer_rmap2%forward2%entry(i_domain_shifted)
       i_index_fix = integer_rmap2%index%entry(i_domain_shifted)
       integer_rmap2%inverse(i_range_1_fix, i_range_2_fix)%entry(i_index_fix) &
            = i_domain
    end if

    ! remove the particle from the forward map (with the side effect
    ! of fixing the forward map for the shifted value)
    call integer_varray_remove_entry(integer_rmap2%forward1, i_domain)
    call integer_varray_remove_entry(integer_rmap2%forward2, i_domain)
    call integer_varray_remove_entry(integer_rmap2%index, i_domain)

    ! remove the inverse map
    i_index_fix = integer_varray_n_entry(integer_rmap2%inverse(i_range_1_old, &
         i_range_2_old))
    i_domain_fix = integer_rmap2%inverse(i_range_1_old, i_range_2_old)%entry(&
         i_index_fix)
    call integer_varray_remove_entry( &
         integer_rmap2%inverse(i_range_1_old, i_range_2_old), i_index_old)

    if (i_index_fix /= i_index_old) then
       ! fix index map
       integer_rmap2%index%entry(i_domain_fix) = i_index_old
    end if

  end subroutine integer_rmap2_remove

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the data is consistent.
  subroutine integer_rmap2_check(integer_rmap2, name, n_domain, n_range_1, &
       n_range_2, continue_on_error)

    !> Structure to check.
    type(integer_rmap2_t) :: integer_rmap2
    !> Check name.
    character(len=*), intent(in) :: name
    !> Size of domain.
    integer, intent(in) :: n_domain
    !> Size of first range.
    integer, intent(in) :: n_range_1
    !> Size of second range.
    integer, intent(in) :: n_range_2
    !> Whether to continue despite error.
    logical, intent(in) :: continue_on_error

    integer :: i_domain, i_range_1, i_range_2, i_index

    if (.not. allocated(integer_rmap2%inverse)) then
       return
    end if

    if ((n_domain /= integer_varray_n_entry(integer_rmap2%forward1)) &
         .or. (n_domain /= integer_varray_n_entry(integer_rmap2%forward2)) &
         .or. (n_domain /= integer_varray_n_entry(integer_rmap2%index)) &
         .or. (n_range_1 /= size(integer_rmap2%inverse, 1)) &
         .or. (n_range_2 /= size(integer_rmap2%inverse, 2))) then
       write(0,*) 'ERROR integer_rmap2 A:', name
       write(0,*) 'n_domain', n_domain
       write(0,*) 'n_range_1', n_range_1
       write(0,*) 'n_range_2', n_range_2
       write(0,*) 'integer_varray_n_entry(integer_rmap2%forward1)', &
            integer_varray_n_entry(integer_rmap2%forward1)
       write(0,*) 'integer_varray_n_entry(integer_rmap2%forward2)', &
            integer_varray_n_entry(integer_rmap2%forward2)
       write(0,*) 'integer_varray_n_entry(integer_rmap2%index)', &
            integer_varray_n_entry(integer_rmap2%index)
       write(0,*) 'size(integer_rmap2%inverse, 1)', &
            size(integer_rmap2%inverse, 1)
       write(0,*) 'size(integer_rmap2%inverse, 2)', &
            size(integer_rmap2%inverse, 2)
       call assert(786992107, continue_on_error)
    end if

    do i_domain = 1,n_domain
       i_range_1 = integer_rmap2%forward1%entry(i_domain)
       i_range_2 = integer_rmap2%forward2%entry(i_domain)
       if ((i_range_1 < 1) .or. (i_range_1 > n_range_1) &
            .or. (i_range_2 < 1) .or. (i_range_2 > n_range_2)) then
          write(0,*) 'ERROR integer_rmap2 B:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range_1', i_range_1
          write(0,*) 'i_range_2', i_range_2
          write(0,*) 'n_range_1', n_range_1
          write(0,*) 'n_range_2', n_range_2
          call assert(723392756, continue_on_error)
       end if

       i_index = integer_rmap2%index%entry(i_domain)
       if ((i_index < 1) .or. (i_index &
            > integer_varray_n_entry(integer_rmap2%inverse(i_range_1, &
            i_range_2)))) then
          write(0,*) 'ERROR integer_rmap2 C:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range_1', i_range_1
          write(0,*) 'i_range_2', i_range_2
          write(0,*) 'i_index', i_index
          write(0,*) 'integer_varray_n_entry(' &
               // 'integer_rmap2%inverse(i_range_1, i_range_2))', &
               integer_varray_n_entry(integer_rmap2%inverse(i_range_1, &
               i_range_2))
          call assert(317458796, continue_on_error)
       end if
       if (i_domain &
            /= integer_rmap2%inverse(i_range_1, i_range_2)%entry(i_index)) then
          write(0,*) 'ERROR integer_rmap2 D:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range_1', i_range_1
          write(0,*) 'i_range_2', i_range_2
          write(0,*) 'i_index', i_index
          write(0,*) 'integer_rmap2%inverse(i_range_1, ' &
               // 'i_range_2)%entry(i_index)', &
               integer_rmap2%inverse(i_range_1, i_range_2)%entry(i_index)
          call assert(662733308, continue_on_error)
       end if
    end do

    do i_range_1 = 1,n_range_1
       do i_range_2 = 1,n_range_2
          do i_index = 1,integer_varray_n_entry( &
               integer_rmap2%inverse(i_range_1, i_range_2))
             i_domain &
                  = integer_rmap2%inverse(i_range_1, i_range_2)%entry(i_index)
             if ((i_domain < 1) .or. (i_domain > n_domain)) then
                write(0,*) 'ERROR integer_rmap2 E:', name
                write(0,*) 'i_range_1', i_range_1
                write(0,*) 'i_range_2', i_range_2
                write(0,*) 'i_index', i_index
                write(0,*) 'i_domain', i_domain
                write(0,*) 'n_domain', n_domain
                call assert(639449827, continue_on_error)
             end if
             if ((i_range_1 /= integer_rmap2%forward1%entry(i_domain)) &
                  .or. (i_range_2 /= integer_rmap2%forward2%entry(i_domain)) &
                  .or. (i_index /= integer_rmap2%index%entry(i_domain))) then
                write(0,*) 'ERROR integer_rmap2 F:', name
                write(0,*) 'i_domain', i_domain
                write(0,*) 'i_range_1', i_range_1
                write(0,*) 'i_range_2', i_range_2
                write(0,*) 'integer_rmap2%forward1%entry(i_domain)', &
                     integer_rmap2%forward1%entry(i_domain)
                write(0,*) 'integer_rmap2%forward2%entry(i_domain)', &
                     integer_rmap2%forward2%entry(i_domain)
                write(0,*) 'i_index', i_index
                write(0,*) 'integer_rmap2%index%entry(i_domain)', &
                     integer_rmap2%index%entry(i_domain)
                call assert(636832060, continue_on_error)
             end if
          end do
       end do
    end do

  end subroutine integer_rmap2_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_integer_rmap2(val)

    !> Value to pack.
    type(integer_rmap2_t), intent(in) :: val

    integer :: i_1, i_2, total_size
    logical :: is_allocated

    total_size = 0
    is_allocated = allocated(val%inverse)
    total_size = total_size + pmc_mpi_pack_size_logical(is_allocated)
    if (is_allocated) then
       total_size = total_size &
            + pmc_mpi_pack_size_integer(size(val%inverse, 1))
       total_size = total_size &
            + pmc_mpi_pack_size_integer(size(val%inverse, 2))
       do i_1 = 1,size(val%inverse, 1)
          do i_2 = 1,size(val%inverse, 2)
             total_size = total_size &
                  + pmc_mpi_pack_size_integer_varray(val%inverse(i_1, i_2))
          end do
       end do
    end if
    total_size = total_size + pmc_mpi_pack_size_integer_varray(val%forward1)
    total_size = total_size + pmc_mpi_pack_size_integer_varray(val%forward2)
    total_size = total_size + pmc_mpi_pack_size_integer_varray(val%index)
    pmc_mpi_pack_size_integer_rmap2 = total_size

  end function pmc_mpi_pack_size_integer_rmap2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_integer_rmap2(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(integer_rmap2_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i_1, i_2
    logical :: is_allocated

    prev_position = position
    is_allocated = allocated(val%inverse)
    call pmc_mpi_pack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       call pmc_mpi_pack_integer(buffer, position, size(val%inverse, 1))
       call pmc_mpi_pack_integer(buffer, position, size(val%inverse, 2))
       do i_1 = 1,size(val%inverse, 1)
          do i_2 = 1,size(val%inverse, 2)
             call pmc_mpi_pack_integer_varray(buffer, position, &
                  val%inverse(i_1, i_2))
          end do
       end do
    end if
    call pmc_mpi_pack_integer_varray(buffer, position, val%forward1)
    call pmc_mpi_pack_integer_varray(buffer, position, val%forward2)
    call pmc_mpi_pack_integer_varray(buffer, position, val%index)
    call assert(283629348, &
         position - prev_position <= pmc_mpi_pack_size_integer_rmap2(val))
#endif

  end subroutine pmc_mpi_pack_integer_rmap2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_integer_rmap2(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(integer_rmap2_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i_1, i_2, n_1, n_2
    logical :: is_allocated

    prev_position = position
    call pmc_mpi_unpack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       call pmc_mpi_unpack_integer(buffer, position, n_1)
       call pmc_mpi_unpack_integer(buffer, position, n_2)
       call integer_rmap2_set_ranges(val, n_1, n_2)
       do i_1 = 1,size(val%inverse, 1)
          do i_2 = 1,size(val%inverse, 2)
             call pmc_mpi_unpack_integer_varray(buffer, position, &
                  val%inverse(i_1, i_2))
          end do
       end do
    else
       if (allocated(val%inverse)) then
          deallocate(val%inverse)
       end if
    end if
    call pmc_mpi_unpack_integer_varray(buffer, position, val%forward1)
    call pmc_mpi_unpack_integer_varray(buffer, position, val%forward2)
    call pmc_mpi_unpack_integer_varray(buffer, position, val%index)
    call assert(796602256, &
         position - prev_position <= pmc_mpi_pack_size_integer_rmap2(val))
#endif

  end subroutine pmc_mpi_unpack_integer_rmap2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_integer_rmap2
