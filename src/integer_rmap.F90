! Copyright (C) 2011-2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_integer_rmap module.

!> The integer_rmap_t structure and assocated subroutines.
module pmc_integer_rmap

  use pmc_integer_varray
  use pmc_util
  use pmc_mpi

  !> A map from integers to integers, together with its multi-valued
  !> inverse.
  !!
  !! The forward map takes integer \f$i\f$ in the domain
  !! 1,...,n_domain to an integer \f$j\f$ in the range
  !! 1,...,n_range. This is stored with <tt>j =
  !! integer_rmap%%forward%%entry(i)</tt>. This map will generally not be
  !! one-to-one or onto.
  !!
  !! The inverse map is multi-valued, with
  !! <tt>integer_rmap%%inverse(j)</tt> containing all the inverses of
  !! \f$j\f$. The entries in the inverse map are given by
  !! <tt>inverse_rmap%%index</tt>. The relationships between
  !! the forward and reverse maps are as follows.
  !!
  !! Given \f$i\f$, let:
  !! <pre>
  !! j = integer_rmap%%forward%%entry(i)
  !! k = integer_rmap%%index%%entry(i)
  !! </pre>
  !! Then:
  !! <pre>
  !! integer_rmap%%inverse(j)%%entry(k) == i
  !! </pre>
  !!
  !! Alternatively, given \f$j\f$ and \f$k\f$, let:
  !! <pre>
  !! i = integer_rmap%%inverse(j)%%entry(k)
  !! </pre>
  !! Then:
  !! <pre>
  !! integer_rmap%%forward%%entry(i) == j
  !! integer_rmap%%index%%entry(i) == k
  !! </pre>
  type integer_rmap_t
     !> Forward map (single valued).
     type(integer_varray_t) :: forward
     !> Inverse map (multi-valued).
     type(integer_varray_t), allocatable :: inverse(:)
     !> Forward map to inverse map entries (single valued).
     type(integer_varray_t) :: index
  end type integer_rmap_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Sets the maximum range of the forward map.
  elemental subroutine integer_rmap_set_range(integer_rmap, n_range)

    !> Structure to initialize.
    type(integer_rmap_t), intent(out) :: integer_rmap
    !> Size of range space.
    integer, intent(in) :: n_range

    if (allocated(integer_rmap%inverse)) then
       deallocate(integer_rmap%inverse)
    end if
    allocate(integer_rmap%inverse(n_range))

  end subroutine integer_rmap_set_range

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets an integer_rmap to have no mappings.
  elemental subroutine integer_rmap_zero(integer_rmap)

    !> Structure to zero.
    type(integer_rmap_t), intent(inout) :: integer_rmap

    call integer_varray_zero(integer_rmap%forward)
    if (allocated(integer_rmap%inverse)) then
       call integer_varray_zero(integer_rmap%inverse)
    end if
    call integer_varray_zero(integer_rmap%index)

  end subroutine integer_rmap_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the map value of the next free domain value to \c i_range.
  subroutine integer_rmap_append(integer_rmap, i_range)

    !> Map to append to.
    type(integer_rmap_t), intent(inout) :: integer_rmap
    !> Range value.
    integer, intent(in) :: i_range

    call assert(200002696, allocated(integer_rmap%inverse))
    call assert(549740445, i_range >= 1)
    call assert(145872613, i_range <= size(integer_rmap%inverse))

    ! grow map by one element
    call integer_varray_append(integer_rmap%forward, i_range)
    call integer_varray_append(integer_rmap%inverse(i_range), &
         integer_varray_n_entry(integer_rmap%forward))
    call integer_varray_append(integer_rmap%index, &
         integer_varray_n_entry(integer_rmap%inverse(i_range)))

  end subroutine integer_rmap_append

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Change the map value of \c i_domain to \c i_range.
  subroutine integer_rmap_change(integer_rmap, i_domain, i_range)

    !> Map to change.
    type(integer_rmap_t), intent(inout) :: integer_rmap
    !> Domain value.
    integer, intent(in) :: i_domain
    !> Range value.
    integer, intent(in) :: i_range

    integer :: i_range_old, i_index_old, i_domain_shifted

    call assert(483243151, allocated(integer_rmap%inverse))
    call assert(709581778, i_domain >= 1)
    call assert(494892311, &
         i_domain <= integer_varray_n_entry(integer_rmap%forward))

    call assert(590911054, i_range >= 1)
    call assert(859774512, i_range <= size(integer_rmap%inverse))

    i_range_old = integer_rmap%forward%entry(i_domain)
    if (i_range_old == i_range) return
    i_index_old = integer_rmap%index%entry(i_domain)

    ! remove the old inverse map
    call integer_varray_remove_entry(integer_rmap%inverse(i_range_old), &
         i_index_old)
    if (i_index_old &
         <= integer_varray_n_entry(integer_rmap%inverse(i_range_old))) then
       ! the removed entry wasn't the last one, so the last entry
       ! was moved and needs fixing
       i_domain_shifted = integer_rmap%inverse(i_range_old)%entry(i_index_old)
       integer_rmap%index%entry(i_domain_shifted) = i_index_old
    end if

    ! set the new map and inverse
    integer_rmap%forward%entry(i_domain) = i_range
    call integer_varray_append(integer_rmap%inverse(i_range), i_domain)
    integer_rmap%index%entry(i_domain) &
         = integer_varray_n_entry(integer_rmap%inverse(i_range))

  end subroutine integer_rmap_change

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Replace the map at the given \c i_domain with the map value of
  !> the last entry, and delete the last entry.
  subroutine integer_rmap_remove(integer_rmap, i_domain)

    !> Map to remove from.
    type(integer_rmap_t), intent(inout) :: integer_rmap
    !> Domain value to replace.
    integer, intent(in) :: i_domain

    integer :: i_range_old, i_index_old, i_domain_shifted, i_range_fix
    integer :: i_index_fix, i_domain_fix

    call assert(790455753, allocated(integer_rmap%inverse))
    call assert(745161821, i_domain >= 1)
    call assert(143043782, i_domain &
         <= integer_varray_n_entry(integer_rmap%forward))

    ! Deleting particles shifts the end particles into the empty slots
    ! in the aero_particle_array and the aero_sorted forward and
    ! reverse indexes. All must be fixed in the right order to
    ! maintain consistency.

    i_range_old = integer_rmap%forward%entry(i_domain)
    i_index_old = integer_rmap%index%entry(i_domain)

    ! old loc of shifted value
    i_domain_shifted = integer_varray_n_entry(integer_rmap%forward)
    if (i_domain_shifted /= i_domain) then
       i_range_fix = integer_rmap%forward%entry(i_domain_shifted)
       i_index_fix = integer_rmap%index%entry(i_domain_shifted)
       integer_rmap%inverse(i_range_fix)%entry(i_index_fix) = i_domain
    end if

    ! remove the particle from the forward map (with the side effect
    ! of fixing the forward map for the shifted value)
    call integer_varray_remove_entry(integer_rmap%forward, i_domain)
    call integer_varray_remove_entry(integer_rmap%index, i_domain)

    ! remove the inverse map
    i_index_fix = integer_varray_n_entry(integer_rmap%inverse(i_range_old))
    i_domain_fix = integer_rmap%inverse(i_range_old)%entry(i_index_fix)
    call integer_varray_remove_entry(integer_rmap%inverse(i_range_old), &
         i_index_old)

    if (i_index_fix /= i_index_old) then
       ! fix index map
       integer_rmap%index%entry(i_domain_fix) = i_index_old
    end if

  end subroutine integer_rmap_remove

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check that the data is consistent.
  subroutine integer_rmap_check(integer_rmap, name, n_domain, n_range, &
       continue_on_error)

    !> Structure to check.
    type(integer_rmap_t) :: integer_rmap
    !> Check name.
    character(len=*), intent(in) :: name
    !> Number of domain items.
    integer, intent(in) :: n_domain
    !> Number of image items.
    integer, intent(in) :: n_range
    !> Whether to continue despite error.
    logical, intent(in) :: continue_on_error

    integer :: i_domain, i_range, i_index

    if (.not. allocated(integer_rmap%inverse)) then
       return
    end if

    if ((n_domain /= integer_varray_n_entry(integer_rmap%forward)) &
         .or. (n_domain /= integer_varray_n_entry(integer_rmap%index)) &
         .or. (n_range /= size(integer_rmap%inverse))) then
       write(0,*) 'ERROR integer_rmap A:', name
       write(0,*) 'n_domain', n_domain
       write(0,*) 'n_range', n_range
       write(0,*) 'integer_varray_n_entry(integer_rmap%forward)', &
            integer_varray_n_entry(integer_rmap%forward)
       write(0,*) 'integer_varray_n_entry(integer_rmap%index)', &
            integer_varray_n_entry(integer_rmap%index)
       write(0,*) 'size(integer_rmap%inverse)', size(integer_rmap%inverse)
       call assert(973643016, continue_on_error)
    end if

    do i_domain = 1,n_domain
       i_range = integer_rmap%forward%entry(i_domain)
       if ((i_range < 1) .or. (i_range > n_range)) then
          write(0,*) 'ERROR integer_rmap B:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range', i_range
          write(0,*) 'n_range', n_range
          call assert(798857945, continue_on_error)
       end if

       i_index = integer_rmap%index%entry(i_domain)
       if ((i_index < 1) &
            .or. (i_index &
            > integer_varray_n_entry(integer_rmap%inverse(i_range)))) then
          write(0,*) 'ERROR integer_rmap C:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range', i_range
          write(0,*) 'i_index', i_index
          write(0,*) 'integer_varray_n_entry(integer_rmap%inverse(i_range))', &
               integer_varray_n_entry(integer_rmap%inverse(i_range))
          call assert(823748734, continue_on_error)
       end if
       if (i_domain /= integer_rmap%inverse(i_range)%entry(i_index)) then
          write(0,*) 'ERROR integer_rmap D:', name
          write(0,*) 'i_domain', i_domain
          write(0,*) 'i_range', i_range
          write(0,*) 'i_index', i_index
          write(0,*) 'integer_rmap%inverse(i_range)%entry(i_index)', &
               integer_rmap%inverse(i_range)%entry(i_index)
          call assert(735205557, continue_on_error)
       end if
    end do

    do i_range = 1,n_range
       do i_index = 1,integer_varray_n_entry(integer_rmap%inverse(i_range))
          i_domain = integer_rmap%inverse(i_range)%entry(i_index)
          if ((i_domain < 1) .or. (i_domain > n_domain)) then
             write(0,*) 'ERROR integer_rmap E:', name
             write(0,*) 'i_range', i_range
             write(0,*) 'i_index', i_index
             write(0,*) 'i_domain', i_domain
             write(0,*) 'n_domain', n_domain
             call assert(502643520, continue_on_error)
          end if
          if ((i_range /= integer_rmap%forward%entry(i_domain)) &
               .or. (i_index /= integer_rmap%index%entry(i_domain))) then
             write(0,*) 'ERROR integer_rmap F:', name
             write(0,*) 'i_domain', i_domain
             write(0,*) 'i_range', i_range
             write(0,*) 'integer_rmap%forward%entry(i_domain)', &
                  integer_rmap%forward%entry(i_domain)
             write(0,*) 'i_index', i_index
             write(0,*) 'integer_rmap%index%entry(i_domain)', &
                  integer_rmap%index%entry(i_domain)
             call assert(544747928, continue_on_error)
          end if
       end do
    end do

  end subroutine integer_rmap_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_integer_rmap(val)

    !> Value to pack.
    type(integer_rmap_t), intent(in) :: val

    integer :: i, total_size
    logical :: is_allocated

    total_size = 0
    is_allocated = allocated(val%inverse)
    total_size = total_size + pmc_mpi_pack_size_logical(is_allocated)
    if (is_allocated) then
       total_size = total_size + pmc_mpi_pack_size_integer(size(val%inverse))
       do i = 1,size(val%inverse)
          total_size = total_size &
               + pmc_mpi_pack_size_integer_varray(val%inverse(i))
       end do
    end if
    total_size = total_size + pmc_mpi_pack_size_integer_varray(val%forward)
    total_size = total_size + pmc_mpi_pack_size_integer_varray(val%index)
    pmc_mpi_pack_size_integer_rmap = total_size

  end function pmc_mpi_pack_size_integer_rmap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_integer_rmap(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(integer_rmap_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i
    logical :: is_allocated

    prev_position = position
    is_allocated = allocated(val%inverse)
    call pmc_mpi_pack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       call pmc_mpi_pack_integer(buffer, position, size(val%inverse))
       do i = 1,size(val%inverse)
          call pmc_mpi_pack_integer_varray(buffer, position, val%inverse(i))
       end do
    end if
    call pmc_mpi_pack_integer_varray(buffer, position, val%forward)
    call pmc_mpi_pack_integer_varray(buffer, position, val%index)
    call assert(533568488, &
         position - prev_position <= pmc_mpi_pack_size_integer_rmap(val))
#endif

  end subroutine pmc_mpi_pack_integer_rmap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_integer_rmap(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(integer_rmap_t), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n
    logical :: is_allocated

    prev_position = position
    call pmc_mpi_unpack_logical(buffer, position, is_allocated)
    if (is_allocated) then
       call pmc_mpi_unpack_integer(buffer, position, n)
       call integer_rmap_set_range(val, n)
       do i = 1,size(val%inverse)
          call pmc_mpi_unpack_integer_varray(buffer, position, val%inverse(i))
       end do
    else
       if (allocated(val%inverse)) then
          deallocate(val%inverse)
       end if
    end if
    call pmc_mpi_unpack_integer_varray(buffer, position, val%forward)
    call pmc_mpi_unpack_integer_varray(buffer, position, val%index)
    call assert(663161025, &
         position - prev_position <= pmc_mpi_pack_size_integer_rmap(val))
#endif

  end subroutine pmc_mpi_unpack_integer_rmap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_integer_rmap
