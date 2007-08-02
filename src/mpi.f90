! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Wrapper functions for MPI. All of these functions can be called
! irrespective of whether MPI support was compiled in or not. If MPI
! support is not enabled then they do the obvious trivial thing
! (normally nothing).

module pmc_mpi

#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_check_ierr(ierr)

    ! Dies if ierr is not ok.

    integer, intent(in) :: ierr         ! MPI status code

#ifdef PMC_USE_MPI
    if (ierr /= MPI_SUCCESS) then
       call pmc_mpi_abort(1)
    end if
#endif

  end subroutine pmc_mpi_check_ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_init()

    ! Initialize MPI.

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_init(ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_abort(status)

    ! Abort the program.

    integer, intent(in) :: status       ! Status flag to abort with

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_abort(MPI_COMM_WORLD, status, ierr)
#else
    call exit(status)
#endif

  end subroutine pmc_mpi_abort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_finalize()

    ! Shut down MPI.

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_finalize(ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_rank()

    ! Returns the rank of the current process.

#ifdef PMC_USE_MPI
    use mpi

    integer :: rank, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_rank = rank
#else
    pmc_mpi_rank = 0
#endif

  end function pmc_mpi_rank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_size()

    ! Returns the total number of processes.

#ifdef PMC_USE_MPI
    use mpi

    integer :: size, ierr

    call mpi_comm_size(MPI_COMM_WORLD, size, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_size = size
#else
    pmc_mpi_size = 1
#endif

  end function pmc_mpi_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_broadcast_integer(val)

    ! Broadcast the given value from process 0 to all other processes.

    integer, intent(inout) :: val ! value to broadcast

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, 1, MPI_INTEGER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_broadcast_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_broadcast_string(val)

    ! Broadcast the given value from process 0 to all other processes.

    character(len=*), intent(inout) :: val ! value to broadcast

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, len(val), MPI_CHARACTER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_broadcast_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_broadcast_char_array(val)

    ! Broadcast the given value from process 0 to all other processes.

    character, intent(inout) :: val(:) ! value to broadcast

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, size(val), MPI_CHARACTER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_broadcast_char_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_integer_size(val)

    ! Determines the number of bytes required to pack the given value.

    integer, intent(in) :: val          ! value to pack

    pmc_mpi_pack_integer_size = 4

  end function pmc_mpi_pack_integer_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_real_size(val)

    ! Determines the number of bytes required to pack the given value.

    real*8, intent(in) :: val           ! value to pack

    pmc_mpi_pack_real_size = 8

  end function pmc_mpi_pack_real_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_string_size(val)

    ! Determines the number of bytes required to pack the given value.

    character(len=*), intent(in) :: val ! value to pack

    pmc_mpi_pack_string_size = len_trim(val) + pmc_mpi_pack_integer_size(1)

  end function pmc_mpi_pack_string_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_logical_size(val)

    ! Determines the number of bytes required to pack the given value.

    use pmc_util

    logical, intent(in) :: val          ! value to pack

    integer, save :: saved_size = 0
    character :: buffer(20)

#ifdef PMC_USE_MPI
    integer :: ierr

    if (saved_size == 0) then
       ! first time through, compute the size
       call mpi_pack(val, 1, MPI_LOGICAL, buffer, size(buffer), &
            saved_size, MPI_COMM_WORLD, ierr)
       call pmc_mpi_check_ierr(ierr)
       call assert(saved_size < size(buffer))
    end if
#endif

    pmc_mpi_pack_logical_size = saved_size

  end function pmc_mpi_pack_logical_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_integer_array_size(val)

    ! Determines the number of bytes required to pack the given value.

    integer, intent(in) :: val(:)       ! value to pack

    pmc_mpi_pack_integer_array_size = pmc_mpi_pack_integer_size(1) &
         + size(val) * pmc_mpi_pack_integer_size(1)

  end function pmc_mpi_pack_integer_array_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_real_array_size(val)

    ! Determines the number of bytes required to pack the given value.

    real*8, intent(in) :: val(:)        ! value to pack

    pmc_mpi_pack_real_array_size = pmc_mpi_pack_integer_size(1) &
         + size(val) * pmc_mpi_pack_real_size(1d0)

  end function pmc_mpi_pack_real_array_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_string_array_size(val)

    ! Determines the number of bytes required to pack the given value.

    character(len=*), intent(in) :: val(:) ! value to pack

    integer :: i, total_size

    total_size = pmc_mpi_pack_integer_size(1)
    do i = 1,size(val)
       total_size = total_size + pmc_mpi_pack_string_size(val(i))
    end do
    pmc_mpi_pack_string_array_size = total_size

  end function pmc_mpi_pack_string_array_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_real_array_2d_size(val)

    ! Determines the number of bytes required to pack the given value.

    real*8, intent(in) :: val(:,:)      ! value to pack

    pmc_mpi_pack_real_array_2d_size = 2 * pmc_mpi_pack_integer_size(1) &
         + size(val) * pmc_mpi_pack_real_size(1d0)

  end function pmc_mpi_pack_real_array_2d_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_integer(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    integer, intent(in) :: val          ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_INTEGER, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_integer_size(val))
#endif

  end subroutine pmc_mpi_pack_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_real(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    real*8, intent(in) :: val           ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_REAL8, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_real_size(val))
#endif

  end subroutine pmc_mpi_pack_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_string(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    character(len=*), intent(in) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, length, ierr

    prev_position = position
    length = len_trim(val)
    call pmc_mpi_pack_integer(buffer, position, length)
    call mpi_pack(val, length, MPI_CHARACTER, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_string_size(val))
#endif

  end subroutine pmc_mpi_pack_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_logical(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    logical, intent(in) :: val           ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_LOGICAL, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_logical_size(val))
#endif

  end subroutine pmc_mpi_pack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_integer_array(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    integer, intent(in) :: val(:)       ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    n = size(val)
    call pmc_mpi_pack_integer(buffer, position, n)
    call mpi_pack(val, n, MPI_INTEGER, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position &
         == pmc_mpi_pack_integer_array_size(val))
#endif

  end subroutine pmc_mpi_pack_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_real_array(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    real*8, intent(in) :: val(:)        ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    n = size(val)
    call pmc_mpi_pack_integer(buffer, position, n)
    call mpi_pack(val, n, MPI_REAL8, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_real_array_size(val))
#endif

  end subroutine pmc_mpi_pack_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_string_array(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    character(len=*), intent(in) :: val(:) ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    n = size(val)
    call pmc_mpi_pack_integer(buffer, position, n)
    do i = 1,n
       call pmc_mpi_pack_string(buffer, position, val(i))
    end do
    call assert(position - prev_position == pmc_mpi_pack_string_array_size(val))
#endif

  end subroutine pmc_mpi_pack_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_real_array_2d(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    real*8, intent(in) :: val(:,:)      ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, n1, n2, ierr

    prev_position = position
    n1 = size(val, 1)
    n2 = size(val, 2)
    call pmc_mpi_pack_integer(buffer, position, n1)
    call pmc_mpi_pack_integer(buffer, position, n2)
    call mpi_pack(val, n1*n2, MPI_REAL8, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position &
         == pmc_mpi_pack_real_array_2d_size(val))
#endif

  end subroutine pmc_mpi_pack_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_integer(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    integer, intent(out) :: val         ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_INTEGER, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_integer_size(val))
#endif

  end subroutine pmc_mpi_unpack_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_real(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    real*8, intent(out) :: val          ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_REAL8, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_real_size(val))
#endif

  end subroutine pmc_mpi_unpack_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_string(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    character(len=*), intent(out) :: val ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, length, ierr

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, length)
    if (length > len(val)) then
       call pmc_mpi_abort(3251)
    end if
    val = ''
    call mpi_unpack(buffer, size(buffer), position, val, length, &
         MPI_CHARACTER, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_string_size(val))
#endif

  end subroutine pmc_mpi_unpack_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_logical(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    logical, intent(out) :: val          ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_LOGICAL, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_logical_size(val))
#endif

  end subroutine pmc_mpi_unpack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_integer_array(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    integer, pointer :: val(:)          ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    allocate(val(n))
    call mpi_unpack(buffer, size(buffer), position, val, n, MPI_INTEGER, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position &
         == pmc_mpi_pack_integer_array_size(val))
#endif

  end subroutine pmc_mpi_unpack_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_real_array(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    real*8, pointer :: val(:)           ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    allocate(val(n))
    call mpi_unpack(buffer, size(buffer), position, val, n, MPI_REAL8, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position == pmc_mpi_pack_real_array_size(val))
#endif

  end subroutine pmc_mpi_unpack_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_string_array(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    character(len=*), pointer :: val(:) ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    allocate(val(n))
    do i = 1,n
       call pmc_mpi_unpack_string(buffer, position, val(i))
    end do
    call assert(position - prev_position == pmc_mpi_pack_string_array_size(val))
#endif

  end subroutine pmc_mpi_unpack_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_real_array_2d(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    real*8, pointer :: val(:,:)         ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, n1, n2, ierr

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n1)
    call pmc_mpi_unpack_integer(buffer, position, n2)
    allocate(val(n1,n2))
    call mpi_unpack(buffer, size(buffer), position, val, n1*n2, MPI_REAL8, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(position - prev_position &
         == pmc_mpi_pack_real_array_2d_size(val))
#endif

  end subroutine pmc_mpi_unpack_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_average_real(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    real*8, intent(in) :: val           ! value to average
    real*8, intent(out) :: val_avg      ! result

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_reduce(val, val_avg, 1, MPI_REAL8, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    if (pmc_mpi_rank() == 0) then
       val_avg = val_avg / dble(pmc_mpi_size())
    end if
#else
    val_avg = val
#endif

  end subroutine pmc_mpi_reduce_average_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_average_real_array(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_util

    real*8, intent(in) :: val(:)        ! value to average
    real*8, intent(out) :: val_avg(:)   ! result

#ifdef PMC_USE_MPI
    integer :: ierr

    call assert(size(val) == size(val_avg))
    call mpi_reduce(val, val_avg, size(val), MPI_REAL8, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    if (pmc_mpi_rank() == 0) then
       val_avg = val_avg / dble(pmc_mpi_size())
    end if
#else
    val_avg = val
#endif

  end subroutine pmc_mpi_reduce_average_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_average_real_array_2d(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_util

    real*8, intent(in) :: val(:,:)      ! value to average
    real*8, intent(out) :: val_avg(:,:) ! result

#ifdef PMC_USE_MPI
    integer :: ierr

    call assert(size(val,1) == size(val_avg,1))
    call assert(size(val,2) == size(val_avg,2))
    call mpi_reduce(val, val_avg, size(val), MPI_REAL8, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    if (pmc_mpi_rank() == 0) then
       val_avg = val_avg / dble(pmc_mpi_size())
    end if
#else
    val_avg = val
#endif

  end subroutine pmc_mpi_reduce_average_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mpi
