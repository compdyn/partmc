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

  logical function pmc_mpi_support()

    ! Whether MPI support is compiled in.

#ifdef PMC_USE_MPI
    pmc_mpi_support = .true.
#else
    pmc_mpi_support = .false.
#endif

  end function pmc_mpi_support

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

  subroutine pmc_mpi_barrier()

    ! Synchronize all processors.

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_barrier

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

  subroutine pmc_mpi_bcast_integer(val)

    ! Broadcast the given value from process 0 to all other processes.

    integer, intent(inout) :: val ! value to broadcast

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, 1, MPI_INTEGER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_bcast_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_bcast_string(val)

    ! Broadcast the given value from process 0 to all other processes.

    character(len=*), intent(inout) :: val ! value to broadcast

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, len(val), MPI_CHARACTER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_bcast_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_bcast_packed(val)

    ! Broadcast the given value from process 0 to all other processes.

    character, intent(inout) :: val(:) ! value to broadcast

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, size(val), MPI_CHARACTER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_bcast_packed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_integer(val)

    ! Determines the number of bytes required to pack the given value.

    integer, intent(in) :: val          ! value to pack

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_integer, ierr)
    call pmc_mpi_check_ierr(ierr)
#else
    pmc_mpi_pack_size_integer = 0
#endif

  end function pmc_mpi_pack_size_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_real(val)

    ! Determines the number of bytes required to pack the given value.

    real*8, intent(in) :: val           ! value to pack

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(1, MPI_REAL8, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_real, ierr)
    call pmc_mpi_check_ierr(ierr)
#else
    pmc_mpi_pack_size_real = 0
#endif

  end function pmc_mpi_pack_size_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_string(val)

    ! Determines the number of bytes required to pack the given value.

    character(len=*), intent(in) :: val ! value to pack

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(len_trim(val), MPI_CHARACTER, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_string, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_string = pmc_mpi_pack_size_string &
         + pmc_mpi_pack_size_integer(len_trim(val))
#else
    pmc_mpi_pack_size_string = 0
#endif

  end function pmc_mpi_pack_size_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_logical(val)

    ! Determines the number of bytes required to pack the given value.

    logical, intent(in) :: val          ! value to pack

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(1, MPI_LOGICAL, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_logical, ierr)
    call pmc_mpi_check_ierr(ierr)
#else
    pmc_mpi_pack_size_logical = 0
#endif

  end function pmc_mpi_pack_size_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_complex(val)

    ! Determines the number of bytes required to pack the given value.

    complex*16, intent(in) :: val       ! value to pack

    pmc_mpi_pack_size_complex = 16

  end function pmc_mpi_pack_size_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_integer_array(val)

    ! Determines the number of bytes required to pack the given value.

    integer, intent(in) :: val(:)       ! value to pack

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(size(val), MPI_INTEGER, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_integer_array, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_integer_array = pmc_mpi_pack_size_integer_array &
         + pmc_mpi_pack_size_integer(size(val))
#else
    pmc_mpi_pack_size_integer_array = 0
#endif

  end function pmc_mpi_pack_size_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_real_array(val)

    ! Determines the number of bytes required to pack the given value.

    real*8, intent(in) :: val(:)        ! value to pack

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(size(val), MPI_REAL8, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_real_array, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_real_array = pmc_mpi_pack_size_real_array &
         + pmc_mpi_pack_size_integer(size(val))
#else
    pmc_mpi_pack_size_real_array = 0
#endif

  end function pmc_mpi_pack_size_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_string_array(val)

    ! Determines the number of bytes required to pack the given value.

    character(len=*), intent(in) :: val(:) ! value to pack

    integer :: i, total_size

    total_size = pmc_mpi_pack_size_integer(size(val))
    do i = 1,size(val)
       total_size = total_size + pmc_mpi_pack_size_string(val(i))
    end do
    pmc_mpi_pack_size_string_array = total_size

  end function pmc_mpi_pack_size_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pmc_mpi_pack_size_real_array_2d(val)

    ! Determines the number of bytes required to pack the given value.

    real*8, intent(in) :: val(:,:)      ! value to pack

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(size(val), MPI_REAL8, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_real_array_2d, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_real_array_2d = pmc_mpi_pack_size_real_array_2d &
         + pmc_mpi_pack_size_integer(size(val,1)) &
         + pmc_mpi_pack_size_integer(size(val,2))
#else
    pmc_mpi_pack_size_real_array_2d = 0
#endif

  end function pmc_mpi_pack_size_real_array_2d

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
    call assert(913495993, position - prev_position == pmc_mpi_pack_size_integer(val))
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
    call assert(395354132, position - prev_position == pmc_mpi_pack_size_real(val))
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
    call assert(607212018, position - prev_position == pmc_mpi_pack_size_string(val))
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
    call assert(104535200, position - prev_position == pmc_mpi_pack_size_logical(val))
#endif

  end subroutine pmc_mpi_pack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_pack_complex(buffer, position, val)

    ! Packs the given value into the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    complex*16, intent(in) :: val       ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_COMPLEX16, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(640416372, position - prev_position == pmc_mpi_pack_size_complex(val))
#endif

  end subroutine pmc_mpi_pack_complex

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
    call assert(698601296, position - prev_position &
         == pmc_mpi_pack_size_integer_array(val))
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
    call assert(825718791, position - prev_position == pmc_mpi_pack_size_real_array(val))
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
    call assert(630900704, position - prev_position == pmc_mpi_pack_size_string_array(val))
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
    call assert(567349745, position - prev_position &
         == pmc_mpi_pack_size_real_array_2d(val))
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
    call assert(890243339, position - prev_position == pmc_mpi_pack_size_integer(val))
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
    call assert(570771632, position - prev_position == pmc_mpi_pack_size_real(val))
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
    call assert(503378058, position - prev_position == pmc_mpi_pack_size_string(val))
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
    call assert(694750528, position - prev_position == pmc_mpi_pack_size_logical(val))
#endif

  end subroutine pmc_mpi_unpack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_unpack_complex(buffer, position, val)

    ! Unpacks the given value from the buffer, advancing position.

    use pmc_util

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position  ! current buffer position
    complex*16, intent(out) :: val      ! value to pack

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_COMPLEX16, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(969672634, position - prev_position == pmc_mpi_pack_size_complex(val))
#endif

  end subroutine pmc_mpi_unpack_complex

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
    call assert(565840919, position - prev_position &
         == pmc_mpi_pack_size_integer_array(val))
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
    call assert(782875761, position - prev_position == pmc_mpi_pack_size_real_array(val))
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
    call assert(320065648, position - prev_position == pmc_mpi_pack_size_string_array(val))
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
    call assert(781681739, position - prev_position &
         == pmc_mpi_pack_size_real_array_2d(val))
#endif

  end subroutine pmc_mpi_unpack_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_real(val, val_avg)

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

  end subroutine pmc_mpi_reduce_avg_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_real_array(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_util

    real*8, intent(in) :: val(:)        ! value to average
    real*8, intent(out) :: val_avg(:)   ! result

#ifdef PMC_USE_MPI
    integer :: ierr

    call assert(915136121, size(val) == size(val_avg))
    call mpi_reduce(val, val_avg, size(val), MPI_REAL8, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    if (pmc_mpi_rank() == 0) then
       val_avg = val_avg / dble(pmc_mpi_size())
    end if
#else
    val_avg = val
#endif

  end subroutine pmc_mpi_reduce_avg_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_reduce_avg_real_array_2d(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_util

    real*8, intent(in) :: val(:,:)      ! value to average
    real*8, intent(out) :: val_avg(:,:) ! result

#ifdef PMC_USE_MPI
    integer :: ierr

    call assert(131229046, size(val,1) == size(val_avg,1))
    call assert(992122167, size(val,2) == size(val_avg,2))
    call mpi_reduce(val, val_avg, size(val), MPI_REAL8, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    if (pmc_mpi_rank() == 0) then
       val_avg = val_avg / dble(pmc_mpi_size())
    end if
#else
    val_avg = val
#endif

  end subroutine pmc_mpi_reduce_avg_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_allreduce_average_real(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    real*8, intent(in) :: val           ! value to average
    real*8, intent(out) :: val_avg      ! result

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_allreduce(val, val_avg, 1, MPI_REAL8, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    val_avg = val_avg / dble(pmc_mpi_size())
#else
    val_avg = val
#endif

  end subroutine pmc_mpi_allreduce_average_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_allreduce_average_real_array(val, val_avg)

    ! Computes the average of val across all processes, storing the
    ! result in val_avg on the root process.

    use pmc_util

    real*8, intent(in) :: val(:)        ! value to average
    real*8, intent(out) :: val_avg(:)   ! result

#ifdef PMC_USE_MPI
    integer :: ierr

    call assert(948533359, size(val) == size(val_avg))
    call mpi_allreduce(val, val_avg, size(val), MPI_REAL8, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    val_avg = val_avg / dble(pmc_mpi_size())
#else
    val_avg = val
#endif

  end subroutine pmc_mpi_allreduce_average_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mpi
