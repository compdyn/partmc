! Copyright (C) 2007-2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_mpi module.

!> Wrapper functions for MPI.
!!
!! All of these functions can be called irrespective of whether MPI
!! support was compiled in or not. If MPI support is not enabled then
!! they do the obvious trivial thing (normally nothing).
module pmc_mpi

  use pmc_util
  
#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Whether MPI support is compiled in.
  logical function pmc_mpi_support()

#ifdef PMC_USE_MPI
    pmc_mpi_support = .true.
#else
    pmc_mpi_support = .false.
#endif

  end function pmc_mpi_support

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dies if \c ierr is not ok.
  subroutine pmc_mpi_check_ierr(ierr)

    !> MPI status code.
    integer, intent(in) :: ierr

#ifdef PMC_USE_MPI
    if (ierr /= MPI_SUCCESS) then
       call pmc_mpi_abort(1)
    end if
#endif

  end subroutine pmc_mpi_check_ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize MPI.
  subroutine pmc_mpi_init()

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_init(ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Abort the program.
  subroutine pmc_mpi_abort(status)

    !> Status flag to abort with.
    integer, intent(in) :: status

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_abort(MPI_COMM_WORLD, status, ierr)
#else
    call exit(status)
#endif

  end subroutine pmc_mpi_abort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Shut down MPI.
  subroutine pmc_mpi_finalize()

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_finalize(ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Synchronize all processors.
  subroutine pmc_mpi_barrier()

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_barrier

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the rank of the current process.
  integer function pmc_mpi_rank()

#ifdef PMC_USE_MPI
    integer :: rank, ierr

    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_rank = rank
#else
    pmc_mpi_rank = 0
#endif

  end function pmc_mpi_rank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the total number of processes.
  integer function pmc_mpi_size()

#ifdef PMC_USE_MPI
    integer :: size, ierr

    call mpi_comm_size(MPI_COMM_WORLD, size, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_size = size
#else
    pmc_mpi_size = 1
#endif

  end function pmc_mpi_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Broadcast the given value from process 0 to all other processes.
  subroutine pmc_mpi_bcast_integer(val)

    !> Value to broadcast.
    integer, intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, 1, MPI_INTEGER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_bcast_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Broadcast the given value from process 0 to all other processes.
  subroutine pmc_mpi_bcast_string(val)

    !> Value to broadcast.
    character(len=*), intent(inout) :: val

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, len(val), MPI_CHARACTER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_bcast_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Broadcast the given value from process 0 to all other processes.
  subroutine pmc_mpi_bcast_packed(val)

    !> Value to broadcast.
    character, intent(inout) :: val(:)

#ifdef PMC_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast(val, size(val), MPI_CHARACTER, root, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#endif

  end subroutine pmc_mpi_bcast_packed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_integer(val)

    !> Value to pack.
    integer, intent(in) :: val

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_integer, ierr)
    call pmc_mpi_check_ierr(ierr)
!DEBUG
    pmc_mpi_pack_size_integer = 4
!DEBUG
#else
    pmc_mpi_pack_size_integer = 0
#endif

  end function pmc_mpi_pack_size_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_real(val)

    !> Value to pack.
    real*8, intent(in) :: val

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(1, MPI_REAL8, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_real, ierr)
    call pmc_mpi_check_ierr(ierr)
!DEBUG
    pmc_mpi_pack_size_real = 8
!DEBUG
#else
    pmc_mpi_pack_size_real = 0
#endif

  end function pmc_mpi_pack_size_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_string(val)

    !> Value to pack.
    character(len=*), intent(in) :: val

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(len_trim(val), MPI_CHARACTER, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_string, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_string = pmc_mpi_pack_size_string &
         + pmc_mpi_pack_size_integer(len_trim(val))
!DEBUG
    pmc_mpi_pack_size_string = len_trim(val) &
         + pmc_mpi_pack_size_integer(len_trim(val))
!DEBUG
#else
    pmc_mpi_pack_size_string = 0
#endif

  end function pmc_mpi_pack_size_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_logical(val)

    !> Value to pack.
    logical, intent(in) :: val

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(1, MPI_LOGICAL, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_logical, ierr)
    call pmc_mpi_check_ierr(ierr)
!DEBUG
    pmc_mpi_pack_size_logical = 4
!DEBUG
#else
    pmc_mpi_pack_size_logical = 0
#endif

  end function pmc_mpi_pack_size_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_complex(val)

    !> Value to pack.
    complex*16, intent(in) :: val

    pmc_mpi_pack_size_complex = 16

  end function pmc_mpi_pack_size_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_integer_array(val)

    !> Value to pack.
    integer, intent(in) :: val(:)

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(size(val), MPI_INTEGER, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_integer_array, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_integer_array = pmc_mpi_pack_size_integer_array &
         + pmc_mpi_pack_size_integer(size(val))
!DEBUG
    pmc_mpi_pack_size_integer_array = size(val) * 4 &
         + pmc_mpi_pack_size_integer(size(val))
!DEBUG
#else
    pmc_mpi_pack_size_integer_array = 0
#endif

  end function pmc_mpi_pack_size_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_real_array(val)

    !> Value to pack.
    real*8, intent(in) :: val(:)

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(size(val), MPI_REAL8, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_real_array, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_real_array = pmc_mpi_pack_size_real_array &
         + pmc_mpi_pack_size_integer(size(val))
!DEBUG
    pmc_mpi_pack_size_real_array = size(val) * 8 &
         + pmc_mpi_pack_size_integer(size(val))
!DEBUG
#else
    pmc_mpi_pack_size_real_array = 0
#endif

  end function pmc_mpi_pack_size_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_string_array(val)

    !> Value to pack.
    character(len=*), intent(in) :: val(:)

    integer :: i, total_size

    total_size = pmc_mpi_pack_size_integer(size(val))
    do i = 1,size(val)
       total_size = total_size + pmc_mpi_pack_size_string(val(i))
    end do
    pmc_mpi_pack_size_string_array = total_size

  end function pmc_mpi_pack_size_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_real_array_2d(val)

    !> Value to pack.
    real*8, intent(in) :: val(:,:)

    integer :: ierr

#ifdef PMC_USE_MPI
    call mpi_pack_size(size(val), MPI_REAL8, MPI_COMM_WORLD, &
         pmc_mpi_pack_size_real_array_2d, ierr)
    call pmc_mpi_check_ierr(ierr)
    pmc_mpi_pack_size_real_array_2d = pmc_mpi_pack_size_real_array_2d &
         + pmc_mpi_pack_size_integer(size(val,1)) &
         + pmc_mpi_pack_size_integer(size(val,2))
!DEBUG
    pmc_mpi_pack_size_real_array_2d = size(val) * 8 &
         + pmc_mpi_pack_size_integer(size(val,1)) &
         + pmc_mpi_pack_size_integer(size(val,2))
!DEBUG
#else
    pmc_mpi_pack_size_real_array_2d = 0
#endif

  end function pmc_mpi_pack_size_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_integer(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    integer, intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_INTEGER, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(913495993, &
         position - prev_position == pmc_mpi_pack_size_integer(val))
#endif

  end subroutine pmc_mpi_pack_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_real(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    real*8, intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_REAL8, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(395354132, &
         position - prev_position == pmc_mpi_pack_size_real(val))
#endif

  end subroutine pmc_mpi_pack_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_string(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    character(len=*), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, length, ierr

    prev_position = position
    length = len_trim(val)
    call pmc_mpi_pack_integer(buffer, position, length)
    call mpi_pack(val, length, MPI_CHARACTER, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(607212018, &
         position - prev_position == pmc_mpi_pack_size_string(val))
#endif

  end subroutine pmc_mpi_pack_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_logical(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    logical, intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_LOGICAL, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(104535200, &
         position - prev_position == pmc_mpi_pack_size_logical(val))
#endif

  end subroutine pmc_mpi_pack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_complex(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    complex*16, intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_pack(val, 1, MPI_COMPLEX16, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(640416372, &
         position - prev_position == pmc_mpi_pack_size_complex(val))
#endif

  end subroutine pmc_mpi_pack_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_integer_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    integer, intent(in) :: val(:)

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    n = size(val)
    call pmc_mpi_pack_integer(buffer, position, n)
    call mpi_pack(val, n, MPI_INTEGER, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(698601296, &
         position - prev_position == pmc_mpi_pack_size_integer_array(val))
#endif

  end subroutine pmc_mpi_pack_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_real_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    real*8, intent(in) :: val(:)

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    n = size(val)
    call pmc_mpi_pack_integer(buffer, position, n)
    call mpi_pack(val, n, MPI_REAL8, buffer, size(buffer), &
         position, MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(825718791, &
         position - prev_position == pmc_mpi_pack_size_real_array(val))
#endif

  end subroutine pmc_mpi_pack_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_string_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    character(len=*), intent(in) :: val(:)

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    n = size(val)
    call pmc_mpi_pack_integer(buffer, position, n)
    do i = 1,n
       call pmc_mpi_pack_string(buffer, position, val(i))
    end do
    call assert(630900704, &
         position - prev_position == pmc_mpi_pack_size_string_array(val))
#endif

  end subroutine pmc_mpi_pack_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_real_array_2d(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    real*8, intent(in) :: val(:,:)

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
    call assert(567349745, &
         position - prev_position == pmc_mpi_pack_size_real_array_2d(val))
#endif

  end subroutine pmc_mpi_pack_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_integer(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    integer, intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_INTEGER, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(890243339, &
         position - prev_position == pmc_mpi_pack_size_integer(val))
#endif

  end subroutine pmc_mpi_unpack_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_real(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    real*8, intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_REAL8, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(570771632, &
         position - prev_position == pmc_mpi_pack_size_real(val))
#endif

  end subroutine pmc_mpi_unpack_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_string(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    character(len=*), intent(out) :: val

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
    call assert(503378058, &
         position - prev_position == pmc_mpi_pack_size_string(val))
#endif

  end subroutine pmc_mpi_unpack_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_logical(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    logical, intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_LOGICAL, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(694750528, &
         position - prev_position == pmc_mpi_pack_size_logical(val))
#endif

  end subroutine pmc_mpi_unpack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_complex(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    complex*16, intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position, ierr

    prev_position = position
    call mpi_unpack(buffer, size(buffer), position, val, 1, MPI_COMPLEX16, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(969672634, &
         position - prev_position == pmc_mpi_pack_size_complex(val))
#endif

  end subroutine pmc_mpi_unpack_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_integer_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    integer, pointer :: val(:)

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    deallocate(val)
    allocate(val(n))
    call mpi_unpack(buffer, size(buffer), position, val, n, MPI_INTEGER, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(565840919, &
         position - prev_position == pmc_mpi_pack_size_integer_array(val))
#endif

  end subroutine pmc_mpi_unpack_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_real_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    real*8, pointer :: val(:)

#ifdef PMC_USE_MPI
    integer :: prev_position, n, ierr

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    deallocate(val)
    allocate(val(n))
    call mpi_unpack(buffer, size(buffer), position, val, n, MPI_REAL8, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(782875761, &
         position - prev_position == pmc_mpi_pack_size_real_array(val))
#endif

  end subroutine pmc_mpi_unpack_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_string_array(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    character(len=*), pointer :: val(:)

#ifdef PMC_USE_MPI
    integer :: prev_position, i, n

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n)
    deallocate(val)
    allocate(val(n))
    do i = 1,n
       call pmc_mpi_unpack_string(buffer, position, val(i))
    end do
    call assert(320065648, &
         position - prev_position == pmc_mpi_pack_size_string_array(val))
#endif

  end subroutine pmc_mpi_unpack_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_real_array_2d(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    real*8, pointer :: val(:,:)

#ifdef PMC_USE_MPI
    integer :: prev_position, n1, n2, ierr

    prev_position = position
    call pmc_mpi_unpack_integer(buffer, position, n1)
    call pmc_mpi_unpack_integer(buffer, position, n2)
    deallocate(val)
    allocate(val(n1,n2))
    call mpi_unpack(buffer, size(buffer), position, val, n1*n2, MPI_REAL8, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
    call assert(781681739, position - prev_position &
         == pmc_mpi_pack_size_real_array_2d(val))
#endif

  end subroutine pmc_mpi_unpack_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_real(val, val_avg)

    !> Value to average.
    real*8, intent(in) :: val
    !> Result.
    real*8, intent(out) :: val_avg

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the sum of \c val across all processes, storing the
  !> result in \c val_sum on the root process.
  subroutine pmc_mpi_reduce_sum_integer(val, val_sum)

    !> Value to sum.
    integer, intent(in) :: val
    !> Result.
    integer, intent(out) :: val_sum

#ifdef PMC_USE_MPI
    integer :: ierr

    call mpi_reduce(val, val_sum, 1, MPI_INTEGER, MPI_SUM, 0, &
         MPI_COMM_WORLD, ierr)
    call pmc_mpi_check_ierr(ierr)
#else
    val_sum = val
#endif

  end subroutine pmc_mpi_reduce_sum_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_real_array(val, val_avg)

    !> Value to average.
    real*8, intent(in) :: val(:)
    !> Result.
    real*8, intent(out) :: val_avg(:)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_reduce_avg_real_array_2d(val, val_avg)

    !> Value to average.
    real*8, intent(in) :: val(:,:)
    !> Result.
    real*8, intent(out) :: val_avg(:,:)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_allreduce_average_real(val, val_avg)

    !> Value to average.
    real*8, intent(in) :: val
    !> Result.
    real*8, intent(out) :: val_avg

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the average of val across all processes, storing the
  !> result in val_avg on the root process.
  subroutine pmc_mpi_allreduce_average_real_array(val, val_avg)

    !> Value to average.
    real*8, intent(in) :: val(:)
    !> Result.
    real*8, intent(out) :: val_avg(:)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_mpi
