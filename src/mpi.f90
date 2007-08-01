! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Wrapper functions for MPI. All of these functions can be called
! irrespective of whether MPI support was compiled in or not. If MPI
! support is not enabled then they do the obvious trivial thing
! (normally nothing).

module pmc_mpi
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_check_ierr(ierr)

    ! Dies if ierr is not ok.

    integer, intent(in) :: ierr         ! MPI status code

#ifdef PMC_USE_MPI
    if (ierr /= MPI_SUCCESS) then
       mpi_abort(MPI_COMM_WORLD, 1)
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
    call mpi_abort(MPI_COMM_WORLD, status)
#else
    call exit(status)
#endif

  end subroutine pmc_mpi_abort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_mpi_finalize()

    ! Shut down MPI.

#ifdef PMC_USE_MPI
    call mpi_finalize()
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

end module pmc_mpi
