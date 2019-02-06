! Copyright (C) 2019 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_solver_stats module

!> The solver_stats_t type and associated subroutines
module pmc_solver_stats

  use pmc_constants,                  only : i_kind, dp

  implicit none
  private

  public :: solver_stats_t

  !> Solver statistics
  !!
  !! Holds information related to a solver run
  type :: solver_stats_t
    !> Status code
    integer(kind=i_kind) :: status_code
    !> Integration start time [s]
    real(kind=dp) :: start_time__s
    !> Integration end time [s]
    real(kind=dp) :: end_time__s
    !> Number of steps
    integer(kind=i_kind) :: num_steps
    !> Right-hand side evaluations
    integer(kind=i_kind) :: RHS_evals
    !> Linear solver setups
    integer(kind=i_kind) :: LS_setups
    !> Error test failures
    integer(kind=i_kind) :: error_test_fails
    !> Non-Linear solver iterations
    integer(kind=i_kind) :: NLS_iters
    !> Non-Linear solver convergence failures
    integer(kind=i_kind) :: NLS_convergence_fails
    !> Direct Linear Solver Jacobian evaluations
    integer(kind=i_kind) :: DLS_Jac_evals
    !> Direct Linear Solver right-hand size evaluations
    integer(kind=i_kind) :: DLS_RHS_evals
    !> Last time step [s]
    real(kind=dp) :: last_time_step__s
    !> Next time step [s]
    real(kind=dp) :: next_time_step__s
  contains
    !> Print the solver statistics
    procedure :: print => do_print
    !> Assignment
    procedure :: assignValue
    generic :: assignment(=) => assignValue
  end type solver_stats_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Print the solver statistics
  subroutine do_print( this, file_unit )

    !> Solver statistics
    class(solver_stats_t), intent(in) :: this
    !> File unit to output to
    integer(kind=i_kind), optional :: file_unit

    integer(kind=i_kind) :: f_unit = 6

    if( present( file_unit ) ) f_unit = file_unit

    write(f_unit,*) "Status code:                 ", this%status_code
    write(f_unit,*) "Integration start time [s]:  ", this%start_time__s
    write(f_unit,*) "Integration end time [s]:    ", this%end_time__s
    write(f_unit,*) "Number of steps:             ", this%num_steps
    write(f_unit,*) "Right-hand side evals:       ", this%RHS_evals
    write(f_unit,*) "Linear solver setups:        ", this%LS_setups
    write(f_unit,*) "Error test failures:         ", this%error_test_fails
    write(f_unit,*) "Non-Linear solver iterations:", this%NLS_iters
    write(f_unit,*) "Non-Linear convergence fails:", this%NLS_convergence_fails
    write(f_unit,*) "DLS Jacobian evals:          ", this%DLS_Jac_evals
    write(f_unit,*) "DLS Right-hand side evals:   ", this%DLS_RHS_evals
    write(f_unit,*) "Last time step [s]:          ", this%last_time_step__s
    write(f_unit,*) "Next time step [s]:          ", this%next_time_step__s

  end subroutine do_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Assign a value to all members of solver stats
  subroutine assignValue( this, new_value )

    !> Solver statistics
    class(solver_stats_t), intent(out) :: this
    !> Value to assign
    integer(kind=i_kind), intent(in) :: new_value

    this%status_code           = new_value
    this%start_time__s         = real( new_value, kind=dp )
    this%end_time__s           = real( new_value, kind=dp )
    this%num_steps             = new_value
    this%RHS_evals             = new_value
    this%LS_setups             = new_value
    this%error_test_fails      = new_value
    this%NLS_iters             = new_value
    this%NLS_convergence_fails = new_value
    this%DLS_Jac_evals         = new_value
    this%DLS_RHS_evals         = new_value
    this%last_time_step__s     = real( new_value, kind=dp )
    this%next_time_step__s     = real( new_value, kind=dp )

  end subroutine assignValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_solver_stats
