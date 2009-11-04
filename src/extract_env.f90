! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_env program.

!> Read NetCDF output files and write out the environment variables in
!> text format.
program extract_env

  use netcdf

  integer, parameter :: dp = kind(0.d0)
  integer, parameter :: out_unit = 64

  character(len=1000) :: in_prefix, in_filename, out_filename
  integer :: ncid
  integer :: varid_time, varid_temp, varid_rh
  integer :: varid_pres, varid_height
  real(kind=dp) :: time, temp, rh, pres, height
  integer :: ios, i_time, status, n_time

  ! process commandline arguments
  if (command_argument_count() .ne. 2) then
     write(6,*) 'Usage: extract_env <netcdf_state_prefix> <output_filename>'
     stop 2
  endif
  call get_command_argument(1, in_prefix)
  call get_command_argument(2, out_filename)

  ! open output file
  open(unit=out_unit, file=out_filename, status='replace', iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     stop 1
  end if

  ! write information
  write(*,'(a,a)') "Output file: ", trim(out_filename)
  write(*,'(a)') "  Each row of output is one time."
  write(*,'(a)') "  The columns of output are:"
  write(*,'(a)') "    column 1: time (s)"
  write(*,'(a)') "    column 2: temperature (K)"
  write(*,'(a)') "    column 3: relative_humidity (1)"
  write(*,'(a)') "    column 4: pressure (Pa)"
  write(*,'(a)') "    column 5: mixing height (m)"

  ! read NetCDF files
  i_time = 0
  n_time = 0
  do while (.true.)
     i_time = i_time + 1
     write(in_filename,'(a,i8.8,a)') trim(in_prefix), i_time, ".nc"
     status = nf90_open(in_filename, NF90_NOWRITE, ncid)
     if (status /= NF90_NOERR) then
        exit
     end if
     n_time = i_time

     call nc_check(nf90_inq_varid(ncid, "time", varid_time))
     call nc_check(nf90_get_var(ncid, varid_time, time))

     call nc_check(nf90_inq_varid(ncid, "temperature", varid_temp))
     call nc_check(nf90_get_var(ncid, varid_temp, temp))

     call nc_check(nf90_inq_varid(ncid, "relative_humidity", varid_rh))
     call nc_check(nf90_get_var(ncid, varid_rh, rh))

     call nc_check(nf90_inq_varid(ncid, "pressure", varid_pres))
     call nc_check(nf90_get_var(ncid, varid_pres, pres))

     call nc_check(nf90_inq_varid(ncid, "height", varid_height))
     call nc_check(nf90_get_var(ncid, varid_height, height))

     call nc_check(nf90_close(ncid))

     ! output data
     write(out_unit, '(5e30.15e3)') time, temp, rh, pres, height
  end do

  if (n_time == 0) then
     write(*,'(a,a)') 'ERROR: no input file found matching: ', &
          trim(in_filename)
     stop 1
  end if

  close(out_unit)

contains

  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       stop 1
    end if

  end subroutine nc_check

#ifdef DEFINE_LOCAL_COMMAND_ARGUMENT_COUNT
  integer function command_argument_count()
    command_argument_count = iargc()
  end function command_argument_count
#endif
#ifdef DEFINE_LOCAL_GET_COMMAND_ARGUMENT
  subroutine get_command_argument(i, arg)
    integer, intent(in) :: i
    character(len=*), intent(out) :: arg
    call getarg(i, arg)
  end subroutine get_command_argument
#endif
  
end program extract_env
