! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Read NetCDF state files and write out the environment variables in
! text format.

program extract_env

  use netcdf

  integer, parameter :: out_unit = 64

  character(len=1000) :: in_prefix, in_filename, out_filename
  integer :: ncid
  integer :: varid_time, varid_temp, varid_rh
  integer :: varid_pres, varid_height
  real*8 :: time, temp, rh, pres, height
  integer :: ios, i_time, status, n_time

  ! process commandline arguments
  if (iargc() .ne. 2) then
     write(6,*) 'Usage: extract_env <netcdf_state_prefix> <output_filename>'
     call exit(2)
  endif
  call getarg(1, in_prefix)
  call getarg(2, out_filename)

  ! open output file
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if

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
     call exit(1)
  end if

  close(out_unit)

contains

  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

end program extract_env
