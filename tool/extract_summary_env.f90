! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Read a NetCDF summary file and write out the gas concenrations in
! text format.

program extract_summary_env

  use netcdf

  integer, parameter :: out_unit = 64

  character(len=1000) :: in_filename, out_filename
  integer :: ncid
  integer :: dimid_time, dimid_env
  integer :: varid_time, varid_env
  integer :: varid_env_state
  integer :: n_time, n_env
  character(len=1000) :: tmp_str, env_names
  real*8, allocatable :: time(:)
  real*8, allocatable :: env_state(:,:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_env

  ! process commandline arguments
  if (iargc() .ne. 2) then
     write(6,*) 'Usage: extract_summary_env <netcdf_filename> <output_filename>'
     call exit(2)
  endif
  call getarg(1, in_filename)
  call getarg(2, out_filename)

  ! read NetCDF file
  call nc_check(nf90_open(in_filename, NF90_NOWRITE, ncid))

  call nc_check(nf90_inq_dimid(ncid, "time", dimid_time))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_time, &
       tmp_str, n_time))
  allocate(time(n_time))
  call nc_check(nf90_inq_varid(ncid, "time", varid_time))
  call nc_check(nf90_get_var(ncid, varid_time, time))
  write(*,*) "n_time:", n_time
  write(*,*) "min time:", minval(time)
  write(*,*) "max time:", maxval(time)

  call nc_check(nf90_inq_dimid(ncid, "env", dimid_env))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_env, &
       tmp_str, n_env))
  call nc_check(nf90_inq_varid(ncid, "env", varid_env))
  call nc_check(nf90_get_att(ncid, varid_env, &
       "names", env_names))
  write(*,*) "n_env:", n_env
  write(*,*) "env_names: ", trim(env_names)
  if (n_env /= 4) then
     write(0,*) 'ERROR: invalid size of env dimension'
     call exit(1)
  end if

  call nc_check(nf90_inq_varid(ncid, "env_state", varid_env_state))
  call nc_check(nf90_Inquire_Variable(ncid, varid_env_state, tmp_str, &
       xtype, ndims, dimids, nAtts))
  if ((ndims /= 2) &
       .or. (dimids(1) /= dimid_env) &
       .or. (dimids(2) /= dimid_time)) then
     write(*,*) "ERROR: unexpected env_state dimids"
     call exit(1)
  end if
  allocate(env_state(n_env, n_time))
  call nc_check(nf90_get_var(ncid, varid_env_state, env_state))

  call nc_check(nf90_close(ncid))

  ! output data
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if
  do i_time = 1,n_time
     write(out_unit, '(e30.15e3)', advance='no') time(i_time)
     do i_env = 1,n_env
        write(out_unit, '(e30.15e3)', advance='no') env_state(i_env, i_time)
     end do
     write(out_unit, '(a)') ''
  end do
  close(out_unit)

  deallocate(time)
  deallocate(env_state)

contains

  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

end program extract_summary_env
