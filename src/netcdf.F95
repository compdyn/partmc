! Copyright (C) 2007-2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_netcdf module.

!> Wrapper functions for NetCDF.
module pmc_netcdf

  use netcdf
  use pmc_util

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the status of a NetCDF function call.
  subroutine pmc_nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       call die_msg(291021908, nf90_strerror(status))
    end if

  end subroutine pmc_nc_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the status of a NetCDF function call and prints the given
  !> error message on failure.
  subroutine pmc_nc_check_msg(status, error_msg)

    !> Status return value.
    integer, intent(in) :: status
    !> Error message in case of failure.
    character(len=*), intent(in) :: error_msg

    if (status /= NF90_NOERR) then
       call die_msg(701841139, trim(error_msg) &
            // " : " // trim(nf90_strerror(status)))
    end if

  end subroutine pmc_nc_check_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open a NetCDF file for reading.
  subroutine pmc_nc_open_read(filename, ncid)

    !> Filename of NetCDF file to open.
    character(len=*), intent(in) :: filename
    !> NetCDF file ID, in data mode.
    integer, intent(out) :: ncid

    call pmc_nc_check_msg(nf90_open(filename, NF90_NOWRITE, ncid), &
         "opening " // trim(filename))

  end subroutine pmc_nc_open_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close a NetCDF file.
  subroutine pmc_nc_close(ncid)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_check(nf90_close(ncid))

  end subroutine pmc_nc_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single real from a NetCDF file.
  subroutine pmc_nc_read_real(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real(kind=dp), intent(out) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = 0d0
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))
    
  end subroutine pmc_nc_read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single integer from a NetCDF file.
  subroutine pmc_nc_read_integer(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(out) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = 0
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))
    
  end subroutine pmc_nc_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple real array from a NetCDF file.
  subroutine pmc_nc_read_real_1d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    real(kind=dp), intent(out) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = 0d0
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))
    
  end subroutine pmc_nc_read_real_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple integer array from a NetCDF file.
  subroutine pmc_nc_read_integer_1d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    integer, intent(out) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = 0
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))
    
  end subroutine pmc_nc_read_integer_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple real 2D array from a NetCDF file.
  subroutine pmc_nc_read_real_2d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    real(kind=dp), intent(out) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = 0d0
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))
    
  end subroutine pmc_nc_read_real_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple integer 2D array from a NetCDF file.
  subroutine pmc_nc_read_integer_2d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    integer, intent(out) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = 0
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))
    
  end subroutine pmc_nc_read_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple logical 2D array from a NetCDF file.
  subroutine pmc_nc_read_logical_2d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    logical, intent(out) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    ! integer kind that can represent -1, 0, 1
    integer, parameter :: short = selected_int_kind(0)

    integer :: varid, status
    logical :: use_must_be_present
    integer(kind=short) :: int_var(size(var,1), size(var,2))

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = .false.
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, int_var), &
         "getting variable " // trim(name))
    where (int_var == 0_short)
       var = .false.
    elsewhere
       var = .true.
    end where
    
  end subroutine pmc_nc_read_logical_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write attributes for a variable to a NetCDF file.
  subroutine pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
       description)

    !> NetCDF file ID, in define mode.
    integer, intent(in) :: ncid
    !> Variable ID to write attributes for.
    integer, intent(in) :: varid
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    if (present(unit)) then
       call pmc_nc_check(nf90_put_att(ncid, varid, "unit", unit))
    end if
    if (present(long_name)) then
       call pmc_nc_check(nf90_put_att(ncid, varid, "long_name", long_name))
    end if
    if (present(standard_name)) then
       call pmc_nc_check(nf90_put_att(ncid, varid, "standard_name", &
            standard_name))
    end if
    if (present(description)) then
       call pmc_nc_check(nf90_put_att(ncid, varid, "description", &
            description))
    end if
    
  end subroutine pmc_nc_write_atts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a single real to a NetCDF file.
  subroutine pmc_nc_write_real(ncid, var, name, unit, long_name, &
       standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real(kind=dp), intent(in) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, dimids(0)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    call pmc_nc_check(nf90_put_var(ncid, varid, var))
    
  end subroutine pmc_nc_write_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a single integer to a NetCDF file.
  subroutine pmc_nc_write_integer(ncid, var, name, unit, long_name, &
       standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, dimids(0)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    call pmc_nc_check(nf90_put_var(ncid, varid, var))
    
  end subroutine pmc_nc_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple real array to a NetCDF file.
  subroutine pmc_nc_write_real_1d(ncid, var, name, dimids, unit, long_name, &
       standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real(kind=dp), intent(in) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(1)
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1 /)
    count = (/ size(var, 1) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_real_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple integer array to a NetCDF file.
  subroutine pmc_nc_write_integer_1d(ncid, var, name, dimids, unit, &
       long_name, standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(1)
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1 /)
    count = (/ size(var, 1) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_integer_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple real 2D array to a NetCDF file.
  subroutine pmc_nc_write_real_2d(ncid, var, name, dimids, unit, long_name, &
       standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real(kind=dp), intent(in) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(2)
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(2), count(2)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1, 1 /)
    count = (/ size(var, 1), size(var, 2) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_real_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple integer 2D array to a NetCDF file.
  subroutine pmc_nc_write_integer_2d(ncid, var, name, dimids, unit, &
       long_name, standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(2)
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(2), count(2)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1, 1 /)
    count = (/ size(var, 1), size(var, 2) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple logical 2D array to a NetCDF file.
  subroutine pmc_nc_write_logical_2d(ncid, var, name, dimids, unit, &
       long_name, standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    logical, intent(in) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(2)
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    ! integer kind that can represent -1, 0, 1
    integer, parameter :: short = selected_int_kind(0)

    integer :: varid, start(2), count(2)
    integer(kind=short) :: int_var(size(var,1), size(var,2))

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_BYTE, dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    where (var)
       int_var = 1_short
    elsewhere
       int_var = 0_short
    end where

    start = (/ 1, 1 /)
    count = (/ size(var, 1), size(var, 2) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, int_var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_logical_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_netcdf
