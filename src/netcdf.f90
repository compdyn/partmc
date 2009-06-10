! Copyright (C) 2007-2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_netcdf module.

!> Wrapper functions for NetCDF.
module pmc_netcdf

  use netcdf
  use pmc_util

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the status of a NetCDF function call.
  subroutine pmc_nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       call die_msg(291021908, nf90_strerror(status))
    end if

  end subroutine pmc_nc_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check the status of a NetCDF function call and prints the given
  !> error message on failure.
  subroutine pmc_nc_check_msg(status, error_msg)

    !> Status return value.
    integer, intent(in) :: status
    !> Error message in case of failure.
    character(len=*), intent(in) :: error_msg

    if (status /= NF90_NOERR) then
       call die_msg(291021908, trim(error_msg) &
            // " : " // trim(nf90_strerror(status)))
    end if

  end subroutine pmc_nc_check_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open a NetCDF file for reading.
  subroutine pmc_nc_open_read(filename, ncid)

    !> Filename of NetCDF file to open.
    character(len=*), intent(in) :: filename
    !> NetCDF file ID, in data mode.
    integer, intent(out) :: ncid

    call pmc_nc_check_msg(nf90_open(filename, NF90_NOWRITE, ncid), &
         "opening " // trim(filename))

  end subroutine pmc_nc_open_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close a NetCDF file.
  subroutine pmc_nc_close(ncid)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_check(nf90_close(ncid))

  end subroutine pmc_nc_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single real from a NetCDF file.
  subroutine pmc_nc_read_real(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real*8, intent(out) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(out) :: unit

    integer :: varid

    call pmc_nc_check(nf90_inq_varid(ncid, name, varid))
    call pmc_nc_check(nf90_get_var(ncid, varid, var))
    call pmc_nc_check(nf90_get_att(ncid, varid, "unit", unit))
    
  end subroutine pmc_nc_read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a single integer from a NetCDF file.
  subroutine pmc_nc_read_integer(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(out) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(out) :: unit

    integer :: varid

    call pmc_nc_check(nf90_inq_varid(ncid, name, varid))
    call pmc_nc_check(nf90_get_var(ncid, varid, var))
    call pmc_nc_check(nf90_get_att(ncid, varid, "unit", unit))
    
  end subroutine pmc_nc_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple real array from a NetCDF file.
  subroutine pmc_nc_read_real_1d(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    real*8, intent(out) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(out) :: unit

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_inq_varid(ncid, name, varid))
    call pmc_nc_check(nf90_get_var(ncid, varid, var))
    call pmc_nc_check(nf90_get_att(ncid, varid, "unit", unit))
    
  end subroutine pmc_nc_read_real_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple integer array from a NetCDF file.
  subroutine pmc_nc_read_integer_1d(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    integer, intent(out) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(out) :: unit

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_inq_varid(ncid, name, varid))
    call pmc_nc_check(nf90_get_var(ncid, varid, var))
    call pmc_nc_check(nf90_get_att(ncid, varid, "unit", unit))
    
  end subroutine pmc_nc_read_integer_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple real 2D array from a NetCDF file.
  subroutine pmc_nc_read_real_2d(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    real*8, intent(out) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(out) :: unit

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_inq_varid(ncid, name, varid))
    call pmc_nc_check(nf90_get_var(ncid, varid, var))
    call pmc_nc_check(nf90_get_att(ncid, varid, "unit", unit))
    
  end subroutine pmc_nc_read_real_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple integer 2D array from a NetCDF file.
  subroutine pmc_nc_read_integer_2d(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read, must be correctly sized.
    integer, intent(out) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(out) :: unit

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_inq_varid(ncid, name, varid))
    call pmc_nc_check(nf90_get_var(ncid, varid, var))
    call pmc_nc_check(nf90_get_att(ncid, varid, "unit", unit))
    
  end subroutine pmc_nc_read_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a single real to a NetCDF file.
  subroutine pmc_nc_write_real(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real*8, intent(in) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(in) :: unit

    integer :: varid, dimids(0)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimids, varid))
    call pmc_nc_check(nf90_put_att(ncid, varid, "unit", unit))
    call pmc_nc_check(nf90_enddef(ncid))

    call pmc_nc_check(nf90_put_var(ncid, varid, var))
    
  end subroutine pmc_nc_write_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a single integer to a NetCDF file.
  subroutine pmc_nc_write_integer(ncid, var, name, unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(in) :: unit

    integer :: varid, dimids(0)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, dimids, varid))
    call pmc_nc_check(nf90_put_att(ncid, varid, "unit", unit))
    call pmc_nc_check(nf90_enddef(ncid))

    call pmc_nc_check(nf90_put_var(ncid, varid, var))
    
  end subroutine pmc_nc_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple real array to a NetCDF file.
  subroutine pmc_nc_write_real_1d(ncid, var, name, unit, dimids)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real*8, intent(in) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(in) :: unit
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(1)

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimids, varid))
    call pmc_nc_check(nf90_put_att(ncid, varid, "unit", unit))
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1 /)
    count = (/ size(var, 1) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_real_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple integer array to a NetCDF file.
  subroutine pmc_nc_write_integer_1d(ncid, var, name, unit, dimids)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(in) :: unit
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(1)

    integer :: varid, start(1), count(1)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, dimids, varid))
    call pmc_nc_check(nf90_put_att(ncid, varid, "unit", unit))
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1 /)
    count = (/ size(var, 1) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_integer_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple real 2D array to a NetCDF file.
  subroutine pmc_nc_write_real_2d(ncid, var, name, unit, dimids)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real*8, intent(in) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(in) :: unit
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(2)

    integer :: varid, start(2), count(2)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, dimids, varid))
    call pmc_nc_check(nf90_put_att(ncid, varid, "unit", unit))
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1, 1 /)
    count = (/ size(var, 1), size(var, 2) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_real_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple integer 2D array to a NetCDF file.
  subroutine pmc_nc_write_integer_2d(ncid, var, name, unit, dimids)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), intent(in) :: unit
    !> NetCDF dimension IDs of the variable
    integer, intent(in) :: dimids(2)

    integer :: varid, start(2), count(2)

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, dimids, varid))
    call pmc_nc_check(nf90_put_att(ncid, varid, "unit", unit))
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1, 1 /)
    count = (/ size(var, 1), size(var, 2) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))
    
  end subroutine pmc_nc_write_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_netcdf
