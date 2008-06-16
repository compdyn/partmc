! Copyright (C) 2007, 2008 Matthew West
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

  !> Write a simple array to a NetCDF file.
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

  !> Write a simple array to a NetCDF file.
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

  !> Write a simple array to a NetCDF file.
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

  !> Write a simple array to a NetCDF file.
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
