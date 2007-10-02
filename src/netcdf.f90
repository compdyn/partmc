! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! NetCDF support functions.

module pmc_netcdf

  use netcdf
  use pmc_util

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pmc_nc_check(status)

    ! Check the status of a NetCDF function call.

    integer, intent(in) :: status       ! status return value

    if (status /= NF90_NOERR) then
       call die_msg(291021908, nf90_strerror(status))
    end if

  end subroutine pmc_nc_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_netcdf
