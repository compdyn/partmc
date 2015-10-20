! Copyright (C) 2007-2010, 2012 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_netcdf module.

!> Wrapper functions for NetCDF. These all take a NetCDF \c ncid in
!> data mode and return with it again in data mode. Shifting to define
!> mode is handled internally within each subroutine.
module pmc_netcdf

  use netcdf
  use pmc_util
  use pmc_rand

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
         "opening " // trim(filename) // " for reading")

  end subroutine pmc_nc_open_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Open a NetCDF file for writing.
  subroutine pmc_nc_open_write(filename, ncid)

    !> Filename of NetCDF file to open.
    character(len=*), intent(in) :: filename
    !> NetCDF file ID, in data mode, returns in data mode.
    integer, intent(out) :: ncid

    call pmc_nc_check_msg(nf90_create(filename, NF90_CLOBBER, ncid), &
         "opening " // trim(filename) // " for writing")
    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine pmc_nc_open_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close a NetCDF file.
  subroutine pmc_nc_close(ncid)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid

    call pmc_nc_check_msg(nf90_close(ncid), "closing NetCDF file")

  end subroutine pmc_nc_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write basic information to a NetCDF file.
  subroutine pmc_nc_write_info(ncid, uuid, source, write_rank, write_n_proc)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> UUID for this data set.
    character(len=PMC_UUID_LEN), intent(in) :: uuid
    !> Source name for this data.
    character(len=*), intent(in) :: source
    !> Rank to write into file.
    integer, intent(in), optional :: write_rank
    !> Number of processes to write into file.
    integer, intent(in), optional :: write_n_proc

    character(len=(len_trim(source) + 100)) :: history
    integer :: use_rank, use_n_proc

    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
         trim(source) // " output file"))
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "source", source))
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "UUID", uuid))
    call iso8601_date_and_time(history)
    history((len_trim(history)+1):) = (" created by " // trim(source))
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "history", history))
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.4"))
    call pmc_nc_check(nf90_enddef(ncid))

    if (present(write_rank)) then
       use_rank = write_rank
    else
       use_rank = pmc_mpi_rank()
    end if
    if (present(write_n_proc)) then
       use_n_proc = write_n_proc
    else
       use_n_proc = pmc_mpi_size()
    end if

#ifdef PMC_USE_MPI
    call pmc_nc_write_integer(ncid, use_rank + 1, "process", &
         description="the process number (starting from 1) " &
         // "that output this data file")
    call pmc_nc_write_integer(ncid, use_n_proc, "total_processes", &
         description="total number of processes")
#endif

  end subroutine pmc_nc_write_info

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
    !> Data to read.
    real(kind=dp), intent(inout), allocatable :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status, dimids(1), size1
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok, set it to empty
       var = [real(kind=dp)::]
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_variable(ncid, varid, dimids=dimids), &
         "determining size of variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_dimension(ncid, dimids(1), len=size1), &
         "determining size of dimension number " &
         // trim(integer_to_string(dimids(1))))
    call ensure_real_array_size(var, size1)
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))

  end subroutine pmc_nc_read_real_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple integer array from a NetCDF file.
  subroutine pmc_nc_read_integer_1d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read.
    integer, intent(inout), allocatable :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status, dimids(1), size1
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = [integer::]
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_variable(ncid, varid, dimids=dimids), &
         "determining size of variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_dimension(ncid, dimids(1), len=size1), &
         "determining size of dimension number " &
         // trim(integer_to_string(dimids(1))))
    call ensure_integer_array_size(var, size1)
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))

  end subroutine pmc_nc_read_integer_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple real 2D array from a NetCDF file.
  subroutine pmc_nc_read_real_2d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read.
    real(kind=dp), intent(inout), allocatable :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status, dimids(2), size1, size2
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = reshape([real(kind=dp)::], [0, 0])
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_variable(ncid, varid, dimids=dimids), &
         "determining size of variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_dimension(ncid, dimids(1), len=size1), &
         "determining size of dimension number " &
         // trim(integer_to_string(dimids(1))))
    call pmc_nc_check_msg(nf90_inquire_dimension(ncid, dimids(2), len=size2), &
         "determining size of dimension number " &
         // trim(integer_to_string(dimids(2))))
    call ensure_real_array_2d_size(var, size1, size2)
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))

  end subroutine pmc_nc_read_real_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read a simple integer 2D array from a NetCDF file.
  subroutine pmc_nc_read_integer_2d(ncid, var, name, must_be_present)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to read.
    integer, intent(inout), allocatable :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Whether the variable must be present in the NetCDF file
    !> (default .true.).
    logical, optional, intent(in) :: must_be_present

    integer :: varid, status, dimids(2), size1, size2
    logical :: use_must_be_present

    if (present(must_be_present)) then
       use_must_be_present = must_be_present
    else
       use_must_be_present = .true.
    end if
    status = nf90_inq_varid(ncid, name, varid)
    if ((.not. use_must_be_present) .and. (status == NF90_ENOTVAR)) then
       ! variable was not present, but that's ok
       var = reshape([integer::], [0, 0])
       return
    end if
    call pmc_nc_check_msg(status, "inquiring variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_variable(ncid, varid, dimids=dimids), &
         "determining size of variable " // trim(name))
    call pmc_nc_check_msg(nf90_inquire_dimension(ncid, dimids(1), len=size1), &
         "determining size of dimension number " &
         // trim(integer_to_string(dimids(1))))
    call pmc_nc_check_msg(nf90_inquire_dimension(ncid, dimids(2), len=size2), &
         "determining size of dimension number " &
         // trim(integer_to_string(dimids(2))))
    call ensure_integer_array_2d_size(var, size1, size2)
    call pmc_nc_check_msg(nf90_get_var(ncid, varid, var), &
         "getting variable " // trim(name))

  end subroutine pmc_nc_read_integer_2d

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

  !> Create a dimension if necessary, or check its size if it already
  !> exists. In any case return the \c dimid.
  subroutine pmc_nc_ensure_dim(ncid, dim_name, dimid, dim_size, array_dim)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> NetCDF dimension name for the variable.
    character(len=*), intent(in) :: dim_name
    !> NetCDF dimension ID.
    integer, intent(out) :: dimid
    !> Size of dimension.
    integer, intent(in) :: dim_size
    !> Dimension within data array that this NetCDF dim corresponds to.
    integer, intent(in) :: array_dim

    integer :: status, check_dim_size
    character(len=NF90_MAX_NAME) :: check_name

    status = nf90_inq_dimid(ncid, dim_name, dimid)
    if (status == NF90_NOERR) then
       call pmc_nc_check(nf90_Inquire_Dimension(ncid, dimid, check_name, &
            check_dim_size))
       call assert_msg(657263912, check_dim_size == dim_size, &
            "dim " // trim(integer_to_string(array_dim)) // " size " &
            // trim(integer_to_string(dim_size)) &
            // " of data array does not match size " &
            // trim(integer_to_string(dim_size)) // " of '" &
            // trim(dim_name) // "' dim")
    else
       ! could not determine dimid
       if (status /= NF90_EBADDIM) then
          ! the problem was not a missing dimension
          call pmc_nc_check(status)
       end if
       ! the problem was a missing dimension, so make it
       call pmc_nc_check(nf90_redef(ncid))
       call pmc_nc_check(nf90_def_dim(ncid, dim_name, dim_size, dimid))
       call pmc_nc_check(nf90_enddef(ncid))
    end if

  end subroutine pmc_nc_ensure_dim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write a simple real array to a NetCDF file.
  subroutine pmc_nc_write_real_1d(ncid, var, name, dimids, dim_name, unit, &
       long_name, standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real(kind=dp), intent(in) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable (either this or dim_name
    !> must be present).
    integer, optional, intent(in) :: dimids(1)
    !> NetCDF dimension name for the variable (either this or dimids
    !> must be present).
    character(len=*), optional, intent(in) :: dim_name
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(1), count(1), use_dimids(1)

    if (present(dimids)) then
       use_dimids = dimids
    elseif (present(dim_name)) then
       call pmc_nc_ensure_dim(ncid, dim_name, use_dimids(1), size(var), 1)
    else
       call die_msg(891890123, "either dimids or dim_name must be present")
    end if
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, use_dimids, varid))
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
  subroutine pmc_nc_write_integer_1d(ncid, var, name, dimids, dim_name, unit, &
       long_name, standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var(:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable (either this or dim_name
    !> must be present).
    integer, optional, intent(in) :: dimids(1)
    !> NetCDF dimension name for the variable (either this or dimids
    !> must be present).
    character(len=*), optional, intent(in) :: dim_name
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(1), count(1), use_dimids(1)

    if (present(dimids)) then
       use_dimids = dimids
    elseif (present(dim_name)) then
       call pmc_nc_ensure_dim(ncid, dim_name, use_dimids(1), size(var), 1)
    else
       call die_msg(464170526, "either dimids or dim_name must be present")
    end if
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, use_dimids, varid))
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
  subroutine pmc_nc_write_real_2d(ncid, var, name, dimids, dim_name_1, &
       dim_name_2, unit, long_name, standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    real(kind=dp), intent(in) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable (either \c dimids or both
    !> \c dim_name_1 and \c dim_name_2 must be present).
    integer, optional, intent(in) :: dimids(2)
    !> First NetCDF dimension name for the variable (either \c dimids
    !> or both \c dim_name_1 and \c dim_name 2 must be present).
    character(len=*), optional, intent(in) :: dim_name_1
    !> Second NetCDF dimension name for the variable (either \c dimids
    !> or both \c dim_name_1 and \c dim_name 2 must be present).
    character(len=*), optional, intent(in) :: dim_name_2
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(2), count(2), use_dimids(2)

    if (present(dimids)) then
       use_dimids = dimids
    elseif (present(dim_name_1) .and. present(dim_name_2)) then
       call pmc_nc_ensure_dim(ncid, dim_name_1, use_dimids(1), size(var, 1), 1)
       call pmc_nc_ensure_dim(ncid, dim_name_2, use_dimids(2), size(var, 2), 2)
    else
       call die_msg(959111259, &
            "either dimids or both dim_name_1 and dim_name_2 must be present")
    end if
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_DOUBLE, use_dimids, varid))
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
  subroutine pmc_nc_write_integer_2d(ncid, var, name, dimids, dim_name_1, &
       dim_name_2, unit, long_name, standard_name, description)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Data to write.
    integer, intent(in) :: var(:,:)
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension IDs of the variable (either \c dimids or both
    !> \c dim_name_1 and \c dim_name_2 must be present).
    integer, optional, intent(in) :: dimids(2)
    !> First NetCDF dimension name for the variable (either \c dimids
    !> or both \c dim_name_1 and \c dim_name 2 must be present).
    character(len=*), optional, intent(in) :: dim_name_1
    !> Second NetCDF dimension name for the variable (either \c dimids
    !> or both \c dim_name_1 and \c dim_name 2 must be present).
    character(len=*), optional, intent(in) :: dim_name_2
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit
    !> Long name of variable.
    character(len=*), optional, intent(in) :: long_name
    !> Standard name of variable.
    character(len=*), optional, intent(in) :: standard_name
    !> Description of variable.
    character(len=*), optional, intent(in) :: description

    integer :: varid, start(2), count(2), use_dimids(2)

    if (present(dimids)) then
       use_dimids = dimids
    elseif (present(dim_name_1) .and. present(dim_name_2)) then
       call pmc_nc_ensure_dim(ncid, dim_name_1, use_dimids(1), size(var, 1), 1)
       call pmc_nc_ensure_dim(ncid, dim_name_2, use_dimids(2), size(var, 2), 2)
    else
       call die_msg(669381383, &
            "either dimids or both dim_name_1 and dim_name_2 must be present")
    end if
    call pmc_nc_check(nf90_redef(ncid))
    call pmc_nc_check(nf90_def_var(ncid, name, NF90_INT, use_dimids, varid))
    call pmc_nc_write_atts(ncid, varid, unit, long_name, standard_name, &
         description)
    call pmc_nc_check(nf90_enddef(ncid))

    start = (/ 1, 1 /)
    count = (/ size(var, 1), size(var, 2) /)
    call pmc_nc_check(nf90_put_var(ncid, varid, var, &
         start = start, count = count))

  end subroutine pmc_nc_write_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_netcdf
