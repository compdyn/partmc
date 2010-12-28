! Copyright (C) 2009-2010 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_gas program.

!> Read NetCDF output files and write out the gas mixing ratios in text
!> format.
program extract_gas

  use netcdf

  integer, parameter :: dp = kind(0.d0)
  integer, parameter :: out_unit = 64

  character(len=1000) :: in_prefix, in_filename, out_filename
  integer :: ncid
  integer :: dimid_gas_species
  integer :: varid_time, varid_gas_species
  integer :: varid_gas_mixing_ratio
  integer :: n_gas_species
  character(len=1000) :: tmp_str, gas_species_names, remaining_species
  real(kind=dp) :: time
  real(kind=dp), allocatable :: gas_mixing_ratio(:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec, status, n_time, i
  character(len=36) :: uuid, run_uuid

  ! process commandline arguments
  if (command_argument_count() .ne. 2) then
     write(6,*) 'Usage: extract_gas <netcdf_state_prefix> <output_filename>'
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
  write(*,'(a)') "    column  1: time (s)"

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

     ! read and check uuid
     call nc_check_msg(nf90_get_att(ncid, NF90_GLOBAL, "UUID", uuid), &
          "getting global attribute 'UUID'")
     if (i_time == 1) then
        run_uuid = uuid
     else
        if (run_uuid /= uuid) then
           write(0,*) 'ERROR: UUID mismatch at: ' // trim(in_filename)
           stop 1
        end if
     end if

     ! read time
     call nc_check_msg(nf90_inq_varid(ncid, "time", varid_time), &
          "getting variable ID for 'time'")
     call nc_check_msg(nf90_get_var(ncid, varid_time, time), &
          "getting variable 'time'")

     ! read gas_species
     call nc_check_msg(nf90_inq_dimid(ncid, "gas_species", &
          dimid_gas_species), "getting dimension ID for 'gas_species'")
     call nc_check_msg(nf90_Inquire_Dimension(ncid, dimid_gas_species, &
          tmp_str, n_gas_species), "inquiring dimension 'gas_species'")
     call nc_check_msg(nf90_inq_varid(ncid, "gas_species", &
          varid_gas_species), "getting variable ID for 'gas_species'")
     call nc_check_msg(nf90_get_att(ncid, varid_gas_species, &
          "names", gas_species_names), &
          "getting attribute 'names' for variable 'gas_species'")
     if (i_time == 1) then
        remaining_species = gas_species_names
        do i_spec = 1,n_gas_species
           if (i_spec < n_gas_species) then
              i = index(remaining_species, ',')
           else
              i = index(remaining_species, ' ')
           end if
           write(*,'(a,i2,a,a,a)') '    column ', i_spec + 1, ': ', &
                remaining_species(1:(i-1)), ' (ppb)'
           remaining_species = remaining_species((i+1):)
        end do
     end if
     
     ! read gas_mixing_ratio
     call nc_check_msg(nf90_inq_varid(ncid, "gas_mixing_ratio", &
          varid_gas_mixing_ratio), &
          "getting variable ID for 'gas_mixing_ratio'")
     call nc_check_msg(nf90_Inquire_Variable(ncid, varid_gas_mixing_ratio, &
          tmp_str, xtype, ndims, dimids, nAtts), &
          "inquiring variable 'gas_mixing_ratio'")
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_gas_species)) then
        write(0,*) "ERROR: unexpected gas_mixing_ratio dimids"
        stop 1
     end if
     allocate(gas_mixing_ratio(n_gas_species))
     call nc_check_msg(nf90_get_var(ncid, varid_gas_mixing_ratio, &
          gas_mixing_ratio), "getting variable 'gas_mixing_ratio'")
     
     call nc_check_msg(nf90_close(ncid), &
          "closing file " // trim(in_filename))

     ! output data
     write(out_unit, '(e30.15e3)', advance='no') time
     do i_spec = 1,n_gas_species
        write(out_unit, '(e30.15e3)', advance='no') &
             gas_mixing_ratio(i_spec)
     end do
     write(out_unit, '(a)') ''

     deallocate(gas_mixing_ratio)
  end do

  if (n_time == 0) then
     write(*,'(a,a)') 'ERROR: no input file found matching: ', &
          trim(in_filename)
     stop 1
  end if

  close(out_unit)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check return status of NetCDF function calls.
  subroutine nc_check_msg(status, error_msg)

    !> Status return value.
    integer, intent(in) :: status
    !> Error message in case of failure.
    character(len=*), intent(in) :: error_msg

    if (status /= NF90_NOERR) then
       write(0,*) trim(error_msg) // " : " // trim(nf90_strerror(status))
       stop 1
    end if

  end subroutine nc_check_msg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEFINE_LOCAL_COMMAND_ARGUMENT
  integer function command_argument_count()
    command_argument_count = iargc()
  end function command_argument_count
  subroutine get_command_argument(i, arg)
    integer, intent(in) :: i
    character(len=*), intent(out) :: arg
    call getarg(i, arg)
  end subroutine get_command_argument
#endif
  
end program extract_gas
