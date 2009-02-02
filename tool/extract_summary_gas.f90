! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Read a NetCDF summary file and write out the gas concenrations in
! text format.

program extract_summary_gas

  use netcdf

  integer, parameter :: out_unit = 64

  character(len=1000) :: in_filename, out_filename
  integer :: ncid
  integer :: dimid_time, dimid_gas_species
  integer :: varid_time, varid_gas_species
  integer :: varid_gas
  integer :: n_time, n_gas_species
  character(len=1000) :: tmp_str, gas_species_names
  real*8, allocatable :: time(:)
  real*8, allocatable :: gas(:,:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec
  real*8 :: val

  ! process commandline arguments
  if (iargc() .ne. 2) then
     write(6,*) 'Usage: extract_summary_gas <netcdf_filename> <output_filename>'
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

  call nc_check(nf90_inq_dimid(ncid, "gas_species", dimid_gas_species))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_gas_species, &
       tmp_str, n_gas_species))
  call nc_check(nf90_inq_varid(ncid, "gas_species", varid_gas_species))
  call nc_check(nf90_get_att(ncid, varid_gas_species, &
       "names", gas_species_names))
  write(*,*) "n_gas_species:", n_gas_species
  write(*,*) "gas_species_names: ", trim(gas_species_names)

  call nc_check(nf90_inq_varid(ncid, "gas", varid_gas))
  call nc_check(nf90_Inquire_Variable(ncid, varid_gas, tmp_str, &
       xtype, ndims, dimids, nAtts))
  if ((ndims /= 2) &
       .or. (dimids(1) /= dimid_gas_species) &
       .or. (dimids(2) /= dimid_time)) then
     write(*,*) "ERROR: unexpected gas dimids"
     call exit(1)
  end if
  allocate(gas(n_gas_species, n_time))
  call nc_check(nf90_get_var(ncid, varid_gas, gas))

  call nc_check(nf90_close(ncid))

  ! output data
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if
  do i_time = 1,n_time
     write(out_unit, '(e30.15e6)', advance='no') time(i_time)
     do i_spec = 1,n_gas_species
        write(out_unit, '(e30.15e6)', advance='no') gas(i_spec, i_time)
     end do
     write(out_unit, '(a)') ''
  end do
  close(out_unit)

  deallocate(time)
  deallocate(gas)

contains

  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

end program extract_summary_gas
