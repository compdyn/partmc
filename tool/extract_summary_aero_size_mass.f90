! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Read a NetCDF summary file and write out the size- and time-resolved
! mass concentration in text format.

program extract_summary_aero_size_mass

  use netcdf

  integer, parameter :: out_unit = 64

  character(len=1000) :: in_filename, out_filename
  integer :: ncid
  integer :: dimid_time, dimid_radius, dimid_aero_species, dimid_unit
  integer :: varid_time, varid_radius, varid_aero_species, varid_unit
  integer :: varid_aero
  integer :: n_time, n_radius, n_aero_species, n_unit
  character(len=1000) :: tmp_str, aero_species_names, data_units
  real*8, allocatable :: time(:)
  real*8, allocatable :: radius(:)
  real*8, allocatable :: aero(:,:,:,:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_radius, i_spec, i_unit
  real*8 :: val

  ! process commandline arguments
  if (iargc() .ne. 2) then
     write(6,*) 'Usage: extract_summary_size_mass <netcdf_filename> <output_filename>'
     call exit(2)
  endif
  call getarg(1, in_filename)
  call getarg(2, out_filename)

  ! write information
  write(*,*) "Output file array A has:"
  write(*,*) "  A(1, j+1) = radius(j) (m)"
  write(*,*) "  A(i+1, 1) = time(i) (s)"
  write(*,*) "  A(i+1, j+1) = mass concentration at time(i) and radius(j) (kg/m^3)"

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

  call nc_check(nf90_inq_dimid(ncid, "radius", dimid_radius))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_radius, &
       tmp_str, n_radius))
  allocate(radius(n_radius))
  call nc_check(nf90_inq_varid(ncid, "radius", varid_radius))
  call nc_check(nf90_get_var(ncid, varid_radius, radius))
  write(*,*) "n_radius:", n_radius
  write(*,*) "min radius:", minval(radius)
  write(*,*) "max radius:", maxval(radius)

  call nc_check(nf90_inq_dimid(ncid, "aero_species", dimid_aero_species))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_species, &
       tmp_str, n_aero_species))
  call nc_check(nf90_inq_varid(ncid, "aero_species", varid_aero_species))
  call nc_check(nf90_get_att(ncid, varid_aero_species, &
       "names", aero_species_names))
  write(*,*) "n_aero_species:", n_aero_species
  write(*,*) "aero_species_names: ", trim(aero_species_names)

  call nc_check(nf90_inq_dimid(ncid, "unit", dimid_unit))
  call nc_check(nf90_Inquire_Dimension(ncid, dimid_unit, &
       tmp_str, n_unit))
  if (n_unit /= 4) then
     write(*,*) "ERROR: unexpected number of units"
     call exit(1)
  end if
  call nc_check(nf90_inq_varid(ncid, "unit", varid_unit))
  call nc_check(nf90_get_att(ncid, varid_unit, &
       "data_units", data_units))
  write(*,*) "n_unit:", n_unit
  write(*,*) "data_units: ", trim(data_units)

  call nc_check(nf90_inq_varid(ncid, "aero", varid_aero))
  call nc_check(nf90_Inquire_Variable(ncid, varid_aero, tmp_str, &
       xtype, ndims, dimids, nAtts))
  if ((ndims /= 4) &
       .or. (dimids(1) /= dimid_radius) &
       .or. (dimids(2) /= dimid_aero_species) &
       .or. (dimids(3) /= dimid_unit) &
       .or. (dimids(4) /= dimid_time)) then
     write(*,*) "ERROR: unexpected aero dimids"
     call exit(1)
  end if
  allocate(aero(n_radius, n_aero_species, n_unit, n_time))
  call nc_check(nf90_get_var(ncid, varid_aero, aero))

  call nc_check(nf90_close(ncid))

  ! output data
  i_unit = 3 ! mass
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if
  write(out_unit, '(e30.15e3)', advance='no') 0d0
  do i_radius = 1,n_radius
     write(out_unit, '(e30.15e3)', advance='no') radius(i_radius)
  end do
  write(out_unit, '(a)') ''
  do i_time = 1,n_time
     write(out_unit, '(e30.15e3)', advance='no') time(i_time)
     do i_radius = 1,n_radius
        val = 0d0
        do i_spec = 1,n_aero_species
           val = val + aero(i_radius, i_spec, i_unit, i_time)
        end do
        write(out_unit, '(e30.15e3)', advance='no') val
     end do
     write(out_unit, '(a)') ''
  end do
  close(out_unit)

  deallocate(time)
  deallocate(radius)
  deallocate(aero)

contains

  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

end program extract_summary_aero_size_mass
