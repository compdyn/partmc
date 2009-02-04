! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Read NetCDF state files and write out the aerosol concentrations in
! text format.

program extract_state_aero_size_num

  use netcdf

  integer, parameter :: out_unit = 64

  character(len=1000) :: in_prefix, in_filename, out_filename
  integer :: ncid
  integer :: dimid_aero_species, dimid_aero_particle
  integer :: varid_time, varid_aero_species
  integer :: varid_aero_comp_mass, varid_aero_density
  integer :: varid_aero_comp_vol
  integer :: n_aero_species, n_aero_particle
  character(len=1000) :: tmp_str, aero_species_names
  real*8 :: time
  real*8, allocatable :: aero_comp_mass(:,:)
  real*8, allocatable :: aero_density(:)
  real*8, allocatable :: aero_comp_vol(:)
  real*8, allocatable :: aero_dist(:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec, i_part, status
  integer :: n_bin, i_bin
  real*8 :: r_min, r_max, radius, volume

  ! process commandline arguments
  if (iargc() .ne. 5) then
     write(6,*) 'Usage: extract_state_aero_size_num <r_min> <r_max> <n_bin> <netcdf_state_prefix> <output_filename>'
     call exit(2)
  endif
  call getarg(1, tmp_str)
  r_min = string_to_real(tmp_str)
  call getarg(2, tmp_str)
  r_max = string_to_real(tmp_str)
  call getarg(3, tmp_str)
  n_bin = string_to_integer(tmp_str)
  call getarg(4, in_prefix)
  call getarg(5, out_filename)

  ! open output file
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if
  write(out_unit, '(e30.15e3)', advance='no') 0d0
  do i_bin = 1,n_bin
     radius = exp(dble(i_bin - 1) / dble(n_bin - 1) &
          * (log(r_max) - log(r_min)) + log(r_min))
     write(out_unit, '(e30.15e3)', advance='no') radius
  end do
  write(out_unit, '(a)') ''

  ! write information
  write(*,*) "Output file array A has:"
  write(*,*) "  A(1, j+1) = radius(j) (m)"
  write(*,*) "  A(i+1, 1) = time(i) (s)"
  write(*,*) "  A(i+1, j+1) = mass concentration at time(i) and radius(j) (kg/m^3)"

  ! process NetCDF files
  allocate(aero_dist(n_bin))
  i_time = 0
  do while (.true.)
     i_time = i_time + 1
     write(in_filename,'(a,i8.8,a)') trim(in_prefix), i_time, ".nc"
     status = nf90_open(in_filename, NF90_NOWRITE, ncid)
     if (status /= NF90_NOERR) then
        exit
     end if

     ! read time
     call nc_check(nf90_inq_varid(ncid, "time", varid_time))
     call nc_check(nf90_get_var(ncid, varid_time, time))

     ! read aero_species
     call nc_check(nf90_inq_dimid(ncid, "aero_species", dimid_aero_species))
     call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_species, &
          tmp_str, n_aero_species))
     call nc_check(nf90_inq_varid(ncid, "aero_species", varid_aero_species))
     call nc_check(nf90_get_att(ncid, varid_aero_species, &
          "names", aero_species_names))
     if (i_time == 1) then
        write(*,*) "n_aero_species:", n_aero_species
        write(*,*) "aero_species_names: ", trim(aero_species_names)
     end if
     
     ! read aero_particle dimension
     call nc_check(nf90_inq_dimid(ncid, "aero_particle", dimid_aero_particle))
     call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_particle, &
          tmp_str, n_aero_particle))
     
     ! read aero_comp_mass
     call nc_check(nf90_inq_varid(ncid, "aero_comp_mass", &
          varid_aero_comp_mass))
     call nc_check(nf90_Inquire_Variable(ncid, varid_aero_comp_mass, &
          tmp_str, xtype, ndims, dimids, nAtts))
     if ((ndims /= 2) &
          .or. (dimids(1) /= dimid_aero_particle) &
          .or. (dimids(2) /= dimid_aero_species)) then
        write(*,*) "ERROR: unexpected aero_comp_mass dimids"
        call exit(1)
     end if
     allocate(aero_comp_mass(n_aero_particle, n_aero_species))
     call nc_check(nf90_get_var(ncid, varid_aero_comp_mass, &
          aero_comp_mass))
     
     ! read aero_density
     call nc_check(nf90_inq_varid(ncid, "aero_density", &
          varid_aero_density))
     call nc_check(nf90_Inquire_Variable(ncid, varid_aero_density, &
          tmp_str, xtype, ndims, dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_species)) then
        write(*,*) "ERROR: unexpected aero_density dimids"
        call exit(1)
     end if
     allocate(aero_density(n_aero_species))
     call nc_check(nf90_get_var(ncid, varid_aero_density, &
          aero_density))
     
     ! read aero_comp_vol
     call nc_check(nf90_inq_varid(ncid, "aero_comp_vol", &
          varid_aero_comp_vol))
     call nc_check(nf90_Inquire_Variable(ncid, varid_aero_comp_vol, &
          tmp_str, xtype, ndims, dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_particle)) then
        write(*,*) "ERROR: unexpected aero_comp_vol dimids"
        call exit(1)
     end if
     allocate(aero_comp_vol(n_aero_particle))
     call nc_check(nf90_get_var(ncid, varid_aero_comp_vol, &
          aero_comp_vol))
     
     call nc_check(nf90_close(ncid))

     ! compute distribution
     aero_dist = 0d0
     do i_part = 1,n_aero_particle
        volume = sum(aero_comp_mass(i_part,:) / aero_density)
        radius = (volume / (4d0 / 3d0 &
             * 3.14159265358979323846d0))**(1d0/3d0)
        i_bin = ceiling((log(radius) - log(r_min)) &
             / (log(r_max) - log(r_min)) * dble(n_bin - 1) + 0.5d0)
        aero_dist(i_bin) = aero_dist(i_bin) &
             + sum(aero_comp_mass(i_part,:)) / aero_comp_vol(i_part)
     end do

     ! output data
     write(out_unit, '(e30.15e3)', advance='no') time
     do i_bin = 1,n_bin
        write(out_unit, '(e30.15e3)', advance='no') &
             aero_dist(i_bin)
     end do
     write(out_unit, '(a)') ''

     deallocate(aero_comp_mass)
     deallocate(aero_density)
     deallocate(aero_comp_vol)
  end do

  close(out_unit)
  deallocate(aero_dist)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check return status of NetCDF function calls.
  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to a real.
  real*8 function string_to_real(string)

    !> String to convert.
    character(len=*), intent(in) :: string
    
    real*8 :: val
    integer :: ios

    read(string, '(e40.0)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', trim(string), &
            ' to real: IOSTAT = ', ios
       call exit(2)
    end if
    string_to_real = val

  end function string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Convert a string to an integer.
  integer function string_to_integer(string)

    !> String to convert.
    character(len=*), intent(in) :: string
    
    integer :: val
    integer :: ios

    read(string, '(i20)', iostat=ios) val
    if (ios /= 0) then
       write(0,'(a,a,a,i3)') 'Error converting ', trim(string), &
            ' to integer: IOSTAT = ', ios
       call exit(1)
    end if
    string_to_integer = val

  end function string_to_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program extract_state_aero_size_num
