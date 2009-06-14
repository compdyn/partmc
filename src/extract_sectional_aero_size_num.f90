! Copyright (C) 2009 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Read NetCDF sectional files and write out the aerosol concentrations
! in text format.

program extract_sectional_aero_size_num

  use netcdf

  integer, parameter :: MAX_N_TIME = 10000
  integer, parameter :: out_unit = 64

  character(len=1000) :: in_prefix, in_filename, out_filename
  integer :: ncid
  integer :: dimid_aero_radius
  integer :: varid_time, varid_aero_number_concentration
  integer :: varid_aero_radius, varid_aero_radius_widths
  integer :: n_bin
  character(len=1000) :: tmp_str
  real*8 :: time
  real*8, allocatable :: aero_dist(:,:)
  real*8, allocatable :: aero_radius(:)
  real*8, allocatable :: aero_radius_widths(:)
  real*8, allocatable :: aero_number_concentration(:)
  real*8, allocatable :: save_aero_radius(:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec, i_part, status
  integer :: i_bin, n_time, new_n_bin
  real*8 :: dlnr

  ! process commandline arguments
  if (iargc() .ne. 2) then
     write(6,*) 'Usage: extract_sectional_aero_size_num' &
          // ' <netcdf_state_prefix> <output_filename>'
     call exit(2)
  endif
  call getarg(1, in_prefix)
  call getarg(2, out_filename)

  ! process NetCDF files
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
     if (n_time >= MAX_N_TIME) then
        write(0,*) 'ERROR: can only process up to MAX_N_TIME times: ', &
             MAX_N_TIME
        call exit(1)
     end if

     ! read time
     call nc_check(nf90_inq_varid(ncid, "time", varid_time))
     call nc_check(nf90_get_var(ncid, varid_time, time))

     ! read aero_radius dimension
     call nc_check(nf90_inq_dimid(ncid, "aero_radius", dimid_aero_radius))
     call nc_check(nf90_Inquire_Dimension(ncid, dimid_aero_radius, &
          tmp_str, new_n_bin))

     ! allocate aero_dist on first time
     if (i_time == 1) then
        n_bin = new_n_bin
        allocate(aero_dist(n_bin, MAX_N_TIME))
        aero_dist = 0d0
     else
        if (new_n_bin /= n_bin) then
           write(*,*) "ERROR: n_bin changed"
           call exit(1)
        end if
     end if
        
     ! read aero_radius variable
     call nc_check(nf90_inq_varid(ncid, "aero_radius", &
          varid_aero_radius))
     call nc_check(nf90_Inquire_Variable(ncid, varid_aero_radius, &
          tmp_str, xtype, ndims, dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_radius)) then
        write(*,*) "ERROR: unexpected aero_radius dimids"
        call exit(1)
     end if
     allocate(aero_radius(n_bin))
     call nc_check(nf90_get_var(ncid, varid_aero_radius, &
          aero_radius))
     if (i_time == 1) then
        allocate(save_aero_radius(n_bin))
        save_aero_radius = aero_radius
     end if

     ! read aero_radius_widths variable
     call nc_check(nf90_inq_varid(ncid, "aero_radius_widths", &
          varid_aero_radius_widths))
     call nc_check(nf90_Inquire_Variable(ncid, varid_aero_radius_widths, &
          tmp_str, xtype, ndims, dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_radius)) then
        write(*,*) "ERROR: unexpected aero_radius_widths dimids"
        call exit(1)
     end if
     allocate(aero_radius_widths(n_bin))
     call nc_check(nf90_get_var(ncid, varid_aero_radius_widths, &
          aero_radius_widths))

     ! read aero_number_concentration
     call nc_check(nf90_inq_varid(ncid, "aero_number_concentration", &
          varid_aero_number_concentration))
     call nc_check(nf90_Inquire_Variable(ncid, &
          varid_aero_number_concentration, tmp_str, xtype, ndims, &
          dimids, nAtts))
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_radius)) then
        write(*,*) "ERROR: unexpected aero_number_concentration dimids"
        call exit(1)
     end if
     allocate(aero_number_concentration(n_bin))
     call nc_check(nf90_get_var(ncid, varid_aero_number_concentration, &
          aero_number_concentration))
     
     call nc_check(nf90_close(ncid))

     ! compute distribution
     dlnr = aero_radius_widths(1)
     aero_dist(:, i_time) = aero_number_concentration

     deallocate(aero_radius)
     deallocate(aero_radius_widths)
     deallocate(aero_number_concentration)
  end do

  if (n_time == 0) then
     write(*,'(a,a)') 'ERROR: no input file found matching: ', &
          trim(in_filename)
     call exit(1)
  end if

  ! write information
  write(*,*) "Output file array A has:"
  write(*,*) "  A(i, 1) = radius(i) (m)"
  write(*,*) "  A(i, j+1) = number concentration at radius(i) and" &
       // " time(j) (#/m^3)"
  write(*,*) "Radius bins have logarithmic width:"
  write(*,*) "  d(ln(r)) = ln(radius(i+1)/radius(i)) =", dlnr

  ! open output file
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     call exit(1)
  end if

  ! output data
  do i_bin = 1,n_bin
     write(out_unit, '(e30.15e3)', advance='no') &
          save_aero_radius(i_bin)
     do i_time = 1,n_time
        write(out_unit, '(e30.15e3)', advance='no') &
             aero_dist(i_bin, i_time)
     end do
     write(out_unit, '(a)') ''
  end do

  close(out_unit)
  deallocate(save_aero_radius)
  deallocate(aero_dist)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Check return status of NetCDF function calls.
  subroutine nc_check(status)

    !> Status return value.
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write(0,*) nf90_strerror(status)
       call exit(1)
    end if

  end subroutine nc_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program extract_sectional_aero_size_num
