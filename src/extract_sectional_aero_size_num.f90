! Copyright (C) 2009-2010 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The extract_sectional_aero_size_num program.

!> Read NetCDF sectional output files and write out the aerosol number
!> size distributions in text format.
program extract_sectional_aero_size_num

  use netcdf

  integer, parameter :: dp = kind(0.d0)
  integer, parameter :: MAX_N_TIME = 10000
  integer, parameter :: out_unit = 64

  character(len=1000) :: in_prefix, in_filename, out_filename
  integer :: ncid
  integer :: dimid_aero_diam
  integer :: varid_time, varid_aero_number_concentration
  integer :: varid_aero_diam, varid_aero_diam_widths
  integer :: n_bin
  character(len=1000) :: tmp_str
  real(kind=dp) :: time
  real(kind=dp), allocatable :: aero_dist(:,:)
  real(kind=dp), allocatable :: aero_diam(:)
  real(kind=dp), allocatable :: aero_diam_widths(:)
  real(kind=dp), allocatable :: aero_number_concentration(:)
  real(kind=dp), allocatable :: save_aero_diam(:)
  integer :: xtype, ndims, nAtts
  integer, dimension(nf90_max_var_dims) :: dimids
  integer :: ios, i_time, i_spec, i_part, status
  integer :: i_bin, n_time, new_n_bin
  real(kind=dp) :: log_width
  character(len=36) :: uuid, run_uuid

  ! process commandline arguments
  if (command_argument_count() .ne. 2) then
     write(6,*) 'Usage: extract_sectional_aero_size_num' &
          // ' <netcdf_state_prefix> <output_filename>'
     stop 2
  endif
  call get_command_argument(1, in_prefix)
  call get_command_argument(2, out_filename)

  ! process NetCDF files
  i_time = 0
  n_time = 0
  n_bin = 0 ! HACK to shut up gfortran warning
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
        stop 1
     end if

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

     ! read aero_diam dimension
     call nc_check_msg(nf90_inq_dimid(ncid, "aero_diam", &
          dimid_aero_diam), "getting dimension ID for 'aero_diam'")
     call nc_check_msg(nf90_Inquire_Dimension(ncid, dimid_aero_diam, &
          tmp_str, new_n_bin), "inquiring dimension 'aero_diam'")

     ! allocate aero_dist on first time
     if (i_time == 1) then
        n_bin = new_n_bin
        allocate(aero_dist(n_bin, MAX_N_TIME))
        aero_dist = 0d0
     else
        if (new_n_bin /= n_bin) then
           write(*,*) "ERROR: n_bin changed"
           stop 1
        end if
     end if
        
     ! read aero_diam variable
     call nc_check_msg(nf90_inq_varid(ncid, "aero_diam", &
          varid_aero_diam), &
          "getting variable ID for 'aero_diam'")
     call nc_check_msg(nf90_Inquire_Variable(ncid, varid_aero_diam, &
          tmp_str, xtype, ndims, dimids, nAtts), &
          "inquiring variable 'aero_diam'")
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_diam)) then
        write(*,*) "ERROR: unexpected aero_diam dimids"
        stop 1
     end if
     allocate(aero_diam(n_bin))
     call nc_check_msg(nf90_get_var(ncid, varid_aero_diam, &
          aero_diam), "getting variable 'aero_diam'")
     if (i_time == 1) then
        allocate(save_aero_diam(n_bin))
        save_aero_diam = aero_diam
     end if

     ! read aero_diam_widths variable
     call nc_check_msg(nf90_inq_varid(ncid, "aero_diam_widths", &
          varid_aero_diam_widths), &
          "getting variable ID 'aero_diam_widths'")
     call nc_check_msg(nf90_Inquire_Variable(ncid, varid_aero_diam_widths, &
          tmp_str, xtype, ndims, dimids, nAtts), &
          "inquiring variable 'aero_diam_widths'")
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_diam)) then
        write(*,*) "ERROR: unexpected aero_diam_widths dimids"
        stop 1
     end if
     allocate(aero_diam_widths(n_bin))
     call nc_check_msg(nf90_get_var(ncid, varid_aero_diam_widths, &
          aero_diam_widths), "getting variable 'aero_diam_widths'")

     ! read aero_number_concentration
     call nc_check_msg(nf90_inq_varid(ncid, "aero_number_concentration", &
          varid_aero_number_concentration), &
          "getting variable ID for 'aero_number_concentration'")
     call nc_check_msg(nf90_Inquire_Variable(ncid, &
          varid_aero_number_concentration, tmp_str, xtype, ndims, &
          dimids, nAtts), "inquiring variable 'aero_number_concentration'")
     if ((ndims /= 1) &
          .or. (dimids(1) /= dimid_aero_diam)) then
        write(*,*) "ERROR: unexpected aero_number_concentration dimids"
        stop 1
     end if
     allocate(aero_number_concentration(n_bin))
     call nc_check_msg(nf90_get_var(ncid, varid_aero_number_concentration, &
          aero_number_concentration), &
          "getting variable 'aero_number_concentration'")
     
     call nc_check_msg(nf90_close(ncid), &
          "closing file " // trim(in_filename))

     ! compute distribution
     log_width = aero_diam_widths(1)
     aero_dist(:, i_time) = aero_number_concentration

     deallocate(aero_diam)
     deallocate(aero_diam_widths)
     deallocate(aero_number_concentration)
  end do

  if (n_time == 0) then
     write(*,'(a,a)') 'ERROR: no input file found matching: ', &
          trim(in_filename)
     stop 1
  end if

  ! write information
  write(*,*) "Output file array A has:"
  write(*,*) "  A(i, 1) = diameter(i) (m)"
  write(*,*) "  A(i, j+1) = number concentration at diameter(i) and" &
       // " time(j) (#/m^3)"
  write(*,*) "Diameter bins have logarithmic width:"
  write(*,*) "  log_width = ln(diameter(i+1)) - ln(diameter(i)) =", log_width

  ! open output file
  open(unit=out_unit, file=out_filename, status='replace', iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     stop 1
  end if

  ! output data
  do i_bin = 1,n_bin
     write(out_unit, '(e30.15e3)', advance='no') &
          save_aero_diam(i_bin)
     do i_time = 1,n_time
        write(out_unit, '(e30.15e3)', advance='no') &
             aero_dist(i_bin, i_time)
     end do
     write(out_unit, '(a)') ''
  end do

  close(out_unit)
  deallocate(save_aero_diam)
  deallocate(aero_dist)

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
  
end program extract_sectional_aero_size_num
