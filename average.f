! -*- mode: f90; -*-
! Copyright (C) 2005,2006 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

program average

  character*100 dataname
  character*100 format_string_in
  character*100 format_string_out
  integer n
  integer m
  real*8 array(1000,100)
  integer i

  n = 160
  m = 4
  dataname = "state_1300_moments.d"
  format_string_in = "(e20.0,e20.10,e20.0,e20.10)"
  format_string_out = "(e20.10,e20.10,e20.10,e20.10)"
  call do_average(dataname, format_string_in, format_string_out, n, m)

  n = 160
  m = 8
  dataname = "state_1300_moments_comp.d"
  format_string_in = "(e20.0,e20.10,e20.0,e20.0,e20.0,e20.10,e20.10,e20.10)"
  format_string_out = "(e20.10,e20.10,e20.10,e20.10,e20.10,e20.10,e20.10,e20.10)"
  call do_average(dataname, format_string_in, format_string_out, n, m)

  n = 20
  m = 2
  dataname = "state_1300_composition.d"
  format_string_in = "(e20.0,e20.0)"
  format_string_out = "(e20.10,e20.10)"
  call do_average(dataname, format_string_in, format_string_out, n, m)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine do_average(dataname, format_string_in, format_string_out, &
       n, m)

    character*100, intent(in) :: dataname
    character*100, intent(in) :: format_string_in
    character*100, intent(in) :: format_string_out
    integer, intent(in) :: n, m

    real*8 array(n,m)

    call calc_average(dataname, format_string_in, n, m, array)
    call write_data(dataname, format_string_out, n, m, array)
    
  end subroutine do_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine write_data(dataname, format_string, n, m, array)

    character*100, intent(in) :: dataname
    character*100, intent(in) :: format_string
    integer, intent(in) :: n, m
    real*8, intent(out) :: array(n,m)

    character*100 :: filename
    integer i

    write(filename,'(a3,a30)') 'av/', dataname
    open(30,file=filename)
    
    do i = 1,n
       write(30,format_string) array(i,1:m)
    end do

    close(unit=30)

  end subroutine write_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_average(dataname, format_string, n, m, array)

    character*100, intent(in) :: dataname
    character*100, intent(in) :: format_string
    integer, intent(in) :: n, m
    real*8, intent(out) :: array(n,m)

    integer, parameter :: n_files = 100

    character*100 :: filename
    real*8 :: new_array(n,m)
    integer r

    array = 0d0
    do r = 0,(n_files - 1)
       write(filename,'(a1,i2.2,a1,a30)') 'r', r, '/', dataname
       call read_data(filename, format_string, n, m, new_array)
       array = array + new_array
    end do

    array = array / dble(n_files)

  end subroutine calc_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_data(filename, format_string, n, m, array)

    character*100, intent(in) :: filename
    character*100, intent(in) :: format_string
    integer, intent(in) :: n, m
    real*8, intent(out) :: array(n,m)

    open(30,file=filename)
    
    do i = 1,n
       read(30,format_string) array(i,1:m)
    end do

    close(unit=30)

  end subroutine read_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program average
