! -*- mode: f90; -*-
! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Routines to read data out of .spec files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_read_spec

  type spec_file
     character(len=300) :: name  ! filename
     integer :: unit             ! attached unit
     integer :: line_num         ! current line number
  end type spec_file

  logical, parameter :: DEBUG_OUTPUT = .false.

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_spec(spec, filename, unit)

    type(spec_file), intent(out) :: spec     ! spec file
    character(len=*), intent(in) :: filename ! name of file to open
    integer, intent(in) :: unit              ! unit number to use

    integer ios

    spec%name = trim(filename)
    spec%unit = unit
    open(unit=spec%unit, status='old', file=spec%name, iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: unable to open file ', spec%name, ': ', ios
       call exit(1)
    end if

  end subroutine open_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine close_spec(spec)

    type(spec_file), intent(in) :: spec     ! spec file

    close(spec%unit)

  end subroutine close_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_init_dist(spec, dist_type, dist_args)

    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(out) :: dist_type ! type of init distribution
    real*8, intent(out) :: dist_args(:)     ! distribution parameters

    call read_string(spec, 'dist_type', dist_type)
    if (trim(dist_type) == 'log_normal') then
       call read_real(spec, 'dist_mean_diam', dist_args(1))
       call read_real(spec, 'dist_std_dev', dist_args(2))
    elseif (trim(dist_type) == 'exp') then
       call read_real(spec, 'dist_mean_vol', dist_args(1))
    else
       write(*,'(a,a,a,a,a,i3)') 'ERROR: Unknown distribution type ', &
            trim(dist_type), ' in file ', trim(spec%name), &
            ' at line ', spec%line_num
       call exit(1)
    end if

  end subroutine read_init_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_material(spec, mat)

    use mod_material

    type(spec_file), intent(inout) :: spec  ! spec file
    type(material), intent(out) :: mat      ! material data

    integer :: n_spec

    call read_integer(spec, 'n_spec', n_spec)
    call allocate_material(mat, n_spec)
    call read_integer(spec, 'i_water', mat%i_water)
    call read_real_array(spec, n_spec, 'rho', mat%rho)
    call read_integer_array(spec, n_spec, 'nu', mat%nu)
    call read_real_array(spec, n_spec, 'eps', mat%eps)
    call read_real_array(spec, n_spec, 'M_w', mat%M_w)

  end subroutine read_material

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_environ(spec, env)

    use mod_environ

    type(spec_file), intent(inout) :: spec  ! spec file
    type(environ), intent(out) :: env       ! environment data

    integer :: n_temps

    call read_integer(spec, 'n_temps', n_temps)
    call allocate_environ_temps(env, n_temps)
    call read_real_array(spec, n_temps, 'temp_times', env%temp_times)
    call read_real_array(spec, n_temps, 'temps', env%temps)
    call read_real(spec, 'RH', env%RH)
    call read_real(spec, 'pressure', env%p)
    
  end subroutine read_environ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_string(spec, name, var)

    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(in) :: name    ! name that should start line
    character(len=*), intent(out) :: var    ! variable to store data

    character(len=300) :: line, rest
    integer :: ios

    call read_next_data_line(spec, line)
    call check_name(spec, line, name, rest)
    read(rest, '(a)', iostat=ios) var
    call check_read_iostat(spec, ios, 'string')
    if (DEBUG_OUTPUT) then
       write(*,*) 'string: ', trim(name), ' = ', trim(var)
    end if
    
  end subroutine read_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_integer(spec, name, var)

    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(in) :: name    ! name that should start line
    integer, intent(out) :: var             ! variable to store data

    character(len=300) :: line, rest
    integer :: ios

    call read_next_data_line(spec, line)
    call check_name(spec, line, name, rest)
    read(rest, '(i20)', iostat=ios) var
    call check_read_iostat(spec, ios, 'integer')
    if (DEBUG_OUTPUT) then
       write(*,*) 'integer: ', trim(name), ' = ', var
    end if

  end subroutine read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_real(spec, name, var)

    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(in) :: name    ! name that should start line
    real*8, intent(out) :: var              ! variable to store data

    character(len=300) :: line, rest
    integer :: ios

    call read_next_data_line(spec, line)
    call check_name(spec, line, name, rest)
    read(rest, '(f20.0)', iostat=ios) var
    call check_read_iostat(spec, ios, 'real')
    if (DEBUG_OUTPUT) then
       write(*,*) 'real: ', trim(name), ' = ', var
    end if

  end subroutine read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_logical(spec, name, var)

    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(in) :: name    ! name that should start line
    logical, intent(out) :: var             ! variable to store data

    character(len=300) :: line, rest
    character(len=20) :: str_var
    integer :: ios

    call read_next_data_line(spec, line)
    call check_name(spec, line, name, rest)
    read(rest, '(a)', iostat=ios) str_var
    call check_read_iostat(spec, ios, 'string')
    if ((trim(str_var) == 'yes') &
         .or. (trim(str_var) == 'y') &
         .or. (trim(str_var) == 'true') &
         .or. (trim(str_var) == 't') &
         .or. (trim(str_var) == '1')) then
       var = .true.
    elseif ((trim(str_var) == 'no') &
         .or. (trim(str_var) == 'n') &
         .or. (trim(str_var) == 'false') &
         .or. (trim(str_var) == 'f') &
         .or. (trim(str_var) == '0')) then
       var = .false.
    else
       write(*,'(a,a,a,i3)') 'ERROR: reading logical from file ', &
            spec%name, ' at line ', spec%line_num
       call exit(1)
    end if
    if (DEBUG_OUTPUT) then
       write(*,*) 'logical: ', trim(name), ' = ', var
    end if

  end subroutine read_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_integer_array(spec, len, name, var)

    type(spec_file), intent(inout) :: spec  ! spec file
    integer, intent(in) :: len              ! array length
    character(len=*), intent(in) :: name    ! name that should start line
    integer, intent(out) :: var(len)        ! variable to store data

    character(len=300) :: line, rest, temp
    integer :: ios, ind, i

    call read_next_data_line(spec, line)
    call check_name(spec, line, name, rest)
    do i = 1,len
       ind = index(rest, ' ')
       temp = rest(1:ind)
       if (len_trim(temp) == 0) then
          write(*,'(a,a,a,i3)') 'ERROR: insufficient data in file ', &
               trim(spec%name), ' on line ', spec%line_num
          call exit(1)
       end if
       rest = rest((ind+1):len_trim(rest))
       read(temp, '(i20)', iostat=ios) var(i)
       call check_read_iostat(spec, ios, 'integer array')
    end do
    if (len_trim(rest) > 0) then
       write(*,'(a,a,a,i3)') 'ERROR: too much data in file ', &
            trim(spec%name), ' on line ', spec%line_num
       call exit(1)
    end if
    if (DEBUG_OUTPUT) then
       write(*,*) 'integer array: ', trim(name), ' = ', var
    end if
    
  end subroutine read_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_real_array(spec, len, name, var)

    type(spec_file), intent(inout) :: spec  ! spec file
    integer, intent(in) :: len              ! array length
    character(len=*), intent(in) :: name    ! name that should start line
    real*8, intent(out) :: var(len)         ! variable to store data

    character(len=300) :: line, rest, temp
    integer :: ios, ind, i

    call read_next_data_line(spec, line)
    call check_name(spec, line, name, rest)
    do i = 1,len
       ind = index(rest, ' ')
       temp = rest(1:ind)
       if (len_trim(temp) == 0) then
          write(*,'(a,a,a,i3)') 'ERROR: insufficient data in file ', &
               trim(spec%name), ' on line ', spec%line_num
          call exit(1)
       end if
       rest = rest((ind+1):len_trim(rest))
       read(temp, '(f20.0)', iostat=ios) var(i)
       call check_read_iostat(spec, ios, 'real array')
    end do
    if (len_trim(rest) > 0) then
       write(*,'(a,a,a,i3)') 'ERROR: too much data in file ', &
            trim(spec%name), ' on line ', spec%line_num
       call exit(1)
    end if
    if (DEBUG_OUTPUT) then
       write(*,*) 'real array: ', trim(name), ' = ', var
    end if

  end subroutine read_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_read_iostat(spec, ios, type)

    type(spec_file), intent(inout) :: spec  ! spec file
    integer, intent(in) :: ios              ! iostat result
    character(len=*), intent(in) :: type    ! type being read during error

    if (ios /= 0) then
       write(*,'(a,a,a,a,i3,a,i4)') 'ERROR: reading ', trim(type), &
            ' from file ', trim(spec%name), ' at line ', spec%line_num, &
            ': IOSTAT = ', ios
       call exit(1)
    end if

  end subroutine check_read_iostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_next_data_line(spec, line)

    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(out) :: line   ! complete line read

    logical :: done

    done = .false.
    do while (.not. done)
       call read_line(spec, line)
       call strip_comment(line)
       if (len_trim(line) > 0) then
          done = .true.
       end if
    end do

  end subroutine read_next_data_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_line(spec, line)

    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(out) :: line   ! complete line read

    integer ios

    read(unit=spec%unit, fmt='(a)', end=100, iostat=ios) line
    spec%line_num = spec%line_num + 1
    if (ios /= 0) then
       write(*,*) 'ERROR: reading from ', trim(spec%name), &
            ' at line ', spec%line_num, ': IOSTAT = ', ios
       call exit(1)
    end if

    return

100 write(*,*) 'ERROR: end of file ', trim(spec%name), &
         ' unexpected at line ', spec%line_num
    call exit(1)
    
  end subroutine read_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine strip_comment(line)
    
    character(len=*), intent(inout) :: line ! complete input line
    
    integer hash_index
    
    hash_index = index(line, '#')
    if (hash_index > 0) then
       line = line(1:(hash_index - 1))
    end if
    
  end subroutine strip_comment
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine check_name(spec, line, name, rest)
    
    type(spec_file), intent(inout) :: spec  ! spec file
    character(len=*), intent(in) :: line   ! complete input line
    character(len=*), intent(in) :: name   ! name that should start line
    character(len=*), intent(out) :: rest  ! remainder of line without name
    
    integer name_len

    name_len = len_trim(name)
    if (index(line, name(1:name_len)) /= 1) then
       write(*,'(a,i3,a,a,a,a)') 'ERROR: line ', spec%line_num, &
            ' of input file ', trim(spec%name), &
            ' must begin with: ', trim(name)
       call exit(1)
    end if
    rest = line((name_len+2):len_trim(line))
    
  end subroutine check_name
  
end module mod_read_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
