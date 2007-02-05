! -*- mode: f90; -*-
! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Generic driver for Monte Carlo simulations using input specification.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program run_mc

  use mod_bin
  use mod_array
  use mod_array_hybrid
  use mod_init_dist
  use mod_condensation
  use mod_kernel_sedi
  use mod_kernel_golovin
  use mod_material
  use mod_environ
  use mod_mc_fix_hybrid

  integer, parameter :: in_unit = 32
  integer, parameter :: MAX_DIST_ARGS = 10

  character(len=300) :: in_name
  integer :: MM, M, M_new, i_loop
  character(len=300) :: out_file_name
  integer, allocatable :: MH(:), bin_n(:)
  type(bin_p), allocatable ::  VH(:)
  real*8, allocatable :: bin_v(:), n_den(:), bin_g(:), bin_gs(:,:), V(:,:)
  real*8 :: dlnr
  integer :: line_num, i, ios

  character(len=300) :: output_name ! name of output files
  integer :: n_loop             ! number of Monte Carlo loops
  real*8 :: N_0                 ! particle concentration (#/m^3)
  character(len=100) :: kernel_name ! coagulation kernel name

  real*8 :: t_max               ! total simulation time (s)
  real*8 :: del_t               ! timestep (s)
  real*8 :: t_print             ! output interval (0 disables) (s)
  real*8 :: t_state             ! state output interval (0 disables) (s)
  real*8 :: t_progress          ! progress printing interval (0 disables) (s)

  integer :: n_spec             ! number of species
  type(material) :: mat         ! material data
  integer :: n_temps            ! number of temperatures
  type(environ) :: env          ! environment data

  integer :: n_init_dist        ! number of initial distributions
  integer, allocatable :: dist_n_part(:) ! distribution particle numbers
  character(len=100), allocatable :: dist_types(:) ! distribution names
  real*8, allocatable :: dist_args(:,:) ! distribution arguments
  real*8, allocatable :: dist_vol_frac(:,:) ! distribution composition

  integer :: n_bin              ! number of bins
  real*8 :: v_min               ! volume of smallest bin (m^3)
  integer :: scal               ! scale factor (integer)

  integer :: rand_init          ! random initialization (0 to 
  logical :: do_coagulation     ! whether to do coagulation (yes/no)
  logical :: do_condensation    ! whether to do condensation (yes/no)
  logical :: do_restart         ! whether to restart from stored state (yes/no)
      character(len=300) :: restart_name ! filename to restart from
  
  ! check there is exactly one commandline argument
  if (iargc() .ne. 1) then
     write(6,*) 'Usage: run_mc <filename.d>'
     call exit(2)
  endif
  
  ! get and check first commandline argument (must be "filename.spec")
  call getarg(1, in_name)
  i = len_trim(in_name)
  if (in_name((i-4):i) /= '.spec') then
     write(6,*) 'ERROR: input filename must end in .spec'
     call exit(2)
  end if

  ! open the input file
  open(unit=in_unit, status='old', file=in_name, iostat=ios)
  if (ios /= 0) then
     write(*,*) 'ERROR: unable to open file ', in_name, ': ', ios
     call exit(1)
  end if

  line_num = 0

  call read_string(in_name, in_unit, line_num, 'output_name', output_name)
  call read_integer(in_name, in_unit, line_num, 'n_loop', n_loop)
  call read_real(in_name, in_unit, line_num, 'N_0', N_0)
  call read_string(in_name, in_unit, line_num, 'kernel', kernel_name)

  call read_real(in_name, in_unit, line_num, 't_max', t_max)
  call read_real(in_name, in_unit, line_num, 'del_t', del_t)
  call read_real(in_name, in_unit, line_num, 't_print', t_print)
  call read_real(in_name, in_unit, line_num, 't_state', t_state)
  call read_real(in_name, in_unit, line_num, 't_progress', t_progress)

  call read_integer(in_name, in_unit, line_num, 'n_spec', n_spec)
  call allocate_material(mat, n_spec)
  call read_integer(in_name, in_unit, line_num, 'i_water', mat%i_water)
  call read_real_array(in_name, in_unit, line_num, n_spec, 'rho', mat%rho)
  call read_integer_array(in_name, in_unit, line_num, n_spec, 'nu', mat%nu)
  call read_real_array(in_name, in_unit, line_num, n_spec, 'eps', mat%eps)
  call read_real_array(in_name, in_unit, line_num, n_spec, 'M_w', mat%M_w)

  call read_integer(in_name, in_unit, line_num, 'n_temps', n_temps)
  call allocate_environ_temps(env, n_temps)
  call read_real_array(in_name, in_unit, line_num, n_temps, &
       'temp_times', env%temp_times)
  call read_real_array(in_name, in_unit, line_num, n_temps, &
       'temps', env%temps)
  call read_real(in_name, in_unit, line_num, 'RH', env%RH)
  call read_real(in_name, in_unit, line_num, 'pressure', env%p)

  call read_integer(in_name, in_unit, line_num, 'n_init_dist', n_init_dist)
  allocate(dist_n_part(n_init_dist))
  allocate(dist_vol_frac(n_init_dist, n_spec))
  allocate(dist_types(n_init_dist))
  allocate(dist_args(n_init_dist, MAX_DIST_ARGS))

  do i = 1,n_init_dist
     call read_integer(in_name, in_unit, line_num, 'n_p', dist_n_part(i))
     call read_real_array(in_name, in_unit, line_num, n_spec, &
          'vol_frac', dist_vol_frac(i,:))
     call read_string(in_name, in_unit, line_num, 'dist_type', dist_types(i))
     if (trim(dist_types(i)) == 'log_normal') then
        call read_real(in_name, in_unit, line_num, 'dist_mean_diam', &
             dist_args(i,1))
        call read_real(in_name, in_unit, line_num, 'dist_std_dev', &
             dist_args(i,2))
     elseif (trim(dist_types(i)) == 'exp') then
        call read_real(in_name, in_unit, line_num, 'dist_mean_vol', &
             dist_args(i,1))
     else
        write(*,'(a,a,a,a,a,i3)') 'ERROR: Unknown distribution type ', &
             trim(dist_types(i)), ' in file ', trim(in_name), &
             ' at line ', line_num
        call exit(1)
     end if
  end do

  call read_integer(in_name, in_unit, line_num, 'n_bin', n_bin)
  call read_real(in_name, in_unit, line_num, 'v_min', v_min)
  call read_integer(in_name, in_unit, line_num, 'scal', scal)

  call read_integer(in_name, in_unit, line_num, 'rand_init', rand_init)
  call read_logical(in_name, in_unit, line_num, &
       'do_coagulation', do_coagulation)
  call read_logical(in_name, in_unit, line_num, &
       'do_condensation', do_condensation)
  call read_logical(in_name, in_unit, line_num, &
       'do_restart', do_restart)
  call read_string(in_name, in_unit, line_num, 'restart_name', restart_name)

  close(unit=in_unit)


  MM = sum(dist_n_part)
  allocate(MH(n_bin), VH(n_bin), V(MM,n_spec), bin_v(n_bin), n_den(n_bin))
  allocate(bin_g(n_bin), bin_gs(n_bin,n_spec), bin_n(n_bin))
  call init_hybrid(n_spec, MH, VH)

  write(out_file_name, '(a,a,a)') 'out_', trim(output_name), '.d'
  open(30,file=out_file_name)
  call print_header(n_loop, n_bin, n_spec, nint(t_max / t_print) + 1)

  if (rand_init /= 0) then
     call srand(rand_init)
  else
     call srand(time())
  end if

  do i_loop = 1,n_loop
     
     call make_bin_grid(n_bin, scal, v_min, bin_v, dlnr)
     call zero_v(MM, n_spec, V)

     M = 1
     do i = 1,n_init_dist
        if (trim(dist_types(i)) == 'log_normal') then
           call init_log_normal(dist_args(i,1), dist_args(i,2), n_bin, &
                bin_v, n_den)
        elseif (trim(dist_types(i)) == 'exp') then
           call init_exp(dist_args(i,1), n_bin, bin_v, n_den)
        else
           write(*,*) 'ERROR: unknown distribution type: ', trim(dist_types(i))
           call exit(1)
        end if
        call dist_to_n(dist_n_part(i), dlnr, n_bin, bin_v, n_den, bin_n)
        call compute_volumes(n_bin, n_spec, dist_vol_frac(i,:), &
             MM, M, M + dist_n_part(i) - 1, bin_n, bin_v, dlnr, V, M_new)
        M = M + M_new
     end do

     env%V_comp = dble(M) / N_0

     ! call equlibriate_particle for each particle in V
     if (do_condensation) then
        do i = 1,M
           call equilibriate_particle(n_spec, V(i,:), env, mat)
        end do
     end if

     if (trim(kernel_name) == 'sedi') then
        call mc_fix_hybrid(MM, M, n_spec, V, n_bin, MH, VH, &
             bin_v, bin_g, bin_gs, bin_n, dlnr, &
             kernel_sedi, t_max, t_print, t_state, t_progress, del_t, &
             do_coagulation, do_condensation, do_restart, restart_name, &
             i_loop, env, mat)
     elseif (trim(kernel_name) == 'golovin') then
        call mc_fix_hybrid(MM, M, n_spec, V, n_bin, MH, VH, &
             bin_v, bin_g, bin_gs, bin_n, dlnr, &
             kernel_golovin, t_max, t_print, t_state, t_progress, del_t, &
             do_coagulation, do_condensation, do_restart, restart_name, &
             i_loop, env, mat)
     else
        write(*,*) 'ERROR: Unknown kernel type; ', trim(kernel_name)
        call exit(1)
     end if

  enddo

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_string(in_name, in_unit, line_num, name, var)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    integer, intent(inout) :: line_num      ! current line number
    character(len=*), intent(in) :: name    ! name that should start line
    character(len=*), intent(out) :: var    ! variable to store data

    character(len=300) :: line, rest
    integer :: ios

    call read_next_data_line(in_name, in_unit, line, line_num)
    call check_name(line, name, in_name, line_num, rest)
    read(rest, '(a)', iostat=ios) var
    call check_read_iostat(ios, in_name, line_num, 'string')
! DEBUG
    write(*,*) 'string: ', trim(name), ' = ', trim(var)
! DEBUG
    
  end subroutine read_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_integer(in_name, in_unit, line_num, name, var)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    integer, intent(inout) :: line_num      ! current line number
    character(len=*), intent(in) :: name    ! name that should start line
    integer, intent(out) :: var             ! variable to store data

    character(len=300) :: line, rest
    integer :: ios

    call read_next_data_line(in_name, in_unit, line, line_num)
    call check_name(line, name, in_name, line_num, rest)
    read(rest, '(i20)', iostat=ios) var
    call check_read_iostat(ios, in_name, line_num, 'integer')
! DEBUG
    write(*,*) 'integer: ', trim(name), ' = ', var
! DEBUG

  end subroutine read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_real(in_name, in_unit, line_num, name, var)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    integer, intent(inout) :: line_num      ! current line number
    character(len=*), intent(in) :: name    ! name that should start line
    real*8, intent(out) :: var              ! variable to store data

    character(len=300) :: line, rest
    integer :: ios

    call read_next_data_line(in_name, in_unit, line, line_num)
    call check_name(line, name, in_name, line_num, rest)
    read(rest, '(f20.0)', iostat=ios) var
    call check_read_iostat(ios, in_name, line_num, 'real')
! DEBUG
    write(*,*) 'real: ', trim(name), ' = ', var
! DEBUG

  end subroutine read_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_logical(in_name, in_unit, line_num, name, var)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    integer, intent(inout) :: line_num      ! current line number
    character(len=*), intent(in) :: name    ! name that should start line
    logical, intent(out) :: var             ! variable to store data

    character(len=300) :: line, rest
    character(len=20) :: str_var
    integer :: ios

    call read_next_data_line(in_name, in_unit, line, line_num)
    call check_name(line, name, in_name, line_num, rest)
    read(rest, '(a)', iostat=ios) str_var
    call check_read_iostat(ios, in_name, line_num, 'string')
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
       write(*,'(a,a,a,i3)') 'ERROR: reading logical from file ', in_name, &
            ' at line ', line_num
       call exit(1)
    end if
! DEBUG
    write(*,*) 'logical: ', trim(name), ' = ', var
! DEBUG

  end subroutine read_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_integer_array(in_name, in_unit, line_num, len, name, var)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    integer, intent(inout) :: line_num      ! current line number
    integer, intent(in) :: len              ! array length
    character(len=*), intent(in) :: name    ! name that should start line
    integer, intent(out) :: var(len)        ! variable to store data

    character(len=300) :: line, rest, temp
    integer :: ios, ind, i

    call read_next_data_line(in_name, in_unit, line, line_num)
    call check_name(line, name, in_name, line_num, rest)
    do i = 1,len
       ind = index(rest, ' ')
       temp = rest(1:ind)
       if (len_trim(temp) == 0) then
          write(*,'(a,a,a,i3)') 'ERROR: insufficient data in file ', &
               trim(in_name), ' on line ', line_num
          call exit(1)
       end if
       rest = rest((ind+1):len_trim(rest))
       read(temp, '(i20)', iostat=ios) var(i)
       call check_read_iostat(ios, in_name, line_num, 'integer array')
    end do
    if (len_trim(rest) > 0) then
       write(*,'(a,a,a,i3)') 'ERROR: too much data in file ', &
            trim(in_name), ' on line ', line_num
       call exit(1)
    end if
! DEBUG
    write(*,*) 'integer array: ', trim(name), ' = ', var
! DEBUG
    
  end subroutine read_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_real_array(in_name, in_unit, line_num, len, name, var)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    integer, intent(inout) :: line_num      ! current line number
    integer, intent(in) :: len              ! array length
    character(len=*), intent(in) :: name    ! name that should start line
    real*8, intent(out) :: var(len)        ! variable to store data

    character(len=300) :: line, rest, temp
    integer :: ios, ind, i

    call read_next_data_line(in_name, in_unit, line, line_num)
    call check_name(line, name, in_name, line_num, rest)
    do i = 1,len
       ind = index(rest, ' ')
       temp = rest(1:ind)
       if (len_trim(temp) == 0) then
          write(*,'(a,a,a,i3)') 'ERROR: insufficient data in file ', &
               trim(in_name), ' on line ', line_num
          call exit(1)
       end if
       rest = rest((ind+1):len_trim(rest))
       read(temp, '(f20.0)', iostat=ios) var(i)
       call check_read_iostat(ios, in_name, line_num, 'real array')
    end do
    if (len_trim(rest) > 0) then
       write(*,'(a,a,a,i3)') 'ERROR: too much data in file ', &
            trim(in_name), ' on line ', line_num
       call exit(1)
    end if
! DEBUG
    write(*,*) 'real array: ', trim(name), ' = ', var
! DEBUG

  end subroutine read_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_read_iostat(ios, in_name, line_num, type)

    integer, intent(in) :: ios              ! iostat result
    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: line_num         ! current line number
    character(len=*), intent(in) :: type    ! type being read during error

    if (ios /= 0) then
       write(*,'(a,a,a,a,i3,a,i4)') 'ERROR: reading ', trim(type), &
            ' from file ', trim(in_name), ' at line ', line_num, &
            ': IOSTAT = ', ios
       call exit(1)
    end if

  end subroutine check_read_iostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_next_data_line(in_name, in_unit, line, line_num)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    character(len=*), intent(out) :: line   ! complete line read from in_unit
    integer, intent(inout) :: line_num      ! current line number

    logical :: done

    done = .false.
    do while (.not. done)
       call read_line(in_name, in_unit, line, line_num)
       call strip_comment(line)
       if (len_trim(line) > 0) then
          done = .true.
       end if
    end do

  end subroutine read_next_data_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_line(in_name, in_unit, line, line_num)

    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: in_unit          ! unit number to read from
    character(len=*), intent(out) :: line   ! complete line read from in_unit
    integer, intent(inout) :: line_num      ! current line number

    integer ios

    read(unit=in_unit, fmt='(a)', end=100, iostat=ios) line
    line_num = line_num + 1
    if (ios /= 0) then
       write(*,*) 'ERROR: reading from ', trim(in_name), &
            ' at line ', line_num, ': IOSTAT = ', ios
       call exit(1)
    end if

    return

100 write(*,*) 'ERROR: end of file ', trim(in_name), &
         ' unexpected at line ', line_num
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
  
  subroutine check_name(line, name, in_name, line_num, rest)
    
    character(len=*), intent(in) :: line   ! complete input line
    character(len=*), intent(in) :: name   ! name that should start line
    character(len=*), intent(in) :: in_name ! filename of input file
    integer, intent(in) :: line_num        ! number of line in input file
    character(len=*), intent(out) :: rest  ! remainder of line without name
    
    integer name_len

    name_len = len_trim(name)
    if (index(line, name(1:name_len)) /= 1) then
       write(*,'(a,i3,a,a,a,a)') 'ERROR: line ', line_num, &
            ' of input file ', trim(in_name), &
            ' must begin with: ', trim(name)
       call exit(1)
    end if
    rest = line((name_len+2):len_trim(line))
    
  end subroutine check_name
  
end program run_mc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
