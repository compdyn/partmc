! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Driver to compute equilibrium radius.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program equilib

  use mod_inout
  use mod_aero_data
  use mod_environ
  use mod_util
  use mod_condensation

  character(len=*), parameter :: aero_filename = "aerosol.dat"

  type(inout_file_t) :: file
  type(environ) :: env
  type(aero_data_t) :: aero_data
  real*8, allocatable :: V(:)
  character(len=300) :: tmp_str
  integer :: i

  open(unit=40, file="temp.dat")
  write(40, *) "Realling running!"
  close(40)

  ! read aerosol data
  call inout_open_read(aero_filename, file)
  call inout_read_aero_data(file, aero_data)
  call inout_close(file)

  ! check command line arguments
  if (iargc() .ne. (1 + aero_data%n_spec)) then
     write(6,*) 'Usage: partmc <temp> <RH> <species volumes>'
     write(6,'(a,i3,a)') 'Need ', (aero_data%n_spec - 1), ' species volumes'
     call exit(2)
  end if
  
  ! get and check commandline arguments
  call getarg(1, tmp_str)
  env%T = string_to_real(tmp_str)

  call getarg(2, tmp_str)
  env%RH = string_to_real(tmp_str)

  allocate(V(aero_data%n_spec))
  V = 0d0
  do i = 2,aero_data%n_spec
     call getarg(i + 1, tmp_str)
     V(i) = string_to_real(tmp_str)
  end do

  call equilibriate_particle(aero_data%n_spec, V, env, aero_data)

  write(*,*) 'temp = ', env%T, ' K'
  write(*,*) 'RH = ', env%RH
  write(*,*) ''

  do i = 1,aero_data%n_spec
     write(*,*) aero_data%name(i), V(i), ' m^{-3}'
  end do

  write(*,*) ''
  write(*,*) 'Equilibrium water volume: ', V(aero_data%i_water), ' m^{-3}'
  write(*,*) 'Dry volume: ', sum(V) - V(aero_data%i_water), ' m^{-3}'
  write(*,*) 'Wet volume: ', sum(V), ' m^{-3}'

end program equilib
