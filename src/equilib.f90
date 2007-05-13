! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Driver to compute equilibrium radius.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program equilib

  use mod_read_spec
  use mod_material
  use mod_environ
  use mod_util
  use mod_condensation

  character(len=*), parameter :: aero_file = "aerosol.dat"

  type(spec_file) :: aero_spec
  type(environ) :: env
  type(material) :: mat
  real*8, allocatable :: V(:)
  character(len=300) :: tmp_str
  integer :: i

  ! read aerosol data
  call open_spec(aero_spec, aero_file)
  call read_material_from_file(aero_spec, mat)
  call close_spec(aero_spec)

  ! check command line arguments
  if (iargc() .ne. (1 + mat%n_spec)) then
     write(6,*) 'Usage: partmc <temp> <RH> <species volumes>'
     write(6,'(a,i3,a)') 'Need ', (mat%n_spec - 1), ' species volumes'
     call exit(2)
  end if
  
  ! get and check commandline arguments
  call getarg(1, tmp_str)
  env%T = string_to_real(tmp_str)

  call getarg(2, tmp_str)
  env%RH = string_to_real(tmp_str)

  allocate(V(mat%n_spec))
  do i = 2,mat%n_spec
     call getarg(i + 1, tmp_str)
     V(i) = string_to_real(tmp_str)
  end do

  call equilibriate_particle(mat%n_spec, V, env, mat)

  do i = 1,mat%n_spec
     write(*,*) mat%name(i), V(i)
  end do

  write(*,*) ''
  write(*,*) 'Equilibrium water volume: ', V(mat%i_water)
  write(*,*) 'Dry volume: ', sum(V) - V(mat%i_water)
  write(*,*) 'West volume: ', sum(V)

end program equilib
