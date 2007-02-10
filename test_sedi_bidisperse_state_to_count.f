! -*- mode: f90; -*-
! Copyright (C) 2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Process state files to produce number of small particles and volume
! of big particle.

program test_sedi_bidisperse_state_to_count
  
  use mod_environ
  use mod_material
  use mod_state
  use mod_array

  integer, parameter :: out_unit = 33   ! output unit number
  integer, parameter :: state_unit = 34 ! state file unit number
  character(len=*), parameter :: out_name = "counts_sedi_bidisperse_mc.d"
  character(len=*), parameter :: state_prefix = "state_sedi_bidisperse_mc_"
  integer, parameter :: n_time = 600    ! number of state files
  integer, parameter :: time_inc = 10   ! increment for state files
  
  character(len=1000) :: state_name     ! name of state file to read
  integer, allocatable :: MH(:)         ! number of particles per bin
  type(bin_p), allocatable :: VH(:)     ! particle volumes (m^3)
  integer :: n_bin                      ! number of bins
  integer :: n_spec                     ! number of species
  real*8, allocatable :: bin_v(:)       ! volume of particles in bins
  type(environ) :: env                  ! environment state
  real*8 :: time                        ! current time (s)
  real*8 :: dlnr                        ! bin scale factor

  integer i, j, k, n_small
  real*8 v_big

  open(unit=out_unit, file=out_name)

  do i = 0,n_time,time_inc
     write(state_name, "(a,i8.8,a)") state_prefix, i, ".d"
     call read_state_header(state_unit, state_name, n_bin, n_spec)

     if (i == 0) then
        allocate(MH(n_bin))
        allocate(VH(n_bin))
        allocate(bin_v(n_bin))
        call init_array(n_spec, MH, VH)
     else
        call zero_array(n_spec, MH, VH)
     end if
     call read_state_bins(state_unit, state_name, n_bin, bin_v, dlnr)
     call read_state(state_unit, state_name, n_bin, n_spec, MH, VH, env, time)

     ! if there is only one particle, assume it is big
     if (sum(MH) == 1) then
        n_small = 0
        j = 1
     else
        do j = 1,n_bin
           if (MH(j) > 0) then
              n_small = MH(j)
              exit
           end if
        end do
     end if

     v_big = 0d0
     do j = (j + 1),n_bin
        do k = 1,MH(j)
           v_big = v_big + particle_volume(VH(j)%p(k,:))
        end do
     end do

     write(*,'(a8,a14,a14)') &
          't', 'n_small', 'v_big'
     write(*,'(f8.1,e14.5,e14.5)') &
          time, dble(n_small) / env%V_comp / dlnr, v_big / env%V_comp / dlnr
     write(out_unit,'(e20.10,e20.10,e20.10)') &
          time, dble(n_small) / env%V_comp / dlnr, v_big / env%V_comp / dlnr
  end do

  close(out_unit)

end program test_sedi_bidisperse_state_to_count
