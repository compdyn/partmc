! Copyright (C) 2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

program partmc

  use pmc_condense
  use pmc_util

  integer, parameter :: unit_D = 20
  integer, parameter :: unit_H = 21
  integer, parameter :: unit_delta_star = 22
  integer, parameter :: unit_Ddot = 23
  integer, parameter :: unit_Hdot_i = 24
  integer, parameter :: unit_Hdot_env = 25

  integer, parameter :: i_diam_start = 10

  type(condense_rates_inputs_t) :: inputs
  type(condense_rates_outputs_t) :: outputs
  real(kind=dp) :: diams(81), rhs(31)
  integer :: i_diam, i_rh

  inputs%T = 288d0
  inputs%Tdot = 0d0
  inputs%p = 1d5
  inputs%V_comp = 1d-6
  inputs%D_dry = 1d-6
  inputs%kappa = 0.5d0

  call linspace(0.98d0, 1.04d0, rhs)
  call logspace(inputs%D_dry, 1d-4, diams)

  open(unit=unit_D, file="D.txt", status='replace')
  open(unit=unit_H, file="H.txt", status='replace')
  open(unit=unit_delta_star, file="delta_star.txt", status='replace')
  open(unit=unit_Ddot, file="Ddot.txt", status='replace')
  open(unit=unit_Hdot_i, file="Hdot_i.txt", status='replace')
  open(unit=unit_Hdot_env, file="Hdot_env.txt", status='replace')

  do i_diam = i_diam_start,size(diams)
     do i_rh = 1,size(rhs)
        inputs%D = diams(i_diam)
        inputs%H = rhs(i_rh)
        if (i_rh == 1) then
           write(unit_D, '(e35.20e3)') inputs%D
        end if
        if (i_diam == i_diam_start) then
           write(unit_H, '(e35.20e3)') inputs%H
        end if
        call condense_rates(inputs, outputs)
        write(unit_delta_star, '(e35.20e3)', advance='no') outputs%delta_star
        write(unit_Ddot, '(e35.20e3)', advance='no') outputs%Ddot
        write(unit_Hdot_i, '(e35.20e3)', advance='no') outputs%Hdot_i
        write(unit_Hdot_env, '(e35.20e3)', advance='no') outputs%Hdot_env
     end do
     write(unit_delta_star, *) ''
     write(unit_Ddot, *) ''
     write(unit_Hdot_i, *) ''
     write(unit_Hdot_env, *) ''
  end do

  close(unit_D)
  close(unit_H)
  close(unit_delta_star)
  close(unit_Ddot)
  close(unit_Hdot_i)
  close(unit_Hdot_env)

end program partmc
