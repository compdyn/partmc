! Copyright (C) 2010, 2011 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Compute the nucleation of SO4 aerosol particles from gas-phase H2SO4
! using a continuous ODE model.
!
! This is hard-coded to certain values and is only designed for use in
! the automated tests of PartMC.
!
! \f[ A = {\rm aerosol number conc} (# / m^3) \f]
! \f[ G = {\rm gas concentration} ({\rm molecules} / m^3) \f]
!
! \f[ \dot{A} = K G^2 \f]
! \f[ G + A {\rm nucleate-vol} {\rm so4-molec-dens} = G_0 \f]
!
! \f[ \dot{G} + \dot{A} {\rm nucleate-vol} {\rm so4-molec-dens} = 0 \f]
! \f[ A = (G_0 - G) / {\rm nucleate-vol} / {\rm so4-molec-dens} \f]
! \f[ \dot{G} = - {\rm nucleate-vol} {\rm so4-molec-dens}
! {\rm nucleate-coeff} G^2 \f]

program test_nucleate_ode
  
  use pmc_util
  use pmc_constants

  !> Conversion from (molecules m^{-3}) to (ppb).
  real(kind=dp), parameter :: conc_to_ppb = 3.977d-17

  !> Sulfate density (kg m^{-3}).
  real(kind=dp), parameter :: so4_dens = 1800d0
  !> Sulfate molar weight (kg mol^{-1}).
  real(kind=dp), parameter :: so4_molar_weight = 96d-3

  !> Coefficent for nucleation rate (m^3 s^{-1}).
  real(kind=dp), parameter :: nucleate_coeff = 1d-18
  !> Diameter of nucleated aerosol particles (m).
  real(kind=dp), parameter :: nucleate_diam = 1d-9

  !> Initial gas-phase H2SO4 mixing ratio (ppb).
  real(kind=dp), parameter :: init_h2so4_mix_rat = 1d-2
  !> Total simulation time (s).
  real(kind=dp), parameter :: t_max = 3600d0
  !> Timestep (s).
  real(kind=dp), parameter :: del_t = 1d0
  !> How often to print progress (s).
  real(kind=dp), parameter :: t_progress = 600d0
  !> How often to print output (s).
  real(kind=dp), parameter :: t_output = 60d0
  !> Output unit number for gas mixing ratio.
  integer, parameter :: unit_gas = 33
  !> Output unit number for aerosol number.
  integer, parameter :: unit_aero_time = 34
  !> Output filename for gas mixing ratio.
  character(len=*), parameter :: name_gas = "out/nucleate_ode_gas.txt"
  !> Output filename for aerosol number.
  character(len=*), parameter :: name_aero_time &
       = "out/nucleate_ode_aero_time.txt"
  !> Output filename for aerosol mass.

  real(kind=dp) :: init_h2so4_conc ! molecules / m^3
  real(kind=dp) :: nucleate_vol    ! m^3
  real(kind=dp) :: time            ! s
  real(kind=dp) :: h2so4_conc      ! molecules / m^3
  real(kind=dp) :: h2so4_mix_rat   ! ppb
  real(kind=dp) :: aero_num_conc   ! particles / m^3
  real(kind=dp) :: aero_mass_conc  ! kg / m^3
  real(kind=dp) :: so4_molec_dens  ! molecules / m^3
  integer :: n_step, i_step

  init_h2so4_conc = init_h2so4_mix_rat / conc_to_ppb
  nucleate_vol = const%pi / 6d0 * nucleate_diam**3
  so4_molec_dens = so4_dens / so4_molar_weight * const%avagadro

  open(unit=unit_gas, file=name_gas)
  open(unit=unit_aero_time, file=name_aero_time)
  time = 0d0
  h2so4_conc = init_h2so4_conc
  n_step = nint(t_max / del_t) + 1
  aero_num_conc = 0d0

  h2so4_mix_rat = h2so4_conc * conc_to_ppb
  aero_mass_conc = aero_num_conc * nucleate_vol * so4_dens
  write(*,'(a20,a20,a20)') &
       'time', 'aero_mass_conc', 'h2so4_mix_rat'
  write(*,'(e20.10,e20.10,e20.10)') &
       time, aero_mass_conc, h2so4_mix_rat
  write(unit_gas,'(e20.10,e20.10)') &
       time, h2so4_mix_rat
  write(unit_aero_time,'(e20.10,e20.10,e20.10,e20.10,e20.10)') &
       time, aero_num_conc, aero_mass_conc, 0d0, 0d0

  do i_step = 2,n_step
     time = dble(i_step - 1) * del_t
     call nucleate_step(h2so4_conc, del_t)
     aero_num_conc = (init_h2so4_conc - h2so4_conc) / nucleate_vol &
          / so4_molec_dens
     h2so4_mix_rat = h2so4_conc * conc_to_ppb
     aero_mass_conc = aero_num_conc * nucleate_vol * so4_dens
     if (mod(i_step - 1, nint(t_progress / del_t)) .eq. 0) then
        write(*,'(a20,a20,a20)') &
             'time', 'aero_mass_conc', 'h2so4_mix_rat'
        write(*,'(e20.10,e20.10,e20.10)') &
             time, aero_mass_conc, h2so4_mix_rat
     end if
     if (mod(i_step - 1, nint(t_output / del_t)) .eq. 0) then
        write(unit_gas,'(e20.10,e20.10)') &
             time, h2so4_mix_rat
        write(unit_aero_time,'(e20.10,e20.10,e20.10,e20.10,e20.10)') &
             time, aero_num_conc, aero_mass_conc, 0d0, 0d0
     end if
  end do

  close(unit_gas)
  close(unit_aero_time)
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  real(kind=dp) function nucleate_f(h2so4_conc)

    !> Concentration of gas-phase H2SO4 (molecules m^{-3}).
    real(kind=dp), intent(in) :: h2so4_conc

    nucleate_f = - nucleate_vol * so4_molec_dens * nucleate_coeff &
         * h2so4_conc**2
    
  end function nucleate_f
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine nucleate_step(h2so4_conc, del_t)
    
    !> Concentration of gas-phase H2SO4 (molecules m^{-3}).
    real(kind=dp), intent(inout) :: h2so4_conc
    !> Time-step (s).
    real(kind=dp), intent(in) :: del_t

    real(kind=dp) :: k1, k2, k3, k4
    
    ! integrate ODE with Runge-Kutta-4
    
    k1 = del_t * nucleate_f(h2so4_conc)
    k2 = del_t * nucleate_f(h2so4_conc + k1/2d0)
    k3 = del_t * nucleate_f(h2so4_conc + k2/2d0)
    k4 = del_t * nucleate_f(h2so4_conc + k3)
    
    h2so4_conc = h2so4_conc + k1/6d0 + k2/3d0 + k3/3d0 + k4/6d0
    
  end subroutine nucleate_step
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program test_nucleate_ode
