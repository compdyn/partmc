! Copyright (C) 2010 Matthew West
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
! \f[ \dot{G} = - {\rm nucleate-vol} {\rm so4-molec-dens} {\rm nucleate-coeff} G^2 \f]

program test_nucleate_ode
  
  use pmc_util
  use pmc_constants

  !> Conversion from (molecules m^{-3}) to (ppb).
  real(kind=dp), parameter :: conc_to_ppb = 3.977d-5
  !> Sulfate molecular density (molecules m^{-3}).
  real(kind=dp), parameter :: so4_molec_dens = 1.129151585625d28

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
  real(kind=dp), parameter :: t_progress = 60d0
  !> Output unit number.
  integer, parameter :: out_unit = 33
  !> Output filename.
  character(len=*), parameter :: out_name = "out/nucleate_ode.txt"

  real(kind=dp) :: init_h2so4_conc, nucleate_vol
  real(kind=dp) :: time, h2so4_conc, aero_conc
  integer :: n_step, i_step

  init_h2so4_conc = init_h2so4_mix_rat / conc_to_ppb
  nucleate_vol = const%pi / 6d0 * nucleate_diam**3

  open(unit=out_unit, file=out_name)
  time = 0d0
  h2so4_conc = init_h2so4_conc
  n_step = nint(t_max / del_t) + 1
  aero_conc = 0d0
  write(*,'(a8,a14,a14)') 'time', 'aero_conc', 'h2so4_conc'
  write(*,'(f8.1,e14.5,e14.5)') time, aero_conc, h2so4_conc
  write(out_unit,'(e20.10,e20.10,e20.10)') &
       time, aero_conc, h2so4_conc * conc_to_ppb
  do i_step = 2,n_step
     time = dble(i_step - 1) * del_t
     call nucleate_step(h2so4_conc, del_t)
     aero_conc = (init_h2so4_conc - h2so4_conc) / nucleate_vol &
          / so4_molec_dens
     if (mod(i_step - 1, nint(t_progress / del_t)) .eq. 0) then
        write(*,'(a8,a14,a14)') 'time', 'aero_conc', 'h2so4_conc'
        write(*,'(f8.1,e14.5,e14.5)') time, aero_conc, h2so4_conc
        write(out_unit,'(e20.10,e20.10,e20.10)') &
             time, aero_conc, h2so4_conc * conc_to_ppb
     end if
  end do

  close(out_unit)
  
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
