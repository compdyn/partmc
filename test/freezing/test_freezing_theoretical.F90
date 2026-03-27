! ---
! This program is used to calculate the theoretical frozen fraction of
! mono-disperse INPs composed of two species in internal mixing at a
! constant temperature during the immersion freezing process.
! All INPs are assumed to be immersed in supercooled water.
! Calculation of J_het is based on the ABIFM method.

! References:
! Tang, W. (2024). Particle-resolved simulations of immersion freezing with multi-species ice nucleating particles
! (Master's thesis, University of Illinois at Urbana-Champaign). https://hdl.handle.net/2142/124611

! Knopf, D. A., & Alpert, P. A. (2013). A water activity based model of
! heterogeneous ice nucleation kinetics for freezing of water and aqueous
! solution droplets. Faraday discussions, 165, 513-534. Doi:10.1039/C3FD00035D.


program theoretical_freezing

  use pmc_env_state
  use pmc_constants
  implicit none

  ! Temperature (unit: Kelvin)
  real(kind=dp), parameter :: temperature = 243.15d0
  ! Melting point (unit: Kelvin)
  real(kind=dp), parameter :: T0 = 273.15d0
  ! Time duration (unit: s)
  real(kind=dp), parameter :: total_time = 600.0d0
  ! Output time interval (unit: s)
  real(kind=dp), parameter :: out_dt = 60.0d0
  ! ABIFM parameters (m and c) for species 1.
  real(kind=dp), parameter :: abifm_m_1 = 17.62106d0, abifm_c_1 = 1.42411d0
  ! ABIFM parameters (m and c) for species 2.
  real(kind=dp), parameter :: abifm_m_2 = 54.48075d0, abifm_c_2 = -10.66873d0
  ! The ratio of total surface for species 1 and 2 (for ABIFM only).
  real(kind=dp), parameter :: s_ratio_1 = 0.5d0, s_ratio_2 = 0.5d0
  ! INAS parameters (a and b) for mineral dust. (for singular only)
  real(kind=dp), parameter :: inas_a = -0.517d0, inas_b = 8.934d0
  ! Constant freezing rate. (for const only)
  real(kind=dp), parameter :: freezing_rate = 0.01123456789d0
  ! Dry diameter of INPs (unit: m)
  real(kind=dp), parameter :: Dp_dry = 1d-6
  character(len=100) :: out_filename, imf_scheme
  integer, parameter :: out_unit = 65
  integer, parameter :: IMMERSION_FREEZING_SCHEME_CONST = 0
  integer, parameter :: IMMERSION_FREEZING_SCHEME_SINGULAR = 1
  integer, parameter :: IMMERSION_FREEZING_SCHEME_ABIFM = 2
  real(kind=dp) :: time, Jhet_1, Jhet_2, Jhet_mean, Phi_mean, frozen_fraction
  real(kind=dp) :: ns
  integer :: ios, immersion_freezing_scheme_type

  call getarg(1, imf_scheme)
  if (trim(imf_scheme) == 'const') then
     immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_CONST
  elseif (trim(imf_scheme) == 'singular') then
     immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_SINGULAR
  elseif (trim(imf_scheme) == 'ABIFM') then
     immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_ABIFM
  else
     print*, "Unknown immersion freezing scheme: " // trim(imf_scheme)
     immersion_freezing_scheme_type = -1
  endif
  out_filename = "out/freezing_theoretical_" // trim(imf_scheme) // &
       "_data.txt"
  ! open output
  open(unit=out_unit, file=out_filename, iostat=ios)
  if (ios /= 0) then
     write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
     stop 1
  end if

  if (immersion_freezing_scheme_type == &
       IMMERSION_FREEZING_SCHEME_ABIFM) then
     call compute_Jhet_ABIFM(temperature, abifm_m_1, abifm_c_1, Jhet_1)
     call compute_Jhet_ABIFM(temperature, abifm_m_2, abifm_c_2, Jhet_2)
     Jhet_mean = Jhet_1 * s_ratio_1 + Jhet_2 * s_ratio_2

     time = 0d0
     do while(time .le. total_time)
        Phi_mean = Jhet_mean * time
        frozen_fraction = 1d0 - exp(- const%pi * Dp_dry**2 * Phi_mean)
        write(out_unit,'(e20.10)') frozen_fraction
        time = time + out_dt
     end do
  elseif (immersion_freezing_scheme_type == &
       IMMERSION_FREEZING_SCHEME_SINGULAR) then
     time = 0d0
     do while(time .le. total_time)
        if (time .eq. 0d0) then
           frozen_fraction = 0.0d0
        else
           call compute_INAS(temperature, inas_a, inas_b, ns)
           frozen_fraction = 1d0 - exp(- const%pi * Dp_dry**2 * ns)
        end if
        write(out_unit,'(e20.10)') frozen_fraction
        time = time + out_dt
     end do
  elseif (immersion_freezing_scheme_type == &
       IMMERSION_FREEZING_SCHEME_CONST) then
     time = 0d0
     do while(time .le. total_time)
        frozen_fraction = 1d0 - exp(-freezing_rate * time)
        write(out_unit,'(e20.10)') frozen_fraction
        time = time + out_dt
     end do
  end if

  close(out_unit)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate J_het value using ABIFM method (Knopf et al., 2013)
  subroutine compute_Jhet_ABIFM(T, abifm_m, abifm_c, Jhet)

    real(kind=dp), intent(in) :: T, abifm_m, abifm_c
    real(kind=dp), intent(out) :: Jhet
    real(kind=dp) :: es, ei, a_w_ice
    es = env_state_saturated_vapor_pressure_wrt_water(T)
    ei = env_state_saturated_vapor_pressure_wrt_ice(T)
    a_w_ice = ei / es
    Jhet = 10d0 ** (abifm_m * (1d0 - a_w_ice) + abifm_c) * 10000d0

  end subroutine compute_Jhet_ABIFM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the INAS density (Niemand et al., 2012)
  subroutine compute_INAS(T, inas_a, inas_b, ns)

    real(kind=dp), intent(in) :: T, inas_a, inas_b
    real(kind=dp), intent(out) :: ns
    ns = exp(inas_a * (T - T0) + inas_b)

  end subroutine compute_INAS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program theoretical_freezing
