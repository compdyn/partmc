! ---
program theoretical_freezing
    use pmc_env_state
    implicit none
    integer, parameter :: out_unit = 65
    real(kind=dp), parameter :: temperature = 253.15
    real(kind=dp), parameter :: total_time = 600
    real(kind=dp), parameter :: out_dt = 10
    real(kind=dp), parameter :: abifm_m_1 = 17.62106, abifm_c_1 = 1.42411
    real(kind=dp), parameter :: abifm_m_2 = 54.48075, abifm_c_2 = -10.66873
    real(kind=dp), parameter :: s_ratio_1 = 0.5, s_ratio_2 = 0.5
    real(kind=dp), parameter :: Dp_dry = 1d-6
    character(len=*), parameter :: out_filename = "out/freezing_theoretical_data.txt"
    real(kind=dp) :: time, Jhet_1, Jhet_2, Jhet_mean, Phi_mean, frozen_fraction
    integer :: ios
    ! open output
    open(unit=out_unit, file=out_filename, iostat=ios)
    if (ios /= 0) then
        write(0,'(a,a,a,i4)') 'ERROR: unable to open file ', &
          trim(out_filename), ' for writing: ', ios
        stop 1
    end if
    call compute_Jhet_ABIFM(temperature, abifm_m_1, abifm_c_1, Jhet_1)
    call compute_Jhet_ABIFM(temperature, abifm_m_2, abifm_c_2, Jhet_2)
    Jhet_mean = Jhet_1 * s_ratio_1 + Jhet_2 * s_ratio_2

    time = 0
    do while(time .le. total_time)
        Phi_mean = Jhet_mean * time
        frozen_fraction = 1 - exp(- const%pi * Dp_dry**2 * Phi_mean)
        write(out_unit,'(e20.10)') frozen_fraction
        time = time + out_dt
    end do

    close(out_unit)

    contains

    subroutine compute_Jhet_ABIFM(T, abifm_m, abifm_c, Jhet)
        use pmc_env_state
        use pmc_constants
        implicit none
        real(kind=dp), intent(in) :: T, abifm_m, abifm_c
        real(kind=dp), intent(out) :: Jhet
        real(kind=dp) :: es, ei, a_w_ice
        call env_state_saturated_vapor_pressure_water_2(T, es)
        call env_state_saturated_vapor_pressure_ice_2(T, ei)
        a_w_ice = ei / es
        Jhet = 10 ** (abifm_m  * (1 - a_w_ice) + abifm_c) * 10000

    end subroutine compute_Jhet_ABIFM
end program theoretical_freezing


