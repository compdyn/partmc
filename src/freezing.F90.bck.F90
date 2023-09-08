! Copyright (C) 2021 University of Illinois at Urbana-Champaign
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!tate_saturated_vapor_pressure_ice The pmc_freezing module.

module pmc_freezing
    use pmc_aero_state
    use pmc_env_state
    use pmc_aero_data
    use pmc_util
    use pmc_aero_particle
    use pmc_constants
    use pmc_rand
    implicit none

    integer :: freeze_module_run_time = 0

contains

    subroutine freeze(aero_state, aero_data, env_state_initial, &
        env_state_final, del_t, do_freezing_CNT, freezing_rate, &
        do_coating, coating_spec, coating_ratio)

        !> Aerosol state.
        type(aero_state_t), intent(inout) :: aero_state
        !> Aerosol data.
        type(aero_data_t), intent(in) :: aero_data
        !> Environment state at the start of the timestep.
        type(env_state_t), intent(in) :: env_state_initial
        !> Environment state at the end of the timestep. The rel_humid
        !> value will be ignored and overwritten with a new value.
        type(env_state_t), intent(inout) :: env_state_final
        !> Total time to integrate.
        real(kind=dp), intent(in) :: del_t

        integer :: i_part
        real(kind=dp) :: tmp
        logical :: do_freezing_CNT, do_coating
        real(kind=dp) :: freezing_rate
        real(kind=dp) :: a_w_ice, pis, pvs
        !real(kind=dp) :: abifm_m, abifm_c
        real(kind=dp) :: p_freeze, p_frozen
        real(kind=dp), allocatable :: H2O_masses(:), total_masses(:), &
            H2O_frac(:)
        character(len=*), intent(in) :: coating_spec
        real(kind=dp) :: coating_ratio
        real(kind=dp) :: rand
        integer :: clock_start, clock_end

        call system_clock(clock_start)
        total_masses = aero_state_masses(aero_state, aero_data)
        H2O_masses = aero_state_masses(aero_state, aero_data, include=(/"H2O"/))
        H2O_frac = H2O_masses / total_masses
        !print*, sum(H2O_frac) / aero_state_n_part(aero_state), &
        !    env_state_final%rel_humid, sum(H2O_masses)
        !print*, aero_data%abifm_m
        !print*, aero_data%abifm_c
        call env_state_saturated_vapor_pressure_water(env_state_final, pvs)
        call env_state_saturated_vapor_pressure_ice(env_state_final, pis)
        a_w_ice = pis / pvs

        do i_part = 1, aero_state_n_part(aero_state)
            if (H2O_frac(i_part) < 1e-2) then
                cycle
            end if
            rand = pmc_random()
            !p = 1 - exp(-J_het * immersed_surface_area[i] * del_t)
!p = 1 - exp(-.00123456789 * del_t)
            !print*, "Begin CNT"
            if (do_freezing_CNT) then
                call ABIFM(i_part, p_freeze, aero_state, aero_data, a_w_ice, del_t, &
                        do_coating, coating_spec, coating_ratio)
            else
                p_freeze = 1 - exp(freezing_rate * del_t)
            end if
            !print*, "freezing_rate:", freezing_rate, "; p=", p

            !aero_state%apa%particle(i_part)%P_frozen = p_freeze
            p_frozen = aero_state%apa%particle(i_part)%P_frozen
            aero_state%apa%particle(i_part)%P_frozen = 1 - (1 - p_frozen) &
                * (1 - p_freeze)

            if (rand < p_freeze) then
                aero_state%apa%particle(i_part)%frozen = .TRUE.
            end if
        end do
        call system_clock(clock_end)
        freeze_module_run_time = freeze_module_run_time + clock_end - clock_start

    end subroutine freeze

    subroutine unfreeze(aero_state, aero_data, env_state_initial, &
        env_state_final)
        !> Aerosol state.
        type(aero_state_t), intent(inout) :: aero_state
        !> Aerosol data.
        type(aero_data_t), intent(in) :: aero_data
        !> Environment state at the start of the timestep.
        type(env_state_t), intent(in) :: env_state_initial
        !> Environment state at the end of the timestep. The rel_humid
        !> value will be ignored and overwritten with a new value.
        type(env_state_t), intent(inout) :: env_state_final
        integer :: i_part 

        if (env_state_final%temp > const%water_freeze_temp) then
            do i_part = 1, aero_state_n_part(aero_state)
                aero_state%apa%particle(i_part)%P_frozen = 0d0   
                aero_state%apa%particle(i_part)%frozen = .FALSE.
            end do
        end if

    end subroutine unfreeze

    subroutine ABIFM(i_part, P_freezing, aero_state, aero_data, a_w_ice, del_t, &
            do_coating, coating_spec, coating_ratio)
        implicit none
        integer :: i_part, i_spec
        type(aero_data_t), intent(in) :: aero_data
        !type(env_state_t), intent(in) :: env_state
        real(kind=dp), intent(out) :: P_freezing
        type(aero_state_t), intent(inout) :: aero_state

        real(kind=dp) :: aerosol_diameter
        real(kind=dp) :: immersed_surface_area
        real(kind=dp) :: total_vol
        real(kind=dp) :: surface_ratio
        real(kind=dp) :: a_w_ice
        real(kind=dp) :: abifm_m, abifm_c
        real(kind=dp) :: del_t
        real(kind=dp) :: j_het, j_het_x_aera
        real(kind=dp) :: pvs, pis
        logical :: do_coating
        character(len=*), intent(in) :: coating_spec
        real(kind=dp) :: coating_ratio
        integer :: i_coat_spec
        real(kind=dp) :: coat_spec_vol
        !real(kind=dp) :: aero_particle_diameter

        aerosol_diameter =  aero_particle_dry_diameter(aero_state%apa%particle(i_part), aero_data)
       
        immersed_surface_area = const%pi * aerosol_diameter **2

        if (do_coating) then
            i_coat_spec = string_array_find(aero_data%name, coating_spec)    
        endif
        !print*, i_coat_spec, aero_data%name(i_coat_spec)
        total_vol = 0d0
        do i_spec = 1,aero_data_n_spec(aero_data)
            if (i_spec == aero_data%i_water) then
                cycle
            end if
            total_vol = total_vol + aero_state%apa%particle(i_part)%vol(i_spec)
        end do
        if (do_coating) then
            coat_spec_vol = aero_state%apa%particle(i_part)%vol(i_coat_spec)
        endif

        j_het_x_aera = 0d0
        !print*, "---------------"
        do i_spec = 1,aero_data_n_spec(aero_data)
            if (i_spec == aero_data%i_water) then
                cycle
            end if
            !call env_state_saturated_vapor_pressure_water(env_state, pvs)
            !call env_state_saturated_vapor_pressure_ice(env_state, pis)
            !a_w_ice = pis / pvs
            !print*, "pis =", pis, "pvs =", pvs, "a_w_ice = ", a_w_ice
            abifm_m = aero_data%abifm_m(i_spec)
            abifm_c = aero_data%abifm_c(i_spec)

            ! Coating with Illite
            !abifm_m = 54.48075
            !abifm_c = -10.66873

            j_het = 10 ** (abifm_m  * (1 - a_w_ice) + abifm_c) * 10000
            !print*, "area = ", immersed_surface_area, "aw = ", a_w_ice, "rate = ", &
            !    j_het, "abifm_m = ", abifm_m, "abifm_c = ", abifm_c
            if (do_coating) then
                if (i_spec == i_coat_spec) then
                    surface_ratio = coating_ratio
                else
                    surface_ratio = aero_state%apa%particle(i_part)%vol(i_spec)&
                        / (total_vol - coat_spec_vol) * (1 - coating_ratio)
                endif
            else
                surface_ratio = aero_state%apa%particle(i_part)%vol(i_spec) / total_vol
            endif
            !print*, aero_data%name(i_spec), surface_ratio
            !print*, i_spec, surface_ratio
            j_het_x_aera = j_het_x_aera + j_het * immersed_surface_area * &
                surface_ratio
            

        end do
            !print*, j_het
            !immersed_surface_area = 1.8150015899967208e-12
            !P_freezing = 1 - exp(-j_het * immersed_surface_area  * del_t)
        P_freezing = 1 - exp(-j_het_x_aera * del_t)
            !print*, "area = ", immersed_surface_area,  "aw = ", a_w_ice, "rate = ",& 
            !    j_het, "P = ", P_freezing

    end subroutine ABIFM

end module pmc_freezing
