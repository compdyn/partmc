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
    interface ABIFM
        module procedure ABIFM_particle
        module procedure ABIFM_max
    end interface

contains

    subroutine freeze(aero_state, aero_data, env_state_initial, &
        env_state_final, del_t, do_freezing_CNT, freezing_rate, &
        do_coating, coating_spec, coating_ratio)

        implicit none
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

        integer :: i_part, i_bin, i_class, n_bins, n_class
        real(kind=dp) :: tmp
        logical :: do_freezing_CNT, do_coating
        real(kind=dp) :: freezing_rate
        real(kind=dp) :: a_w_ice, pis, pvs
        real(kind=dp) :: abifm_m, abifm_c
        real(kind=dp) :: p_freeze, p_frozen

        real(kind=dp) :: p_freeze_max, radius_max, diameter_max
        logical :: freeze_thisTime

        integer :: n_entry, ind, spec_bin = 13
        integer :: k_th, n_parts_in_bin
        real(kind=dp) :: aerosol_dry_radius, aerosol_radius, rand
        real(kind=dp), allocatable :: H2O_masses(:), total_masses(:), &
            H2O_frac(:)
        character(len=*), intent(in) :: coating_spec
        real(kind=dp) :: coating_ratio
        integer :: i_spec_max
        real(kind=dp) :: j_het_max

        integer :: loop_count = 0
        integer :: clock_start, clock_end

        call system_clock(clock_start)
        call aero_state_sort(aero_state, aero_data)

        total_masses = aero_state_masses(aero_state, aero_data)
        H2O_masses = aero_state_masses(aero_state, aero_data, include=(/"H2O"/))
        H2O_frac = H2O_masses / total_masses
        call env_state_saturated_vapor_pressure_water(env_state_final, pvs)
        call env_state_saturated_vapor_pressure_ice(env_state_final, pis)
        a_w_ice = pis / pvs
        call ABIFM_max_spec(aero_data, aero_state, a_w_ice, i_spec_max, j_het_max)
        !print*, i_spec_max, aero_data%name(i_spec_max)
        !print*, size(aero_state%awa%weight, 2)

        n_bins = aero_sorted_n_bin(aero_state%aero_sorted)
        n_class = aero_sorted_n_class(aero_state%aero_sorted)
        loop_count = 0

        loop_bins: do i_bin = 1, n_bins
            loop_classes: do i_class = 1, n_class
                !print*, "i_bin = ", i_bin, "i_class = ", i_class
                n_parts_in_bin = integer_varray_n_entry(aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
                !print*, "pass"
                radius_max = aero_state%aero_sorted%bin_grid%edges(i_bin + 1)
                diameter_max = radius_max * 2
                if (do_freezing_CNT) then
                    !print*, "Before ABIFM_max"
                    call ABIFM_max(diameter_max, p_freeze_max, aero_data, j_het_max, del_t)
                    !print*, "After ABIFM_max"
                else
                    p_freeze_max = 1 - exp(freezing_rate * del_t)
                endif
                !print*, "diameter_max = ", diameter_max

                k_th = n_parts_in_bin + 1
                loop_choosed_particles: do while(.TRUE.)
                    !print*, "k_th = ", k_th
                    rand = pmc_random_geometric(p_freeze_max)
                    k_th = k_th - rand
                    if (k_th <= 0) then
                        EXIT loop_choosed_particles
                    endif
                    loop_count = loop_count + 1
                    i_part = aero_state%aero_sorted%size_class%inverse(i_bin, i_class)%entry(k_th)
                    !print*, i_bin, k_th, i_part, aero_state_n_part(aero_state), rand
                    if (aero_state%apa%particle(i_part)%frozen) then
                        cycle
                    end if
                    if (H2O_frac(i_part) < 1e-2) then
                        cycle
                    end if
                    if (do_freezing_CNT) then
                        !print*, "Before ABIFM_particle"
                        call ABIFM_particle(i_part, p_freeze, aero_state, aero_data, a_w_ice, del_t)
                        !print*, "After ABIFM_particle"
                        if (p_freeze > p_freeze_max) then
                            print*, "Warning! p_freeze > p_freeze_max. &
                            p_freeze = ", p_freeze, "; p_freeze_max = "&
                            , p_freeze_max
                        endif
                        rand = pmc_random()
                        !print*, "p_freeze = ", p_freeze, "; p_freeze_max = ", &
                            !p_freeze_max, "; ratio = ",  p_freeze / p_freeze_max
                        if (rand < p_freeze / p_freeze_max) then
                            aero_state%apa%particle(i_part)%frozen = .TRUE.
                            !aero_state%apa%particle(i_part)%P_frozen = 1 - (1 - p_frozen) &
                            !    * (1 - p_freeze)
                        endif 
                    else
                         aero_state%apa%particle(i_part)%frozen = .TRUE.
                    endif
                enddo loop_choosed_particles
            enddo loop_classes
         enddo loop_bins
         !print*, "Total particles=", aero_state_n_part(aero_state), "; loop=", loop_count

        call system_clock(clock_end)
        !print*, clock_end - clock_start
        freeze_module_run_time = freeze_module_run_time + clock_end - clock_start
        !print*, freeze_module_run_time

    end subroutine freeze

    subroutine unfreeze(aero_state, aero_data, env_state_initial, &
        env_state_final)
        implicit none
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
                !aero_state%apa%particle(i_part)%P_frozen = 0d0   
                aero_state%apa%particle(i_part)%frozen = .FALSE.
            end do
        end if

    end subroutine unfreeze
    subroutine ABIFM_particle(i_part, P_freezing, aero_state, aero_data, &
            a_w_ice, del_t)
        !integer :: i_part
        !type(aero_data_t), intent(in) :: aero_data
        !type(env_state_t), intent(in) :: env_state
        !real(kind=dp), intent(out) :: P_freezing
        !type(aero_state_t), intent(inout) :: aero_state

        !real(kind=dp) :: aerosol_diameter
        !real(kind=dp) :: immersed_surface_area
        !real(kind=dp) :: a_w_ice
        !real(kind=dp) :: abifm_m, abifm_c
        !real(kind=dp) :: del_t
        !real(kind=dp) :: j_het
        !real(kind=dp) :: pvs, pis

        !!real(kind=dp) :: aero_particle_diameter

        !aerosol_diameter =  aero_particle_dry_diameter(aero_state%apa%particle(i_part), aero_data)
       
        !immersed_surface_area = const%pi * aerosol_diameter **2
        !call env_state_saturated_vapor_pressure_water(env_state, pvs)
        !call env_state_saturated_vapor_pressure_ice(env_state, pis)
        !a_w_ice = pis / pvs
        !!print*, "pis =", pis, "pvs =", pvs, "a_w_ice = ", a_w_ice
        !j_het = 10 ** (abifm_m  * (1 - a_w_ice) + abifm_c)
        !!print*, j_het
        !P_freezing = 1 - exp(-j_het * immersed_surface_area  * del_t)
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
        

        !real(kind=dp) :: aero_particle_diameter
        !print*, "i_part = ", i_part
        !print*, "n_part = ", aero_state_n_part(aero_state)
        aerosol_diameter =  aero_particle_dry_diameter(aero_state%apa%particle(i_part), aero_data)
        !print*, "aerosol_diameter = ", aerosol_diameter
       
        immersed_surface_area = const%pi * aerosol_diameter **2

        
        total_vol = 0d0
        do i_spec = 1,aero_data_n_spec(aero_data)
            if (i_spec == aero_data%i_water) then
                cycle
            end if
            total_vol = total_vol + aero_state%apa%particle(i_part)%vol(i_spec)
        end do

        j_het_x_aera = 0d0
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
            surface_ratio = aero_state%apa%particle(i_part)%vol(i_spec) / total_vol
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

    end subroutine ABIFM_particle

    subroutine ABIFM_max_spec(aero_data, aero_state, a_w_ice, i_spec_max, j_het_max)
        implicit none
        !type(env_state_t), intent(in) :: env_state
        type(aero_data_t), intent(in) :: aero_data
        type(aero_state_t), intent(in) :: aero_state
        real(kind=dp) :: a_w_ice
        real(kind=dp) :: abifm_m, abifm_c
        real(kind=dp) :: j_het, j_het_max
        !real(kind=dp) :: pvs, pis
        integer :: i_spec
        integer, intent(out) :: i_spec_max
               
        j_het_max = -9999.0
        do i_spec = 1,aero_data_n_spec(aero_data)
            if (i_spec == aero_data%i_water) then
                cycle
            end if
            abifm_m = aero_data%abifm_m(i_spec)
            abifm_c = aero_data%abifm_c(i_spec)
            !print*, "pis =", pis, "pvs =", pvs, "a_w_ice = ", a_w_ice
            j_het = 10 ** (abifm_m  * (1 - a_w_ice) + abifm_c) * 10000
            if (j_het > j_het_max) then
                j_het_max = j_het
                i_spec_max = i_spec
            end if
        end do

    end subroutine ABIFM_max_spec

    subroutine ABIFM_max(diameter_max, P_freezing, aero_data, j_het_max, del_t)
        implicit none
        !type(env_state_t), intent(in) :: env_state
        real(kind=dp), intent(out) :: P_freezing
        type(aero_data_t), intent(in) :: aero_data
        real(kind=dp) :: diameter_max
        real(kind=dp) :: immersed_surface_area
        !real(kind=dp) :: a_w_ice
        !real(kind=dp) :: abifm_m, abifm_c
        real(kind=dp) :: del_t
        real(kind=dp) :: j_het_max
        !real(kind=dp) :: pvs, pis
        integer :: i_spec
               
        immersed_surface_area = const%pi * diameter_max **2
        P_freezing = 1 - exp(-j_het_max * immersed_surface_area  * del_t)

    end subroutine ABIFM_max
    
end module pmc_freezing
