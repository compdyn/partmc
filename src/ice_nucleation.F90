! Copyright (C) 2024 University of Illinois at Urbana-Champaign
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!tate_saturated_vapor_pressure_ice The pmc_ice_nucleation module.

module pmc_ice_nucleation
    use pmc_aero_state
    use pmc_env_state
    use pmc_aero_data
    use pmc_util
    use pmc_aero_particle
    use pmc_constants
    use pmc_rand

    implicit none

    integer, parameter :: IMMERSION_FREEZING_SCHEME_INVALID = 0
    integer, parameter :: IMMERSION_FREEZING_SCHEME_CONST = 1
    integer, parameter :: IMMERSION_FREEZING_SCHEME_SINGULAR = 2
    integer, parameter :: IMMERSION_FREEZING_SCHEME_ABIFM = 3

    !> Used to record the runtime of the ice nucleation module.
    integer :: freeze_module_run_time = 0
    !> True: using the binned-tau leaping algorithm for time-dependent scheme.
    logical :: do_speedup = .True.
    !> False: using the naive algorithm for time-dependent scheme.
    !logical :: do_speedup = .False.

    interface ABIFM
        module procedure ABIFM_particle
        module procedure ABIFM_max
    end interface

contains

    !> Main subroutine for immersion freezing simulation.
    subroutine immersion_freezing(aero_state, aero_data, env_state_initial, &
        env_state_final, del_t, immersion_freezing_scheme_type, freezing_rate)
        !do_coating, coating_spec, coating_ratio)
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
        !> Flag for coating effect.
        !logical ::  do_coating
        !> Immersion freezing scheme (e.g., ABIFM, sigular, constant rate ...).
        !character(len=*), intent(in) :: immersion_freezing_scheme
        integer, intent(in) :: immersion_freezing_scheme_type
        !> Freezing rate (only used for the constant rate scheme).
        real(kind=dp) :: freezing_rate
        !> Coating species (only used for the coating effect).
        !character(len=*), intent(in) :: coating_spec
        !> Surface ratio coated (only used for the coating effect).
        !real(kind=dp) :: coating_ratio
        !> Start and end time (used for timing).
        integer :: clock_start, clock_end

        !> Record the start time.
        call system_clock(clock_start)

        !> Call the immersion freezing subroutine according to the immersion
        !> freezing scheme.
        if ((immersion_freezing_scheme_type .eq. IMMERSION_FREEZING_SCHEME_ABIFM) &
            .OR. (immersion_freezing_scheme_type .eq. IMMERSION_FREEZING_SCHEME_CONST)) then
            if (do_speedup) then
                call immersion_freezing_time_dependent(aero_state, aero_data, env_state_initial, &
                    env_state_final, del_t, immersion_freezing_scheme_type, freezing_rate)!, &
                    !do_coating, coating_spec, coating_ratio)
            else
                call immersion_freezing_time_dependent_naive(aero_state, aero_data, env_state_initial, &
                    env_state_final, del_t, immersion_freezing_scheme_type, freezing_rate)!, &
                    !do_coating, coating_spec, coating_ratio)
            end if
        else if (immersion_freezing_scheme_type .eq. &
                IMMERSION_FREEZING_SCHEME_SINGULAR) then
            call immersion_freezing_singular(aero_state, aero_data, env_state_initial, &
                env_state_final)
        else
            call assert_msg(121370299, .false., &
                'Error type of immersion freezing scheme')
        endif
        !> Record the end time.
        call system_clock(clock_end)
        !> Add to the total time.
        freeze_module_run_time = freeze_module_run_time + clock_end - clock_start

    end subroutine immersion_freezing

    !> Initialization for the sigular scheme, sampling the freezing temperature
    !> for each particles.
    subroutine singular_initialize(aero_state, aero_data)
        implicit none
         !> Aerosol state.
        type(aero_state_t), intent(inout) :: aero_state
        !> Aerosol data.
        type(aero_data_t), intent(in) :: aero_data
        integer :: i_part
        real(kind=dp) :: a_INAS, b_INAS, p, S, T0, temp
        real(kind=dp) :: aerosol_diameter

        T0 = const%water_freeze_temp
        do i_part = 1, aero_state_n_part(aero_state)
            a_INAS = -0.517 
            b_INAS = 8.934
            aerosol_diameter = aero_particle_dry_diameter(aero_state%apa%particle(i_part), aero_data)
            S = const%pi * aerosol_diameter **2
            p = pmc_random()
            temp = (log(1 - p) + exp(-S * exp(-a_INAS * T0 + b_INAS))) / (-S)
            aero_state%apa%particle(i_part)%imf_temperature = T0 + (log(temp) &
                     - b_INAS) / a_INAS
            !print*, aerosol_diameter, aero_state%apa%particle(i_part)%imf_temperature
        end do
    end subroutine singular_initialize

    !> Simulation for singular scheme, deciding whether to freeze for each
    !> particle. Run in each time step.
    subroutine immersion_freezing_singular(aero_state, aero_data, env_state_initial, env_state_final)
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
        real(kind=dp), allocatable :: H2O_masses(:), total_masses(:), &
            H2O_frac(:)
        integer :: i_part, i_bin, i_class, n_bins, n_class
        
        allocate(total_masses(aero_state_n_part(aero_state)))
        allocate(H2O_masses(aero_state_n_part(aero_state)))
        allocate(H2O_frac(aero_state_n_part(aero_state)))

        total_masses = aero_state_masses(aero_state, aero_data)
        H2O_masses = aero_state_masses(aero_state, aero_data, include=(/"H2O"/))
        H2O_frac = H2O_masses / total_masses
        do i_part = 1, aero_state_n_part(aero_state)
            if (aero_state%apa%particle(i_part)%frozen) then
                cycle
            end if
            if (H2O_frac(i_part) < 1e-2) then
                cycle
            end if
            if (env_state_final%temp .le. &
                    aero_state%apa%particle(i_part)%imf_temperature) then
                aero_state%apa%particle(i_part)%frozen = .TRUE.
                aero_state%apa%particle(i_part)%den_ice = &
                        const%reference_ice_density
                aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
            end if
        end do
        deallocate(total_masses)
        deallocate(H2O_masses)
        deallocate(H2O_frac)

    end subroutine immersion_freezing_singular

    !> Simulation for time-dependent scheme (e.g., ABIFM, constant rate),
    !> deciding whether to freeze for each particle. Run in each time step.   
    !> This subroutine applys the binned-tau leaping algorithm for speeding up.
    subroutine immersion_freezing_time_dependent(aero_state, aero_data, env_state_initial, &
        env_state_final, del_t, immersion_freezing_scheme_type, freezing_rate)!, &
        !do_coating, coating_spec, coating_ratio)

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
        !logical :: do_freezing_CNT, do_coating
        !logical ::  do_coating
        !character(len=*), intent(in) :: immersion_freezing_scheme
        integer, intent(in) :: immersion_freezing_scheme_type
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
        !character(len=*), intent(in) :: coating_spec
        !real(kind=dp) :: coating_ratio
        integer :: i_spec_max
        real(kind=dp) :: j_het_max
        integer :: rand_geometric

        integer :: loop_count = 0
        
        allocate(total_masses(aero_state_n_part(aero_state)))
        allocate(H2O_masses(aero_state_n_part(aero_state)))
        allocate(H2O_frac(aero_state_n_part(aero_state)))

        call aero_state_sort(aero_state, aero_data)

        total_masses = aero_state_masses(aero_state, aero_data)
        H2O_masses = aero_state_masses(aero_state, aero_data, include=(/"H2O"/))
        H2O_frac = H2O_masses / total_masses
        call env_state_saturated_vapor_pressure_water(env_state_final, pvs)
        call env_state_saturated_vapor_pressure_ice(env_state_final, pis)
        a_w_ice = pis / pvs
        if (immersion_freezing_scheme_type .eq. IMMERSION_FREEZING_SCHEME_ABIFM) then
            call ABIFM_max_spec(aero_data, aero_state, a_w_ice, i_spec_max, j_het_max)
        endif

        n_bins = aero_sorted_n_bin(aero_state%aero_sorted)
        n_class = aero_sorted_n_class(aero_state%aero_sorted)
        loop_count = 0

        loop_bins: do i_bin = 1, n_bins
            loop_classes: do i_class = 1, n_class
                n_parts_in_bin = integer_varray_n_entry(aero_state%aero_sorted%size_class%inverse(i_bin, i_class))
                radius_max = aero_state%aero_sorted%bin_grid%edges(i_bin + 1)
                diameter_max = radius_max * 2
                if (immersion_freezing_scheme_type .eq. &
                        IMMERSION_FREEZING_SCHEME_ABIFM) then
                    call ABIFM_max(diameter_max, p_freeze_max, aero_data, j_het_max, del_t)
                else if (immersion_freezing_scheme_type .eq. &
                        IMMERSION_FREEZING_SCHEME_CONST) then
                    p_freeze_max = 1 - exp(freezing_rate * del_t)
                endif

                k_th = n_parts_in_bin + 1
                loop_choosed_particles: do while(.TRUE.)
                    rand_geometric = pmc_random_geometric(p_freeze_max)
                    k_th = k_th - rand_geometric
                    if (k_th <= 0) then
                        EXIT loop_choosed_particles
                    endif
                    loop_count = loop_count + 1
                    i_part = aero_state%aero_sorted%size_class%inverse(i_bin, i_class)%entry(k_th)
                    if (aero_state%apa%particle(i_part)%frozen) then
                        cycle
                    end if
                    if (H2O_frac(i_part) < 1e-2) then
                        cycle
                    end if
                    if (immersion_freezing_scheme_type .eq. &
                            IMMERSION_FREEZING_SCHEME_ABIFM) then
                        call ABIFM_particle(i_part, p_freeze, aero_state, aero_data, a_w_ice, del_t)
                        if (p_freeze > p_freeze_max) then
                            print*, "Warning! p_freeze > p_freeze_max. "&
                            ,"p_freeze = ", p_freeze, "; p_freeze_max = "&
                            , p_freeze_max
                        endif
                        rand = pmc_random()
                        if (rand < p_freeze / p_freeze_max) then
                            !print*, i_part
                            aero_state%apa%particle(i_part)%frozen = .TRUE.
                            aero_state%apa%particle(i_part)%den_ice = &
                                const%reference_ice_density
                            aero_state%apa%particle(i_part)%ice_shape_phi = 1d0

                        endif 
                    else
                         aero_state%apa%particle(i_part)%frozen = .TRUE.
                         aero_state%apa%particle(i_part)%den_ice = &
                                const%reference_ice_density
                        aero_state%apa%particle(i_part)%ice_shape_phi = 1d0

                    endif
                enddo loop_choosed_particles
            enddo loop_classes
         enddo loop_bins

         deallocate(total_masses)
         deallocate(H2O_masses)
         deallocate(H2O_frac)
       
    end subroutine immersion_freezing_time_dependent

    !> Simulation for time-dependent scheme (e.g., ABIFM, constant rate),
    !> deciding whether to freeze for each particle. Run in each time step.   
    !> This subroutine applys the naive algorithm for reference.
    subroutine immersion_freezing_time_dependent_naive(aero_state, aero_data, &
        env_state_initial, env_state_final, del_t, &
        immersion_freezing_scheme_type, &
        freezing_rate)!, do_coating, coating_spec, coating_ratio)

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
        !logical :: do_coating
        !character(len=*), intent(in) :: immersion_freezing_scheme
        integer, intent(in) :: immersion_freezing_scheme_type
        real(kind=dp) :: freezing_rate
        real(kind=dp) :: a_w_ice, pis, pvs
        !real(kind=dp) :: abifm_m, abifm_c
        real(kind=dp) :: p_freeze, p_frozen
        real(kind=dp), allocatable :: H2O_masses(:), total_masses(:), &
            H2O_frac(:)
        !character(len=*), intent(in) :: coating_spec
        !real(kind=dp) :: coating_ratio
        real(kind=dp) :: rand
        integer :: clock_start, clock_end

        call system_clock(clock_start)

        allocate(total_masses(aero_state_n_part(aero_state)))
        allocate(H2O_masses(aero_state_n_part(aero_state)))
        allocate(H2O_frac(aero_state_n_part(aero_state)))

        total_masses = aero_state_masses(aero_state, aero_data)
        H2O_masses = aero_state_masses(aero_state, aero_data, include=(/"H2O"/))
        H2O_frac = H2O_masses / total_masses
        
        call env_state_saturated_vapor_pressure_water(env_state_final, pvs)
        call env_state_saturated_vapor_pressure_ice(env_state_final, pis)
        a_w_ice = pis / pvs

        do i_part = 1, aero_state_n_part(aero_state)
            !print*, i_part, H2O_frac(i_part)
            if (aero_state%apa%particle(i_part)%frozen) then
                cycle
            end if
            if (H2O_frac(i_part) < 1e-2) then
                cycle
            end if
            rand = pmc_random()
            
            if (immersion_freezing_scheme_type .eq. &
                    IMMERSION_FREEZING_SCHEME_ABIFM) then
                call ABIFM_particle(i_part, p_freeze, aero_state, aero_data, a_w_ice, del_t) !do_coating, coating_spec, coating_ratio)
            else if (immersion_freezing_scheme_type .eq. &
                    IMMERSION_FREEZING_SCHEME_CONST) then
                p_freeze = 1 - exp(freezing_rate * del_t)
            end if
            p_frozen = aero_state%apa%particle(i_part)%P_frozen
            aero_state%apa%particle(i_part)%P_frozen = 1 - (1 - p_frozen) &
                * (1 - p_freeze)

            if (rand < p_freeze) then
                aero_state%apa%particle(i_part)%frozen = .TRUE.
                aero_state%apa%particle(i_part)%den_ice = &
                        const%reference_ice_density
                aero_state%apa%particle(i_part)%ice_shape_phi = 1d0
            end if
        end do

        deallocate(total_masses)
        deallocate(H2O_masses)
        deallocate(H2O_frac)

        call system_clock(clock_end)
        freeze_module_run_time = freeze_module_run_time + clock_end - clock_start

    end subroutine immersion_freezing_time_dependent_naive

    !> This subroutine simulates the melting process.(Set frozen = .False.
    !> for each particle once the temperature is higher than water freezing
    !> tempearture)
    subroutine melting(aero_state, aero_data, env_state_initial, &
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
                aero_state%apa%particle(i_part)%frozen = .FALSE.
                aero_state%apa%particle(i_part)%den_ice = -9999d0
                aero_state%apa%particle(i_part)%ice_shape_phi = -9999d0
            end do
        end if

    end subroutine melting

    !> Calculating the freezing probability for one particle using ABIFM 
    !> method (Knopf et al.,2013)

    subroutine ABIFM_particle(i_part, P_freezing, aero_state, aero_data, &
            a_w_ice, del_t)
        
        implicit none
        integer :: i_part, i_spec
        type(aero_data_t), intent(in) :: aero_data
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
        

       
        aerosol_diameter =  aero_particle_dry_diameter(aero_state%apa%particle(i_part), aero_data)
       
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
            
            abifm_m = aero_data%abifm_m(i_spec)
            abifm_c = aero_data%abifm_c(i_spec)

           

            j_het = 10 ** (abifm_m  * (1 - a_w_ice) + abifm_c) * 10000
            surface_ratio = aero_state%apa%particle(i_part)%vol(i_spec) / total_vol
            j_het_x_aera = j_het_x_aera + j_het * immersed_surface_area * &
                surface_ratio
            

        end do
            
        P_freezing = 1 - exp(-j_het_x_aera * del_t)

    end subroutine ABIFM_particle

    !> Calculating the heterogeneous ice nucleation rate coefficient 
    !> for each species, and fining the species having the largest rate.
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

    !> Calculating the maximum freezing probability for particles in 
    !> one bin using ABIFM method (Knopf et al.,2013). Only used by
    !> the binned-tau leaping algorithm.
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

    !> Read the specification for a kernel type from a spec file and
    !> generate it.
    subroutine spec_file_read_immersion_freezing_scheme_type(file, immersion_freezing_scheme_type)

        !> Spec file.
        type(spec_file_t), intent(inout) :: file
        !> Kernel type.
        integer, intent(out) :: immersion_freezing_scheme_type

        character(len=SPEC_LINE_MAX_VAR_LEN) :: imf_scheme

        !> \page input_format_immersion_freezing_scheme Input File Format: Immersion freezing scheme
        !!
        !! The immersion freezing scheme is specified by the parameter:
        !!   - \b immersion_freezing_scheme (string): the type of immersion freezing scheme
        !!     must be one of: \c sedi for the gravitational sedimentation
        !!     kernel; \c additive for the additive kernel; \c constant
        !!     for the constant kernel; \c brown for the Brownian kernel,
        !!     or \c zero for no immersion freezing
        !!
        !! If \c immersion_freezing_scheme is \c additive, the kernel coefficient needs to be
        !! provided using the \c additive_kernel_coeff parameter
        !!
        !! See also:
        !!   - \ref spec_file_format --- the input file text format

        call spec_file_read_string(file, 'immersion_freezing_scheme', imf_scheme)
        if (trim(imf_scheme) == 'const') then
           immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_CONST
        elseif (trim(imf_scheme) == 'singular') then
           immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_SINGULAR
        elseif (trim(imf_scheme) == 'ABIFM') then
           immersion_freezing_scheme_type = IMMERSION_FREEZING_SCHEME_ABIFM
        else
           call spec_file_die_msg(920761229, file, &
                "Unknown immersion freezing scheme: " // trim(imf_scheme))
        end if

    end subroutine spec_file_read_immersion_freezing_scheme_type

end module pmc_ice_nucleation
