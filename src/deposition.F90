! Copyright (C) 2012 Jeffrey H Curtis
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_deposition module.

!> Aerosol deposition.
module pmc_deposition

  use pmc_aero_state
  use pmc_aero_data
  use pmc_constants
  use pmc_aero_particle
  use pmc_rand
  use pmc_env_state

  !> Land type
  integer, parameter :: EVERGREEN_NEEDLE = 1
  !>
  integer, parameter :: EVERGREEN_BROAD = 2
  !>
  integer, parameter :: DECIDUOUS_NEEDLE = 3
  !>
  integer, parameter :: DECIDUOUS_BROAD = 4
  !>
  integer, parameter :: FOREST_MIXED = 5
  !>
  integer, parameter :: GRASS = 6
  !>
  integer, parameter :: CROPLAND = 7
  !>
  integer, parameter :: DESERT = 8
  !>
  integer, parameter :: URBAN = 15

  !> Seasonal categories
  !>
  integer, parameter :: SUMMER = 1
  !>
  integer, parameter :: AUTUMN = 2
  !>
  integer, parameter :: HARVEST = 3
  !>
  integer, parameter :: WINTER = 4
  !>
  integer, parameter :: SPRING = 5

  !> Stability class
  integer, parameter :: PAS_A = 1
  !>
  integer, parameter :: PAS_B = 2
  !>
  integer, parameter :: PAS_C = 3
  !>
  integer, parameter :: PAS_D = 4
  !>
  integer, parameter :: PAS_E = 5
  !>
  integer, parameter :: PAS_F = 6

  !>
  real(kind=dp), parameter, dimension(5,15) :: a_ref =  reshape(&
       [[2.0, 2.0, 2.0, 2.0, 2.0], [5.0, 5.0, 5.0, 5.0, 5.0], &
       [2.0, 2.0, 5.0, 5.0, 2.0], [5.0, 5.0, 10.0, 10.0, 5.0], &
       [5.0, 5.0, 5.0, 5.0, 5.0], [2.0, 2.0, 5.0, 5.0, 2.0], &
       [2.0, 2.0, 5.0, 5.0, 2.0], [0., 0., 0., 0., 0.], [0., 0., 0., 0.,0.], &
       [10.0, 10.0, 10.0, 10.0, 10.0], [10.0, 10.0, 10.0, 10.0, 10.0], &
       [0., 0., 0., 0., 0.], [0., 0., 0., 0., 0.], [0., 0., 0., 0., 0.], &
       [10.0, 10.0, 10.0, 10.0, 10.0]],[5,15])
  !>
  real(kind=dp), parameter, dimension(15) :: alpha_ref = [ 1.0, 0.6, 1.1, &
       0.8, 0.8, 1.2, 1.2, 50.0, 50.0, 1.3, 2.0, 50.0, 100.0, 100.0, 1.5]
  !>
  real(kind=dp), parameter, dimension(15) :: gamma_ref = [ 0.56, 0.58, 0.56, &
       0.56, 0.56, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.54, 0.50, 0.50, 0.56]

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Remove particles from an aero state by dry deposition.
  subroutine dry_dep_aero_state(aero_state, aero_data, env_state, aer_res_a, &
       ustar, gamma, A, alpha, dt)

    !> Aerosol states.
    type(aero_state_t), intent(inout) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment states.
    type(env_state_t), intent(in) :: env_state
    !> Aerodynamic resistance (s/m).
    real(kind=dp), intent(in) :: aer_res_a
    !> Friction velocity (m/s).
    real(kind=dp), intent(in) :: ustar
    !> Parameter for collection efficiency from Brownian diffusion.
    real(kind=dp), intent(in) :: gamma
    !> Characteristic radius (mm).
    real(kind=dp), intent(in) :: A
    !> Parameter for collection efficiency from impaction.
    real(kind=dp), intent(in) :: alpha
    !> Timestep (s).
    real(kind=dp), intent(in) :: dt

    ! Deposition velocity array
    real(kind=dp), allocatable :: vd(:)
    ! Particle removal probability array
    real(kind=dp), allocatable :: remove_prob(:)

    allocate(vd(aero_state%apa%n_part))
    allocate(remove_prob(aero_state%apa%n_part))

    vd = 0.0d0
    remove_prob = 0.0d0

    ! Compute deposition velocity array
    call compute_dep_vel(aero_state, aero_data, aer_res_a, ustar, &
         env_state, gamma, A, alpha, vd)
    ! Compute the removal probability array
    call compute_dep_prob(vd, dt, env_state%height, remove_prob)
    ! Remove particles based on probabilities
    call aero_particles_remove_by_dep(aero_state, remove_prob)

    deallocate(vd)
    deallocate(remove_prob)

  end subroutine dry_dep_aero_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute deposition velocities.
  subroutine compute_dep_vel(aero_state, aero_data, aer_res_a, ustar, &
       env_state, gamma, A, alpha, vd)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerodynamic resistance (s/m)
    real(kind=dp), intent(in) :: aer_res_a
    !> Friction velocity (m/s).
    real(kind=dp), intent(in) :: ustar
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Parameter for collection efficiency from Brownian diffusion.
    real(kind=dp), intent(in) :: gamma
    !> Characteristic radius of large collectors (m).
    real(kind=dp), intent(in) :: A
    !> Parameter for collection efficiency from impaction.
    real(kind=dp), intent(in) :: alpha
    !> Dry deposition velocity (m/s).
    real(kind=dp), intent(inout) :: vd(:)

    integer :: i_part
    real(kind=dp) :: diameter
    real(kind=dp) :: density
    real(kind=dp) :: rs
    real(kind=dp) :: vs
    real(kind=dp) :: rho_air
    real(kind=dp) :: viscosd
    real(kind=dp) :: viscosk
    real(kind=dp) :: lambda

    ! Compute values that do not vary per particle
    rho_air = env_state_air_den(env_state)
    viscosd = env_state_air_dynamic_viscosity(env_state)
    viscosk = viscosd / rho_air
    lambda = env_state_air_mean_free_path(env_state)

    do i_part = 1,aero_state%apa%n_part
       diameter = aero_particle_diameter(aero_state%apa%particle(i_part))
       density = aero_particle_density(aero_state%apa%particle(i_part), &
            aero_data)
       vs = calculate_vs(diameter, density, viscosd, lambda)
       rs = calculate_rs(diameter, vs, env_state%temp, ustar, gamma, A, &
           alpha, lambda)
       vd(i_part) = calculate_vd(aer_res_a, rs, vs)
    end do

  end subroutine compute_dep_vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the probability of each particle being removed by dry deposition.
  subroutine compute_dep_prob(vd, dt, dz, remove_prob)

    !> Dry deposition velocities (m/s).
    real(kind=dp), intent(in) :: vd(:)
    !> Time step (s).
    real(kind=dp), intent(in) :: dt
    !> Delta z of depositing layer (m).
    real(kind=dp), intent(in) :: dz
    !> Probability array of removal.
    real(kind=dp), intent(inout) :: remove_prob(:)

    integer :: i_part

    do i_part = 1,size(vd)
       remove_prob(i_part) = (dt / dz) * vd(i_part)
    end do

  end subroutine compute_dep_prob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests the particle array to see if the particle is to be removed by
  !> dry deposition based on the provided probability.
  subroutine aero_particles_remove_by_dep(aero_state, remove_prob)

    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Probabilities of each particle being removed.
    real(kind=dp), intent(in) :: remove_prob(:)

    type(aero_info_t) :: aero_info
    integer :: i_part

    do i_part = aero_state%apa%n_part,1,-1
         if (pmc_random() < remove_prob(i_part)) then
            call aero_info_allocate(aero_info)
            aero_info%id = aero_state%apa%particle(i_part)%id
            aero_info%action = AERO_INFO_DEPOSITION
            aero_info%other_id = 0
            call aero_state_remove_particle_with_info(aero_state, i_part, &
                 aero_info)
            call aero_info_deallocate(aero_info)
         end if
    end do

  end subroutine aero_particles_remove_by_dep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the particle settling velocity.
  !> Seinfeld and Pandis eq. 19.18.
  real(kind=dp) function calculate_vs(diameter, density, viscosd, lambda)

     !> Particle diameter.
     real(kind=dp), intent(in) :: diameter
     !> Particle density.
     real(kind=dp), intent(in) :: density
     !>
     real(kind=dp), intent(in) :: viscosd
     !> Mean free path of air (m).
     real(kind=dp), intent(in) :: lambda

     real(kind=dp) :: C_c

     C_c = slip_correction_factor(diameter, lambda)

     calculate_vs = ((diameter)**2.0d0 * density * const%grav * C_c) &
          / (18.0d0 * viscosd)

  end function calculate_vs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate surface resistance of a given particle.
  real(kind=dp) function calculate_rs(diameter, vs, temperature, ustar, &
       gamma, A, alpha, lambda)

    !> Particle diameter (m).
    real(kind=dp), intent(in) :: diameter
    !> Settling velocity (m/s).
    real(kind=dp), intent(in) :: vs
    !> Temperature (K).
    real(kind=dp), intent(in) :: temperature
    !> Friction velocity (m/s).
    real(kind=dp), intent(in) :: ustar
    !> Parameter for collection efficiency from Brownian diffusion.
    real(kind=dp), intent(in) :: gamma
    !> Characteristic radius of large collectors (m).
    real(kind=dp), intent(in) :: A
    !> Parameter for collection efficiency from impaction.
    real(kind=dp), intent(in) :: alpha
    !> Mean free path of air (m).
    real(kind=dp), intent(in) :: lambda

    real(kind=dp) :: eps_0
    real(kind=dp) :: diff
    real(kind=dp) :: nu
    real(kind=dp) :: EB, EIM, EIN, R1
    real(kind=dp) :: beta
    real(kind=dp) :: Sc
    real(kind=dp) :: st

    ! Compute Schmidt number
    diff = compute_brownian_diff(diameter, temperature, lambda)
    nu = 1.45d-5
    Sc = nu / diff

    ! Compute Stokes number
    St = calculate_st(vs, ustar, A)

    ! EB: Collection efficiency from Brownian diffusion
    ! Equation 6 from Zhang et al. (2001)
    EB = Sc**(-gamma)

    ! EIM: Collection efficiency from impaction
    ! Equation 7b from Zhang et al. (2001)
    beta = 2.0d0
    EIM = (St / (alpha + St))**(beta)

    ! EIN: Collection efficiency from interception
    ! Equation 8 from Zhang et al. (2001)
    EIN = .5 * (diameter / A)**2.0d0

    ! R1: Correction factor for sticking
    ! Equation 9 from Zhang et al. (2001)
    R1 = exp(-St**.5)

    ! Empirical constant
    eps_0 = 3.0d0

    ! Equation 5 from Zhang et al. (2001)
    calculate_rs = 1.0d0 / (eps_0 * ustar * (EB + EIM + EIN) * R1)

  end function calculate_rs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the dry deposition velocity.
  !> Seinfeld and Pandis eq. 19.7.
  real(kind=dp) function calculate_vd(ra, rs, vs)

    !> Aerodynamic resistance (s/m)
    real(kind=dp), intent(in) :: ra
    !> Surface resistance (s/m).
    real(kind=dp), intent(in) :: rs
    !> Settling velocity (m/s).
    real(kind=dp), intent(in) :: vs

    real(kind=dp) :: tmp

    tmp = ra + rs + ra * rs * vs

    calculate_vd = (1.0d0 / tmp) + vs

  end function calculate_vd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the Cunningham slip correction factor for a single particle.
  !> Seinfeld and Pandis eq. 9.34.
  real(kind=dp) function slip_correction_factor(diameter, lambda)

     !> Particle diameter (m).
     real(kind=dp), intent(in) :: diameter
     !> Mean free path of air (m).
     real(kind=dp), intent(in) :: lambda

     slip_correction_factor = 1.0d0 + ((2.0d0 * lambda) / diameter) * ( 1.257 &
         + .4 * exp(-(1.1 * diameter) / (2.0d0 * lambda)))

  end function slip_correction_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the Stokes number.
  !> Zhang et al (2001) for vegetative surface.
  real(kind=dp) function calculate_st(vs, ustar, A)

   !> Settling velocity (m/s).
   real(kind=dp), intent(in) :: vs
   !> Friction velocity (m/s).
   real(kind=dp), intent(in) :: ustar
   !> Characteristic radius of large collectors (m).
   real(kind=dp), intent(in) :: A

   calculate_st = (vs * ustar) / (const%grav * A)

  end function calculate_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the Brownian diffusivity of a particle.
  !> Seinfeld and Pandis eq. 9.73.
  real(kind=dp) function compute_brownian_diff(diameter, temperature, lambda)

    !> Particle diameter (m).
    real(kind=dp), intent(in) :: diameter
    !> Temperature (K).
    real(kind=dp), intent(in) :: temperature
    !> Mean free path of air (m).
    real(kind=dp), intent(in) :: lambda

    real(kind=dp) :: C_c

    C_c = slip_correction_factor(diameter, lambda)

    compute_brownian_diff = (const%boltzmann * temperature * C_c) / &
        (3.0d0 * const%pi * const%air_dyn_visc * diameter)

  end function compute_brownian_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  subroutine spec_file_read_stability_class(file, times, classes)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Times (s).
    real(kind=dp), pointer :: times(:)
    !> Rates (s^{-1}).
    real(kind=dp), pointer :: classes(:)

    integer :: n_lines, i, n_time, i_time
    character(len=SPEC_LINE_MAX_VAR_LEN), pointer :: species_name(:)
    real(kind=dp), pointer :: species_data(:,:)

    allocate(species_name(0))
    allocate(species_data(0,0))
    deallocate(times)
    deallocate(classes)
    call spec_file_read_real_named_array(file, 0, species_name, &
         species_data)
    n_time = size(species_data,2)
    allocate(times(n_time))
    allocate(classes(n_time))
    do i_time = 1,n_time
       times(i_time) = species_data(1,i_time)
       classes(i_time) = species_data(2,i_time)
    end do

  end subroutine spec_file_read_stability_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the friction velocity.
  real(kind=dp) function dry_dep_ustar()

     dry_dep_ustar = .10d0

  end function dry_dep_ustar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Computes the aerodynamic resistance.
  real(kind=dp) function dry_dep_aero_resistance()

    dry_dep_aero_resistance  = 0.0d0

  end function dry_dep_aero_resistance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Assigns values for dry deposition given a land type.
  !> Based on Zhang et al 2001 Table 2.
  subroutine dry_dep_land_parameters(gamma, A, alpha, land_type)

    !>
    real(kind=dp), intent(inout) :: gamma
    !>
    real(kind=dp), intent(inout) :: A
    !>
    real(kind=dp), intent(inout) :: alpha
    !>
    integer, intent(in) :: land_type

    integer :: SC

    SC = 1

    A = a_ref(SC,GRASS)
    gamma = gamma_ref(GRASS)
    alpha = alpha_ref(GRASS)

  end subroutine dry_dep_land_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_deposition
