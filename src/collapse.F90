! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_collapse module.

!> Aerosol particle restructuring.
module pmc_collapse

  use pmc_env_state
  use pmc_aero_data
  use pmc_aero_state
  use pmc_gas_data
  use pmc_gas_state
  use pmc_mpi
#ifdef PMC_USE_MPI
  use mpi
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  !> Checks for collapse of fractal black carbon particles based on coating.
  subroutine collapse(env_state, aero_data, aero_state, gas_data, gas_state)

    !> Environmental state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state

    integer :: i_part
    integer :: i_spec
    real(kind=dp), allocatable :: kelvin_length_species(:)
    integer, allocatable :: condensable_index(:,:)
    real(kind=dp), allocatable :: surface_tension(:)
    real(kind=dp), allocatable :: vapor_pressure(:)
    integer :: bc_index, aero_index, gas_index
    real(kind=dp) :: F, zeta, chi, x_f, Q

    ! FIXME: Most of these parameters will be specified in input files
    n_condensable = 2
    allocate(condensable_index(n_condensable,2))
    condensable_index(1,1) = gas_data_spec_by_name(gas_data, "H2SO4")
    condensable_index(1,2) = aero_data_spec_by_name(aero_data, "SO4")
    condensable_index(2,1) = gas_data_spec_by_name(gas_data, "ARO1")
    condensable_index(2,2) = aero_data_spec_by_name(aero_data, "ARO1")
    allocate(surface_tension(n_condensable))
    surface_tension(1) = 72.0d0 / 1000d0 ! N/m
    surface_tension(2) = 72.0d0 / 1000d0 ! N/m
    allocate(vapor_pressure(n_condensable))
    vapor_pressure = 101.25 ! Pa


    ! Compute kelvin length of each species
    allocate(kelvin_length_species(n_condensable))
    do i_spec = 1,n_condensable
      gas_index = condensable_index(i_spec,1)
      aero_index = condensable_index(i_spec,2)
      kelvin_length_species(i_spec) = kelvin_length(surface_tension(i_spec), &
         aero_data%density(aero_index), aero_data%molec_weight(aero_index), &
         env_state%temp)
      vapor_pressure(i_spec) = env_state%pressure &
           * gas_state%mix_rat(gas_index) / 1d9
    end do

    bc_index = aero_data_spec_by_name(aero_data, "BC")
    do i_part = 1,aero_state_n_part(aero_state)
       bc_mass = aero_particle_species_mass(aero_state%apa%particle( &
                     i_part), bc_index, aero_data)
       if (bc_mass > 0d0) then
          do i_spec = 1,n_condensable
             aero_index = condensable_index(i_spec,2)

             x_f = (aero_state%apa%particle(i_part)%vol(aero_index) * &
                  aero_data%density(aero_index) &
                  / aero_data%molec_weight(aero_index)) &
                  / aero_particle_moles(aero_state%apa%particle(i_part), &
                  aero_data)
             zeta = super_sat(vapor_pressure(condensable_index(i_spec,1)), &
                gas_state%sat_vapor_pressure(condensable_index(i_spec,1)), &
                x_f)

             chi = coating_distribution(kelvin_length_species(i_spec), zeta, &
                  aero_state%apa%particle(i_part)%fractal%prime_radius)
             if (chi < 0.6) then ! uniform condensation
                F = 0.1d0
             else
                F = 1.0d0
             end if
             N = fractal_vol_to_num_of_monomers( &
                  aero_state%apa%particle(i_part)%fractal, &
                  aero_state%apa%particle(i_part)%vol(bc_index))
             Q = (aero_particle_species_mass(aero_state%apa%particle( &
                  i_part), aero_index, aero_data) / bc_mass) &
                  * ((2d0 * aero_data%density(bc_index)) &
                  / (N * aero_data%density(aero_index)))
             theta = fill_angle(F, Q)
             ! Test for collapse based on filling angle
             if (theta > 45.0d0) then
                aero_state%apa%particle(i_part)%fractal%frac_dim = 2.8d0
             end if
          end do
       end if
    end do

  end subroutine collapse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determine characteristic Kelvin length.
  !!
  !! Based on Equation 6 in Chen et al [2018].
  real(kind=dp) function kelvin_length(surface_tension, density, molec_weight, &
       temperature)

    !> Surface tension of species.
    real(kind=dp) :: surface_tension
    !> Density of aerosol species.
    real(kind=dp) :: density
    !> Molecular weight of aerosol species.
    real(kind=dp) :: molec_weight
    !> Temperature (K).
    real(kind=dp) :: temperature

    kelvin_length = (2.0d0 * surface_tension * (molec_weight / density)) / &
       (const%univ_gas_const * temperature)

  end function kelvin_length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the supersaturation of a condensable species.
  real(kind=dp) function super_sat(vapor_pressure, sat_vapor_pressure, &
       mole_fraction)

    !> Vapor pressure of species.
    real(kind=dp) :: vapor_pressure
    !> Saturation vapor pressure of species.
    real(kind=dp) :: sat_vapor_pressure
    !> Mole fraction of species.
    real(kind=dp) :: mole_fraction

    super_sat = (vapor_pressure / (sat_vapor_pressure * mole_fraction)) - 1d0

  end function super_sat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  !!
  !! Based on Equation 12 of Chen et al [2018].
  real(kind=dp) function coating_distribution(kelvin_length, super_sat, &
       prime_radius)

    !> Characteristic Kelvin length.
    real(kind=dp) :: kelvin_length
    !> Vapor supersaturation.
    real(kind=dp) :: super_sat
    !> Radius of primary particle.
    real(kind=dp) :: prime_radius

    coating_distribution = kelvin_length / (prime_radius * super_sat)

  end function coating_distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solves for fill angle.
  !!
  !! Based on Equations 14 and 24 of Chen et al [2018].
  real(kind=dp) function fill_angle(F,Q)

    !>
    real(kind=dp) :: F
    !>
    real(kind=dp) :: Q

    real(kind=dp), parameter :: NEWTON_REL_TOL = 1d-14
    integer, parameter :: NEWTON_MAX_STEPS = 10

    ! Set initial guess
    x = 0.01
    do newton_step = 1,NEWTON_MAX_STEPS
       f = (2.*cos(x)**3 - 3.0*cos(x)**2 + 1 - F*Q)
       df = -6.0*cos(x)**2*sin(x) + 3*sin(2*x)
       x = x - f / df
       if (abs(f / df) / (abs(x + f / df) + abs(x)) &
            < NEWTON_REL_TOL) exit
    end do

    call assert_msg(589404526, newton_step < NEWTON_MAX_STEPS, &
         "collapse Newton loop failed to converge")

    fill_angle = x

  end function fill_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_collapse
