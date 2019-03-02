! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_collapse module.

!> Aerosol particle restructuring.
!!
!! References
!!   - C.&nbsp;Chen,O.&nbsp;Y.&nbsp;Enekwizu, X.&nbsp;Fan, 
!!     C.&nbsp;D.&nbsp;Dobrzanski, E.&nbsp;V.&nbsp;Ivanova, Y.&nbsp;Ma,
!!     and A.&nbsp;F.&nbsp;Khalizov (2018) Single Parameter for Predicting the
!!     Morphology of Atmospheric Black Carbon, <i>Environmental Science \&
!!     Technology </i>, 52(24), 14169-14179. DOI: <a
!!     href="http://dx.doi.org/10.1021/acs.est.8b04201"> 
!!     10.1021/acs.est.8b04201 </a>
module pmc_collapse

  use pmc_constants
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
    integer, allocatable :: condensable_gas_index(:)
    integer, allocatable :: condensable_aero_index(:)
    real(kind=dp), allocatable :: surface_tension(:)
    real(kind=dp), allocatable :: vapor_pressure(:)
    integer :: bc_index, i_aero, i_gas 
    real(kind=dp) :: F, zeta, chi, x_f, Q

    ! FIXME: Most of these parameters will be specified in input files
    n_condensable = 2 ! Number of condensing species
    allocate(condensable_gas_index(n_condensable))
    allocate(condensable_aero_index(n_condensable))
    condensable_gas_index(1) = gas_data_spec_by_name(gas_data, "H2SO4")
    condensable_aero_index(1) = aero_data_spec_by_name(aero_data, "SO4")
    condensable_gas_index(2) = gas_data_spec_by_name(gas_data, "ARO1")
    condensable_aero_index(2) = aero_data_spec_by_name(aero_data, "ARO1")
    allocate(surface_tension(n_condensable))
    surface_tension(1) = 72.0d0 / 1000d0 ! N/m
    surface_tension(2) = 72.0d0 / 1000d0 ! N/m

    ! Compute kelvin length of each species
    allocate(vapor_pressure(n_condensable))
    allocate(kelvin_length_species(n_condensable))
    do i_cond_spec = 1,n_condensable
      i_gas = condensable_gas_index(i_cond_spec)
      i_aero = condensable_aero_index(i_cond_spec)
      kelvin_length_species(i_cond_spec) = kelvin_length(surface_tension( &
           i_cond_spec), aero_data%density(i_aero), &
           aero_data%molec_weight(i_aero), env_state%temp)
      vapor_pressure(i_cond_spec) = env_state%pressure &
           * gas_state%mix_rat(i_gas) / 1d9
    end do

    bc_index = aero_data_spec_by_name(aero_data, "BC")
    do i_part = 1,aero_state_n_part(aero_state)
       bc_mass = aero_particle_species_mass(aero_state%apa%particle( &
                     i_part), bc_index, aero_data)
       if (bc_mass > 0d0) then
          do i_cond_spec = 1,n_condensable
             i_aero = condensable_aero_index(i_cond_spec)
             i_gas = condensable_gas_index(i_cond_spec)
             x_f = (aero_state%apa%particle(i_part)%vol(i_aero) * &
                  aero_data%density(i_aero) &
                  / aero_data%molec_weight(i_aero)) &
                  / aero_particle_moles(aero_state%apa%particle(i_part), &
                  aero_data)
             zeta = super_sat(vapor_pressure(i_cond_spec), &
                gas_state%sat_vapor_pressure(i_gas), &
                x_f)

             chi = coating_distribution(kelvin_length_species(i_cond_spec), &
                  zeta, aero_state%apa%particle(i_part)%fractal%prime_radius)
             if (chi < 0.6) then ! uniform condensation
                F = 0.1d0
             else
                F = 1.0d0
             end if
             ! Possible to have particles that are too small
             N = max(fractal_vol_to_num_of_monomers( &
                  aero_state%apa%particle(i_part)%fractal, &
                  aero_state%apa%particle(i_part)%vol(bc_index)),1d0)
             Q = (aero_particle_species_mass(aero_state%apa%particle( &
                  i_part), i_aero, aero_data) / bc_mass) &
                  * ((2d0 * aero_data%density(bc_index)) &
                  / (N * aero_data%density(i_aero)))
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
  real(kind=dp) function kelvin_length(surface_tension, density, &
       molec_weight, temperature)

    !> Surface tension of species.
    real(kind=dp), intent(in) :: surface_tension
    !> Density of aerosol species.
    real(kind=dp), intent(in) :: density
    !> Molecular weight of aerosol species.
    real(kind=dp), intent(in) :: molec_weight
    !> Temperature (K).
    real(kind=dp), intent(in) :: temperature

    kelvin_length = (2.0d0 * surface_tension * (molec_weight / density)) / &
       (const%univ_gas_const * temperature)

  end function kelvin_length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the supersaturation of a condensable species.
  real(kind=dp) function super_sat(vapor_pressure, sat_vapor_pressure, &
       mole_fraction)

    !> Vapor pressure of species.
    real(kind=dp), intent(in) :: vapor_pressure
    !> Saturation vapor pressure of species.
    real(kind=dp), intent(in) :: sat_vapor_pressure
    !> Mole fraction of species.
    real(kind=dp), intent(in) :: mole_fraction

    super_sat = (vapor_pressure / (sat_vapor_pressure * mole_fraction)) - 1d0

  end function super_sat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  !!
  !! Based on Equation 12 of Chen et al [2018].
  real(kind=dp) function coating_distribution(kelvin_length, super_sat, &
       prime_radius)

    !> Characteristic Kelvin length.
    real(kind=dp), intent(in) :: kelvin_length
    !> Vapor supersaturation.
    real(kind=dp), intent(in) :: super_sat
    !> Radius of primary particle.
    real(kind=dp), intent(in) :: prime_radius

    coating_distribution = kelvin_length / (prime_radius * super_sat)

  end function coating_distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Solves for fill angle.
  !!
  !! Based on Equations 14 and 24 of Chen et al [2018].
  real(kind=dp) function fill_angle(F_val,Q_val)

    !>
    real(kind=dp), intent(in) :: F_val
    !>
    real(kind=dp), intent(in) :: Q_val

    real(kind=dp) :: f, df, x
    real(kind=dp), parameter :: NEWTON_REL_TOL = 1d-3
    integer, parameter :: NEWTON_MAX_STEPS = 1000

    ! Set initial guess
    x = 0d0 
    do newton_step = 1,NEWTON_MAX_STEPS
       f = (2.*cos(x)**3 - 3.0*cos(x)**2 + 1 - F_val *Q_val)
       df = -6.0*cos(x)**2*sin(x) + 3*sin(2*x)
       x = x - f / df
       if (abs(f / df) / (abs(x + f / df) + abs(x)) &
            < NEWTON_REL_TOL) exit
    end do

    call assert_msg(589404526, newton_step < NEWTON_MAX_STEPS, &
         "collapse Newton loop failed to converge")

    fill_angle = x * 180d0 / const%pi

  end function fill_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_collapse
