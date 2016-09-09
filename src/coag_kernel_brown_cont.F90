! Copyright (C) 2012 Jian Tian
! Copyright (C) 2005-2011 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coag_kernel_brown_cont module.
!!
!! The coagulation kernel is based on Eq. 6 in
!! S. Vemury and S. E. Pratsinis, Self-preserving size distributions
!! of agglomerates, Journal of Aerosol Science, Vol. 26, No. 2,
!! pp. 175-185, 1995.

!> Brownian coagulation kernel in continuum regime based on
!> Vemury and Pratsinis [1995].
module pmc_coag_kernel_brown_cont

  use pmc_env_state
  use pmc_constants
  use pmc_util
  use pmc_aero_particle
  use pmc_aero_data

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the Brownian coagulation kernel in continuum regime.
  !!
  !! Use Eq. 6 of Vemury and Pratsinis [1995].
  subroutine kernel_brown_cont(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)

    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Kernel k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: k

    real(kind=dp) :: v1, v2, d1, d2

    v1 = aero_particle_volume(aero_particle_1)
    v2 = aero_particle_volume(aero_particle_2)
    d1 = aero_particle_density(aero_particle_1, aero_data)
    d2 = aero_particle_density(aero_particle_2, aero_data)

    call kernel_brown_cont_helper(v1, d1, v2, d2, aero_data, &
         env_state%temp, k)

  end subroutine kernel_brown_cont

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the minimum and maximum Brownian coagulation kernel in continuum
  !> regime based on Vemury and Pratsinis [1995].
  !!
  !! Finds the minimum and maximum kernel values between particles of
  !! volumes v1 and v2, by sampling over possible densities.
  subroutine kernel_brown_cont_minmax(v1, v2, aero_data, env_state, &
       k_min, k_max)

    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Minimum kernel value (m^3/s).
    real(kind=dp), intent(out) :: k_min
    !> Maximum kernel value (m^3/s).
    real(kind=dp), intent(out) :: k_max

    !> Number of density sample points.
    integer, parameter :: n_sample = 3

    real(kind=dp) :: d1, d2, d_min, d_max, k
    integer :: i, j
    logical :: first

    d_min = minval(aero_data%density)
    d_max = maxval(aero_data%density)

    first = .true.
    do i = 1,n_sample
       do j = 1,n_sample
          d1 = interp_linear_disc(d_min, d_max, n_sample, i)
          d2 = interp_linear_disc(d_min, d_max, n_sample, j)
          call kernel_brown_cont_helper(v1, d1, v2, d2, aero_data, &
               env_state%temp, k)
          if (first) then
             first = .false.
             k_min = k
             k_max = k
          else
             k_min = min(k_min, k)
             k_max = max(k_max, k)
          end if
       end do
    end do

  end subroutine kernel_brown_cont_minmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Helper function that does the actual coagulation kernel computation.
  !!
  !! Helper function. Do not call directly. Instead use kernel_brown_cont().
  !!
  !! Use Eq. 6 of Vemury and Pratsinis [1995].
  subroutine kernel_brown_cont_helper(v1, d1, v2, d2, aero_data, &
       temp, bckernel)

    !> Volume of first particle (m^3).
    real(kind=dp), intent(in) :: v1
    !> Density of first particle (kg/m^3).
    real(kind=dp), intent(in) :: d1
    !> Volume of second particle (m^3).
    real(kind=dp), intent(in) :: v2
    !> Density of second particle (kg/m^3).
    real(kind=dp), intent(in) :: d2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Kernel k(a,b) (m^3/s).
    real(kind=dp), intent(out) :: bckernel

    real(kind=dp) :: N_i, N_j, N_i_inv_df, N_j_inv_df

    N_i = aero_data_vol_to_num_of_monomers(aero_data, v1)
    N_j = aero_data_vol_to_num_of_monomers(aero_data, v2)
    N_i_inv_df = N_i**(1d0 / aero_data%fractal%frac_dim)
    N_j_inv_df = N_j**(1d0 / aero_data%fractal%frac_dim)
    bckernel = 2d0 * const%boltzmann * temp / 3d0 / const%air_dyn_visc &
         * (1d0 / N_i_inv_df + 1d0 / N_j_inv_df) * (N_i_inv_df + N_j_inv_df)

  end subroutine kernel_brown_cont_helper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coag_kernel_brown_cont
