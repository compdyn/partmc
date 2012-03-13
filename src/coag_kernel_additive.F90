! Copyright (C) 2005-2011 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_coag_kernel_additive module.

!> Additive coagulation kernel.
module pmc_coag_kernel_additive

  use pmc_bin_grid
  use pmc_env_state
  use pmc_util
  use pmc_constants
  use pmc_constants
  use pmc_aero_binned
  use pmc_aero_data
  use pmc_aero_dist
  use pmc_aero_data
  use pmc_aero_particle
  
  !> Scaling coefficient for constant kernel.
  real(kind=dp), parameter :: beta_1 = 1000d0

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Additive coagulation kernel.
  subroutine kernel_additive(aero_particle_1, aero_particle_2, &
       aero_data, env_state, k)

    !> First particle.
    type(aero_particle_t), intent(in) :: aero_particle_1
    !> Second particle.
    type(aero_particle_t), intent(in) :: aero_particle_2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel.
    real(kind=dp), intent(out) :: k

    k = beta_1 * (aero_particle_volume(aero_particle_1) &
         + aero_particle_volume(aero_particle_2))

  end subroutine kernel_additive
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Minimum and maximum values of the additive kernel.
  subroutine kernel_additive_minmax(v1, v2, aero_data, env_state, k_min, k_max)

    !> Volume of first particle.
    real(kind=dp), intent(in) :: v1
    !> Volume of second particle.
    real(kind=dp), intent(in) :: v2
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Coagulation kernel minimum value.
    real(kind=dp), intent(out) :: k_min
    !> Coagulation kernel maximum value.
    real(kind=dp), intent(out) :: k_max
    
    k_min = beta_1 * (v1 + v2)
    k_max = k_min
    
  end subroutine kernel_additive_minmax
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exact solution with the additive coagulation kernel and
  !> exponential initial condition.
  !!
  !! Given input paramaters \f$R\f$ and \f$N_0\f$ we let the mean
  !! volume be \f$v_\mu = \frac{4\pi}{3} R^3\f$ and define the
  !! rescaled times \f$\tau = N_0 v_\mu \beta_1 t\f$ and \f$T = 1 -
  !! e^{-\tau}\f$, where \f$\beta_1\f$ is the fixed kernel scaling
  !! parameter. Then the solution is
  !! \f[
  !!     n(D,t) \ {\rm d}\ln D
  !!     = \frac{\pi}{2} D^3 
  !!       \frac{N_0}{v} \frac{1 - T}{\sqrt{T}}
  !!       \exp\left(-(1 + T) \frac{v}{v_\mu}\right)
  !!       I_1\left(2 \frac{v}{v_\mu} \sqrt{T}\right) {\rm d}\ln D
  !! \f]
  !! where \f$I_1(x)\f$ is the <a
  !! href="http://en.wikipedia.org/wiki/Bessel_function">modified
  !! Bessel function of the first kind</a> and \f$v = \frac{\pi}{6}
  !! D^3\f$.
  !!
  !! For small \f$x\f$ we have \f$I_1(x) \approx \frac{x}{2}\f$, so
  !! this solution has initial condition
  !! \f[
  !!     n(D,t) \ {\rm d}\ln D
  !!     = \frac{\pi}{2} D^3 \frac{N_0}{v_\mu}
  !!     \exp\left(-\frac{v}{v_\mu}\right) {\rm d}\ln D
  !! \f]
  subroutine soln_additive_exp(bin_grid, aero_data, time, num_conc, &
       radius_at_mean_vol, env_state, aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Current time.
    real(kind=dp), intent(in) :: time
    !> Particle number concentration (#/m^3).
    real(kind=dp), intent(in) :: num_conc
    !> Mean init radius (m).
    real(kind=dp), intent(in) :: radius_at_mean_vol
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Output state.
    type(aero_binned_t), intent(inout) :: aero_binned
    
    real(kind=dp) :: tau, T, rat_v, nn, b, x, mean_vol
    integer :: k
    
    mean_vol = rad2vol(radius_at_mean_vol, aero_data%fractal)
    if (time .eq. 0d0) then
       do k = 1,bin_grid%n_bin
          aero_binned%num_conc(k) = const%pi/2d0 &
               * (2d0 * bin_grid%center_radius(k))**3 * num_conc / mean_vol &
               * exp(-(rad2vol(bin_grid%center_radius(k), aero_data%fractal) &
               / mean_vol))
       end do
    else
       tau = num_conc * mean_vol * beta_1 * time
       T = 1d0 - exp(-tau)
       do k = 1,bin_grid%n_bin
          rat_v = rad2vol(bin_grid%center_radius(k), aero_data%fractal) / mean_vol
          x = 2d0 * rat_v * sqrt(T)
          if (x .lt. 500d0) then
             call bessi1(x, b)
             nn = num_conc / rad2vol(bin_grid%center_radius(k), aero_data%fractal) &
                  * (1d0 - T) / sqrt(T) * exp(-((1d0 + T) * rat_v)) * b
          else
             ! For very large volumes we can use the asymptotic
             ! approximation I_1(x) \approx e^x / sqrt(2 pi x) and
             ! simplify the result to avoid the overflow from
             ! multiplying a huge bessel function result by a very
             ! tiny exponential.
             nn = num_conc / rad2vol(bin_grid%center_radius(k), aero_data%fractal) &
                  * (1d0 - T) / sqrt(T) &
                  * exp((2d0*sqrt(T) - T - 1d0) * rat_v) &
                  / sqrt(4d0 * const%pi * rat_v * sqrt(T))
          end if
          aero_binned%num_conc(k) = const%pi/2d0 &
               * (2d0 * bin_grid%center_radius(k))**3 * nn
       end do
    end if

    aero_binned%vol_conc = 0d0
    do k = 1,bin_grid%n_bin
       aero_binned%vol_conc(k,1) = rad2vol(bin_grid%center_radius(k), &
            aero_data%fractal) * aero_binned%num_conc(k)
    end do
    
  end subroutine soln_additive_exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Modified Bessel function of the first kind \f$ I_1(x) \f$.
  !!
  !! This looks like it was taken from Numerical Recipes.
  !!
  !! FIXME: Where did this code come from? What license does it have?
  subroutine bessi1(x, r)

    !> Function argument.
    real(kind=dp), intent(in) :: x
    !> Function value.
    real(kind=dp), intent(out) :: r
    
    real(kind=dp) ax
    real(kind=dp) p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
    data p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0, &
         0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
    data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1, &
         -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1, &
         -0.2895312d-1,0.1787654d-1,-0.420059d-2/
    
    if (abs(x) .lt. 3.75d0) then
       y = (x / 3.75d0)**2
       r = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    elseif (abs(x) .lt. 500d0) then
       ax = abs(x)
       y = 3.75d0 / ax
       r = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+ &
            y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
       if (x .lt. 0d0) r = -r
    else
       
    end if
    
  end subroutine bessi1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_coag_kernel_additive
