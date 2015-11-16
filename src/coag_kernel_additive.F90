! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
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

    call aero_binned_set_sizes(aero_binned, bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))

    mean_vol = aero_data_rad2vol(aero_data, radius_at_mean_vol)
    if (time .eq. 0d0) then
       do k = 1,bin_grid_size(bin_grid)
          aero_binned%num_conc(k) = const%pi/2d0 &
               * (2d0 * bin_grid%centers(k))**3 * num_conc / mean_vol &
               * exp(-(aero_data_rad2vol(aero_data, bin_grid%centers(k)) &
               / mean_vol))
       end do
    else
       tau = num_conc * mean_vol * beta_1 * time
       T = 1d0 - exp(-tau)
       do k = 1,bin_grid_size(bin_grid)
          rat_v = aero_data_rad2vol(aero_data, bin_grid%centers(k)) / mean_vol
          x = 2d0 * rat_v * sqrt(T)
          if (x .lt. 500d0) then
             call bessi1(x, b)
             nn = num_conc / aero_data_rad2vol(aero_data, &
                  bin_grid%centers(k)) &
                  * (1d0 - T) / sqrt(T) * exp(-((1d0 + T) * rat_v)) * b
          else
             ! For very large volumes we can use the asymptotic
             ! approximation I_1(x) \approx e^x / sqrt(2 pi x) and
             ! simplify the result to avoid the overflow from
             ! multiplying a huge bessel function result by a very
             ! tiny exponential.
             nn = num_conc / aero_data_rad2vol(aero_data, &
                  bin_grid%centers(k)) &
                  * (1d0 - T) / sqrt(T) &
                  * exp((2d0*sqrt(T) - T - 1d0) * rat_v) &
                  / sqrt(4d0 * const%pi * rat_v * sqrt(T))
          end if
          aero_binned%num_conc(k) = const%pi/2d0 &
               * (2d0 * bin_grid%centers(k))**3 * nn
       end do
    end if

    aero_binned%vol_conc = 0d0
    do k = 1,bin_grid_size(bin_grid)
       aero_binned%vol_conc(k,1) = aero_data_rad2vol(aero_data, &
            bin_grid%centers(k)) * aero_binned%num_conc(k)
    end do

  end subroutine soln_additive_exp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Modified Bessel function of the first kind \f$ I_1(x) \f$.
  subroutine bessi1(x, r)

    !> Function argument.
    real(kind=dp), intent(in) :: x
    !> Function value.
    real(kind=dp), intent(out) :: r

    call calci1 (x, r, 1 )

  end subroutine bessi1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculates modified Bessel functions of the first kind \f$ I_1(x) \f$.
  subroutine calci1 ( arg, result, jint )

    !*************************************************************************
    !
    !! CALCI1 computes various I1 Bessel functions.
    !
    !  Discussion:
    !
    !    This routine computes modified Bessel functioons of the first kind
    !    and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
    !    arguments X.
    !
    !    The main computation evaluates slightly modified forms of
    !    minimax approximations generated by Blair and Edwards, Chalk
    !    River (Atomic Energy of Canada Limited) Report AECL-4928,
    !    October, 1974.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody, Laura Stoltz.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
    !    the argument must be less than XMAX.
    !
    !    Output, real ( kind = 8 ) RESULT, the value of the function,
    !    which depends on the input value of JINT:
    !    1, RESULT = I1(x);
    !    2, RESULT = exp(-x) * I1(x);
    !
    !    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
    !    1, I1(x);
    !    2, exp(-x) * I1(x);
    !

    real ( kind = 8 ) a
    real ( kind = 8 ) arg
    real ( kind = 8 ) b
    real ( kind = 8 ) exp40
    real ( kind = 8 ) forty
    integer ( kind = 4 ) j
    integer ( kind = 4 ) jint
    real ( kind = 8 ) one5
    real ( kind = 8 ) p(15)
    real ( kind = 8 ) pbar
    real ( kind = 8 ) pp(8)
    real ( kind = 8 ) q(5)
    real ( kind = 8 ) qq(6)
    real ( kind = 8 ) rec15
    real ( kind = 8 ) result
    real ( kind = 8 ) sump
    real ( kind = 8 ) sumq
    real ( kind = 8 ) two25
    real ( kind = 8 ) x
    real ( kind = 8 ) xinf
    real ( kind = 8 ) xmax
    real ( kind = 8 ) xsmall
    real ( kind = 8 ) xx
    !
    !  Mathematical constants
    !
    data one5 / 15.0d0 /
    data exp40 / 2.353852668370199854d17 /
    data forty / 40.0d0 /
    data rec15 / 6.6666666666666666666d-2 /
    data two25 / 225.0d0 /
    !
    !  Machine-dependent constants
    !
    data xsmall /5.55d-17/
    data xinf /1.79d308/
    data xmax /713.987d0/
    !
    !  Coefficients for XSMALL <= ABS(ARG) < 15.0
    !
    data p/-1.9705291802535139930d-19,-6.5245515583151902910d-16, &
           -1.1928788903603238754d-12,-1.4831904935994647675d-09, &
           -1.3466829827635152875d-06,-9.1746443287817501309d-04, &
           -4.7207090827310162436d-01,-1.8225946631657315931d+02, &
           -5.1894091982308017540d+04,-1.0588550724769347106d+07, &
           -1.4828267606612366099d+09,-1.3357437682275493024d+11, &
           -6.9876779648010090070d+12,-1.7732037840791591320d+14, &
           -1.4577180278143463643d+15/
    data q/-4.0076864679904189921d+03, 7.4810580356655069138d+06, &
           -8.0059518998619764991d+09, 4.8544714258273622913d+12, &
           -1.3218168307321442305d+15/
    !
    !  Coefficients for 15.0 <= ABS(ARG)
    !
    data pp/-6.0437159056137600000d-02, 4.5748122901933459000d-01, &
            -4.2843766903304806403d-01, 9.7356000150886612134d-02, &
            -3.2457723974465568321d-03,-3.6395264712121795296d-04, &
             1.6258661867440836395d-05,-3.6347578404608223492d-07/
    data qq/-3.8806586721556593450d+00, 3.2593714889036996297d+00, &
            -8.5017476463217924408d-01, 7.4212010813186530069d-02, &
            -2.2835624489492512649d-03, 3.7510433111922824643d-05/
    data pbar/3.98437500d-01/

    x = abs ( arg )
    !
    !  Return for ABS(ARG) < XSMALL.
    !
    if ( x < xsmall ) then

       result = 0.5D+00 * x
    !
    !  XSMALL <= ABS(ARG) < 15.0.
    !
    else if ( x < one5 ) then

       xx = x * x
       sump = p(1)
       do j = 2, 15
          sump = sump * xx + p(j)
       end do
       xx = xx - two25

       sumq = (((( &
            xx + q(1) ) &
            * xx + q(2) ) &
            * xx + q(3) ) &
            * xx + q(4) ) &
            * xx + q(5)

       result = ( sump / sumq ) * x

       if ( jint == 2 ) then
          result = result * exp ( -x )
       end if

    else if ( jint == 1 .and. xmax < x ) then

       result = xinf

    else
       !
       !  15.0 <= ABS(ARG).
       !
       xx = 1.0D+00 / x - rec15

       sump = (((((( &
                   pp(1) &
            * xx + pp(2) ) &
            * xx + pp(3) ) &
            * xx + pp(4) ) &
            * xx + pp(5) ) &
            * xx + pp(6) ) &
            * xx + pp(7) ) &
            * xx + pp(8)

       sumq = ((((( &
              xx + qq(1) ) &
            * xx + qq(2) ) &
            * xx + qq(3) ) &
            * xx + qq(4) ) &
            * xx + qq(5) ) &
            * xx + qq(6)

       result = sump / sumq

    if ( jint /= 1 ) then
       result = ( result + pbar ) / sqrt ( x )
    else
       !
       !  Calculation reformulated to avoid premature overflow.
       !
       if ( xmax - one5 < x ) then
          a = exp ( x - forty )
          b = exp40
       else
          a = exp ( x )
          b = 1.0D+00
       end if

       result = ( ( result * a + pbar * a ) / sqrt ( x ) ) * b

       end if
    end if

    if ( arg < 0.0D+00 ) then
       result = -result
    end if

    return

  end subroutine calci1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_coag_kernel_additive
