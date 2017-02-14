! Copyright (C) 2012-2017 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_stats module.

!> The \c stats_t type and associated subroutines.
module pmc_stats

  use pmc_util
  use pmc_netcdf

  !> Structure for online computation of mean and variance.
  type stats_t
     !> Number of samples included so far.
     integer :: n = 0
     !> Current mean estimate.
     real(kind=dp) :: mean
     !> Current variance estimate.
     real(kind=dp) :: var
  end type stats_t

  !> Structure for online computation of 1D arrays of mean and variance.
  type stats_1d_t
     !> Number of samples included so far.
     integer, allocatable :: n(:)
     !> Current mean estimates.
     real(kind=dp), allocatable :: mean(:)
     !> Current variance estimates.
     real(kind=dp), allocatable :: var(:)
  end type stats_1d_t

  !> Structure for online computation of 2D arrays of mean and variance.
  type stats_2d_t
     !> Number of samples included so far.
     integer, allocatable :: n(:, :)
     !> Current mean estimates.
     real(kind=dp), allocatable :: mean(:, :)
     !> Current variance estimates.
     real(kind=dp), allocatable :: var(:, :)
  end type stats_2d_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clear data statistics collected so far.
  subroutine stats_clear(stats)

    !> Statistics structure to clear.
    type(stats_t), intent(inout) :: stats

    stats%n = 0

  end subroutine stats_clear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clear data statistics collected so far.
  subroutine stats_1d_clear(stats)

    !> Statistics structure to clear.
    type(stats_1d_t), intent(inout) :: stats

    if (allocated(stats%n)) then
       deallocate(stats%n)
       deallocate(stats%mean)
       deallocate(stats%var)
    end if

  end subroutine stats_1d_clear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Clear data statistics collected so far.
  subroutine stats_2d_clear(stats)

    !> Statistics structure to clear.
    type(stats_2d_t), intent(inout) :: stats

    if (allocated(stats%n)) then
       deallocate(stats%n)
       deallocate(stats%mean)
       deallocate(stats%var)
    end if

  end subroutine stats_2d_clear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new data value to a \c stats_t structure.
  subroutine stats_add(stats, data)

    !> Statistics structure to add to.
    type(stats_t), intent(inout) :: stats
    !> Data value to add.
    real(kind=dp), intent(in) :: data

    stats%n = stats%n + 1
    call update_mean_var(stats%mean, stats%var, data, stats%n)

  end subroutine stats_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add all new data values to a \c stats_1d_t structure.
  subroutine stats_1d_add(stats, data)

    !> Statistics structure to add to.
    type(stats_1d_t), intent(inout) :: stats
    !> Data values to add.
    real(kind=dp), intent(in) :: data(:)

    integer :: i

    if (.not. allocated(stats%n)) then
       allocate(stats%n(size(data)))
       allocate(stats%mean(size(data)))
       allocate(stats%var(size(data)))
       stats%n = 0
    end if

    call assert_msg(851913829, size(stats%n) == size(data), &
         "size mismatch between existing n and newly added data")
    call assert_msg(933021218, size(stats%mean) == size(data), &
         "size mismatch between existing mean and newly added data")
    call assert_msg(170092922, size(stats%var) == size(data), &
         "size mismatch between existing var and newly added data")
    stats%n = stats%n + 1
    do i = 1,size(data)
       call update_mean_var(stats%mean(i), stats%var(i), data(i), stats%n(i))
    end do

  end subroutine stats_1d_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a new single data value to a \c stats_1d_t structure.
  subroutine stats_1d_add_entry(stats, data, i)

    !> Statistics structure to add to.
    type(stats_1d_t), intent(inout) :: stats
    !> Data value to add.
    real(kind=dp), intent(in) :: data
    !> Index of data value to add.
    integer, intent(in) :: i

    call assert_msg(802003511, i >= 1, "cannot use a non-positive index")
    call ensure_integer_array_size(stats%n, i, only_grow=.true.)
    call ensure_real_array_size(stats%mean, i, only_grow=.true.)
    call ensure_real_array_size(stats%var, i, only_grow=.true.)

    stats%n(i) = stats%n(i) + 1
    call update_mean_var(stats%mean(i), stats%var(i), data, stats%n(i))

  end subroutine stats_1d_add_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add all new data values to a \c stats_2d_t structure.
  subroutine stats_2d_add(stats, data)

    !> Statistics structure to add to.
    type(stats_2d_t), intent(inout) :: stats
    !> Data values to add.
    real(kind=dp), intent(in) :: data(:, :)

    integer :: i, j

    if (.not. allocated(stats%n)) then
       allocate(stats%n(size(data, 1), size(data, 2)))
       allocate(stats%mean(size(data, 1), size(data, 2)))
       allocate(stats%var(size(data, 1), size(data, 2)))
       stats%n = 0
    else
       call assert_msg(563659350, all(shape(stats%n) == shape(data)), &
            "size mismatch between existing n and newly added data")
       call assert_msg(274220601, all(shape(stats%mean) == shape(data)), &
            "size mismatch between existing mean and newly added data")
       call assert_msg(939179895, all(shape(stats%var) == shape(data)), &
            "size mismatch between existing var and newly added data")
    end if
    stats%n = stats%n + 1
    do i = 1,size(data, 1)
       do j = 1,size(data, 2)
          call update_mean_var(stats%mean(i, j), stats%var(i, j), data(i, j), &
               stats%n(i, j))
       end do
    end do

  end subroutine stats_2d_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a row of new data values to a \c stats_2d_t structure.
  subroutine stats_2d_add_row(stats, data, i)

    !> Statistics structure to add to.
    type(stats_2d_t), intent(inout) :: stats
    !> Data values to add.
    real(kind=dp), intent(in) :: data(:)
    !> Row of data value to add.
    integer, intent(in) :: i

    integer :: j

    call assert_msg(549391523, i >= 1, "cannot use a non-positive row")
    if (allocated(stats%n)) then
       call assert_msg(286470660, size(stats%n, 2) == size(data), &
            "size mismatch between existing n and newly added data")
       call assert_msg(901102174, size(stats%mean, 2) == size(data), &
            "size mismatch between existing mean and newly added data")
       call assert_msg(993806885, size(stats%var, 2) == size(data), &
            "size mismatch between existing var and newly added data")
    end if
    call ensure_integer_array_2d_size(stats%n, i, size(stats%n, 2), &
         only_grow=.true.)
    call ensure_real_array_2d_size(stats%mean, i, size(stats%mean, 2), &
         only_grow=.true.)
    call ensure_real_array_2d_size(stats%var, i, size(stats%var, 2), &
         only_grow=.true.)

    do j = 1,size(stats%n, 2)
       stats%n(i, j) = stats%n(i, j) + 1
       call update_mean_var(stats%mean(i, j), stats%var(i, j), data(j), &
            stats%n(i, j))
    end do

  end subroutine stats_2d_add_row

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a column of new data values to a \c stats_2d_t structure.
  subroutine stats_2d_add_col(stats, data, j)

    !> Statistics structure to add to.
    type(stats_2d_t), intent(inout) :: stats
    !> Data values to add.
    real(kind=dp), intent(in) :: data(:)
    !> Column of data value to add.
    integer, intent(in) :: j

    integer :: i

    call assert_msg(944718123, j >= 1, "cannot use a non-positive column")
    if (allocated(stats%n)) then
       call assert_msg(170275693, size(stats%n, 1) == size(data), &
            "size mismatch between existing n and newly added data")
       call assert_msg(659257452, size(stats%mean, 1) == size(data), &
            "size mismatch between existing mean and newly added data")
       call assert_msg(279552980, size(stats%var, 1) == size(data), &
            "size mismatch between existing var and newly added data")
    end if
    call ensure_integer_array_2d_size(stats%n, size(stats%n, 1), j, &
         only_grow=.true.)
    call ensure_real_array_2d_size(stats%mean, size(stats%mean, 1), j, &
         only_grow=.true.)
    call ensure_real_array_2d_size(stats%var, size(stats%var, 1), j, &
         only_grow=.true.)

    do i = 1,size(stats%n, 1)
       stats%n(i, j) = stats%n(i, j) + 1
       call update_mean_var(stats%mean(i, j), stats%var(i, j), data(i), &
            stats%n(i, j))
    end do

  end subroutine stats_2d_add_col

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Add a single new data value to a \c stats_2d_t structure.
  subroutine stats_2d_add_entry(stats, data, i, j)

    !> Statistics structure to add to.
    type(stats_2d_t), intent(inout) :: stats
    !> Data values to add.
    real(kind=dp), intent(in) :: data
    !> First index of data value to add.
    integer, intent(in) :: i
    !> Second index of data value to add.
    integer, intent(in) :: j

    call assert_msg(522548347, (i >= 1) .and. (j >= 1), &
         "cannot use non-positive indexes")
    call ensure_integer_array_2d_size(stats%n, i, j)
    call ensure_real_array_2d_size(stats%mean, i, j)
    call ensure_real_array_2d_size(stats%var, i, j)

    stats%n(i, j) = stats%n(i, j) + 1
    call update_mean_var(stats%mean(i, j), stats%var(i, j), data, &
         stats%n(i, j))

  end subroutine stats_2d_add_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the 95% confidence interval offset from the mean.
  function stats_conf_95_offset(stats)

    !> Statistics structure to find confidence interval for.
    type(stats_t), intent(in) :: stats

    !> Return offset so that [mean - offset, mean + offset] is the 95%
    !> confidence interval.
    real(kind=dp) :: stats_conf_95_offset

    stats_conf_95_offset = conf_95_offset(stats%var, stats%n)

  end function stats_conf_95_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the 95% confidence interval offset from the mean.
  function stats_1d_conf_95_offset(stats)

    !> Statistics structure to find confidence interval for.
    type(stats_1d_t), intent(in) :: stats

    !> Return offset so that [mean - offset, mean + offset] is the 95%
    !> confidence interval.
    real(kind=dp) :: stats_1d_conf_95_offset(size(stats%n))

    integer :: i

    do i = 1,size(stats%n)
       stats_1d_conf_95_offset(i) = conf_95_offset(stats%var(i), stats%n(i))
    end do

  end function stats_1d_conf_95_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute the 95% confidence interval offset from the mean.
  function stats_2d_conf_95_offset(stats)

    !> Statistics structure to find confidence interval for.
    type(stats_2d_t), intent(in) :: stats

    !> Return offset so that [mean - offset, mean + offset] is the 95%
    !> confidence interval.
    real(kind=dp) :: stats_2d_conf_95_offset(size(stats%n, 1), &
         size(stats%n, 2))

    integer :: i, j

    do i = 1,size(stats%n, 1)
       do j = 1,size(stats%n, 2)
          stats_2d_conf_95_offset(i, j) = conf_95_offset(stats%var(i, j), &
               stats%n(i, j))
       end do
    end do

  end function stats_2d_conf_95_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a running average and variance.
  !!
  !! Given a sequence of data <tt>x(i)</tt> for <tt>i = 1,...,n</tt>,
  !! this should be called like
  !! <pre>
  !! do i = 1,n
  !!   call update_mean_var(mean, var, x(i), i)
  !! end do
  !! </pre>
  !! After each call the variables \c mean and \c var will be the
  !! sample mean and sample variance of the sequence elements up to \c
  !! i.
  !!
  !! This computes the sample mean and sample variance using a
  !! recurrence. The initial sample mean is \f$m_1 = x_1\f$ and the
  !! initial sample variance is \f$v_1 = 0\f$ for \f$n = 1\f$, and
  !! then for \f$n \ge 2\f$ we use the mean update
  !! \f[
  !!     m_n = m_{n-1} + \frac{x_n - m_{n-1}}{n}
  !! \f]
  !! and the variance update
  !! \f[
  !!     v_n = \frac{n - 2}{n - 1} v_{n-1}
  !!           + \frac{(x_n - m_{n-1})^2}{n}.
  !! \f]
  !! Then \f$m_n\f$ and \f$v_n\f$ are the sample mean and sample
  !! variance for \f$\{x_1,\ldots,x_n\}\f$ for each \f$n\f$.
  !!
  !! The derivation of these formulas begins with the definitions for
  !! the running total
  !! \f[
  !!     t_n = \sum_{i=1}^n x_i
  !! \f]
  !! and running sum of square differences
  !! \f[
  !!     s_n = \sum_{i=1}^n (x_i - m_n)^2.
  !! \f]
  !! Then the running mean is \f$m_n = \frac{t_n}{n}\f$, the running
  !! population variance is \f$p_n = \frac{1}{n} \left( s_n - n m_n^2
  !! \right)\f$, and the running sample variance \f$v_n = \frac{1}{n -
  !! 1} \left( s_n - n m_n^2 \right)\f$.
  !!
  !! We can then compute the mean update above, and observing that
  !! \f[
  !!     s_n = \sum_{i=1}^n x_i^2 - n m_n^2
  !! \f]
  !! we can compute the sum-of-square-dfferences update identity
  !! \f[
  !!     s_n = s_{n-1} + \frac{n - 1}{n} (x_n - m_{n-1})^2.
  !! \f]
  !! The algorithm then follows immediately. The population variance
  !! update is given by
  !! \f[
  !!     p_n = \frac{n-1}{n} p_{n-1} + \frac{(x_n - m_n)(x_n - m_{n-1})}{n}
  !!         = \frac{n-1}{n} p_{n-1} + \frac{n-1}{n^2} (x_n - m_{n-1})^2.
  !! \f]
  !!
  !! This algorithm (in a form where \f$m_n\f$ and \f$s_n\f$ are
  !! tracked) originally appeared in:
  !!
  !! B. P. Welford [1962] "Note on a Method for Calculating Corrected
  !! Sums of Squares and Products", Technometrics 4(3), 419-420.
  !!
  !! Numerical tests performed by M. West on 2012-04-12 seem to
  !! indicate that there is no substantial difference between tracking
  !! \f$s_n\f$ versus \f$v_n\f$.
  !!
  !! The same method (tracking \f$m_n\f$ and \f$s_n\f$) is presented
  !! on page 232 in Section 4.2.2 of Knuth:
  !!
  !! D. E. Knuth [1988] "The Art of Computer Programming, Volume 2:
  !! Seminumerical Algorithms", third edition, Addison Wesley Longman,
  !! ISBN 0-201-89684-2.
  !!
  !! An analysis of the error introduced by different variance
  !! computation methods is given in:
  !!
  !! T. F. Chan, G. H. Golub, and R. J. LeVeque [1983] "Algorithms for
  !! Computing the Sample Variance: Analysis and Recommendations", The
  !! American Statistician 37(3), 242-247.
  !!
  !! The relative error in \f$s_n\f$ of Welford's method (tracking
  !! \f$m_n\f$ and \f$s_n\f$) is of order \f$n \kappa \epsilon\f$,
  !! where \f$\epsilon\f$ is the machine precision and \f$\kappa\f$ is
  !! the condition number for the problem, which is given by
  !! \f[
  !!     \kappa = \frac{\|x\|_2}{\sqrt{s_n}}
  !!            = \sqrt{1 + \frac{m_n^2}{p_n}}
  !!            \approx \frac{m_n}{p_n}.
  !! \f]
  !!
  !! This analysis was apparently first given in:
  !!
  !! T. F. C. Chan and J. G. Lewis [1978] "Rounding error analysis of
  !! algorithms for computing means and standard deviations",
  !! Technical Report No. 284, The Johns Hopkins University,
  !! Department of Mathematical Sciences.
  subroutine update_mean_var(mean, var, data, n)

    !> Mean value to update (on entry \f$m_{n-1}\f$, on exit \f$m_n\f$).
    real(kind=dp), intent(inout) :: mean
    !> Variance value to update (on entry \f$v_{n-1}\f$, on exit \f$v_n\f$).
    real(kind=dp), intent(inout) :: var
    !> Data value \f$x_n\f$.
    real(kind=dp), intent(in) :: data
    !> Number \f$n\f$ of this data value.
    integer, intent(in) :: n

    real(kind=dp) :: data_diff

    call assert_msg(376972566, n >= 1, &
         "must have number n >= 1, not: " // trim(integer_to_string(n)))
    if (n == 1) then
       mean = data
       var = 0d0
    else
       data_diff = data - mean
       mean = mean + data_diff / real(n, kind=dp)
       var = real(n - 2, kind=dp) / real(n - 1, kind=dp) * var &
            + data_diff**2 / real(n, kind=dp)
    end if

  end subroutine update_mean_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a fairly tight upper-bound on the Student's t coefficient
  !> for the 95% confidence interval.
  !!
  !! The number of degrees of freedom is one less than \c n_sample. If
  !! a set of \f$n\f$ numbers has sample mean \f$\mu\f$ and sample
  !! standard deviation \f$\sigma\f$, then the 95% confidence interval
  !! for the mean is \f$[\mu - r\sigma/\sqrt{n}, \mu +
  !! r\sigma/\sqrt{n}]\f$, where <tt>r =
  !! student_t_95_coeff(n_sample)</tt>.
  !!
  !! The method used here was written by MW on 2011-05-01, based on
  !! the following empirical observation. If \f$f(\nu) =
  !! t_{0.95,\nu}\f$ is the function we want, where \f$\nu = n - 1\f$
  !! is the number of degrees-of-freedom, then set \f$g(\nu) = f(\nu)
  !! - L\f$, where \f$L = \Phi^{-1}(0.975)\f$ is the limiting value
  !! given by the Gaussian CDF \f$\Phi\f$. We observe numerically that
  !! \f$g'(\nu) < -1\f$ and \f$g'(\nu) \to -1\f$ as \f$\nu \to
  !! \infty\f$. Thus \f$g(\nu)\f$ is well-approximated by \f$A/\nu\f$
  !! for some \f$A\f$. Furthermore, if \f$g(\nu^*) = A^*/\nu\f$, then
  !! \f$g(\nu) < A^*/\nu\f$ for \f$\nu > \nu^*\f$. We thus have
  !! \f$f(\nu) \le (f(\nu^*) - L) (\nu^* / \nu) + L\f$ for \f$\nu \ge
  !! \nu^*\f$. By using a sequence of known \f$(\nu^*, f(\nu^*))\f$
  !! pairs we can thus construct a fairly tight upper bound.
  !!
  !! This implementation has an error of below 0.1% for all values of
  !! \c n_sample.
  real(kind=dp) function student_t_95_coeff(n_sample)

    !> Number of samples.
    integer, intent(in) :: n_sample

    real(kind=dp), parameter :: limit = 1.959963984540054d0
    real(kind=dp), parameter, dimension(15) :: values &
         = (/ 12.7062047364d0, 4.30265272991d0, 3.18244630528d0, &
         2.7764451052d0, 2.57058183661d0, 2.44691184879d0, 2.36462425101d0, &
         2.30600413503d0, 2.26215716274d0, 2.22813885196d0, 2.20098516008d0, &
         2.17881282966d0, 2.16036865646d0, 2.14478668792d0, 2.13144954556d0 /)

    integer :: n_dof

    n_dof = n_sample - 1
    call assert(359779741, n_dof >= 1)
    if (n_dof <= 15) then
       student_t_95_coeff = values(n_dof)
    elseif (n_dof <= 20) then
       student_t_95_coeff = (2.11990529922d0 - limit) * 16d0 &
            / real(n_dof, kind=dp) + limit
    else
       student_t_95_coeff = (2.07961384473d0 - limit) * 21d0 &
            / real(n_dof, kind=dp) + limit
    end if

  end function student_t_95_coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> 95% confidence interval offset from mean.
  !!
  !! If \c mean and \c var are the sample mean and sample variance of
  !! \c n data values, then
  !! <pre>
  !! offset = conf_95_offset(var, n)
  !! </pre>
  !! means that the 95% confidence interval for the mean is
  !! <tt>[mean - offset, mean + offset]</tt>.
  !!
  !! If \c n_sample is one or less then zero is returned.
  function conf_95_offset(var, n_sample)

    !> Sample variance of data.
    real(kind=dp), intent(in) :: var
    !> Number of samples.
    integer, intent(in) :: n_sample

    !> Return offset from mean for the 95% confidence interval.
    real(kind=dp) :: conf_95_offset

    if (n_sample <= 1) then
       conf_95_offset = 0d0
       return
    end if

    conf_95_offset = student_t_95_coeff(n_sample) * sqrt(var) &
         / sqrt(real(n_sample, kind=dp))

  end function conf_95_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write statistics (mean and 95% conf. int.) to a NetCDF file.
  subroutine stats_output_netcdf(stats, ncid, name, unit)

    !> Statistics structure to write.
    type(stats_t), intent(in) :: stats
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit

    call pmc_nc_write_real(ncid, stats%mean, name, unit=unit)
    call pmc_nc_write_real(ncid, stats_conf_95_offset(stats), &
         trim(name) // "_ci_offset", unit=unit)

  end subroutine stats_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write statistics (mean and 95% conf. int.) to a NetCDF file.
  subroutine stats_1d_output_netcdf(stats, ncid, name, dim_name, unit)

    !> Statistics structure to write.
    type(stats_1d_t), intent(in) :: stats
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> NetCDF dimension name for the variable.
    character(len=*), optional, intent(in) :: dim_name
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit

    call pmc_nc_write_real_1d(ncid, stats%mean, name, dim_name=dim_name, &
         unit=unit)
    call pmc_nc_write_real_1d(ncid, stats_1d_conf_95_offset(stats), &
         trim(name) // "_ci_offset", dim_name=dim_name, unit=unit)

  end subroutine stats_1d_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write statistics (mean and 95% conf. int.) to a NetCDF file.
  subroutine stats_2d_output_netcdf(stats, ncid, name, dim_name_1, &
       dim_name_2, unit)

    !> Statistics structure to write.
    type(stats_2d_t), intent(in) :: stats
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Variable name in NetCDF file.
    character(len=*), intent(in) :: name
    !> First NetCDF dimension name for the variable.
    character(len=*), optional, intent(in) :: dim_name_1
    !> Second NetCDF dimension name for the variable.
    character(len=*), optional, intent(in) :: dim_name_2
    !> Unit of variable.
    character(len=*), optional, intent(in) :: unit

    call pmc_nc_write_real_2d(ncid, stats%mean, name, dim_name_1=dim_name_1, &
         dim_name_2=dim_name_2, unit=unit)
    call pmc_nc_write_real_2d(ncid, stats_2d_conf_95_offset(stats), &
         trim(name) // "_ci_offset", dim_name_1=dim_name_1, &
         dim_name_2=dim_name_2, unit=unit)

  end subroutine stats_2d_output_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write statistics (mean and 95% conf. int.) to a text file.
  !!
  !! The format has three columns:
  !!     `dim mean ci_offset`
  !! where \c dim is the optional dimension argument, \c mean is the mean value
  !! and \c ci_offset is the 95% confidence interval offset, so the 95% CI is
  !! `mean - ci_offset, mean + ci_offset`.
  subroutine stats_1d_output_text(stats, filename, dim)

    !> Statistics structure to write.
    type(stats_1d_t), intent(in) :: stats
    !> Filename to write to.
    character(len=*), intent(in) :: filename
    !> Dimension array (independent variable).
    real(kind=dp), intent(in) :: dim(:)

    real(kind=dp), allocatable :: data(:,:)

    if (size(dim) /= size(stats%n)) then
       call die_msg(460147728, 'dim size ' // integer_to_string(size(dim)) &
            // ' does not match stats size ' &
            // integer_to_string(size(stats%n)))
    end if

    allocate(data(size(stats%n), 3))
    data(:,1) = dim
    data(:,2) = stats%mean
    data(:,3) = stats_1d_conf_95_offset(stats)
    call savetxt_2d(filename, data)

  end subroutine stats_1d_output_text

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_stats
