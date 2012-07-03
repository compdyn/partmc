! Copyright (C) 2009-2012 Matthew West
! Copyright (C) 2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The test_fractal_radii_conversion program.

!> Do radii conversions between geometric and mobility
!> equivalent radius and volume under various fractal parameters
!> to check if the results are within tolerance.
program test_fractal_radii_conversion

  use pmc_util
  use pmc_constants
  use pmc_fractal

  real(kind=dp), parameter :: temp = 300d0
  real(kind=dp), parameter :: pressure = 1d5

  type(fractal_t) :: fractal
  real(kind=dp) :: d_f(7), f(2), R0(4), Rme(100), Rgeo(100)
  integer :: i_df, i_f, i_R0, i_part
  logical :: Rme_match, Rgeo_match

  call fractal_set_spherical(fractal)

  d_f = (/1.8d0, 2d0, 2.2d0, 2.4d0, 2.6d0, 2.8d0, 3d0/)
  f = (/1d0, 1.43d0/)
  R0 = (/1d-10, 1d-9, 1d-8, 1d-7/)

  Rme_match = .true.
  Rgeo_match = .true.
  do i_R0 = 1,size(R0)
     fractal%prime_radius = R0(i_R0)
     Rme = logspace(fractal%prime_radius, fractal%prime_radius * 1d3, 100)
     Rgeo = logspace(fractal%prime_radius, fractal%prime_radius * 1d3, 100)
     do i_f = 1,size(f)
        fractal%vol_fill_factor = f(i_f)
        do i_df = 1,size(d_f)
           fractal%frac_dim = d_f(i_df)
           do i_part = 1,size(Rme)
              if (.not. almost_equal(Rme(i_part), &
                   fractal_vol_to_mobility_rad(fractal, &
                   fractal_mobility_rad_to_vol(fractal, Rme(i_part), &
                   temp, pressure), temp, pressure))) then
                 Rme_match = .false.
                 write(*, '(a,e12.3,3x,a,f5.2,3x,a,f5.2,3x,a,e12.3)') &
                      'Mobility equivalent radii mismatch at: Rme = ', &
                      Rme(i_part), 'df = ', d_f(i_df), 'f = ', f(i_f), &
                      'R0 = ', R0(i_R0)
              end if
              if (.not. almost_equal(Rgeo(i_part), &
                   fractal_vol2rad(fractal, fractal_rad2vol(fractal, Rgeo(i_part))))) then
                 Rgeo_match = .false.
                 write(*, '(a,e12.3,3x,a,f5.2,3x,a,f5.2,3x,a,e12.3)') &
                      'Geometric radii mismatch at: Rgeo = ', &
                      Rgeo(i_part), 'df = ', d_f(i_df), 'f = ', f(i_f), &
                      'R0 = ', R0(i_R0)
              end if
           end do
        end do
     end do
  end do

  call assert_msg(435301873, Rme_match, &
       "Test failed: Mobility equivalent radius does not match.")
  write(*,'(a)') 'Mobility equivalent radius values match ' &
       // 'within the given relative tolerance.'

  call assert_msg(812341709, Rgeo_match, &
       "Test failed: Geometric radius does not match.")
  write(*,'(a)') 'Geometric radius values match ' &
       // 'within the given relative tolerance.'

end program test_fractal_radii_conversion
