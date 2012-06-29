! Copyright (C) 2012 Jian Tian
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The test_radii_conversion program.

!> Do radii conversions between geometric and mobility
!> equivalent radius and volume under various fractal parameters
!> to check if the results are within tolerance.
program test_radii_conversion
  
  use pmc_util
  use pmc_constants
  use pmc_fractal

  type(fractal_t) :: fractal  
  real(kind=dp) :: d_f(7), f(2), R0(4), Rme(100), Rgeo(100), vol(100)
  integer :: i, j, k, m, vcount_Rme, vcount_Rgeo
  real(kind=dp), allocatable :: Rme_err(:), Rgeo_err(:)
  real(kind=dp), parameter :: tk = 300d0
  real(kind=dp), parameter :: press = 1d+5

  call fractal_allocate(fractal)

  d_f = (/1.8d0, 2d0, 2.2d0, 2.4d0, 2.6d0, 2.8d0, 3d0/)
  f = (/1d0, 1.43d0/)
  R0 = (/1d-10, 1d-9, 1d-8, 1d-7/)
  
  vcount_Rme = 0
  vcount_Rgeo = 0
  do i = 1, size(R0)
     fractal%prime_radius = R0(i)
     call logspace(fractal%prime_radius, fractal%prime_radius*1d+3, Rme)
     call logspace(fractal%prime_radius, fractal%prime_radius*1d+3, Rgeo)
     do j = 1, size(f)
        fractal%vol_fill_factor = f(j)
        do k = 1, size(d_f)
           fractal%frac_dim = d_f(k)
           do m = 1, size(Rme)
              if (.not. almost_equal(Rme(m), vol2Rme(Rme2vol(Rme(m), tk, press, fractal) &
                   , tk, press, fractal))) then
                 vcount_Rme = vcount_Rme + 1
                 write(*, '(a, e12.3, 3x, a, f5.2, 3x, a, f5.2, 3x, a, e12.3)') &
                      'Mobility equivalent radius values mismatch at: Rme = ', &
                      Rme(m), 'df = ', d_f(k), 'f = ', f(j), 'R0 = ', R0(i)
              end if
              if (.not. almost_equal(Rgeo(m), vol2rad(rad2vol(Rgeo(m), fractal), fractal))) then
                 vcount_Rgeo = vcount_Rgeo + 1
                 write(*, '(a, e12.3, 3x, a, f5.2, 3x, a, f5.2, 3x, a, e12.3)') &
                      'Geometric radius values mismatch at: Rgeo = ', &
                      Rgeo(m), 'df = ', d_f(k), 'f = ', f(j), 'R0 = ', R0(i)
              end if
           end do
        end do
     end do
  end do

  allocate (Rme_err(vcount_Rme))
  allocate (Rgeo_err(vcount_Rgeo))

  if (size(Rme_err) == 0) then
     write(*,'(a)') 'Mobility equivalent radius values match within the given relative tolerance.'
  end if

  if (size(Rgeo_err) == 0) then
     write(*,'(a)') 'Geometric radius values match within the given relative tolerance.'
  end if

  call fractal_deallocate(fractal)

end program test_radii_conversion
