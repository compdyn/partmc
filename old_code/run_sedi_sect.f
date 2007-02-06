! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Copyright (C) Andreas Bott
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Sectional code based on coad1d.f by Andreas Bott
! http://www.meteo.uni-bonn.de/mitarbeiter/ABott/

program run_sedi_sect
  
  use mod_bin
  use mod_array
  use mod_kernel_sedi
  use mod_util
  use mod_init_dist
  use mod_environ
  use mod_material
  
  integer n, scal, isw
  real*8 eps, u, rho, emin, N_0, V_0
  real*8 t_max, t_print, t_progress, del_t
  parameter (n = 220)           ! number of bins
  parameter (scal = 4)          ! bin mesh scale factor
  parameter (isw = 1)           ! kernel (0 = long, 1 = hall, 2 = golovin)
  parameter (eps = 100d0)       ! epsilon (cm^2/s^3)
  parameter (u = 2.5d0)         ! u' (m/s)
  parameter (rho = 1000d0)      ! particle density (kg/m^3)
  parameter (emin = 1d-15)      ! 
  parameter (N_0 = 1d9)         ! particle number concentration (#/m^3)
  parameter (V_0 = 4.1886d-15)  ! mean volume of initial distribution (m^3)
  parameter (t_max = 600d0)     ! total simulation time (seconds)
  parameter (t_print = 60d0)    ! interval between printing (s)
  parameter (t_progress = 10d0) ! interval between progress (s)
  parameter (del_t = 1d0)       ! timestep (s)
  
  real*8 dlnr
  real*8 c(n,n)
  integer ima(n,n)
  real*8 g(n), r(n), e(n)
  real*8 k_bin(n,n), ck(n,n), ec(n,n)
  real*8 taug(n), taup(n), taul(n), tauu(n)
  real*8 prod(n), ploss(n)
  real*8 bin_v(n)
  real*8 n_den(n), bin_g_den(n), bin_n_den(n)
  real*8 time, last_print_time, last_progress_time
  
  integer i, j, i_time, num_t
  logical do_print, do_progress
  
  type(environ) :: env
  type(material) :: mat
  
  real*8 pi
  parameter (pi = 3.14159265358979323846d0)
  
  call allocate_material(mat, 1)
  mat%i_water = 1
  mat%rho = (/ rho /)
  mat%nu = (/ 0 /)
  mat%eps = (/ 0d0 /)
  mat%M_w = (/ 18d-3 /)
  
  env%T = 288d0        ! (K)
  env%RH = 0.999d0     ! (1)
  env%p = 1d5          ! (Pa)
  env%dTdt = -0.01d0   ! (K s^{-1})
  
  ! g   : spectral mass distribution (mg/cm^3)
  ! e   : droplet mass grid (mg)
  ! r   : droplet radius grid (um)
  ! dlnr: constant grid distance of logarithmic grid 
  
  ! mass and radius grid
  call make_bin_grid(n, scal, rho, bin_v, dlnr)
  do i = 1,n
     r(i) = vol2rad(bin_v(i)) * 1d6         ! radius in m to um
     e(i) = bin_v(i) * rho * 1d6   ! volume in m^3 to mass in mg
  enddo
  
  ! initial mass distribution
  call init_exp(V_0, n, bin_v, n_den)
  do i = 1,n
     g(i) = n_den(i) * bin_v(i) * rho * N_0
  enddo
  
  call courant(n, dlnr, scal, c, ima, g, r, e)
  
  ! precompute kernel values for all pairs of bins
  call bin_kernel(n, bin_v, kernel_sedi, env, k_bin)
  call smooth_bin_kernel(n, k_bin, ck)
  do i = 1,n
     do j = 1,n
        ck(i,j) = ck(i,j) * 1d6  ! m^3/s to cm^3/s
     enddo
  enddo
  
  ! multiply kernel with constant timestep and logarithmic grid distance
  do i = 1,n
     do j = 1,n
        ck(i,j) = ck(i,j) * del_t * dlnr
     enddo
  enddo
  
  ! output file
  open(30, file = 'out_sedi_sect.d')
  call print_header(1, n, 1, nint(t_max / t_print) + 1)
  
  ! initialize time
  last_progress_time = 0d0
  time = 0d0
  
  ! initial output
  call check_event(time, del_t, t_print, last_print_time, do_print)
  if (do_print) then
     do i = 1,n
        bin_g_den(i) = g(i) / rho
        bin_n_den(i) = bin_g_den(i) / bin_v(i)
     enddo
     call print_info_density(0d0, n, 1, bin_v, bin_g_den, &
          bin_g_den, bin_n_den, env, mat)
  endif
  
  ! main time-stepping loop
  num_t = nint(t_max / del_t)
  do i_time = 1, num_t
     time = t_max * dble(i_time) / dble(num_t)
     
     call coad(n, del_t, taug, taup, taul, tauu, prod, ploss, &
          c, ima, g, r, e, ck, ec)
     
     ! print output
     call check_event(time, del_t, t_print, last_print_time, &
          do_print)
     if (do_print) then
        do i = 1,n
           bin_g_den(i) = g(i) / rho
           bin_n_den(i) = bin_g_den(i) / bin_v(i)
        enddo
        call print_info_density(0d0, n, 1, bin_v, bin_g_den, &
             bin_g_den, bin_n_den, env, mat)
     endif
     
     ! print progress to stdout
     call check_event(time, del_t, t_progress, last_progress_time, &
          do_progress)
     if (do_progress) then
        write(6,'(a6,a8)') 'step', 'time'
        write(6,'(i6,f8.1)') i_time, time
     endif
  enddo
  
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine coad(n, dt, taug, taup, taul, tauu, prod, ploss, &
       c, ima, g, r, e, ck, ec)
    
    ! collision subroutine, exponential approach
    
    integer n
    real*8 dt
    real*8 taug(n)
    real*8 taup(n)
    real*8 taul(n)
    real*8 tauu(n)
    real*8 prod(n)
    real*8 ploss(n)
    real*8 c(n,n)
    integer ima(n,n)
    real*8 g(n)
    real*8 r(n)
    real*8 e(n)
    real*8 ck(n,n)
    real*8 ec(n,n)
    
    real*8 gmin
    parameter (gmin = 1d-60)
    
    integer i, i0, i1, j, k, kp
    real*8 x0, gsi, gsj, gsk, gk, x1, flux
    
    do i = 1,n
       prod(i) = 0d0
       ploss(i) = 0d0
    enddo
    
    ! lower and upper integration limit i0,i1
    do i0 = 1,(n - 1)
       if (g(i0) .gt. gmin) goto 2000
    enddo
2000 continue
    do i1 = (n - 1),1,-1
       if (g(i1) .gt. gmin) goto 2010
    enddo
2010 continue
    
    do i = i0,i1
       do j = i,i1
          k = ima(i,j)
          kp = k + 1
          
          x0 = ck(i,j) * g(i) * g(j)
          x0 = min(x0, g(i) * e(j))
          
          if (j .ne. k) x0 = min(x0, g(j) * e(i))
          gsi = x0 / e(j)
          gsj = x0 / e(i)
          gsk = gsi + gsj
          
          ! loss from positions i, j
          ploss(i) = ploss(i) + gsi
          ploss(j) = ploss(j) + gsj
          g(i) = g(i) - gsi
          g(j) = g(j) - gsj
          gk = g(k) + gsk
          
          if (gk .gt. gmin) then
             x1 = log(g(kp) / gk + 1d-60)
             flux = gsk / x1 * (exp(0.5d0 * x1) &
                  - exp(x1 * (0.5d0 - c(i,j))))
             flux = min(flux, gk)
             g(k) = gk - flux
             g(kp) = g(kp) + flux
             ! gain for positions i, j
             prod(k) =  prod(k) + gsk - flux           
             prod(kp) = prod(kp) + flux
          endif
       enddo
    enddo
    
  end subroutine coad
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine courant(n, dlnr, scal, c, ima, g, r, e)
    
    integer n
    real*8 dlnr
    integer scal
    real*8 c(n,n)
    integer ima(n,n)
    real*8 g(n)
    real*8 r(n)
    real*8 e(n)
    
    integer i, j, k, kk
    real*8 x0
    
    do i = 1,n
       do j = i,n
          x0 = e(i) + e(j)
          do k = j,n
             if ((e(k) .ge. x0) .and. (e(k-1) .lt. x0)) then
                if (c(i,j) .lt. 1d0 - 1d-08) then
                   kk = k - 1
                   c(i,j) = log(x0 / e(k-1)) / (3d0 * dlnr)
                else
                   c(i,j) = 0d0
                   kk = k
                endif
                ima(i,j) = min(n - 1, kk)
                goto 2000
             endif
          enddo
2000      continue
          c(j,i) = c(i,j)
          ima(j,i) = ima(i,j)
       enddo
    enddo
    
  end subroutine courant
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine smooth_bin_kernel(n_bin, k, k_smooth)
    
    ! Smooths kernel values for bin pairs, and halves the self-rate.
    
    integer, intent(in) :: n_bin                !  number of bins
    real*8, intent(in) :: k(n_bin,n_bin)        !  kernel values
    real*8, intent(out) :: k_smooth(n_bin,n_bin) !  smoothed kernel values
    
    integer i, j, im, ip, jm, jp
    
    do i = 1,n_bin
       do j = 1,n_bin
          jm = max0(j - 1, 1)
          im = max0(i - 1, 1)
          jp = min0(j + 1, n_bin)
          ip = min0(i + 1, n_bin)
          k_smooth(i,j) = 0.125d0 * (k(i,jm) + k(im,j) &
               + k(ip,j) + k(i,jp)) &
               + 0.5d0 * k(i,j)
          if (i .eq. j) then
             k_smooth(i,j) = 0.5d0 * k_smooth(i,j)
          endif
       enddo
    enddo
    
  end subroutine smooth_bin_kernel
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program run_sedi_sect
