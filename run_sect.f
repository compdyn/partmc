! -*- mode: f90; -*-
! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Copyright (C) Andreas Bott
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Sectional code based on coad1d.f by Andreas Bott
! http://www.meteo.uni-bonn.de/mitarbeiter/ABott/
! A. Bott, A flux method for the numerical solution of the stochastic
! collection equation, J. Atmos. Sci. 55, 2284-2293, 1998.

module mod_run_sect
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run_sect(n_bin, bin_v, dlnr, n_den, num_conc, kernel, &
       t_max, del_t, t_output, t_progress, output_unit, mat, env)

    ! Run a sectional simulation.
  
    use mod_bin
    use mod_array
    use mod_kernel_sedi
    use mod_util
    use mod_init_dist
    use mod_environ
    use mod_material

    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(in) :: dlnr          ! bin scale factor
    real*8, intent(in) :: n_den(n_bin)  ! initial number density
    real*8, intent(in) :: num_conc      ! particle concentration (#/m^3)
    
    real*8, intent(in) :: t_max         ! final time (seconds)
    real*8, intent(in) :: del_t         ! timestep for coagulation
    real*8, intent(in) :: t_output      ! output interval (0 disables) (s)
    real*8, intent(in) :: t_progress    ! progress interval (0 disables) (s)
    integer, intent(in) :: output_unit  ! unit number to output to

    type(environ), intent(inout) :: env ! environment state
    type(material), intent(in) :: mat   ! material properties
    
    real*8 c(n_bin,n_bin)
    integer ima(n_bin,n_bin)
    real*8 g(n_bin), r(n_bin), e(n_bin)
    real*8 k_bin(n_bin,n_bin), ck(n_bin,n_bin), ec(n_bin,n_bin)
    real*8 taug(n_bin), taup(n_bin), taul(n_bin), tauu(n_bin)
    real*8 prod(n_bin), ploss(n_bin)
    real*8 bin_g_den(n_bin), bin_n_den(n_bin)
    real*8 time, last_output_time, last_progress_time
    
    integer i, j, i_time, num_t
    logical do_output, do_progress
  
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_environ
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(environ), intent(in) :: env 
         real*8, intent(out) :: k
       end subroutine kernel
    end interface

    ! g   : spectral mass distribution (mg/cm^3)
    ! e   : droplet mass grid (mg)
    ! r   : droplet radius grid (um)
    ! dlnr: constant grid distance of logarithmic grid 
    
    ! mass and radius grid
    do i = 1,n_bin
       r(i) = vol2rad(bin_v(i)) * 1d6       ! radius in m to um
       e(i) = bin_v(i) * mat%rho(1) * 1d6   ! volume in m^3 to mass in mg
    end do
    
    ! initial mass distribution
    do i = 1,n_bin
       g(i) = n_den(i) * bin_v(i) * mat%rho(1) * num_conc
       if (g(i) .le. 1d-80) g(i) = 0d0    ! avoid problem with gnuplot
    end do
    
    call courant(n_bin, dlnr, e, ima, c)
    
    ! precompute kernel values for all pairs of bins
    call bin_kernel(n_bin, bin_v, kernel_sedi, env, k_bin)
    call smooth_bin_kernel(n_bin, k_bin, ck)
    do i = 1,n_bin
       do j = 1,n_bin
          ck(i,j) = ck(i,j) * 1d6  ! m^3/s to cm^3/s
       end do
    end do
    
    ! multiply kernel with constant timestep and logarithmic grid distance
    do i = 1,n_bin
       do j = 1,n_bin
          ck(i,j) = ck(i,j) * del_t * dlnr
       end do
    end do
    
    ! initialize time
    last_progress_time = 0d0
    time = 0d0
    call init_environ(env, time)
    
    ! initial output
    call check_event(time, del_t, t_output, last_output_time, do_output)
    if (do_output) then
       do i = 1,n_bin
          bin_g_den(i) = g(i) / mat%rho(1)
          bin_n_den(i) = bin_g_den(i) / bin_v(i)
       end do
       call output_info_density(output_unit, 0d0, n_bin, 1, bin_v, bin_g_den, &
            bin_g_den, bin_n_den, env, mat, 1)
    end if
    
    ! main time-stepping loop
    num_t = nint(t_max / del_t)
    do i_time = 1, num_t
       time = t_max * dble(i_time) / dble(num_t)
       call update_environ(env, time)
       
       call coad(n_bin, del_t, taug, taup, taul, tauu, prod, ploss, &
            c, ima, g, r, e, ck, ec)
       
       ! print output
       call check_event(time, del_t, t_output, last_output_time, &
            do_output)
       if (do_output) then
          do i = 1,n_bin
             bin_g_den(i) = g(i) / mat%rho(1)
             bin_n_den(i) = bin_g_den(i) / bin_v(i)
          end do
          call output_info_density(output_unit, 0d0, n_bin, 1, bin_v, &
               bin_g_den, bin_g_den, bin_n_den, env, mat, 1)
       end if
       
       ! print progress to stdout
       call check_event(time, del_t, t_progress, last_progress_time, &
            do_progress)
       if (do_progress) then
          write(6,'(a6,a8)') 'step', 'time'
          write(6,'(i6,f8.1)') i_time, time
       end if
    end do

  end subroutine run_sect
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine coad(n_bin, dt, taug, taup, taul, tauu, prod, ploss, &
       c, ima, g, r, e, ck, ec)
    
    ! Collision subroutine, exponential approach.
    
    integer n_bin
    real*8 dt
    real*8 taug(n_bin)
    real*8 taup(n_bin)
    real*8 taul(n_bin)
    real*8 tauu(n_bin)
    real*8 prod(n_bin)
    real*8 ploss(n_bin)
    real*8 c(n_bin,n_bin)
    integer ima(n_bin,n_bin)
    real*8 g(n_bin)
    real*8 r(n_bin)
    real*8 e(n_bin)
    real*8 ck(n_bin,n_bin)
    real*8 ec(n_bin,n_bin)
    
    real*8, parameter :: gmin = 1d-60
    
    integer i, i0, i1, j, k, kp
    real*8 x0, gsi, gsj, gsk, gk, x1, flux

    do i = 1,n_bin
       prod(i) = 0d0
       ploss(i) = 0d0
    end do
    
    ! lower and upper integration limit i0,i1
    do i0 = 1,(n_bin - 1)
       if (g(i0) .gt. gmin) goto 2000
    end do
2000 continue
    do i1 = (n_bin - 1),1,-1
       if (g(i1) .gt. gmin) goto 2010
    end do
2010 continue
    
    do i = i0,i1
       do j = i,i1
          k = ima(i,j) ! k = 0 means that i + j goes nowhere
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

          if (k > 0) then ! do we have a valid bin for the coagulation result?
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
             end if
          end if
       end do
    end do
    
  end subroutine coad
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine courant(n_bin, dlnr, e, ima, c)

    ! Determines the Courant number for each bin pair.
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    real*8, intent(in) :: e(n_bin)      ! droplet mass grid (mg)
    integer, intent(out) :: ima(n_bin,n_bin) ! i + j goes in bin ima(i,j)
    real*8, intent(out) :: c(n_bin,n_bin) ! Courant number for bin pairs
    
    integer i, j, k, kk
    real*8 x0

    c = 0d0 ! added to avoid uninitialized access errors
    ima = 0 ! ima(i,j) = 0 means that particles i + j go nowhere
    do i = 1,n_bin
       do j = i,n_bin
          x0 = e(i) + e(j)
          ! FIXME: replace below with particle_in_bin()
          do k = j,n_bin
             if ((e(k) .ge. x0) .and. (e(k-1) .lt. x0)) then
                if (c(i,j) .lt. 1d0 - 1d-08) then
                   kk = k - 1
                   c(i,j) = log(x0 / e(k-1)) / (3d0 * dlnr)
                else
                   c(i,j) = 0d0
                   kk = k
                end if
                ima(i,j) = min(n_bin - 1, kk)
                goto 2000
             end if
          end do
2000      continue
          c(j,i) = c(i,j)
          ima(j,i) = ima(i,j)
       end do
    end do
    
  end subroutine courant
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine smooth_bin_kernel(n_bin, k, k_smooth)
    
    ! Smooths kernel values for bin pairs, and halves the self-rate.
    
    integer, intent(in) :: n_bin        ! number of bins
    real*8, intent(in) :: k(n_bin,n_bin) ! kernel values
    real*8, intent(out) :: k_smooth(n_bin,n_bin) ! smoothed kernel values
    
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
          end if
       end do
    end do
    
  end subroutine smooth_bin_kernel
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module mod_run_sect
