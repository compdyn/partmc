! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Copyright (C) Andreas Bott
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Sectional code based on coad1d.f by Andreas Bott
! http://www.meteo.uni-bonn.de/mitarbeiter/ABott/
! Released under the GPL to Nicole Riemer (personal communication)
! A. Bott, A flux method for the numerical solution of the stochastic
! collection equation, J. Atmos. Sci. 55, 2284-2293, 1998.

module mod_run_sect

  use mod_inout

  type run_sect_opt_t
    real*8 :: t_max                     ! final time (s)
    real*8 :: del_t                     ! timestep for coagulation (s)
    real*8 :: t_output                  ! output interval (0 disables) (s)
    real*8 :: t_progress                ! progress interval (0 disables) (s)
    logical :: do_coagulation           ! whether to do coagulation
  end type run_sect_opt_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run_sect(bin_grid, gas_data, aero_data, aero_dist, env, &
       kernel, sect_opt, summary_file)

    ! Run a sectional simulation.
  
    use mod_bin_grid
    use mod_aero_binned
    use mod_kernel_sedi
    use mod_util
    use mod_aero_dist
    use mod_env
    use mod_aero_data
    use mod_kernel
    use mod_output_summary
    use mod_gas_data
    use mod_gas_state

    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_dist_t), intent(inout) :: aero_dist ! aerosol distribution
    type(env_t), intent(inout) :: env   ! environment state
    type(run_sect_opt_t), intent(in) :: sect_opt ! options
    type(inout_file_t), intent(inout) :: summary_file ! summary output file
    
    real*8 c(bin_grid%n_bin,bin_grid%n_bin)
    integer ima(bin_grid%n_bin,bin_grid%n_bin)
    real*8 g(bin_grid%n_bin), r(bin_grid%n_bin), e(bin_grid%n_bin)
    real*8 k_bin(bin_grid%n_bin,bin_grid%n_bin)
    real*8 ck(bin_grid%n_bin,bin_grid%n_bin), ec(bin_grid%n_bin,bin_grid%n_bin)
    real*8 taug(bin_grid%n_bin), taup(bin_grid%n_bin)
    real*8 taul(bin_grid%n_bin), tauu(bin_grid%n_bin)
    real*8 prod(bin_grid%n_bin), ploss(bin_grid%n_bin)
    real*8 time, last_output_time, last_progress_time
    type(aero_binned_t) :: aero_binned
    type(gas_state_t) :: gas_state
    
    integer i, j, i_time, num_t
    logical do_output, do_progress
  
    interface
       subroutine kernel(v1, v2, env, k)
         use mod_env
         real*8, intent(in) :: v1
         real*8, intent(in) :: v2
         type(env_t), intent(in) :: env   
         real*8, intent(out) :: k
       end subroutine kernel
    end interface

    ! g   : spectral mass distribution (mg/cm^3)
    ! e   : droplet mass grid (mg)
    ! r   : droplet radius grid (um)
    ! dlnr: constant grid distance of logarithmic grid 

    if (aero_data%n_spec /= 1) then
       write(0,*) 'ERROR: run_sect() can only handle one aerosol species at present'
       call exit(1)
    end if

    ! output data structure
    call aero_binned_alloc(bin_grid%n_bin, aero_data%n_spec, aero_binned)
    aero_binned%vol_den = 0d0
    call gas_state_alloc(gas_data%n_spec, gas_state)
    
    ! mass and radius grid
    do i = 1,bin_grid%n_bin
       r(i) = vol2rad(bin_grid%v(i)) * 1d6           ! radius in m to um
       e(i) = bin_grid%v(i) * aero_data%rho(1) * 1d6 ! vol in m^3 to mass in mg
    end do
    
    ! initial mass distribution
    call aero_dist_add_to_binned(bin_grid, aero_dist, aero_binned)
    ! avoid problem with gnuplot
    where (aero_binned%num_den .le. 1d-80) aero_binned%num_den = 0d0
    where (aero_binned%vol_den .le. 1d-80) aero_binned%vol_den = 0d0
    
    call courant(bin_grid%n_bin, bin_grid%dlnr, e, ima, c)
    
    ! precompute kernel values for all pairs of bins
    call bin_kernel(bin_grid%n_bin, bin_grid%v, kernel_sedi, env, k_bin)
    call smooth_bin_kernel(bin_grid%n_bin, k_bin, ck)
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          ck(i,j) = ck(i,j) * 1d6  ! m^3/s to cm^3/s
       end do
    end do
    
    ! multiply kernel with constant timestep and logarithmic grid distance
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          ck(i,j) = ck(i,j) * sect_opt%del_t * bin_grid%dlnr
       end do
    end do
    
    ! initialize time
    last_progress_time = 0d0
    time = 0d0
    call init_environ(env, time)
    
    ! initial output
    call check_event(time, sect_opt%del_t, sect_opt%t_output, &
         last_output_time, do_output)
    if (do_output) then
       call output_summary(summary_file, 0d0, &
            bin_grid, aero_data, aero_binned, gas_data, gas_state, env, 1)
    end if
    
    ! main time-stepping loop
    num_t = nint(sect_opt%t_max / sect_opt%del_t)
    do i_time = 1, num_t

       if (sect_opt%do_coagulation) then
          g = aero_binned%vol_den(:,1) * aero_data%rho(1)
          call coad(bin_grid%n_bin, sect_opt%del_t, taug, taup, taul, &
               tauu, prod, ploss, c, ima, g, r, e, ck, ec)
          aero_binned%vol_den(:,1) = g / aero_data%rho(1)
          aero_binned%num_den = aero_binned%vol_den(:,1) / bin_grid%v
          ! avoid problem with gnuplot
          where (aero_binned%num_den .le. 1d-80) aero_binned%num_den = 0d0
          where (aero_binned%vol_den .le. 1d-80) aero_binned%vol_den = 0d0
       end if

       time = sect_opt%t_max * dble(i_time) / dble(num_t)

       call update_environ(env, time)
       call environ_update_gas_state(env, sect_opt%del_t, gas_data, gas_state)
       call environ_update_aero_binned(env, sect_opt%del_t, bin_grid, &
            aero_data, aero_binned)
       
       ! print output
       call check_event(time, sect_opt%del_t, sect_opt%t_output, &
            last_output_time, do_output)
       if (do_output) then
          call output_summary(summary_file, time, &
               bin_grid, aero_data, aero_binned, gas_data, gas_state, env, 1)
       end if
       
       ! print progress to stdout
       call check_event(time, sect_opt%del_t, sect_opt%t_progress, &
            last_progress_time, do_progress)
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

    use mod_util
    
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
          ! FIXME: should use particle_in_bin(), but that is actually
          ! slightly different than what was always done here
          k = find_1d(n_bin, e, x0)
          if (k < n_bin) then
             k = k + 1
             if (c(i,j) .lt. 1d0 - 1d-08) then
                kk = k - 1
                c(i,j) = log(x0 / e(k-1)) / (3d0 * dlnr)
             else
                c(i,j) = 0d0
                kk = k
             end if
             ima(i,j) = min(n_bin - 1, kk)
          end if
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
