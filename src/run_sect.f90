! Copyright (C) 2005-2010 Nicole Riemer and Matthew West
! Copyright (C) Andreas Bott
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_run_sect module.

!> 1D sectional simulation.
!!
!! Sectional code based on \c coad1d.f by Andreas Bott
!!     - http://www.meteo.uni-bonn.de/mitarbeiter/ABott/
!!     - Released under the GPL to Nicole Riemer (personal communication)
!!     - A. Bott, A flux method for the numerical solution of the
!!       stochastic collection equation, J. Atmos. Sci. 55, 2284-2293,
!!       1998.
module pmc_run_sect

  use pmc_bin_grid
  use pmc_aero_binned
  use pmc_util
  use pmc_aero_dist
  use pmc_env_data
  use pmc_env_state
  use pmc_aero_data
  use pmc_kernel
  use pmc_output
  use pmc_gas_data
  use pmc_gas_state

  !> Options to control the operation of run_sect().
  type run_sect_opt_t
     !> Final time (s).
    real(kind=dp) :: t_max
    !> Timestep for coagulation (s).
    real(kind=dp) :: del_t
    !> Output interval (0 disables) (s).
    real(kind=dp) :: t_output
    !> Progress interval (0 disables) (s).
    real(kind=dp) :: t_progress
    !> Whether to do coagulation.
    logical :: do_coagulation
    !> Output prefix.
     character(len=300) :: prefix
    !> Type of coagulation kernel.
    integer :: kernel_type
     !> UUID of the simulation.
     character(len=PMC_UUID_LEN) :: uuid
  end type run_sect_opt_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run a sectional simulation.
  subroutine run_sect(bin_grid, gas_data, aero_data, aero_dist, &
       env_data, env_state, sect_opt)
  
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol distribution.
    type(aero_dist_t), intent(inout) :: aero_dist
    !> Environment data.
    type(env_data_t), intent(inout) :: env_data
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Options.
    type(run_sect_opt_t), intent(in) :: sect_opt
    
    real(kind=dp) c(bin_grid%n_bin,bin_grid%n_bin)
    integer ima(bin_grid%n_bin,bin_grid%n_bin)
    real(kind=dp) g(bin_grid%n_bin), r(bin_grid%n_bin), e(bin_grid%n_bin)
    real(kind=dp) k_bin(bin_grid%n_bin,bin_grid%n_bin)
    real(kind=dp) ck(bin_grid%n_bin,bin_grid%n_bin)
    real(kind=dp) ec(bin_grid%n_bin,bin_grid%n_bin)
    real(kind=dp) taug(bin_grid%n_bin), taup(bin_grid%n_bin)
    real(kind=dp) taul(bin_grid%n_bin), tauu(bin_grid%n_bin)
    real(kind=dp) prod(bin_grid%n_bin), ploss(bin_grid%n_bin)
    real(kind=dp) time, last_output_time, last_progress_time
    type(env_state_t) :: old_env_state
    type(aero_binned_t) :: aero_binned
    type(gas_state_t) :: gas_state
    
    integer i, j, i_time, num_t, i_summary
    logical do_output, do_progress
  
    ! g         : spectral mass distribution (mg/cm^3)
    ! e         : droplet mass grid (mg)
    ! r         : droplet radius grid (um)
    ! log_width : constant grid distance of logarithmic grid 

    if (aero_data%n_spec /= 1) then
       call die_msg(844211192, &
            'run_sect() can only use one aerosol species')
    end if

    ! output data structure
    call aero_binned_allocate_size(aero_binned, bin_grid%n_bin, &
         aero_data%n_spec)
    aero_binned%vol_conc = 0d0
    call gas_state_allocate_size(gas_state, gas_data%n_spec)
    
    ! mass and radius grid
    do i = 1,bin_grid%n_bin
       r(i) = bin_grid%center_radius(i) * 1d6 ! radius in m to um
       e(i) = rad2vol(bin_grid%center_radius(i)) &
            * aero_data%density(1) * 1d6 ! vol in m^3 to mass in mg
    end do
    
    ! initial mass distribution
    call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
         aero_dist)
    
    call courant(bin_grid%n_bin, bin_grid%log_width, e, ima, c)
    
    ! initialize time
    last_progress_time = 0d0
    time = 0d0
    i_summary = 1
    
    ! precompute kernel values for all pairs of bins
    call bin_kernel(bin_grid%n_bin, bin_grid%center_radius, aero_data, &
         sect_opt%kernel_type, env_state, k_bin)
    call smooth_bin_kernel(bin_grid%n_bin, k_bin, ck)
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          ck(i,j) = ck(i,j) * 1d6  ! m^3/s to cm^3/s
       end do
    end do
    
    ! multiply kernel with constant timestep and logarithmic grid distance
    do i = 1,bin_grid%n_bin
       do j = 1,bin_grid%n_bin
          ck(i,j) = ck(i,j) * sect_opt%del_t * bin_grid%log_width
       end do
    end do
    
    ! initial output
    call check_event(time, sect_opt%del_t, sect_opt%t_output, &
         last_output_time, do_output)
    if (do_output) then
       call output_sectional(sect_opt%prefix, bin_grid, aero_data, &
            aero_binned, gas_data, gas_state, env_state, i_summary, &
            time, sect_opt%t_output, sect_opt%uuid)
    end if
    
    ! main time-stepping loop
    num_t = nint(sect_opt%t_max / sect_opt%del_t)
    call env_state_allocate(old_env_state)
    do i_time = 1, num_t

       if (sect_opt%do_coagulation) then
          g = aero_binned%vol_conc(:,1) * aero_data%density(1)
          call coad(bin_grid%n_bin, sect_opt%del_t, taug, taup, taul, &
               tauu, prod, ploss, c, ima, g, r, e, ck, ec)
          aero_binned%vol_conc(:,1) = g / aero_data%density(1)
          aero_binned%num_conc = aero_binned%vol_conc(:,1) &
               / rad2vol(bin_grid%center_radius)
       end if

       time = sect_opt%t_max * real(i_time, kind=dp) / real(num_t, kind=dp)

       call env_state_copy(env_state, old_env_state)
       call env_data_update_state(env_data, env_state, time, &
            update_rel_humid = .true.)
       call env_state_update_gas_state(env_state, sect_opt%del_t, &
            old_env_state, gas_data, gas_state)
       call env_state_update_aero_binned(env_state, sect_opt%del_t, &
            old_env_state, bin_grid, aero_data, aero_binned)
       
       ! print output
       call check_event(time, sect_opt%del_t, sect_opt%t_output, &
            last_output_time, do_output)
       if (do_output) then
          i_summary = i_summary + 1
          call output_sectional(sect_opt%prefix, bin_grid, aero_data, &
               aero_binned, gas_data, gas_state, env_state, i_summary, &
               time, sect_opt%t_output, sect_opt%uuid)
       end if
       
       ! print progress to stdout
       call check_event(time, sect_opt%del_t, sect_opt%t_progress, &
            last_progress_time, do_progress)
       if (do_progress) then
          write(*,'(a6,a8)') 'step', 'time'
          write(*,'(i6,f8.1)') i_time, time
       end if
    end do

    call aero_binned_deallocate(aero_binned)
    call gas_state_deallocate(gas_state)

  end subroutine run_sect
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Collision subroutine, exponential approach.
  subroutine coad(n_bin, dt, taug, taup, taul, tauu, prod, ploss, &
       c, ima, g, r, e, ck, ec)
    
    integer n_bin
    real(kind=dp) dt
    real(kind=dp) taug(n_bin)
    real(kind=dp) taup(n_bin)
    real(kind=dp) taul(n_bin)
    real(kind=dp) tauu(n_bin)
    real(kind=dp) prod(n_bin)
    real(kind=dp) ploss(n_bin)
    real(kind=dp) c(n_bin,n_bin)
    integer ima(n_bin,n_bin)
    real(kind=dp) g(n_bin)
    real(kind=dp) r(n_bin)
    real(kind=dp) e(n_bin)
    real(kind=dp) ck(n_bin,n_bin)
    real(kind=dp) ec(n_bin,n_bin)
    
    real(kind=dp), parameter :: gmin = 1d-60
    
    integer i, i0, i1, j, k, kp
    real(kind=dp) x0, gsi, gsj, gsk, gk, x1, flux

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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the Courant number for each bin pair.
  subroutine courant(n_bin, log_width, e, ima, c)

    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Bin scale factor.
    real(kind=dp), intent(in) :: log_width
    !> Droplet mass grid (mg).
    real(kind=dp), intent(in) :: e(n_bin)
    !> i + j goes in bin ima(i,j).
    integer, intent(out) :: ima(n_bin,n_bin)
    !> Courant number for bin pairs.
    real(kind=dp), intent(out) :: c(n_bin,n_bin)
    
    integer i, j, k, kk
    real(kind=dp) x0

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
                c(i,j) = log(x0 / e(k-1)) / (3d0 * log_width)
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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Smooths kernel values for bin pairs, and halves the self-rate.
  subroutine smooth_bin_kernel(n_bin, k, k_smooth)
    
    !> Number of bins.
    integer, intent(in) :: n_bin
    !> Kernel values.
    real(kind=dp), intent(in) :: k(n_bin,n_bin)
    !> Smoothed kernel values.
    real(kind=dp), intent(out) :: k_smooth(n_bin,n_bin)
    
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
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_run_sect
