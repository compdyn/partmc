! Copyright (C) 2005-2012 Nicole Riemer and Matthew West
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
  use pmc_scenario
  use pmc_env_state
  use pmc_aero_data
  use pmc_coag_kernel
  use pmc_output
  use pmc_gas_data
  use pmc_gas_state
  use pmc_tchem_interface
  use pmc_mosaic

  !> Options controlling the operation of run_sect().
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
    !> Whether to do TChem chemistry.
    logical :: do_tchem
    !> Whether to do MOSAIC chemistry.
    logical :: do_mosaic
    !> Output prefix.
     character(len=300) :: prefix
    !> Type of coagulation kernel.
    integer :: coag_kernel_type
     !> UUID of the simulation.
     character(len=PMC_UUID_LEN) :: uuid
    !> NetCDF file to read the initial aero_binned state from.
    character(len=PMC_MAX_FILENAME_LEN) :: aero_init_filename
  end type run_sect_opt_t

  !> Bott (1998) single-moment flux coagulation (the historical scheme):
  !! advects the volume distribution and derives number from it.
  integer, parameter :: SECT_COAG_BOTT = 1
  !> Jacobson (2002) semi-implicit two-moment coagulation: evolves the
  !! per-bin number and per-species volume as independent moments. This is
  !! the single-type algorithm used by MOSAIC (jacobson2002_singletype_coag
  !! in module_mosaic_coag1d).
  integer, parameter :: SECT_COAG_TWO_MOMENT = 2

  !> Coagulation scheme used by run_sect(), selected at compile time.
  !! Set to SECT_COAG_TWO_MOMENT to use the two-moment method instead of the
  !! historical Bott flux scheme.
  integer, parameter :: SECT_COAG_METHOD = 2 !SECT_COAG_BOTT

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Run a sectional simulation.
  subroutine run_sect(bin_grid, gas_data, gas_state, aero_data, aero_dist, &
       scenario, env_state, run_sect_opt)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state (initial value on entry; updated during the run).
    type(gas_state_t), intent(inout) :: gas_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol distribution.
    type(aero_dist_t), intent(inout) :: aero_dist
    !> Environment data.
    type(scenario_t), intent(inout) :: scenario
    !> Environment state.
    type(env_state_t), intent(inout) :: env_state
    !> Options.
    type(run_sect_opt_t), intent(in) :: run_sect_opt

    real(kind=dp) c(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    integer ima(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) g(bin_grid_size(bin_grid)), r(bin_grid_size(bin_grid))
    real(kind=dp) gs(bin_grid_size(bin_grid),aero_data_n_spec(aero_data))
    real(kind=dp) bin_vol_frac(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))
    real(kind=dp) e(bin_grid_size(bin_grid))
    real(kind=dp) k_bin(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) ck(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) ec(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) taug(bin_grid_size(bin_grid)), taup(bin_grid_size(bin_grid))
    real(kind=dp) taul(bin_grid_size(bin_grid)), tauu(bin_grid_size(bin_grid))
    real(kind=dp) prod(bin_grid_size(bin_grid)), ploss(bin_grid_size(bin_grid))
    real(kind=dp) time, last_output_time, last_progress_time
    real(kind=dp) bin_vol_tot
    type(env_state_t) :: old_env_state
    type(aero_binned_t) :: aero_binned

    integer i, j, i_time, num_t, i_summary, n_spec
    integer ncid
    logical do_output, do_progress

    call check_time_multiple("t_max", run_sect_opt%t_max, &
         "del_t", run_sect_opt%del_t)
    call check_time_multiple("t_output", run_sect_opt%t_output, &
         "del_t", run_sect_opt%del_t)
    call check_time_multiple("t_progress", run_sect_opt%t_progress, &
         "del_t", run_sect_opt%del_t)

    ! g         : total spectral volume distribution (m^3/m^3)
    ! gs        : per-species spectral volume distribution (m^3/m^3)
    ! e         : single-particle volume grid (m^3)
    ! r         : droplet radius grid (um)
    ! log_width : constant grid distance of logarithmic grid
    !
    ! The transport is done in particle volume rather than mass because
    ! the destination bin of a coagulation depends on the (composition-
    ! independent) particle volume. For a single species this is
    ! algebraically identical to the original mass-based formulation.

    n_spec = aero_data_n_spec(aero_data)

    ! gas_state arrives already sized (zeroed when no gas-phase chemistry,
    ! or populated from gas_init when do_tchem is set)

    ! volume and radius grid
    do i = 1,bin_grid_size(bin_grid)
       r(i) = bin_grid%centers(i) * 1d6 ! radius in m to um
       e(i) = aero_data_rad2vol(aero_data, bin_grid%centers(i)) ! vol in m^3
    end do

    ! initial mass distribution: build it from the aero_dist ...
    call aero_binned_set_sizes(aero_binned, bin_grid_size(bin_grid), n_spec)
    call aero_binned_add_aero_dist(aero_binned, bin_grid, aero_data, &
         aero_dist)

    ! ... or read the aero_binned state directly from a NetCDF file
    !call pmc_nc_open_read(run_sect_opt%aero_init_filename, ncid)
    !call aero_binned_input_netcdf(aero_binned, ncid, bin_grid, &
    !     aero_data)
    !call pmc_nc_close(ncid)


    call courant(bin_grid_size(bin_grid), bin_grid%widths(1), e, ima, c)

    ! initialize MOSAIC (allocates its data structures; matches run_part)
    if (run_sect_opt%do_mosaic) then
#ifdef PMC_USE_MOSAIC
       call mosaic_init(env_state, aero_data, run_sect_opt%del_t, .false.)
#endif
    end if

    ! initialize time
    last_progress_time = 0d0
    time = 0d0
    i_summary = 1

    ! initial output
    call check_event(time, run_sect_opt%del_t, run_sect_opt%t_output, &
         last_output_time, do_output)
    if (do_output) then
       call output_sectional(run_sect_opt%prefix, bin_grid, aero_data, &
            aero_binned, gas_data, gas_state, env_state, i_summary, &
            time, run_sect_opt%t_output, run_sect_opt%uuid)
    end if

    ! main time-stepping loop
    num_t = nint(run_sect_opt%t_max / run_sect_opt%del_t)
    do i_time = 1, num_t

       if (run_sect_opt%do_coagulation) then
          if (SECT_COAG_METHOD == SECT_COAG_TWO_MOMENT) then
             ! Jacobson (2002) semi-implicit two-moment scheme: evolves the
             ! per-bin number and per-species volume as independent moments.
             call coag_two_moment(bin_grid, aero_data, env_state, &
                  run_sect_opt%coag_kernel_type, run_sect_opt%del_t, &
                  aero_binned)
          else
             ! per-bin mean composition (volume fractions) for the kernel,
             ! falling back to pure species 1 in empty bins
             do i = 1,bin_grid_size(bin_grid)
                bin_vol_tot = sum(aero_binned%vol_conc(i,:))
                if (bin_vol_tot > 0d0) then
                   bin_vol_frac(i,:) = aero_binned%vol_conc(i,:) / bin_vol_tot
                else
                   bin_vol_frac(i,:) = 0d0
                   bin_vol_frac(i,1) = 1d0
                end if
             end do

             ! recompute kernel for the current per-bin mean density
             call bin_kernel(bin_grid_size(bin_grid), bin_grid%centers, &
                  aero_data, run_sect_opt%coag_kernel_type, env_state, &
                  bin_vol_frac, k_bin)
             call smooth_bin_kernel(bin_grid_size(bin_grid), k_bin, ck)
             ! multiply kernel with constant timestep and logarithmic grid
             ! distance (kernel and volume grid are both in SI, so no unit
             ! conversion is needed)
             do i = 1,bin_grid_size(bin_grid)
                do j = 1,bin_grid_size(bin_grid)
                   ck(i,j) = ck(i,j) * run_sect_opt%del_t * bin_grid%widths(i)
                end do
             end do

             g = sum(aero_binned%vol_conc, dim=2)
             gs = aero_binned%vol_conc
             call coad(bin_grid_size(bin_grid), n_spec, run_sect_opt%del_t, &
                  taug, taup, taul, tauu, prod, ploss, c, ima, g, gs, r, e, &
                  ck, ec)
             aero_binned%vol_conc = gs
             aero_binned%num_conc = g &
                  / aero_data_rad2vol(aero_data, bin_grid%centers)
          end if
       end if

       time = run_sect_opt%t_max * real(i_time, kind=dp) &
            / real(num_t, kind=dp)

       old_env_state = env_state
       call scenario_update_env_state(scenario, env_state, time)
       call scenario_update_gas_state(scenario, run_sect_opt%del_t, &
            env_state, old_env_state, gas_data, gas_state)
       call scenario_update_aero_binned(scenario, run_sect_opt%del_t, &
            env_state, old_env_state, bin_grid, aero_data, aero_binned)

       if (run_sect_opt%do_tchem) then
#ifdef PMC_USE_TCHEM
          call pmc_tchem_interface_solve_sect(env_state, aero_data, bin_grid, &
               aero_binned, gas_data, gas_state, run_sect_opt%del_t)
#endif
       end if

       if (run_sect_opt%do_mosaic) then
#ifdef PMC_USE_MOSAIC
          call mosaic_timestep_sect(env_state, aero_data, bin_grid, &
               aero_binned, gas_data, gas_state)
#endif
       end if

       ! print output
       call check_event(time, run_sect_opt%del_t, run_sect_opt%t_output, &
            last_output_time, do_output)
       if (do_output) then
          i_summary = i_summary + 1
          call output_sectional(run_sect_opt%prefix, bin_grid, aero_data, &
               aero_binned, gas_data, gas_state, env_state, i_summary, &
               time, run_sect_opt%t_output, run_sect_opt%uuid)
       end if

       ! print progress to stdout
       call check_event(time, run_sect_opt%del_t, run_sect_opt%t_progress, &
            last_progress_time, do_progress)
       if (do_progress) then
          write(*,'(a6,a8)') 'step', 'time'
          write(*,'(i6,f8.1)') i_time, time
       end if
    end do

    if (run_sect_opt%do_mosaic) then
#ifdef PMC_USE_MOSAIC
       call mosaic_cleanup()
#endif
    end if

  end subroutine run_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the specification for a run_sect simulation from a spec file.
  subroutine spec_file_read_run_sect(file, run_sect_opt, aero_data, &
       bin_grid, gas_data, gas_state_init, env_state, aero_dist_init, scenario)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Options controlling the operation of run_sect().
    type(run_sect_opt_t), intent(inout) :: run_sect_opt
    !> Aerosol data.
    type(aero_data_t), intent(out) :: aero_data
    !> Bin grid.
    type(bin_grid_t), intent(out) :: bin_grid
    !> Initial aerosol state.
    type(aero_dist_t), intent(out) :: aero_dist_init
    !> Scenario data.
    type(scenario_t), intent(out) :: scenario
    !> Environmental state.
    type(env_state_t), intent(out) :: env_state
    !> Gas data.
    type(gas_data_t), intent(out) :: gas_data
    !> Initial gas state.
    type(gas_state_t), intent(out) :: gas_state_init

    character(len=PMC_MAX_FILENAME_LEN) :: sub_filename
    character(len=PMC_MAX_FILENAME_LEN) :: tchem_gas_filename, &
         tchem_aero_filename, tchem_numerics_filename
    integer :: n_grid_cells
    type(spec_file_t) :: sub_file

    call spec_file_read_string(file, 'output_prefix', run_sect_opt%prefix)

    call spec_file_read_real(file, 't_max', run_sect_opt%t_max)
    call spec_file_read_real(file, 'del_t', run_sect_opt%del_t)
    call spec_file_read_real(file, 't_output', run_sect_opt%t_output)
    call spec_file_read_real(file, 't_progress', run_sect_opt%t_progress)

    call spec_file_read_logical(file, 'do_tchem', run_sect_opt%do_tchem)
    if (run_sect_opt%do_tchem) then
#ifdef PMC_USE_TCHEM
       call spec_file_read_string(file, 'tchem_gas_config', &
            tchem_gas_filename)
       call spec_file_read_string(file, 'tchem_aero_config', &
            tchem_aero_filename)
       call spec_file_read_string(file, 'tchem_numerics_config', &
            tchem_numerics_filename)
#endif
    end if

    call spec_file_read_logical(file, 'do_mosaic', run_sect_opt%do_mosaic)
    if (run_sect_opt%do_mosaic .and. (.not. mosaic_support())) then
       call spec_file_die_msg(318492655, file, &
            'cannot use MOSAIC, support is not compiled in')
    end if
    if (run_sect_opt%do_mosaic .and. run_sect_opt%do_tchem) then
       call spec_file_die_msg(318492656, file, &
            'cannot use MOSAIC and TChem simultaneously')
    end if

    if (run_sect_opt%do_tchem) then
#ifdef PMC_USE_TCHEM
       n_grid_cells = 1
       call pmc_tchem_initialize(tchem_gas_filename, tchem_aero_filename, &
            tchem_numerics_filename, gas_data, aero_data, n_grid_cells)
#endif
    end if

    call spec_file_read_radius_bin_grid(file, bin_grid)

    ! gas and aerosol material data: when do_tchem is set these are provided
    ! by TChem (pmc_tchem_initialize, above) and are not read here
    if (.not. run_sect_opt%do_tchem) then
       call spec_file_read_string(file, 'gas_data', sub_filename)
       call spec_file_open(sub_filename, sub_file)
       call spec_file_read_gas_data(sub_file, gas_data)
       call spec_file_close(sub_file)
    end if

    ! initial gas state: only read when there is gas-phase chemistry (TChem or
    ! MOSAIC); otherwise start from a zero gas state
    if (run_sect_opt%do_tchem .or. run_sect_opt%do_mosaic) then
       call spec_file_read_string(file, 'gas_init', sub_filename)
       call spec_file_open(sub_filename, sub_file)
       call spec_file_read_gas_state(sub_file, gas_data, gas_state_init)
       call spec_file_close(sub_file)
    else
       call gas_state_set_size(gas_state_init, gas_data_n_spec(gas_data))
    end if

    if (.not. run_sect_opt%do_tchem) then
       call spec_file_read_string(file, 'aerosol_data', sub_filename)
       call spec_file_open(sub_filename, sub_file)
       call spec_file_read_aero_data(sub_file, aero_data)
       call spec_file_close(sub_file)
    end if

    call spec_file_read_fractal(file, aero_data%fractal)

    call spec_file_read_string(file, 'aerosol_init', sub_filename)
    call spec_file_open(sub_filename, sub_file)
    call spec_file_read_aero_dist(sub_file, aero_data, .false., aero_dist_init)
    call spec_file_close(sub_file)

    ! NetCDF file holding the initial aero_binned state, used when the
    ! aero_binned is read directly rather than built from aerosol_init.
    ! Disabled for now: the initial state is built from aerosol_init (the
    ! aero_dist) in run_sect(); see the matching commented-out read there.
    !call spec_file_read_string(file, 'aerosol_init_netcdf', &
    !     run_sect_opt%aero_init_filename)

    call spec_file_read_scenario(file, gas_data, aero_data, .false., scenario)
    call spec_file_read_env_state(file, env_state)

    call spec_file_read_logical(file, 'do_coagulation', &
         run_sect_opt%do_coagulation)
    if (run_sect_opt%do_coagulation) then
       call spec_file_read_coag_kernel_type(file, &
            run_sect_opt%coag_kernel_type)
       if (run_sect_opt%coag_kernel_type == COAG_KERNEL_TYPE_ADDITIVE) then
          call spec_file_read_real(file, 'additive_kernel_coeff', &
               env_state%additive_kernel_coefficient)
       end if
    else
       run_sect_opt%coag_kernel_type = COAG_KERNEL_TYPE_INVALID
    end if

    call spec_file_close(file)

    ! finished reading .spec data, now do the run

    call pmc_srand(0, 0)

    call uuid4_str(run_sect_opt%uuid)

    call scenario_init_env_state(scenario, env_state, 0d0)

  end subroutine spec_file_read_run_sect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Collision subroutine, exponential approach.
  !!
  !! Transports total particle volume \c g between bins with the Bott
  !! (1998) flux scheme, and carries the per-species volume \c gs along
  !! with it. Mass (volume) lost from a bin is removed in proportion to
  !! that bin's current composition; mass gained in a destination bin
  !! is added with its mixed composition and the advective flux to the
  !! next bin carries the destination bin's (post-gain) composition.
  !! For \c n_spec == 1 this reduces exactly to the original scheme.
  subroutine coad(n_bin, n_spec, dt, taug, taup, taul, tauu, prod, ploss, &
       c, ima, g, gs, r, e, ck, ec)

    integer n_bin
    integer n_spec
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
    real(kind=dp) gs(n_bin,n_spec)
    real(kind=dp) r(n_bin)
    real(kind=dp) e(n_bin)
    real(kind=dp) ck(n_bin,n_bin)
    real(kind=dp) ec(n_bin,n_bin)

    real(kind=dp), parameter :: gmin = 1d-60

    integer i, i0, i1, j, k, kp, i_spec
    real(kind=dp) x0, gsi, gsj, gsk, gk, x1, flux
    real(kind=dp) gain_s(n_spec), loss_i_s, loss_j_s, transfer_s

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

          ! loss from positions i, j (split by each donor bin's
          ! composition) and accumulate the per-species gain
          do i_spec = 1, n_spec
             loss_i_s = merge(gsi * gs(i,i_spec) / g(i), 0d0, g(i) > 0d0)
             loss_j_s = merge(gsj * gs(j,i_spec) / g(j), 0d0, g(j) > 0d0)
             gs(i,i_spec) = gs(i,i_spec) - loss_i_s
             gs(j,i_spec) = gs(j,i_spec) - loss_j_s
             gain_s(i_spec) = loss_i_s + loss_j_s
          end do
          ploss(i) = ploss(i) + gsi
          ploss(j) = ploss(j) + gsj
          g(i) = g(i) - gsi
          g(j) = g(j) - gsj

          if (k > 0) then ! do we have a valid bin for the coagulation result?
             gk = g(k) + gsk

             if (gk .gt. gmin) then
                ! add the gained volume into bin k, mixing composition
                do i_spec = 1, n_spec
                   gs(k,i_spec) = gs(k,i_spec) + gain_s(i_spec)
                end do
                g(k) = gk

                x1 = log(g(kp) / gk + 1d-60)
                flux = gsk / x1 * (exp(0.5d0 * x1) &
                     - exp(x1 * (0.5d0 - c(i,j))))
                flux = min(flux, gk)

                ! advect the flux from bin k to bin kp carrying bin k's
                ! current (post-gain) composition
                do i_spec = 1, n_spec
                   transfer_s = merge(flux * gs(k,i_spec) / g(k), 0d0, &
                        g(k) > 0d0)
                   gs(k,i_spec) = gs(k,i_spec) - transfer_s
                   gs(kp,i_spec) = gs(kp,i_spec) + transfer_s
                end do
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

  !> Two-moment coagulation, semi-implicit scheme of Jacobson (2002).
  !!
  !! Unlike the single-moment Bott flux scheme in coad() (which advects the
  !! volume distribution and derives number from it), this evolves the per-bin
  !! number concentration and the per-species volume concentration as
  !! independent moments. It is the single-type algorithm used by MOSAIC
  !! (jacobson2002_singletype_coag in module_mosaic_coag1d), adapted here to
  !! PartMC's aero_binned_t and coagulation kernel.
  !!
  !! For each donor pair (i,j) the combined single-particle volume
  !! e(i) + e(j) is split between the two grid bins k and k+1 that bracket it,
  !! with the volume-conserving weights f_ijk of Jacobson (2002) eqn 8. Total
  !! per-species volume is conserved to round-off; number decreases as
  !! particles merge. The update is semi-implicit (unconditionally stable), and
  !! the step is sub-cycled so that no bin loses more than a set fraction of its
  !! number per sub-step, which keeps the moments accurate for large del_t.
  !!
  !! \reference Jacobson, M. Z., Analysis of aerosol interactions with numerical
  !! techniques for solving coagulation..., J. Geophys. Res., 107(D19), 4366,
  !! 2002.
  subroutine coag_two_moment(bin_grid, aero_data, env_state, &
       coag_kernel_type, del_t, aero_binned)

    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Type of coagulation kernel.
    integer, intent(in) :: coag_kernel_type
    !> Timestep (s).
    real(kind=dp), intent(in) :: del_t
    !> Binned aerosol distribution (updated in place).
    type(aero_binned_t), intent(inout) :: aero_binned

    !> Maximum fractional number loss from a bin per sub-step.
    real(kind=dp), parameter :: frac_loss_limit = 0.5d0
    !> Cap on the number of sub-steps.
    integer, parameter :: max_nsubstep = 1000

    integer :: n_bin, n_spec, i, j, k, kp1, l, isubstep, nsubstep
    real(kind=dp) :: bin_vol_tot, vol_ipj, f_tmp, tmpa, tmpb, del_t_sub
    real(kind=dp) :: t1_num, t3_num, t1_vol_denom
    real(kind=dp) :: e(bin_grid_size(bin_grid))
    real(kind=dp) :: bin_vol_frac(bin_grid_size(bin_grid), &
         aero_data_n_spec(aero_data))
    real(kind=dp) :: k_bin(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) :: beta(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) :: cnum(bin_grid_size(bin_grid))
    real(kind=dp) :: cnum_old(bin_grid_size(bin_grid))
    real(kind=dp) :: cvol(bin_grid_size(bin_grid),aero_data_n_spec(aero_data))
    real(kind=dp) :: t1_vol(aero_data_n_spec(aero_data))
    real(kind=dp) :: f_ieqk(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) :: f_klo(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    real(kind=dp) :: f_khi(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    integer :: klo_of_ij(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    integer :: khi_of_ij(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    integer :: ilo_of_jk(bin_grid_size(bin_grid),bin_grid_size(bin_grid))
    integer :: ihi_of_jk(bin_grid_size(bin_grid),bin_grid_size(bin_grid))

    n_bin = bin_grid_size(bin_grid)
    n_spec = aero_data_n_spec(aero_data)
    if (n_bin < 1) return

    ! single-particle volume at each bin center (m^3)
    do i = 1,n_bin
       e(i) = aero_data_rad2vol(aero_data, bin_grid%centers(i))
    end do

    ! per-bin mean composition (volume fractions) for the kernel, falling
    ! back to pure species 1 in empty bins
    do i = 1,n_bin
       bin_vol_tot = sum(aero_binned%vol_conc(i,:))
       if (bin_vol_tot > 0d0) then
          bin_vol_frac(i,:) = aero_binned%vol_conc(i,:) / bin_vol_tot
       else
          bin_vol_frac(i,:) = 0d0
          bin_vol_frac(i,1) = 1d0
       end if
    end do

    ! coagulation kernel beta(i,j) (m^3/s) for the current per-bin composition
    call bin_kernel(n_bin, bin_grid%centers, aero_data, coag_kernel_type, &
         env_state, bin_vol_frac, k_bin)

    ! Work in actual concentrations (per m^3), i.e. per-log-width densities
    ! times the bin widths, because the coagulation rate is quadratic in
    ! concentration. The widths cancel when converting back at the end.
    do i = 1,n_bin
       cnum(i) = aero_binned%num_conc(i) * bin_grid%widths(i)
       cvol(i,:) = aero_binned%vol_conc(i,:) * bin_grid%widths(i)
    end do

    ! f_ijk of Jacobson (2002) eqn 8: for each donor pair (i,j) the combined
    ! volume e(i)+e(j) is split between bins klo and khi = klo+1. Only the two
    ! non-zero k are stored per (i,j).
    f_ieqk = 0d0
    f_klo = 0d0
    f_khi = 0d0
    klo_of_ij = 0
    khi_of_ij = 0
    do i = 1,n_bin
       do j = 1,n_bin
          vol_ipj = e(i) + e(j)
          if (vol_ipj >= e(n_bin)) then
             if (i == n_bin) f_ieqk(i,j) = 1d0
             klo_of_ij(i,j) = n_bin
             f_klo(i,j) = 1d0
          else
             do k = 1,n_bin - 1
                kp1 = k + 1
                if ((vol_ipj >= e(k)) .and. (vol_ipj < e(kp1))) then
                   f_tmp = ((e(kp1) - vol_ipj) / (e(kp1) - e(k))) &
                        * e(k) / vol_ipj
                   if (i == k) f_ieqk(i,j) = f_tmp
                   klo_of_ij(i,j) = k
                   f_klo(i,j) = f_tmp
                   khi_of_ij(i,j) = kp1
                   f_khi(i,j) = 1d0 - f_tmp
                end if
             end do
          end if
       end do
    end do

    ! for each (j,k), the range of donor i for which i + j --> k is possible
    do k = 1,n_bin
       do j = 1,k
          ilo_of_jk(j,k) = n_bin + 1
          ihi_of_jk(j,k) = 0
          do i = 1,k - 1
             if ((klo_of_ij(i,j) == k) .or. (khi_of_ij(i,j) == k)) then
                ilo_of_jk(j,k) = min(ilo_of_jk(j,k), i)
                ihi_of_jk(j,k) = max(ihi_of_jk(j,k), i)
             end if
          end do
       end do
    end do

    ! choose the number of sub-steps so that no bin loses more than
    ! frac_loss_limit of its number per sub-step (the scheme is stable for any
    ! step, but sub-cycling keeps the two moments accurate for large del_t)
    tmpa = 0d0
    do j = 1,n_bin
       tmpb = 0d0
       do i = 1,n_bin
          tmpb = tmpb + k_bin(i,j) * cnum(i)
       end do
       tmpb = tmpb - 0.5d0 * k_bin(j,j) * cnum(j)
       tmpa = max(tmpa, tmpb)
    end do
    nsubstep = 1
    if (tmpa * del_t > frac_loss_limit) then
       nsubstep = min(max_nsubstep, int(tmpa * del_t / frac_loss_limit) + 1)
    end if
    del_t_sub = del_t / real(nsubstep, kind=dp)

    ! kernel times the sub-step (dimensionless collision volume)
    beta = k_bin * del_t_sub

    ! semi-implicit solve over the sub-steps
    do isubstep = 1,nsubstep
       cnum_old = cnum

       do k = 1,n_bin
          ! T3 term: fractional loss out of bin k (eqns 7 and 9). The f_ieqk
          ! part that stays in k is excluded.
          t3_num = 0d0
          if (k < n_bin) then
             do j = 1,n_bin
                t3_num = t3_num &
                     + (1d0 - f_ieqk(k,j)) * beta(k,j) * cnum_old(j)
             end do
          end if

          ! T1 term: production into bin k from smaller donor pairs
          t1_num = 0d0
          t1_vol = 0d0
          if (k > 1) then
             do j = 1,k
                do i = ilo_of_jk(j,k), ihi_of_jk(j,k)
                   if (klo_of_ij(i,j) == k) then
                      f_tmp = f_klo(i,j)
                   else if (khi_of_ij(i,j) == k) then
                      f_tmp = f_khi(i,j)
                   else
                      cycle
                   end if
                   tmpa = f_tmp * beta(j,i) * cnum_old(j)
                   t1_num = t1_num + tmpa * cnum(i) * e(i)
                   do l = 1,n_spec
                      t1_vol(l) = t1_vol(l) + tmpa * cvol(i,l)
                   end do
                end do
             end do
             t1_num = t1_num / e(k)
          end if

          ! semi-implicit update (Jacobson 2002 eqns 7 and 9); the volume loss
          ! rate t3 is identical for number and volume
          t1_vol_denom = 1d0 + t3_num
          cnum(k) = (cnum(k) + t1_num) / (1d0 + t3_num)
          do l = 1,n_spec
             cvol(k,l) = (cvol(k,l) + t1_vol(l)) / t1_vol_denom
          end do
       end do
    end do

    ! convert actual concentrations back to per-log-width densities
    do i = 1,n_bin
       aero_binned%num_conc(i) = max(0d0, cnum(i)) / bin_grid%widths(i)
       aero_binned%vol_conc(i,:) = max(0d0, cvol(i,:)) / bin_grid%widths(i)
    end do

  end subroutine coag_two_moment

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
          ! this is basically the same as particle_in_bin(), but that
          ! is actually slightly different than what was always done
          ! here

          ! MW 2011-04-28: I think the above comment no longer
          ! applies, and we can make this just
          ! bin_grid_find(). FIXME.
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
