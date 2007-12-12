! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Data processing.

module pmc_output_processed

  use pmc_bin_grid
  use pmc_aero_data
  use pmc_aero_state
  use pmc_aero_binned
  use pmc_gas_data
  use pmc_gas_state
  use pmc_env_state
  use pmc_util
  use pmc_process_spec
  use pmc_netcdf
  use netcdf

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_processed_open(prefix, i_loop, ncid)

    ! Open the processed state output file.

    character(len=*), intent(in) :: prefix ! prefix of files to write
    integer, intent(in) :: i_loop       ! current loop number
    integer, intent(out) :: ncid        ! new NetCDF file ID, in data mode

    character(len=len(prefix)+20) :: filename
    character(len=500) :: history

    write(filename, '(a,a,i4.4,a)') trim(prefix), '_', i_loop, '.nc'
    call pmc_nc_check(nf90_create(filename, NF90_CLOBBER, ncid))

    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
         "PartMC output file"))
    call iso8601_date_and_time(history)
    history((len_trim(history)+1):) = " created by PartMC"
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "history", history))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine output_processed_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_processed_close(ncid)

    ! Close the processed state output file.

    integer, intent(out) :: ncid        ! new NetCDF file ID, in data mode

    call pmc_nc_check(nf90_close(ncid))

  end subroutine output_processed_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_processed(ncid, process_spec_list, bin_grid, &
       aero_data, aero_state, gas_data, gas_state, env, index, time, &
       del_t, i_loop)

    ! Write the current processed state.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(process_spec_t), intent(in) :: process_spec_list(:) ! processings specs
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    type(env_state_t), intent(in) :: env      ! environment state
    integer, intent(in) :: index        ! filename index
    real*8, intent(in) :: time          ! current time (s)
    real*8, intent(in) :: del_t         ! current output time-step (s)
    integer, intent(in) :: i_loop       ! current loop number

    integer :: i

    call process_time(ncid, time, index, del_t)
    do i = 1,size(process_spec_list)
       if (process_spec_list(i)%type == "env") then
          call process_env(ncid, process_spec_list(i)%suffix, &
               time, index, env)
       elseif (process_spec_list(i)%type == "gas") then
          call process_gas(ncid, process_spec_list(i)%suffix, &
               time, index, gas_data, gas_state)
       elseif (process_spec_list(i)%type == "aero") then
          call process_aero(ncid, process_spec_list(i)%suffix, &
               time, index, bin_grid, aero_data, aero_state)
       elseif ((process_spec_list(i)%type == "kappa") &
          .or. (process_spec_list(i)%type == "comp") &
          .or. (process_spec_list(i)%type == "n_orig_part") &
          .or. (process_spec_list(i)%type == "optic_absorb") &
          .or. (process_spec_list(i)%type == "optic_scatter") &
          .or. (process_spec_list(i)%type == "optic_extinct")) then
          call process_hist_new(ncid, time, index, bin_grid, &
               env, aero_data, aero_state, process_spec_list(i))
       else
          call die(450985234)
       end if
    end do

  end subroutine output_processed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_processed_binned(ncid, process_spec_list, &
       bin_grid, aero_data, aero_binned, gas_data, gas_state, env, &
       index, time, del_t)

    ! Write the current binned data.
    
    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(process_spec_t), intent(in) :: process_spec_list(:) ! processings specs
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_binned_t), intent(in) :: aero_binned ! binned aerosol data
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state
    type(env_state_t), intent(in) :: env      ! environment state
    integer, intent(in) :: index        ! filename index
    real*8, intent(in) :: time          ! current time (s)
    real*8, intent(in) :: del_t         ! current output time-step (s)

    integer :: i

    call process_time(ncid, time, index, del_t)
    do i = 1,size(process_spec_list)
       if (process_spec_list(i)%type == "env") then
          call process_env(ncid, process_spec_list(i)%suffix, &
               time, index, env)
       elseif (process_spec_list(i)%type == "gas") then
          call process_gas(ncid, process_spec_list(i)%suffix, &
               time, index, gas_data, gas_state)
       elseif (process_spec_list(i)%type == "aero") then
          call output_aero(ncid, process_spec_list(i)%suffix, &
               time, index, bin_grid, aero_data, aero_binned)
       end if
    end do

  end subroutine output_processed_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_dim_time(ncid, dimid_time)

    ! Write the time dimension to the given NetCDF file if it is not
    ! already present and in any case return the associated dimid.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    integer, intent(out) :: dimid_time  ! dimid of the time dimension

    integer :: status, varid_time, varid_time_widths

    status = nf90_inq_dimid(ncid, "time", dimid_time)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "time", NF90_UNLIMITED, dimid_time))
    call pmc_nc_check(nf90_def_var(ncid, "time", NF90_DOUBLE, &
         dimid_time, varid_time))
    call pmc_nc_check(nf90_put_att(ncid, varid_time, "unit", "s"))
    call pmc_nc_check(nf90_def_var(ncid, "time_widths", NF90_DOUBLE, &
         dimid_time, varid_time_widths))
    call pmc_nc_check(nf90_put_att(ncid, varid_time_widths, "unit", "s"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_dim_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_dim_radius(ncid, bin_grid, dimid_radius)

    ! Write the radius dimension to the given NetCDF file if it is not
    ! already present and in any case return the associated dimid.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    integer, intent(out) :: dimid_radius ! dimid of the radius dimension

    integer :: status, i_bin, varid_radius
    integer :: dimid_radius_edges, varid_radius_edges, varid_radius_widths
    real*8 :: radius_centers(bin_grid%n_bin), radius_edges(bin_grid%n_bin + 1)
    real*8 :: radius_widths(bin_grid%n_bin)

    status = nf90_inq_dimid(ncid, "radius", dimid_radius)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "radius", &
         bin_grid%n_bin, dimid_radius))
    call pmc_nc_check(nf90_def_dim(ncid, "radius_edges", &
         bin_grid%n_bin + 1, dimid_radius_edges))
    call pmc_nc_check(nf90_def_var(ncid, "radius", NF90_DOUBLE, &
         dimid_radius, varid_radius))
    call pmc_nc_check(nf90_put_att(ncid, varid_radius, "unit", "m"))
    call pmc_nc_check(nf90_def_var(ncid, "radius_edges", NF90_DOUBLE, &
         dimid_radius_edges, varid_radius_edges))
    call pmc_nc_check(nf90_put_att(ncid, varid_radius_edges, "unit", "m"))
    call pmc_nc_check(nf90_def_var(ncid, "radius_widths", NF90_DOUBLE, &
         dimid_radius, varid_radius_widths))
    call pmc_nc_check(nf90_put_att(ncid, varid_radius_widths, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_bin = 1,bin_grid%n_bin
       radius_centers(i_bin) = vol2rad(bin_grid%v(i_bin))
       radius_widths(i_bin) = bin_grid%dlnr
    end do
    do i_bin = 1,(bin_grid%n_bin + 1)
       radius_edges(i_bin) = vol2rad(bin_edge(bin_grid, i_bin))
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_radius, radius_centers))
    call pmc_nc_check(nf90_put_var(ncid, varid_radius_edges, radius_edges))
    call pmc_nc_check(nf90_put_var(ncid, varid_radius_widths, radius_widths))

  end subroutine ensure_nc_dim_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_dim_aero_species(ncid, aero_data, dimid_aero_species)

    ! Write the aero species dimension to the given NetCDF file if it
    ! is not already present and in any case return the associated
    ! dimid.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(out) :: dimid_aero_species ! dimid of the species dimension

    integer :: status, i_spec
    integer :: varid_aero_species, varid_aero_species_widths
    integer :: aero_species_centers(aero_data%n_spec)
    real*8 :: aero_species_widths(aero_data%n_spec)
    character(len=(AERO_NAME_LEN * aero_data%n_spec)) :: aero_species_names

    status = nf90_inq_dimid(ncid, "aero_species", dimid_aero_species)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "aero_species", &
         aero_data%n_spec, dimid_aero_species))
    aero_species_names = ""
    do i_spec = 1,aero_data%n_spec
       aero_species_names((len_trim(aero_species_names) + 1):) &
            = trim(aero_data%name(i_spec))
       if (i_spec < aero_data%n_spec) then
          aero_species_names((len_trim(aero_species_names) + 1):) = ","
       end if
    end do
    call pmc_nc_check(nf90_def_var(ncid, "aero_species", NF90_INT, &
         dimid_aero_species, varid_aero_species))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_species, "unit", "1"))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_species, "names", &
         aero_species_names))
    call pmc_nc_check(nf90_def_var(ncid, "aero_species_widths", NF90_DOUBLE, &
         dimid_aero_species, varid_aero_species_widths))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero_species_widths, &
         "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_spec = 1,aero_data%n_spec
       aero_species_centers(i_spec) = i_spec
       aero_species_widths(i_spec) = 1d0
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_species, &
         aero_species_centers))
    call pmc_nc_check(nf90_put_var(ncid, varid_aero_species_widths, &
         aero_species_widths))

  end subroutine ensure_nc_dim_aero_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_dim_gas_species(ncid, gas_data, dimid_gas_species)

    ! Write the gas species dimension to the given NetCDF file if it
    ! is not already present and in any case return the associated
    ! dimid.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(gas_data_t), intent(in) :: gas_data ! gas_data structure
    integer, intent(out) :: dimid_gas_species ! dimid of the species dimension

    integer :: status, i_spec
    integer :: varid_gas_species, varid_gas_species_widths
    integer :: gas_species_centers(gas_data%n_spec)
    real*8 :: gas_species_widths(gas_data%n_spec)
    character(len=(GAS_NAME_LEN * gas_data%n_spec)) :: gas_species_names

    status = nf90_inq_dimid(ncid, "gas_species", dimid_gas_species)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "gas_species", &
         gas_data%n_spec, dimid_gas_species))
    gas_species_names = ""
    do i_spec = 1,gas_data%n_spec
       gas_species_names((len_trim(gas_species_names) + 1):) &
            = trim(gas_data%name(i_spec))
       if (i_spec < gas_data%n_spec) then
          gas_species_names((len_trim(gas_species_names) + 1):) = ","
       end if
    end do
    call pmc_nc_check(nf90_def_var(ncid, "gas_species", NF90_INT, &
         dimid_gas_species, varid_gas_species))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species, "unit", "1"))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species, "names", &
         gas_species_names))
    call pmc_nc_check(nf90_def_var(ncid, "gas_species_widths", NF90_DOUBLE, &
         dimid_gas_species, varid_gas_species_widths))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species_widths, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_spec = 1,gas_data%n_spec
       gas_species_centers(i_spec) = i_spec
       gas_species_widths(i_spec) = 1d0
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_gas_species, &
         gas_species_centers))
    call pmc_nc_check(nf90_put_var(ncid, varid_gas_species_widths, &
         gas_species_widths))

  end subroutine ensure_nc_dim_gas_species

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_dim_unit(ncid, dimid_unit)

    ! Write the unit dimension to the given NetCDF file if it is not
    ! already present and in any case return the associated dimid.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    integer, intent(out) :: dimid_unit  ! dimid of the unit dimension

    integer :: status, varid_unit, unit_centers(4)
    character(len=(4*30)) :: unit_names, unit_units

    status = nf90_inq_dimid(ncid, "unit", dimid_unit)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "unit", &
         4, dimid_unit))
    unit_names = "num_den,vol_den,mass_den,mole_den"
    unit_units = "#/m^3,m^3/m^3,kg/m^3,mole/m^3"
    call pmc_nc_check(nf90_def_var(ncid, "unit", NF90_INT, &
         dimid_unit, varid_unit))
    call pmc_nc_check(nf90_put_att(ncid, varid_unit, "unit", "1"))
    call pmc_nc_check(nf90_put_att(ncid, varid_unit, "names", unit_names))
    call pmc_nc_check(nf90_put_att(ncid, varid_unit, "data_units", unit_units))

    call pmc_nc_check(nf90_enddef(ncid))

    unit_centers = (/ 1, 2, 3, 4 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_unit, unit_centers))

  end subroutine ensure_nc_dim_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_dim_env(ncid, dimid_env)

    ! Write the env dimension to the given NetCDF file if it is not
    ! already present and in any case return the associated dimid.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    integer, intent(out) :: dimid_env  ! dimid of the env dimension

    integer :: status, varid_env, env_centers(4)
    character(len=(4*30)) :: env_names, env_units

    status = nf90_inq_dimid(ncid, "env", dimid_env)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "env", &
         4, dimid_env))
    env_names = "temp,rel_humid,pressure,height"
    env_units = "K,1,Pa,m"
    call pmc_nc_check(nf90_def_var(ncid, "env", NF90_INT, &
         dimid_env, varid_env))
    call pmc_nc_check(nf90_put_att(ncid, varid_env, "unit", "1"))
    call pmc_nc_check(nf90_put_att(ncid, varid_env, "names", env_names))
    call pmc_nc_check(nf90_put_att(ncid, varid_env, "data_units", env_units))

    call pmc_nc_check(nf90_enddef(ncid))

    env_centers = (/ 1, 2, 3, 4 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_env, env_centers))

  end subroutine ensure_nc_dim_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_dim_step(ncid, process_spec, dimid_step)

    ! Write the step dimension for the given process_spec to the given
    ! NetCDF file if it is not already present and in any case return
    ! the associated dimid.

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(process_spec_t), intent(in) :: process_spec ! process spec
    integer, intent(out) :: dimid_step  ! dimid of the step dimension

    integer :: status, i_step, varid_step
    integer :: dimid_step_edges, varid_step_edges, varid_step_widths
    real*8 :: step_centers(process_spec%n_step)
    real*8 :: step_edges(process_spec%n_step + 1)
    real*8 :: step_widths(process_spec%n_step)
    character(len=100) :: dim_name, dim_name_edges, dim_name_widths, dim_unit

    if (process_spec%type == "kappa") then
       dim_name = 'critical_supersat'
       dim_unit = '1'
    elseif (process_spec%type == "comp") then
       dim_name = 'composition'
       dim_unit = '1'
    elseif (process_spec%type == "n_orig_part") then
       dim_name = 'n_orig_part'
       dim_unit = '1'
    elseif (process_spec%type == "optic_absorb") then
       dim_name = 'absorb_cross_section_area'
       dim_unit = 'm^2'
    elseif (process_spec%type == "optic_scatter") then
       dim_name = 'scatter_cross_section_area'
       dim_unit = 'm^2'
    elseif (process_spec%type == "optic_extinct") then
       dim_name = 'extinct_cross_section_area'
       dim_unit = 'm^2'
    else
       call die_msg(912387902, &
            "unknown process_spec%type: " // process_spec%type)
    end if
    dim_name_edges = dim_name
    dim_name_edges((len_trim(dim_name_edges)+1):) = "_edges"
    dim_name_widths = dim_name
    dim_name_widths((len_trim(dim_name_widths)+1):) = "_widths"

    status = nf90_inq_dimid(ncid, dim_name, dimid_step)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, dim_name, &
         process_spec%n_step, dimid_step))
    call pmc_nc_check(nf90_def_dim(ncid, dim_name_edges, &
         process_spec%n_step + 1, dimid_step_edges))
    call pmc_nc_check(nf90_def_var(ncid, dim_name, NF90_DOUBLE, &
         dimid_step, varid_step))
    call pmc_nc_check(nf90_put_att(ncid, varid_step, "unit", dim_unit))
    call pmc_nc_check(nf90_def_var(ncid, dim_name_edges, NF90_DOUBLE, &
         dimid_step_edges, varid_step_edges))
    call pmc_nc_check(nf90_put_att(ncid, varid_step_edges, "unit", dim_unit))
    call pmc_nc_check(nf90_def_var(ncid, dim_name_widths, NF90_DOUBLE, &
         dimid_step, varid_step_widths))
    call pmc_nc_check(nf90_put_att(ncid, varid_step_widths, "unit", dim_unit))

    call pmc_nc_check(nf90_enddef(ncid))

    ! FIXME: add a logspace_with_centers() etc, and also use for bin_grid
    if (process_spec%log_scale) then
       call logspace(process_spec%min_val, process_spec%max_val, &
            process_spec%n_step + 1, step_edges)
       do i_step = 1,process_spec%n_step
          step_centers(i_step) = sqrt(step_edges(i_step) &
               * step_edges(i_step + 1))
          step_widths(i_step) = (log10(process_spec%max_val) &
               - log10(process_spec%min_val)) / dble(process_spec%n_step)
       end do
    else
       call linspace(process_spec%min_val, process_spec%max_val, &
            process_spec%n_step + 1, step_edges)
       do i_step = 1,process_spec%n_step
          step_centers(i_step) = (step_edges(i_step) &
               + step_edges(i_step + 1)) / 2d0
          step_widths(i_step) = (process_spec%max_val - process_spec%min_val) &
               / dble(process_spec%n_step)
       end do
    end if
    call pmc_nc_check(nf90_put_var(ncid, varid_step, step_centers))
    call pmc_nc_check(nf90_put_var(ncid, varid_step_edges, step_edges))
    call pmc_nc_check(nf90_put_var(ncid, varid_step_widths, step_widths))

  end subroutine ensure_nc_dim_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_var_env_state(ncid, varid_env_state)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    integer, intent(out) :: varid_env_state ! varid of env_state

    integer :: dimid_time, dimid_env, dimids_env_state(2), status

    status = nf90_inq_varid(ncid, "env_state", varid_env_state)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    call ensure_nc_dim_env(ncid, dimid_env)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_env_state = (/ dimid_env, dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, "env_state", NF90_DOUBLE, &
         dimids_env_state, varid_env_state))
    call pmc_nc_check(nf90_put_att(ncid, varid_env_state, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_var_gas(ncid, gas_data, varid_gas)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(gas_data_t), intent(in) :: gas_data ! gas data
    integer, intent(out) :: varid_gas   ! varid of gas

    type(inout_file_t) :: file
    integer :: i_spec, status

    integer :: dimid_time, dimid_gas_species, dimids_gas(2)

    status = nf90_inq_varid(ncid, "gas", varid_gas)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    call ensure_nc_dim_gas_species(ncid, gas_data, dimid_gas_species)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_gas = (/ dimid_gas_species, dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, "gas", NF90_DOUBLE, &
         dimids_gas, varid_gas))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas, "unit", "ppb"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_var_aero(ncid, bin_grid, aero_data, varid_aero)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(out) :: varid_aero  ! varid of aero

    integer :: dimid_time, dimid_radius, dimid_aero_species, dimid_unit
    integer :: status, dimids_aero(4)

    status = nf90_inq_varid(ncid, "aero", varid_aero)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    call ensure_nc_dim_radius(ncid, bin_grid, dimid_radius)
    call ensure_nc_dim_aero_species(ncid, aero_data, dimid_aero_species)
    call ensure_nc_dim_unit(ncid, dimid_unit)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_aero = (/ dimid_radius, dimid_aero_species, dimid_unit, dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, "aero", NF90_DOUBLE, &
         dimids_aero, varid_aero))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_var_hist(ncid, process_spec, bin_grid, &
       aero_data, varid_hist)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    type(process_spec_t), intent(in) :: process_spec ! process spec
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    integer, intent(out) :: varid_hist  ! varid of hist

    integer :: dimid_time, dimid_step, dimid_radius, dimid_aero_species
    integer :: dimid_unit, status, dimids_hist(5)

    status = nf90_inq_varid(ncid, process_spec%suffix, varid_hist)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    call ensure_nc_dim_step(ncid, process_spec, dimid_step)
    call ensure_nc_dim_radius(ncid, bin_grid, dimid_radius)
    call ensure_nc_dim_aero_species(ncid, aero_data, dimid_aero_species)
    call ensure_nc_dim_unit(ncid, dimid_unit)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_hist = (/ dimid_step, dimid_radius, dimid_aero_species, &
         dimid_unit, dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, process_spec%suffix, NF90_DOUBLE, &
         dimids_hist, varid_hist))
    call pmc_nc_check(nf90_put_att(ncid, varid_hist, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_hist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_time(ncid, time, index, del_t)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    real*8, intent(in) :: del_t         ! output timestep of current time

    integer :: dimid_time, varid_time, varid_time_widths

    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_inq_varid(ncid, "time", varid_time))
    call pmc_nc_check(nf90_put_var(ncid, varid_time, (/time/), &
         start = (/index/), count = (/1/)))
    call pmc_nc_check(nf90_inq_varid(ncid, "time_widths", varid_time_widths))
    call pmc_nc_check(nf90_put_var(ncid, varid_time_widths, (/del_t/), &
         start = (/index/), count = (/1/)))

  end subroutine process_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_env(ncid, suffix, time, index, env)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    character(len=*), intent(in) :: suffix ! suffix of the file
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(env_state_t), intent(in) :: env      ! environment state

    integer :: varid_env_state
    integer :: start(2), count(2)
    real*8 :: data(4)

    call ensure_nc_var_env_state(ncid, varid_env_state)

    data(1) = env%temp
    data(2) = env%rel_humid
    data(3) = env%pressure
    data(4) = env%height

    start = (/ 1, index /)
    count = (/ 4, 1 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_env_state, data, &
         start = start, count = count))

  end subroutine process_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_gas(ncid, suffix, time, index, gas_data, gas_state)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    character(len=*), intent(in) :: suffix ! suffix of the file
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(gas_data_t), intent(in) :: gas_data ! gas data
    type(gas_state_t), intent(in) :: gas_state ! gas state

    integer :: varid_gas
    integer :: start(2), count(2)

    call ensure_nc_var_gas(ncid, gas_data, varid_gas)

    start = (/ 1, index /)
    count = (/ gas_data%n_spec, 1 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_gas, gas_state%conc, &
         start = start, count = count))

  end subroutine process_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_aero(ncid, suffix, time, index, bin_grid, &
       aero_data, aero_binned)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    character(len=*), intent(in) :: suffix ! suffix for the output filename
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_binned_t), intent(in) :: aero_binned ! aero_binned structure

    real*8 :: aero(bin_grid%n_bin, aero_data%n_spec, 4)
    integer :: i_bin, i_spec, varid_aero
    integer :: start(4), count(4)

    call ensure_nc_var_aero(ncid, bin_grid, aero_data, varid_aero)

    do i_bin = 1,bin_grid%n_bin
       do i_spec = 1,aero_data%n_spec
          aero(i_bin, i_spec, 1) = &
               aero_binned%num_den(i_bin) / dble(aero_data%n_spec)
          aero(i_bin, i_spec, 2) = &
               aero_binned%vol_den(i_bin, i_spec)
          aero(i_bin, i_spec, 3) = &
               aero_binned%vol_den(i_bin, i_spec) &
               * aero_data%density(i_spec)
          aero(i_bin, i_spec, 4) = &
               aero_binned%vol_den(i_bin, i_spec) &
               * aero_data%density(i_spec) &
               / aero_data%molec_weight(i_spec)
       end do
    end do

    start = (/ 1, 1, 1, index /)
    count = (/ bin_grid%n_bin, aero_data%n_spec, 4, 1 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_aero, aero, &
         start = start, count = count))

  end subroutine output_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_aero(ncid, suffix, time, index, bin_grid, &
       aero_data, aero_state)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    character(len=*), intent(in) :: suffix ! suffix for the output filename
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(aero_state_t), intent(in) :: aero_state ! aero_state structure

    type(aero_binned_t) :: aero_binned

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)

    call output_aero(ncid, suffix, time, index, bin_grid, &
         aero_data, aero_binned)

    call aero_binned_free(aero_binned)

  end subroutine process_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_hist(ncid, time, index, bin_grid, aero_data, &
       process_spec, hist)

    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(bin_grid_t), intent(in) :: bin_grid ! bin_grid structure
    type(aero_data_t), intent(in) :: aero_data ! aero_data structure
    type(process_spec_t), intent(in) :: process_spec ! process spec
    real*8, intent(in) :: hist(process_spec%n_step, bin_grid%n_bin, &
         aero_data%n_spec, 4)           ! histogram data

    integer :: varid_hist
    integer :: start(5), count(5)

    call ensure_nc_var_hist(ncid, process_spec, bin_grid, aero_data, &
         varid_hist)

    start = (/ 1, 1, 1, 1, index /)
    count = (/ process_spec%n_step, bin_grid%n_bin, aero_data%n_spec, 4, 1 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_hist, hist, &
         start = start, count = count))

  end subroutine output_hist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_hist_new(ncid, time, index, bin_grid, &
       env, aero_data, aero_state, process_spec)

    ! Compute histogram by calling the step_comp() function on each
    ! particle.
    
    integer, intent(in) :: ncid         ! NetCDF file ID, in data mode
    real*8, intent(in) :: time          ! current time (s)
    integer, intent(in) :: index        ! current index
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(env_state_t), intent(in) :: env      ! environment state
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    type(aero_state_t), intent(in) :: aero_state ! aerosol state
    type(process_spec_t), intent(in) :: process_spec ! process spec

    real*8 :: hist(process_spec%n_step, bin_grid%n_bin, aero_data%n_spec, 4)
    integer :: i_step, i_bin, i_part, i
    type(aero_particle_t), pointer :: aero_particle
    real*8 :: step_width, scale, rh, supersat, val
    integer, allocatable :: a_species(:), b_species(:)
    character(len=200) :: error_str

    if (process_spec%log_scale) then
       step_width = (log10(process_spec%max_val) &
            - log10(process_spec%min_val)) / dble(process_spec%n_step)
    else
       step_width = (process_spec%max_val - process_spec%min_val) &
            / dble(process_spec%n_step)
    end if

    if (process_spec%type == "comp") then
       allocate(a_species(size(process_spec%a_species)))
       do i = 1,size(process_spec%a_species)
          a_species(i) = aero_data_spec_by_name(aero_data, &
               process_spec%a_species(i))
          if (a_species(i) == 0) then
             write(error_str, '(a,a)') 'unknown species: ', &
                  process_spec%a_species(i)
             call die_msg(194029329, error_str)
          end if
       end do
       allocate(b_species(size(process_spec%b_species)))
       do i = 1,size(process_spec%b_species)
          b_species(i) = aero_data_spec_by_name(aero_data, &
               process_spec%b_species(i))
          if (b_species(i) == 0) then
             write(error_str, '(a,a)') 'unknown species: ', &
                  process_spec%b_species(i)
             call die_msg(283298183, error_str)
          end if
       end do
    end if

    hist = 0d0
    scale = 1d0 / bin_grid%dlnr / step_width / aero_state%comp_vol
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bins(i_bin)%n_part
          aero_particle => aero_state%bins(i_bin)%particle(i_part)

          if (process_spec%type == "kappa") then
             rh = aero_particle_kappa_rh(aero_particle, aero_data, env)
             supersat = rh - 1d0
             val = supersat
          elseif (process_spec%type == "comp") then
             val = aero_particle_comp(aero_particle, a_species, b_species)
          elseif (process_spec%type == "n_orig_part") then
             val = dble(aero_particle%n_orig_part) + 0.5d0
          elseif (process_spec%type == "optic_absorb") then
             val = aero_particle%absorb_cross_sect
          elseif (process_spec%type == "optic_scatter") then
             val = aero_particle%scatter_cross_sect
          elseif (process_spec%type == "optic_extinct") then
             val = aero_particle%absorb_cross_sect &
                  + aero_particle%scatter_cross_sect
          else
             call die(298420915)
          end if

          if (process_spec%log_scale) then
             i_step = logspace_find(process_spec%min_val, &
                  process_spec%max_val, process_spec%n_step + 1, val)
          else
             i_step = linspace_find(process_spec%min_val, &
                  process_spec%max_val, process_spec%n_step + 1, val)
          end if
          i_step = max(1, i_step)
          i_step = min(process_spec%n_step, i_step)

          hist(i_step, i_bin, :, 1) = hist(i_step, i_bin, :, 1) &
               + 1d0 / dble(aero_data%n_spec) * scale
          hist(i_step, i_bin, :, 2) = hist(i_step, i_bin, :, 2) &
               + aero_particle%vol * scale
          hist(i_step, i_bin, :, 3) = hist(i_step, i_bin, :, 3) &
               + aero_particle%vol * aero_data%density * scale
          hist(i_step, i_bin, :, 4) = hist(i_step, i_bin, :, 4) &
               + aero_particle%vol * aero_data%density &
               / aero_data%molec_weight * scale
       end do
    end do

    if (process_spec%type == "comp") then
       deallocate(a_species, b_species)
    end if

    call output_hist(ncid, time, index, bin_grid, aero_data, &
         process_spec, hist)

  end subroutine process_hist_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_comp(aero_particle, a_species, b_species)

    type(aero_particle_t), intent(in) :: aero_particle ! particle
    integer, intent(in) :: a_species(:) ! first list of species
    integer, intent(in) :: b_species(:) ! second list of species

    integer :: i
    real*8 :: a_vol, b_vol

    ! FIXME: can we just do sum(vol(a_species))?
    a_vol = 0d0
    do i = 1,size(a_species)
       a_vol = a_vol + aero_particle%vol(a_species(i))
    end do
    b_vol = 0d0
    do i = 1,size(b_species)
       b_vol = b_vol + aero_particle%vol(b_species(i))
    end do
    call assert(880038232, a_vol >= 0d0)
    call assert(715496111, b_vol >= 0d0)
    if ((a_vol == 0d0) .and. (b_vol == 0d0)) then
       aero_particle_comp = 0d0
    else
       aero_particle_comp = b_vol / (a_vol + b_vol)
    end if

  end function aero_particle_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_output_processed
