! Copyright (C) 2007, 2008 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_output_processed module.

!> Transform data into arrays for output to a NetCDF file.
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

  !> Open the processed state output file.
  !!
  !! The filename is of the form \c prefix_loop.nc if \c i_loop is
  !! positive, otherwise it is \c prefix.nc.
  subroutine output_processed_open(prefix, i_loop, ncid)

    !> Prefix of files to write.
    character(len=*), intent(in) :: prefix
    !> Current loop number, or 0 to ignore the loop number.
    integer, intent(in) :: i_loop
    !> New NetCDF file ID, in data mode.
    integer, intent(out) :: ncid

    character(len=len(prefix)+20) :: filename
    character(len=500) :: history

    if (i_loop > 0) then
       write(filename, '(a,a,i4.4,a)') trim(prefix), '_', i_loop, '.nc'
    else
       write(filename, '(a,a)') trim(prefix), '.nc'
    end if
    call pmc_nc_check(nf90_create(filename, NF90_CLOBBER, ncid))

    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "title", &
         "PartMC output file"))
    call iso8601_date_and_time(history)
    history((len_trim(history)+1):) = " created by PartMC"
    call pmc_nc_check(nf90_put_att(ncid, NF90_GLOBAL, "history", history))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine output_processed_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Close the processed state output file.
  subroutine output_processed_close(ncid)

    !> New NetCDF file ID, in data mode.
    integer, intent(out) :: ncid

    call pmc_nc_check(nf90_close(ncid))

  end subroutine output_processed_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current processed state.
  subroutine output_processed(ncid, process_spec_list, bin_grid, &
       aero_data, aero_state, gas_data, gas_state, env_state, index, time, &
       del_t, i_loop)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Process specs.
    type(process_spec_t), intent(in) :: process_spec_list(:)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current output time-step (s).
    real*8, intent(in) :: del_t
    !> Current loop number.
    integer, intent(in) :: i_loop

    integer :: i

    call process_time(ncid, time, index, del_t)
    do i = 1,size(process_spec_list)
       if (process_spec_list(i)%type == "env") then
          call process_env(ncid, time, index, env_state, &
               process_spec_list(i))
       elseif (process_spec_list(i)%type == "gas") then
          call process_gas(ncid, time, index, gas_data, gas_state, &
               process_spec_list(i))
       elseif (process_spec_list(i)%type == "aero") then
          call process_aero(ncid, time, index, bin_grid, aero_data, &
               aero_state, process_spec_list(i))
       elseif ((process_spec_list(i)%type == "kappa") &
          .or. (process_spec_list(i)%type == "comp") &
          .or. (process_spec_list(i)%type == "n_orig_part") &
          .or. (process_spec_list(i)%type == "optic_absorb") &
          .or. (process_spec_list(i)%type == "optic_scatter") &
          .or. (process_spec_list(i)%type == "optic_extinct")) then
          call process_hist(ncid, time, index, bin_grid, &
               env_state, aero_data, aero_state, process_spec_list(i))
       else
          call die(450985234)
       end if
    end do

  end subroutine output_processed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the current binned data.
  subroutine output_processed_binned(ncid, process_spec_list, &
       bin_grid, aero_data, aero_binned, gas_data, gas_state, env_state, &
       index, time, del_t)
    
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Process specs.
    type(process_spec_t), intent(in) :: process_spec_list(:)
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Binned aerosol data.
    type(aero_binned_t), intent(in) :: aero_binned
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Filename index.
    integer, intent(in) :: index
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current output time-step (s).
    real*8, intent(in) :: del_t

    integer :: i

    call process_time(ncid, time, index, del_t)
    do i = 1,size(process_spec_list)
       if (process_spec_list(i)%type == "env") then
          call process_env(ncid, time, index, env_state, &
               process_spec_list(i))
       elseif (process_spec_list(i)%type == "gas") then
          call process_gas(ncid, time, index, gas_data, gas_state, &
               process_spec_list(i))
       elseif (process_spec_list(i)%type == "aero") then
          call output_aero(ncid, time, index, bin_grid, aero_data, &
               aero_binned, process_spec_list(i))
       end if
    end do

  end subroutine output_processed_binned
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the time dimension to the given NetCDF file if it is not
  !> already present and in any case return the associated dimid.
  subroutine ensure_nc_dim_time(ncid, dimid_time)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the time dimension.
    integer, intent(out) :: dimid_time

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

  !> Write the radius dimension to the given NetCDF file if it is not
  !> already present and in any case return the associated dimid.
  subroutine ensure_nc_dim_radius(ncid, bin_grid, dimid_radius)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Dimid of the radius dimension.
    integer, intent(out) :: dimid_radius

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
       radius_edges(i_bin) = vol2rad(bin_grid_edge(bin_grid, i_bin))
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_radius, radius_centers))
    call pmc_nc_check(nf90_put_var(ncid, varid_radius_edges, radius_edges))
    call pmc_nc_check(nf90_put_var(ncid, varid_radius_widths, radius_widths))

  end subroutine ensure_nc_dim_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the dry_radius dimension to the given NetCDF file if it is not
  !> already present and in any case return the associated dimid.
  subroutine ensure_nc_dim_radius_dry(ncid, bin_grid, dimid_radius)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Dimid of the radius dimension.
    integer, intent(out) :: dimid_radius

    integer :: status, i_bin, varid_radius
    integer :: dimid_radius_edges, varid_radius_edges, varid_radius_widths
    real*8 :: radius_centers(bin_grid%n_bin), radius_edges(bin_grid%n_bin + 1)
    real*8 :: radius_widths(bin_grid%n_bin)

    status = nf90_inq_dimid(ncid, "dry_radius", dimid_radius)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, "dry_radius", &
         bin_grid%n_bin, dimid_radius))
    call pmc_nc_check(nf90_def_dim(ncid, "dry_radius_edges", &
         bin_grid%n_bin + 1, dimid_radius_edges))
    call pmc_nc_check(nf90_def_var(ncid, "dry_radius", NF90_DOUBLE, &
         dimid_radius, varid_radius))
    call pmc_nc_check(nf90_put_att(ncid, varid_radius, "unit", "m"))
    call pmc_nc_check(nf90_def_var(ncid, "dry_radius_edges", NF90_DOUBLE, &
         dimid_radius_edges, varid_radius_edges))
    call pmc_nc_check(nf90_put_att(ncid, varid_radius_edges, "unit", "m"))
    call pmc_nc_check(nf90_def_var(ncid, "dry_radius_widths", NF90_DOUBLE, &
         dimid_radius, varid_radius_widths))
    call pmc_nc_check(nf90_put_att(ncid, varid_radius_widths, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

    do i_bin = 1,bin_grid%n_bin
       radius_centers(i_bin) = vol2rad(bin_grid%v(i_bin))
       radius_widths(i_bin) = bin_grid%dlnr
    end do
    do i_bin = 1,(bin_grid%n_bin + 1)
       radius_edges(i_bin) = vol2rad(bin_grid_edge(bin_grid, i_bin))
    end do
    call pmc_nc_check(nf90_put_var(ncid, varid_radius, radius_centers))
    call pmc_nc_check(nf90_put_var(ncid, varid_radius_edges, radius_edges))
    call pmc_nc_check(nf90_put_var(ncid, varid_radius_widths, radius_widths))

  end subroutine ensure_nc_dim_radius_dry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Write the aero species dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine ensure_nc_dim_aero_species(ncid, aero_data, dimid_aero_species)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_aero_species

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

  !> Write the gas species dimension to the given NetCDF file if it
  !> is not already present and in any case return the associated
  !> dimid.
  subroutine ensure_nc_dim_gas_species(ncid, gas_data, dimid_gas_species)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Gas_data structure.
    type(gas_data_t), intent(in) :: gas_data
    !> Dimid of the species dimension.
    integer, intent(out) :: dimid_gas_species

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
    call pmc_nc_check(nf90_put_att(ncid, varid_gas_species_widths, &
         "unit", "1"))

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

  !> Write the unit dimension to the given NetCDF file if it is not
  !> already present and in any case return the associated dimid.
  subroutine ensure_nc_dim_unit(ncid, dimid_unit)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the unit dimension.
    integer, intent(out) :: dimid_unit

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

  !> Write the env dimension to the given NetCDF file if it is not
  !> already present and in any case return the associated dimid.
  subroutine ensure_nc_dim_env(ncid, dimid_env)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Dimid of the env dimension.
    integer, intent(out) :: dimid_env

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

  !> Write the step dimension for the given process_spec to the given
  !> NetCDF file if it is not already present and in any case return
  !> the associated dimid.
  subroutine ensure_nc_dim_step(ncid, process_spec, dimid_step)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec
    !> Dimid of the step dimension.
    integer, intent(out) :: dimid_step

    integer :: status, i_step, varid_step
    integer :: dimid_step_edges, varid_step_edges, varid_step_widths
    real*8 :: step_centers(process_spec%n_step)
    real*8 :: step_edges(process_spec%n_step + 1)
    real*8 :: step_widths(process_spec%n_step)
    character(len=PROCESS_SPEC_NAME_LEN) :: dim_name_edges, &
         dim_name_widths, dim_unit

    if (process_spec%type == "kappa") then
       dim_unit = '1'
    elseif (process_spec%type == "comp") then
       dim_unit = '1'
    elseif (process_spec%type == "n_orig_part") then
       dim_unit = '1'
    elseif (process_spec%type == "optic_absorb") then
       dim_unit = 'm^2'
    elseif (process_spec%type == "optic_scatter") then
       dim_unit = 'm^2'
    elseif (process_spec%type == "optic_extinct") then
       dim_unit = 'm^2'
    else
       call die_msg(912387902, &
            "unknown process_spec%type: " // trim(process_spec%type))
    end if
    dim_name_edges = process_spec%dim_name
    dim_name_edges((len_trim(dim_name_edges)+1):) = "_edges"
    dim_name_widths = process_spec%dim_name
    dim_name_widths((len_trim(dim_name_widths)+1):) = "_widths"

    status = nf90_inq_dimid(ncid, process_spec%dim_name, dimid_step)
    if (status == NF90_NOERR) return
    if (status /= NF90_EBADDIM) call pmc_nc_check(status)

    ! dimension not defined, so define now define it
    call pmc_nc_check(nf90_redef(ncid))

    call pmc_nc_check(nf90_def_dim(ncid, process_spec%dim_name, &
         process_spec%n_step, dimid_step))
    call pmc_nc_check(nf90_def_dim(ncid, dim_name_edges, &
         process_spec%n_step + 1, dimid_step_edges))
    call pmc_nc_check(nf90_def_var(ncid, process_spec%dim_name, NF90_DOUBLE, &
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

  subroutine ensure_nc_var_env_state(ncid, process_spec, varid_env_state)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec
    !> Varid of env_state.
    integer, intent(out) :: varid_env_state

    integer :: dimid_time, dimid_env, dimids_env_state(2), status

    status = nf90_inq_varid(ncid, process_spec%name, varid_env_state)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    call ensure_nc_dim_env(ncid, dimid_env)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_env_state = (/ dimid_env, dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, process_spec%name, NF90_DOUBLE, &
         dimids_env_state, varid_env_state))
    call pmc_nc_check(nf90_put_att(ncid, varid_env_state, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_env_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_var_gas(ncid, process_spec, gas_data, varid_gas)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Varid of gas.
    integer, intent(out) :: varid_gas

    type(inout_file_t) :: file
    integer :: i_spec, status

    integer :: dimid_time, dimid_gas_species, dimids_gas(2)

    status = nf90_inq_varid(ncid, process_spec%name, varid_gas)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    call ensure_nc_dim_gas_species(ncid, gas_data, dimid_gas_species)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_gas = (/ dimid_gas_species, dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, process_spec%name, NF90_DOUBLE, &
         dimids_gas, varid_gas))
    call pmc_nc_check(nf90_put_att(ncid, varid_gas, "unit", "ppb"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_var_aero(ncid, process_spec, bin_grid, &
       aero_data, varid_aero)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec
    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Varid of aero.
    integer, intent(out) :: varid_aero

    integer :: dimid_time, dimid_radius, dimid_aero_species, dimid_unit
    integer :: status, dimids_aero(4)

    status = nf90_inq_varid(ncid, process_spec%name, varid_aero)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    if (process_spec%radius == "wet") then
       call ensure_nc_dim_radius(ncid, bin_grid, dimid_radius)
    elseif (process_spec%radius == "dry") then
       call ensure_nc_dim_radius_dry(ncid, bin_grid, dimid_radius)
    else
       call die_msg(484208239, &
            "radius must be 'wet' or 'dry', not: " &
            // trim(process_spec%radius))
    end if
    call ensure_nc_dim_aero_species(ncid, aero_data, dimid_aero_species)
    call ensure_nc_dim_unit(ncid, dimid_unit)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_aero = (/ dimid_radius, dimid_aero_species, dimid_unit, &
         dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, process_spec%name, NF90_DOUBLE, &
         dimids_aero, varid_aero))
    call pmc_nc_check(nf90_put_att(ncid, varid_aero, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ensure_nc_var_hist(ncid, process_spec, bin_grid, &
       aero_data, varid_hist)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec
    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Varid of hist.
    integer, intent(out) :: varid_hist

    integer :: dimid_time, dimid_step, dimid_radius, dimid_aero_species
    integer :: dimid_unit, status, dimids_hist(5)

    status = nf90_inq_varid(ncid, process_spec%name, varid_hist)
    if (status == NF90_NOERR) return
    if (status /= NF90_ENOTVAR) call pmc_nc_check(status)

    ! variable not defined, so define now define it
    call ensure_nc_dim_step(ncid, process_spec, dimid_step)
    if (process_spec%radius == "wet") then
       call ensure_nc_dim_radius(ncid, bin_grid, dimid_radius)
    elseif (process_spec%radius == "dry") then
       call ensure_nc_dim_radius_dry(ncid, bin_grid, dimid_radius)
    else
       call die_msg(894298398, &
            "radius must be 'wet' or 'dry', not: " &
            // trim(process_spec%radius))
    end if
    call ensure_nc_dim_aero_species(ncid, aero_data, dimid_aero_species)
    call ensure_nc_dim_unit(ncid, dimid_unit)
    call ensure_nc_dim_time(ncid, dimid_time)

    call pmc_nc_check(nf90_redef(ncid))

    dimids_hist = (/ dimid_step, dimid_radius, dimid_aero_species, &
         dimid_unit, dimid_time /)
    call pmc_nc_check(nf90_def_var(ncid, process_spec%name, NF90_DOUBLE, &
         dimids_hist, varid_hist))
    call pmc_nc_check(nf90_put_att(ncid, varid_hist, "unit", "1"))

    call pmc_nc_check(nf90_enddef(ncid))

  end subroutine ensure_nc_var_hist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_time(ncid, time, index, del_t)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current index.
    integer, intent(in) :: index
    !> Output timestep of current time.
    real*8, intent(in) :: del_t

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

  subroutine process_env(ncid, time, index, env_state, &
       process_spec)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current index.
    integer, intent(in) :: index
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec

    integer :: varid_env_state
    integer :: start(2), count(2)
    real*8 :: data(4)

    call ensure_nc_var_env_state(ncid, process_spec, varid_env_state)

    data(1) = env_state%temp
    data(2) = env_state%rel_humid
    data(3) = env_state%pressure
    data(4) = env_state%height

    start = (/ 1, index /)
    count = (/ 4, 1 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_env_state, data, &
         start = start, count = count))

  end subroutine process_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_gas(ncid, time, index, gas_data, gas_state, &
       process_spec)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current index.
    integer, intent(in) :: index
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec

    integer :: varid_gas
    integer :: start(2), count(2)

    call ensure_nc_var_gas(ncid, process_spec, gas_data, varid_gas)

    start = (/ 1, index /)
    count = (/ gas_data%n_spec, 1 /)
    call pmc_nc_check(nf90_put_var(ncid, varid_gas, gas_state%conc, &
         start = start, count = count))

  end subroutine process_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_aero(ncid, time, index, bin_grid, &
       aero_data, aero_binned, process_spec)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current index.
    integer, intent(in) :: index
    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero_binned structure.
    type(aero_binned_t), intent(in) :: aero_binned
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec

    real*8 :: aero(bin_grid%n_bin, aero_data%n_spec, 4)
    integer :: i_bin, i_spec, varid_aero
    integer :: start(4), count(4)

    call ensure_nc_var_aero(ncid, process_spec, bin_grid, &
         aero_data, varid_aero)

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

  subroutine process_aero(ncid, time, index, bin_grid, &
       aero_data, aero_state, process_spec)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current index.
    integer, intent(in) :: index
    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Aero_state structure.
    type(aero_state_t), intent(in) :: aero_state
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec

    type(aero_binned_t) :: aero_binned

    call aero_binned_alloc(aero_binned, bin_grid%n_bin, aero_data%n_spec)
    if (process_spec%radius == "wet") then
       call aero_state_to_binned(bin_grid, aero_data, aero_state, aero_binned)
    elseif (process_spec%radius == "dry") then
       call aero_state_to_binned_dry(bin_grid, aero_data, aero_state, &
            aero_binned)
    else
       call die(378223982)
    end if

    call output_aero(ncid, time, index, bin_grid, aero_data, &
         aero_binned, process_spec)

    call aero_binned_free(aero_binned)

  end subroutine process_aero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_hist(ncid, time, index, bin_grid, aero_data, &
       process_spec, hist)

    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current index.
    integer, intent(in) :: index
    !> Bin_grid structure.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aero_data structure.
    type(aero_data_t), intent(in) :: aero_data
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec
    !> Histogram data.
    real*8, intent(in) :: &
         hist(process_spec%n_step, bin_grid%n_bin, aero_data%n_spec, 4)

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

  !> Compute histogram by calling the step_comp() function on each
  !> particle.
  subroutine process_hist(ncid, time, index, bin_grid, &
       env_state, aero_data, aero_state, process_spec)
    
    !> NetCDF file ID, in data mode.
    integer, intent(in) :: ncid
    !> Current time (s).
    real*8, intent(in) :: time
    !> Current index.
    integer, intent(in) :: index
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> Process spec.
    type(process_spec_t), intent(in) :: process_spec

    real*8 :: hist(process_spec%n_step, bin_grid%n_bin, aero_data%n_spec, 4)
    integer :: i_step, i_bin, i_part, i, i_bin_hist
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
                  trim(process_spec%a_species(i))
             call die_msg(194029329, error_str)
          end if
       end do
       allocate(b_species(size(process_spec%b_species)))
       do i = 1,size(process_spec%b_species)
          b_species(i) = aero_data_spec_by_name(aero_data, &
               process_spec%b_species(i))
          if (b_species(i) == 0) then
             write(error_str, '(a,a)') 'unknown species: ', &
                  trim(process_spec%b_species(i))
             call die_msg(283298183, error_str)
          end if
       end do
    end if

    hist = 0d0
    scale = 1d0 / bin_grid%dlnr / step_width / aero_state%comp_vol
    do i_bin = 1,bin_grid%n_bin
       do i_part = 1,aero_state%bin(i_bin)%n_part
          aero_particle => aero_state%bin(i_bin)%particle(i_part)

          if (process_spec%type == "kappa") then
             rh = aero_particle_kappa_rh(aero_particle, aero_data, env_state)
             supersat = rh - 1d0
             val = supersat
          elseif (process_spec%type == "comp") then
             val = aero_particle_comp(aero_particle, aero_data, &
                  process_spec%comp_type, a_species, b_species)
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

          if (process_spec%radius == "wet") then
             i_bin_hist = i_bin
          elseif (process_spec%radius == "dry") then
             i_bin_hist = bin_grid_particle_in_bin(bin_grid, &
                  aero_particle_solute_volume(aero_particle, aero_data))
          else
             call die_msg(842982739, &
                  "radius must be 'wet' or 'dry', not: " &
                  // trim(process_spec%radius))
          end if

          hist(i_step, i_bin_hist, :, 1) = hist(i_step, i_bin_hist, :, 1) &
               + 1d0 / dble(aero_data%n_spec) * scale
          hist(i_step, i_bin_hist, :, 2) = hist(i_step, i_bin_hist, :, 2) &
               + aero_particle%vol * scale
          hist(i_step, i_bin_hist, :, 3) = hist(i_step, i_bin_hist, :, 3) &
               + aero_particle%vol * aero_data%density * scale
          hist(i_step, i_bin_hist, :, 4) = hist(i_step, i_bin_hist, :, 4) &
               + aero_particle%vol * aero_data%density &
               / aero_data%molec_weight * scale
       end do
    end do

    if (process_spec%type == "comp") then
       deallocate(a_species, b_species)
    end if

    call output_hist(ncid, time, index, bin_grid, aero_data, &
         process_spec, hist)

  end subroutine process_hist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real*8 function aero_particle_comp(aero_particle, aero_data, &
       comp_type, a_species, b_species)

    !> Particle.
    type(aero_particle_t), intent(in) :: aero_particle
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Composition type ("volume", "mass", or "mole").
    character(len=PROCESS_SPEC_TYPE_LEN) :: comp_type
    !> First list of species.
    integer, intent(in) :: a_species(:)
    !> Second list of species.
    integer, intent(in) :: b_species(:)

    integer :: i
    real*8 :: a_val, b_val

    a_val = 0d0
    do i = 1,size(a_species)
       if (comp_type == "volume") then
          a_val = a_val + aero_particle%vol(a_species(i))
       elseif (comp_type == "mass") then
          a_val = a_val + aero_particle%vol(a_species(i)) &
               * aero_data%density(a_species(i))
       elseif (comp_type == "mole") then
          a_val = a_val + aero_particle%vol(a_species(i)) &
               * aero_data%density(a_species(i)) &
               / aero_data%molec_weight(a_species(i))
       else
          call die_msg(298305208, "unknown comp_type: " // trim(comp_type))
       end if
    end do
    b_val = 0d0
    do i = 1,size(b_species)
       if (comp_type == "volume") then
          b_val = b_val + aero_particle%vol(b_species(i))
       elseif (comp_type == "mass") then
          b_val = b_val + aero_particle%vol(b_species(i)) &
               * aero_data%density(b_species(i))
       elseif (comp_type == "mole") then
          b_val = b_val + aero_particle%vol(b_species(i)) &
               * aero_data%density(b_species(i)) &
               / aero_data%molec_weight(b_species(i))
       else
          call die_msg(842198113, "unknown comp_type: " // trim(comp_type))
       end if
    end do
    call assert(880038232, a_val >= 0d0)
    call assert(715496111, b_val >= 0d0)
    if ((a_val == 0d0) .and. (b_val == 0d0)) then
       aero_particle_comp = 0d0
    else
       aero_particle_comp = b_val / (a_val + b_val)
    end if

  end function aero_particle_comp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module pmc_output_processed
