! Copyright (C) 2005-2009 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_aero_mode module.

!> The aero_mode_t structure and associated subroutines.
module pmc_aero_mode

  use pmc_bin_grid
  use pmc_util
  use pmc_constants
  use pmc_spec_read
  use pmc_aero_data
  use pmc_mpi
  use pmc_rand
#ifdef PMC_USE_MPI
  use mpi
#endif

  !> Maximum length of an aero_dist mode name.
  integer, parameter :: AERO_MODE_NAME_LEN = 300
  !> Maximum length of an aero_dist mode type.
  integer, parameter :: AERO_MODE_TYPE_LEN = 100

  !> An aerosol size distribution mode.
  !!
  !! Each mode is assumed to be fully internally mixed so that every
  !! particle has the same composition. The \c num_conc array then
  !! stores the number concentration distribution.
  type aero_mode_t
     !> Mode name, used to track particle sources.
     character(len=AERO_MODE_NAME_LEN) :: name
     !> Mode type ("log_normal", "exp", or "mono").
     character(len=AERO_MODE_TYPE_LEN) :: type
     !> Mean radius of mode (m).
     real*8 :: mean_radius
     !> Log base 10 of geometric standard deviation of radius, if necessary (m).
     real*8 :: log10_std_dev_radius
     !> Total number concentration of mode (#/m^3).
     real*8 :: num_conc
     !> Species fractions by volume [length \c aero_data%%n_spec] (1).
     real*8, pointer :: vol_frac(:)
  end type aero_mode_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_mode.
  subroutine aero_mode_allocate(aero_mode)

    !> Aerosol mode.
    type(aero_mode_t), intent(out) :: aero_mode

    allocate(aero_mode%vol_frac(0))

  end subroutine aero_mode_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Allocates an aero_mode of the given size.
  subroutine aero_mode_allocate_size(aero_mode, n_spec)

    !> Aerosol mode.
    type(aero_mode_t), intent(out) :: aero_mode
    !> Number of species.
    integer, intent(in) :: n_spec

    allocate(aero_mode%vol_frac(n_spec))

  end subroutine aero_mode_allocate_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Free all storage.
  subroutine aero_mode_deallocate(aero_mode)

    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode

    deallocate(aero_mode%vol_frac)

  end subroutine aero_mode_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Copy an aero_mode.
  subroutine aero_mode_copy(aero_mode_from, aero_mode_to)

    !> Aerosol mode original.
    type(aero_mode_t), intent(in) :: aero_mode_from
    !> Aerosol mode copy.
    type(aero_mode_t), intent(inout) :: aero_mode_to

    call aero_mode_deallocate(aero_mode_to)
    call aero_mode_allocate_size(aero_mode_to, size(aero_mode_from%vol_frac))
    aero_mode_to%name = aero_mode_from%name
    aero_mode_to%type = aero_mode_from%type
    aero_mode_to%mean_radius = aero_mode_from%mean_radius
    aero_mode_to%log10_std_dev_radius = aero_mode_from%log10_std_dev_radius
    aero_mode_to%num_conc = aero_mode_from%num_conc
    aero_mode_to%vol_frac = aero_mode_from%vol_frac

  end subroutine aero_mode_copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution, normalized so that
  !> sum(num_conc(k) * dlnr) = 1
  subroutine num_conc_log_normal(mean_radius, log_sigma, bin_grid, num_conc)
    
    !> Geometric mean radius (m).
    real*8, intent(in) :: mean_radius
    !> log_10(geom. std dev) (1).
    real*8, intent(in) :: log_sigma
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real*8, intent(out) :: num_conc(bin_grid%n_bin)
    
    integer :: k
    
    do k = 1,bin_grid%n_bin
       num_conc(k) = 1d0 / (sqrt(2d0 * const%pi) * log_sigma) * &
            dexp(-(dlog10(vol2rad(bin_grid%v(k))) &
            - dlog10(mean_radius))**2d0 &
            / (2d0 * log_sigma**2d0)) / dlog(10d0)
    end do
    
    ! The formula above was originally for a distribution in
    ! log_10(r), while we are using log_e(r) for our bin grid. The
    ! division by dlog(10) at the end corrects for this.

    ! Remember that log_e(r) = log_10(r) * log_e(10).
    
  end subroutine num_conc_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Compute a log-normal distribution in volume.
  subroutine vol_conc_log_normal(mean_radius, log_sigma, bin_grid, vol_conc)
    
    !> Geometric mean radius (m).
    real*8, intent(in) :: mean_radius
    !> log_10(geom. std dev) (1).
    real*8, intent(in) :: log_sigma
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume concentration (V(ln(r))d(ln(r))).
    real*8, intent(out) :: vol_conc(bin_grid%n_bin)
    
    real*8 :: num_conc(bin_grid%n_bin)

    call num_conc_log_normal(mean_radius, log_sigma, bin_grid, num_conc)
    vol_conc = num_conc * bin_grid%v
    
  end subroutine vol_conc_log_normal
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume
  !> n(v) = 1 / mean_vol * exp(- v / mean_vol)
  !> Normalized so that sum(num_conc(k) * dlnr) = 1
  subroutine num_conc_exp(mean_radius, bin_grid, num_conc)
    
    !> Mean radius (m).
    real*8, intent(in) :: mean_radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real*8, intent(out) :: num_conc(bin_grid%n_bin)
    
    integer :: k
    real*8 :: mean_vol, num_conc_vol
    
    mean_vol = rad2vol(mean_radius)
    do k = 1,bin_grid%n_bin
       num_conc_vol = 1d0 / mean_vol * exp(-(bin_grid%v(k) / mean_vol))
       call vol_to_lnr(vol2rad(bin_grid%v(k)), num_conc_vol, num_conc(k))
    end do
    
  end subroutine num_conc_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Exponential distribution in volume.
  subroutine vol_conc_exp(mean_radius, bin_grid, vol_conc)
    
    !> Mean radius (m).
    real*8, intent(in) :: mean_radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume concentration (V(ln(r))d(ln(r))).
    real*8, intent(out) :: vol_conc(bin_grid%n_bin)
    
    real*8 :: num_conc(bin_grid%n_bin)

    call num_conc_exp(mean_radius, bin_grid, num_conc)
    vol_conc = num_conc * bin_grid%v
    
  end subroutine vol_conc_exp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution.
  !> Normalized so that sum(num_conc(k) * dlnr) = 1
  subroutine num_conc_mono(radius, bin_grid, num_conc)
    
    !> Radius of each particle (m^3).
    real*8, intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Number concentration (#(ln(r))d(ln(r))).
    real*8, intent(out) :: num_conc(bin_grid%n_bin)
    
    integer :: k

    num_conc = 0d0
    k = bin_grid_particle_in_bin(bin_grid, rad2vol(radius))
    num_conc(k) = 1d0 / bin_grid%dlnr
    
  end subroutine num_conc_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Mono-disperse distribution in volume.
  subroutine vol_conc_mono(radius, bin_grid, vol_conc)
    
    !> Radius of each particle (m^3).
    real*8, intent(in) :: radius
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Volume concentration (V(ln(r))d(ln(r))).
    real*8, intent(out) :: vol_conc(bin_grid%n_bin)
    
    integer :: k

    vol_conc = 0d0
    k = bin_grid_particle_in_bin(bin_grid, rad2vol(radius))
    vol_conc(k) = 1d0 / bin_grid%dlnr * rad2vol(radius)
    
  end subroutine vol_conc_mono
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the binned number concentration for an aero_mode.
  subroutine aero_mode_num_conc(aero_mode, bin_grid, aero_data, &
       num_conc)

    !> Aero mode for which to compute number concentration.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Number concentration (#(ln(r))d(ln(r))).
    real*8, intent(out) :: num_conc(bin_grid%n_bin)

    if (aero_mode%type == "log_normal") then
       call num_conc_log_normal(aero_mode%mean_radius, &
            aero_mode%log10_std_dev_radius, bin_grid, num_conc)
    elseif (aero_mode%type == "exp") then
       call num_conc_exp(aero_mode%mean_radius, bin_grid, num_conc)
    elseif (aero_mode%type == "mono") then
       call num_conc_mono(aero_mode%mean_radius, bin_grid, num_conc)
    else
       call die_msg(719625922, "Unknown aero_mode type")
    end if
    num_conc = num_conc * aero_mode%num_conc

  end subroutine aero_mode_num_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return the binned per-species volume concentration for an
  !> aero_mode.
  subroutine aero_mode_vol_conc(aero_mode, bin_grid, aero_data, &
       vol_conc)

    !> Aero mode for which to compute volume concentration.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Bin grid.
    type(bin_grid_t), intent(in) :: bin_grid
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Volume concentration (V(ln(r))d(ln(r))).
    real*8, intent(out) :: vol_conc(bin_grid%n_bin, aero_data%n_spec)

    integer :: i_spec
    real*8 :: vol_conc_total(bin_grid%n_bin)

    if (aero_mode%type == "log_normal") then
       call vol_conc_log_normal(aero_mode%mean_radius, &
            aero_mode%log10_std_dev_radius, bin_grid, vol_conc_total)
    elseif (aero_mode%type == "exp") then
       call vol_conc_exp(aero_mode%mean_radius, bin_grid, &
            vol_conc_total)
    elseif (aero_mode%type == "mono") then
       call vol_conc_mono(aero_mode%mean_radius, bin_grid, &
            vol_conc_total)
    else
       call die_msg(314169653, "Unknown aero_mode type")
    end if
    vol_conc_total = vol_conc_total * aero_mode%num_conc
    do i_spec = 1,aero_data%n_spec
       vol_conc(:,i_spec) = vol_conc_total &
            * aero_mode%vol_frac(i_spec)
    end do

  end subroutine aero_mode_vol_conc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Return a radius randomly sampled from the mode distribution.
  subroutine aero_mode_sample_radius(aero_mode, radius)

    !> Aero_mode to sample radius from.
    type(aero_mode_t), intent(in) :: aero_mode
    !> Sampled radius (m).
    real*8, intent(out) :: radius

    if (aero_mode%type == "log_normal") then
       radius = 10d0**rand_normal(log10(aero_mode%mean_radius), &
            aero_mode%log10_std_dev_radius)
    elseif (aero_mode%type == "exp") then
       radius = vol2rad(- rad2vol(aero_mode%mean_radius) * log(pmc_random()))
    elseif (aero_mode%type == "mono") then
       radius = aero_mode%mean_radius
    else
       call die_msg(749122931, "Unknown aero_mode type")
    end if

  end subroutine aero_mode_sample_radius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read volume fractions from a data file.
  subroutine spec_read_vol_frac(file, aero_data, vol_frac)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol species volume fractions.
    real*8, intent(out) :: vol_frac(:)

    integer :: n_species, species, i
    character(len=MAX_VAR_LEN) :: read_name
    type(spec_file_t) :: read_file
    character(len=MAX_VAR_LEN), pointer :: species_name(:)
    real*8, pointer :: species_data(:,:)
    real*8 :: tot_vol_frac

    ! read the aerosol data from the specified file
    call spec_read_string(file, 'mass_frac', read_name)
    call spec_read_open(read_name, read_file)
    allocate(species_name(0))
    allocate(species_data(0,0))
    call spec_read_real_named_array(read_file, 0, species_name, species_data)
    call spec_read_close(read_file)

    ! check the data size
    n_species = size(species_data, 1)
    if (n_species < 1) then
       call die_msg(427666881, 'file ' // trim(read_name) &
            // ' must contain at least one line of data')
    end if
    if (size(species_data, 2) /= 1) then
       call die_msg(427666881, 'each line in file ' &
            // trim(read_name) // ' must contain exactly one data value')
    end if

    ! copy over the data
    vol_frac = 0d0
    do i = 1,n_species
       species = aero_data_spec_by_name(aero_data, species_name(i))
       if (species == 0) then
          call die_msg(775942501, 'unknown species ' &
               // trim(species_name(i)) // ' in file ' &
               // trim(read_name))
       end if
       vol_frac(species) = species_data(i,1)
    end do
    deallocate(species_name)
    deallocate(species_data)

    ! convert mass fractions to volume fractions
    vol_frac = vol_frac / aero_data%density
    
    ! normalize
    tot_vol_frac = sum(vol_frac)
    if ((minval(vol_frac) < 0d0) .or. (tot_vol_frac <= 0d0)) then
       call die_msg(356648030, 'vol_frac in ' // trim(read_name) &
            // ' is not positive')
    end if
    vol_frac = vol_frac / tot_vol_frac

  end subroutine spec_read_vol_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read the shape (number concentration profile) of one mode of an aerosol
  !> distribution.
  subroutine spec_read_aero_mode_shape(file, aero_mode)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode

    character(len=MAX_VAR_LEN) :: mode_type

    call spec_read_real(file, 'num_conc', aero_mode%num_conc)
    call spec_read_string(file, 'mode_type', mode_type)
    if (len_trim(mode_type) < AERO_MODE_TYPE_LEN) then
       aero_mode%type = mode_type(1:AERO_MODE_TYPE_LEN)
    else
       call spec_read_die_msg(284789262, file, "mode_type string too long")
    end if
    if (trim(mode_type) == 'log_normal') then
       call spec_read_real(file, 'mean_radius', aero_mode%mean_radius)
       call spec_read_real(file, 'log_std_dev', aero_mode%log10_std_dev_radius)
    elseif (trim(mode_type) == 'exp') then
       call spec_read_real(file, 'mean_radius', aero_mode%mean_radius)
    elseif (trim(mode_type) == 'mono') then
       call spec_read_real(file, 'radius', aero_mode%mean_radius)
    else
       call spec_read_die_msg(729472928, file, "Unknown distribution type")
    end if

  end subroutine spec_read_aero_mode_shape

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Read one mode of an aerosol distribution (number concentration and
  !> volume fractions).
  subroutine spec_read_aero_mode(file, aero_data, aero_mode, eof)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aero_data data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol mode.
    type(aero_mode_t), intent(inout) :: aero_mode
    !> If eof instead of reading data.
    logical :: eof

    character(len=MAX_VAR_LEN) :: tmp_str
    type(spec_line_t) :: line

    call spec_line_allocate(line)
    call spec_read_line(file, line, eof)
    if (.not. eof) then
       call spec_read_check_line_name(file, line, "mode_name")
       call spec_read_check_line_length(file, line, 1)
       call aero_mode_deallocate(aero_mode)
       call aero_mode_allocate_size(aero_mode, aero_data%n_spec)
       tmp_str = line%data(1) ! hack to avoid gfortran warning
       aero_mode%name = tmp_str(1:AERO_MODE_NAME_LEN)
       allocate(aero_mode%vol_frac(aero_data%n_spec))
       call spec_read_vol_frac(file, aero_data, aero_mode%vol_frac)
       call spec_read_aero_mode_shape(file, aero_mode)
    end if
    call spec_line_deallocate(line)

  end subroutine spec_read_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Determines the number of bytes required to pack the given value.
  integer function pmc_mpi_pack_size_aero_mode(val)

    !> Value to pack.
    type(aero_mode_t), intent(in) :: val

    pmc_mpi_pack_size_aero_mode = &
         pmc_mpi_pack_size_string(val%name) &
         + pmc_mpi_pack_size_string(val%type) &
         + pmc_mpi_pack_size_real(val%mean_radius) &
         + pmc_mpi_pack_size_real(val%log10_std_dev_radius) &
         + pmc_mpi_pack_size_real(val%num_conc) &
         + pmc_mpi_pack_size_real_array(val%vol_frac)

  end function pmc_mpi_pack_size_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Packs the given value into the buffer, advancing position.
  subroutine pmc_mpi_pack_aero_mode(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_mode_t), intent(in) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_pack_string(buffer, position, val%name)
    call pmc_mpi_pack_string(buffer, position, val%type)
    call pmc_mpi_pack_real(buffer, position, val%mean_radius)
    call pmc_mpi_pack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_pack_real(buffer, position, val%num_conc)
    call pmc_mpi_pack_real_array(buffer, position, val%vol_frac)
    call assert(579699255, &
         position - prev_position == pmc_mpi_pack_size_aero_mode(val))
#endif

  end subroutine pmc_mpi_pack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Unpacks the given value from the buffer, advancing position.
  subroutine pmc_mpi_unpack_aero_mode(buffer, position, val)

    !> Memory buffer.
    character, intent(inout) :: buffer(:)
    !> Current buffer position.
    integer, intent(inout) :: position
    !> Value to pack.
    type(aero_mode_t), intent(out) :: val

#ifdef PMC_USE_MPI
    integer :: prev_position

    prev_position = position
    call pmc_mpi_unpack_string(buffer, position, val%name)
    call pmc_mpi_unpack_string(buffer, position, val%type)
    call pmc_mpi_unpack_real(buffer, position, val%mean_radius)
    call pmc_mpi_unpack_real(buffer, position, val%log10_std_dev_radius)
    call pmc_mpi_unpack_real(buffer, position, val%num_conc)
    call pmc_mpi_unpack_real_array(buffer, position, val%vol_frac)
    call assert(874467577, &
         position - prev_position == pmc_mpi_pack_size_aero_mode(val))
#endif

  end subroutine pmc_mpi_unpack_aero_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_aero_mode
