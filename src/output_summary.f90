! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

module mod_output_summary

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_summary_open(output_unit, output_file, n_loop, n_bin, &
       n_spec, n_time)

    ! Open an output file for writing.

    integer, intent(in) :: output_unit  ! unit number to output to
    character(len=300) :: output_file   ! output filename
    integer, intent(in) :: n_loop       ! number of loops
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: n_time       ! number of times

    open(unit=output_unit, file=output_file)

    write(output_unit,'(a10,i10)') 'n_loop', n_loop
    write(output_unit,'(a10,i10)') 'n_bin', n_bin
    write(output_unit,'(a10,i10)') 'n_time', n_time
    write(output_unit,'(a10,i10)') 'n_spec', n_spec
    
  end subroutine output_summary_open
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_summary(output_unit, time, bin_grid, aero_data, &
       bin_v, bin_vs, bin_n, env, i_loop, comp_vol)

    ! Write the current binned data to the output file. This version
    ! of the function takes absolute number and absolute volume
    ! per-bin (as produced by a particle-resolved code, for example).
    
    use mod_bin
    use mod_aero_data
    use mod_environ
    
    integer, intent(in) :: output_unit  ! unit number to output to
    real*8, intent(in) :: time          ! simulation time
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    real*8, intent(in) :: bin_v(bin_grid%n_bin) ! volume in bins (m^3)
    real*8, intent(in) :: bin_vs(bin_grid%n_bin,aero_data%n_spec) ! volume per
                                        ! bin and species (m^3)
    integer, intent(in) :: bin_n(bin_grid%n_bin) ! number in bins (dim-less)
    type(environ), intent(in) :: env    ! environment state
    integer, intent(in) :: i_loop       ! current loop number
    real*8, intent(in) :: comp_vol      ! FIXME: temporary hack until
                                        ! aero_dist stores densities
    
    real*8 :: bin_v_den(bin_grid%n_bin), bin_n_den(bin_grid%n_bin)
    real*8 :: bin_vs_den(bin_grid%n_bin,aero_data%n_spec)
    
    bin_v_den = bin_v / comp_vol / bin_grid%dlnr
    bin_vs_den = bin_vs / comp_vol / bin_grid%dlnr
    bin_n_den = dble(bin_n) / comp_vol / bin_grid%dlnr
    call output_summary_density(output_unit, time, bin_grid, aero_data, &
         bin_v_den, bin_vs_den, bin_n_den, env, i_loop)
    
  end subroutine output_summary
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_summary_density(output_unit, time, bin_grid, &
       aero_data, bin_v_den, bin_vs_den, bin_n_den, env, i_loop)

    ! Write the current binned data to the output file. This version
    ! of the function takes number and volume densities (as produced
    ! by a sectional code, for example).

    use mod_bin
    use mod_aero_data
    use mod_environ
    use mod_util
    
    integer, intent(in) :: output_unit  ! unit number to output to
    real*8, intent(in) :: time          ! simulation time
    type(bin_grid_t), intent(in) :: bin_grid ! bin grid
    type(aero_data_t), intent(in) :: aero_data ! aerosol data
    real*8, intent(in) :: bin_v_den(bin_grid%n_bin) ! volume density in bins (1)
    real*8, intent(in) :: bin_vs_den(bin_grid%n_bin,aero_data%n_spec) ! volume
                                        ! density per bin and per species (1)
    real*8, intent(in) :: bin_n_den(bin_grid%n_bin) ! num dens in bins (1/m^3)
    type(environ), intent(in) :: env    ! environment state
    integer, intent(in) :: i_loop       ! current loop number
    
    integer k

    write(output_unit,'(a10,i20)') 'loop_num', i_loop
    write(output_unit,'(a10,e20.10)') 'time(s)', time
    write(output_unit,'(a10,e20.10)') 'temp(K)', env%T
    write(output_unit,'(a10,e20.10)') 'RH(1)', env%RH
    write(output_unit,'(a1,a9,a20,a20,a20,a30)') '#', 'bin_num', &
         'radius(m)', 'tot num (#/m^3)', 'tot vol (m^3/m^3)', &
         'vol per species (m^3/m^3)'
    do k = 1,bin_grid%n_bin
       write(output_unit, '(i10,23e20.10)') k, vol2rad(bin_grid%v(k)), &
            bin_n_den(k), bin_v_den(k), bin_vs_den(k,:)
    end do
    
  end subroutine output_summary_density
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_output_summary
