! Copyright (C) 2005-2007 Nicole Riemer and Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

module mod_output

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_open(output_unit, output_file, n_loop, n_bin, &
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
    
  end subroutine output_open
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_info(output_unit, time, n_bin, n_spec, &
       bin_v, bin_g, bin_gs, bin_n, dlnr, env, aero_data, i_loop)

    ! Write the current binned data to the output file. This version
    ! of the function takes absolute number and absolute volume
    ! per-bin (as produced by a particle-resolved code, for example).
    
    use mod_aero_data
    use mod_environ
    
    integer, intent(in) :: output_unit  ! unit number to output to
    real*8, intent(in) :: time          ! simulation time
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(in) :: bin_g(n_bin)  ! volume in bins (m^3)
    real*8, intent(in) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(in) :: bin_n(n_bin) ! number in bins (dimensionless)
    real*8, intent(in) :: dlnr          ! bin scale factor
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    integer, intent(in) :: i_loop       ! current loop number
    
    real*8 bin_g_den(n_bin), bin_gs_den(n_bin,n_spec)
    real*8 bin_n_den(n_bin)
    
    bin_g_den = bin_g / env%V_comp / dlnr
    bin_gs_den = bin_gs / env%V_comp / dlnr
    bin_n_den = dble(bin_n) / env%V_comp / dlnr
    call output_info_density(output_unit, time, n_bin, n_spec, bin_v, &
         bin_g_den, bin_gs_den, bin_n_den, env, aero_data, i_loop)
    
  end subroutine output_info
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine output_info_density(output_unit, time, n_bin, n_spec, bin_v, &
       bin_g_den, bin_gs_den, bin_n_den, env, aero_data, i_loop)

    ! Write the current binned data to the output file. This version
    ! of the function takes number and volume densities (as produced
    ! by a sectional code, for example).
    
    use mod_aero_data
    use mod_environ
    use mod_util
    
    integer, intent(in) :: output_unit  ! unit number to output to
    real*8, intent(in) :: time          ! simulation time
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(in) :: n_spec       ! number of species
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(in) :: bin_g_den(n_bin) ! volume density in bins (1)
    real*8, intent(in) :: bin_gs_den(n_bin,n_spec) ! spec vol den in bins (1)
    real*8, intent(in) :: bin_n_den(n_bin) ! number density in bins (1/m^3)
    type(environ), intent(in) :: env    ! environment state
    type(aero_data_t), intent(in) :: aero_data   ! aerosol data
    integer, intent(in) :: i_loop       ! current loop number
    
    integer k

    write(output_unit,'(a10,i20)') 'loop_num', i_loop
    write(output_unit,'(a10,e20.10)') 'time(s)', time
    write(output_unit,'(a10,e20.10)') 'temp(K)', env%T
    write(output_unit,'(a10,e20.10)') 'RH(1)', env%RH
    write(output_unit,'(a1,a9,a20,a20,a20,a30)') '#', 'bin_num', &
         'radius(m)', 'tot num (#/m^3)', 'tot vol (m^3/m^3)', &
         'vol per species (m^3/m^3)'
    do k = 1,n_bin
       write(output_unit, '(i10,23e20.10)') k, vol2rad(bin_v(k)), &
            bin_n_den(k), bin_g_den(k), bin_gs_den(k,:)
    end do
    
  end subroutine output_info_density
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_output
