! Copyright (C) 2007 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! Interface to MOSAIC

module mod_mosaic
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine singlestep_mosaic(M, n_spec, n_bin, MH, VH, &
       bin_v, bin_g, bin_gs, bin_n, dlnr, t, del_t, env, mat)

    use mod_constants
    use mod_util
    use mod_array
    use mod_bin 
    use mod_condensation
    use mod_environ
    use mod_material
    use mod_state
    use mod_gas
    use module_data_mosaic_kind ! NO_AUTO_DEP
    use module_data_mosaic_aero ! NO_AUTO_DEP
    use module_data_mosaic_gas ! NO_AUTO_DEP
    use module_data_mosaic_main ! NO_AUTO_DEP
    
    integer, intent(inout) :: M         ! actual number of particles
    integer, intent(in) :: n_spec       ! number of species
    integer, intent(in) :: n_bin        ! number of bins
    integer, intent(inout) :: MH(n_bin) ! number of particles per bin
    type(bin_p), intent(inout) :: VH(n_bin) ! particle volumes (m^3)
    
    real*8, intent(in) :: bin_v(n_bin)  ! volume of particles in bins (m^3)
    real*8, intent(out) :: bin_g(n_bin) ! volume in bins  
    real*8, intent(out) :: bin_gs(n_bin,n_spec) ! species volume in bins
    integer, intent(out) :: bin_n(n_bin) ! number in bins
    real*8, intent(in) :: dlnr          ! bin scale factor
    
    real*8, intent(in) :: t             ! current time (s)
    real*8, intent(in) :: del_t         ! timestep for coagulation
    
    type(environ), intent(inout) :: env ! environment state
    type(material), intent(in) :: mat   ! material properties

    ! MOSAIC function interfaces
    interface
       subroutine UpdateTime()
       end subroutine UpdateTime
       subroutine SolarZenithAngle()
       end subroutine SolarZenithAngle
       subroutine IntegrateChemistry()
       end subroutine IntegrateChemistry
       subroutine DoMassBalance()
       end subroutine DoMassBalance
    end interface

    ! local variables
    real*8 :: t_utc ! time (s since 00:00 UTC)
    integer :: i_bin, i_num, i_spec, i_mosaic

    t_utc = env%start_time + t

    ! set MOSAIC data

    ! beginning date
    tbeg_mo = 0 ! FIXME: determine these from env%start_day
    tbeg_dd = 0
    tbeg_hh = 0

    ! beginning time
    tbeg_hh = 0 ! FIXME: determine these from env%start_time
    tbeg_mm = 0
    tbeg_ss = 0

    ! run time
    trun_dd = 0 ! hopefully we don't need to set these
    trun_hh = 0
    trun_mm = 0
    trun_ss = 0

    ! transport timestep (min)
    dt_min = 0d0
    ! aerosol optics timestep (min)
    dt_aeroptic_min = 0d0

    ! geographic location
    rlon = env%longitude
    rlat = env%latitude
    zalt_m = env%altitude

    ! environmental parameters
    RH = env%RH
    te = env%T
    pr_atm = env%p / const%atm

    ! aerosol data
    nbin_a = M
    i_mosaic = 0
    do i_bin = 1,n_bin
       do i_num = 1,MH(i_bin)
          do i_spec = 1,n_spec
             i_mosaic = i_mosaic + 1
             aer(i_spec, jtotal, i_mosaic) = VH(i_bin)%p(i_num, i_spec)
          end do
       end do
    end do

    ! parameters
    mmode = 1               ! 1 = time integration, 2 = parametric analysis
    mgas = 1                ! 1 = gas chem on, 0 = gas chem off
    maer = 1                ! 1 = aer chem on, 0 = aer chem off
    mcld = 0                ! 1 = cld chem on, 0 = cld chem off
    maeroptic = 1           ! 1 = aer_optical on, 0 = aer_optical off
    mshellcore = 1          ! 0 = no shellcore, 1 = core is BC only
                            ! 2 = core is BC and DUST
    msolar = 1              ! 1 = diurnally varying phot, 2 = fixed phot
    mphoto = 2              ! 1 = Rick's param, 2 = Yang's param
    mGAS_AER_XFER = 1       ! 1 = gas-aerosol partitioning, 0 = do not partition
    mDYNAMIC_SOLVER = 1     ! 1 = astem, 2 = lsodes
    alpha_ASTEM = 0.5d0     ! solver parameter. range: 0.01 - 1.0
    rtol_eqb_ASTEM = 0.01d0 ! relative eqb tolerance. range: 0.01 - 0.03
    ptol_mol_ASTEM = 0.01d0 ! percent mol tolerance.  range: 0.01 - 1.0

    call UpdateTime
    if(msolar.eq.1)then
       call SolarZenithAngle
    end if
    call IntegrateChemistry
    call DoMassBalance

  end subroutine singlestep_mosaic

end module mod_mosaic
