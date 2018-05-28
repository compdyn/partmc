include 'f90/TRACERS_AEROSOLS_SOA.F90'

PROGRAM KPP_ROOT_Driver

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize
  use TRACERS_SOA, only: soa_bins, soa_cond, soa_tr, &
                         soa_init, soa_calc, &
                         isoa_apin, isoa_isop, isoa_oh, isoa_o3, isoa_no3, &
                         isoa_isopp1, isoa_isopp2, isoa_apinp1, isoa_apinp2

      KPP_REAL :: T, DVAL(NSPEC)
      KPP_REAL :: RSTATE(20)
      INTEGER :: i, v
! EQSAM
      real, parameter    :: eqsam_conv = 1.e12 / avog ! unit conversion factor
                                                      ! from molecules cm-3
                                                      ! to umol m-3
      integer, parameter :: nca=11     ! number of input variables
      integer, parameter :: nco=36     ! number of output variables
      integer            :: iopt=1     ! aerosol state and history
                                       ! 1 selects metastable (wet)
                                       ! 2 selects solid (dry)
      integer, parameter :: loop=1     ! only one timestep
      integer, parameter :: imax=1     ! only one timestep
      integer, parameter :: ipunit=99  ! diagnostic output printing (not used)
      real, dimension(imax, nca) :: yi ! input array
      real, dimension(imax, nco) :: yo ! output array
! SOA
      real,parameter :: molec2ug=1.d12/avog ! multiply with MW later!!
      real :: voc2nox_denom
  
!~~~> Initialization 

      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO
     
      CALL Initialize()

! SOA
      soa_cond%dt=dt
      soa_tr%ivoc(isoa_apin)=ind_Terpenes
      soa_tr%ivoc(isoa_isop)=ind_Isoprene
      soa_tr%ioxidant(isoa_oh)=ind_OH
      soa_tr%ioxidant(isoa_o3)=ind_O3
      soa_tr%ioxidant(isoa_no3)=ind_NO3
      soa_tr%igas(isoa_isopp1)=ind_isopp1g
      soa_tr%igas(isoa_isopp2)=ind_isopp2g
      soa_tr%igas(isoa_apinp1)=ind_apinp1g
      soa_tr%igas(isoa_apinp2)=ind_apinp2g
      soa_tr%iaer(isoa_isopp1)=ind_isopp1a
      soa_tr%iaer(isoa_isopp2)=ind_isopp2a
      soa_tr%iaer(isoa_apinp1)=ind_apinp1a
      soa_tr%iaer(isoa_apinp2)=ind_apinp2a
      call soa_init(nvar)

      CALL InitSaveData()

!~~~> Time loop
      T = TSTART
kron: DO WHILE (T < TEND)

        TIME = T
        CALL GetMass( C, DVAL )
        WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,       &
                   ( TRIM(SPC_NAMES(MONITOR(i))),           &
                     C(MONITOR(i))/CFACTOR, i=1,NMONITOR )
        CALL SaveData()
        CALL Update_SUN() 
        CALL Update_RCONST()

! SOA, before chemistry. VOCs in ug m-3, oxidants in molecules cm-3
        soa_tr%voc(isoa_apin)=c(ind_Terpenes)*molec2ug*136.d0
        soa_tr%voc(isoa_isop)=c(ind_Isoprene)*molec2ug*68.d0
        soa_tr%oxidant(isoa_oh)=c(ind_OH)
        soa_tr%oxidant(isoa_o3)=c(ind_O3)
        soa_tr%oxidant(isoa_no3)=c(ind_NO3)

! KPP chemistry
        CALL INTEGRATE( TIN = T, TOUT = T+DT, RSTATUS_U = RSTATE, &
        ICNTRL_U = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) )
        T = RSTATE(1)

! EQSAM
        yi = 0.d0
        yo = 0.d0
        yi(1, 1) = temp          ! temperature (K)
        yi(1, 2) = RH * 1.e-2    ! relative humidity (0-1)
        yi(1, 3) = (c(ind_nh3) + c(ind_nh4)) * eqsam_conv ! total ammonia/um (umol m-3)
        yi(1, 4) = c(ind_so4) * eqsam_conv ! total sulfate (umol m-3)
        yi(1, 5) = (c(ind_hno3) + c(ind_no3p)) * eqsam_conv ! total nitrate (umol m-3)
        yi(1, 11) = press * 1.e-2 ! pressure (hPa)
        call eqsam_v03d(yi,yo,nca,nco,iopt,loop,imax,ipunit)
        c(ind_hno3) = yo(1, 9) / eqsam_conv ! HNO3
        c(ind_nh3) = yo(1, 10) / eqsam_conv ! NH3
        c(ind_nh4) = yo(1, 19) / eqsam_conv ! total NH4 (aq+s)
        c(ind_no3p) = yo(1, 20) / eqsam_conv ! total NO3p (aq+s)
        c(ind_so4) = yo(1, 21) / eqsam_conv ! total SO4 (aq+s)

! SOA, after chemistry
        soa_cond%rr(isoa_isop,isoa_oh)=2.55d-11*exp(410.d0/temp) ! Isoprene + OH
        soa_cond%rr(isoa_isop,isoa_o3)=1.23d-14*exp(-2013.d0/temp) ! Isoprene + O3
!        soa_cond%rr(isoa_isop,isoa_no3)=7.8d-13 ! Isoprene + NO3
!        soa_cond%rr(isoa_apin,isoa_oh)=2.51d-11*exp(444.d0/temp) ! Terpenes + OH
        soa_cond%rr(isoa_apin,isoa_o3)=1.40d-14*exp(-732.d0/temp) ! Terpenes + O3
!        soa_cond%rr(isoa_apin,isoa_no3)=5.77d-13*exp(841.d0/temp) ! Terpenes + NO3
        do v = 1, soa_bins
          soa_tr%gas(v)=c(soa_tr%igas(v))*molec2ug*12.d0*1.6d0
          soa_tr%aer(v)=c(soa_tr%iaer(v))*molec2ug*12.d0*1.6d0
        enddo
        voc2nox_denom=4.2d-12*exp(180.d0/temp)*c(ind_NO)+&  ! XO2 + NO
                      3.5d-12*exp(250.d0/temp)*c(ind_HO2)+& ! XO2 + HO2
                      1.7d-14*exp(1300.d0/temp)*c(ind_XO2)  ! XO2 + XO2
        if (voc2nox_denom == 0.d0) then
          soa_cond%voc2nox=0.d0
        else
          soa_cond%voc2nox=4.2d-12*exp(180.d0/temp)*c(ind_NO)/&  ! XO2 + NO
                           voc2nox_denom
        endif
        soa_cond%temp=temp
        soa_cond%nvoa=molec2ug*(c(ind_SO4)*96.d0 + &
                                c(ind_NH4)*19.d0 + &
                                c(ind_NO3p)*63.d0)
        call soa_calc
        do v = 1, soa_bins
          c(soa_tr%igas(v))=soa_tr%gas(v)/(molec2ug*12.d0*1.6d0)
          c(soa_tr%iaer(v))=soa_tr%aer(v)/(molec2ug*12.d0*1.6d0)
        enddo

      END DO kron
!~~~> End Time loop

      CALL GetMass( C, DVAL )
      WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,     &
               ( TRIM(SPC_NAMES(MONITOR(i))),           &
                 C(MONITOR(i))/CFACTOR, i=1,NMONITOR ) 
      TIME = T
      CALL SaveData()
      CALL CloseSaveData()

991   FORMAT(F6.1,'%. T=',E10.3,2X,400(A,'=',E11.4,'; '))

END PROGRAM KPP_ROOT_Driver

