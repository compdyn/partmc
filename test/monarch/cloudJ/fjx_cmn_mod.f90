!------------------------------------------------------------------------------
!    'cmn_fjx_mod.f90'  for fast-JX code v 7.4+ (prather 8/15)
!         note that module and file begin with 'cmn_"
!         small changes: LREF=51 instead of hardwired, JVMAP replaces JMAP
!         MASFAC param added, also cloud params - see below
!------------------------------------------------------------------------------
!
! NB - ALL of these common variables are set paramters,
!    They are NOT to be used as variables for a local solution
!    Thus this entire set is 'in' only after it is initialized
!-----------------------------------------------------------------------
!
! !INTERFACE:
!
      MODULE FJX_CMN_MOD

      implicit none
      public

!-----------------------------------------------------------------------
      ! Turn on/off RRTM-G absorption cross sections
      logical, parameter::  LRRTMG= .true.
!-----------------------------------------------------------------------

      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
      integer, parameter ::  JXL_=100, JXL1_=JXL_+1
      ! JXL2_: 2*JXL_ + 2 = mx no.levels in basic FJX grid (mid-level)
      integer, parameter ::  JXL2_=2*JXL_+2
      ! W_   = dim = no. of Fast-J Wavelength bins:  currenly only 18, TROP-ONLY is done by zeroing FL fluxes
      integer, parameter ::  W_=18
      ! S_   = dim = number of wavelength bins INLCUDING the Solar-J extensions (RRTMG value = 27)
      integer, parameter ::  S_=27
      ! X_   = dim = max no. of X-section data sets (input data)
      integer, parameter ::  X_=72
      ! A_   = dim = max no. of Aerosol Mie sets (input data) not including clouds and SSA
      integer, parameter ::  A_=40
      ! SSA_   = dim = no. of strat sulfate aerosol types (input data)
      integer, parameter ::  SSA_=18
      ! C_   = dim = no. of cld-data sets (input data) ///currently just 4xLiq and 2xIce
      integer, parameter ::  C_=6
      ! N_  = no. of levels in Mie scattering arrays
      !     = 2*NC+1 = 4*(L_+1) + 1`+ 2*sum(JADDLV)
      integer, parameter ::  N_=601
      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      integer, parameter ::  M_=4
      ! M2_ = 2*M_ = 8, replaces MFIT
      integer, parameter ::  M2_=2*M_

!-----------------------------------------------------------------------
      ! 4 Gauss pts = 8-stream
      real*8, DIMENSION(M_), parameter  ::  &
                         EMU = [.06943184420297d0, .33000947820757d0, &
                                .66999052179243d0, .93056815579703d0]
      real*8, DIMENSION(M_), parameter  :: &
                         WT  = [.17392742256873d0, .32607257743127d0, &
                                .32607257743127d0, .17392742256873d0]
!-----------------------------------------------------------------------

      ! MASFAC: Conversion factor for pressure to column density
      real*8, parameter   ::  &
                         MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
      ! ZZHT: scale height (cm) used above top of CTM ZHL(LPAR+1)
      real*8, parameter   :: ZZHT = 5.d5
      ! RAD: Radius of Earth (cm)
      real*8, parameter   :: RAD = 6375.d5
      ! ATAU: heating rate (factor increase from one layer to the next)
      real*8, parameter   :: ATAU = 1.120d0
      ! ATAU0: minimum heating rate
      real*8, parameter   :: ATAU0 = 0.010d0
      ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
      integer, parameter  :: JTAUMX = (N_ - 4*JXL_)/2

      character*25, dimension(8), parameter :: TITCLD =  &
      ['clear sky - no clouds    ', &
       'avg cloud cover          ', &
       'avg cloud cover^3/2      ', &
       'ICAs - avg direct beam   ', &
       'ICAs - random N ICAs     ', &
       'QCAs - midpt of bins     ', &
       'QCAs - avg clouds in bins', &
       'ICAs - use all ICAs***   ']

!---- Variables in file 'FJX_spec.dat' (RD_XXX)
      ! WL: Centres of wavelength bins - 'effective wavelength'  (nm)
      real*8  WL(S_)
      ! WBIN: Boundaries of wavelength bins                  (microns)
      real*8  WBIN(S_+1)
      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      ! FW: Solar flux in W/m2
      ! FP: PAR quantum action spectrum
      real*8  FL(S_),FW(S_),FP(S_)
      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      real*8  QRAYL(S_)
      ! SJSUB:  intended for breakdown of the super-bins (1:27) into smaller sub-bins.
      real*8  SJSUB(S_,15)
      real*8  QO2(W_,3)   ! QO2: O2 cross-sections
      real*8  QO3(W_,3)   ! QO3: O3 cross-sections
      real*8  Q1D(W_,3)   ! Q1D: O3 => O(1D) quantum yield

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      real*8  QQQ(W_,3,X_)
      ! TQQ: Temperature for supplied cross sections
      real*8  TQQ(3,X_)
      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      integer LQQ(X_)

      ! TITLEJX: Title (short & long) for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*6  TITLEJX(X_)
      CHARACTER*16 TITLEJL(X_)
      ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1  SQQ(X_)

!---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)
      ! TITLAA: Aerosol Mie Titles
      character*12  TITLAA(A_)
      ! QAA: Aerosol scattering phase functions
      real*8  QAA(5,A_)
      ! WAA: 5 Wavelengths for the supplied phase functions
      real*8  WAA(5,A_)
      ! PAA: Phase function: first 8 terms of expansion
      real*8  PAA(8,5,A_)
      ! RAA: Effective radius associated with aerosol type
      real*8  RAA(A_)
      ! SAA: Single scattering albedo
      real*8  SAA(5,A_)
      ! DAA: density (g/cm^3)
      real*8  DAA(A_)
      ! NAA: Number of categories for scattering phase functions
      integer NAA

!---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)
      ! NCC: Number of categories for cloud scattering phase functions
      integer NCC
      ! TITLCC: Cloud type titles
      character*12  TITLCC(C_)
      ! RCC: Effective radius associated with cloud type
      real*8  RCC(C_)
      ! GCC: Effective geometric cross section
      real*8  GCC(C_)
      ! DCC: density (g/cm^3)
      real*8  DCC(C_)
      ! QCC: Cloud Q-ext
      real*8  QCC(S_,C_)
      ! WCC: Wavelengths for supplied phase functions
      real*8  WCC(S_,C_)
      ! SCC: Single scattering albedo
      real*8  SCC(S_,C_)
      ! PCC: Phase function: first 8 terms of expansion
      real*8  PCC(8,S_,C_)

!---- Variables in file 'FJX_scat-ssa.dat' (RD_SSA)
      ! NSS: Number of categories for Strat Sulf Aerosol scattering phase functions
      integer NSS
      ! TITLSS: Cloud type titles
      character*12  TITLSS(SSA_)
      ! RSS: Effective radius associated with cloud type
      real*8  RSS(SSA_)
      ! GSS: Effective geometric cross section
      real*8  GSS(SSA_)
      ! DSS: density (g/cm^3)
      real*8  DSS(SSA_)
      ! TSS: temperature (K)
      real*8  TSS(SSA_)
      ! WSS: weight percent sulfuric acid (%)
      real*8  WSS(SSA_)
      ! QSS: Q-ext      ----begin wavelength dependent quantities
      real*8  QSS(S_,SSA_)
      ! SSS: Single scattering albedo
      real*8  SSS(S_,SSA_)
      ! PSS: Phase function: first 8 terms of expansion
      real*8  PSS(8,S_,SSA_)

!---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)
      ! WMM: U Michigan aerosol wavelengths
      real*8  WMM(6)
      ! UMAER: U Michigan aerosol data sets
      real*8  UMAER(3,6,21,33)

!---- Variables in file 'atmos_std.dat' (RD_PROF)
      integer, parameter ::  LREF=51   ! layer dim. in reference profiles
      integer, parameter ::  JREF=18   ! latitude dim. in reference profiles

!----- T and O3, H2O, CH4, reference profiles added underscore _ because TREF used in RRTMG_SW
      real*8, DIMENSION(LREF,JREF,12) :: T_REF, O_REF, H2O_REF, CH4_REF
      integer NJX,NW1,NW2,NS1,NS2

!-----------------NEW for FJX72 parameters for cloud grid now here------
      integer, parameter :: &
            LPAR= 57, LWEPAR=34  &   !this can be set by CTM code
           ,L_=LPAR, L1_=L_+1 &   ! L_ = number of CTM layers
           ,L2_=2*L_+2 &        ! no. levels in the Fast-JX grid that
                       ! includes both layer edges and layer mid-points
           ,JVL_=LPAR &  ! vertical(levels) dim for J-values sent to CTM
           ,JVN_=101 &  ! max no. of J-values
           ,AN_=25     ! # FJX aerosols in layer (needs NDX for each)

!-----------------------------------------------------------------------
      ! variables used to map fast-JX J's onto CTM J's
!-----------------------------------------------------------------------
      real*8  JFACTA(JVN_)  ! multiplication factor for fast-JX calculated J
      integer JIND(JVN_)    ! index arrays that map Jvalue(j) onto rates
      integer NRATJ         ! number of Photolysis reactions in CTM chemistry, NRATJ <= JVN_
      integer NWBIN         ! used to shut off SFlux (and Fast-J calc) for strat wavelengths
                            ! only values that invoke action are 8 & 12 (mimic W_=8 or 12)
      integer NSBIN         ! likewise used to zero FL, if NSBIN = 18 then FL(19:27)=0.0

      character*6 JVMAP(JVN_) !label of J-value used to match w/FJX J's
      character*50 JLABEL(JVN_) ! label of J-value used in the chem model

! Cloud Cover parameters
!-----------------------------------------------------------------------
! NB CBIN_ was set at 20, but with NRG=6 groups, 10 gives a more reasonable number of ICAs
      integer, parameter :: CBIN_ = 10     ! # of quantized cloud fraction bins
      integer, parameter :: ICA_ = 20000   ! Max # of indep colm atmospheres
      integer, parameter :: NQD_ = 4       ! # of cloud fraction bins (4)

      real*8,  parameter ::  CPI    = 3.141592653589793d0
      real*8,  parameter ::  C2PI   = 2.d0*CPI
      real*8,  parameter ::  CPI180 = CPI/180.d0
      real*8,  parameter ::  G0     = 9.80665d0
      real*8,  parameter ::  G100 = 100.d0/G0
!-------data to set up the random number sequence for use in cloud-JX
      integer, parameter :: NRAN_ = 10007  ! dimension for random number
      real*4   RAN4(NRAN_)      ! Random number set

! Solar J parameters - for Cloud-J v7.4 just do one sub-bin per RRTMg superbin
!-----------------------------------------------------------------------
      integer, parameter:: W_r = S_-W_              ! Cloud-J fix, W_r = 82 for RRTMg in rrsw_fasj_cmn.f90
      real*8, dimension(L1_,W_r)  ::  TAUG_RRTMG    ! Cloud-J fix, is defined in rrsw_fasj_cmn.f90

      integer,  parameter, dimension(S_) ::   NGC =        &
        [  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,          &
           1,  1,  1,  1,  1,  1,  1,  1,  1,  1,          &
           1,  1,  1,  1,  1,  1,  1  ]
! actual RRMTMg values
!      integer,  parameter, dimension(S_) ::   NGC =       &
!       [   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,         &
!           1,  1,  1,  1,  1,  1,  1,  5, 10,  2,         &
!          10, 10,  8,  8, 12,  6, 12]


      END MODULE FJX_CMN_MOD
