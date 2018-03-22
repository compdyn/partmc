subroutine eqsam_v03d(yi,yo,nca,nco,iopt,loop,imax,ipunit,in) 
!
implicit none
!___________________________________________________________________________________________________________________________________
!      Written by Swen Metzger 3/11/99. Modified 2002, 2003.
!
!      Department of Atmospheric Chemistry, Max-Planck-Institute for Chemistry.
!      email: metzger@mpch-mainz.mpg.de
!      http://www.mpch-mainz.mpg.de/~metzger
!
!      COPYRIGHT 1999-2003
!
!      purpose
!      -------
!      EQSAM (EQuilibrium Simplified Aerosol Model) is a new and computationally efficient thermodynamic
!      aerosol composition model that allows to calculate the gas/aerosol equilibrium partitioning,
!      including aerosol water, sufficiently fast and accurate for global (or even regional) modeling.
!      EQSAM is based on a number of parameterizations, including single solute molalities and activity 
!      coefficients (AC). The thermodynamic framework (domains and subdomains, internally mixed aerosols) 
!      is the same as of more sophisticated thermodynamic equilibrium models (EQMs), e.g. of ISORROPIA 
!      (Nenes et al., 1998). Details are given in the references below (and the references therein).
!
!      The main assumption on which EQSAM/EQMs are based is thermodynamical and chemical equilibrium. 
!      From this assumption it directly follows that the aerosol water activity (aw) equals the ambient 
!      relative humidity (RH), if the water vapor pressure is sufficiently larger than the partial vapor
!      pressure of the aerosol compounds. This is approximately true for tropospheric aerosols. Given the 
!      large amount of water vapor present, water vapor and aerosol water equilibrate relatively faster 
!      compared to all other aerosol compounds. This is subsequently also true for single aerosol compounds.
!      The water activity of single solutes must also equal RH under this assumption. Therefore, the so 
!      called ZSR-relation is (and can be) used to calculate the aerosol associated water mass (simply
!      from the sum of all water mass fractions that are derived from measured single solute molalities). 
!
!      In contrast to other EQMs, EQSAM utilizes the fact that the RH fixes the water activity 
!      (under the above assumptions) and the consequence that any changes in RH also causes changes in 
!      the aerosol water mass and, hence, aerosol activity (including activity coefficients). Thus, an decrease
!      (increase) in RH decrease (increases) the aerosol water mass (and water activity). This can change the
!      aerosol composition, e.g. due to condensation (evaporation/crystallization), because the vapor pressure 
!      above the aerosol reduces (increases). In turn, a vapor pressure reduction (increase) due to changes
!      in the aerosol composition is compensated by an associated condensation (evaporation) of water vapor 
!      to maintain the aerosol molality to remain constant (because aw=RH). Furthermore, the aerosol water 
!      mainly depends on the aerosol mass and the type of solute, so that parameterizations of single solute 
!      molalities and activity coefficients can be defined, only depending on the type of solute and RH.
!      The advantage of using such parameterizations is that the entire aerosol equilibrium composition 
!      can be solved analytically, i.e. non-iteratively, which considerably reduces the amount of CPU time 
!      that is usually need for aerosol thermodynamic calculations (especially if an EQM is incorporated in
!      an aerosol dynamical model that is in turn embedded in a high resolution regional or global model).
!
!      However, EQSAM should still be regarded as a starting point for further developments. There is still 
!      room for improvements. For instance, this code is not yet numerically optimized (vectorized) and a 
!      number of improvements with respect to an explicit treatment of additional equilibrium reactions,
!      missing (or only implicit) dissociation, and a basic parameterization of the water uptake. 
!      
!      Note that EQSAM was originally developed to calculate the gas/aerosol equilibrium partitioning of the 
!      ammonium-sulfate-nitrate-water system for climate models, excluding solid compounds. 
!      This version (eqsam_v03d.f90) is extended with respect to sea salt. Solids/hysteresis are treated in a 
!      simplified manner. Results of a box model comparison with ISORROPIA will be available from the web page.
!      Please also note that the water uptake is based on additional (unpublished) parameterizations for single 
!      solute molalities, which are derived from tabulated measurements used in ISORROPIA. Note further that 
!      this extended version (eqsam_v03d.f90) is not yet published. A publication is in progress.
!
! ToDo:
!     Split ion-pairs into ions for water parameterizations (since info is actually available)
!     Include uptake/dissociation of NH3, HNO3, HCl (mainly to get pH right at near neutral conditions)
!     Extension to K+,Ca++,Mg++, CO2/(CO3)2--/HCO3-,SOA,etc.. (maybe not)
!     Vectorization. Translation of hardcoded formulas in array syntax.
!     I/O Interface and program structure clean up.
!     EQSAM info webpage.
!
! Version History:
!
!  eqsam_v03d.f90 (MPI-CH, June 2003): 
!   - gama parameterizations now according to Metzger 2002 (JGR Appendix)
!   - improved pH calculations (still restricted to strong acids)
!   - removed bug that lead to too high nitrate formation at dry and cold regions (UT/LS) 
!   - removed bug in solid/hysteresis calculations 
!     (both bugs introduced in eqsam_v03b.f90 by cleaning up eqsam_v02a.f90)
!   
!  eqsam_v03c.f90 (MPI-CH, April 2003):
!   - more accurate paramterizations of single solute molalities (Na, Cl species)
!   - cleanded up RHD subdomain structure
!   - improved water uptake (Na, Cl species)
!
!  eqsam_v03b.f90 (MPI-CH, March 2003): 
!                 System extended to HCl,Cl-/Na+. 
!                 Parameterization (fit) of additional HNO3 uptake removed.
!                 Instead, complete analytical solution of equilibrium reactions, based on the AC-RH relationship.
!  eqsam_v03.f90  (IMAU, October 1999): 
!                 Test version (included in TM3).
!  eqsam_v02a.f90 (IMAU, April 2000):
!                 Box model version.
!  eqsam_v02.f90  (IMAU, October 1999):
!                 TM3 version.
!                 Version including solids and additional HNO3 uptake on acidic aerosols (parameterized).
!  eqsam_v01b.f90 (MPI-CH, January 2003):
!                 Same as eqsam_v01a.f90 (additional lines though uncommented for test purposes only).
!  eqsam_v01a.f90 (IMAU, April 2000):
!                 Box model version.
!  eqsam_v01.f90  (IMAU, October 1999):
!                 TM3 version.
!                 First and most basic version (without solids) for better vectorization (for global modeling).
!                 System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, H2O 
!                 based on equilibrium / internal mixture assumption / aw=rh / ZSR-relation
!                 parameterization of activcity coefficients (AC), i.e. an AC-RH relationship
!
!      
!      interface
!      ---------
!      call  eqsam_v03d(yi,yo,nca,nco,iopt,loop,imax,ipunit,in)
!
!      yi = input  array (imax, nca)
!      yo = output array (imax, nco)
!      imax = max loop (e.g. time steps)
!      nca >= 11
!      nc0 >= 35
!      iopt = 1 metastable 
!      iopt = 2 solids 
!      iopt = 3 hysteresis (metastable/solids) for online calculations
!      iopt = 31 hysteresis lower branch 
!      iopt = 32 hysteresis upper branch 
!      ipunit = I/O unit (can be skipped)
!      in = array        (can be skipped)
!         
!      method
!      ------
!      equilibrium / internal mixture assumption / aw=rh
!      System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, HCl,Cl-/Na+, H2O 
!              (K+,Ca++,Mg++)
!      external
!      --------
!      program    eqmd.f90    (driver only needed for the box model version)
!      subroutine gribio.f90  (provides diagnostics output in grib/binary/ascii format)
!      
!      references
!      ---------
!      Swen Metzger Ph.D Thesis, University Utrecht, 2000.
!         http://www.library.uu.nl/digiarchief/dip/diss/1930853/inhoud.htm
!
!      Metzger, S. M., F. J. Dentener, J. Lelieveld, and S. N. Pandis, 
!         GAS/AEROSOL PARTITIONING I: A COMPUTATIONALLY EFFICIENT MODEL, 
!         J Geophys. Res., 107, D16, 10.1029/2001JD001102, 2002
!         http://www.agu.org/journals/jd/jd0216/2001JD001102/index.html
!      Metzger, S. M., F. J. Dentener, A. Jeuken, and M. Krol, J. Lelieveld, 
!         GAS/AEROSOL PARTITIONING II: GLOBAL MODELING RESULTS, 
!         J Geophys. Res., 107, D16, 10.1029/2001JD001103, 2002.
!         http://www.agu.org/journals/jd/jd0216/2001JD001103/index.html
!___________________________________________________________________________________________________________________________________
real,parameter                                :: RH_HIST_DW=1.50                                    ! mean value for mixture of wet (2) and dry (1) gridboxes (needed for HYSTERESIS)
real,parameter                                :: T0=298.15,T1=298.0,AVO=6.03e23,R=82.0567e-6,     & ! in cu.m*atm/deg/mole
                                                 r_kcal  = 1.986E-3                                 ! Ideal gas constant [kcal K-1.mole-1]
real,parameter                                :: RHMAX=0.99,RHMIN=0.0001                            ! restrict to max / min RH
real,parameter                                :: MWNH4=18.,MWSO4=96.,MWNO3=62.,MWCl=35.5            ! mole mass of species considered
real,parameter                                :: MWNa=23.0,MWCa=40.1,MWN=14.0, MWS=32.1
real,parameter                                :: MWH20=55.51*18.01,ZERO=0.0
real,parameter                                :: GF1=0.25,GF2=0.50,GF3=0.40,GF4=1.00,K=2.           ! exponents of AC-RH functions
!______________________________________________
integer,parameter                             :: NPAIR=10 
!
integer                                       :: ii,il,IHYST
integer,intent(in)                            :: nca,nco,imax,loop,ipunit
integer,intent(inout)                         :: iopt
!______________________________________________
integer,dimension(6),intent(in)               :: in
!______________________________________________
real                                          :: T0T,TT,RH,PX,RHD,KAN,KAC,ZIONIC,RH_HIST,GAMA,GG,GF,GFN
real                                          :: X00,X01,X02,X03,X04,X05,X08,X09,X10,X11
real                                          :: X0,X1,X2,X3,X4,X5,X6,XK10,XK6
real                                          :: ZFLAG,ZKAN,ZKAC,PH,COEF,HPLUS,AKW,XKW,MOLAL
real                                          :: TNH4,TSO4,TNO3,TNa,TCl,TPo,TCa,TMg
real                                          :: PNH4,PSO4,PNO3,PCl,PNa,GNO3,GNH3,GSO4,GHCl
real                                          :: ASO4,ANO3,ANH4,ACl,ANa,SNH4,SSO4,SNO3,SCl,SNa
real                                          :: WH2O,PM,PMs,PMt,RINC,DON,RATIONS,GR,NO3P,NH4P
!_______________________________________________
real,dimension(imax,nca),intent(in)           :: yi
real,dimension(imax,nco),intent(out)          :: yo
real,dimension(8)                             :: w1,w2
real,dimension(8)                             :: RHDA,RHDE,RHDX,RHDZ    ! RHD / MRHD arrays for different aerosol types
real,dimension(NPAIR)                         :: M0,MW,NW,ZW            ! arrays of ion pairs
!
! salt solutes:
!   1 = NACl,  2 = (NA)2SO4, 3 = NANO3,  4 = (NH4)2SO4,  5 = NH4NO3, 6 = NH4CL,   7 = 2H-SO4
!   8 = NH4HSO4,   9 = NAHSO4, 10 = (NH4)3H(SO4)2
!
DATA MW(1:NPAIR)/ 58.5, 142.0, 88.0, 132.0, 80.0, 53.5, 98.0, 115.0, 120.0, 247.0/ ! mole mass of the salt solute
DATA NW(1:NPAIR)/  2.0,   2.5,  2.5,   2.5,  3.5,  1.0,  4.5,   2.0,   2.0,   2.5/ ! square of max. dissocation number (not consistent)
DATA ZW(1:NPAIR)/ 0.67,   1.0,  1.0,   1.0,  1.0,  1.0,  0.5,   1.0,   1.0,   1.0/ ! exponents of water activity functions
!
DATA RHDA(1:8)/0.32840, 0.4906, 0.6183, 0.7997, 0.67500, 0.5000, 0.4000, 0.0000/ ! RHD / MRHD values as of ISORROPIA / SCAPE (T=298.15K)
DATA RHDE(1:8)/-1860.0, -431.0, 852.00, 80.000, 262.000, 3951.0, 384.00, 0.0000/ ! Temp. coeff.
!___________________________________________________________________________________________________________________________________
IHYST=2
IF(IOPT.EQ.31) THEN      ! SOLID HYSTORY
   IHYST=1
   IOPT=3
ELSEIF(IOPT.EQ.32) THEN  ! WET   HISTORY
   IHYST=2
   IOPT=3
ENDIF
!-------------------------------------------------------------------------------
!DLW:010406: Commented out print block.
!-------------------------------------------------------------------------------
!write(ipunit,*)'eqsam_v03d ...'
!print*,'                                                          '
!print*,'              EQuilibrium Simplified Aerosol Model (EQSAM)'
!print*,'                        for global modeling               ' 
!print*,'                                 by                       '
!print*,'                         Swen Metzger, MPI-CH             '
!print*,'                         Copyright 1999-2003              '
!print*,'                    >> metzger@mpch-mainz.mpg.de <<       '
!print*,'                     last change:  04. June, 2003         '
!print*,'                           (version 3.0d)                 '
!print*,'                   gas/aerosol calculations assuming      '
!print*,'                  System: NH3,NH4+/H2SO4+,HSO4-,SO4--     '
!print*,'                      HNO3,NO3-, HCl,Cl-/Na+, H2O         '
!if(iopt.eq.1) then
!print*,'                         metastable aeorsols              '
!elseif(iopt.eq.2) then
!print*,'                            solid aeorsols                '
!elseif(iopt.eq.3) then
!print*,'                             hysteresis                   '
!print*,'                         (metastable/solids)              '
!if(IHYST.eq.1) then
!print*,'                            solid hystory                 '
!elseif(IHYST.eq.2) then
!print*,'                             wet hystory                  '
!endif
!endif
!print*,'                                                          '
!print*,'loop over ',loop,' data sets'
!print*,'   '
!-------------------------------------------------------------------------------
!DLW:010406: End comment out print block.
!-------------------------------------------------------------------------------
!___________________________________________________________________________________________________________________________________
yo=0.;w1=0.;w2=0.                                                                                   ! init/reset
!___________________________________________________________________________________________________________________________________
do il=1,loop

! get old relative humidity to calculate aerosol hysteresis (online only)

   RH_HIST = 2.                                                                                     ! WET HISTORY (DEFAULT)
   IF(IHYST.EQ.1.OR.IOPT.EQ.2)  RH_HIST = 1.                                                        ! SET TO SOLIDS

!  meteorology
   TT = yi(il,1)                    ! T                      [K]
   RH = yi(il,2)                    ! RH                     [0-1]
   PX = yi(il,11)                   ! p                      [hPa]
!
! gas+aerosol:
   w1(1) = yi(il,6)                 ! Na+ (ss  + xsod) (a)   [umol/m^3 air]
   w1(2) = yi(il,4)                 ! H2SO4    + SO4-- (p)   [umol/m^3 air]
   w1(3) = yi(il,3)                 ! NH3  (g) + NH4+  (p)   [umol/m^3 air]
   w1(4) = yi(il,5)                 ! HNO3 (g) + NO3-  (p)   [umol/m^3 air]
   w1(5) = yi(il,7)                 ! HCl  (g) + Cl-   (p)   [umol/m^3 air]
   w1(6) = yi(il, 8)                ! K+   (p) from Dust     [umol/m^3 air]
   w1(7) = yi(il, 9)                ! Ca++ (p) from Dust     [umol/m^3 air]
   w1(8) = yi(il,10)                ! Mg++ (p) from Dust     [umol/m^3 air]
!______________________________________________

   zflag=1.

   w1=w1*1.0e-6                     ! [mol/m^3 air]

   TNa   = w1(1)                    ! total input sodium   (g+p) 
   TSO4  = w1(2)                    ! total input sulfate  (g+p) 
   TNH4  = w1(3)                    ! total input ammonium (g+p)
   TNO3  = w1(4)                    ! total input nitrate  (g+p) 
   TCl   = w1(5)                    ! total input chloride (g+p) 
   TPo   = w1(6)                    ! total input potasium (g+p) 
   TCa   = w1(7)                    ! total input calcium  (g+p)
   TMg   = w1(8)                    ! total input magnesium(g+p)

! SULFATE RICH

      if((w1(1)+w1(3)+w1(6)+2.*(w1(7)+w1(8))).le.(2.*w1(2))) then
          zflag=3.
      endif

! SULFATE VERY RICH CASE if (NH4+Na+K+2(Ca+Mg))/SO4 < 1

      if((w1(1)+w1(3)+w1(6)+2.*(w1(7)+w1(8))).le.w1(2)) then
          zflag=4.
      endif

! SULFATE NEUTRAL CASE

      if((w1(1)+w1(3)+w1(6)+2.*(w1(7)+w1(8))).gt.(2.*w1(2))) then
          zflag=2.
      endif

! SULFATE POOR AND CATION POOR CASE

      if((w1(1)+w1(6)+2.*(w1(7)+w1(8))).gt.(2.*w1(2))) then       
          zflag=1.
      endif

      IF ( RH .LT. RHMIN ) RH=RHMIN
      IF ( RH .GT. RHMAX ) RH=RHMAX

! CALCULATE TEMPERATURE DEPENDENCY FOR SOME RHDs

      RHDX(:)=RHDA(:)*exp(RHDE(:)*(1./TT-1./T0))
      RHDZ(:)=RHDX(:)
      
! ACCOUNT FOR VARIOUS AMMOMIUM/SODIUM SULFATE SALTS ACCORDING TO MEAN VALUE AS OF ISORROPIA

      GG=2.0                          ! (Na)2SO4 / (NH4)2SO4 IS THE PREFFERED SPECIES FOR SULFATE DEFICIENT CASES
      IF(ZFLAG.EQ.3.) THEN
         IF(RH.LE.RHDZ(7)) THEN       ! ACCOUNT FOR MIXTURE OF (NH4)2SO4(s) & NH4HSO4(s) & (NH4)3H(SO4)2(s) 
            GG=1.677                  !                        (Na)2SO4 &  NaHSO4
!           GG=1.5
         ELSEIF(RH.GT.RHDZ(7).AND.RH.LE.RHDZ(5)) THEN ! MAINLY (Na)2SO4 / (NH4)2SO4(s) & (NH4)3H(SO4)2(s)
            GG=1.75
!           GG=1.5
         ELSEIF(RH.GE.RHDZ(5)) THEN   ! (NH4)2SO4(S) & NH4HSO4(S) & SO4-- & HSO4-
            GG=1.5                    !  (Na)2SO4 &  NaHSO4
         ENDIF
      ENDIF
      IF(ZFLAG.EQ.4.) GG=1.0          ! IF SO4 NEUTRALIZED, THEN ONLY AS NaHSO4 / NH4HSO4(S) OR  HSO4- / H2SO4

      RHD=RH
      IF(IOPT.EQ.2.OR.RH_HIST.LT.RH_HIST_DW) THEN   ! GET RHD FOR SOLIDS / HYSTERESIS
!
! GET LOWEST DELIQUESCENCE RELATIVE HUMIDITIES ACCORDING TO THE CONCENTRATION DOMAIN (APROXIMATION) 
! BASED ON RHD / MRHD ISORROPIA/SCAPE
!
      w2(:)=1.
      do ii=1,8
         if(w1(ii).le.1.e-12) w2(ii)=0.  ! skip compound in RHD calculation if value is concentration is zero or rather small
      enddo

! GET LOWEST RHD ACCORDING TO THE CONCENTRATION DOMAIN

! zflag=1. (cation rich)  ...
! 1. sea salt      aerosol          : RHDX(1)=MgCl2
! 2. mineral dust  aerosol          : RHDX(2)=Ca(NO3)2
!
! zflag=2. (sulfate neutral) ...
! 3. ammonium + nitrate             : RHDX(3)= NH4NO3
! 4. ammonium + sulfate             : RHDX(4)=(NH4)2SO4        
! 5. ammonium + sulfate mixed salt  : RHDX(5)=(NH4)3H(SO4)2, (NH4)2SO4        
! 6. ammonium + nitrate  + sulfate  : RHDX(6)=(NH4)2SO4, NH4NO3, NA2SO4, NH4CL
!
! zflag=3. (sulfate poor) ...
! 7. ammonium + sulfate  (1:1,1.5)  : RHDX(7)= NH4HSO4
!
! zflag=4. (sulfate very poor) ...
! 8. sulfuric acid                  : RHDX(8)= H2SO4       

!WRITE(IPUNIT,*)'zflag = ', ZFLAG     !DLW
!WRITE(IPUNIT,*)'GG    = ', GG        !DLW
!WRITE(IPUNIT,*)'w1(1:8) = ', W1(1:8) !DLW
!WRITE(IPUNIT,*)'w2(1:8) = ', W2(1:8) !DLW

   IF(ZFLAG.EQ.1.)THEN

      RHD=W2(1)+W2(5)                     ! Na+  dependency
      IF(RHD.EQ.0.)  RHDX(1)=1. 
      RHD=W2(6)+W2(7)+W2(8)               ! K+/Ca++/Mg++ dependency (incl. ss)
      IF(RHD.EQ.0.)  RHDX(2)=1. 

      RHD=MINVAL(RHDX(1:2))

   ELSEIF(ZFLAG.EQ.2.)THEN

      RHD=W2(3)*W2(4)                     ! NH4+ & NO3-  dependency
      IF(RHD.EQ.0.)  RHDX(3)=1. 
      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(GG.NE.2.)   RHD=0.               ! account only for pure (NH4)2SO4
      IF(RHD.EQ.0.)  RHDX(4)=1. 
      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(RHD.EQ.0.)  RHDX(5)=1. 
      RHD=W2(2)+W2(3)+W2(4)+W2(5)         ! (NH4)2SO4, NH4NO3, NA2SO4, NH4CL dependency
      IF(RHD.EQ.0.)  RHDX(6)=1. 

!     RHD=MINVAL(RHDX(3:4))
      RHD=MINVAL(RHDX(3:6))

   ELSEIF(ZFLAG.EQ.3.)THEN

      RHD=W2(2)+W2(3)                     ! NH4+ & SO4-- dependency
      IF(RHD.EQ.0.)  RHDX(7)=1. 
      RHD=RHDX(7)

   ELSEIF(ZFLAG.EQ.4.)THEN

      RHD=W2(2)                           ! H2SO4 dependency (assume no dry aerosol)
      IF(RHD.EQ.0.)  RHDX(8)=1. 

      RHD=RHDX(8)

   ENDIF ! ZFLAG
   ! WRITE(IPUNIT,*)'RHDX(1:8) = ', RHDX(1:8)   !DLW
   ENDIF ! SOLIDS

! GET WATER ACTIVITIES ACCORDING TO METZGER, 2000.
! FUNCTION DERIVED FROM ZSR RELATIONSHIP DATA (AS USED IN ISORROPIA)

      M0(:) = ((NW(:)*MWH20/MW(:)*(1./RH-1.)))**ZW(:)

! CALCULATE TEMPERATURE DEPENDENT EQUILIBRIUM CONSTANTS

      T0T=T0/TT
      COEF=1.0+LOG(T0T)-T0T

! EQUILIBRIUM CONSTANT NH4NO3(s) <==> NH3(g) + HNO3(g) [atm^2] (ISORROPIA)

      XK10 = 5.746e-17
      XK10= XK10 * EXP(-74.38*(T0T-1.0) + 6.120*COEF)
      KAN = XK10/(R*TT)/(R*TT)

! EQUILIBRIUM CONSTANT  NH4CL(s) <==> NH3(g) + HCL(g) [atm^2] (ISORROPIA)

      XK6  = 1.086e-16
      XK6 = XK6 * EXP(-71.00*(T0T-1.0) + 2.400*COEF)
      KAC = XK6/(R*TT)/(R*TT)

!
! CALCULATE AUTODISSOCIATION CONSTANT (KW) FOR WATER H2O <==> H(aq) + OH(aq) [mol^2/kg^2] (ISORROPIA)

      XKW  = 1.010e-14
      XKW = XKW *EXP(-22.52*(T0T-1.0) + 26.920*COEF)

! GET MEAN MOLAL IONIC ACTIVITY COEFF ACCORDING TO METZGER, 2002.

      GAMA=0.0
      IF(RH.GE.RHD)                GAMA=(RH**ZFLAG/(1000./ZFLAG*(1.-RH)+ZFLAG))
      GAMA = GAMA**GF1            ! ONLY GAMA TYPE OF NH4NO3, NaCl, etc. NEEDED SO FAR

      GAMA=0.0
      GFN=K*K                      ! K=2, i.e. condensation of 2 water molecules per 1 mole ion pair
      GF=GFN*GF1                   ! = GFN[=Nw=4] * GF1[=(1*1^1+1*1^1)/2/Nw=1/4] = 1
                                   ! ONLY GAMA TYPE OF NH4NO3, NH4Cl, etc. needed so far

      IF(RH.GE.RHD) GAMA=RH**GF/((GFN*MWH20*(1./RH-1.)))**GF1

      GAMA = min(GAMA,1.0)        ! FOCUS ON 0-1 SCALE
      GAMA = max(GAMA,0.0)
      GAMA = (1.-GAMA)**K          ! transplate into aqueous phase equillibrium and account for 
                                   ! enhanced uptake of aerosol precursor gases with increasing RH
                                   ! (to match the results of ISORROPIA)

! CALCULATE RHD DEPENDENT EQ: IF RH <  RHD => NH4NO3(s) <==> NH3 (g) + HNO3(g) (ISORROPIA)
!                             IF RH >> RHD => HNO3  (g)   -> NO3 (aq)

      X00  = MAX(ZERO,MIN(TNa,GG*TSO4))       ! MAX SODIUM   SULFATE
      X0   = MAX(ZERO,MIN(TNH4,GG*TSO4-X00))  ! MAX AMMOMIUM SULFATE
      X01  = MAX(ZERO,MIN(TNa-X00, TNO3))     ! MAX SODIUM   NITRATE
      X1   = MAX(ZERO,MIN(TNH4-X0,TNO3-X01))  ! MAX AMMOMIUM NITRATE
!
      X02  = MAX(ZERO,MIN(TNa-X01-X00,TCl))   ! MAX SODIUM   CHLORIDE
      X03  = MAX(ZERO,MIN(TNH4-X0-X1,TCl-X02))! MAX AMMOMIUM CHLORIDE

      X2   = MAX(TNH4-X1-X0-X03,ZERO)         ! INTERIM RESIDUAL NH3
      X3   = MAX(TNO3-X1-X01,ZERO)            ! INTERIM RESIDUAL HNO3
      X04  = MAX(TSO4-(X0+X00)/GG,ZERO)       ! INTERIM RESIDUAL H2SO4
      X05  = MAX(TCl-X03-X02,ZERO)            ! INTERIM RESIDUAL HCl
!     X06  = MAX(TNa-X02-X01-X00,ZERO)        ! INTERIM RESIDUAL Na (should be zero for electro-neutrality in input data)
!
      ZKAN=2.
      IF(RH.GE.RHD) ZKAN=ZKAN*GAMA

      X4   = X2 + X3
      X5   = SQRT(X4*X4+KAN*ZKAN*ZKAN)
      X6   = 0.5*(-X4+X5)
      X6   = MIN(X1,X6)
      
      GHCl = X05                              ! INTERIM RESIDUAl HCl
      GNH3 = X2 + X6                          ! INTERIM RESIDUAl NH3
      GNO3 = X3 + X6                          ! RESIDUAl HNO3
      GSO4 = X04                              ! RESIDUAl H2SO4
      PNa  = X02 + X01 + X00                  ! RESIDUAl Na (neutralized)
      
      ZKAC=2.
      IF(RH.GE.RHD) ZKAC=ZKAC*GAMA

      X08   = GNH3 + GHCl
      X09   = SQRT(X08*X08+KAC*ZKAC*ZKAC)
      X10   = 0.5*(-X08+X09)
      X11   = MIN(X03,X10)

      GHCl = GHCl + X11                       ! RESIDUAL HCl
      GNH3 = GNH3 + X11                       ! RESIDUAL NH3

! GO SAVE ...

      IF(GHCl.LT.0.)   GHCl=0.
      IF(GSO4.LT.0.)   GSO4=0.
      IF(GNH3.LT.0.)   GNH3=0.
      IF(GNO3.LT.0.)   GNO3=0.
      IF(PNa.LT.0.)    PNa=0.
      IF(GSO4.GT.TSO4) GSO4=TSO4
      IF(GNH3.GT.TNH4) GNH3=TNH4
      IF(GNO3.GT.TNO3) GNO3=TNO3
      IF(GHCl.GT.TCl)  GHCl=TCl
      IF(PNa.GT.TNa)   PNa=TNa
!     IF(PNa.LT.TNa)   print*,il,' PNa.LT.TNa => no electro-neutrality in input data! ',PNa,TNa


! DEFINE AQUEOUSE PHASE (NO SOLID NH4NO3 IF NO3/SO4>1, TEN BRINK, ET AL., 1996, ATMOS ENV, 24, 4251-4261)

!     IF(TSO4.EQ.ZERO.AND.TNO3.GT.ZERO.OR.TNO3/TSO4.GE.1.) RHD=RH

!     IF(IOPT.EQ.2.AND.RH.LT.RHD.OR.IOPT.EQ.2.AND.RH_HIST.LT.RH_HIST_DW) THEN        ! SOLIDS / HYSTERESIS
      IF(RH_HIST.EQ.1.AND.RH.LT.RHD) THEN        ! SOLIDS / HYSTERESIS

       ! EVERYTHING DRY, ONLY H2SO4 (GSO4) REMAINS IN THE AQUEOUSE PHASE

         ANH4 = 0.
         ASO4 = 0.
         ANO3 = 0.
         ACl  = 0.
         ANa  = 0.

      ELSE  !  SUPERSATURATED SOLUTIONS NO SOLID FORMATION

         ASO4 = TSO4 - GSO4
         ANH4 = TNH4 - GNH3
         ANO3 = TNO3 - GNO3
         ACl  = TCl  - GHCl
         ANa  = PNa

      ENDIF ! SOLIDS / HYSTERESIS

! CALCULATE AEROSOL WATER [kg/m^3(air)]
!
! salt solutes:
!   1 = NACl,  2 = (NA)2SO4, 3 = NANO3,  4 = (NH4)2SO4,  5 = NH4NO3, 6 = NH4CL,   7 = 2H-SO4
!   8 = NH4HSO4,   9 = NAHSO4, 10 = (NH4)3H(SO4)2
!
      IF(ZFLAG.EQ.1.) WH2O = ASO4/M0( 2) + ANO3/M0(3) + ACl/M0(6)
      IF(ZFLAG.EQ.2.) WH2O = ASO4/M0( 9) + ANO3/M0(5) + ACl/M0(6)
      IF(ZFLAG.EQ.3.) WH2O = ASO4/M0( 8) + ANO3/M0(5) + ACl/M0(6)
      IF(ZFLAG.EQ.4.) WH2O = ASO4/M0( 8) + GSO4/M0(7)

! CALCULATE AQUEOUS PHASE PROPERTIES

!     PH    = 9999.
      PH    = 7.
      MOLAL = 0.
      HPLUS = 0.
      ZIONIC= 0.

      IF(WH2O.GT.0.) THEN            

         ! CALCULATE AUTODISSOCIATION CONSTANT (KW) FOR WATER

         AKW=XKW*RH*WH2O*WH2O                         ! H2O <==> H+ + OH- with kw [mol^2/kg^2]
         AKW=AKW**0.5                                 ! [OH-] = [H+] [mol]
         AKW=MAX( AKW, 1.0E-30 )                      ! DLW, 11-14-06, to prevent division by 0 in HPLUS below.

         ! Calculate hydrogen molality [mol/kg], i.e. H+ of the ions:
         !                                   Na+, NH4+, NO3-, Cl-, SO4--, HH-SO4- [mol/kg(water)]
         !                                   with [OH-] = kw/[H+]

         HPLUS = (-ANa/WH2O-ANH4/WH2O+ANO3/WH2O+ACl/WH2O+GG*ASO4/WH2O+GG*GSO4/WH2O+ & 
            SQRT(( ANa/WH2O+ANH4/WH2O-ANO3/WH2O-ACl/WH2O-GG*ASO4/WH2O-GG*GSO4/WH2O)**2+XKW/AKW*WH2O))/2.

         ! Calculate pH

         PH=-ALOG10(HPLUS)                             ! aerosol pH 

         ! Calculate ionic strength [mol/kg]

         ZIONIC=0.5*(ANa+ANH4+ANO3+ACl+ASO4*GG*GG+GSO4*GG*GG+XKW/AKW*WH2O*WH2O)
         ZIONIC=ZIONIC/WH2O                                           ! ionic strength [mol/kg]
!        ZIONIC=min(ZIONIC,200.0)                                     ! limit for output
!        ZIONIC=max(ZIONIC,0.0)

      ENDIF ! AQUEOUS PHASE 
!
!-------------------------------------------------------
! calculate diagnostic output consistent with other EQMs ...
!
      ASO4 = ASO4 + GSO4                                     ! assuming H2SO4 remains aqueous

      TNa   = TNa  * 1.e6                                    ! total input sodium   (g+p)  [umol/m^3]
      TSO4  = TSO4 * 1.e6                                    ! total input sulfate  (g+p)  [umol/m^3]
      TNH4  = TNH4 * 1.e6                                    ! total input ammonium (g+p)  [umol/m^3]
      TNO3  = TNO3 * 1.e6                                    ! total input nitrate  (g+p)  [umol/m^3]
      TCl   = TCl  * 1.e6                                    ! total input chloride (g+p)  [umol/m^3]
      TPo   = TPo  * 1.e6                                    ! total input potasium (g+p)  [umol/m^3]
      TCa   = TCa  * 1.e6                                    ! total input calcium  (g+p)  [umol/m^3]
      TMg   = TMg  * 1.e6                                    ! total input magnesium(g+p)  [umol/m^3]
!
! residual gas:
      GNH3 = GNH3 * 1.e6                                     ! residual NH3
      GSO4 = GSO4 * 1.e6                                     ! residual H2SO4
      GNO3 = GNO3 * 1.e6                                     ! residual HNO3
      GHCl = GHCl * 1.e6                                     ! residual HCl

! total particulate matter (neutralized)
      PNH4=TNH4-GNH3                                         ! particulate ammonium   [umol/m^3]
      PNO3=TNO3-GNO3                                         ! particulate nitrate    [umol/m^3]
      PCl =TCl -GHCl                                         ! particulate chloride   [umol/m^3]
      PNa =TNa                                               ! particulate sodium     [umol/m^3]
      PSO4=TSO4                                              ! particulate sulfate    [umol/m^3]

! liquid matter
      ASO4 = ASO4 * 1.e6                                     ! aqueous phase sulfate  [umol/m^3] 
      ANH4 = ANH4 * 1.e6                                     ! aqueous phase ammonium [umol/m^3]
      ANO3 = ANO3 * 1.e6                                     ! aqueous phase nitrate  [umol/m^3]
      ACl  = ACl  * 1.e6                                     ! aqueous phase chloride [umol/m^3]
      ANa  = ANa  * 1.e6                                     ! aqueous phase sodium   [umol/m^3]

! solid matter
      SNH4=PNH4-ANH4                                         ! solid phase ammonium   [umol/m^3]
      SSO4=PSO4-ASO4                                         ! solid phase sulfate    [umol/m^3]
      SNO3=PNO3-ANO3                                         ! solid phase nitrate    [umol/m^3]
      SCl =PCl -ACl                                          ! solid phase chloride   [umol/m^3]
      SNa =PNa -ANa                                          ! solid phase sodium     [umol/m^3]

! GO SAVE ...

      IF(SNH4.LT.0.)   SNH4=0.
      IF(SSO4.LT.0.)   SSO4=0.
      IF(SNO3.LT.0.)   SNO3=0.
      IF(SCl.LT.0.)    SCl=0.
      IF(SNa.LT.0.)    SNa=0.

      PM=SNH4+SSO4+SNO3+SNH4+SCl+SNa+ANH4+ASO4+ANO3+ACl+ANa  ! total PM [umol/m^3]
      PMs=SNH4*MWNH4+SSO4*MWSO4+SNO3*MWNO3+SCl*MWCl+SNa*MWNa ! dry particulate matter (PM)     [ug/m^3]
      PMt=PMs+ANH4*MWNH4+ASO4*MWSO4+ANO3*MWNO3+ACl*MWCl+  &
          ANa*MWNa                                           ! total (dry + wet) PM, excl. H20 [ug/m^3]

      WH2O = WH2O * 1.e9                                     ! convert aerosol water from [kg/m^3] to [ug/m^3]
      IF(WH2O.LT.1.e-3) WH2O=0.

! UPDATE HISTORY RH FOR HYSTERESIS (ONLINE CALCULATIONS ONLY)

      RH_HIST=2.                                             ! wet
      IF(WH2O.EQ.0.) RH_HIST=1.                              ! dry

      RINC = 1.
      IF(PMt.GT.0.)   RINC = (WH2O/PMt+1)**(1./3.)           ! approx. radius increase due to water uptake
      IF(RINC.EQ.0.)  RINC = 1.

      RATIONS = 0.
      IF(PSO4.GT.0.) RATIONS = PNO3/PSO4                     ! nitrate / sulfate mol ratio

      GR = 0.
      IF(GNO3.GT.0.) GR = GNH3/GNO3                          ! gas ratio = residual NH3 / residual HNO3   [-]

      DON = 0.
      IF((PNO3+2.*PSO4).GT.0.) DON = 100.*PNH4/(PNO3+2.*PSO4)! degree of neutralization by ammonia : ammonium / total nitrate + sulfate  [%]

      NO3P = 0.
      IF(TNO3.GT.0.) NO3P = 100.*PNO3/TNO3                   ! nitrate  partitioning = nitrate  / total nitrate    [%]

      NH4P = 0.
      IF(TNH4.GT.0.) NH4P = 100.*PNH4/TNH4                   ! ammonium partitioning = ammonium / total ammonium   [%]
!
! store aerosol species for diagnostic output:
!______________________________________________________________
! Input values:
      yo(il, 1) = TT   - 273.15                                       ! T                                    [degC]
      yo(il, 2) = RH   * 100.00                                       ! RH                                      [%]
      yo(il, 3) = TNH4                                                ! total input ammonium (g+p)       [umol/m^3]
      yo(il, 4) = TSO4                                                ! total input sulfate  (g+p)       [umol/m^3]
      yo(il, 5) = TNO3                                                ! total input nitrate  (g+p)       [umol/m^3]
      yo(il, 6) = TNa                                                 ! total input sodium   (p)         [umol/m^3]
      yo(il,33) = TCl                                                 ! total input chloride (g+p)       [umol/m^3]
      yo(il, 7) = TPo                                                 ! total input potasium (p)         [umol/m^3]
      yo(il,34) = TCa                                                 ! total input calcium  (p)         [umol/m^3]
      yo(il,35) = TMg                                                 ! total input magnesium(p)         [umol/m^3]
      yo(il,25) = PX                                                  ! atmospheric pressure                  [hPa]
! Output values:
      yo(il, 8) = GHCL                                                ! residual HCl   (g)               [umol/m^3]
      yo(il, 9) = GNO3                                                ! residual HNO3  (g)               [umol/m^3]
      yo(il,10) = GNH3                                                ! residual NH3   (g)               [umol/m^3]
      yo(il,11) = GSO4                                                ! residual H2SO4 (aq)              [umol/m^3]
      yo(il,12) = WH2O                                                ! aerosol Water  (aq)                [ug/m^3]
      yo(il,13) = PH                                                  ! aerosol pH                            [log]
      yo(il,14) = ZFLAG                                               ! concnetration domain [1=SP,2=SN,3=SR,4=SVR]
      yo(il,15) = PM                                                  ! total particulate matter         [umol/m^3]
      yo(il,16) = SNH4                                                ! solid ammonium (s)               [umol/m^3]
      yo(il,17) = SNO3                                                ! solid nitrate  (s)               [umol/m^3]
      yo(il,18) = SSO4                                                ! solid sulfate  (s)               [umol/m^3]
      yo(il,19) = PNH4                                                ! particulate ammonium (p=a+s)     [umol/m^3]
      yo(il,20) = PNO3                                                ! particulate nitrate  (p=a+s)     [umol/m^3]
      yo(il,21) = PSO4                                                ! particulate sulfate  (p=a+s)     [umol/m^3]
      yo(il,22) = RATIONS                                             ! mol ratio Nitrate/Sulfate (p)           [-]
      yo(il,23) = GAMA                                                ! activity coefficient (e.g. NH4NO3)           [-]
      yo(il,24) = ZIONIC                                              ! ionic strength (aq)                [mol/kg]
      yo(il,26) = PMt                                                 ! total PM (liquids & solids)        [ug/m^3]
      yo(il,27) = PMs                                                 ! total PM (solid)                   [ug/m^3]
      yo(il,28) = RINC                                                ! radius increase (H2O/PMt+1)**(1/3)      [-]
      yo(il,29) = SCl                                                 ! solid chloride (s)               [umol/m^3]
      yo(il,30) = SNa                                                 ! solid sodium (s)                 [umol/m^3]
      yo(il,31) = PCl                                                 ! particulate chloride (p=a+s)     [umol/m^3]
      yo(il,32) = PNa                                                 ! particulate sodium (p=a+s)       [umol/m^3]
      yo(il,36) = RHD                                                 ! RH of deliquescence
enddo
!
end subroutine eqsam_v03d

