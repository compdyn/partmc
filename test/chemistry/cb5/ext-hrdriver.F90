!OJORBA3
       MODULE MODULE_EBI_SOLVER_EXT
! DRIVER NMMB-EBI SOLVER CB05_EXT CMAQ
!
      USE MODULE_DM_PARALLEL
      USE MODULE_BSC_CHEM_DATA
!OJORBA3      USE GC_SPC_mod

      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE
!
!-----------------------------------------------------------------------
!
      CONTAINS
!OJORBA2
C***********************************************************************
C   Portions of Models-3/CMAQ software were developed or based on      *
C   information from various groups: Federal Government employees,     *
C   contractors working on a United States Government contract, and    *
C   non-Federal sources (including research institutions).  These      *
C   research institutions have given the Government permission to      *
C   use, prepare derivative works, and distribute copies of their      *
C   work in Models-3/CMAQ to the public and to permit others to do     *
C   so.  EPA therefore grants similar permissions for use of the       *
C   Models-3/CMAQ software, but users are requested to provide copies  *
C   of derivative works to the Government without restrictions as to   *
C   use by others.  Users are responsible for acquiring their own      *
C   copies of commercial software associated with Models-3/CMAQ and    *
C   for complying with vendor requirements.  Software copyrights by    *
C   the MCNC Environmental Modeling Center are used with their         *
C   permissions subject to the above restrictions.                     *
C***********************************************************************

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/work/rep/CCTM/src/chem/ebi_cb05cl_ae5/hrdriver.F,v 1.2 2008/08/05 16:01:39 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

!OJORBA2       SUBROUTINE CHEM( CGRID, JDATE, JTIME, TSTEP )
      SUBROUTINE RUN_EBI_EXT(NTIMESTEP,DT,NPHS,TEMPERATURE
     &          ,START_HOUR
     &          ,CHEM,NUM_GAS_TOTAL
     &          ,P8W,RHO_PHY,WATER,NUM_WATER,P_QV
     &          ,JPNL
     &          ,IDS,IDE,JDS,JDE,LM
     &          ,IMS,IME,JMS,JME,KMS,KME
     &          ,ITS,ITE,JTS,JTE)
!OJORBA2
C**********************************************************************
C
C  FUNCTION: Driver subroutine for Euler Backward Iterative solver
C
C  PRECONDITIONS: For the CB05CL family of mechanisms
C
C  KEY SUBROUTINES/FUNCTIONS CALLED:  HRINIT, PHOT, HRCALCKS, HRSOLVER
C
C  REVISION HISTORY: Created by EBI solver program, Jan. 31, 2008
C                       Based on the algorithm in "Test of Two Numerical
C                       Schemes for Use in Atmospheric Transport-Chemistry
C                       Models", O. Hertel, R. Berkowicz, J. Christensen,
C                       and O. Hov, Atm Env., Vol. 27A, No. 16, 1993.
C                       Original MEBI code developed by Ho-Chun Huang,
C                       SUNY, Albany -- "On the performance of numerical
C                       solvers for a chemistry submodel in three-dimensional
C                       air quality models 1. Box model simulations",
C                       H. Huang and J.S. Chang, JGR, Vol 106, No. D17, 2001.
C                       This version replaces Huang and Chang use of numerical
C                       solutions with analytical solutions derived in
C                       Hertel et al.
C
C**********************************************************************

!OJORBA2      USE HGRD_DEFN             ! horizontal domain specifications
!OJORBA2      USE VGRD_DEFN             ! vertical layer specifications
      USE EXT_HRDATA
!OJORBA2
      USE EXT_CONST
      USE EXT_RXCM
!OJORBA2
      IMPLICIT NONE

C..Includes:
!OJORBA2      INCLUDE SUBST_IOPARMS   ! Io/api parameters
!OJORBA2      INCLUDE SUBST_IOFDESC   ! Io/api file descriptions
!OJORBA2      INCLUDE SUBST_IODECL    ! Io/api declarations
!OJORBA2      INCLUDE SUBST_FILES_ID  ! CMAQ files
!OJORBA2      INCLUDE SUBST_CONST     ! CMAQ constants
!OJORBA2      INCLUDE SUBST_GC_SPC    ! Gas chem species names and MWs
!OJORBA2      INCLUDE SUBST_RXCMMN    ! Mechanism reaction common block
!OJORBA2      INCLUDE SUBST_GC_EMIS   ! Gas chem emissions name and mapping tables

!OJORBA2#ifdef emis_chem
!OJORBA2      INCLUDE SUBST_EMPR_CH   ! Emissions processing in chem
!OJORBA2#else
!OJORBA2      INCLUDE SUBST_EMPR_VD   ! Emissions processing in vdif
!OJORBA2#endif

!OJORBA2      INCLUDE SUBST_PACTL_ID  ! Process analysis control parameters
!OJORBA2
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM
     &                      ,IMS,IME,JMS,JME,KMS,KME
     &                      ,ITS,ITE,JTS,JTE
     &                      ,NPHS,NTIMESTEP,START_HOUR
     &                      ,NUM_WATER,P_QV
!
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE),INTENT(IN) :: JPNL
!
      REAL,INTENT(IN) :: DT

      REAL,DIMENSION(ITS:ITE,JTS:JTE,1:LM,1:NUM_WATER),INTENT(IN) ::
     & WATER

      INTEGER,INTENT(IN) :: NUM_GAS_TOTAL

      REAL,DIMENSION(ITS:ITE,JTS:JTE,1:LM,NUM_GAS_TOTAL),
     & INTENT(INOUT) :: CHEM

      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: TEMPERATURE

      REAL,DIMENSION(ITS:ITE,1:LM+1,JTS:JTE),INTENT(IN) :: P8W

      REAL,DIMENSION(ITS:ITE,1:LM,JTS:JTE),INTENT(IN) :: RHO_PHY
!
!-----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
!-----------------------------------------------------------------------
!OJORBA2
C..Arguments:
!OJORBA2      REAL, POINTER :: CGRID( :,:,:,: )  ! Species concentrations
      REAL,DIMENSION(ITS:ITE,JTS:JTE
     &    ,1:LM,1:72) :: CGRID    !!!OJORBA1 N_GC_SPCD=73 !
!OJORBA2
      INTEGER JDATE           ! Current date (YYYYDDD)
      INTEGER JTIME           ! Current time (HHMMSS)
      INTEGER TSTEP( 2 )      ! Time step vector (HHMMSS)

C..Parameters:
      REAL, PARAMETER :: PA2ATM = 1.0 / STDATMPA    ! Pascal to atm conv fac

C..External Functions:
!OJORBA2      INTEGER INDEX1          ! Looks up name in a list
!OJORBA2      INTEGER JUNIT           ! Gets logical device number
!OJORBA2      INTEGER SEC2TIME        ! Returns time interval from seconds
!OJORBA2      INTEGER TIME2SEC        ! Returns seconds in time interval

C..Saved Local Variables:

!OJORBA2      CHARACTER( 16 ), SAVE :: PNAME = 'HRDRIVER'     ! Program name

      INTEGER, SAVE :: ISTFL            ! Unit no. of iteration stat output file
      LOGICAL, SAVE :: LFIRST = .TRUE.  ! Flag for first call to this subroutine

      REAL, SAVE :: MAOMV               ! Mol Wt of air over Mol Wt of water

C..Scratch Local Variables:
      CHARACTER( 132 ) :: MSG           ! Message text
      CHARACTER(  16 ) :: VNAME         ! Name of I/O API data variable

      INTEGER C, E, L, R, S   ! Loop indices

      INTEGER AVGEBI          ! Average no. of EBI iterations
      INTEGER DELT_SEC        ! EBI max time step in seconds
      INTEGER ESP             ! Loop index for emissions species
      INTEGER ITMSTEP         ! Chemistry integration interval (sec)
      INTEGER LEV             ! Layer index
      INTEGER MIDDATE         ! Date at time step midpoint
      INTEGER MIDTIME         ! Time at time step midpoint
      INTEGER MNEBI           ! Min no. of EBI iterations
      INTEGER MXEBI           ! Max no. of EBI iterations
      INTEGER NDARK           ! Number of layer 1 cells in darkness
      INTEGER NPH             ! Index for number of phot. rxns in PHOT
      INTEGER SPC             ! Species loop index
      INTEGER STATUS          ! Status code
      INTEGER VAR             ! Variable number on I/O API file

      LOGICAL LSUNLIGHT       ! Flag for sunlight

      REAL ATMPRES            ! Cell pressure
      REAL CHEMSTEP           ! Chemistry integration interval (min)
      REAL H2O                ! Cell H2O mixing ratio (ppmV)
      REAL SUMEBI             ! Sum of EBI iterations
      REAL TEMP               ! Cell Temperature

!OJORBA2      REAL PRES(    NCOLS, NROWS, NLAYS )        ! Cell pressure (Pa)
!OJORBA2      REAL QV(      NCOLS, NROWS, NLAYS )        ! Cell water vapor (Kg/Kg air)
!OJORBA2      REAL TA(      NCOLS, NROWS, NLAYS )        ! Cell temperature (K)
!OJORBA2      REAL RJIN( NPHOTAB )                       ! J-values for a cell
!OJORBA2      REAL RJ( NCOLS, NROWS, NLAYS, NPHOTAB )    ! J-values for each cell
      REAL PRES(ITS:ITE,JTS:JTE,1:LM)        ! Cell pressure (Pa)
      REAL QV(ITS:ITE,JTS:JTE,1:LM)  ! Cell water vapor (Kg/Kgair)
      REAL TA(ITS:ITE,JTS:JTE,1:LM)        ! Cell temperature (K)
      REAL RJIN( NPHOTAB )                       ! J-values for a cell
      REAL RJ(ITS:ITE,JTS:JTE,1:LM, NPHOTAB ) ! J-values for each cell

      INTEGER KFLIP
!OJORBA2

      INTEGER     GXOFF, GYOFF            ! global origin offset from file
C for INTERPX
      INTEGER, SAVE :: STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3


C**********************************************************************

!OJORBA2      IF( N_GC_SPC .EQ. 0 ) RETURN
      ITMSTEP = NPHS * DT
!OJORBA2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  On first call, call routines to set-up for EBI solver
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IF( LFIRST ) THEN

!OJORBA2         LOGDEV = INIT3( )
         LOGDEV = 6
!OJORBA2
!OJORBA2         IF( MECHNAME .NE. 'CB05CL_AE5_AQ' ) THEN
!OJORBA2             MSG = 'This version of the EBI solver can only be used with'
!OJORBA2     &            // ' the CB05CL_AE5_AQ chemical mechanism'
!OJORBA2             CALL M3EXIT( PNAME, 0, 0, MSG, XSTAT1 )
!OJORBA2         ENDIF

!OJORBA2         IF( INDEX( MECHNAME, 'AE' ) .NE. 0 ) THEN
!OJORBA2           L_AE_VRSN = .TRUE.
!OJORBA2         ELSE
           L_AE_VRSN = .FALSE.
!OJORBA2         ENDIF

!OJORBA2         IF( INDEX( MECHNAME, 'AQ' ) .NE. 0 ) THEN
!OJORBA2           L_AQ_VRSN = .TRUE.
!OJORBA2         ELSE
           L_AQ_VRSN = .FALSE.
!OJORBA2         ENDIF

!OJORBA2         IF( LIRR ) THEN
!OJORBA2            MSG = 'IRR Analysis not allowed with EBI solver'
!OJORBA2            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
!OJORBA2         ENDIF
!OJORBA3
         CALL EXT_HRINIT
!OJORBA3
!OJORBA2         ITMSTEP = TIME2SEC( TSTEP( 2 ) )
         CHEMSTEP = FLOAT( ITMSTEP ) / 60.0
!OJORBA2         WRITE( LOGDEV, 92000 ) CHEMSTEP, DELTAT

!OJORBA2         WRITE( LOGDEV, 92020 )
!OJORBA2         DO SPC = 1, N_GC_SPC
!OJORBA2            WRITE( LOGDEV, 92040 ) GC_SPC( SPC ), RTOL( SPC )
!OJORBA2         ENDDO

         MAOMV =  MWAIR / MWWAT

c..If emissions processing requested stop
!OJORBA2         IF( EMISCH ) THEN

!OJORBA2            MSG = 'ERROR: EBI solver not configured to '//
!OJORBA2     &            'process emissions in chemistry'
!OJORBA2            CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )

!OJORBA2         ENDIF   ! End if doing emissions


!OJORBA2#ifdef hrstats
!OJORBA2         ISTFL = JUNIT()
!OJORBA2         OPEN( UNIT=ISTFL, FILE='iterstat.dat' )
!OJORBA2         WRITE( ISTFL, 94020 )
!OJORBA2#endif

!OJORBA2         CALL SUBHFILE ( MET_CRO_3D, GXOFF, GYOFF,
!OJORBA2     &                   STRTCOLMC3, ENDCOLMC3, STRTROWMC3, ENDROWMC3 )

         LFIRST = .FALSE.

      ENDIF      ! First time
!OJORBA2
      RKI(:)=0.
      RXRAT(:)=0.
      YC(:)=0.
      YC0(:)=0.
      YCP(:)=0.
      PROD(:)=0.
      LOSS(:)=0.
      PNEG(:)=0.

      DO C = 1,N_GC_SPC_CHEM  !=72 LOOP OVER SPECIES
        CGRID(:,:,:,C) = CHEM(:,:,:,EBI_SPC(C)) !MAX(CHEM(:,:,:,EBI_SPC(C)),1.0E-30)
      ENDDO
!OJORBA2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C  Set date and time to center of time step, get necessary physical
C  data, and get photolysis rates
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!OJORBA2      MIDDATE = JDATE
!OJORBA2      MIDTIME = JTIME
!OJORBA2      ITMSTEP = TIME2SEC( TSTEP( 2 ) )
      CHEMSTEP = FLOAT( ITMSTEP ) / 60.0D+00
!OJORBA2      CALL NEXTIME( MIDDATE, MIDTIME, SEC2TIME( ITMSTEP / 2 ) )

C.. Compute number of time step loops and step size for EBI solver
      DELT_SEC = DELTAT * 60.0 + 0.1
      IF( DELT_SEC .GE. ITMSTEP ) THEN
         N_EBI_STEPS = 1
         EBI_TMSTEP = FLOAT( ITMSTEP ) / 60.0
      ELSE
         IF( MOD( ITMSTEP, DELT_SEC ) .EQ. 0 ) THEN
            N_EBI_STEPS = ITMSTEP / DELT_SEC
         ELSE
            N_EBI_STEPS = ITMSTEP / DELT_SEC + 1
         ENDIF
         EBI_TMSTEP =  FLOAT( ITMSTEP ) / FLOAT( N_EBI_STEPS ) / 60.0
      ENDIF

      N_INR_STEPS = 1

C.. Get ambient temperature in K

!OJORBA2      VNAME = 'TA'
!OJORBA2      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
!OJORBA2     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
!OJORBA2     &                    MIDDATE, MIDTIME, TA ) ) THEN
!OJORBA2         MSG = 'Could not read TA from MET_CRO_3D'
!OJORBA2         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
!OJORBA2      ENDIF

C.. Get specific humidity in Kg H2O / Kg air
!OJORBA2      VNAME = 'QV'
!OJORBA2      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
!OJORBA2     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
!OJORBA2     &                    MIDDATE, MIDTIME, QV ) ) THEN
!OJORBA2         MSG = 'Could not read QV from MET_CRO_3D'
!OJORBA2         CALL M3EXIT( PNAME, JDATE, JTIME, MSG, XSTAT1 )
!OJORBA2      ENDIF

C.. Get pressure in Pascals
!OJORBA2      VNAME = 'PRES'
!OJORBA2      IF ( .NOT. INTERPX( MET_CRO_3D, VNAME, PNAME,
!OJORBA2     &                    STRTCOLMC3,ENDCOLMC3, STRTROWMC3,ENDROWMC3, 1,NLAYS,
!OJORBA2     &                    MIDDATE, MIDTIME, PRES ) ) THEN
!OJORBA2         MSG = 'Could not read PRES from MET_CRO_3D'
!OJORBA2         CALL M3EXIT ( PNAME, JDATE, JTIME, MSG, XSTAT1 )
!OJORBA2      ENDIF

C.. Get photolysis rates in /min
!OJORBA2      CALL PHOT ( MIDDATE, MIDTIME, JDATE, JTIME, NDARK, RJ )

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Top of loop over cells
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!OJORBA2#ifdef hrstats
!OJORBA2      MNEBI = 1000
!OJORBA2      MXEBI = 0
!OJORBA2      SUMEBI = 0.0
!OJORBA2#endif

!OJORBA2      DO L = 1, NLAYS
!OJORBA2         DO R = 1, MY_NROWS
!OJORBA2            DO C = 1, MY_NCOLS
!OJORBA13      DO L = 1,LM-2 !OJORBA12
      DO L = 1,JPNL(ITS,JTS) 
      DO R = JTS,JTE !OJORBA1 1, MY_NROWS
         DO C = ITS,ITE !OJORBA1 1, MY_NCOLS
!OJORBA12            DO L = 1,JPNL(C,R) !OJORBA1 1, NLAYS
!OJORBA12        KFLIP=JPNL(C,R)+1-L+(LM-JPNL(C,R)) ! The flipped index will thus
            KFLIP=LM+1-L
!       be for the NMMB arrays
!OJORBA2

c..Load ICs
               DO S = 1, N_GC_SPC_CHEM
!OJORBA2                  YC( S ) = MAX( CGRID( C, R, L, S ), 1.0E-30 )
                  YC( S ) = MAX( CGRID( C, R, KFLIP, S ), 1.0E-30 )
!OJORBA1
               ENDDO


c..Set physical quantities
!OJORBA2               TEMP = TA( C, R, L )
!OJORBA2               ATMPRES = PA2ATM * PRES( C, R, L )
!OJORBA2               H2O  = MAX ( QV( C, R, L ) * MAOMV *  1.0E+06, 0.0 )
               TEMP = TEMPERATURE(C,R,KFLIP) !TA( C, R, L )
               ATMPRES = p8w(C,L,R) * PA2ATM !P8W in Pa
               H2O  = MAX(WATER(C,R,kflip,P_QV) * MAOMV * 1.0e+06,0.0)
!OJORBA2
c..Get rate constants
               LSUNLIGHT = .FALSE.
!OJORBA2
! NEED TO MATCH RJ NPHOTAB WITH ph_... vars
!OJORBA2               DO NPH = 1, NPHOTAB
!OJORBA2                  RJIN( NPH ) = RJ( C, R, L, NPH )
!OJORBA2                  IF( RJ( C, R, L, NPH ) .GT. 0.0 ) LSUNLIGHT = .TRUE.
!OJORBA2               ENDDO
        RJIN(:)    = 0.
        RJIN(1)    = ph_no2(C,L,R)                  !NO2_SAPRC99
        RJIN(2)    = ph_o33p(C,L,R)                 !O3_O3P_IUPAC04
        RJIN(3)    = ph_o31d(C,L,R)                 !O3_O1D_IUPAC04
        RJIN(4)    = ph_no3o(C,L,R)                 !NO3NO2_SAPRC99
        RJIN(5)    = ph_no3o2(C,L,R)                !NO3NO_SAPRC99
        RJIN(6)    = ph_hno2(C,L,R)                 !HONO_IUPAC04
        RJIN(7)    = ph_h2o2(C,L,R)
        RJIN(8)    = ph_hno4(C,L,R) !HNO4=HO2NO2
        RJIN(9)    = ph_hno3(C,L,R)
        RJIN(10)   = ph_n2o5(C,L,R)
        RJIN(11)   = ph_ntr(C,L,R) !!!!!
        RJIN(12)   = ph_mepx(C,L,R) !ROOH
!                   = ph_mepx(C,L,R) !MEPX
        RJIN(13)   = ph_ch2or(C,L,R)
        RJIN(14)   = ph_ch2om(C,L,R)
        RJIN(15)   = ph_ald2(C,L,R) !!!CCHO_R_SAPRC99  ????
        RJIN(16)   = ph_pan(C,L,R) !!!!!
        RJIN(17)   = ph_pacd(C,L,R) !!!!!
        RJIN(18)   = ph_aldx(C,L,R) !!!!!C2CHO_SAPRC99
!                     ph_pan(C,L,R) !PANX
        RJIN(19)   = ph_mgly(C,L,R) !!!!!
        RJIN(20)   = ph_ispd(C,L,R) !!!!!
        RJIN(21)   = 0.             !  <CL2_IUPAC04>
        RJIN(22)   = 0.             !  <HOCL_IUPAC04>
        RJIN(23)   = 0.             !  <FMCL_IUPAC04>

       DO NPH=1,23
        IF (RJIN(NPH) .GT. 0.0) LSUNLIGHT = .TRUE.
       ENDDO
!OJORBA2
!OJORBA3
               CALL EXT_HRCALCKS( NPHOTAB, LSUNLIGHT, RJIN, TEMP,
     &                        ATMPRES, H2O, RKI )
!OJORBA3

c..Call EBI solver
               N_EBI_IT = 0
!OJORBA3
               CALL EXT_HRSOLVER( JDATE, JTIME, C, R, L )
!OJORBA3
!OJORBA2#ifdef hrstats
!OJORBA2               MXEBI  = MAX( MXEBI, N_EBI_IT )
!OJORBA2               MNEBI  = MIN( MNEBI, N_EBI_IT )
!OJORBA2               SUMEBI = SUMEBI + FLOAT( N_EBI_IT )
!OJORBA2#endif


c..Update concentration array
               DO S = 1, N_GC_SPC_CHEM
!OJORBA2                  CGRID( C, R, L, S ) = YC( S )
                  CGRID( C, R, KFLIP, S ) = YC( S )
!OJORBA2
               ENDDO

            ENDDO
         ENDDO
      ENDDO

!OJORBA2
      DO C = 1,N_GC_SPC_CHEM  !=72 LOOP OVER SPECIES
        CHEM(:,:,:,EBI_SPC(C)) = CGRID(:,:,:,C) 
      ENDDO
!OJORBA2

!OJORBA2#ifdef hrstats
!OJORBA2      AVGEBI = SUMEBI / FLOAT( NCOLS * NROWS * NLAYS )
!OJORBA2      WRITE( ISTFL, 94040 ) JDATE, JTIME, MNEBI, AVGEBI, MXEBI
!OJORBA2#endif

      RETURN

C*********************** FORMAT STATEMENTS ****************************

92000 FORMAT( / 10X, 'Euler Backward Iterative Parameters -'
     &        / 10X, 'Chemistry Integration Time Interval (min):', F12.4,
     &        / 10X, 'EBI maximum time step (min):              ', F12.4 )

92020 FORMAT( //10X, 'Species convergence tolerances:' )

92040 FORMAT(   10X, A16, 2X, 1PE12.2 )

92060 FORMAT( / 10X, 'Emissions Processing in Chemistry ...'
     &        / 10X, 'Number of Emissions Layers:         ', I3
     &        / 10X, 'out of total Number of Model Layers:', I3 )


94020 FORMAT( 'DATE      TIME ', 'MNEBI AVEBI MXEBI' )

94040 FORMAT( I7, 1X, I6, 1X, 3( I5, 1X ) )
!OJORBA2      END
      END SUBROUTINE RUN_EBI_EXT

      END MODULE MODULE_EBI_SOLVER_EXT
!OJORBA2
