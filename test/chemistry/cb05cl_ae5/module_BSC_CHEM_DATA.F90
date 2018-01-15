!-----------------------------------------------------------------------
!
                        module module_bsc_chem_data
!
!-----------------------------------------------------------------------
!
      implicit none
!      
      SAVE
!-----------------------------------------------------------------------
!***  DEFINITION OF CLOUD CHEMISTRY MAPPING INDEXES AND ARRAYS
!-----------------------------------------------------------------------
        INTEGER, PARAMETER :: N_GC_SPC = 53
        INTEGER, PARAMETER :: N_GC_SPC_CHEM = 72
        CHARACTER*16 :: GC_SPC(100)  !If more than 100 chem species increase
        REAL :: GC_MOLWT(100)

contains

  subroutine init_bsc_chem_data()        
!*******************************************************
! INITIALIZATION OF MAPPING ARRAYS FOR CLOUD CHEMISTRY *
!*******************************************************

! FIXME This is for MECH==2 or MECH==22 - Is this right for CB05?

      GC_SPC( 1) = 'NO2'
      GC_SPC( 2) = 'NO'
      GC_SPC( 3) = 'O3'
      GC_SPC( 4) = 'NO3'
      GC_SPC( 5) = 'N2O5'
      GC_SPC( 6) = 'HNO3'
      GC_SPC( 7) = 'HNO2' !'HONO'
      GC_SPC( 8) = 'HNO4' !'PNA'
      GC_SPC( 9) = 'H2O2'
      GC_SPC(10) = 'NTR'
      GC_SPC(11) = 'ROOH'
      GC_SPC(12) = 'FORMALDEHYDE' !'FORM'
      GC_SPC(13) = 'ACETALDEHYDE' !'ALD2'
      GC_SPC(14) = 'GENERIC_ALDEHYDE' !'ALDX'
      GC_SPC(15) = 'ETHANE' !'PAR'
      GC_SPC(16) = 'CO'
      GC_SPC(17) = 'METHYLHYDROPEROX' !'MEPX'
      GC_SPC(18) = 'METHANOL' !'MEOH'
      GC_SPC(19) = 'FORMIC_ACID' !'FACD'
      GC_SPC(20) = 'PAN'
      GC_SPC(21) = 'PEROXYACETIC_ACI' !'PACD'
      GC_SPC(22) = 'ACETIC_ACID' !'AACD'
      GC_SPC(23) = 'PPN' !'PANX'
      GC_SPC(24) = 'ETHENE' !'OLE'
      GC_SPC(25) = 'ETHENE' !'ETH'
      GC_SPC(26) = 'ETHENE' !'IOLE'
      GC_SPC(27) = 'TOLUENE' !'TOL'
      GC_SPC(28) = '2-CRESOL' !'CRES'
      GC_SPC(29) = 'OPEN'
      GC_SPC(30) = 'METHYL_GLYOXAL' !'MGLY'
      GC_SPC(31) = 'O-XYLENE' !'XYL'
      GC_SPC(32) = 'ISOPRENE' !'ISOP'
      GC_SPC(33) = 'ISPD'
      GC_SPC(34) = 'PINENE' !'TERP'
      GC_SPC(35) = 'SO2'
      GC_SPC(36) = 'H2SO4' !'SULF'
      GC_SPC(37) = 'ETHANOL' !'ETOH'
      GC_SPC(38) = 'ETHANE' !'ETHA'
!NON-ADVECTED SPECIES
      GC_SPC(39) = 'O'
      GC_SPC(40) = 'O1D'
      GC_SPC(41) = 'OH'
      GC_SPC(42) = 'HO2'
      GC_SPC(43) = 'XO2'
      GC_SPC(44) = 'XO2N'
      GC_SPC(45) = 'MEO2'
      GC_SPC(46) = 'HCO3'
      GC_SPC(47) = 'C2O3'
      GC_SPC(48) = 'CXO3'
      GC_SPC(49) = 'ROR'
      GC_SPC(50) = 'TO2'
      GC_SPC(51) = 'CRO'
      GC_SPC(52) = 'CH4'
!OJORBA2
      GC_SPC(53) = 'NO2E'
!OJORBA2

      GC_MOLWT( 1) = 46.0
      GC_MOLWT( 2) = 30.0
      GC_MOLWT( 3) = 48.0
      GC_MOLWT( 4) = 62.0
      GC_MOLWT( 5) = 108.0
      GC_MOLWT( 6) = 63.0
      GC_MOLWT( 7) = 47.0
      GC_MOLWT( 8) = 79.0
      GC_MOLWT( 9) = 34.0
      GC_MOLWT(10) = 130.0
      GC_MOLWT(11) = 62.0
      GC_MOLWT(12) = 30.0
      GC_MOLWT(13) = 44.0
      GC_MOLWT(14) = 44.0
      GC_MOLWT(15) = 14.0
      GC_MOLWT(16) = 28.0
      GC_MOLWT(17) = 48.0
      GC_MOLWT(18) = 32.0
      GC_MOLWT(19) = 46.0
      GC_MOLWT(20) = 121.0
      GC_MOLWT(21) = 76.0
      GC_MOLWT(22) = 60.0
      GC_MOLWT(23) = 121.0
      GC_MOLWT(24) = 27.0
      GC_MOLWT(25) = 28.0
      GC_MOLWT(26) = 48.0
      GC_MOLWT(27) = 92.0
      GC_MOLWT(28) = 108.0
      GC_MOLWT(29) = 84.0
      GC_MOLWT(30) = 72.0
      GC_MOLWT(31) = 106.0
      GC_MOLWT(32) = 68.0
      GC_MOLWT(33) = 70.0
      GC_MOLWT(34) = 136.0
      GC_MOLWT(35) = 64.0
      GC_MOLWT(36) = 98.0
      GC_MOLWT(37) = 46.0
      GC_MOLWT(38) = 30.0
!NON-ADVECTED SPECIES
      GC_MOLWT(39) = 16.0
      GC_MOLWT(40) = 16.0
      GC_MOLWT(41) = 17.0
      GC_MOLWT(42) = 33.0
      GC_MOLWT(43) = 1.0
      GC_MOLWT(44) = 1.0
      GC_MOLWT(45) = 47.0
      GC_MOLWT(46) = 63.0
      GC_MOLWT(47) = 75.0
      GC_MOLWT(48) = 75.0
      GC_MOLWT(49) = 31.0
      GC_MOLWT(50) = 109.0
      GC_MOLWT(51) = 139.0
      GC_MOLWT(52) = 18.0
!OJORBA2
      GC_MOLWT(53) = 46.0
!OJORBA2

  end subroutine init_bsc_chem_data

end module MODULE_BSC_CHEM_DATA
