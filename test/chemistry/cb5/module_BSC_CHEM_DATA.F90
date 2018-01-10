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
        INTEGER, PARAMETER :: N_GC_SPC = 32
        INTEGER, PARAMETER :: N_GC_SPC_CHEM = 72
        CHARACTER*16 :: GC_SPC(100)  !If more than 100 chem species increase
        REAL :: GC_MOLWT(100)

contains

  subroutine init_bsc_chem_data()        
!*******************************************************
! INITIALIZATION OF MAPPING ARRAYS FOR CLOUD CHEMISTRY *
!*******************************************************

        GC_SPC(  1)= 'O1D             '
        GC_SPC(  2)= 'H2O2            '
        GC_SPC(  3)= 'PAN             '
        GC_SPC(  4)= 'CRO             '
        GC_SPC(  5)= 'TOL             '
        GC_SPC(  6)= 'N2O5            '
        GC_SPC(  7)= 'XYL             '
        GC_SPC(  8)= 'XO2N            '
        GC_SPC(  9)= 'HONO            '
        GC_SPC( 10)= 'PNA             '
        GC_SPC( 11)= 'TO2             '
        GC_SPC( 12)= 'HNO3            '
        GC_SPC( 13)= 'ROR             '
        GC_SPC( 14)= 'CRES            '
        GC_SPC( 15)= 'MGLY            '
        GC_SPC( 16)= 'CO              '
        GC_SPC( 17)= 'ETH             '
        GC_SPC( 18)= 'XO2             '
        GC_SPC( 19)= 'OPEN            '
        GC_SPC( 20)= 'PAR             '
        GC_SPC( 21)= 'FORM            '
        GC_SPC( 22)= 'ISOP            '
        GC_SPC( 23)= 'OLE             '
        GC_SPC( 24)= 'ALD2            '
        GC_SPC( 25)= 'O3              '
        GC_SPC( 26)= 'NO2             '
        GC_SPC( 27)= 'OH              '
        GC_SPC( 28)= 'HO2             '
        GC_SPC( 29)= 'O               '
        GC_SPC( 30)= 'NO3             '
        GC_SPC( 31)= 'NO              '
        GC_SPC( 32)= 'C2O3            '

        GC_MOLWT(  1) = 16.0 
        GC_MOLWT(  2) = 34.0 
        GC_MOLWT(  3) = 121.0 
        GC_MOLWT(  4) = 139.0 
        GC_MOLWT(  5) = 92.0 
        GC_MOLWT(  6) = 108.0 
        GC_MOLWT(  7) = 106.0 
        GC_MOLWT(  8) = 1.0 
        GC_MOLWT(  9) = 47.0 
        GC_MOLWT( 10) = 79.0 
        GC_MOLWT( 11) = 109.0 
        GC_MOLWT( 12) = 63.0 
        GC_MOLWT( 13) = 31.0 
        GC_MOLWT( 14) = 108.0 
        GC_MOLWT( 15) = 72.0 
        GC_MOLWT( 16) = 28.0 
        GC_MOLWT( 17) = 28.0 
        GC_MOLWT( 18) = 1.0 
        GC_MOLWT( 19) = 84.0 
        GC_MOLWT( 20) = 14.0 
        GC_MOLWT( 21) = 30.0 
        GC_MOLWT( 22) = 68.0 
        GC_MOLWT( 23) = 27.0 
        GC_MOLWT( 24) = 44.0 
        GC_MOLWT( 25) = 48.0 
        GC_MOLWT( 26) = 46.0 
        GC_MOLWT( 27) = 17.0 
        GC_MOLWT( 28) = 33.0 
        GC_MOLWT( 29) = 16.0 
        GC_MOLWT( 30) = 62.0 
        GC_MOLWT( 31) = 30.0 
        GC_MOLWT( 32) = 75.0 

  end subroutine init_bsc_chem_data

end module MODULE_BSC_CHEM_DATA
