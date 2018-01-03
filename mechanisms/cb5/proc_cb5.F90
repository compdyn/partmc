! Copyright (C) 2017 Matthew Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_proc_cb5 program

!> A tool to generate CB-5 input files for phlex-chem from the ebi-cb05
!! solver in the MONARCH model.
program pmc_proc_cb5

  use pmc_constants,                only : i_kind, dp
  use pmc_util,                     only : die_msg, warn_msg, to_string, string_t

#ifdef PMC_USE_JSON
  use json_module
#endif

  implicit none

  type(string_t), allocatable :: spec_name(:)


  ! Set up the species name list
  allocate(  spec_name(71)%string

  spec_name(1)%string = "NO2"
  spec_name(2)%string = "NO"
  spec_name(3)%string = "O"
  spec_name(4)%string = "O3"
  spec_name(5)%string = "NO3"
  spec_name(6)%string = "O1D"
  spec_name(7)%string = "OH"
  spec_name(8)%string = "HO2"
  spec_name(9)%string = "N2O5"
  spec_name(10)%string = "HNO3"
  spec_name(11)%string = "HONO"
  spec_name(12)%string = "PNA"
  spec_name(13)%string = "H2O2"
  spec_name(14)%string = "XO2"
  spec_name(15)%string = "XO2N"
  spec_name(16)%string = "NTR"
  spec_name(17)%string = "ROOH"
  spec_name(18)%string = "FORM"
  spec_name(19)%string = "ALD2"
  spec_name(20)%string = "ALDX"
  spec_name(21)%string = "PAR"
  spec_name(22)%string = "CO"
  spec_name(23)%string = "MEO2"
  spec_name(24)%string = "MEPX"
  spec_name(25)%string = "MEOH"
  spec_name(26)%string = "HCO3"
  spec_name(27)%string = "FACD"
  spec_name(28)%string = "C2O3"
  spec_name(29)%string = "PAN"
  spec_name(30)%string = "PACD"
  spec_name(31)%string = "AACD"
  spec_name(32)%string = "CXO3"
  spec_name(33)%string = "PANX"
  spec_name(34)%string = "ROR"
  spec_name(35)%string = "OLE"
  spec_name(36)%string = "ETH"
  spec_name(37)%string = "IOLE"
  spec_name(38)%string = "TOL"
  spec_name(39)%string = "CRES"
  spec_name(40)%string = "TO2"
  spec_name(41)%string = "TOLRO2"
  spec_name(42)%string = "OPEN"
  spec_name(43)%string = "CRO"
  spec_name(44)%string = "MGLY"
  spec_name(45)%string = "XYL"
  spec_name(46)%string = "XYLRO2"
  spec_name(47)%string = "ISOP"
  spec_name(48)%string = "ISPD"
  spec_name(49)%string = "ISOPRXN"
  spec_name(50)%string = "TERP"
  spec_name(51)%string = "TRPRXN"
  spec_name(52)%string = "SO2"
  spec_name(53)%string = "SULF"
  spec_name(54)%string = "SULRXN"
  spec_name(55)%string = "ETOH"
  spec_name(56)%string = "ETHA"
  spec_name(57)%string = "CL2"
  spec_name(58)%string = "CL"
  spec_name(59)%string = "HOCL"
  spec_name(60)%string = "CLO"
  spec_name(61)%string = "FMCL"
  spec_name(62)%string = "HCL"
  spec_name(63)%string = "TOLNRXN"
  spec_name(64)%string = "TOLHRXN"
  spec_name(65)%string = "XYLNRXN"
  spec_name(66)%string = "XYLHRXN"
  spec_name(67)%string = "BENZENE"
  spec_name(68)%string = "BENZRO2"
  spec_name(69)%string = "BNZNRXN"
  spec_name(70)%string = "BNZHRXN"
  spec_name(71)%string = "SESQ"

  

end program pmc_proc_cb5
