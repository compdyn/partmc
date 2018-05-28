!===============================================================================
module TRACERS_SOA
!@sum Module for calculating the secondary organic aerosol formation, based on
!@+ Tsigaridis and Kanakidou (2007) and references therein.
!@auth Kostas Tsigaridis (kostas.tsigaridis@nasa.gov)
!===============================================================================
!===============================================================================
implicit none
private
public :: soa_bins,soa_tracers,soa_conditions,soa_init,soa_calc, &
          isoa_apin,isoa_isop,isoa_oh,isoa_o3,isoa_no3,&
          isoa_isopp1,isoa_isopp2,isoa_apinp1,isoa_apinp2
!-------------------------------------------------------------------------------
!@param soa_precursors number of soa precursors
integer, parameter :: soa_precursors=2
!@param isoa_apin,isoa_isop indices of soa precursors
integer, parameter :: isoa_apin=1,isoa_isop=2
!@param soa_oxidants number of oxidants that react with soa precursors
integer, parameter :: soa_oxidants=3
!@param isoa_oh,isoa_o3,isoa_no3 indices of oxidants
integer, parameter :: isoa_oh=1,isoa_o3=2,isoa_no3=3
!@param soa_products number of soa products
integer, parameter :: soa_products=2
!@param isoa_isopp1,isoa_isopp2,isoa_apinp1,isoa_apinp2 indices of soa tracers
integer, parameter :: isoa_isopp1=1,isoa_isopp2=2, &
                      isoa_apinp1=3,isoa_apinp2=4
!@param soa_bins total number of semi-volatile soa compounds,
!@+ =soa_precursors*soa_products
integer, parameter :: soa_bins=soa_precursors*soa_products
!-------------------------------------------------------------------------------
type soa_tracers
!@var voc gas-phase concentration of precursor tracer (molecules cm-3)
real*8, dimension(soa_precursors) :: voc
!@var oxidant gas-phase concentration of oxidant tracer (molecules cm-3)
real*8, dimension(soa_oxidants) :: oxidant
!@var gas gas-phase concentration of tracer (ug m-3)
!@var aer aerosol-phase concentration of tracer (ug m-3)
real*8, dimension(soa_bins) :: gas, aer
!@var ivoc index of voc tracers of the host model
!@var ioxidant index of oxidant tracers of the host model
!@var igas index of gas-phase soa tracers of the host model
!@var iaer index of aerosol-phase soa tracers of the host model
integer, dimension(soa_precursors) :: ivoc
integer, dimension(soa_oxidants) :: ioxidant
integer, dimension(soa_bins) :: igas, iaer
!@var ivocinv inverse ivoc
!@var ioxidantinv inverse ioxidant
!@var igasinv inverse igas
!@var iaerinv inverse iaer
integer, allocatable, dimension(:) :: ivocinv, ioxidantinv, igasinv, iaerinv
!@var chem_prod amount of chemically produced tracer in the gas-phase
!@var partition amount of tracer that partitioned TO the aerosol phase (if
!@+ positive, evaporated FROM the aerosol phase, if negative
real*8, dimension(soa_bins) :: chem_prod, partition
end type soa_tracers
!-------------------------------------------------------------------------------
type soa_properties
!@var amass_lownox Reference mass-based partitioning yields for low NOx
!@var amass_highnox Reference mass-based partitioning yields for high NOx
!@var amass Mass-based partitioning yields for current NOx conditions
!@var sat saturation concentration (C*) of the soa bins (ug m-3). The lower
!@+ the saturation concentration, the less volatile the bin.
!@var dH enthalpy of vaporization of the soa bins (KJ mol-1). The higher
!@+ the enthalpy of vaporization, the strongest the dependency of saturation
!@+ concentration to temperature.
!@var Tref reference temperature for which sat corresponds to
real*8, dimension(soa_bins) :: amass_lownox,amass_highnox,amass,sat,dH,Tref
!@var soa index of parent voc
integer, dimension(soa_bins) :: parent
end type soa_properties
!-------------------------------------------------------------------------------
type soa_conditions
!@var dt timestep (seconds)
!@var voc2nox Ratio of ((XO2+NO) + (XO2+HO2) + (XO2+XO2))/(XO2+NO), for
!@+ determining whether conditions are high or low NOx (dimensionless)
!@var temp temperature (K)
!@var nvoa concentration of non-volatile organic aerosol (ug m-3)
real*8 :: dt, voc2nox, temp, nvoa
!@var rr reaction rates (m3 molecules-1 s-1)
real*8, dimension(soa_precursors, soa_oxidants) :: rr
end type soa_conditions
!-------------------------------------------------------------------------------
!@var soa_tr soa precursor and oxidant concentrations (molecules cm-3),
!@+ soa tracers concentrations (ug m-3), indices and diagnostics
!@+ (ug m-3 timestep-1)
type(soa_tracers), public :: soa_tr
!@var soa_prop soa bin properties and parent voc maps
type(soa_properties) :: soa_prop
!@var soa_cond current contitions of the atmosphere
type(soa_conditions), public :: soa_cond
!-------------------------------------------------------------------------------
contains
!===============================================================================
!===============================================================================
subroutine soa_init(ntm_host)
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
!@var ntm_host number of tracers in the host model
integer, intent(in) :: ntm_host
integer :: i
!-------------------------------------------------------------------------------
allocate(soa_tr%ivocinv(ntm_host))
allocate(soa_tr%ioxidantinv(ntm_host))
allocate(soa_tr%igasinv(ntm_host))
allocate(soa_tr%iaerinv(ntm_host))
soa_tr%ivocinv=0
soa_tr%ioxidantinv=0
soa_tr%igasinv=0
soa_tr%iaerinv=0
do i=1,soa_precursors
  soa_tr%ivocinv(soa_tr%ivoc(i))=i
enddo
do i=1,soa_oxidants
  soa_tr%ioxidantinv(soa_tr%ioxidant(i))=i
enddo
do i=1,soa_bins
  soa_tr%igasinv(soa_tr%igas(i))=i
  soa_tr%iaerinv(soa_tr%iaer(i))=i
enddo
soa_prop%parent(isoa_apinp1)=isoa_apin
soa_prop%parent(isoa_apinp2)=isoa_apin
soa_prop%parent(isoa_isopp1)=isoa_isop
soa_prop%parent(isoa_isopp2)=isoa_isop
soa_prop%amass_lownox(isoa_apinp1)=0.192d0
soa_prop%amass_lownox(isoa_apinp2)=0.215d0
soa_prop%amass_lownox(isoa_isopp1)=0.0288d0
soa_prop%amass_lownox(isoa_isopp2)=0.232d0
soa_prop%amass_highnox(isoa_apinp1)=0.0138d0
soa_prop%amass_highnox(isoa_apinp2)=0.461d0
soa_prop%amass_highnox(isoa_isopp1)=soa_prop%amass_lownox(isoa_isopp1)*&
                                    soa_prop%amass_highnox(isoa_apinp1)/&
                                    soa_prop%amass_lownox(isoa_apinp1)
soa_prop%amass_highnox(isoa_isopp2)=soa_prop%amass_lownox(isoa_isopp2)*&
                                    soa_prop%amass_highnox(isoa_apinp2)/&
                                    soa_prop%amass_lownox(isoa_apinp2)
soa_prop%sat(isoa_apinp1)=15.7d0 ! Presto et al., 2005
soa_prop%sat(isoa_apinp2)=385.d0 ! Presto et al., 2005
soa_prop%sat(isoa_isopp1)=1.d0/1.62d0 ! Henze and Seinfeld, 2006
soa_prop%sat(isoa_isopp2)=1.d0/0.00862d0 ! Henze and Seinfeld, 2006
soa_prop%dH(isoa_isopp1)=42.d0 ! Chung and Seinfeld, JGR, 2002; Henze and Seinfeld, GRL, 2006; Kleindienst et al., GRL, 2007
soa_prop%dH(isoa_isopp2)=42.d0 ! Chung and Seinfeld, JGR, 2002; Henze and Seinfeld, GRL, 2006; Kleindienst et al., GRL, 2007
soa_prop%dH(isoa_apinp1)=72.9d0
soa_prop%dH(isoa_apinp2)=72.9d0
soa_prop%Tref(isoa_apinp1)=295.d0 ! Presto et al., 2005
soa_prop%Tref(isoa_apinp2)=295.d0 ! Presto et al., 2005
soa_prop%Tref(isoa_isopp1)=295.d0 ! Henze and Seinfeld, 2006
soa_prop%Tref(isoa_isopp2)=295.d0 ! Henze and Seinfeld, 2006
!===============================================================================
end subroutine soa_init
!===============================================================================
!===============================================================================
subroutine soa_calc
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
soa_prop%amass=soa_cond%voc2nox*soa_prop%amass_highnox+&
               (1.d0-soa_cond%voc2nox)*soa_prop%amass_lownox
call soa_gas
call soa_aer
!===============================================================================
end subroutine soa_calc
!===============================================================================
!===============================================================================
subroutine soa_gas
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
real*8, dimension(soa_bins) :: tr0_gas
real*8, dimension(soa_precursors) :: prod, dest
integer :: i, j
!-------------------------------------------------------------------------------
tr0_gas=soa_tr%gas
prod=0.d0
dest=0.d0
do i=1,soa_precursors
  do j=1,soa_oxidants ! accumulated over all oxidants, assuming same products
    prod(i)=prod(i)+soa_cond%rr(i,j)*soa_tr%voc(i)*soa_tr%oxidant(j)
  enddo
enddo
do i=1,soa_bins
  soa_tr%gas(i)=(tr0_gas(i)+&
                soa_prop%amass(i)*prod(soa_prop%parent(i))*soa_cond%dt)/&
                (1.d0+dest(soa_prop%parent(i))*soa_cond%dt)
enddo
!===============================================================================
end subroutine soa_gas
!===============================================================================
!===============================================================================
subroutine soa_aer
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!CONTINUE FROM HERE: Calculate Kp before doing the partitioning, based on composition changes
call soa_partition
!===============================================================================
end subroutine soa_aer
!===============================================================================
!===============================================================================
subroutine soa_partition
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
integer :: iter ! iteration count
real*8, dimension(soa_bins) :: ksi,soa_tot ! see Donahue et al., 2006
!@var Ccurr value of C* at current temperature (ug m-3)
real*8, dimension(soa_bins) :: Ccurr
!@param maxit maximum number of iterations allowed
integer, parameter :: maxit=10000 !!!!!!!!! this is EXTREMELY HIGH!!!!!!!!!
!@var Mo total mass of condensed OA at equilibrium (ug m-3)
!@var Mo_guess total mass of condensed OA while searching for Mo (ug m-3)
real*8 :: Mo, Mo_guess
!-------------------------------------------------------------------------------
soa_tot=soa_tr%gas+soa_tr%aer
if (sum(soa_tot)+soa_cond%nvoa == 0.d0) return
call clausius_clapeyron(Ccurr)
Mo=soa_cond%nvoa
Mo_guess=max(tiny(Mo),Mo)
iter=0
do ! loop indefinitely until a solution is found, or iter > maxit
  iter=iter+1
  ksi=1.d0/(1.d0+Ccurr/Mo_guess)
  Mo=soa_cond%nvoa+sum(soa_tot*ksi)
  if (Mo == Mo_guess) then
!    print *,'Solution found after ',iter,' iterations'
!    print *,'ksi=',ksi
!    print *,'temp=',soa_cond%temp
!    print *,'Ccurr=',Ccurr
    exit
  else
    if (iter > maxit) then
      print *,'iter>maxit, breaking operation. Mo= ',Mo,' Mo_guess= ',Mo_guess
      exit
    endif
  endif
  Mo_guess=Mo
enddo
soa_tr%partition=soa_tot*ksi-soa_tr%aer ! positive means TO the aerosol phase
soa_tr%gas=soa_tot*(1.d0-ksi)
soa_tr%aer=soa_tot*ksi
!===============================================================================
end subroutine soa_partition
!===============================================================================
!===============================================================================
subroutine clausius_clapeyron(C)
!===============================================================================
use modelE_Global, only: gasc
implicit none
!-------------------------------------------------------------------------------
real*8, dimension(soa_bins) :: C
!-------------------------------------------------------------------------------
C=soa_prop%sat*soa_prop%Tref/soa_cond%temp*& ! 1.d3 converts KJ to J
  exp(1.d3*soa_prop%dH/gasc*(1.d0/soa_prop%Tref-1.d0/soa_cond%temp))
!===============================================================================
end subroutine clausius_clapeyron
!===============================================================================
!===============================================================================
end module TRACERS_SOA
!===============================================================================
