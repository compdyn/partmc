program box_model

  use pmc_camp_core
  use pmc_camp_state
  use pmc_chem_spec_data
  use pmc_constants

  !! [NO2 photolysis modules]
  use pmc_mechanism_data
  use pmc_rxn_data
  use pmc_rxn_photolysis
  use pmc_rxn_factory

  !! [NO2 photolysis modules]

  implicit none

  type(camp_core_t), pointer :: camp_core
  type(camp_state_t), pointer :: camp_state

  integer(kind=i_kind) :: idx_O3, idx_NO2
  type(chem_spec_data_t), pointer :: chem_spec_data

  integer(kind=i_kind) :: i_time

  !! [NO2 photolysis variables]
  integer(kind=i_kind) :: i_rxn
  character(len=:), allocatable :: photo_label
  type(mechanism_data_t), pointer :: mechanism
  type(rxn_factory_t) :: rxn_factory
  class(rxn_data_t), pointer :: photo_rxn
  type(rxn_update_data_photolysis_t) :: NO2_photolysis

  !! [NO2 photolysis variables]

  camp_core => camp_core_t( "my_config_file.json" )
  call camp_core%initialize( )

  if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
    write(*,*) "Something's gone wrong!"
    stop 3
  end if

  idx_O3 = chem_spec_data%gas_state_id( "O3" )
  idx_NO2 = chem_spec_data%gas_state_id( "NO2" )
  if( idx_O3.eq.0 .or. idx_NO2.eq.0 ) then
    write(*,*) "Missing species!"
    stop 3
  end if

  !! [Find NO2 photolysis]
  if( .not.camp_core%get_mechanism( "my simple mechanism", mechanism ) ) then
    write(*,*) "Missing mechanism!"
    stop 3
  end if

  do i_rxn = 1, mechanism%size( )
    photo_rxn => mechanism%get_rxn( i_rxn )
    select type( photo_rxn )
      class is( rxn_photolysis_t )
        if( photo_rxn%property_set%get_string( "my photo label", photo_label ) ) then
          if( photo_label .eq. "NO2 photolysis" ) then
            call rxn_factory%initialize_update_data( photo_rxn, NO2_photolysis )
          end if
        end if
    end select
  end do
  !! [Find NO2 photolysis]

  call camp_core%solver_initialize( )
  camp_state => camp_core%new_state( )

  camp_state%env_state%temp     = 275.4    ! Temperature in K
  camp_state%env_state%pressure = 101532.2 ! Pressure in Pa
  call camp_state%update_env_state( )

  camp_state%state_var( idx_O3   ) = 13.3 ! [O3] in ppm
  camp_state%state_var( idx_NO2 ) = 3.3  ! [NO2] in ppm

  !! [Set NO2 photolysis]
  call NO2_photolysis%set_rate( 12.2d0 ) ! rate in s-1
  call camp_core%update_rxn_data( NO2_photolysis )
  !! [Set NO2 photolysis]

  write(*,*) "time", "O3", "NO2"
  do i_time = 1, 100
    call camp_core%solve( camp_state, 1.0d-2 ) ! time step in s
    write(*,*) i_time*1.0e-2, &
               camp_state%state_var( idx_O3 ), &
               camp_state%state_var( idx_NO2 )
  end do

  deallocate( camp_core )
  deallocate( camp_state )

end program box_model
