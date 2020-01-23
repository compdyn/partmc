program box_model

  use pmc_camp_core
  use pmc_camp_state
  use pmc_chem_spec_data
  use pmc_constants
  use pmc_mechanism_data
  use pmc_rxn_data
  use pmc_rxn_photolysis
  use pmc_rxn_factory

  !! [MPI modules]
#ifdef USE_MPI
  use mpi
  use pmc_mpi
#endif
  !! [MPI modules]

  implicit none

  type(camp_core_t), pointer :: camp_core
  type(camp_state_t), pointer :: camp_state

  integer(kind=i_kind) :: idx_O3, idx_NO, idx_NO2, idx_O2
  type(chem_spec_data_t), pointer :: chem_spec_data

  character(len=*), parameter :: fmt_hdr = "(A10,',',A10,',',A10,',',A10,',',A10)"
  character(len=*), parameter :: fmt_dat = "(ES10.4,',',ES10.4,',',ES10.4,',',ES10.4,',',ES10.4)"

  integer(kind=i_kind) :: i_time

  integer(kind=i_kind) :: i_rxn
  character(len=:), allocatable :: photo_label
  type(mechanism_data_t), pointer :: mechanism
  type(rxn_factory_t) :: rxn_factory
  class(rxn_data_t), pointer :: photo_rxn
  type(rxn_update_data_photolysis_t) :: NO2_photolysis

  !! [MPI variables]
#ifdef USE_MPI
  integer(kind=i_kind) :: pos, pack_size
  character, allocatable :: buffer(:)
#endif
  !! [MPI variables]

  !! [wrap initialization]
#ifdef USE_MPI
  call pmc_mpi_init( )

  if( pmc_mpi_rank( ) .eq. 0 ) then
#endif

    camp_core => camp_core_t( "my_config_file.json" )
    call camp_core%initialize( )
    !! [wrap initialization]

    if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
      write(*,*) "Something's gone wrong!"
      stop 3
    end if

    idx_O3  = chem_spec_data%gas_state_id( "O3"  )
    idx_NO  = chem_spec_data%gas_state_id( "NO"  )
    idx_NO2 = chem_spec_data%gas_state_id( "NO2" )
    idx_O2  = chem_spec_data%gas_state_id( "O2"  )
    if( idx_O3.eq.0 .or. idx_NO2.eq.0 .or.idx_O2.eq.0 ) then
      write(*,*) "Missing species!"
      stop 3
    end if

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

    !! [get pack size]
    if( .not.NO2_photolysis%is_attached( ) ) then
      write(*,*) "Missing NO2 photolysis reaction!"
      stop 3
    end if

#ifdef USE_MPI
    ! Pack the core and the NO2 photolysis update data object
    pack_size = camp_core%pack_size( )               + &
                NO2_photolysis%pack_size( )
    allocate( buffer( pack_size ) )
    !! [get pack size]

    !! [pack objects]
    pos = 0
    call camp_core%bin_pack(      buffer, pos )
    call NO2_photolysis%bin_pack( buffer, pos )

  end if ! primary process
  !! [pack objects]

  !! [pass indices]
  call pmc_mpi_bcast_integer( idx_O3  )
  call pmc_mpi_bcast_integer( idx_NO  )
  call pmc_mpi_bcast_integer( idx_NO2 )
  call pmc_mpi_bcast_integer( idx_O2  )
  !! [pass indices]

  !! [pass the buffer]
  call pmc_mpi_bcast_integer( pack_size )

  if( pmc_mpi_rank( ) .gt. 0 ) then
    allocate( buffer( pack_size ) )
  end if

  call pmc_mpi_bcast_packed( buffer )
  !! [pass the buffer]

  !! [unpack the objects]
  if( pmc_mpi_rank( ) .gt. 0 ) then

    camp_core => camp_core_t( )
    pos = 0
    call camp_core%bin_unpack(      buffer, pos )
    call NO2_photolysis%bin_unpack( buffer, pos )

  end if

  deallocate( buffer )
#endif
  !! [unpack the objects]

  call camp_core%solver_initialize( )
  camp_state => camp_core%new_state( )

  camp_state%env_state%temp     = 275.4    ! Temperature in K
  camp_state%env_state%pressure = 101532.2 ! Pressure in Pa
  call camp_state%update_env_state( )

  camp_state%state_var( idx_O3   ) = 0.13  ! [O3] in ppm
  camp_state%state_var( idx_NO   ) = 0.02  ! [NO] in ppm
  camp_state%state_var( idx_NO2  ) = 0.053 ! [NO2] in ppm
  camp_state%state_var( idx_O2   ) = 2.1e5 ! [O2] in ppm

  !! [Set NO2 photolysis]
  call NO2_photolysis%set_rate( 12.2d0 ) ! rate in s-1
  call camp_core%update_rxn_data( NO2_photolysis )
  !! [Set NO2 photolysis]

  !! [output]
#ifdef USE_MPI
  if( pmc_mpi_rank( ) .eq. 1 ) then
#endif

    write(*,fmt_hdr) "time", "O3", "NO", "NO2", "O2"
    do i_time = 1, 100
      call camp_core%solve( camp_state, 1.0d-15 ) ! time step in s
      write(*,fmt_dat) i_time*1.0e-15, &
                       camp_state%state_var( idx_O3  ), &
                       camp_state%state_var( idx_NO  ), &
                       camp_state%state_var( idx_NO2 ), &
                       camp_state%state_var( idx_O2  )
    end do

#ifdef USE_MPI
  end if

  call pmc_mpi_finalize( )
#endif
  !! [output]

  deallocate( camp_core )
  deallocate( camp_state )

end program box_model
