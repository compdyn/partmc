! Copyright (C) 2017-2019 Matt Dawson
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_camp_box_model program

!> Driver for the \ref camp_chem "CAMP" box model
program pmc_camp_box_model

  use pmc_camp_box_model_data
  use pmc_util

  implicit none

  ! New-line character
  character(len=*), parameter :: new_line = char(10)
  ! Output file unit
  integer(kind=i_kind), parameter :: OUTPUT_FILE_UNIT = 7

  ! Config file list path
  character(len=300) :: config_file_arg
  character(len=:), allocatable :: config_file
  ! Output file name
  character(len=300) :: output_file_arg
  character(len=:), allocatable :: output_file

  ! The box model data
  type(camp_box_model_data_t), pointer :: box_model

  ! Get the configuration and output file names
  if ( command_argument_count( ) /= 2 ) then
    write(*,*) "Usage: camp_box_model box_model_config.json output_file.txt"
    call die_msg( 695622653, "Incorrect number of command line arguments" )
  end if
  call get_command_argument( 1, config_file_arg )
  call get_command_argument( 2, output_file_arg )
  config_file = trim( config_file_arg )
  output_file = trim( output_file_arg )

  ! Create a new box model
  box_model => camp_box_model_data_t( config_file )

  ! Open the output file
  open( unit = OUTPUT_FILE_UNIT, file = output_file, status = "replace", &
        action = "write" )

  ! Run the camp-chem box model
  call box_model%run( OUTPUT_FILE_UNIT )

  ! Close the output file
  close( OUTPUT_FILE_UNIT )

  ! Free the box model
  deallocate( box_model )

end program pmc_camp_box_model
