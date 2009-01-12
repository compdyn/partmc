
add_executable(partmc partmc.f90 aero_state.f90 aero_binned.f90
		bin_grid.f90 condensation.f90 constants.f90
		env_data.f90 env_state.f90 aero_dist.f90
		kernel_golovin.f90 kernel_sedi.f90 kernel_constant.f90
		kernel_brown.f90 kernel_zero.f90 aero_data.f90
		run_exact.f90 run_mc.f90 util.f90 run_sect.f90
		output_state.f90 output_state_netcdf.f90 mosaic.f90
		gas_data.f90 gas_state.f90 coagulation.f90 kernel.f90
		output_processed.f90 inout.f90 rand.f90
		aero_particle.f90 aero_particle_array.f90 mpi.f90
		process_spec.f90 netcdf.f90)

SET_TARGET_PROPERTIES(partmc PROPERTIES COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(partmc PROPERTIES LINK_FLAGS "${EXTRA_LINK_FLAGS}")
IF(${EXTRA_DEFINITIONS})
  SET_TARGET_PROPERTIES(partmc PROPERTIES COMPILE_DEFINITIONS ${EXTRA_DEFINITIONS})
ENDIF(${EXTRA_DEFINITIONS})
TARGET_LINK_LIBRARIES(partmc ${EXTRA_LIBS})
