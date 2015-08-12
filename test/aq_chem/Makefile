F77       = gfortran
F90       = gfortran
C         = gcc
EXECS	  = test_aq_chem
COMP      = GCC
GCC       = -cpp -pedantic -Wall -fbounds-check -Wuninitialized -Wmaybe-uninitialized -I/usr/local/include
INTEL     = 
PGI       = 
OPTFLAGS  = $($(COMP))
OBJS 	  = constants.o util.o mpi.o rand.o netcdf.o spec_line.o spec_file.o env_state.o aq_integrate.o aq_spec_data.o aq_state.o aq_rxn_line.o aq_rxn_data.o aq_rxn_file.o aq_mech_data.o aq_integrate_f.o aq_chem.o test_aq_chem.o

.PHONY: clean

default: $(EXECS)

all: $(EXECS)

test_aq_chem: $(OBJS) 
	$(F90) $(OPTFLAGS) -o test_aq_chem $(OBJS) -L/usr/local/lib -lnetcdff -lsundials_cvodes -lsundials_nvecserial

test_aq_chem.o:  test_aq_chem.f90
	$(F90) $(OPTFLAGS) -c test_aq_chem.f90
aq_chem.o:  aq_chem.f90
	$(F90) $(OPTFLAGS) -c aq_chem.f90
aq_rxn_data.o:  aq_rxn_data.f90
	$(F90) $(OPTFLAGS) -c aq_rxn_data.f90
aq_mech_data.o:  aq_mech_data.f90
	$(F90) $(OPTFLAGS) -c aq_mech_data.f90
aq_rxn_file.o:  aq_rxn_file.f90
	$(F90) $(OPTFLAGS) -c aq_rxn_file.f90
aq_rxn_line.o:  aq_rxn_line.f90
	$(F90) $(OPTFLAGS) -c aq_rxn_line.f90
aq_state.o:  aq_state.f90
	$(F90) $(OPTFLAGS) -c aq_state.f90
aq_spec_data.o:  aq_spec_data.f90
	$(F90) $(OPTFLAGS) -c aq_spec_data.f90
aq_integrate_f.o:  aq_integrate_f.f90
	$(F90) $(OPTFLAGS) -c aq_integrate_f.f90
aq_integrate.o:  aq_integrate.c
	$(C) -c aq_integrate.c
env_state.o:  pmc_src/env_state.f90
	$(F90) $(OPTFLAGS) -c pmc_src/env_state.f90
spec_file.o:  pmc_src/spec_file.f90
	$(F90) $(OPTFLAGS) -c pmc_src/spec_file.f90
spec_line.o:  pmc_src/spec_line.f90
	$(F90) $(OPTFLAGS) -c pmc_src/spec_line.f90
netcdf.o:  pmc_src/netcdf.f90
	$(F90) $(OPTFLAGS) -c pmc_src/netcdf.f90
mpi.o:  pmc_src/mpi.f90
	$(F90) $(OPTFLAGS) -c pmc_src/mpi.f90
rand.o:  pmc_src/rand.f90
	$(F90) $(OPTFLAGS) -c pmc_src/rand.f90
util.o:  pmc_src/util.f90
	$(F90) $(OPTFLAGS) -c pmc_src/util.f90
constants.o:  pmc_src/constants.f90
	$(F90) $(OPTFLAGS) -c pmc_src/constants.f90

clean:
	/bin/rm -f *.o *.mod *~ PI* $(EXECS) rdb.* startup.*

