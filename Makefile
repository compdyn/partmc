.SUFFIXES:
.SUFFIXES: .f90 .o .deps .mod

VERSION = 1.1.0
DIST_NAME = partmc-$(VERSION)
DATE := $(shell date +"%Y-%m-%d")

# set to "yes" if building as a developer, otherwise "no"
DEV_BUILD = yes

FC = gfortran
# -O              optimize
# -g              debugging
# -pg             profiling (must also be used on LDFLAGS)
# -fbounds-check  check array accesses
# -Wno-unused     disable reporting of unused variables
FFLAGS = -g -Jsrc -Isrc -x f95-cpp-input -fimplicit-none -W -Wall -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused $(MOSAIC_MODDIR) -fbounds-check -Wp,-DPMC_USE_MOSAIC $(NETCDF_MODDIR)
LDFLAGS = $(MOSAIC_LIBDIR) $(NETCDF_LIBDIR)

MOSAIC_MODDIR = -I/usr/include
MOSAIC_LIBDIR = -L/usr/lib
MOSAIC_LIB = -lmosaic

NETCDF_MODDIR = -I/usr/include
NETCDF_LIBDIR = -L/usr/lib
NETCDF_LIB = -lnetcdf

-include Makefile.local

PROGS := src/partmc test/bidisperse/bidisperse_ode equilib/equilib	\
	test/poisson/poisson_sample

CLEAN_DIRS = test/bidisperse/out test/emission/out test/golovin/out	\
        test/mosaic/out test/poisson/out test/sedi/out urban_plume/out

OTHER := src/aero_state src/aero_binned src/bin_grid src/condensation	\
	src/constants src/env_data src/env src/aero_dist		\
	src/kernel_golovin src/kernel_sedi src/kernel_constant		\
	src/kernel_brown src/kernel_zero src/aero_data src/run_exact	\
	src/run_mc src/util src/run_sect src/output_state src/mosaic	\
	src/gas_data src/gas_state src/coagulation src/kernel		\
	src/output_processed src/inout src/rand_poisson			\
	src/aero_particle src/aero_particle_array src/mpi		\
	src/process_spec src/netcdf

EXTRA_DIST := dust_salt.sh dust_salt_part1.spec dust_salt_part2.spec	\
	golovin.sh golovin_exact.spec golovin_mc.spec			\
	sedi_bidisperse.sh sedi_bidisperse_mc.spec sedi_exp.sh		\
	sedi_exp_mc.spec sedi_exp_sect.spec

partmc_OBJS := src/partmc.o src/bin_grid.o src/aero_state.o		\
	src/aero_dist.o src/condensation.o src/kernel_sedi.o		\
	src/kernel_golovin.o src/kernel_constant.o src/kernel_brown.o	\
	src/kernel_zero.o src/aero_data.o src/env_data.o src/env.o	\
	src/run_mc.o src/gas_data.o src/gas_state.o src/run_exact.o	\
	src/run_sect.o src/util.o src/constants.o src/output_state.o	\
	src/mosaic.o src/coagulation.o src/kernel.o			\
	src/output_processed.o src/inout.o src/aero_binned.o		\
	src/rand_poisson.o src/aero_particle.o				\
	src/aero_particle_array.o src/mpi.o src/process_spec.o		\
	src/netcdf.o
process_state_OBJS := src/process_state.o src/bin_grid.o	\
	src/env_data.o src/env.o src/aero_data.o src/aero_state.o	\
	src/output_state.o src/util.o src/constants.o src/gas_data.o	\
	src/gas_state.o src/inout.o src/aero_particle.o			\
	src/aero_particle_array.o src/mpi.o src/aero_dist.o		\
	src/aero_binned.o src/rand_poisson.o src/process_state_hist.o	\
	src/mosaic.o src/process_spec.o src/process.o src/netcdf.o
bidisperse_ode_OBJS := test/bidisperse/bidisperse_ode.o			\
	src/kernel_sedi.o src/env_data.o src/env.o src/constants.o	\
	src/aero_data.o src/util.o src/gas_data.o src/gas_state.o	\
	src/aero_state.o src/bin_grid.o src/inout.o src/aero_dist.o	\
	src/aero_binned.o src/rand_poisson.o src/aero_particle.o	\
	src/aero_particle_array.o src/mpi.o
equilib_OBJS := equilib/equilib.o src/aero_data.o src/env_data.o	\
	src/env.o src/condensation.o src/util.o src/aero_state.o	\
	src/constants.o src/gas_data.o src/gas_state.o src/bin_grid.o	\
	src/aero_dist.o src/inout.o src/aero_binned.o			\
	src/rand_poisson.o src/aero_particle.o				\
	src/aero_particle_array.o src/mpi.o
poisson_sample_OBJS := test/poisson/poisson_sample.o src/util.o	\
	src/rand_poisson.o src/constants.o

ALL_FILES = $(PROGS) $(OTHER)
ALL_SOURCE = $(ALL_FILES:%=%.f90)
ALL_OBJS = $(ALL_FILES:%=%.o)
ALL_DEPS = $(ALL_FILES:%=%.deps)

.PHONY: all

ifeq ($(DEV_BUILD),yes)
# developers should rebuild Makefile.deps and TAGS
all: TAGS $(PROGS)

TAGS: $(ALL_SOURCE)
	etags $(ALL_SOURCE)
else
# non-developers should only build the programs
all: $(PROGS)
endif

-include $(ALL_FILES:%=%.deps)
%.deps: %.f90 tool/f90_mod_deps.py
	tool/f90_mod_deps.py -o $@ -d "(pmc_.*)" -D "src/\1.mod" -m "(.*)" -M "src/\1.mod" $<

src/%.o src/pmc_%.mod: src/%.f90 src/%.deps
	$(FC) $(FFLAGS) -c -o $(patsubst %.f90,%.o,$<) $<
test/%.o: test/%.f90 test/%.deps
	$(FC) $(FFLAGS) -c -o $(patsubst %.f90,%.o,$<) $<
equilib/%.o: equilib/%.f90 equilib/%.deps
	$(FC) $(FFLAGS) -c -o $(patsubst %.f90,%.o,$<) $<

src/partmc: $(partmc_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(partmc_OBJS) $(MOSAIC_LIB) $(NETCDF_LIB)
src/process_state: $(process_state_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(process_state_OBJS) $(MOSAIC_LIB) $(NETCDF_LIB)
equilib/equilib: $(equilib_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(equilib_OBJS)
test/bidisperse/bidisperse_ode: $(bidisperse_ode_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(bidisperse_ode_OBJS)
test/poisson/poisson_sample: $(poisson_sample_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(poisson_sample_OBJS)

.PHONY: clean
clean:
	rm -f TAGS $(PROGS) $(ALL_OBJS) src/*.mod

.PHONY: cleanall
cleanall: clean
	find . -name *~ -exec rm {} \;
	find . -name *.pyc -exec rm {} \;
	find . -name .gdb_history -exec rm {} \;
	find . -name gmon.out -exec rm {} \;
	find . -name gprof_* -exec rm {} \;
	rm -rf test/bidisperse/out/* test/emission/out/*		\
               test/golovin/out/* test/mosaic/out/* test/poisson/out/*	\
               test/sedi/out/* urban_plume/out/*

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<

.PHONY: dist
dist: Makefile.deps TAGS
	grep -q "Version $(VERSION)" README
	grep -q "Released $(DATE)" README
	grep -q "$(VERSION) - $(DATE)" README
	grep -q "DEV_BUILD = no" Makefile
	mkdir $(DIST_NAME)
	cp $(DIST_FILES) $(ALL_SOURCE) $(DIST_NAME)
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -r $(DIST_NAME)
