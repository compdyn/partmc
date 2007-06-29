.SUFFIXES:
.SUFFIXES: .f90 .o .deps .mod

VERSION = 1.0.0
DIST_NAME = partmc-$(VERSION)
DATE := $(shell date +"%Y-%m-%d")

# set to "yes" if building as a developer, otherwise "no"
DEV_BUILD = yes

# run "make FC=pgf95" or "make FC=pgf90" o use the Portland Grou
#  compiler instead
FC = gfortran

MOSAIC_DIR = $(HOME)/proj/mosaic/trunk/compile/
MOSAIC_LIBDIR = $(MOSAIC_DIR)
MOSAIC_MODDIR = $(MOSAIC_DIR)

ifeq ($(FC),gfortran)
    # -O              optimize
    # -g              debugging
    # -pg             profiling
    # -fbounds-check  check array accesses
    # -Wno-unused     disable reporting of unused variables
  FFLAGS = -g -Jsrc -Isrc -x f95-cpp-input -fimplicit-none -W -Wall -Wunused-labels -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused -I$(MOSAIC_MODDIR) -fbounds-check
  LDFLAGS = -L$(MOSAIC_LIBDIR)
endif
ifeq ($(FC),pgf95)
    # -Mbounds      array bounds checking
    # -Mdclchk      check for undeclared variables
  FFLAGS = -O -Mpreprocess -DUSE_F95_RAND -I$(MOSAIC_MODDIR)
  LDFLAGS = -L$(MOSAIC_LIBDIR)
endif
ifeq ($(FC),pgf90)
  FFLAGS = -O -Mpreprocess -DUSE_F95_RAND -I$(MOSAIC_MODDIR)
  LDFLAGS = -L$(MOSAIC_LIBDIR)
endif

-include Makefile.local

PROGS := src/process_summary src/process_full src/process_average	\
	src/partmc test/sedi_bidisperse_ode				\
	test/sedi_bidisperse_output_to_count src/equilib

OTHER := src/aero_state src/bin src/condensation src/constants		\
	src/environ src/aero_dist src/kernel_golovin src/kernel_sedi	\
	src/kernel_constant src/kernel_brown src/aero_data		\
	src/run_exact src/run_mc src/util src/run_sect			\
	src/output_state src/read_spec src/mosaic src/gas_data		\
	src/gas_state src/coagulation src/kernel src/output_summary

EXTRA_DIST := dust_salt.sh dust_salt_part1.spec dust_salt_part2.spec	\
	golovin.sh golovin_exact.spec golovin_mc.spec			\
	sedi_bidisperse.sh sedi_bidisperse_mc.spec sedi_exp.sh		\
	sedi_exp_mc.spec sedi_exp_sect.spec

process_summary_OBJS := src/process_summary.o src/util.o src/constants.o
process_full_OBJS := src/process_full.o src/bin.o src/environ.o		\
	src/aero_data.o src/aero_state.o src/output_state.o src/util.o	\
	src/constants.o src/gas_data.o src/gas_state.o
process_average_OBJS := src/process_average.o
partmc_OBJS := src/partmc.o src/read_spec.o src/bin.o src/aero_state.o	\
	src/aero_dist.o src/condensation.o src/kernel_sedi.o		\
	src/kernel_golovin.o src/kernel_constant.o src/kernel_brown.o	\
	src/aero_data.o src/environ.o src/run_mc.o src/gas_data.o	\
	src/gas_state.o src/run_exact.o src/run_sect.o src/util.o	\
	src/constants.o src/output_state.o src/mosaic.o			\
	src/coagulation.o src/kernel.o src/output_summary.o
sedi_bidisperse_ode_OBJS := test/sedi_bidisperse_ode.o			\
	src/kernel_sedi.o src/environ.o src/constants.o			\
	src/aero_data.o src/util.o src/gas_data.o src/gas_state.o	\
	src/aero_state.o src/bin.o
sedi_bidisperse_output_to_count_OBJS :=				\
	test/sedi_bidisperse_output_to_count.o src/environ.o	\
	src/aero_data.o src/output_state.o src/aero_state.o	\
	src/constants.o src/util.o src/bin.o src/gas_data.o	\
	src/gas_state.o
equilib_OBJS := src/equilib.o src/aero_data.o src/environ.o		\
	src/condensation.o src/read_spec.o src/util.o src/aero_state.o	\
	src/constants.o src/gas_data.o src/gas_state.o src/bin.o	\
	src/aero_dist.o

ALL_FILES = $(PROGS) $(OTHER)
ALL_SOURCE = $(ALL_FILES:%=%.f90)
ALL_OBJS = $(ALL_FILES:%=%.o)
ALL_DEPS = $(ALL_FILES:%=%.deps)

.PHONY: all

ifeq ($(DEV_BUILD),yes)
# developers should rebuild Makefile.deps and TAGS
all: TAGS $(PROGS)

# centralized dependencies
#Makefile.deps: make_mod_deps.py $(ALL_SOURCE)
#	./make_mod_deps.py -o $@ $(ALL_SOURCE)

TAGS: $(ALL_SOURCE)
	etags $(ALL_SOURCE)
else
# non-developers should only build the programs
all: $(PROGS)
endif

# centralized dependencies
# also need to remove the %.deps dependency for each source file
# if we use a centralized dependency file. We can make each source file
# depend on Makefile.deps, but then we are forced to do a complete rebuild
# if an file changes.
#-include Makefile.deps

# we can also do per-sourcefile deps, instead of a single Makefile.deps
-include $(ALL_FILES:%=%.deps)
%.deps: %.f90 make_mod_deps.py
	./make_mod_deps.py -o $@ $<

src/%.o src/mod_%.mod: src/%.f90 src/%.deps
	$(FC) $(FFLAGS) -c -o $(patsubst %.f90,%.o,$<) $<
test/%.o test/mod_%.mod: test/%.f90 test/%.deps
	$(FC) $(FFLAGS) -c -o $(patsubst %.f90,%.o,$<) $<

src/process_summary: $(process_summary_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(process_summary_OBJS)
src/process_full: $(process_full_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(process_full_OBJS)
src/process_average: $(process_average_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(process_average_OBJS)
src/partmc: $(partmc_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(partmc_OBJS) -lmosaic
test/sedi_bidisperse_ode: $(sedi_bidisperse_ode_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(sedi_bidisperse_ode_OBJS)
test/sedi_bidisperse_output_to_count: $(sedi_bidisperse_output_to_count_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(sedi_bidisperse_output_to_count_OBJS)
src/equilib: $(equilib_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(equilib_OBJS)

.PHONY: clean
clean:
	rm -f $(PROGS) $(ALL_OBJS)

.PHONY: cleanall
cleanall: clean
	rm -f *~ src/*~ test/*~ src/*.mod

.PHONY: distclean
distclean: cleanall
	rm -f test/*.d test/*.pdf test/*.eps test/gmon.out test/gprof_*

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
