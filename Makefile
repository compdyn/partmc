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

MOSAIC_LIBDIR = /Users/mwest/t/mos/MOSAIC.24/compile/
MOSAIC_MODDIR = /Users/mwest/t/mos/MOSAIC.24/compile/

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

PROGS := src/process_out src/process_state src/process_average	\
	src/partmc test/sedi_bidisperse_ode			\
	test/sedi_bidisperse_state_to_count src/equilib

OTHER := src/array src/bin src/condensation src/constants src/environ	\
	src/init_dist src/kernel_golovin src/kernel_sedi		\
	src/kernel_constant src/kernel_brown src/material		\
	src/run_exact src/run_mc src/util src/run_sect src/state	\
	src/read_spec src/mosaic src/gas

EXTRA_DIST := dust_salt.sh dust_salt_part1.spec dust_salt_part2.spec	\
	golovin.sh golovin_exact.spec golovin_mc.spec			\
	sedi_bidisperse.sh sedi_bidisperse_mc.spec sedi_exp.sh		\
	sedi_exp_mc.spec sedi_exp_sect.spec

process_out_OBJS := src/process_out.o src/util.o src/constants.o
process_state_OBJS := src/process_state.o src/bin.o src/environ.o	\
	src/material.o src/array.o src/state.o src/util.o		\
	src/constants.o
process_average_OBJS := src/process_average.o
partmc_OBJS := src/partmc.o src/read_spec.o src/bin.o src/array.o	\
	src/init_dist.o src/condensation.o src/kernel_sedi.o		\
	src/kernel_golovin.o src/kernel_constant.o src/kernel_brown.o	\
	src/material.o src/environ.o src/run_mc.o src/gas.o		\
	src/run_exact.o src/run_sect.o src/util.o src/constants.o	\
	src/state.o src/mosaic.o
sedi_bidisperse_ode_OBJS := test/sedi_bidisperse_ode.o		\
	src/kernel_sedi.o src/environ.o src/constants.o src/material.o	\
	src/util.o
sedi_bidisperse_state_to_count_OBJS :=			\
	test/sedi_bidisperse_state_to_count.o src/environ.o	\
	src/material.o src/state.o src/array.o src/constants.o	\
	src/util.o src/bin.o
equilib_OBJS := src/equilib.o src/material.o src/environ.o		\
	src/condensation.o src/read_spec.o src/util.o src/array.o	\
	src/constants.o src/gas.o src/bin.o

ALL_FILES = $(PROGS) $(OTHER)
ALL_SOURCE = $(patsubst %,%.f90,$(ALL_FILES))
ALL_OBJS = $(patsubst %,%.o,$(ALL_FILES))
ALL_DEPS = $(patsubst %,%.deps,$(ALL_FILES))

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
-include $(patsubst %,%.deps,$(ALL_FILES))
%.deps: %.f90 make_mod_deps.py
	./make_mod_deps.py -o $@ $<

src/%.o src/mod_%.mod: src/%.f90 src/%.deps
	$(FC) $(FFLAGS) -c -o $(patsubst %.f90,%.o,$<) $<
test/%.o test/mod_%.mod: test/%.f90 test/%.deps
	$(FC) $(FFLAGS) -c -o $(patsubst %.f90,%.o,$<) $<

src/process_out: $(process_out_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(process_out_OBJS)
src/process_state: $(process_state_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(process_state_OBJS)
src/process_average: $(process_average_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(process_average_OBJS)
src/partmc: $(partmc_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(partmc_OBJS) -lmosaic
test/sedi_bidisperse_ode: $(sedi_bidisperse_ode_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(sedi_bidisperse_ode_OBJS)
test/sedi_bidisperse_state_to_count: $(sedi_bidisperse_state_to_count_OBJS)
	$(FC) $(LDFLAGS) -o $@ $(sedi_bidisperse_state_to_count_OBJS)
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
