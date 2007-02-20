
VERSION = 1.0.0
DIST_NAME = partbox-$(VERSION)
DATE := $(shell date +"%Y-%m-%d")

# run "make FC=pgf95" to use the Portland Group compiler instead
FC = gfortran

ifeq ($(FC),gfortran)
    # -O              optimize
    # -g              debugging
    # -pg             profiling
    # -fbounds-check  check array accesses
    # -Wno-unused     disable reporting of unused variables
  FFLAGS = -O -g -ffree-form -x f95-cpp-input -fimplicit-none -W -Wall -Wunused-labels -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused
  LDFLAGS = 
endif
ifeq ($(FC),pgf95)
    # -Mbounds      array bounds checking
    # -Mdclchk      check for undeclared variables
  FFLAGS = -O -Mfree -Mpreprocess -DUSE_F95_RAND
  LDFLAGS =
endif

PROGS = process_out process_state test_sedi_bidisperse_ode		\
	process_average partbox test_sedi_bidisperse_state_to_count

OTHER = array bin condensation constants environ init_dist	\
	kernel_golovin kernel_sedi kernel_constant kernel_brown	\
	material run_exact run_mc util run_sect state read_spec

DIST_FILES = Makefile Makefile.deps makedeps.py TODO COPYING README

TEST_FILES = test_dust_salt.sh test_dust_salt_part1.spec		\
             test_dust_salt_part2.spec test_golovin.sh			\
             test_golovin_exact.spec test_golovin_mc.spec		\
             test_sedi_bidisperse.sh test_sedi_bidisperse_mc.spec	\
             test_sedi_bidisperse_ode.f					\
             test_sedi_bidisperse_state_to_count.f test_sedi_exp.sh	\
             test_sedi_exp_mc.spec test_sedi_exp_sect.spec

ALL_FILES = $(PROGS) $(OTHER)
ALL_SOURCE = $(patsubst %,%.f,$(ALL_FILES))
ALL_OBJS = $(patsubst %,%.o,$(ALL_FILES))

all: Makefile.deps TAGS $(PROGS)

Makefile.deps: $(ALL_SOURCE)
	./makedeps.py --progs $(PROGS) --other $(OTHER)

-include Makefile.deps

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o : %.mod

clean:
	rm -f $(PROGS) *.o *.mod TAGS

cleanall: clean
	rm -f *~ *.d *.pdf *.eps gmon.out gprof_*

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<

dist: Makefile.deps
	grep -q "Version $(VERSION)" README
	grep -q "Released $(DATE)" README
	grep -q "$(VERSION) - $(DATE)" README
	mkdir $(DIST_NAME)
	cp $(DIST_FILES) $(TEST_FILES) $(ALL_SOURCE) $(DIST_NAME)
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -r $(DIST_NAME)

TAGS: $(ALL_SOURCE)
	etags $(ALL_SOURCE)
