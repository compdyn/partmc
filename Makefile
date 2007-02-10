
VERSION = 1.0.0
DIST_NAME = partbox-$(VERSION)

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

PROGS = process_out process_state test_sedi_bidisperse_ode average	\
	partbox test_sedi_bidisperse_state_to_count

OTHER = array bin condensation constants environ init_dist	\
	kernel_golovin kernel_sedi kernel_constant kernel_brown	\
	material run_exact run_mc util run_sect state read_spec

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
	mkdir $(DIST_NAME)
	cp Makefile Makefile.deps makedeps.py TODO COPYING README $(ALL_SOURCE) $(DIST_NAME)
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -r $(DIST_NAME)

TAGS: $(ALL_SOURCE)
	etags $(ALL_SOURCE)
