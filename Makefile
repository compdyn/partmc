
VERSION = 0.9.0
DIST_NAME = partmc-$(VERSION)

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

PROGS = process_out process_state run_golovin_exact			\
	run_golovin_fix_hybrid run_constant_exact			\
	run_constant_fix_hybrid run_sedi_fix_hybrid run_sedi_ode	\
	run_sedi_sect run_brown_fix_hybrid average

OTHER = array array_hybrid bin condensation constants environ	\
	init_dist kernel_golovin kernel_sedi kernel_constant	\
	kernel_brown material mc_exact mc_fix_hybrid util state

FILES = $(PROGS) $(OTHER)

all: TAGS $(PROGS)

-include Makefile.deps

deps: $(patsubst %,%.f,$(FILES))
	./makedeps.py --progs $(PROGS) --other $(OTHER)

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o : %.mod

clean:
	rm -f $(PROGS) *.o *.mod TAGS

distclean: clean
	rm -f Makefile.deps

cleanall: clean
	rm -f *~ *.d gmon.out gprof_*

check:
	ftnchek-3.3.1/ftnchek *.f

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<

dist: Makefile.deps
	mkdir $(DIST_NAME)
	cp Makefile Makefile.deps makedeps.py TODO COPYING README $(patsubst %,%.f,$(FILES)) $(DIST_NAME)
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -r $(DIST_NAME)

TAGS:
	etags $(patsubst %,%.f,$(FILES))
