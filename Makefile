
VERSION = 0.9.0
DIST_NAME = partmc-$(VERSION)

F77 = gfortran

ifeq ($(F77),gfortran)
    # -O              optimize
    # -g              debugging
    # -pg             profiling
    # -fbounds-check  check array accesses
    # -Wno-unused     disable reporting of unused variables
  FFLAGS = -O -ffree-form -x f95-cpp-input -fimplicit-none -W -Wall -Wunused-labels -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused
  LDFLAGS = 
endif
ifeq ($(F77),pgf95)
    # -Mbounds      array bounds checking
    # -Mdclchk      check for undeclared variables
  FFLAGS = -O -Mfree -Mpreprocess -DUSE_F95_RAND
  LDFLAGS =
endif

PROGS = process_out process_state run_golovin_exact			\
	run_golovin_fix_hybrid run_constant_exact			\
	run_constant_fix_hybrid run_sedi_fix_hybrid run_sedi_ode	\
	run_sedi_sect average

OTHER = array array_hybrid bin condensation constants environ	\
	init_dist kernel_golovin kernel_sedi kernel_constant	\
	kernel_brown material mc_exact mc_fix_hybrid util state

FILES = $(PROGS) $(OTHER)

all: TAGS $(PROGS)

include $(patsubst %,%.deps,$(FILES))

%.deps: %.f
	./makedeps.py $(patsubst %.f,%,$<)

%.o: %.f
	$(F77) $(FFLAGS) -c -o $@ $<

%.o : %.mod

# rules for building all programs
define set_program_build
$(1): $(1).o
	$$(F77) $$(LDFLAGS) -o $$@ $(1).o $$($(1)_deps)
endef
$(foreach prog,$(PROGS),$(eval $(call set_program_build,$(prog))))

clean:
	rm -f $(PROGS) *.o *.mod *.deps TAGS

cleanall: clean
	rm -f *~ *.d gmon.out gprof_*

check:
	ftnchek-3.3.1/ftnchek *.f

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<

dist:
	mkdir $(DIST_NAME)
	cp makemake.py Makefile COPYING README $(patsubst %,%.f,$(FILES)) $(DIST_NAME)
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -r $(DIST_NAME)

TAGS:
	etags $(patsubst %,%.f,$(FILES))
