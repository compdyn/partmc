
VERSION = 1.0.0
DIST_NAME = hpmc-$(VERSION)

# useful flags:
#   -O              optimize
#   -g              debugging
#   -pg             profiling
#   -fbounds-check  check array accesses
# FIXME: remove -Wno-unused to start reporting unused variables again
FFLAGS = -O -fimplicit-none -W -Wall -Wunused-labels -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused
LDFLAGS = 

F77 = gfortran

PROGS = process_out run_golovin_adapt run_golovin_exact run_golovin_fix run_golovin_var run_sedi_adapt run_sedi_fix run_sedi_fix_hybrid run_sedi_fix_split run_sedi_fix_super run_sedi_ode run_sedi_sect run_sedi_var condensation_plot

# temporary hack
FREEFORM = condensation constants environ material

run_golovin_adapt.o: array.o init_dist.o mc_adapt.o bin.o
run_sedi_fix_hybrid.o: array.o init_dist.o mc_fix_hybrid.o kernel_sedi.o condensation.o environ.o material.o constants.o bin.o array_hybrid.o util.o
condensation_plot.o: condensation.o array.o array_hybrid.o bin.o environ.o material.o util.o constants.o
array.o: bin.o
array_hybrid.o: array.o bin.o util.o
condensation.o: array.o array_hybrid.o bin.o environ.o material.o util.o constants.o
environ.o: constants.o
mc_fix_hybrid.o: array.o array_hybrid.o condensation.o environ.o material.o bin.o util.o constants.o

# temporary hack
freeflag = $(if $(findstring $(1),$(FREEFORM)),-ffree-form,-ffixed-form)

%.o: %.f
	$(F77) $(FFLAGS) $(call freeflag,$(basename $<)) -c -o $@ $<

all: TAGS $(PROGS)

%.o : %.mod

clean:
	rm -f $(PROGS) *.o *.mod

cleanall: clean
	rm -f *~ *.d gmon.out gprof_*

check:
	ftnchek-3.3.1/ftnchek *.f

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<

dist:
	mkdir $(DIST_NAME)
	cp Makefile $(ALL_SOURCE) $(DIST_NAME)
	tar czf $(DIST_NAME).tar.gz $(DIST_NAME)
	rm -r $(DIST_NAME)

TAGS:
	etags $(patsubst %,%.f,$(FILES))

