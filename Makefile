
#
# Auto-generated Makefile --- DO NOT EDIT
#

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

PROGS = process_out process_state run_golovin_adapt run_golovin_exact run_golovin_fix run_golovin_var run_constant_exact run_constant_fix_hybrid run_sedi_adapt run_sedi_fix run_sedi_fix_hybrid run_sedi_fix_split run_sedi_fix_super run_sedi_ode run_sedi_sect run_sedi_var condensation_plot average

OTHER = array array_hybrid array_split array_super bin condensation constants environ init_dist kernel_golovin kernel_sedi kernel_constant material mc_adapt mc_exact mc_fix mc_fix_hybrid mc_fix_split mc_fix_super mc_var util state

FILES = $(PROGS) $(OTHER)

# temporary hack
FREEFORM = condensation constants environ material process_state state average

process_state.o: bin.o environ.o material.o array_hybrid.o state.o util.o constants.o array.o
run_golovin_adapt.o: array.o init_dist.o mc_adapt.o bin.o material.o environ.o constants.o
run_constant_exact.o: bin.o mc_exact.o kernel_constant.o array.o environ.o material.o constants.o
run_constant_fix_hybrid.o: bin.o array.o init_dist.o mc_fix_hybrid.o kernel_constant.o condensation.o environ.o material.o constants.o util.o array_hybrid.o state.o
run_sedi_fix_hybrid.o: bin.o array.o init_dist.o mc_fix_hybrid.o kernel_sedi.o condensation.o environ.o material.o constants.o util.o array_hybrid.o state.o
condensation_plot.o: condensation.o array.o array_hybrid.o bin.o environ.o material.o util.o constants.o
array.o: bin.o material.o environ.o constants.o
array_hybrid.o: array.o bin.o util.o material.o environ.o constants.o
condensation.o: array.o array_hybrid.o bin.o environ.o material.o util.o constants.o
environ.o: constants.o material.o
mc_exact.o: array.o environ.o material.o bin.o constants.o
mc_fix_hybrid.o: array.o array_hybrid.o condensation.o environ.o material.o state.o bin.o util.o constants.o
state.o: environ.o constants.o material.o

# temporary hack
freeflag = $(if $(findstring $(1),$(FREEFORM)),-ffree-form,-ffixed-form)

all: TAGS $(PROGS)

%.o: %.f
	$(F77) $(FFLAGS) $(call freeflag,$(basename $<)) -c -o $@ $<

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

make:
	./makemake.py > Makefile.new
	mv Makefile.new Makefile

process_out: process_out.o 
	$(F77) $(LDFLAGS) -o $@ process_out.o 

process_state: process_state.o bin.o environ.o material.o array_hybrid.o state.o util.o constants.o array.o
	$(F77) $(LDFLAGS) -o $@ process_state.o bin.o environ.o material.o array_hybrid.o state.o util.o constants.o array.o

run_golovin_adapt: run_golovin_adapt.o array.o init_dist.o mc_adapt.o bin.o material.o environ.o constants.o
	$(F77) $(LDFLAGS) -o $@ run_golovin_adapt.o array.o init_dist.o mc_adapt.o bin.o material.o environ.o constants.o

run_golovin_exact: run_golovin_exact.o 
	$(F77) $(LDFLAGS) -o $@ run_golovin_exact.o 

run_golovin_fix: run_golovin_fix.o 
	$(F77) $(LDFLAGS) -o $@ run_golovin_fix.o 

run_golovin_var: run_golovin_var.o 
	$(F77) $(LDFLAGS) -o $@ run_golovin_var.o 

run_constant_exact: run_constant_exact.o bin.o mc_exact.o kernel_constant.o array.o environ.o material.o constants.o
	$(F77) $(LDFLAGS) -o $@ run_constant_exact.o bin.o mc_exact.o kernel_constant.o array.o environ.o material.o constants.o

run_constant_fix_hybrid: run_constant_fix_hybrid.o bin.o array.o init_dist.o mc_fix_hybrid.o kernel_constant.o condensation.o environ.o material.o constants.o util.o array_hybrid.o state.o
	$(F77) $(LDFLAGS) -o $@ run_constant_fix_hybrid.o bin.o array.o init_dist.o mc_fix_hybrid.o kernel_constant.o condensation.o environ.o material.o constants.o util.o array_hybrid.o state.o

run_sedi_adapt: run_sedi_adapt.o 
	$(F77) $(LDFLAGS) -o $@ run_sedi_adapt.o 

run_sedi_fix: run_sedi_fix.o 
	$(F77) $(LDFLAGS) -o $@ run_sedi_fix.o 

run_sedi_fix_hybrid: run_sedi_fix_hybrid.o bin.o array.o init_dist.o mc_fix_hybrid.o kernel_sedi.o condensation.o environ.o material.o constants.o util.o array_hybrid.o state.o
	$(F77) $(LDFLAGS) -o $@ run_sedi_fix_hybrid.o bin.o array.o init_dist.o mc_fix_hybrid.o kernel_sedi.o condensation.o environ.o material.o constants.o util.o array_hybrid.o state.o

run_sedi_fix_split: run_sedi_fix_split.o 
	$(F77) $(LDFLAGS) -o $@ run_sedi_fix_split.o 

run_sedi_fix_super: run_sedi_fix_super.o 
	$(F77) $(LDFLAGS) -o $@ run_sedi_fix_super.o 

run_sedi_ode: run_sedi_ode.o 
	$(F77) $(LDFLAGS) -o $@ run_sedi_ode.o 

run_sedi_sect: run_sedi_sect.o 
	$(F77) $(LDFLAGS) -o $@ run_sedi_sect.o 

run_sedi_var: run_sedi_var.o 
	$(F77) $(LDFLAGS) -o $@ run_sedi_var.o 

condensation_plot: condensation_plot.o condensation.o array.o array_hybrid.o bin.o environ.o material.o util.o constants.o
	$(F77) $(LDFLAGS) -o $@ condensation_plot.o condensation.o array.o array_hybrid.o bin.o environ.o material.o util.o constants.o

average: average.o 
	$(F77) $(LDFLAGS) -o $@ average.o 

