VERSION = 1.0.0
DIST_NAME = hpmc-$(VERSION)

# useful flags:
#   -O              optimize
#   -g              debugging
#   -pg             profiling
#   -pedantic       strict F77
#   -fbounds-check  check array accesses
#FFLAGS = -O -fcase-preserve -W -Wall -Wimplicit -Wsurprising -Wunused -Wuninitialized
# for gfortran:
#FFLAGS = -O -fimplicit-none -W -Wall -Wunused -Wconversion -Wunderflow -Wunused-labels
FFLAGS = -O -fimplicit-none -w
LDFLAGS = 

F77 = gfortran

PROGS = \
	process_out \
	run_golovin_adapt \
	run_golovin_exact \
	run_golovin_fix \
	run_golovin_var \
	run_sedi_adapt \
	run_sedi_fix \
	run_sedi_fix_hybrid \
	run_sedi_fix_split \
	run_sedi_fix_super \
	run_sedi_ode \
	run_sedi_sect \
	run_sedi_var

process_out_objs = \
	process_out.o
run_golovin_adapt_objs = \
	run_golovin_adapt.o \
	mc_adapt.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
run_golovin_exact_objs = \
	run_golovin_exact.o \
	mc_exact.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
run_golovin_fix_objs = \
	run_golovin_fix.o \
	mc_fix.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
run_golovin_var_objs = \
	run_golovin_var.o \
	mc_var.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
run_sedi_adapt_objs = \
	run_sedi_adapt.o \
	mc_adapt.o \
	kernel_sedi.o \
	array.o \
	init_dist.o
run_sedi_fix_objs = \
	run_sedi_fix.o \
	mc_fix.o \
	kernel_sedi.o \
	array.o \
	init_dist.o
run_sedi_fix_hybrid_objs = \
	run_sedi_fix_hybrid.o \
	mc_fix_hybrid.o \
	kernel_sedi.o \
	array.o \
	array_hybrid.o \
	bin.o \
	init_dist.o
run_sedi_fix_split_objs = \
	run_sedi_fix_split.o \
	mc_fix_split.o \
	kernel_sedi.o \
	array.o \
	array_split.o \
	init_dist.o
run_sedi_fix_super_objs = \
	run_sedi_fix_super.o \
	mc_fix_super.o \
	kernel_sedi.o \
	array.o \
	array_super.o \
	bin.o \
	init_dist.o \
	util.o
run_sedi_ode_objs = \
	kernel_sedi.o \
	run_sedi_ode.o
run_sedi_sect_objs = \
	kernel_sedi.o \
	array.o \
	run_sedi_sect.o
run_sedi_var_objs = \
	run_sedi_var.o \
	mc_var.o \
	kernel_sedi.o \
	array.o \
	init_dist.o

ALL_OBJS = $(foreach PROG,$(PROGS),$($(PROG)_objs))
ALL_SOURCE = $(ALL_OBJS:.o=.f)

all: $(PROGS)

process_out: $(process_out_objs)
	$(F77) $(LDFLAGS) -o $@ $(process_out_objs)

run_golovin_adapt: $(run_golovin_adapt_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_golovin_adapt_objs)

run_golovin_exact: $(run_golovin_exact_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_golovin_exact_objs)

run_golovin_fix: $(run_golovin_fix_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_golovin_fix_objs)

run_golovin_var: $(run_golovin_var_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_golovin_var_objs)

run_sedi_adapt: $(run_sedi_adapt_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_adapt_objs)

run_sedi_fix: $(run_sedi_fix_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_objs)

run_sedi_fix_hybrid: $(run_sedi_fix_hybrid_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_hybrid_objs)

run_sedi_fix_split: $(run_sedi_fix_split_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_split_objs)

run_sedi_fix_super: $(run_sedi_fix_super_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_super_objs)

run_sedi_ode: $(run_sedi_ode_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_ode_objs)

run_sedi_sect: $(run_sedi_sect_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_sect_objs)

run_sedi_var: $(run_sedi_var_objs)
	$(F77) $(LDFLAGS) -o $@ $(run_sedi_var_objs)

%.o: %.f
	$(F77) $(FFLAGS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o

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
