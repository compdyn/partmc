# useful options:
#   -g                   debugging
#   -pg                  profiling
#   -pedantic            strict F77
#   -fbounds-check       check array accesses
#   -malign-double       align real*8 on 64-bit boundaries
#   -funroll-all-loops   unroll "do" and "do while" loops
FFLAGS = -O2 -fcase-preserve -W -Wall -Wimplicit -Wsurprising -Wunused -Wuninitialized
LDFLAGS = 

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
	run_sedi_ode \
	run_sedi_sect \
	run_sedi_var

PROCESS_OUT_OBJS = \
	process_out.o
RUN_GOLOVIN_ADAPT_OBJS = \
	run_golovin_adapt.o \
	mc_adapt.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
RUN_GOLOVIN_EXACT_OBJS = \
	run_golovin_exact.o \
	mc_exact.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
RUN_GOLOVIN_FIX_OBJS = \
	run_golovin_fix.o \
	mc_fix.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
RUN_GOLOVIN_VAR_OBJS = \
	run_golovin_var.o \
	mc_var.o \
	kernel_golovin.o \
	array.o \
	init_dist.o
RUN_SEDI_ADAPT_OBJS = \
	run_sedi_adapt.o \
	mc_adapt.o \
	kernel_sedi.o \
	array.o \
	init_dist.o
RUN_SEDI_FIX_OBJS = \
	run_sedi_fix.o \
	mc_fix.o \
	kernel_sedi.o \
	array.o \
	init_dist.o
RUN_SEDI_FIX_HYBRID_OBJS = \
	run_sedi_fix_hybrid.o \
	mc_fix_hybrid.o \
	kernel_sedi.o \
	array.o \
	array_hybrid.o \
	bin.o \
	init_dist.o
RUN_SEDI_FIX_SPLIT_OBJS = \
	run_sedi_fix_split.o \
	mc_fix_split.o \
	kernel_sedi.o \
	array.o \
	array_split.o \
	init_dist.o
RUN_SEDI_ODE_OBJS = \
	kernel_sedi.o \
	run_sedi_ode.o
RUN_SEDI_SECT_OBJS = \
	kernel_sedi.o \
	run_sedi_sect.o
RUN_SEDI_VAR_OBJS = \
	run_sedi_var.o \
	mc_var.o \
	kernel_sedi.o \
	array.o \
	init_dist.o

all: $(PROGS)

process_out: $(PROCESS_OUT_OBJS)
	g77 $(LDFLAGS) -o $@ $(PROCESS_OUT_OBJS)

run_golovin_adapt: $(RUN_GOLOVIN_ADAPT_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_GOLOVIN_ADAPT_OBJS)

run_golovin_exact: $(RUN_GOLOVIN_EXACT_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_GOLOVIN_EXACT_OBJS)

run_golovin_fix: $(RUN_GOLOVIN_FIX_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_GOLOVIN_FIX_OBJS)

run_golovin_var: $(RUN_GOLOVIN_VAR_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_GOLOVIN_VAR_OBJS)

run_sedi_adapt: $(RUN_SEDI_ADAPT_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_SEDI_ADAPT_OBJS)

run_sedi_fix: $(RUN_SEDI_FIX_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_SEDI_FIX_OBJS)

run_sedi_fix_hybrid: $(RUN_SEDI_FIX_HYBRID_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_SEDI_FIX_HYBRID_OBJS)

run_sedi_fix_split: $(RUN_SEDI_FIX_SPLIT_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_SEDI_FIX_SPLIT_OBJS)

run_sedi_ode: $(RUN_SEDI_ODE_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_SEDI_ODE_OBJS)

run_sedi_sect: $(RUN_SEDI_SECT_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_SEDI_SECT_OBJS)

run_sedi_var: $(RUN_SEDI_VAR_OBJS)
	g77 $(LDFLAGS) -o $@ $(RUN_SEDI_VAR_OBJS)

%.o: %.f
	g77 $(FFLAGS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o

cleanall: clean
	rm -f *~ *.d gmon.out gprof_*

check:
	ftnchek-3.3.1/ftnchek *.f

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<
