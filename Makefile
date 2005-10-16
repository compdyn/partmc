
OPTS = -O2 -fcase-preserve -W -Wall -Wimplicit -Wsurprising -Wunused -Wuninitialized

RUN_SEDI_ADAPT_OBJS = \
	run_sedi_adapt.o \
	mc_adapt.o \
	kernel_sedi.o \
	particle_array.o
RUN_SEDI_FIX_OBJS = \
	run_sedi_fix.o \
	mc_fix.o \
	kernel_sedi.o \
	particle_array.o
RUN_SEDI_VAR_OBJS = \
	run_sedi_var.o \
	mc_var.o \
	kernel_sedi.o \
	particle_array.o
RUN_GOLOVIN_ADAPT_OBJS = \
	run_golovin_adapt.o \
	mc_adapt.o \
	kernel_golovin.o \
	particle_array.o
RUN_GOLOVIN_EXACT_OBJS = \
	run_golovin_exact.o \
	mc_exact.o \
	kernel_golovin.o \
	particle_array.o
RUN_GOLOVIN_FIX_OBJS = \
	run_golovin_fix.o \
	mc_fix.o \
	kernel_golovin.o \
	particle_array.o
RUN_GOLOVIN_VAR_OBJS = \
	run_golovin_var.o \
	mc_var.o \
	kernel_golovin.o \
	particle_array.o
PROCESS_OUT_OBJS = \
	process_out.o

PROGS = \
	run_sedi_adapt \
	run_sedi_fix \
	run_sedi_var \
	run_golovin_adapt \
	run_golovin_exact \
	run_golovin_fix \
	run_golovin_var \
	process_out

all: $(PROGS)

run_sedi_adapt: $(RUN_SEDI_ADAPT_OBJS)
	g77 $(OPTS) -o $@ $(RUN_SEDI_ADAPT_OBJS)

run_sedi_fix: $(RUN_SEDI_FIX_OBJS)
	g77 $(OPTS) -o $@ $(RUN_SEDI_FIX_OBJS)

run_sedi_var: $(RUN_SEDI_VAR_OBJS)
	g77 $(OPTS) -o $@ $(RUN_SEDI_VAR_OBJS)

run_golovin_adapt: $(RUN_GOLOVIN_ADAPT_OBJS)
	g77 $(OPTS) -o $@ $(RUN_GOLOVIN_ADAPT_OBJS)

run_golovin_exact: $(RUN_GOLOVIN_EXACT_OBJS)
	g77 $(OPTS) -o $@ $(RUN_GOLOVIN_EXACT_OBJS)

run_golovin_fix: $(RUN_GOLOVIN_FIX_OBJS)
	g77 $(OPTS) -o $@ $(RUN_GOLOVIN_FIX_OBJS)

run_golovin_var: $(RUN_GOLOVIN_VAR_OBJS)
	g77 $(OPTS) -o $@ $(RUN_GOLOVIN_VAR_OBJS)

process_out: $(PROCESS_OUT_OBJS)
	g77 $(OPTS) -o $@ $(PROCESS_OUT_OBJS)

%.o: %.f
	g77 $(OPTS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o

cleanall: clean
	rm -f *~ *.d

check:
	ftnchek-3.3.1/ftnchek -declare *.f

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<
