
OPTS = -O -fcase-preserve -W -Wall -Wimplicit -Wsurprising

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
UMSP_MC_OBJS = umsp_mc.o

PROGS = \
	run_sedi_adapt \
	run_sedi_fix \
	run_sedi_var \
	run_golovin_adapt \
	run_golovin_fix \
	run_golovin_var \
	umsp_mc

all: $(PROGS)

run_sedi_adapt: $(RUN_SEDI_ADAPT_OBJS)
	g77 $(OPTS) -o $@ $(RUN_SEDI_ADAPT_OBJS)

run_sedi_fix: $(RUN_SEDI_FIX_OBJS)
	g77 $(OPTS) -o $@ $(RUN_SEDI_FIX_OBJS)

run_sedi_var: $(RUN_SEDI_VAR_OBJS)
	g77 $(OPTS) -o $@ $(RUN_SEDI_VAR_OBJS)

run_golovin_adapt: $(RUN_GOLOVIN_ADAPT_OBJS)
	g77 $(OPTS) -o $@ $(RUN_GOLOVIN_ADAPT_OBJS)

run_golovin_fix: $(RUN_GOLOVIN_FIX_OBJS)
	g77 $(OPTS) -o $@ $(RUN_GOLOVIN_FIX_OBJS)

run_golovin_var: $(RUN_GOLOVIN_VAR_OBJS)
	g77 $(OPTS) -o $@ $(RUN_GOLOVIN_VAR_OBJS)

umsp_mc: $(UMSP_MC_OBJS)
	g77 $(OPTS) -o $@ $(UMSP_MC_OBJS)

%.o: %.f
	g77 $(OPTS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o

cleanall: clean
	rm -f *~ *.d
