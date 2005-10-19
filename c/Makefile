
CC = gcc
CPPFLAGS = 
CFLAGS = -W -Wall
LDFLAGS = 

SRCS = \
	kernel_golovin.c \
	mc_adapt.c \
	mc_exact.c \
	particle_array.c \
	run_golovin_adapt.c \
	run_golovin_exact.c

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
	$(CC) $(LDFLAGS) -o $@ $(RUN_SEDI_ADAPT_OBJS)

run_sedi_fix: $(RUN_SEDI_FIX_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(RUN_SEDI_FIX_OBJS)

run_sedi_var: $(RUN_SEDI_VAR_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(RUN_SEDI_VAR_OBJS)

run_golovin_adapt: $(RUN_GOLOVIN_ADAPT_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(RUN_GOLOVIN_ADAPT_OBJS)

run_golovin_exact: $(RUN_GOLOVIN_EXACT_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(RUN_GOLOVIN_EXACT_OBJS)

run_golovin_fix: $(RUN_GOLOVIN_FIX_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(RUN_GOLOVIN_FIX_OBJS)

run_golovin_var: $(RUN_GOLOVIN_VAR_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(RUN_GOLOVIN_VAR_OBJS)

process_out: $(PROCESS_OUT_OBJS)
	$(CC) $(LDFLAGS) -o $@ $(PROCESS_OUT_OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o

cleanall: clean
	rm -f *~ *.d

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<

%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

-include $(SRCS:.c=.d)
