VERSION = 1.0.0
DIST_NAME = hpmc-$(VERSION)

# useful flags:
#   -O              optimize
#   -g              debugging
#   -pg             profiling
#   -pedantic       strict F77
#   -fbounds-check  check array accesses
#FFLAGS = -g -O -fcase-preserve -W -Wall -Wimplicit -Wsurprising -Wunused -Wuninitialized
# for gfortran:
FFLAGS = -O -fimplicit-none -W -Wall -Wunused -Wconversion -Wunderflow -Wunused-labels -Wimplicit-interface
#FFLAGS = -g -O -fimplicit-none -w
LDFLAGS = 

#F77 = g77
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
	run_sedi_var \
	condensation_plot

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
	init_dist.o \
        condensation.o
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
condensation_plot_objs = \
	condensation_plot.o \
	condensation.o \
	array_hybrid.o \
	array.o

ALL_OBJS = $(foreach PROG,$(PROGS),$($(PROG)_objs))
ALL_SOURCE = $(ALL_OBJS:.o=.f)

#condensation.o: array_hybrid.o
#condensation_plot.o: condensation.o

all: $(PROGS)

FILES := $(patsubst %.o,%,$(ALL_OBJS))

# These rules assume that USE statements are like "^      use module$"
# with no extra stuff on the line, and that the module names are the
# same as the filenames containing the modules, and that all dependencies
# are expressed via USE statements.

# set variables like filename.deps to the immediate dependencies of filename.f
find_deps = $(strip $(filter-out use,$(shell grep '^ *use ' $(1))))
define set_deps
$(1).deps := $(call find_deps,$(1).f)
endef
$(foreach file,$(FILES),$(eval $(call set_deps,$(file))))

# set variables like filename.alldeps to recursively expanded deps of filename.f
recurse_deps = $(sort $(1) $(foreach f,$(1),$($(f).deps)))
find_alldeps = $(if $(filter-out $(1),$(call recurse_deps,$(1))),$(call find_alldeps,$(call recurse_deps,$(1))),$(1))
define set_alldeps
$(1).alldeps := $(call find_alldeps,$(1))
endef
$(foreach file,$(FILES),$(eval $(call set_alldeps,$(file))))

# set object file lists
define set_objs
$(1).objs := $(addsuffix .o,$($(1).deps))
$(1).allobjs := $(addsuffix .o,$($(1).alldeps))
endef
$(foreach file,$(FILES),$(eval $(call set_objs,$(file))))

# establish direct dependencies for object files
define set_obj_deps
$(1).o: $(1).f $$($(1).objs)
	$$(F77) $$(FFLAGS) -c -o $$@ $$<
endef
$(foreach file,$(FILES),$(eval $(call set_obj_deps,$(file))))

# establish rules for building all programs
define set_program_build
$(1): $$($(1).allobjs)
	$$(F77) $$(LDFLAGS) -o $$@ $$($(1).allobjs)
endef
$(foreach prog,$(PROGS),$(eval $(call set_program_build,$(prog))))

test:
	echo $(condensation_plot.deps)
	echo $(condensation.deps)
	echo $(call recurse_deps,condensation_plot)
	echo $(condensation_plot.alldeps)
	echo $(addsuffix .o,$(condensation_plot.alldeps))
	echo $(condensation_plot.temp)
	echo $(condensation_plot.objs)
	echo $(condensation_plot.allobjs)

#process_out: $(process_out_objs)
#	$(F77) $(LDFLAGS) -o $@ $(process_out_objs)
#
#run_golovin_adapt: $(run_golovin_adapt_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_golovin_adapt_objs)
#
#run_golovin_exact: $(run_golovin_exact_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_golovin_exact_objs)
#
#run_golovin_fix: $(run_golovin_fix_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_golovin_fix_objs)
#
#run_golovin_var: $(run_golovin_var_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_golovin_var_objs)
#
#run_sedi_adapt: $(run_sedi_adapt_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_adapt_objs)
#
#run_sedi_fix: $(run_sedi_fix_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_objs)
#
#run_sedi_fix_hybrid: $(run_sedi_fix_hybrid_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_hybrid_objs)
#
#run_sedi_fix_split: $(run_sedi_fix_split_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_split_objs)
#
#run_sedi_fix_super: $(run_sedi_fix_super_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_fix_super_objs)
#
#run_sedi_ode: $(run_sedi_ode_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_ode_objs)
#
#run_sedi_sect: $(run_sedi_sect_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_sect_objs)
#
#run_sedi_var: $(run_sedi_var_objs)
#	$(F77) $(LDFLAGS) -o $@ $(run_sedi_var_objs)
#
#condensation_plot: $(condensation_plot_objs)
#	$(F77) $(LDFLAGS) -o $@ $(condensation_plot_objs)

%.o : %.mod

#%.o: %.f
#	$(F77) $(FFLAGS) -c -o $@ $<

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
