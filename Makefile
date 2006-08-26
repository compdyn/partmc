VERSION = 1.0.0
DIST_NAME = hpmc-$(VERSION)

# useful flags:
#   -O              optimize
#   -g              debugging
#   -pg             profiling
#   -fbounds-check  check array accesses
FFLAGS = -O -fimplicit-none -W -Wall -Wunused -Wconversion -Wunderflow -Wunused-labels -Wimplicit-interface
LDFLAGS = 

F77 = gfortran

PROGS = process_out run_golovin_adapt run_golovin_exact			\
	run_golovin_fix run_golovin_var run_sedi_adapt run_sedi_fix	\
	run_sedi_fix_hybrid run_sedi_fix_split run_sedi_fix_super	\
	run_sedi_ode run_sedi_sect run_sedi_var condensation_plot

OTHER = array array_hybrid array_split array_super bin condensation	\
	init_dist kernel_golovin kernel_sedi mc_adapt mc_exact mc_fix	\
	mc_fix_hybrid mc_fix_split mc_fix_super mc_var util

FILES := $(PROGS) $(OTHER)

all: TAGS $(PROGS)

# These rules assume that USE statements are like "^ use module$" with
# no extra stuff on the line, and that the module names are the same as
# the filenames containing the modules with "mod_" prepended, and that
# all dependencies are expressed via USE statements.

# set variables like filename.deps to the immediate dependencies of filename.f
# and also set filename.objs to the corresponding .o files
find_deps = $(patsubst mod_%,%,$(strip $(filter-out use,$(shell grep '^ *use ' $(1)))))
define set_deps
$(1).deps := $$(call find_deps,$(1).f)
$(1).objs := $$(addsuffix .o,$$($(1).deps))
endef
$(foreach file,$(FILES),$(eval $(call set_deps,$(file))))

# set variables like filename.alldeps to recursively expanded deps of filename.f
# and also set filename.allobjs to the corresponding .o files
recurse_deps = $(sort $(1) $(foreach f,$(1),$($(f).deps)))
find_alldeps = $(if $(filter-out $(1),$(call recurse_deps,$(1))),$(call find_alldeps,$(call recurse_deps,$(1))),$(1))
define set_alldeps
$(1).alldeps := $$(call find_alldeps,$(1))
$(1).allobjs := $$(addsuffix .o,$$($(1).alldeps))
endef
$(foreach file,$(FILES),$(eval $(call set_alldeps,$(file))))

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
