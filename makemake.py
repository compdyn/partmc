#!/usr/bin/python

import re, sys

f77 = "gfortran"
fflags = "-O -fimplicit-none -W -Wall -Wunused-labels -Wconversion -Wunderflow -Wimplicit-interface -Wno-unused"
ldflags = ""

progs = ["process_out",
	 "process_state",
	 "run_golovin_adapt",
	 "run_golovin_exact",
	 "run_golovin_fix",
	 "run_golovin_var",
	 "run_sedi_adapt",
	 "run_sedi_fix",
	 "run_sedi_fix_hybrid",
	 "run_sedi_fix_split",
	 "run_sedi_fix_super",
	 "run_sedi_ode",
	 "run_sedi_sect",
	 "run_sedi_var",
	 "condensation_plot",
	 ]

other = ["array",
	 "array_hybrid",
	 "array_split",
	 "array_super",
	 "bin",
	 "condensation",
	 "constants",
	 "environ",
	 "init_dist",
	 "kernel_golovin",
	 "kernel_sedi",
	 "material",
	 "mc_adapt",
	 "mc_exact",
	 "mc_fix",
	 "mc_fix_hybrid",
	 "mc_fix_split",
	 "mc_fix_super",
	 "mc_var",
	 "util",
	 "state",
	 ]

free_form = ["condensation",
	     "constants",
	     "environ",
	     "material",
	     "process_state",
	     "state",
	     ]

all_files = progs + other

def get_deps(file):
    deps = []
    f = open(file + ".f")
    use_re = re.compile("^\s*use mod_(\S*)$")
    for line in f:
	match = use_re.search(line)
	if match:	
	    dep = match.group(1)
	    if dep not in deps:
		deps.append(dep)
    f.close()
    return deps

def get_dep_list(files):
    deps = {}
    for f in files:
	deps[f] = get_deps(f)
    return deps

def expand_deps(deps):
    changed = False
    for f in deps:
	for d in deps[f]:
	    if d not in deps.keys():
		print "ERROR: unknown dependency: " + d	
		sys.exit(1)
	    for new_d in deps[d]:
		if new_d not in deps[f]:
		    deps[f].append(new_d)
		    changed = True
    return changed

def full_expand_deps(deps):
    max_expands = 100
    expands = 0
    changed = True
    while changed:
	if expands == max_expands:
	    print "ERROR: expansion failed to terminate"
	    sys.exit(1)
	changed = expand_deps(deps)
	expands = expands + 1

def print_deps(obj, deps):
    if deps:
	print "%s.o: %s" % (obj, " ".join([(d + ".o") for d in deps]))

######################################################################

print """
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
"""

print "PROGS = " + " ".join(progs)
print
print "OTHER = " + " ".join(other)
print
print "FILES = $(PROGS) $(OTHER)"
print

print "# temporary hack"
print "FREEFORM = " + " ".join(free_form)
print

deps = get_dep_list(all_files)
full_expand_deps(deps)
for f in all_files:
    print_deps(f, deps[f])

print """
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
"""

for f in progs:
    dep_objs = " ".join([(d + ".o") for d in deps[f]])
    print "%s: %s.o %s" % (f, f, dep_objs)
    print "\t$(F77) $(LDFLAGS) -o $@ %s.o %s" % (f, dep_objs)
    print

######################################################################
