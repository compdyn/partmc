#!/usr/bin/python

import re, sys

progs = ["process_out",
	 "process_state",
	 "run_golovin_exact",
	 "run_golovin_fix_hybrid",
         "run_constant_exact",
         "run_constant_fix_hybrid",
	 "run_sedi_fix_hybrid",
	 "run_sedi_ode",
	 "run_sedi_sect",
         "average",
	 ]

other = ["array",
	 "array_hybrid",
	 "bin",
	 "condensation",
	 "constants",
	 "environ",
	 "init_dist",
	 "kernel_golovin",
	 "kernel_sedi",
	 "kernel_constant",
         "kernel_brown",
	 "material",
	 "mc_exact",
	 "mc_fix_hybrid",
	 "util",
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
# Auto-generated Makefile --- DO NOT EDIT --- edit makemake.py instead
#

VERSION = 1.1.0
DIST_NAME = partmc-$(VERSION)

F77 = pgf95

ifeq ($(F77),gfortran)
    # -O              optimize
    # -g              debugging
    # -pg             profiling
    # -fbounds-check  check array accesses
    # -Wno-unused     disable reporting of unused variables
  FFLAGS = -O -ffree-form -fimplicit-none \
           -W -Wall -Wunused-labels -Wconversion -Wunderflow \
           -Wimplicit-interface -Wno-unused
  LDFLAGS = 
endif
ifeq ($(F77),pgf95)
    # -Mbounds      array bounds checking
    # -Mdclchk      check for undeclared variables
  FFLAGS = -O -Mfree -Mpreprocess -DUSE_F95_RAND
  LDFLAGS =
endif
"""

print "PROGS = " + " ".join(progs)
print
print "OTHER = " + " ".join(other)
print
print "FILES = $(PROGS) $(OTHER)"
print
print "all: TAGS $(PROGS)"
print

deps = get_dep_list(all_files)
full_expand_deps(deps)
for f in all_files:
    print_deps(f, deps[f])

print """
%.o: %.f
	$(F77) $(FFLAGS) -c -o $@ $<

%.o : %.mod

clean:
	rm -f $(PROGS) *.o *.mod TAGS

cleanall: clean
	rm -f *~ *.d gmon.out gprof_*

check:
	ftnchek-3.3.1/ftnchek *.f

gprof_%: % gmon.out
	gprof -p -q $< gmon.out > gprof_$<

dist:
	mkdir $(DIST_NAME)
	cp makemake.py Makefile COPYING README $(patsubst %,%.f,$(FILES)) $(DIST_NAME)
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
