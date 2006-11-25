#!/usr/bin/python

SOURCE_SUFFIX = ".f"
OBJECT_SUFFIX = ".o"

import re, sys, os

def usage():
    print """makedep.py --progs <prog_files> --other <other_files>

Generates Makefile.deps for the given program files and supporting
other files, based on lines of the form "use mod_<filename>" in the
code."""

def get_deps(filename):
    deps = []
    full_filename = filename + SOURCE_SUFFIX
    f = open(full_filename)
    if not f:
	print "ERROR: unable to open %s%s" % full_filename
	sys.exit(1)
    use_re = re.compile("^\s*use mod_(\S*)$")
    for line in f:
	match = use_re.search(line)
	if match:	
	    dep = match.group(1)
	    if dep not in deps:
		deps.append(dep)
    f.close()
    return deps

def expand_deps(dep_lists):
    changed = False
    for f in dep_lists:
	for d in dep_lists[f]:
	    if d not in dep_lists.keys():
		print "ERROR: unknown dependency: " + d
		sys.exit(1)
	    for new_d in dep_lists[d]:
		if new_d not in dep_lists[f]:
		    dep_lists[f].append(new_d)
		    changed = True
    return changed

def full_expand_deps(dep_lists):
    max_expands = 100
    expands = 0
    changed = True
    while changed:
	if expands == max_expands:
	    print "ERROR: expansion failed to terminate"
	    sys.exit(1)
	changed = expand_deps(dep_lists)
	expands = expands + 1

def get_dep_lists(files):
    dep_lists = {}
    for f in files:
	dep_lists[f] = get_deps(f)
    return dep_lists

def print_deps(progs, other, dep_lists):
    out_filename = "Makefile.deps"
    outf = open(out_filename, "w")
    outf.write("#\n# DO NOT EDIT --- auto-generated file\n#\n\n")
    for f in dep_lists.keys():
	if dep_lists[f]:
	    deps_str = " ".join([(d + ".o") for d in dep_lists[f]])
	    outf.write("%s.o: %s\n" % (f, deps_str))
    for f in progs:
	deps_str = " ".join([(d + ".o") for d in dep_lists[f]])
	outf.write("\n%s: %s.o %s\n\t$(FC) $(LDFLAGS) -o $@ %s.o %s\n"
		   % (f, f, deps_str, f, deps_str))
    outf.close()

def process_args(args):
    reading = 0 # 0 = None, 1 = progs, 2 = other
    progs = []
    other = []
    for arg in args:
	if arg == "--progs":
	    reading = 1
	elif arg == "--other":
	    reading = 2
	else:
	    if reading == 1: # progs
		progs.append(arg)
	    elif reading == 2: # other
		other.append(arg)
	    else:
		usage()
		sys.exit(1)
    return (progs, other)

def main():
    if len(sys.argv) < 2:
	usage()
	sys.exit(1)
    (progs, other) = process_args(sys.argv[1:])
    all_files = progs + other
    dep_lists = get_dep_lists(all_files)
    full_expand_deps(dep_lists)
    print_deps(progs, other, dep_lists)

if __name__ == "__main__":
    main()
