#!/usr/bin/python

SOURCE_SUFFIX = ".f"
OBJECT_SUFFIX = ".o"

import re, sys, os

def usage():
    print """makedep.py <filename>

Reads a Fortran file and looks for lines of the form
     use mod_<modulename>
and outputs a make-format file with dependencies on <modulename>.f.
The output filename is the same as the input filename with the
extension changed to .deps.
"""

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
    extra_files = []
    for f in dep_lists:
	for d in dep_lists[f]:
	    if d not in dep_lists.keys():
		if d not in extra_files:
		    extra_files.append(d)
		continue
	    for new_d in dep_lists[d]:
		if new_d not in dep_lists[f]:
		    dep_lists[f].append(new_d)
		    changed = True
    for f in extra_files:
	dep_lists[f] = get_deps(f)
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

def get_full_deps(f):
    dep_lists = {}
    dep_lists[f] = get_deps(f)
    full_expand_deps(dep_lists)
    return dep_lists[f]

def print_deps(file, obj, deps, full_deps):
    deps_str = " ".join([(d + ".o") for d in deps])
    full_deps_str = " ".join([(d + ".o") for d in full_deps])
    file.write("%s_deps = %s\n" % (obj, full_deps_str))
    if deps:
	file.write("%s.o: %s\n" % (obj, deps_str))

def main():
    if len(sys.argv) != 2:
	usage()
	sys.exit()
    in_root = sys.argv[1]
    deps = get_deps(in_root)
    full_deps = get_full_deps(in_root)
    out_filename = in_root + ".deps"
    f = open(out_filename, "w")
    print_deps(f, in_root, deps, full_deps)
    f.close()

if __name__ == "__main__":
    main()
