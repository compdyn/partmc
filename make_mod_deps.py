#!/usr/bin/python

import re, sys, os, getopt

def usage():
    print """makemoddeps.py -o <output> <files>

Generates dependencies for the given files, based on lines of the form
"use mod_<filename>" in the code."""

def get_deps_and_mods(filename):
    deps = []
    mods = []
    f = open(filename)
    if not f:
	print "ERROR: unable to open %s%s" % filename
	sys.exit(1)
    use_re = re.compile("^\s*use\s+(mod_\S+)\s*$")
    mod_re = re.compile("^\s*module\s+(\S+)\s*$")
    for line in f:
	match = use_re.search(line)
	if match:	
	    dep = "src/" + match.group(1) + ".mod"
	    if dep not in deps:
		deps.append(dep)
	match = mod_re.search(line)
	if match:	
	    mod = "src/" + match.group(1) + ".mod"
	    if mod not in mods:
		mods.append(mod)
    f.close()
    return (deps, mods)

def write_deps(outf, filename, deps, mods):
    filebase, fileext = os.path.splitext(filename)
    outf.write("%s.o: %s\n" % (filebase, " ".join(deps)))
    for mod in mods:
	outf.write("%s: %s\n"
		   % (mod, " ".join(deps)))

def process_args():
    try:
	opts, args = getopt.getopt(sys.argv[1:], "ho:", ["help", "output="])
    except getopt.GetoptError:
	usage()
	sys.exit(1)
    out_filename = None
    for o, a in opts:
	if o in ("-h", "--help"):
	    usage()
	    sys.exit()
	if o in ("-o", "--output"):
	    out_filename = a
    if len(args) < 1:
	usage()
	sys.exit(1)
    if not out_filename:
	print "ERROR: must specify an output filename"
	usage()
	sys.exit(1)
    return (out_filename, args)

def main():
    out_filename, args = process_args()
    outf = open(out_filename, "w")
    outf.write("# DO NOT EDIT --- auto-generated file\n")
    for filename in args:
	(deps, mods) = get_deps_and_mods(filename)
	if deps:
	    write_deps(outf, filename, deps, mods)
    outf.close()

if __name__ == "__main__":
    main()
