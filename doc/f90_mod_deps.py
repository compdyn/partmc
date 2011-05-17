#!/usr/bin/python
# Copyright (C) 2006-2008, 2011 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import re, sys, os, getopt

def usage():
    print """f90_mod_deps.py [options] <file ...>

Generates dependencies for the given Fortran 90 source files based on
module and use statements in them. Options are:

  -h, --help     This help output.

  -o, --output <file>
                 Specify the output file. If unspecified then output
                 is to stdout.

  -d, --dep-re <regexp>
                 Regular expression to match against each module name
                 within use statements. Defaults to matching
                 everything.

  -D, --dep-template <template>
                 Template expression for the dependency produced from
                 each module matched by the dep-re regular expression.

  -m, --mod-re <regexp>
                 Regular expression to match against each module name
                 within module definition statements. Defaults to
                 matching everything.

  -M, --mod-template <template>
                 Template expression for the dependency target
                 produced from each module matched by the mod-re
                 regular expression.

  -g, --graph    Produce output in .dot format suitable for use with
                 graphviz.

  -i, --indirect
                 Allow indirect references (only for graphing), so if
                 A depends on B and C, and B depends on C, then the
                 direct dependency of A on C will dropped as there is
                 an indirect dependency via B.

  -v, --verbose  Turn on verbose debugging output.

For a discussion of managing Fortran 90 dependencies see:
http://tableau.stanford.edu/~mwest/group/Fortran_90_Module_Dependencies

Example:
f90_mod_deps.py --output src/myfile.deps --dep-re "(pmc_.*)" \\
      --dep-template "src/\1.mod" --mod-re "(.*)" \\
      --mod-template "src/\1.mod" src/myfile.f90
"""

# default options
class Opts:
    output = None
    dep_re = "(.*)"
    dep_template = "\\1.mod"
    mod_re = "(.*)"
    mod_template = "\\1.mod"
    output_format = "make"
    indirect = False
    verbose = False

def get_deps_and_mods(filename, opts):
    if opts.verbose:
	sys.stderr.write("Processing %s\n" % filename)
    deps = []
    mods = []
    f = open(filename)
    if not f:
	print "ERROR: unable to open %s%s" % filename
	sys.exit(1)
    use_line_re = re.compile("^\s*use\s+(\S.+)\s*$")
    cont_line_re = re.compile("^(.*)&\s*$")
    mod_line_re = re.compile("^\s*module\s+(\S+)\s*$")
    split_re = re.compile("\s*,\s*")
    dep_re = re.compile(opts.dep_re)
    mod_re = re.compile(opts.mod_re)
    within_use_statement = False
    line_with_use = False
    for line in f:
	match = use_line_re.search(line)
	if match:
	    within_use_statement = True
	    rest_line = match.group(1)
	else:
	    rest_line = line
	if within_use_statement:
	    match = cont_line_re.search(rest_line)
	    if match:
		rest_line = match.group(1)
	    else:
		within_use_statement = False
	    line_items = split_re.split(rest_line.strip())
	    for item in line_items:
		if item:
		    if opts.verbose:
			sys.stderr.write("use: %s\n" % item)
		    match = dep_re.match(item)
		    if match:
			dep = match.expand(opts.dep_template)
			if opts.verbose:
			    sys.stderr.write("matched to: %s\n" % dep)
			if dep not in deps:
			    deps.append(dep)
	else:
	    # not within_use_statement
	    match = mod_line_re.search(line)
	    if match:
		mod_name = match.group(1)
		if opts.verbose:
		    sys.stderr.write("module: %s\n" % mod_name)
		match = mod_re.match(mod_name)
		if match:
		    mod = match.expand(opts.mod_template)
		    if opts.verbose:
			sys.stderr.write("matched to: %s\n" % mod)
		    if mod not in mods:
			mods.append(mod)
    f.close()
    return (deps, mods)

def write_deps_make(outf, filename, deps, mods):
    filebase, fileext = os.path.splitext(filename)
    outf.write("%s.o: %s\n" % (filebase, " ".join(deps)))
    for mod in mods:
	outf.write("%s: %s\n"
		   % (mod, " ".join(deps)))

def write_deps_graphviz(outf, filename, all_deps, has_mods):
    for mod in all_deps.keys():
        if has_mods[mod]:
            props = " peripheries = 1,"
        else:
            props = " peripheries = 2,"
        outf.write(("    node [shape = box,%s"
                    " href = \"\\ref %s.F90\"] %s\n") % (props, mod, mod))
    for mod in all_deps.keys():
        for dep in all_deps[mod]:
            outf.write("    %s -> %s\n" % (mod, dep))

def propagrate_deps(full_deps, opts):
    if opts.verbose:
        sys.stderr.write("propagating full dependencies...\n")
    changes_made = False
    if opts.verbose:
        sys.stderr.write("for modules: %s\n" % " ".join(full_deps.keys()))
    for mod in full_deps.keys():
        new_deps = set()
        for dep in full_deps[mod]:
            new_deps |= full_deps[dep]
        if opts.verbose:
            sys.stderr.write("current dependencies of %s: %s\n"
                             % (mod, " ".join(full_deps[mod])))
            sys.stderr.write("new dependencies of %s: %s\n"
                             % (mod, " ".join(new_deps)))
        if new_deps > full_deps[mod]:
            if opts.verbose:
                sys.stderr.write("adding dependencies to %s: %s\n"
                                 % (mod, " ".join(new_deps - full_deps[mod])))
            changes_made = True
            full_deps[mod] = new_deps
    if opts.verbose:
        if changes_made:
            sys.stderr.write("dependencies were updated\n")
        else:            
            sys.stderr.write("no updates to dependencies\n")
    return changes_made

def remove_with_indirect(all_deps, opts):
    full_deps = {}
    for mod in all_deps.keys():
        full_deps[mod] = set(all_deps[mod])
    while propagrate_deps(full_deps, opts):
        pass
    if opts.verbose:
        sys.stderr.write("computing minimal dependencies...\n")
    for mod in all_deps.keys():
        if opts.verbose:
            sys.stderr.write("computing minimal dependencies for %s...\n" % mod)
        deps = set(all_deps[mod])
        if opts.verbose:
            sys.stderr.write("direct dependencies: %s\n" % " ".join(deps))
        for dep in all_deps[mod]:
            if opts.verbose:
                sys.stderr.write("removing via %s: %s\n"
                                 % (dep, " ".join(full_deps[dep])))
            deps -= full_deps[dep]
        all_deps[mod] = [x for x in deps]
        if opts.verbose:
            sys.stderr.write("minimal dependencies of %s: %s\n"
                             % (mod, " ".join(all_deps[mod])))

def process_args():
    try:
	opts, args = getopt.getopt(sys.argv[1:], "ho:d:D:m:M:vgi",
				   ["help", "output=", "dep-re=",
				    "dep-template=", "mod-re=",
				    "mod-template=", "graph",
                                    "indirect", "verbose"])
    except getopt.GetoptError:
	print "ERROR: invalid commandline options"
	usage()
	sys.exit(1)
    myopts = Opts()
    for o, a in opts:
	if o in ("-h", "--help"):
	    usage()
	    sys.exit()
	if o in ("-o", "--output"):
	    myopts.output = a
	if o in ("-d", "--dep-re"):
	    myopts.dep_re = a
	if o in ("-D", "--dep-template"):
	    myopts.dep_template = a
	if o in ("-m", "--mod-re"):
	    myopts.mod_re = a
	if o in ("-M", "--mod-template"):
	    myopts.mod_template = a
	if o in ("-g", "--graph"):
	    myopts.output_format = "graphviz"
	if o in ("-i", "--indirect"):
	    myopts.indirect = True
	if o in ("-v", "--verbose"):
	    myopts.verbose = True
	    sys.stderr.write("Verbose output on\n")
    if len(args) < 1:
	usage()
	sys.exit(1)
    if myopts.verbose:
	sys.stderr.write("output = %s\n" % myopts.output)
	sys.stderr.write("dep-re = %s\n" % myopts.dep_re)
	sys.stderr.write("dep-template = %s\n" % myopts.dep_template)
	sys.stderr.write("mod-re = %s\n" % myopts.mod_re)
	sys.stderr.write("mod-template = %s\n" % myopts.mod_template)
	sys.stderr.write("output-format = %s\n" % myopts.output_format)
	sys.stderr.write("indirect = %s\n" % myopts.indirect)
	sys.stderr.write("verbose = %s\n" % myopts.verbose)
    return (myopts, args)

def main():
    (opts, filenames) = process_args()
    if opts.output:
	outf = open(opts.output, "w")
	if opts.verbose:
	    sys.stderr.write("Output to %s\n" % opts.output)
    else:
	outf = sys.stdout
	if opts.verbose:
	    sys.stderr.write("Output to STDOUT\n")
    if opts.output_format == "make":
        outf.write("# DO NOT EDIT --- auto-generated file\n")
        for filename in filenames:
            (deps, mods) = get_deps_and_mods(filename, opts)
            if opts.verbose:
                sys.stderr.write("deps: %s\n" % " ".join(deps))
                sys.stderr.write("mods: %s\n" % " ".join(mods))
            if deps:
                write_deps_make(outf, filename, deps, mods)
    elif opts.output_format == "graphviz":
        all_deps = {}
        has_mods = {}
        for filename in filenames:
            (deps, mods) = get_deps_and_mods(filename, opts)
            if opts.verbose:
                sys.stderr.write("processing %s\n" % filename)
                sys.stderr.write("deps: %s\n" % " ".join(deps))
                sys.stderr.write("mods: %s\n" % " ".join(mods))
            dir, file = os.path.split(filename)
            filebase, fileext = os.path.splitext(file)
            if opts.verbose:
                sys.stderr.write("adding deps to hash: ")
                sys.stderr.write("%s: %s\n" % (mod, " ".join(deps)))
            all_deps[filebase] = deps
            if mods:
                has_mods[filebase] = True
            else:
                has_mods[filebase] = False
        if opts.indirect:
            if opts.verbose:
                sys.stderr.write("using indirect dependencies\n")
            remove_with_indirect(all_deps, opts)
        write_deps_graphviz(outf, filename, all_deps, has_mods)
    else:
        raise Exception("unknown output_format: %s" % opts.output_format)
    outf.close()

if __name__ == "__main__":
    main()
