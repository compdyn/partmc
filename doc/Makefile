
VERSION = 2.4.0

GRAPHVIZ_FILE = partmc_modules.gv
GRAPHVIZ_PDF = $(patsubst %.gv,%.pdf,$(GRAPHVIZ_FILE))
DOC_INDEX = html/index.html

ALL_FSRC = $(wildcard ../src/*.F90 ../src/aero_reps/*.F90 ../src/rxns/*.F90 ../src/sub_models/*.F90)
ALL_CSRC = $(wildcard ../src/*.c ../src/aero_reps/*.c ../src/rxns/*.c ../src/sub_models/*.c)
ALL_SRC = $(ALL_FSRC) $(ALL_CSRC)

PARTMC_FSRC = $(filter-out ../src/extract_% ../src/numeric_% ../src/bin_average_% ../src/getopt.%,$(ALL_FSRC))
PARTMC_CSRC = $(ALL_CSRC)
PARTMC_SRC = $(PARTMC_FSRC) $(PARTMC_CSRC)

.PHONY: all
all: doc

.PHONY: doc
doc: README.html ChangeLog.html $(GRAPHVIZ_PDF) $(DOC_INDEX)

README.html: ../README.md markdown2.py
	echo "<html>" > $@
	echo "<head><title>PartMC $(VERSION)</title></head>" >> $@
	echo "<body bgcolor=\"#ffffff\">" >> $@
	./markdown2.py ../README.md >> $@
	echo "</body></html>" >> $@

ChangeLog.html: ../ChangeLog markdown2.py
	echo "<html>" > $@
	echo "<head><title>ChangeLog: PartMC $(VERSION)</title></head>" >> $@
	echo "<body bgcolor=\"#ffffff\">" >> $@
	echo "<h1>ChangeLog: PartMC $(VERSION)</h1>" >> $@
	./markdown2.py ../ChangeLog >> $@
	echo "</body></html>" >> $@

$(DOC_INDEX): $(ALL_SRC) $(GRAPHVIZ_PDF) Doxyfile
	doxygen

$(GRAPHVIZ_FILE): $(PARTMC_FSRC) f90_mod_deps.py
	echo "digraph partmc_modules {" > $@
	echo "    rankdir = TB;" >> $@
	echo "    node [fontsize = 10, height = 0.3, width = 0.5];" >> $@
	echo "    graph [nodesep = 0.2, ranksep = 0.3];" >> $@
	echo "    edge [arrowsize = 0.7];" >> $@
	./f90_mod_deps.py -d "pmc_(.*)" -D "\1" -m "pmc_(.*)" -M "\1" -g -i $(PARTMC_FSRC) >> $@
	echo "}" >> $@

%.png: %.gv
	dot -Tpng $< > $@

%.eps: %.gv
	dot -Tps $< > $@

%.pdf: %.eps
	epstopdf $<

clean:
	rm -rf html latex README.html ChangeLog.html $(GRAPHVIZ_FILE) $(GRAPHVIZ_PDF)
