
GRAPHVIZ_FILE = partmc_modules.gv
GRAPHVIZ_PDF = $(patsubst %.gv,%.pdf,$(GRAPHVIZ_FILE))
DOC_INDEX = html/index.html

ALL_SOURCE = $(wildcard ../src/*.f90)

PARTMC_SOURCE = $(filter-out ../src/extract_% ../src/numeric_%,$(ALL_SOURCE))

.PHONY: all
all: doc

.PHONY: doc
doc: README.html $(GRAPHVIZ_PDF) $(DOC_INDEX)

README.html: ../README markdown2.py
	echo "<html>" > $@
	echo "<head><title>PartMC $(VERSION)</title></head>" >> $@
	echo "<body bgcolor=\"#ffffff\">" >> $@
	./markdown2.py ../README >> $@
	echo "</body></html>" >> $@

$(DOC_INDEX): $(ALL_SOURCE) $(GRAPHVIZ_PDF) Doxyfile
	doxygen

$(GRAPHVIZ_FILE): $(PARTMC_SOURCE) f90_mod_deps.py
	echo "digraph partmc_modules {" > $@
	echo "    rankdir = TB;" >> $@
	echo "    node [fontsize = 10, height = 0.3, width = 0.5];" >> $@
	echo "    graph [nodesep = 0.2, ranksep = 0.3];" >> $@
	echo "    edge [arrowsize = 0.7];" >> $@
	./f90_mod_deps.py -d "pmc_(.*)" -D "\1" -m "pmc_(.*)" -M "\1" -g -i $(PARTMC_SOURCE) >> $@
	echo "}" >> $@

%.png: %.gv
	dot -Tpng $< > $@

%.eps: %.gv
	dot -Tps $< > $@

%.pdf: %.eps
	epstopdf $<

clean:
	rm -rf html README.html $(GRAPHVIZ_FILE) $(GRAPHVIZ_PDF)
