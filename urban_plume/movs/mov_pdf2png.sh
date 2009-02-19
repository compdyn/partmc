#!/bin/bash

for d in aero_2d_all_pdfs aero_2d_all_no_coag_pdfs aero_2d_water_pdfs aero_2d_water_no_coag_pdfs ; do
    ( cd $d
	for f in *.pdf ; do
	    echo $f
	    gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r200 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sOutputFile=../${d/_pdfs/_pngs}/${f/.pdf/.png} $f
	done )
done
