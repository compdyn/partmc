#!/bin/bash

( cd movs
    for d in aero_2d_all_pdfs aero_2d_all_no_coag_pdfs aero_2d_water_pdfs aero_2d_water_no_coag_pdfs ; do
	( cd $d
	    for f in *.pdf ; do
		echo $f
		gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r800 -dTextAlphaBits=1 -dGraphicsAlphaBits=1 -sOutputFile=t.png $f
		# using AlphaBits=4 directly gives hairline cracks
		convert -resize 25% t.png ../${d/_pdfs/_pngs}/${f/.pdf/.png}
		rm t.png
	    done )
    done )
