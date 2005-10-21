#C***************> /dev/null*****************************************
g77 -o sect_coag sect_coag.f
./sect_coag  		                     <<ein
1		collision kernel: 0 (long), 1 (hall), 2 (golovin)
1		time step	(sec)
10		mean initial radius (um)
4.1886		init water content (g/m**3)
100		eps in cm2 s-3
2.5       	u'	(m/s)
ein
mv sectionm sedi_mt
mv sectionn sedi_nt
#C*****************************************************************
