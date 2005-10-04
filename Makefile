
OPTS = -O -W -Wall -Wimplicit -Wsurprising

MC_SEDI_FIX_OBJS = mc_sedi_fix.o coag_sedi.o part_array.o
MC_SEDI_INTER_OBJS = mc_sedi_inter.o coag_sedi.o part_array.o

PROGS = mc_sedi_fix mc_sedi_inter

all: $(PROGS)

mc_sedi_fix: $(MC_SEDI_FIX_OBJS)
	g77 $(OPTS) -o $@ $(MC_SEDI_FIX_OBJS)

mc_sedi_inter: $(MC_SEDI_INTER_OBJS)
	g77 $(OPTS) -o $@ $(MC_SEDI_INTER_OBJS)

%.o: %.f
	g77 $(OPTS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o *~
