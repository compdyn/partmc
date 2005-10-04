
OPTS = -O -W -Wall -Wimplicit -Wsurprising

MC_SEDI_FIX_OBJS = mc_sedi_fix.o coag_sedi.o part_array.o

PROGS = mc_sedi_fix

all: $(PROGS)

mc_sedi_fix: $(MC_SEDI_FIX_OBJS)
	g77 $(OPTS) -o mc_sedi_fix $(MC_SEDI_FIX_OBJS)

%.o: %.f
	g77 $(OPTS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o *~
