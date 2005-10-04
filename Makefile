
OPTS = -g -W -Wall -Wimplicit -Wsurprising

MC_SEDI_FIX_OBJS = mc_sedi_fix.o kernel_sedi.o particle_array.o
MC_SEDI_INTER_OBJS = mc_sedi_inter.o kernel_sedi.o particle_array.o
UMSP_MC_OBJS = umsp_mc.o

PROGS = mc_sedi_fix mc_sedi_inter umsp_mc

all: $(PROGS)

mc_sedi_fix: $(MC_SEDI_FIX_OBJS)
	g77 $(OPTS) -o $@ $(MC_SEDI_FIX_OBJS)

mc_sedi_inter: $(MC_SEDI_INTER_OBJS)
	g77 $(OPTS) -o $@ $(MC_SEDI_INTER_OBJS)

umsp_mc: $(UMSP_MC_OBJS)
	g77 $(OPTS) -o $@ $(UMSP_MC_OBJS)

%.o: %.f
	g77 $(OPTS) -c -o $@ $<

clean:
	rm -f $(PROGS) *.o *~
