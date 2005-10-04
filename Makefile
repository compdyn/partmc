
OPTS = -Wunused

MC_SEDI_FIX_OBJS = mc_sedi_fix.o coag.o

all: mc_sedi_fix

mc_sedi_fix: $(MC_SEDI_FIX_OBJS)
	g77 $(OPTS) -o mc_sedi_fix $(MC_SEDI_FIX_OBJS)

%.o: %.f
	g77 $(OPTS) -o $@ $<
