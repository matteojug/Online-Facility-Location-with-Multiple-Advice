CFLAGS=-O3 -std=c++1z mpi-sw/facloc_JMS/facloc_JMS.c mpi-sw/facloc_JMS/basic_calls.c mpi-sw/facloc_LOCAL/facloc_LOCAL.c
# CFLAGS+= -g
CFLAGS+= -DHUGE=1e15 # mpi def
CFLAGS+= -DINSTANCE_LIMIT=45000 # offline size limit
# CFLAGS+= -DUSE_PTS_GEN # Uncomment for >2 dimensions

.PHONY: clean offline online synth

all: offline online synth

offline:
	g++ $(CFLAGS) solve_offline.cpp -o solve_offline

online:
	g++ $(CFLAGS) solve_online.cpp -o solve_online

synth:
	g++ $(CFLAGS) solve_synth.cpp -o solve_synth
	
clean:
	rm solve_offline solve_online solve_synth