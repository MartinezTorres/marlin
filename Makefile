CFLAGS += -std=gnu++11 -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing=1 -Wswitch-enum -Wundef -pedantic  -Wfatal-errors

CFLAGS += -I./src

CFLAGS += `pkg-config opencv --cflags`
LFLAGS += `pkg-config opencv --libs`

LFLAGS += -lboost_system -lboost_program_options -lboost_serialization
LFLAGS += -lz -lrt
LFLAGS += -lsnappy -lCharLS -lzstd -llz4 -llzo2


CFLAGS += -fopenmp
LFLAGS += -lgomp

CFLAGS += -Ofast

CFLAGS += -g

CFLAGS += -I./ext
LFLAGS += $(wildcard ./ext/*.a)


#CFLAGS += -DNDEBUG
#CFLAGS += -frename-registers -fopenmp
#CFLAGS += -fno-unroll-loops
#CFLAGS += -funroll-all-loops
#CFLAGS += -fno-align-loops
#CFLAGS += -fno-align-labels
#CFLAGS += -fno-tree-vectorize
#CFLAGS += -falign-functions -falign-labels -falign-jumps -falign-loops -frename-registers -finline-functions
#CFLAGS += -fomit-frame-pointer
#CFLAGS += -fmerge-all-constants -fmodulo-sched -fmodulo-sched-allow-regmoves -funsafe-loop-optimizations -floop-unroll-and-jam

CODECS :=  $(patsubst %.cc,%.o,$(wildcard ./src/codecs/*.cc))


CXX = g++
#CXX = clang++-3.3 -D__extern_always_inline=inline -fslp-vectorize
#CXX = icpc -fast -auto-ilp32 -xHost -fopenmp

all: data ./bin/dcc2017 ./bin/dcc2018

.PHONY: data ext show prof clean realclean

.SECONDARY: $(CODECS)

data:
	@$(MAKE) -C rawzor --no-print-directory
	
ext:
	@$(MAKE) -C ext --no-print-directory
	
./src/codecs/%.o: ./src/codecs/%.cc ./src/codecs/%.hpp ./src/util/*.hpp 
	@echo "CREATING $@"
	@$(CXX) -c -o $@ $< $(CFLAGS)

./bin/dcc2018: ./src/dcc2018.cc $(CODECS) ext
	@echo "CREATING $@"
	@$(CXX) -o $@ $< $(CODECS) $(CFLAGS) $(LFLAGS)

./bin/%: ./src/%.cc
	@echo "CREATING $@"
	@$(CXX) -o $@ $< $(CODECS) $(CFLAGS) $(LFLAGS)

show: ./bin/dcc2017
	./bin/dcc2017
	pdflatex out.tex > /dev/null
	evince out.pdf

prof: ./bin/dcc2017
	 valgrind --dsymutil=yes --cache-sim=yes --branch-sim=yes --dump-instr=yes --trace-jump=no --tool=callgrind --callgrind-out-file=callgrind.out ./eval 
	 kcachegrind callgrind.out

clean:
	rm -f $(CODECS) eval out.tex out.aux out.log out.pdf callgrind.out bin/*

realclean:
	@make -C rawzor clean
	@make -C ext clean
	rm -f $(CODECS) eval out.tex out.aux out.log out.pdf callgrind.out bin/*
