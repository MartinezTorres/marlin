CFLAGS += -std=gnu++17 -Wall -Wextra -Wcast-qual -Wcast-align -Wstrict-aliasing=1 -Wswitch-enum -Wundef -pedantic  -Wfatal-errors -Wshadow

CFLAGS += -I./src

CFLAGS += `pkg-config opencv --cflags`
LFLAGS += `pkg-config opencv --libs`

LFLAGS += -lboost_system -lboost_program_options -lboost_serialization
LFLAGS += -lz -lrt -fPIC

LCODECS += -lsnappy -lCharLS -lzstd -llz4 -llzo2


CFLAGS += -fopenmp
LFLAGS += -lgomp

CFLAGS += -g
CFLAGS += -g -Ofast
CFLAGS += -march=native

CFLAGS += -I./inc
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


CXX = g++-7
#CXX = clang++-3.3 -D__extern_always_inline=inline -fslp-vectorize
#CXX = icpc -fast -auto-ilp32 -xHost -fopenmp

all: lib/libmarlin.so ./bin/benchmark

.PHONY: all clean

./lib/lib%.so: ./src/%.cc ./inc/*.h
	@echo "CREATING $@" ext
	@$(CXX)  -shared -o $@ $< $(CFLAGS) $(LFLAGS)

./bin/%: ./src/%.cc ./inc/*.h
	@echo "CREATING $@" ext
	@$(CXX) -o $@ $< $(CFLAGS) $(LFLAGS)

clean:
	rm -f bin/*
	rm -f lib/*
