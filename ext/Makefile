.PHONY: all

all: gipfeli.a FiniteStateEntropy.a
	@true

EXTLIBFLAGS := -Ofast -march=native

FSEFILES := fse_decompress fse_compress fseU16 huf_compress huf_decompress entropy_common
FiniteStateEntropy.a:
	@echo "DOWNLOADING $*"
	@rm -rf FiniteStateEntropy
	@git clone -q https://github.com/Cyan4973/FiniteStateEntropy
	@cd FiniteStateEntropy && git checkout -q 5443d55fed598867ee543a6a83f7cdd1da6d54fe
	@echo "CREATING $@"
	@for f in $(FSEFILES);\
		do cc -c -I./FiniteStateEntropy/lib -O3 -falign-loops=32 -std=c99 $(EXTLIBFLAGS) ./FiniteStateEntropy/lib/$$f.c -o /tmp/$$f.o;\
	done
	@gcc-ar rcs $@ $(patsubst %,/tmp/%.o,$(FSEFILES))

GIPFELIFILES := decompress entropy entropy_code_builder lz77 gipfeli-internal gipfeli_test
gipfeli.a:
	@echo "DOWNLOADING $*"
	@rm -rf gipfeli
	@git clone -q https://github.com/google/gipfeli
	@cd gipfeli && git checkout -q 04fe241e27f6dcfef239afc6c5e3cee0b4d7c333
	@echo "CREATING $@"
	@for f in $(GIPFELIFILES);\
		do g++ -c -I./gipfeli -O3 $(EXTLIBFLAGS) ./gipfeli/$$f.cc -o /tmp/$$f.o;\
	done
	@gcc-ar rcs $@  $(patsubst %,/tmp/%.o,$(GIPFELIFILES))

clean:
	rm -rf gipfeli FiniteStateEntropy gipfeli.a FiniteStateEntropy.a
