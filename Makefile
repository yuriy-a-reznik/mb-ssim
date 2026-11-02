# nbssim -- N-band filterbank SSIM.  Strict C89 build.

CC      ?= cc
CSTD    ?= -std=c89
WARN     = -Wall -Wextra -pedantic -Wno-long-long
OPT     ?= -O2
CFLAGS  += $(CSTD) $(WARN) $(OPT) -Iinclude -MMD -MP
LDLIBS  += -lm

SRC  = src/core.c src/fft.c src/image.c src/png.c src/ssim.c src/msssim.c \
       src/nbssim.c src/nbssim_dft.c src/stats.c src/data.c src/cache.c src/cli.c
OBJ  = $(SRC:.c=.o)

TOOLS = bin/nbsmetric bin/nbsbench bin/nbstune bin/nbsjoint bin/nbsbank bin/nbscachesel
TESTS = bin/nbstest

all: $(TOOLS) $(TESTS)

bin:
	mkdir -p bin

bin/nbsmetric: tools/nbsmetric.o $(OBJ) | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/nbsbench: tools/nbsbench.o $(OBJ) | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/nbstune: tools/nbstune.o $(OBJ) | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/nbsjoint: tools/nbsjoint.o $(OBJ) | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/nbsbank: tools/nbsbank.o $(OBJ) | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/nbscachesel: tools/nbscachesel.o $(OBJ) | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

bin/nbstest: test/nbstest.o $(OBJ) | bin
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

check: bin/nbstest
	./bin/nbstest

clean:
	rm -f $(OBJ) tools/*.o test/*.o src/*.d tools/*.d test/*.d
	rm -rf bin

-include $(OBJ:.o=.d) tools/nbsmetric.d tools/nbsbench.d tools/nbstune.d tools/nbsjoint.d tools/nbsbank.d test/nbstest.d

.PHONY: all check clean
