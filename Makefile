
CFLAGS=-Wall -Werror -Wno-unknown-pragmas -fopenmp
DEBUG?=0
THREADED?=1
LOGGER?=1

ifeq ($(shell uname), Darwin)
CC=gcc-11
else
CC=cc
CFLAGS+=-lm
endif

ifeq ($(DEBUG), 1)
CFLAGS+=-DDEBUG -g -O0 -fsanitize=address -fno-omit-frame-pointer
else
CFLAGS+=-O2
endif

ifeq ($(THREADED), 1)
CFLAGS+=-DTHREADED
endif

ifeq ($(LOGGER), 1)
CFLAGS+=-DLOGGER
endif

PRGS=bfs uy spadd spgemm

all: $(PRGS)

bfs: bfs.c spmat.o
	$(CC) $(CFLAGS) -o $@ $^

uy: uy.c spmat.o
	$(CC) $(CFLAGS) -o $@ $^

spadd: spadd.c spmat.o
	$(CC) $(CFLAGS) -o $@ $^

spgemm: spgemm.c spmat.o
	$(CC) $(CFLAGS) -o $@ $^

spmat.o: spmat.c spmat.h
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: all clean purge

clean:
	@rm -rf *.o *.dSYM

purge: clean
	@rm -rf $(PRGS)
