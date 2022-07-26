
CFLAGS=-Wall -Werror -Wno-unknown-pragmas -fopenmp
D?=0
T?=0
L?=1
V?=0

ifeq ($(shell uname), Darwin)
CC=gcc-11
else
CC=cc
CFLAGS+=-lm
endif

ifeq ($(D), 1)
CFLAGS+=-DDEBUG -g -O0 -fsanitize=address -fno-omit-frame-pointer
else
CFLAGS+=-O2 -fno-tree-vectorize
endif

ifeq ($(T), 1)
CFLAGS+=-DTHREADED
endif

ifeq ($(L), 1)
CFLAGS+=-DLOGGER
endif

ifeq ($(V), 1)
CFLAGS+=-DVERBOSE
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
