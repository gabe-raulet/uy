ifeq ($(shell uname), Darwin)
CC=gcc-11
else
CC=cc
endif

CFLAGS=-Wall -Werror -Wno-unknown-pragmas -fopenmp
DEBUG?=0
THREADED?=0
LOGGER?=1

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

PRGS=bfs

all: $(PRGS)

bfs: bfs.c uy.o
	$(CC) $(CFLAGS) -o $@ $^

uy.o: uy.c uy.h
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: all clean purge

clean:
	@rm -rf *.o *.dSYM

purge: clean
	@rm -rf $(PRGS)
