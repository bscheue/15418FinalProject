DEBUG=0
CC=gcc
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)
LDFLAGS= -lm

CFILES = rutil.c image.c sim.c
HFILES = rutil.h image.h

all: crun-seq

crun-seq: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o crun-seq $(CFILES) $(LDFLAGS)


clean:
	rm -f crun crun-seq crun-omp
