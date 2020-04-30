DEBUG=0
CC=gcc
OMP=-fopenmp -DOMP
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)
LDFLAGS= -lm

CFILES = cycletimer.c instrument.c rutil.c image.c sim.c crun.c
HFILES = cycletimer.h instrument.h rutil.h image.h

all: crun-seq crun-omp

crun-seq: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o crun-seq $(CFILES) fake_omp.c $(LDFLAGS)

crun-omp: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) $(OMP) -o crun-omp $(CFILES) $(LDFLAGS)

test: $(CFILES) $(HFILES)
	./test.sh

benchmark: $(CFILES) $(HFILES)
	./benchmark.sh

image: $(CFILES) $(HFILES)
	./crun-omp -m 70 -r "image.bmp"

clean:
	rm -rf crun-seq crun-omp test
