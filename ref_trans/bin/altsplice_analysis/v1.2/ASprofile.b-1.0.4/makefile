CC=gcc
CFLAGS=-O3 #-g
LDIBS=-lm

extract-as: extract-as.c lib.c
	$(CC) -o extract-as -I ./ $(CFLAGS) extract-as.c lib.c $(LDLIBS)

extract-as-fpkm: extract-as-fpkm.c lib.c
	$(CC) -o extract-as-fpkm -I ./ $(CFLAGS) extract-as-fpkm.c lib.c $(LDLIBS)

all: extract-as extract-as-fpkm

clean:
	rm -f extract-as extract-as-fpkm *.o
