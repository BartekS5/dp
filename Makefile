CC = gcc
CFLAGS = -O2 -Wall

lattice: lattice.c
	$(CC) $(CFLAGS) -o lattice lattice.c

clean:
	rm -f lattice
