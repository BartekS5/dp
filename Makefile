CC = gcc
CFLAGS = -O2 -Wall

run:
	$(CC) $(CFLAGS) lattice.c -o lattice
	./lattice

clean:
	rm -f lattice
