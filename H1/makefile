CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm

HEADERS = H1lattice.h H1potential.h
OBJECTS = H1lattice.o H1potential.o H1main.o
PROGRAM = MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

