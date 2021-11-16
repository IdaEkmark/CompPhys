CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lgsl -lgslcblas

HEADERS = H1lattice.h H1potential.h
OBJECTS = H1lattice.o H1potential.o H1main.o
PROGRAM = ./Executable_files/MD

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

