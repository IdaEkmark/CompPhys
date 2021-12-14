CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lgsl -lgslcblas

OBJECTS = H2main.o
PROGRAM = ./Executable_files/MC

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

