CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lgsl -lgslcblas

OBJECTS = H3main.o
PROGRAM = ./Executable_files/QS

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

