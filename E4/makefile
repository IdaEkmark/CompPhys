CC = gcc
CFLAGS = -O3 -Wall
LIBS = -lm -lgsl -lgslcblas

OBJECTS = E4main.o fft.o
HEADERS = fft.h
PROGRAM = ./Executable_files/BD

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

