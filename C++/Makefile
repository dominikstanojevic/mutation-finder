GCC=        g++
CFLAGS=		-g -w -Ofast -fopenmp -march=native
PROG=		mutation_finder
LIBS=		-lm -lz -lpthread
OBJS=       utils.o sketch.o map.o align.o

all: main.cpp $(OBJS) align.h map.h sketch.h utils.h
	$(GCC) $(CFLAGS) main.cpp -o $(PROG) $(OBJS)

clean:
	rm *.o
	rm $(PROG)

utils.o :utils.cpp utils.h
	$(GCC) $(CFLAGS) -c utils.cpp

sketch.o :sketch.cpp sketch.h
	$(GCC) $(CFLAGS) -c sketch.cpp

map.o :map.cpp map.h
	$(GCC) $(CFLAGS) -c map.cpp

align.o :align.cpp align.h
	$(GCC) $(CFLAGS) -c align.cpp