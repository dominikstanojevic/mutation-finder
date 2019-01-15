GCC=        g++
CFLAGS=		-g -w -O2
PROG=		mutation_finder
LIBS=		-lm -lz -lpthread
OBJS=       utils.o sketch.o map.o align.o

all:main.cpp $(OBJS)
	$(GCC) $(CFLAGS) main.cpp -o $(PROG) $(OBJS)

clean:
	rm *.o
	rm $(PROG)

utils:utils.cpp utils.h
	$(GCC) $(CFLAGS) -c utils.cpp

sketch:sketch.cpp sketch.h
	$(GCC) $(CFLAGS) -c sketch.cpp

map:map.cpp map.h
	$(GCC) $(CFLAGS) -c map.cpp

align:align.cpp align.h
	$(GCC) $(CFLAGS) -c align.cpp

utils.o: utils.h
sketch.o: sketch.h
map.o: map.h
align.o: align.h