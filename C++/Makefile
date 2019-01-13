CFLAGS=		-g -w -O2
PROG=		mutation_finder
LIBS=		-lm -lz -lpthread
OBJS=       utils.o sketch.o

all:main.cpp $(OBJS)
	g++ $(CFLAGS) main.cpp -o $(PROG) $(OBJS)

clean:
	rm *.o
	rm $(PROG)

utils:utils.cpp utils.h
	g++ $(CFLAGS) -c utils.cpp

sketch:sketch.cpp sketch.h
	g++ $(CFLAGS) -c sketch.cpp

utils.o: utils.h
sketch.o: sketch.h