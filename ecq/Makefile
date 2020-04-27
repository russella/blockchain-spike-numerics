CFLAGS = -std=c++11 -g -Wall -O3

all: ecq

ecq: disttools.o ecq.o 
	g++ -o ecq $^

ecq-pg: disttools.o ecq-pg.o 
	g++ -o ecq-pg $^

disttools.o: disttools.cpp disttools.h
	g++ -c -o $@  $< $(CFLAGS)

ecq.o:	ecq.cpp disttools.h
	g++ -c -o $@ $< $(CFLAGS)

ecq-pg.o:	ecq-pg.cpp disttools.h
	g++ -c -o $@ $< $(CFLAGS)

clean:
	rm -rf *.o ecq
