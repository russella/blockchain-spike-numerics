CFLAGS = -g -Wall

all: pos posthr

pos : pos.o barriertools.o
	g++ -o pos $?

posthr : posthr.o barriertools.o
	g++ -o posthr $?

barriertools.o : barriertools.cpp barriertools.h
	g++ -c -o $@  $< $(CFLAGS)

pos.o:	pos.cpp
	g++ -c -o $@ $< $(CFLAGS)

posthr.o:  posthr.cpp
	g++ -c -o $@ $< $(CFLAGS)

clean:
	rm -rf *.o pos
