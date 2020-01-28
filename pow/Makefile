CFLAGS = -std=c++11 -g -Wall -O2

all: pow powthr

pow : pow.o barriertools.o
	g++ -o pow $?

powthr : powthr.o barriertools.o
	g++ -o powthr $?

barriertools.o : barriertools.cpp barriertools.h
	g++ -c -o $@  $< $(CFLAGS)

pow.o:	pow.cpp
	g++ -c -o $@ $< $(CFLAGS)

powthr.o:  powthr.cpp
	g++ -c -o $@ $< $(CFLAGS)

clean:
	rm -rf *.o pow
