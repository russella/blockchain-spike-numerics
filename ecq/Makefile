CFLAGS = -std=c++11 -g -Wall -O2

all: ecq

ecq : ecq.o disttools.o
	g++ -o ecq $?

disttools.o : disttools.cpp disttools.h
	g++ -c -o $@  $< $(CFLAGS)

ecq.o:	ecq.cpp
	g++ -c -o $@ $< $(CFLAGS)

clean:
	rm -rf *.o ecq
