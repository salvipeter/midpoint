all: example

INCLUDES=-I../libgeom
LIBS=-L../libgeom/release -lgeom
CXXFLAGS=-Wall -pedantic -std=c++17 $(INCLUDES)

example: example.o midpoint.o
	g++ -o $@ $^ $(LIBS)

c0coons.o: midpoint.cc midpoint.hh
