all: example

INCLUDES=-I../libgeom
LIBS=-L../libgeom/release -lgeom
CXXFLAGS=-Wall -pedantic -std=c++17 -O3 $(INCLUDES)

example: example.o midpoint.o regular-domain.o
	g++ -o $@ $^ $(LIBS)

midpoint.o: midpoint.cc midpoint.hh

regular-domain.o: regular-domain.cc regular-domain.hh
