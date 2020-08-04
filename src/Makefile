CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROG=NaTorsion
HEADER=PDBParser.hpp pstream.h

all: ${PROG}

NaTorsion: NaTorsion.cpp NaTorsion.hpp GeometryTools.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: ${PROG}
	cp ${PROG} ../script

clean:
	rm ${PROG}
