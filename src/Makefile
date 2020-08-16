CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROG=NaTorsion MissingRNAatom BPtorsion cssr
HEADER=PDBParser.hpp pstream.h GeometryTools.hpp

all: ${PROG}

NaTorsion: NaTorsion.cpp NaTorsion.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

BPtorsion: BPtorsion.cpp BPtorsion.hpp StringTools.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cssr: cssr.cpp cssr.hpp BPstat.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

MissingRNAatom: MissingRNAatom.cpp MissingRNAatom.hpp Superpose.hpp IdealRNA.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: ${PROG}
	cp ${PROG} ../script

clean:
	rm ${PROG}
