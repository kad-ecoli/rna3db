CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROG=NaTorsion MissingRNAatom BPtorsion cssr newChainID
HEADER=PDBParser.hpp pstream.h GeometryTools.hpp

all: ${PROG}

newChainID: newChainID.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

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
