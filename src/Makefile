CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROG=NaTorsion NaTorsion2 MissingRNAatom BPtorsion cssr newChainID NAbuilder AllDist MissingResidue
HEADER=PDBParser.hpp pstream.h GeometryTools.hpp

all: ${PROG}

AllDist: AllDist.cpp AllDist.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

newChainID: newChainID.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NAbuilder: NAbuilder.cpp NAbuilder.hpp pseudo_torsion_prob.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NaTorsion: NaTorsion.cpp NaTorsion.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NaTorsion2: NaTorsion2.cpp NaTorsion2.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

BPtorsion: BPtorsion.cpp BPtorsion.hpp StringTools.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cssr: cssr.cpp cssr.hpp BPstat.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

MissingRNAatom: MissingRNAatom.cpp MissingRNAatom.hpp Superpose.hpp IdealRNA.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

MissingResidue: MissingResidue.cpp MissingResidue.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: ${PROG}
	cp ${PROG} ../script

clean:
	rm ${PROG}
