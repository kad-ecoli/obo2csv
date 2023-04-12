CC=g++
CFLAGS=-O3
LDFLAGS=-static

all: obo2csv backpropagate curate_GAF

obo2csv: obo2csv.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

backpropagate: backpropagate.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

curate_GAF: curate_GAF.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
clean:
	rm obo2csv backpropagate curate_GAF
