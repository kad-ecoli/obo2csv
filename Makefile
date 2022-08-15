CC=g++
CFLAGS=-O3
LDFLAGS=-static

obo2csv: obo2csv.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

