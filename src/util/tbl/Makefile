#  the following compilers have been tested: 
# g++ (gcc version 3.3.1)  www.gnu.org/software/gcc/gcc.html
# KCC                       www.kai.com
#
# email questions to;
# georgio@vuse.vanderbilt.edu
# cff@vuse.vanderbilt.edu
# look for "which g++" in order to find the compiler
# test the g++ version: "g++ -v"
# the version should be above 2.95.2 in order to compile
# use "setenv CPP path_to_g++"
# set the flags: "setenv CPPFLAGS "-O2 -static"
#

#CPP=/usr/bin/g++
#CPP_FLAGS=-static -O2 

EXE1 = $(ATSP)/bin/tables

SRC = _inp.cc _prn.cc _tbl.cc _tt.cc _tr.cc _to.cc _te.cc\
      appendOdd.cc

O1 = _tbl.o _prn.o _inp.o _tt.o _tr.o _to.o _te.o \
     appendOdd.o

all:  ${EXE1}
${EXE} : ${OBJ}
	${CPP} ${CPP_LD} -o ${EXE} ${OBJ}

${EXE1} : ${O1}
	${CPP} ${CPP_LD} -o ${EXE1} ${O1} -lm

%.o : %.cc
	${CPP} ${CPP_FLAGS} -c ${SRC}

clean:
	rm -f ${EXE} ${O1} core 

tar:
	tar -cvf tbl.tar ${SRC} *.hh Makefile
