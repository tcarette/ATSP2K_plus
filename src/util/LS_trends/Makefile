#  the following compilers have been tested: 
# g++ (gcc version 2.95.2)  www.gnu.org/software/gcc/gcc.html
# KCC                       www.kai.com
#
# email questions to;
# georgio@vuse.vanderbilt.edu
# cff@vuse.vanderbilt.edu
# look for "which g++" in order to find the compiler
# test the g++ version: "g++ -v"
# the version should be above 2.95.2 in order to compile
# use "setenv CPP_COMP path_to_g++"
# set the flags: "setenv CPPFLAGS "-O2 -static"
#


EXED=${ATSP}/bin
EXE = ${EXED}/LS_trends
#EXE = LS_trends

SRC = LS_trends.cc rp.cc repr.cc ls.cc 

OBJ = LS_trends.o  rp.o  repr.o  ls.o

all:  ${EXE}

${EXE} : ${OBJ}
	${CPP} ${CPP_LD} -o ${EXE} ${OBJ} -lm

%.o : %.cc
	${CPP} ${CPP_FLAGS} -c ${SRC}

clean:
	rm -f ${EXE} ${OBJ} core 

tar:
	tar -cvf LS.tar ${SRC} *.hh Makefile
