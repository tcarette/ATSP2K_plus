#
# Makefile for MCHF Atomic Structure Package
# Created by C. Froese Fischer, December 1988
#  ....  Make the ANGular library

LIB = ../../bin

#  ....  Compiler for SUN, DEC, SGI, ..
#FC = f77

#  ....  Compiler for IBM AIX
#FC = xlf

#  .... Compiler for Portland Group
#FC = pgf77  -O2 -byteswapio

#  ....  Flags to be added for compilation SGI
#FFLAGS =  -O  -static 

# ....  Flags to be added for compilation on SUN, IBM AIX
#FFLAGS = -O
#PFLAGS = -O2 

#  .... Flags to be added for compilaton on DEC ALpha
#FFLAGS = -O4 -r4 -non_shared -static -align dcommon  -convert big_endian
#PFLAGS = -O2 -r4 -non_shared -static -align dcommon  -convert big_endian

#  ....  Compiler for IBM AIX
FC = g95

#    ... Definitions for the ANGular library

ANGSRC  = analy1.f analy2.f cfgin2.f cfgn1.f cfgo1.f cfgtst.f \
          cfp.f cfpd.f cfpf.f cfpp.f find.f gstate.f initt.f \
	  mumdad.f ntab1.f orth.f qsort.f save.f \
	  tensor.f tritst.f vijout.f

ANGOBJS = analy1.o analy2.o cfgin2.o cfgn1.o cfgo1.o cfgtst.o \
          cfp.o cfpd.o cfpf.o cfpp.o find.o gstate.o initt.o \
	  mumdad.o ntab1.o orth.o qsort.o save.o \
	  tensor.o tritst.o vijout.o 

# libang.a -- make the archive library

libang.a : $(ANGOBJS) NJGRAF.o
	ar r libang.a $(ANGOBJS) NJGRAF.o
	mv libang.a $(LIB)

