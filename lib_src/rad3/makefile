#
# Makefile for MCHF Atomic Structure Package
# Created by C. Froese Fischer, December 1988
#  ....  Make the Rad3 library

#  ....  Location of libraries and binaries
LIB = ../../bin

#  ....  Compiler for SUN, DEC, SGI, ..
#FC = f77

#  ....  Compiler for IBM AIX
#FC = xlf

#  .... Compiler for Portland Group
#FC = pgf77 -O2 -byteswapio

#  ....  Flags to be added for compilation SGI
#FFLAGS =  -O  -static 

# ....  Flags to be added for compilation on SUN, IBM AIX
#FFLAGS = -O

#  .... Flags to be added for compilaton on DEC ALpha
#FFLAGS = -O4 -r4 -non_shared -static -align dcommon  -convert big_endian

#  ....  Compiler for IBM AIX
FC = g95

#  ....  Make the RADial library for 30 orbitals

#   ...  Definitions for the RAD Libraries

RADSRC =  dyk.f ecore.f grad.f hl.f hlc.f quadr.f quads.f\
	  rlshft.f ykf.f zeta.f zk.f

RADOBJS = dyk.o ecore.o grad.o hl.o hlc.o quadr.o quads.o\
	  rlshft.o ykf.o zeta.o zk.o

#
librad3.a : $(RADOBJS)
	 ar r librad3.a $(RADOBJS)
	 mv librad3.a $(LIB)

