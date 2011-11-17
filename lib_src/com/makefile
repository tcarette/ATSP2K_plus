#
# Makefile for MCHF Atomic Structure Package
# Created by C. Froese Fischer, December 1988
#  ....  Make the COMmon library

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

COMSRC  = acurat.f bisect.f bwint.f ca.f cb.f eptr.f factrl.f fk.f gk.f\
	  gracah.f hnorm.f hwf.f inita.f initr.f intact.f lineqn.f\
	  lval.f pack.f reform.f rk.f rme.f sn.f symmeq.f\
          tinvit.f tql2.f trbak1.f tred1.f tred2.f vk.f

COMOBJS = acurat.o bisect.o bwint.o ca.o cb.o eptr.o factrl.o fk.o gk.o\
          gracah.o hnorm.o hwf.o inita.o initr.o intact.o lineqn.o\
          lval.o pack.o reform.o rk.o rme.o sn.o symmeq.o\
          tinvit.o tql2.o trbak1.o tred1.o tred2.o vk.o


# make the archive library

libcom.a : $(COMOBJS)
	ar r libcom.a $(COMOBJS)
	mv libcom.a $(LIB)


