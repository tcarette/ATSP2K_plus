#
# Makefile for MCHF Atomic Structure Package
# Created by C. Froese Fischer, December 1988
#

#  ....   Libraries to be searched
LIB_ANG_COM = $(LIB)/libang.a  -lcom
LIB_COM_RAD3 = -lcom -lrad3 -lcom
LIB_ANG_COM_RAD3 = $(LIB)/libang.a -lcom -lrad3 -lcom
LIB_COM_RAD6 = -lcom -lrad6 -lcom
LIB_ANG_COM_RAD6 = $(LIB)/libang.a -lcom -lrad6 -lcom

#  ....  Location of libraries and binaries
LIB = ../bin

#  ....  PGI Compiler
#FC = pgf77 -O2 -byteswapio

#  ....  Compiler for SUN, DEC, SGI, ..
#FC = f77

#  ....  Compiler for IBM AIX
#FC = xlf

#  ....  Flags to be added for compilation SGI
#FFLAGS =  -O  -static 
#PFLAGS = -O  -static 

# ....  Flags to be added for compilation on SUN, IBM AIX
#FFLAGS = -O2
#PFLAGS = -O2 
#PFLAGS = -g 

#  .... Flags to be added for compilaton on DEC ALpha
#FFLAGS = -O4 -r4 -non_shared -static -align dcommon -assume noacc -convert big_endian
#PFLAGS = -O2 -r4 -non_shared -static -align dcommon -assume noacc -convert big_endian

#  ....  Compiler g95
FC = g95 -fsloppy-char

UTILITIES = Comp Condens Levels Lines Plotw Printw Relabel Terms
PROGRAMS = Auto Breit Cmchf Ci Gencl HF Hfs Iso Lsjtr Lstr Mchf Mltpol Nonh Photo
install:  ${UTILITIES} ${PROGRAMS}
#
#  ...  Make the UTILITIES
#

Comp : COMP.o
	$(FC) -o Comp $(FFLAGS)  COMP.o -L$(LIB) $(LIB_ANG_COM)
	mv Comp $(LIB)

Condens : CONDENS.o
	$(FC) -o Condens $(FFLAGS) CONDENS.o 
	mv Condens $(LIB)

Levels : LEVELS.o
	$(FC) -o Levels $(FFLAGS) LEVELS.o 
	mv Levels $(LIB)

Lines : LINES.o
	$(FC) -o Lines $(FFLAGS) LINES.o 
	mv Lines $(LIB)

Plotw:  PLOTW.o
	$(FC) -o Plotw PLOTW.o
	mv Plotw $(LIB)

Printw: PRINTW.o
	$(FC) -o Printw PRINTW.o
	mv Printw $(LIB)
	
Relabel : RELABEL.o
	$(FC) -o Relabel $(FFLAGS) RELABEL.o 
	mv Relabel $(LIB)

Terms : TERMS.o
	$(FC) -o Terms $(FFLAGS) TERMS.o 
	mv Terms $(LIB)

#
#  .... Make the programs. Because global optimizaiton on long
#  .... programs causes swap space problems, these will be 
#  .... compiled here with a lower level of optimization
#


Auto :  AUTO.o
	$(FC) -o Auto $(PFLAGS) AUTO.o -L$(LIB) $(LIB_COM_RAD3)
	mv Auto $(LIB)

Breit : BREIT.o
	$(FC) -o Breit $(PFLAGS) BREIT.o -L$(LIB) $(LIB_ANG_COM)
	mv Breit $(LIB)

Cmchf : CMCHF.o
	$(FC) -o Cmchf $(PFLAGS) CMCHF.o -L$(LIB) $(LIB_COM_RAD3)
	mv Cmchf $(LIB)
	
Ci :    CI.o
	$(FC) -o Ci $(PFLAGS) CI.o -L$(LIB) $(LIB_COM_RAD3)
	mv Ci $(LIB)

Gencl:  GENCL.o
	$(FC) -o Gencl $(FFLAGS) GENCL.o
	mv Gencl $(LIB)

HF :    HF.o
	$(FC) -o HF $(PFLAGS) HF.o
	mv HF $(LIB)

Hfs :   HFS.o
	$(FC) -o Hfs $(PFLAGS) HFS.o -L$(LIB) $(LIB_ANG_COM_RAD3)
	mv Hfs $(LIB)

Iso :   ISO.o
	$(FC) -o Iso $(PFLAGS) ISO.o -L$(LIB) $(LIB_COM_RAD3)
	mv Iso $(LIB)

Lsjtr : LSJTR.o
	$(FC) -o Lsjtr $(PFLAGS) LSJTR.o -L$(LIB) $(LIB_ANG_COM_RAD6)
	mv Lsjtr $(LIB)

Lstr :  LSTR.o
	$(FC) -o Lstr $(FFLAGS) LSTR.o -L$(LIB) $(LIB_COM_RAD6)
	mv Lstr $(LIB)

Mchf :  MCHF.o
	$(FC) -o Mchf $(PFLAGS) MCHF.o -L$(LIB) $(LIB_COM_RAD3)
	mv Mchf $(LIB)

Mltpol: MLTPOL.o
	$(FC) -o Mltpol $(PFLAGS) MLTPOL.o -L$(LIB) $(LIB_ANG_COM)
	mv Mltpol $(LIB)

Nonh :  NONH.o
	$(FC) -o Nonh $(PFLAGS) NONH.o -L$(LIB) $(LIB_ANG_COM)
	mv Nonh $(LIB)

Photo:  PHOTO.o
	$(FC) -o Photo $(PFLAGS) PHOTO.o -L$(LIB) $(LIB_COM_RAD6)
	mv Photo $(LIB)

