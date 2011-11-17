################################################################################
#YTRANS
# Last modifications : September 14, 1997
# tests to do : E2/M1
#               lhs = rhs
#*************************************************************************
#
#  This program evaluates the gf values and transition probabilities
#  of electric and magnetic transitions in the LS coupling scheme
#  (length and velocity forms for E1/E2) or LSJ intermediate coupling
#  scheme (length form only).
#
#  M. Godefroid       Laboratoire de Chimie Physique Moleculaire CP160/09
#                     Universite Libre de Bruxelles, Belgium
#  A. Hibbert         Department of Applied Mathematics
#                     Queen's University, Belfast, Northern Ireland
#  C. Froese Fischer  Department of Computer Science
#                     Vanderbilt University, Nashville, U.S.A.
#  G. Gaigalas        Institute of Theoretical Physiscs
#                     and Astronomy, Vilnius, LITHUANIA
#
#  References:
#  ----------
#  1. Hibbert et al,                Comput. Phys. Commun. 51(1988)285
#  2. Froese Fischer et al,         Comput. Phys. Commun. 64(1991)486-500
#  3. Froese Fischer and Godefroid, Comput. Phys. Commun. 64(1991)501-519
#
################################################################################



################################################################################
# makefile for nonh
# Targets:
# . default: default target is trans
# . clean:   remove object files
# . backup: makes a tar file, the date is incorporate in filename
# . TRANS.f: will concantenate all source files into TRANS.f
################################################################################

EXED = ${ATSP}/bin

TAR = trans_`date '+%y%m%d_%H'`.bak.tar   # tar -file : year_month_date_hour
# executable targets
EXEF = $(EXED)/trans

# libraries
LIBDIR = ${ATSP}/lib
LIB_ANG_COM_RAD = -lang -lcom -lrad


# Object files needed for the system
OBJ = trans.o analy2.o almult.o calcul.o cfgin2.o find.o fline.o      \
      grad2.o gstate.o initm2.o orth.o radint.o rconst.o  readw2.o    \
      rmetr.o setup.o trp.o nontrans.o transition.o savenon.o probab.o \
      valnum.o  take.o 
OBJF = $(OBJ) eigvec.o rdegvc.o 

#makinfg
install : $(EXEU) $(EXEF)
$(EXEF): $(OBJF)
	$(FC) -o $(EXEF) $(FC_LD) -L$(LIBDIR) $(OBJF) $(LIB_ANG_COM_RAD)


backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core trans
TRANS.f  : $(SRCF)
	cat $(SRC) > TRANS.f
.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

