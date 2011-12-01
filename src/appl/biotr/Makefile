################################################################################
#
#  This program evaluates the gf values and transition probabilities
#  of electric and magnetic transitions in the LS coupling scheme
#  (length and velocity forms for E1/E2) or LSJ Breit-Pauli intermediate
#  coupling scheme (length form only).
#
#  M. R. Godefroid    Laboratoire de Chimie Physique Moleculaire
#                     Universite Libre de Bruxelles, Belgium
#  A. Hibbert         Department of Applied Mathematics
#                     Queen's University, Belfast, Northern Ireland
#  C. Froese Fischer  Department of Computer Science
#                     Vanderbilt University, Nashville, U.S.A.
#  P. Jonsson         Department of Physics,
#                     Lund Institute of Technology, Lund, Sweden
#  J. Olsen           Department of Theoretical Chemistry,
#                     Chemical Center, Lund, Sweden
#
#
#  References:
#  ----------
#  1. A. Hibbert et al, Comput. Phys. Commun. 51(1988)285
#  2. C. Froese Fischer et al, Comput. Phys. Commun. 64(1991)486-500
#  3. C. Froese Fischer and M. Godefroid, Comput. Phys. Commun. 64(1991)
#     501-519
#  4. P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)
#  5. J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and
#     C. Froese Fischer, Phys. Rev. A, submitted. (#AY5043)
#
#                      C O P Y R I G H T   2000
################################################################################


################################################################################
# makefile for biotr
# Targets:
# . default: default target is biotr
# . clean:   remove object files 
# . backup: makes a tar file, the date is incorporated into the filename
# . BIOTR.f: will concantenate all source files into BITOR.f
################################################################################
# directory where the executable biotr will be stored 
EXED = ${ATSP}/bin

# executable :: biotr
EXE = $(EXED)/biotr

# make a tar file of all sources + Makefile
TAR = biotr_`date '+%y%m%d_%H'`.bak.tar 

#library directory
LIBDIR = ${ATSP}/lib

# library use
RCA =  -lang -lrad -lcom  

# targets
install: $(EXE) 

# objects
OBJ1 = biotrn.o\
	almult.o analy1.o analy2.o angmom.o bndinv.o brkt.o calcul.o \
	cfgin2.o cfgn1.o cfgo1.o citra.o copvec.o find.o fline.o  \
	gets.o grad2.o gstate.o gtrac1.o gtraca.o gtracb.o gtracxvn.o ielsum.o \
	ifnmnx.o initm2.o inprod.o invmat.o lmatrix.o lmatrix1.o lmatrix2.o \
	lulu.o matml4.o nonh1.o nontrans.o nortbp.o nortbpn.o orth.o orthogg.o \
	outls.o ovlset.o pamtmt.o polint.o probab.o prsym.o radint.o rasin.o \
	rconst.o readw2.o rmetr.o rscheck.o savels.o savenon.o scalve.o \
	setsupras.o setup.o setupm.o setvec.o slope.o ti1tv.o tiini.o \
	transition.o trp.o trpmat.o ulla.o vecsum.o wrtmat.o biotrn_mem.o 
OBJ = $(OBJ1) eigvec.o rdegvc.o biotr.o 

# making
$(EXE):  $(OBJ)
	$(FC) -o $(EXE) $(OBJ) $(FC_LD)  -L$(LIBDIR) $(RCA) \
                 ${LAPACK_DIR} ${LAPACK_LIB}
backup : 
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core biotr

# source files
BIOTRS = \
biotr.f       biotrn.f \
almult.f      analy1.f      analy2.f      angmom.f      biotrn_mem.f  \
bndinv.f      brkt.f        calcul.f      cfgin2.f      cfgn1.f       \
cfgo1.f       citra.f       copvec.f      eigvec.f      find.f        \
fline.f       gets.f        grad2.f       gstate.f      gtrac1.f      \
gtraca.f      gtracb.f      gtracxvn.f    ielsum.f      ifnmnx.f      \
initm2.f      inprod.f      invmat.f      lmatrix.f     lmatrix1.f    \
lmatrix2.f    lulu.f        matml4.f      nonh1.f       nontrans.f    \
nortbp.f      nortbpn.f     orth.f        orthogg.f     outls.f       \
ovlset.f      pamtmt.f      polint.f      probab.f      prsym.f       \
radint.f      rasin.f       rconst.f      rdegvc.f      readw2.f      \
rmetr.f       rscheck.f     savels.f \
savenon.f     scalve.f      setsupras.f   setup.f       setupm.f      \
setvec.f      slope.f       ti1tv.f       tiini.f   \
transition.f  trp.f         trpmat.f      ulla.f        vecsum.f      \
wrtmat.f

BIOTR: $(BIOTRS)
	cat $(BIOTRS) > BIOTR

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

