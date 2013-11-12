################################################################################
#
#  This program evaluates ...
#
#  M. R. Godefroid    Laboratoire de Chimie Physique Moleculaire
#                     Universite Libre de Bruxelles, Belgium
#  P. Jonsson         ????
#                     ???
#
#  References:
#  ----------
#
#                      C O P Y R I G H T   2008
################################################################################


################################################################################
# makefile for biomat
# Targets:
# . default: default target is biomat
# . clean:   remove object files 
# . backup: makes a tar file, the date is incorporated into the filename
# . BIOMAT.f: will concantenate all source files into BIOMAT.f
################################################################################
# directory where the executable biotr will be stored 
EXED = ${ATSP}/bin

# executable :: biotr
EXE = $(EXED)/biomat

# make a tar file of all sources + Makefile
TAR = biomat_`date '+%y%m%d_%H'`.bak.tar 

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
	transition.o trp.o trpmat.o ulla.o vecsum.o wrtmat.o biotrn_mem.o \
        writerw.o writerl.o brktin.o
OBJ = $(OBJ1) eigvec.o rdegvc.o biomat.o 

# making
$(EXE):  $(OBJ)
	$(FC) -o $(EXE) $(OBJ) $(FC_LD)  -L$(LIBDIR) $(RCA) \
                ${LAPACK_DIR} ${LAPACK_LIB}
backup : 
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core biomat

# source files
BIOMATS = \
biomat.f      biotrn.f \
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
wrtmat.f      writerw.f     writerl.f     brktin.f

BIOMAT: $(BIOMATS)
	cat $(BIOMATS) > BIOMAT

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

