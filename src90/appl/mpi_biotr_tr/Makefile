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

.SUFFIXES: .f90 .mod

# directory where the executable biotr will be stored 
EXED = ${ATSP}/bin
#EXED = .

#FC = xlf95
#FC_FLAGS = -g -qsuffix=f=f90 
LFLAGS = 

# executable :: biotr_tr
EXEF = $(EXED)/biotr_tr_mpi

# library modules
MOD = ../../lib/libmod

# make a tar file of all sources + Makefile
TAR = biotr_tr_mpi_`date '+%y%m%d_%H'`.bak.tar 

#library directory
LIBDIR = ${ATSP}/lib
#LIBDIR = ${ATSP}/lib

# library use
RCA = -lang90 -lrad90 -lcom90  -lmod90  -lmpiu90

# targets
install: $(EXEF) 

# objects
BIOTR =  mpi_biotr_tr.o inp_biotr_tr.o mpi_read_data.o readw2.o \
         initm2.o brkt.o \
         eigvec.o rdegvc.o almult.o slope.o polint.o biotrn.o biotrn_mem.o \
         ifnmnx.o citra.o radint.o mpi_calcul.o probab.o nortbp.o nontrans.o \
         setup_oc.o trp.o transition.o gets.o scalve.o ielsum.o copvec.o \
         invmat.o ulla.o trpmat.o matml4.o setvec.o pamtmt.o vecsum.o \
         wrtmat.o bndinv.o lulu.o inprod.o prsym.o grad2.o fline.o  \
         tiini.o ti1tv.o gtrac1.o gtracxvn.o setupm.o savels.o


OALL = $(BIOTR) 

# making
$(EXEF):  $(BIOTR)
	$(FC_MPI) -o $(EXEF) $(FC_MPIFAGS) ${OALL}  -L$(LIBDIR) $(RCA) 

backup : 
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core biotr_tr *.mod

SRC = \
biotr_I.f biotr_C.f mpi_biotr.f inp_biotr_tr.f read_data.f readw2.f  \
initm2.f brkt.f eigvec.f rdegvc.f almult.f slope.f polint.f biotrn.f \
biotrn_mem.f ifnmnx.f citra.f radint.f calcul.f probab.f nortbp.f \
nontrans.f setup_oc.f trp.f transition.f gets.f scalve.f ielsum.f \
copvec.f invmat.f ulla.f trpmat.f matml4.f setvec.f pamtmt.f vecsum.f \
wrtmat.f bndinv.f lulu.f inprod.f prsym.f grad2.f fline.f tiini.f  \
ti1tv.f gtrac1.f gtracxvn.f setupm.f savels.f
.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I ${MOD} -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

