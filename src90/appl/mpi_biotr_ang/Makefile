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

# executable :: biotr
EXEF = $(EXED)/biotr_ang_mpi

# library modules
LMOD = ../../lib/libmod

VO = ${LMOD}/vast_kind_param_M.o

# make a tar file of all sources + Makefile
TAR = biotr_ang_mpi_`date '+%y%m%d_%H'`.bak.tar 

#library directory
LIBDIR = ${ATSP}/lib

# library use
RCA =  -lang90 -lrad90 -lcom90  -lmod90 -lmpiu90

# targets
install: $(EXEF) 

# objects
BIOTR =  \
analy2.o       cfgn1.o     lmatrix1.o       orth.o     rscheck.o  \
mpi_angmom.o       find.o      lmatrix2.o       orthogg.o  savels.o   \
gstate.o       outls.o          setupm.o  mpi_cfgo1.o \
inp_biotr_ang.o     mpi_nonh1.o     rasin.o          setsupras.o  \
mpi_cfgin2.o    lmatrix.o   nortbpn.o   setup.o mpi_save_data.o

BIOTRF = $(BIOTR) mpi_biotr_ang.o 

OALL = mpi_biotr_ang.o $(BIOTR) 

# making
$(EXEF):  $(BIOTRF)
	$(FC_MPI) -o $(EXEF) $(FC_MPILD) ${OALL}  -L$(LIBDIR) $(RCA) 

$(BIOTR): ${OBJI} ${OBJC}

backup : 
	tar -cvf $(TAR) *.f90 Makefile
	#mv $(TAR) ../../backup
clean:
	-rm -f *.o core biotr_ang *.mod


.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I ${LMOD} -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

