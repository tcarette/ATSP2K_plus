
################################################################################
#     ==================================================================
#
#       A CONFIGURATION INTERACTION PROGRAM EITHER NON-RELATIVISTIC
#             OR IN THE BREIT-PAULI APPROXIMATION
#
#       by C. Froese Fischer
#          Vanderbilt University
#          Nashville, TN 37235 USA
#          May, 1983
#
#       Modified August, 1992 to combine BREIT and the sparse
#       matrix Dvdson eigensolver.
#
#       Modified by Gediminas Gaigalas for new angular structure
#       Modified for Unsorted version
#       Modified for restarting [Sep 1994]
#       Modified by Gediminas for new angular libraries and orbit-orbit (1997)
#
#              C O P Y R I G H T   2000
#     ==================================================================
#
#       The PARAMETER values in this program define the following:
#               NOD   - Maximum number of points in the range of a
#                         - function
#               LIMD  - Maximum number of basis vectors
#               LCDIM - Number of array elements in a block
#                       of the direct access file
#
#       Variables for these dimensions are
#               NWF   - Number of functions (or electrons)
#               NO    - Number of points in the range of the function
#               NCFG  - Number of configuration state functions
#               NUME  - Number of eigenvalues/eigenvectors
#               NZE   - Number of non-zero matrix elements
#               leigen- logical array ntermdXntermd for requested eigenvalues
#     ------------------------------------------------------------------
################################################################################



################################################################################
# makefile for bpci
# Targets:
# . default: default target is brci
# . clean:   remove object files 
# . backup: makes a tar file, the date is incorporate in filename
# . BPCI.f: will concantenate all source files into BRCI.f
################################################################################



EXED = ${ATSP}/bin

TAR = bpci_`date '+%y%m%d_%H'`.bak.tar   # tar -file : year_month_date_hour
# executable targets
EXE  = $(EXED)/bpci
HDSK =  $(EXED)/bpci_d

# libraruies
LIBDIR = ${ATSP}/lib
LIB_ANG_COM = -ldvdson -lang -lrad -lcom -lrad  -lnet #${LAPACK}

# objects

OBRC = bpci.o lsjmat.o readw.o
OBRCA = $(OBRC) breval.o breitgg.o acnfig.o advance.o alcmat.o alctab.o \
        genintbr.o initm.o isearch.o intval.o ltriang.o orthog.o outlsj.o \
        savebr.o dssbmv2.o setupgg.o findcol.o dinit.o  brci_eig.o brci_zeeman.o
OBRC_d = bpci_d.o lsjmat.o readw.o
OBRCA_d= $(OBRC_d) breval.o breitgg.o acnfig.o advance.o alcmat_d.o alctab.o \
        genintbr.o initm.o isearch.o intval.o ltriang.o orthog.o outlsj.o \
        savebr.o dssbmv2.o setupgg.o findcol.o dinit.o  brci_eig.o brci_zeeman.o

# making
install: $(EXE) ${HDSK}

$(EXE) : $(OBRCA)
	$(FC) -o $(EXE) $(FC_LD) -L$(LIBDIR) $(OBRCA)   \
	$(LIB_ANG_COM) ${LAPACK_DIR} ${LAPACK_LIB}

$(HDSK) : $(OBRCA_d)
	$(FC) -o $(HDSK) $(FC_LD)  -L$(LIBDIR) $(OBRCA_d)   \
	$(LIB_ANG_COM) ${LAPACK_DIR} ${LAPACK_LIB}

#$(LIBDIR)/libang.a $(LIBDIR)/libcom.a
backup: 
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core brci

SBRC = \
       breval.f breitgg.f acnfig.f advance.f alcmat.f alctab.f  genintbr.f \
       initm.f isearch.f intval.f ltriang.f orthog.f outlsj.f savebr.f \
       dssbmv2.f setupgg.f findcol.f alcmat_d.f  brci_eig.f brci_zeeman.f

SBRCF = $(SBRC) brci_F.f lsjmat_F.f readw_F.f
BPCI.f : $(SRCF)
	cat $(SRCF) > BRCI.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

