################################################################################
#     ==================================================================
#
#       A PROGRAM TO FIND SELECTED EIGENVALUES AND EIGENVECTORS
#       GIVEN FILES FOR THE INTERACTION MATRIX OF NON-RELATIVISTIC
#             THE BREIT-PAULI HAMILTONIANS
#
#       by C. Froese Fischer and G. Tachiev
#          Vanderbilt University
#          Nashville, TN 37235 USA
#          May, 2000
#
#                C O P Y R I G H T   2000
#     ==================================================================
#
#       The PARAMETER values in this program define the following:
#               NOD   - Maximum number of points in the range of a
#                         - function
#               LIMD  - Maximum number of basis vectors
#
#       Variables for these dimensions are
#               NWF   - Number of functions (or electrons)
#               NO    - Number of points in the range of the function
#               NCFG  - Number of configuration state functions
#               NUME  - Number of eigenvalues/eigenvectors
#               NZE   - Number of non-zero matrix elements
#
#     ------------------------------------------------------------------
################################################################################



################################################################################
# makefile for bp_eiv
# Targets:
# . default: default target is biotr
# . clean:   remove object files
# . backup: makes a tar file, the date is incorporate in filename
# . BP_EIV.f: will concantenate all source files into BP_EIV.f
################################################################################

EXED = ${ATSP}/bin

TAR = bp_eiv_`date '+%y%m%d_%H'`.bak.tar 
# executable targets
EXE = $(EXED)/bp_eiv

# libraries
LIBDIR = ${ATSP}/lib
RAD_COM = -lrad  -lcom -lrad 
SPARSE = -ldvdson 
LIB = -ldvdson -lang -lrad -lcom  

#objects
SYSO = p_alloc_${FC_MALLOC}.o
OBJ = ${SYSO} bp_eiv.o brci_eig.o  advance.o alcmat.o \
      flsj_comp.o brci_zeeman.o jmem.o lsjmat.o inp_eig.o \
      hmx_lsj.o inp_case.o alctab.o dssbmv.o dinit.o \
      initm.o readw.o lsp_comp.o term_shift.o inp_data.o

install : $(EXE)
$(EXE): $(OBJ)
	$(FC) -o $(EXE) $(FC_LD)  \
                $(OBJ) -L$(LIBDIR) $(LIB) ${LAPACK_DIR} ${LAPACK_LIB}
backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core bp_eiv

SRC = bp_eiv.f bpci.f  advance.f alcmat.f \
      flsj_comp.f brci_zeeman.f jmem.f lsjmat.f inp_eig.f \
      hmx_lsj.f inp_case.f alctab.f dssbmv.f dinit.f \
      initm.f readw.f lsp_comp.f term_shift.f inp_data.f


BP_EIV.f : $(SRC) 
	cat $(SRC) > BPci.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

