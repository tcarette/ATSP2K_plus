################################################################################
#     ==================================================================
#
#      Program to generate the interaction matrix files in sparse matrix
#      format given files of BREIT-PAULI angular integraion data.
#
#
#       by C. Froese Fischer and G. Tachiev
#          Vanderbilt University
#          Nashville, TN 37235 USA
#          May, 2000
#
#                  C O P Y R I G H T   2000
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
# makefile for bp_mat
# Targets:
# . default: default target is bp_mat
# . clean:   remove object files 
# . backup: makes a tar file, the date is incorporate in filename
# . BP_MAT.f: will concantenate all source files into BP_MAT.f
################################################################################


EXED = ${ATSP}/bin

TAR = bp_mat_`date '+%y%m%d_%H'`.bak.tar   # tar -file : year_month_date_hour
# executable targets
EXE = $(EXED)/bp_mat

# libraries
LIBDIR = ${ATSP}/lib
RAD_COM = -lrad  -lcom -lrad
SPARSE = -ldvdson -lnet
LIB = -lang -lrad -lcom -lrad

# for debugging set FFLAGS to -g
#FC_FLAGS = -g -dp -e0

#objects
OBJ = bp_mat.o bpevalf.o genintbrH.o genlst.o \
          ltriang.o input_bpH.o readw.o alctab.o initm.o

install : $(EXE)
$(EXE): $(OBJ)
	$(FC) -o $(EXE) $(FC_LD) -L$(LIBDIR) \
                $(OBJ) $(LIB) 
backup:
	tar -cvf $(TAR) *.f Make*
	mv $(TAR) ../../backup

SRC = bp_mat.f brevalr.f genintbrH.f genlst.f  ltriang.f \
          input_bpH.f readw.f  alctab.f initm.f
clean:
	-rm -f *.o core bp_mat



BP_MAT.f : $(SRC) 
	cat $(SRC) > BPci.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

