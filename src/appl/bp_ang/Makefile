################################################################################
#     ==================================================================
#
#      THIS PROGRAM EVALUATES THE BREIT-PAULI OPERATORS
#
#     ONE-ELECTRON OPERATOR
#     ELECTROSTATIC INTERACTION
#     SPIN-ORBIT INTERACTION
#     SPIN-OTHER-ORBIT INTERACTION
#     SPIN-SPIN INTERACTION
#
#       by A. Hibbert, R. Glass, and C. Froese Fischer
#       Comput. Phys. Commun. 64 (1991) 455
#
#       Modified for Gediminas g-angular structure
#       Modified for Unsorted version
#       Modified for restarting
#                               [Sep 1994]
#       Modified by G. Gaigalas for new libraries and orbit-orbit (1997)
#       Modified by G. Tachiev for sparse matrix data structures (2000)
#
#
#     ==================================================================
#
#       The PARAMETER values in this program define the following:
#               LSDIM - Number of coeffficients in a block
#               NTERMD- Maximum number of terms
#
#     ------------------------------------------------------------------
################################################################################


################################################################################
# makefile for bp_ang
# Targets:
#   default: default target is bp_ang
#   clean:   remove object files
#   backup: makes a tar file, the date is incorporate in filename
#   BP_ANG.f: will concantenate all source files into BP_ANG.f
################################################################################

# directory where the executable bp_ang will be stored
EXED = ${ATSP}/bin

# executable
EXE = $(EXED)/bp_ang

# tar file of sources + Makefile
TAR = bp_ang_`date '+%y%m%d_%H'`.bak.tar  

# libraries
LIBDIR = ${ATSP}/lib

# library use
LIB = -lang -lrad -lcom -lrad -lang -lcom

#objects
OBJ = bp_ang.o acnfig.o breitgg.o brevala.o genintbr_${FC_MALLOC}.o  \
        input_br.o isearch.o ltriang.o orthog.o savebr.o setupgg.o 

install : $(EXE)
$(EXE): $(OBJ)
	$(FC) -o $(EXE) $(FC_LD) -L$(LIBDIR) \
                $(OBJ) $(LIB) 
backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core bp_ang

SRC = bp_ang.f acnfig.f breitgg.f brevala.f genintbr.f  \
        input_br.f isearch.f ltriang.f orthog.f savebr.f setupgg.f 

BP_ANG.f : $(SRC) 
	cat $(SRC) > BPci.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

