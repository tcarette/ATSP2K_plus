################################################################################
#=======================================================================
#
#      A General Program for Computing the Density and the 
#      natural orbital representation.
#
#      This version of the program assumes orthogonal orbitals
#
#      O. Scharf and M.R. Godefroid  ULB Bruxelles
#      A. Boorgoo                    VUB Bruxelles
#      G. Gedeminas                  VURITPA Vilnius
#
#      This program should be linked with libang.a libcom.a 
#
#=======================================================================
################################################################################


################################################################################
# makefile for density
# Targets:
# . default: default target is density
# . clean:   remove object files 
# . backup: makes a tar file, the date is incorporate in filename
# . DENSITY.f: will concantenate all source files into DENSITY.f
################################################################################

#ATSP is set in install_density script, it is the atsp2k main directory 
EXED = ${ATSP}/bin

# Directory where the libraries reside
LIBDIR = ${ATSP}/lib

EXEF = $(EXED)/density

# Libraries to be searched
RCA = -lang -lcom 

# Object files needed for the system

OBJF =  density.o readwt.o readwfn.o blkdta000.o qquadr.o \
        cfgo1.o spin_angular_density.o unitelement.o savenon.o \
        setup.o cfgn1.o setupm.o
		

SRCF =  density.f readwt.f readwfn.f blkdta000.f qquadr.f \
        cfgo1.f spin_angular_density.f unitelement.f savenon.f \
        setup.f cfgn1.f setupm.f
		
# making
install: $(EXEF)
$(EXEF) : $(OBJF)
	$(FC) -o $(EXEF) $(FC_LD) $(OBJF) -L$(LAPACK_DIR)  -L$(LIBDIR)  $(RCA) $(LAPACK_LIB)
backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core density


DENSITY.f:	$(SRCF)
	cat $(SRCF) >DENSITY.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $(RCA) $(LAPACK_LIB) $< -o $@

