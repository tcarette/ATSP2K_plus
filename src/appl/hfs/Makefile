################################################################################
#=======================================================================
#
#      A General Program for Computing Magnetic Dipole and
#      Electric Quadrupole Hyperfine Constants
#
#      This version of the program assumes orthogonal orbitals
#
#      P. Jonsson and C.G. Wahlstrom, Department of Physics,
#                                     Lund Institute of Technology
#                                     P.O. Box 118,S-221 00 Lund
#                                     Sweden
#      C. Froese Fischer,             Department of Computer Science
#                                     Vanderbilt University,
#                                     Nashville, TN 37235
#                                     USA
#      Comput. Phys. Commun. 74 (1993) 399
#
#      Modified for dynamic memory allocation by P. Jonsson, Nov 91
#
#      Modified to calculate electron densities at the nucleus and
#      level field shifts.                       P. Jonsson, Feb 92
#
#      Modified to calculate f-shells by G. Gaigalas  December 1997
#
#      This program should be linked with libang.a libcom.a librad.a
#                       libdudu.a librecls.a libsqlsf1.a libsqlsf2.a
#
#=======================================================================
################################################################################


################################################################################
# makefile for hfs
# Targets:
# . default: default target is hfs
# . clean:   remove object files 
# . backup: makes a tar file, the date is incorporate in filename
# . HFS.f: will concantenate all source files into HFS.f
################################################################################


EXED = ${ATSP}/bin

# Directory where the libraries reside
LIBDIR = ${ATSP}/lib
EXEF = $(EXED)/hfs
TAR = hfs_`date '+%y%m%d_%H'`.bak.tar   # tar -file : year_month_date_hour

# Libraries to be searched
RCA = -lang -lcom -lrad

# Object files needed for the system
OBJF =  hfs.o readwt.o readwfn.o blkdta000.o lsjfact.o lsj.o qquadr.o \
        multwt.o radial1.o radial2.o orthog.o setup.o setupm.o cdens.o cfgn1.o \
        nonhiper.o orbit.o dipol.o contac.o qdrpol.o volum.o savenon.o

SRCF =  hfs.f readwt.f readwfn.f blkdta000.f lsjfact.f lsj.f qquadr.f \
        multwt.f radial1.f radial2.f orthog.f setup.f setupm.f cdens.f cfgn1.f \
        nonhiper.f orbit.f dipol.f contac.f qdrpol.f volum.f savenon.f

# making
install: $(EXEF)
$(EXEF) : $(OBJF)
	$(FC) -o $(EXEF) $(FC_LD) -L$(LIBDIR) $(OBJF)   $(RCA)
backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup
clean:
	-rm -f *.o core hfs


HFS.f:	$(SRCF)
	cat $(SRCF) >HFS.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

