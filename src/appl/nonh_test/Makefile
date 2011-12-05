################################################################################
# ======================================================================
#
#     GENERAL PROGRAM TO COMPUTE MATRIX ELEMENTS OF THE  NON-
#     RELATIVISTIC HAMILTONIAN UNDER THE FOLLOWING ASSUMPTIONS -
#         (1) LS-COUPLING
#         (2) ORTHO-NORMAL CONFIGURATION STATE FUNCTIONS
#         (3) ALL ORBITALS ARE ORTHOGONAL
#
#     WRITTEN BY -
#     C. FROESE FISCHER, DEP'T OF COMPUTER SCIENCE
#        VANDERBILT UNIVERISTY
#     FEBRUARY, 1994
#
#     MODIFIED:
#     AUGUST,   1994  (A Ynnerman for unsorted lists)
#     DECEMBER, 1995                                 ( f-sell included )
#               1997  (G. Gaigalas, Vilnius for new angular codes)
#     DECEMBER, 1998  (C. Froese Fischer and G. Tachiev - block version)
#
#     C O P Y R I G H T  2000
# ======================================================================
################################################################################


################################################################################
# makefile for nonh
# Targets:
# . default: default target is nonh 
# . clean:   remove object files
# . backup: makes a tar file, the date is incorporate in filename
# . NONH.f: will concantenate all source files into NONH.f
################################################################################
#
#   Features:  computes coefficients of H for:
#   . simultanoeus optimization of multiple terms     
#   . simultanoeus optimization of multiple eigenvalues 
#
################################################################################



EXED = ${ATSP}/bin
#EXED = .
TAR = nonh_`date '+%y%m%d_%H'`.bak.tar   # tar -file : year_month_date_hour
# executable targets
SN = $(EXED)/nonh
SNZ = $(EXED)/nonhz

#libraries
LIBDIR = ${ATSP}/lib
LIBUSE = -lcom -ldvdson -lang  -lrad -lcom 

#obkjects
OBJECT_SN  = nonh.o angmomg.o genint.o orbitals.o\
        isearch.o ltriang.o  savels.o ffgg.o analy_blk.o\
        cfgblk.o

OBJECT_SNZ  = nonhz.o angmomg.o genint.o orbitals.o\
        isearch.o ltriang.o  savels.o ffgg.o analy_blk.o\
        cfgblk.o

#making
install: $(SN) ${SNZ}

$(SN): $(OBJECT_SN)
	$(FC) -o $(SN) $(FC_LD) $(FFLAGS) -L$(LIBDIR) $(OBJECT_SN) \
	$(LIBUSE) 

$(SNZ): $(OBJECT_SNZ)
	$(FC) -o $(SNZ) $(FC_LD) $(FFLAGS) -L$(LIBDIR) $(OBJECT_SNZ) \
	$(LIBUSE)
clean:
	-rm -f *.o core nonh nonhz
backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup


SRCF = nonh.f nonhz.f angmomg.f genint.f orbitals.f\
        isearch.f ltriang.f  savels.f ffgg.f  analy_blk.f\
        cfgblk.f

NONH.f  : $(SRCF)
	cat $(SRCF) > NONH.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

