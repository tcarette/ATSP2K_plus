################################################################################
# makefile for iso
# Targets:
# . default: default target is iso
# . clean:   remove object files 
# . backup: makes a tar file, the date is incorporate in filename
# . ISO.f: will concantenate all source files into ISO.f
################################################################################
#
#     This system assumes all orbitals are orthogonal
#     so that there are no overlap integrals.
################################################################################



EXED = ${ATSP}/bin


TAR = iso_`date '+%y%m%d_%H'`.bak.tar   # tar -file : year_month_date_hour
# executable targets
EXEF = $(EXED)/isotope

# Directory where the libraries reside
LIBDIR = ${ATSP}/lib
LIB_ANG_COM = -lang -lcom 

# Object files needed for the system
ISOB = angmomg.o cfgn1.o cfgo1.o ffgg.o
ISOBF = $(ISOB) iso.o readwt.o 

#making
install: $(EXEF) 
$(EXEF): $(ISOBF)
	${FC} -o $(EXEF) $(FC_LD) -L$(LIBDIR) $(ISOBF) $(LIB_ANG_COM)
backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup

clean:
	-rm -f *.o core iso


SRC = angmomg.f cfgn1.f cfgo1.f ffgg.f
SRCF = $(SRC) iso.f readwt.f

ISO.f :  $(SRCF) 
	  cat $(SRCF) > ISO.f
.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

