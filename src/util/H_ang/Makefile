################################################################################
#
# Read the lists produced by nonh and show the angularly integrated 
# Hamiltonian to a human
#
# Written by Thomas Carette
#                   November, 2011
#
################################################################################
# makefile for show_H_ang
# Targets:
# . default: default target is show_H_ang
# . clean:   remove object files
################################################################################



EXED = ${ATSP}/bin
#EXED = .

#libraries
LIBDIR = ${ATSP}/lib
LIBUSE = -lang -lcom

SN=show_H_ang

#obkjects
OBJECT  = show_H_ang.o genint.o orbitals.o analy_blk.o ltriang.o #\
#        cfgblk.o
#        isearch.o ltriang.o  savels.o ffgg.o 
SRCF  = $(foreach a, $(OBJECT), $(patsubst %.o,%.f,$(a)))

#making
install: $(SN)

$(SN): $(OBJECT)
	$(FC) -o $(SN) $(FC_LD) $(FFLAGS) -L$(LIBDIR) $(OBJECT) \
	$(LIBUSE) 

clean:
	-rm -f *.o $(SN)


.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

