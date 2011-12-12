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

#objects. Priorities according to use statements (TODO: clean that)
OPRIO1 = safeio.o
OPRIO2 = mdata.o
OPRIO3 = show_H_ang.o crash.o
OBJECT  = $(OPRIO1) $(OPRIO2) $(OPRIO3)
SRCF  = $(foreach a, $(OBJECT), $(patsubst %.o,%.f90,$(a)))

#making
install: $(SN)

$(SN): $(OBJECT)
	$(FC) -o $(SN) $(FC_LD) $(FFLAGS) -L$(LIBDIR) $(OBJECT) \
	$(LIBUSE) 
	@mv $(SN) $(EXEC)

clean:
	-rm -f *.o $(SN)

mdata.o: data.f90 spintgrl.f90 unsparse.f90

%.o: %.f90
	$(FC) -c $(FC_F90) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

