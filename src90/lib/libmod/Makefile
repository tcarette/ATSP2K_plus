.SUFFIXES: .f90 .mod

LIBDIR = $(ATSP)/lib
LIBA =  $(LIBDIR)/libmod90.a

RM= /bin/rm -f

VASTO=vast_kind_param_M.o

MODOBJ = \
parameters_biotr_C.o \
ang_C.o  com_C.o   mpi_I.o  rad_C.o  \
ang_I.o  com_I.o  rad_I.o  dvdson_I.o dvdson_C.o \
biotr_ang_C.o \
biotr_tr_C.o \
biotr_ang_I.o \
biotr_tr_I.o


install : $(LIBA)
$(LIBA) : ${VASTO} $(MODOBJ)
	@echo " Building install libmod90.a"
	ar -curs $(LIBA) $?

${MODOBJ} : ${VASTO}

clean :
	-rm -f *.o *.mod

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I ../libmod -o $@

SRC =  \
ang_C.f90  com_C.f90  dvdson_C.f90  mpi_I.f90  rad_C.f90 \
ang_I.f90  com_I.f90  dvdson_I.f90  net_I.f90  rad_I.f90 \
vast_kind_param_M.f90 \
biotr_ang_C.f90 \
biotr_ang_I.f90 \
biotr_tr_C.f90 \
biotr_tr_I.f90

