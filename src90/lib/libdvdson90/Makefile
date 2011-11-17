.SUFFIXES: .f90 .mod

LIBDIR = $(ATSP)/lib
LIBA =  $(LIBDIR)/libdvdson90.a

MODDIR = ../libmod

#VPATH = ../libang90:../librad90:../mod

FC90= xlf90
FFLAGS = -g -qsuffix=f=f90

RM= /bin/rm -f

VASTO = ../libmod/vast_kind_param_M.o

MODOBJ = \
rad_I.o  \
rad_C.o

OBJ = \
addabs.o  dscatr.o  dvdrvr.o  multbc.o  orthnrm.o  setup.o   \
dgathr.o  dssbmv.o  dvdson.o  newvec.o  ovflow.o   tstsel.o


install : $(LIBA)
$(LIBA) : $(OBJ)  
	@echo " Building install libang.a"
	ar -curs $(LIBA) $?

$(OBJ) : ${MODOBJ}
${MODOBJ} : ${VASTO}

clean :
	-rm -f *.o *.mod

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I ../libmod -o $@


	
SRC =  \

