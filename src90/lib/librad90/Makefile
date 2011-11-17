.SUFFIXES: .f90 .mod

LIBDIR = $(ATSP)/lib
LIBA =  $(LIBDIR)/librad90.a

MODDIR = ../libmod

#VPATH = ../libang90:../librad90:../mod

#FC90= xlf90
#FFLAGS = -qsuffix=f=f90

RM= /bin/rm -f

VASTO = ../libmod/vast_kind_param_M.o

MODOBJ = \
rad_I.o  \
rad_C.o

OBJ = \
dyk.o  ecore.o  hl.o   quadr.o  rk.o      tk.o   ykk.o   zk.o \
dzk.o  grad.o   hlc.o  quads.o  rlshft.o  ykf.o  zeta.o 

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

