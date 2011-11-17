.SUFFIXES: .f90 .mod

LIBDIR = $(ATSP)/lib
LIBA =  $(LIBDIR)/libcom90.a

MODDIR = ../libmod

#VPATH = ../libang90:../librad90:../mod

#FC90= xlf90
#FFLAGS = -g -qsuffix=f=f90

RM= /bin/rm -f

VASTO = ../libmod/vast_kind_param_M.o

MODOBJ = \
com_I.o  \
com_C.o

OBJ = \
acurat.o        dfill.o  dpbsl.o   intact.o         tinvit.o     \
bisect.o        dgbco.o  dscal.o   lineqn.o         tql2.o     \
bwint.o         dgbdi.o  dswap.o   lval.o           trbak1.o     \
ca.o            dgbfa.o  eptr.o    pack.o           trbak3.o     \
cb.o            dgbsl.o  factrl.o  pythag.o         tred1.o     \
dgeco.o         gracah.o  qsort.o          tred2.o     \
dgedi.o         hnorm.o       tred3.o     \
dgefa.o         hwf.o      \
dgesl.o         idamax.o     zcb.o     \
dasum.o         dnrm2.o  ifill.o     \
daxpy.o         dpbco.o  imtqlv.o  reform.o    \
dcopy.o         dpbdi.o  inita.o   rme.o    \
ddot.o          dpbfa.o  initr.o   symmeq.o    

install : $(LIBA)
$(LIBA) : $(OBJ)  
	@echo " Building install libang.a"
	ar -curs $(LIBA) $?

$(OBJ) : ${MODOBJ}
${MODOBJ} : ${VASTO}

clean :
	-rm -f *.o 

.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I ../libmod -o $@

SRC =  \
acurat.f90        dgbco.f90   eptr.f90           realloc_LINUX.f90    \
alloc_dec.f90     dgbdi.f90   factrl.f90         realloc_t3e.f90    \
alloc_IBMSP.f90   dgbfa.f90   fk.f90             reform.f90    \
alloc_LINUX.f90     dgbsl.f90   gk.f90             rme.f90    \
alloc_t3e.f90     dgeco.f90   gracah.f90         sn.f90    \
bisect.f90        dgedi.f90   hnorm.f90          symmeq.f90    \
bwint.f90         dgefa.f90   hwf.f90            tinvit.f90    \
ca.f90            dgemv.f90   idamax.f90         tql2.f90    \
cb.f90            dger.f90    ifill.f90          trbak1.f90    \
com_C.f90         dgesl.f90   imtqlv.f90         trbak3.f90    \
com_I.f90         dnrm2.f90   inita.f90          tred1.f90    \
dalloc_dec.f90    dpbco.f90   initr.f90          tred2.f90    \
dalloc_ibmSP.f90  dpbdi.f90   intact.f90         tred3.f90    \
dalloc_sun.f90    dpbfa.f90   lineqn.f90         uk.f90    \
dalloc_t3e.f90    dpbsl.f90   lval.f90           vast_kind_param_M.f90    \
dasum.f90         dscal.f90   pack.f90           vk.f90    \
daxpy.f90         dspevx.f90  pythag.f90         zcb.f90    \
dcopy.f90         dspmv.f90   qsort.f90    \
ddot.f90          dspr2.f90   realloc_dec.f90    \
dfill.f90         dswap.f90   realloc_ibmSP.f90    \

