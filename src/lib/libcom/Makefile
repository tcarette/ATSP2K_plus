LIBDIR = $(ATSP)/lib
LIBA =  $(LIBDIR)/libcom.a
LIBOBJ = $(COMOBJ) $(DYNOBJ) $(SCIOBJ) $(BLASOBJ) $(LINPACK)
 
# object files for library 'com'

COMOBJ=  func.o acurat.o bwint.o ca.o cb.o eptr.o factrl.o fk.o gk.o \
	 gracah.o hnorm.o hwf.o inita.o initr.o intact.o lineqn.o \
	 lval.o pack.o qsort.o reform.o rme.o sn.o symmeq.o \
	 uk.o vk.o zcb.o dfill.o ifill.o
# object files which may be part of the SCILIB
SCIOBJ = bisect.o imtqlv.o pythag.o tinvit.o tql2.o trbak1.o trbak3.o\
         tred1.o tred2.o tred3.o 
BLASOBJ= dasum.o daxpy.o dcopy.o ddot.o dscal.o dswap.o dnrm2.o \
         idamax.o dfill.o ifill.o
LINPACK= dgbco.o dgbdi.o dgbfa.o dgbsl.o dgeco.o dgedi.o dgefa.o \
         dgesl.o dpbco.o dpbdi.o dpbfa.o dpbsl.o dgemv.o dspevx.o \
         dger.o dspr2.o dgemv.o dspmv.o
# object files for the dynamic memory allocation
DYNOBJ = alloc_$(FC_MALLOC).o dalloc_$(FC_MALLOC).o realloc_$(FC_MALLOC).o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building install" # $(LIBA)
	ar -curs $(LIBA) $(LIBOBJ)

clean :
	-rm -f *.o

SRC   = acurat.f bwint.f ca.f cb.f eptr.f factrl.f fk.f gk.f \
	gracah.f hnorm.f hwf.f inita.f initr.f intact.f lineqn.f \
	lval.f pack.f qsort.f reform.f rk.f rme.f sn.f symmeq.f vk.f\
	dfill.f ifill.f
COM.f	: $(SRC)
	cat $(SRC) > COM.f
SCISRC  = bisect.f tinvit.f tql2.f trbak1.f tred1.f tred2.f
BLASSRC    =  daxpy.f dcopy.f ddot.f dscal.f idamax.f
SCI.f: $(SCISRC) $(BLASSRC)
	cat $(SCISRC) > SCI.f
DYNSRC =  alloc.f dalloc.f realloc.f
DYN.f   : $(DYNSRC)
	cat $(DYNSRC) >DYN.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@ func.o

