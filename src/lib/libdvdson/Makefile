LIBDIR = $(ATSP)/lib
LIBA = $(LIBDIR)/libdvdson.a
LIBOBJ = $(ANG1OBJ) $(ANG2OBJ) $(OBJDUDU) $(OBJNORE) $(OBJRECLS) \
         $(OBJF1) $(OBJF2)

LIBOBJ = \
addabs.o   dgathr.o   dscatr.o   dssbmv.o   dvdrvr.o   dvdson.o   multbc.o \
newvec.o   orthnrm.o  ovflow.o   setup.o    tstsel.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building install libang.a"
	ar -curs $(LIBA) $(LIBOBJ)
clean: 
	-rm -f *.o


SRCS    = dvdrvr.f dvdson.f addabs.f dgathr.f dscatr.f \
          dssbmv.f multbc.f newvec.f orthnrm.f ovflow.f \
          setup.f tstsel.f
DVDSON.f : $(SRCS)
	cat $(SRCS) >DVDSON.f

%.o : %.f
	$(FC) -c $(FC_FLAGS) $< -o $@

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

