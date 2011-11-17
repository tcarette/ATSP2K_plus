#FC = $(ATSP_F90)
LIBDIR = $(ATSP)/lib
LIBA = $(LIBDIR)/librad.a
LIBOBJ = $(OBJ)

OBJ =    dyk.o dzk.o ecore.o grad.o hl.o hlc.o \
         quadr.o quads.o rk.o rlshft.o tk.o ykf.o ykk.o \
         zeta.o zk.o zz.o

install : $(LIBA)
$(LIBA) : $(OBJ)
	@echo " Building install libang.a"
	ar -curs $(LIBA) $(LIBOBJ)
clean:
	-rm -f *.o

SRC =    dyk.f dzk.f ecore.f grad.f hl.f hlc.f \
         quadr.f quads.f rk.f rlshft.f tk.f ykf.f ykk.f \
         zeta.f zk.f zz.f

RAD.f : $(SRC)
	cat $(SRC) > RAD.f

%.o : %.f
	$(FC) -c $(FC_FLAGS) $< -o $@

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

