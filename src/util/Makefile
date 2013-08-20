# Flags to be added for compilation
NOFFLAGS = $(ARCH)

OBJ=hf.o lsgen.o rotate2.o
SRC=hf.f lsgen.f rotate2.f
OBJU = \
comp.o        plotw.o       tables.o      condens.o     hf.o         \
lsgenf.o      printw.o      w_format.o    gencl.o       levels.o      \
lslines.o     relabel.o     w_unformat.o  lines.o      \
lstr.o        select.o      rtotate2.o

SRCU = \
comp.f        plotw.f       tables.f      condens.f     hf.f         \
lsgenf.f      printw.f      w_format.f    levels.f      \
lslines.f     relabel.f     w_unformat.f  lines.f       \
lstr.f        select.f      rotate2.f

#SRC77 = gencl.f genclf.f
SRCALL = $(SRC) $(SRC1)

LIBDIR = $(ATSP)/lib
LIBUSE = -lcom -lrad -lang -lcom

install : show_H_ang rotate2 w_unformat w_format hf lsgen
# lsgen 
#          hf  comp condens  levels lines \
#          w_unformat gencl select lsgen
#          rotate2 lsgenf
#         plotw printw relabel w_format FLS_trends T_dependence \
#
#	cp save $(ATSP)/bin
#	cd tbl; make; cd ..
#	cd LS_trends; make; cd ..

show_H_ang:
	@echo "....entering H_ang";\
	cd H_ang/;\
		$(MAKE) $@ ;\
	echo "....leaving H_ang" ;\
	mv show_H_ang $(ATSP)/bin;\
	cd ..

hf: hf.o
	$(FC)  -o hf $(FC_FLAGS) hf.o
	mv hf $(ATSP)/bin

lsgen : lsgen.o
	$(FC) -o lsgen $(NOFFLAGS) -L$(LIBDIR) lsgen.f $(LIBUSE)
	mv lsgen $(ATSP)/bin

lsgenf : lsgenf.o
	$(FC) -o lsgenf $(NOFFLAGS) -L$(LIBDIR) lsgenf.f $(LIBUSE)
	mv lsgenf $(ATSP)/bin

comp: comp.o
	$(FC)  -o comp $(FC_FLAGS)  comp.o 
	mv comp $(ATSP)/bin

condens : condens.o 
	$(FC) -o condens $(FC_FLAGS) -L$(LIBDIR) condens.o  $(LIBUSE)
	mv  condens $(ATSP)/bin

gencl : gencl.o
	${FC} -o gencl $(FC_FLAGS) -L$(LIBDIR) gencl.f $(LIBUSE)
	mv gencl $(ATSP)/bin

FLS_trends : FLS_trends.o
	${FC} -o FLS_trends $(FC_FLAGS) FLS_trends.o
	mv FLS_trends $(ATSP)/bin

T_dependence : T_dependence.o
	${FC} -o T_dependence $(FC_FLAGS) T_dependence.o
	mv T_dependence $(ATSP)/bin

levels : levels.o
	$(FC) -o levels $(FC_FLAGS) -L$(LIBDIR) levels.o $(LIBUSE)
	mv levels $(ATSP)/bin

lines: lines.o
	$(FC)  -o lines $(FC_FLAGS) lines.o
	mv lines $(ATSP)/bin

printw: printw.o
	$(FC)  -o  printw $(FC_FLAGS) printw.o
	mv  printw $(ATSP)/bin

plotw : plotw.o
	$(FC)  -o  plotw $(FC_FLAGS) plotw.o
	mv  plotw $(ATSP)/bin

select: select.o 
	$(FC) -o select $(FC_FLAGS) -L$(LIBDIR) select.o  $(LIBUSE)
	mv  select $(ATSP)/bin

relabel : relabel.o
	$(FC) -o relabel $(FC_FLAGS)  relabel.o
	mv relabel $(ATSP)/bin

w_format : w_format.o
	$(FC) -o w_format $(FC_FLAGS)  w_format.o
	mv w_format $(ATSP)/bin

w_unformat : w_unformat.o
	$(FC) -o w_unformat $(FC_FLAGS)  w_unformat.o
	mv w_unformat $(ATSP)/bin

rotate2 : rotate2.o
	${FC} -o rotate2 $(FC_FLAGS) -L$(LIBDIR) rotate2.f $(LIBUSE)
	mv rotate2 $(ATSP)/bin

clean:
	-rm -f *.o
	@cd show_H_ang
	$(MAKE) $@
	@cd ..

gencl.o : gencl.f
	${FC} -c $(ARCH) gencl.f
lsgen.o : lsgen.f
	${FC} -c $(ARCH) lsgen.f

%.o : %.f
	$(FC) -c $(FC_FLAGS) $< -o $@

%.o : %.f90
	$(FC) -c $(FC_FLAGS) $< -o $@

%.o : %.cc
	$(CPP) -c $(CPP_FLAGS) $< -o $@
