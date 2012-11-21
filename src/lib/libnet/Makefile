LIBDIR = $(ATSP)/lib
LIBA = $(LIBDIR)/libnet.a
NOOPT = $(ARCH)    # NOTE!!!! whn compiling ?lsmchf.f make sure that
                   # there is NO OPTIMIZATION, you may need to do it
                   # explicitly, chech the compiler documentation 
LIBOBJ = \
	dlagtf.o  dlanst.o  dlarnv.o  dlasr.o   dopmtr.o  dstebz.o \
	ilaenv.o  dlae2.o   dlagts.o  dlapy2.o  dlartg.o  dlasrt.o \
	dorg2l.o  dstein.o  lsame.o   dlaebz.o   dlarf.o   dspmv.o\
	dlascl.o  dlassq.o  dorg2r.o  dsteqr.o dlaev2.o  dlansp.o \
	dlarfg.o  dlaset.o  dopgtr.o  dsptrd.o  dsterf.o  xerbla.o \
	ieeeck.o  dlaruv.o  dlamch.o  slamch.o  dger.o    dasum.o  \
        dspr2.o

install : $(LIBA) 
$(LIBA) : $(LIBOBJ)
	@echo " Building install libnet.a"
	ar -curs $(LIBA) $(LIBOBJ)
clean:
	-rm -f *.o

slamch.o: slamch.f
	$(FC) $(NOOPT) -c  $<
dlamch.o: dlamch.f
	$(FC) $(NOOPT) -c  $<

%.o : %.f
	$(FC) -c $(FC_FLAGS) $< -o $@

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

