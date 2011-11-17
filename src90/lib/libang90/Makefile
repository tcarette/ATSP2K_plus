.SUFFIXES: .f90 .mod

LIBDIR = $(ATSP)/lib
LIBA =  $(LIBDIR)/libang90.a

MODDIR = ../libmod

#VPATH = ../libang90:../librad90:../mod

#FC90= xlf90
#FFLAGS = -g -qsuffix=f=f90

RM= /bin/rm -f

VASTO = ../libmod/vast_kind_param_M.o

MODOBJ = \
ang_I.o  \
ang_C.o

OBJ = \
a1.o           frma05.o    nonrelat1.o     rmewd3ls.o  two12.o        \
a1a2a3a4ls.o   frma06.o    nonrelat2.o     rmewdls.o   two13.o        \
a1a2a3a4lsp.o  frma07.o    nonrelat3.o     rmewpls.o   two2.o        \
a1a2ls.o       frma08.o    nonrelat31.o    rumt.o      two2a.o        \
a1a2w3ls.o     frma09.o    nonrelat32.o    rumt67.o    two31.o        \
a1a2w3lsp.o    frma10.o    nonrelat33.o    rwls.o      two32.o        \
a1aw2ls.o      frma11.o    nonrelat4.o     savenon.o   two33a.o        \
a1aw2lsp.o     frma12.o    nonrelat41.o    sixj.o      two33b.o        \
aconst.o       frma13.o    nonrelat5.o     sixj1.o     two41.o        \
analy1.o       frma14.o    nonrelat51.o    sixj2.o     two42.o        \
avera.o        frma15.o    nonrelat52.o    sixj3.o     two51.o        \
averd.o        frma16.o    nonrelat53.o    sixj35.o    two52.o        \
averf.o        frma17.o    ntab1.o         sixj4.o     two53.o        \
averp.o        frma18.o    numter.o        sixj5.o     two54.o        \
awp1g.o        frma19.o    numterf.o       sls.o       two55.o        \
awp1ls.o       gracah1.o   numval.o        soo1.o      two56.o        \
c0t5s.o        hibff.o     oneparticle1.o  soo1111.o   twoparticle3.o        \
c1e0sm.o       hibp31.o    oneparticle2.o  soo1111p.o  twoparticle31.o        \
c1e1sm.o       itls.o      orbitorbit.o    soo1122.o   twoparticle32.o        \
calculation.o  itls2.o     recoupls0.o     soo1122p.o  twoparticle33.o        \
ccc.o          itls3.o     recoupls2.o     soo1212.o   twoparticle4.o        \
cfgo1.o        itrexg.o    recoupls3.o     soo1212p.o  twoparticle5.o        \
cfgtst.o       itrexg2.o   recoupls31.o    sooa.o      vijout.o        \
cle0sm.o       ittk.o      recoupls4.o     soob.o      w1.o        \
coulombls.o    ixjtik.o    rlsp0.o         sooc.o      w1g.o        \
dlsa1.o        izas1.o     rlsp00.o        spinor.o    w1w2ls.o        \
dlsa2.o        jfaze.o     rlsp1.o         ss1.o       w1w2lsp.o        \
dlsa3.o        jthn.o      rlsp2.o         ss1111.o    wa1a2ls.o        \
dlsa4.o        mes.o       rlsp3.o         ss1122.o    wa1a2lsp.o        \
dlsa5.o        nine.o      rlsp31.o        ss1212.o    wap1g.o        \
dlsa6.o        nine0.o     rlsp32.o        ss1221.o    wap1ls.o        \
eile.o         nine11.o    rlsp4a.o        ssa.o       wwls1.o        \
frma.o         nine12.o    rlsp4b.o        ssc.o       wwpls1.o        \
frma01.o       nine13.o    rmeals.o        subls.o        \
frma02.o       ninels.o    rmetr.o         tritst.o        \
frma03.o       niness.o    rmewd1ls.o      two1.o        \
frma04.o       nonrelat.o  rmewd2ls.o      two11.o

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

