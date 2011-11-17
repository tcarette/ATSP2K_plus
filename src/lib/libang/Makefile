#gt

LIBDIR = $(ATSP)/lib
LIBA = $(LIBDIR)/libang.a
LIBOBJ = $(ANG1OBJ) $(ANG2OBJ) $(OBJDUDU) $(OBJNORE) $(OBJRECLS) \
         $(OBJF1) $(OBJF2)

# object files for ANG1: no pointer variables
ANG1OBJ = initt.o ntab1.o tritst.o vijout.o
# obj files for ANG2: contains pointer variables
ANG2OBJ = analy1.o cfgo1.o cfgtst.o
OBJDUDU = \
   oneparticle1.o   oneparticle2.o   soo1.o           soo1111.o \
   soo1111p.o       soo1122.o        soo1122p.o       soo1212.o \
   soo1212p.o       sooa.o           soob.o           sooc.o \
   spinor.o         ss1.o            ss1111.o         ss1122.o \
   ss1212.o         ss1221.o         ssa.o            ssc.o \
   two1.o           two11.o          two12.o          two13.o \
   two2.o           two2a.o          two31.o          two32.o \
   two33a.o         two33b.o         two41.o          two42.o \
   two51.o          two52.o          two53.o          two54.o \
   two55.o          two56.o          twoparticle3.o   twoparticle31.o \
   twoparticle32.o  twoparticle33.o  twoparticle4.o   twoparticle5.o

OBJNORE = \
   calculation.o coulombls.o  eile.o       jfaze.o        nonrelat.o \
   nonrelat1.o   nonrelat2.o  nonrelat3.o  nonrelat31.o   nonrelat32.o \
   nonrelat33.o  nonrelat4.o  nonrelat41.o nonrelat5.o    nonrelat51.o \
   nonrelat52.o  nonrelat53.o orbitorbit.o

OBJRECLS = \
   dlsa1.o       dlsa2.o       dlsa3.o       dlsa4.o       dlsa5.o \
   dlsa6.o       recoupls0.o   recoupls2.o   recoupls3.o   recoupls31.o \
   recoupls4.o   rlsp0.o       rlsp00.o      rlsp1.o       rlsp2.o \
   rlsp3.o       rlsp31.o      rlsp32.o      rlsp4a.o      rlsp4b.o

OBJF1 = \
   a1.o           a1a2a3a4ls.o   a1a2a3a4lsp.o  a1a2ls.o       a1a2w3ls.o \
   a1a2w3lsp.o    a1aw2ls.o      a1aw2lsp.o     aconst.o       awp1g.o \
   awp1ls.o       c0t5s.o        c1e0sm.o       c1e1sm.o       cle0sm.o \
   glcons.o       gracah1.o      hibff.o        hibp31.o       itls.o \
   itls2.o        itls3.o        itrexg.o       itrexg2.o      ittk.o \
   ixjtik.o       izas1.o        jthn.o         mes.o          nine.o \
   nine0.o        nine11.o       nine12.o       nine13.o       ninels.o \
   niness.o       numter.o       numval.o       ribls.o        rmeals.o \
   rmewd1ls.o     rmewd2ls.o     rmewd3ls.o     rmewdls.o      rmewpls.o \
   rumt.o         rwls.o         sixj.o         sixj1.o        sixj2.o \
   sixj3.o        sixj35.o       sixj4.o        sixj5.o        sls.o \
   subls.o        termls.o       w1.o           w1g.o          w1w2ls.o \
   w1w2lsp.o      wa1a2ls.o      wa1a2lsp.o     wap1g.o        wap1ls.o \
   wwls1.o        wwpls1.o

OBJF2 = \
   avera.o    averd.o    averf.o    averp.o    frma.o     frma01.o   frma02.o \
   frma03.o   frma04.o   frma05.o   frma06.o   frma07.o   frma08.o   frma09.o \
   frma10.o   frma11.o   frma12.o   frma13.o   frma14.o   frma15.o   frma16.o \
   frma17.o   frma18.o   frma19.o   numterf.o  rumt67.o   trmf.o     trmf15.o

install : $(LIBA)
$(LIBA) : $(LIBOBJ)
	@echo " Building install libang.a"
	ar -curs $(LIBA) $(LIBOBJ) 
clean:
	rm -f *.o core


# ANG
ANG1SRC = initt.f ntab1.f tritst.f vijout.f
ANG1.f : $(ANG1SRC)
	cat $(ANG1SRC) >ANG1.f
ANG2SRC = analy1.f cfgo1.f cfgtst.f
ANG2.f : $(ANG2SRC)
	cat $(ANG2SRC) >ANG2.f
# DUDU
DUDUSRC = \
oneparticle1.f oneparticle2.f soo1.f soo1111.f soo1111p.f soo1122.f soo1122p.f \
soo1212.f soo1212p.f sooa.f soob.f sooc.f spinor.f ss1.f ss1111.f ss1122.f \
ss1212.f ss1221.f ssa.f ssc.f two1.f two11.f two12.f two13.f two2.f two2a.f \
two31.f two32.f two33a.f two33b.f two41.f two42.f two51.f two52.f two53.f \
two54.f two55.f two56.f twoparticle3.f twoparticle31.f twoparticle32.f \
twoparticle33.f  twoparticle4.f   twoparticle5.f
DUDU.f   : $(DUDUSRC)
	cat  $(DUDUSRC) > DUDU.f
#NORE
NORESRC = \
calculation.f  coulombls.f    eile.f         jfaze.f        nonrelat.f \
nonrelat1.f    nonrelat2.f    nonrelat3.f    nonrelat31.f   nonrelat32.f \
nonrelat33.f   nonrelat4.f    nonrelat41.f   nonrelat5.f    nonrelat51.f \
nonrelat52.f   nonrelat53.f   orbitorbit.f
NORE.f   : $(NORESRC)
	cat  $(NORESRC) > NORE.f
#RECLS
RECLSSRC = \
dlsa1.f       dlsa2.f       dlsa3.f       dlsa4.f       dlsa5.f \
dlsa6.f       recoupls0.f   recoupls2.f   recoupls3.f   recoupls31.f \
recoupls4.f   rlsp0.f       rlsp00.f      rlsp1.f       rlsp2.f \
rlsp3.f       rlsp31.f      rlsp32.f      rlsp4a.f      rlsp4b.f 
RECLS.f   : $(RECLSSRC)
	cat  $(RECLSSRC) > RECLS.f
#SQLSF1SRC 
SQLSF1SRC = \
a1.f           a1a2a3a4ls.f   a1a2a3a4lsp.f  a1a2ls.f       a1a2w3ls.f   \
a1a2w3lsp.f    a1aw2ls.f      a1aw2lsp.f     aconst.f       awp1g.f   \
awp1ls.f       c0t5s.f        c1e0sm.f       c1e1sm.f       cle0sm.f   \
glcons.f       gracah1.f      hibff.f        hibp31.f       itls.f   \
itls2.f        itls3.f        itrexg.f       itrexg2.f      ittk.f   \
ixjtik.f       izas1.f        jthn.f         mes.f          nine.f   \
nine0.f        nine11.f       nine12.f       nine13.f       ninels.f   \
niness.f       numter.f       numval.f       old_nine.f     old_nine0.f   \
old_nine11.f   old_nine12.f   ribls.f        rmeals.f       rmewd1ls.f   \
rmewd2ls.f     rmewd3ls.f     rmewdls.f      rmewpls.f      rumt.f   \
rwls.f         sixj.f         sixj1.f        sixj2.f        sixj3.f   \
sixj35.f       sixj4.f        sixj5.f        sls.f          subls.f   \
termls.f       w1.f           w1g.f          w1w2ls.f       w1w2lsp.f   \
wa1a2ls.f      wa1a2lsp.f     wap1g.f        wap1ls.f       wwls1.f \
wwpls1.f
SQLSF1.f   : $(SQLSF1SRC)
	cat  $(SQLSF1SRC) > SQLSF1.f
#SQLSF2SRC
SQLSF2SRC = \
avera.f    averd.f    averf.f    averp.f    frma.f     frma01.f   frma02.f \
frma03.f   frma04.f   frma05.f   frma06.f   frma07.f   frma08.f   frma09.f \
frma10.f   frma11.f   frma12.f   frma13.f   frma14.f   frma15.f   frma16.f \
frma17.f   frma18.f   frma19.f   numterf.f  rumt67.f   trmf.f     trmf15.f
SQLSF2.f   : $(SQLSF2SRC)
	cat  $(SQLSF2SRC) > SQLSF2.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

