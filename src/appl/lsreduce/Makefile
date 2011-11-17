################################################################################
# ======================================================================
#
#     GENERAL PROGRAM TO DETERMINE WHICH CSFs IN CLIST DO
#     INTERACT WITH THE MULTIREFERENCE IN MRLIST
#
#     WRITTEN BY -
#
#     P. JONSSON, DEP'T OF COMPUTER SCIENCE
#        VANDERBILT UNIVERISTY, JULY 1995
#
#     THE PROGRAM IS A SMALL MODIFICATION OF THE GDUNONH
#     PROGRAM BY
#
#     G. GAIGALAS, INSTITUTE OF THEOERETICAL PHYSICS
#        AND ASTRONOMY, VILNIUS, LITHUANIA
#
#     C. FROESE FISCHER, DEP'T OF COMPUTER SCIENCE
#        VANDERBILT UNIVERISTY
#
#      ON THE INPUT.
#      mrlist;    CONTAINS THE MULTIREFERENCE CONFIGURATIONS
#                 WITH ANGULAR COUPLINGS
#      clist:     CONTAINS THE CFSs FROM LSGEN THAT SHOULD BE ANALYZED
#
#      ON THE OUTPUT
#      clist.out; CONTAINS 1) THE MULTIREFERENCE
#                          2) THE CSFs FROM clist THAT INTERACTS WITH THE
#                             MULTIREFERENCE
#
#     C O P Y R I G H T  2000
# ======================================================================
#
################################################################################



################################################################################
# makefile for lsreduce 
# Targets:
# . default: default target is nonh
# . clean:   remove object files
# . archive: makes a tar file, the date is incorporate in filename
# . LSREDUCE.f: will concantenate all source files into LSREDUCE.f
################################################################################

EXED = ${ATSP}/bin

TAR = lsredu_`date '+%y%m%d_%H'`.bak.tar
# executable targets
EXE = $(EXED)/lsreduce

#libraries
LIBDIR = ${ATSP}/lib

LIB = -lcom -ldvdson -lang  -lrad -lcom

OBJ = lsreduce.o analy_blk.o  cfgblk.o  intcfg.o   merge.o     savels.o \
      angmomg.o    ffgg.o    isearch.o  nonrelat.o  vijout.o      \
      cfg_inf.o    genint.o  ltriang.o  orbitals.o
install : ${EXE}
 
${EXE} : $(OBJ) 
	$(FC) -o ${EXE} $(FC_LD) -L$(LIBDIR) $(OBJ) $(LIB)

clean:
	-rm -f *.o core


backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup/


SRC = lsreduce.f analy_blk.f  cfgblk.f  intcfg.f    ltriang.f   orbitals.f \
      angmomg.f    ffgg.f    isearch.f   merge.f     savels.f \
      cfg_inf.f    genint.f  lsreduce.f  nonrelat.f  vijout.f

LSREDUCE.f :  $(SRC)
	  cat $(SRC) > LSREDUCE.f
.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

