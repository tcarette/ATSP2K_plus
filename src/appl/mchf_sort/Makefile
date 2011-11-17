################################################################################
#
#     ------------------------------------------------------------------
#     3          P R O G R A M   L I S T I N G
#     ------------------------------------------------------------------
#
#
#     All comments in the program listing assume the radial function P
#     is the solution of an equation of the form
#
#      P" + ( 2Z/R - Y - L(L+1)/R**2 - E)P = X + T
#
#     where Y is called a potential function
#           X is called an exchange function, and
#           T includes contributions from off-diagonal energy parameter,
#             interactions between configurations, etc.
#
#     The program uses LOG(Z*R) as independent variable and
#                      P/DSQRT(R) as dependent variable.
#
#     Modified by:
#     G. Tachiev for simultaneous optimization of multiple blocks
#     August, 2000
#                  C O P Y R I G H T   2000
#     ------------------------------------------------------------------
#                   M A I N    P R O G R AM
#     ------------------------------------------------------------------
#
#       The MAIN program controls the overall calculation and  allows
#   a series of cases to be processed as one run.  Each  case  itself
#   may  consist  of  a  series of atoms or ions in an iso-electronic
#   sequence.  In each case, all but the initial  estimates  for  the
#   first  are  obtained  by  scaling  the previous results using the
#   scaling of Sec.  (7-2).   Mixing coefficients are left unchanged.
#
################################################################################



################################################################################
# makefile for mchf
# Targets:
# . default: default target is mchf
# . clean:   remove object files
# . backup: makes a tar file, the date is incorporate in filename
# . MCHF.f: will concantenate all source files into MCHF.f
################################################################################
#
#   Features: 
#   . simultanoeus optimization of multiple terms     
#   . simultanoeus optimization of multiple eigenvalues 
#   . dynmaic memory/disk management
#
################################################################################


EXED = ${ATSP}/bin
#EXED = ${HOME}/bin

TAR = mchf_`date '+%y%m%d_%H'`.bak.tar   # tar -file : year_month_d
#  executable targets
SM = $(EXED)/mchf_sort
CDSK = $(EXED)/mchf_C_sort
HCDSK = $(EXED)/mchf_HC_sort

# libraries
LIBDIR = ${ATSP}/lib  
#RAD_COM = -lrad -lcom 
RAD_COM = -lrad -lcom -lnet
SPARSE = -ldvdson 
LIBUSE = -lcom -ldvdson 
# linker flags

#objects
SYSOBJ =  alcsts_${FC_MALLOC}.o diag_allocate_${FC_MALLOC}.o \
          diag_deallocate_${FC_MALLOC}.o

COBJ = alcsts_${FC_MALLOC}_C_disk.o diag_allocate_${FC_MALLOC}.o \
          diag_deallocate_${FC_MALLOC}.o

HCOBJ = alcsts_${FC_MALLOC}_HC_disk.o diag_allocate_${FC_MALLOC}.o \
          diag_deallocate_${FC_MALLOC}.o

USPOBJ= mchf.o alloc_mem.o data.o de.o diag.o diff.o dinit.o        \
        dssbmv.o e.o eig_out.o eijset.o grange.o methd1.o      \
        nmrvs.o nodec.o orthog.o output.o potl.o quad.o        \
        scale.o scf.o search.o solve.o spintgrl.o      \
        summry.o unpacki.o update.o wavefn.o xch.o             \
        diag_hmx.o supdatc.o prprty.o print_memory.o hrestor.o \
        recalcen.o sortcsf.o

install: de.o $(SM) $(CDSK) $(HCDSK)
de.o : de.f
	$(FC) $(FC_FLAGS) -c de.f
$(SM):  $(USPOBJ) $(SYSOBJ)
	$(FC) -o  $(SM) $(FC_LD) -L$(LIBDIR) $(USPOBJ) $(SYSOBJ)\
         $(SPARSE) $(RAD_COM) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(CDSK):  $(USPOBJ) $(COBJ)
	$(FC) -o $(CDSK) $(FC_LD) -L$(LIBDIR) $(USPOBJ) ${COBJ} \
	$(SPARSE) $(RAD_COM) -L ${LAPACK_DIR} ${LAPACK_LIB}

$(HCDSK):  $(USPOBJ) $(HCOBJ)
	$(FC) -o $(HCDSK) $(FC_LD)  -L$(LIBDIR) $(USPOBJ) ${HCOBJ} \
	$(SPARSE) $(RAD_COM) -L ${LAPACK_DIR} ${LAPACK_LIB}


clean:
	-rm -f *.o mchf

backup:
	tar -cvf $(TAR) *.f Makefile
	mv $(TAR) ../../backup

SRCF = mchf.f alloc_mem.f data.f de.f diag.f diff.f dinit.f        \
        dssbmv.f e.f eig_out.f eijset.f grange.f methd1.f      \
        nmrvs.f nodec.f orthog.f output.f potl.f quad.f        \
        scale.f scf.f search.f solve.f spintgrl.f      \
        summry.f unpacki.f update.f wavefn.f xch.f             \
        diag_hmx.f supdatc.f prprty.f print_memory.f alcsts_${FC_MALLOC}.f \
        diag_allocate_${FC_MALLOC}.f diag_deallocate_${FC_MALLOC}.f \
        hrestor.f

MCHF.f  : $(SRCF)
	cat $(SRCF) > Gtrans.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC) -c $(FC_FLAGS) $<  -I . -o $@
.f.o:
	$(FC) -c $(FC_FLAGS) $< -o $@

