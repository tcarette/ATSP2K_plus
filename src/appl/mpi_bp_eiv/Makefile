EXED = ${ATSP}/bin

TAR = `date '+%y_%m_%d:%H'`_bpci_E.tar   # tar -file : year_month_date_hour
# executable targets
EXE = $(EXED)/bp_eiv_mpi


# libraries
#LIBDIR = /u1/georgio/atsp2K/lib #../../../lib
#LIBDIR = /u2/georgio/SPII/atsp2K/lib #../../../lib
LIBDIR = ${ATSP}/lib

RAD_COM = -lrad  -lcom -lrad  -llapack -lblas -lpgftnrtl
SPARSE = -ldvdson 
LIB =  -ldvdson  -lang -lrad -lcom -lrad -lmpiu  

#objects
OBJ = brci_est.o brci_eig.o mpi_bpciE_${FC_MALLOC}.o  advance.o mpi_alcmat.o \
      flsj_comp.o brci_zeeman.o jmem.o mpi_lsjmat.o inp_eig.o \
      mpi_hmx_lsj.o inp_case.o alctab.o mpi_dssbmv.o dinit.o \
      mpi_initm.o mpi_readw.o lsp_comp.o term_shift.o inp_data.o

install : $(EXE)

$(EXE): $(OBJ)
	$(FC_MPI) -o $(EXE) $(FC_MPILD)   \
                $(OBJ) -L$(LIBDIR) $(LIB) ${LAPACK_DIR} ${LAPACK_LIB}

archive:
	tar -cvf $(TAR) *.f Makefile

clean:
	rm -rf *.o core

SRC = brci_eig.f mpi_bpci.f  advance.f mpi_alcmat.f \
      flsj_comp.f brci_zeeman.f jmem.f mpi_lsjmat.f inp_eig.f \
      mpi_hmx_lsj.f mpi_lsjmat.f inp_case.f alctab.f mpi_dssbmv.f dinit.f \
      mpi_initm.f mpi_readw.f lsp_comp.f term_shift.f inp_data.f

bpciE.f : $(SRC) 
	cat $(SRC) > BPci.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I . -o $@
.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -o $@

