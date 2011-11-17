EXED = ${ATSP}/bin

TAR = `date '+%y_%m_%d:%H'`_bpci_H.tar   # tar -file : year_month_date_hour
# executable targets
EXE = $(EXED)/bp_mat_mpi

# libraries
#LIBDIR = /u1/georgio/atsp2K/lib  #../../../lib
#LIBDIR = /u2/georgio/SPII/atsp2K/lib #../../../lib
LIBDIR = ${ATSP}/lib

RAD_COM = -lrad  -lcom -lrad
SPARSE = -ldvdson -lnet
LIB = -lang -lrad -lcom -lrad -lmpiu

# for debugging set FFLAGS to -g
#FC_FLAGS = -g -dp -e0

#objects
OBJ = mpi_bpciH_${FC_MALLOC}.o mpi_brevalr.o mpi_genintbrH.o mpi_genlst.o \
          ltriang.o input_bpH.o mpi_readw.o  alctab.o mpi_initm.o

install : $(EXE)
$(EXE): $(OBJ)
	$(FC_MPI) -o $(EXE) $(FC_MPILD) ${OBJ} \
              -L$(LIBDIR) $(LIB)  
archive:
	tar -cvf $(TAR) *.f Makefile

clean:
	rm -rf *.o core
SRC = mpi_bpciH.f mpi_brevalr.f mpi_genintbrH.f mpi_genlst.f  ltriang.f \
          input_bpH.f mpi_readw.f  alctab.f mpi_initm.f


bpciH.f : $(SRC) 
	cat $(SRC) > BPci.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I . -o $@
.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -o $@

