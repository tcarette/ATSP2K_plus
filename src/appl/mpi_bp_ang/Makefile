
EXED = ${ATSP}/bin

TAR = `date '+%y_%m_%d:%H'`_bpci_H.tar   # tar -file : year_month_date_hour
# executable targets
EXE = $(EXED)/bp_ang_mpi

# libraries
LIBDIR = ${ATSP}/lib

RAD_COM = -lrad  -lcom -lrad
SPARSE = -ldvdson 
LIB = -lang -lrad -lcom -lrad -lmpiu

#objects
OBJ = acnfig.o mpi_breit_${FC_MALLOC}.o breitgg.o mpi_brevala.o mpi_genintbr.o  \
        input_br.o isearch.o ltriang.o orthog.o savebr.o setupgg.o mpi_cfgo1.o

install : $(EXE)
$(EXE): $(OBJ)
	$(FC_MPI) -o $(EXE) $(FC_MPILD) -L$(LIBDIR) \
                $(OBJ) $(LIB) 
archive:
	tar -cvf $(TAR) *.f Makefile

clean:
	rm -rf *.o core

SRC = acnfig.f mpi_breit.f breitgg.f mpi_brevala.f mpi_genintbr.f  \
        input_br.f isearch.f ltriang.f orthog.f savebr.f setupgg.f mpi_cfgo1.f

MPI_BREIT.f : $(SRC) 
	cat $(SRC) > BPci.f

.f90.o:
	@echo "Compiling source file $< into $@"
	$(FC_MPI) -c $(FC_MPIFLAGS) $<  -I . -o $@
.f.o:
	$(FC_MPI) -c $(FC_MPIFLAGS) $< -o $@

